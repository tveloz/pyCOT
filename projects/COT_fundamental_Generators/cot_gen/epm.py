"""
epm.py — EPM and ESPM computation via the fundamental adjacency graph.

Architecture
------------
All traversal state lives in ERC-index space (FundamentalGraph + DFSState).
No species-level closure is computed during the DFS.

  req(S)  = (∪ req_mask[i]) & ~(∪ prod_mask[i])   for i ∈ S
  prod(S) = ∪ prod_mask[i]                          for i ∈ S

These are maintained incrementally by FundamentalGraph.extend_state().

Mode-1 DFS (req != 0): extend the active ERC-set via:
  a) Complementarity: for each species s in req(S), add a minimal producer
     E_prod via comp_by_species[s].
  b) Synergy: for each (E_i ∈ S, E_j ∉ S) with syn (E_i,E_j)→E_k and
     prod(E_k) ∩ req(S) ≠ ∅, add E_j (and E_k via synergy closure).

  Connectivity is NOT checked: it is structurally guaranteed by both
  combination rules.  Option B scan is O(|S|) per state, not O(|ERCs|).

Mode-2 (req == 0, found SSM): extend outward via fundamental synergy:
  for each (E_i ∈ S, E_j ∉ S) add E_j → drops into Mode-1.
  Scan is again O(|S|) × syn-degree, not O(|ERCs|).

Order assignment:
  SSMs sorted by species-mask popcount.
  order(C) = max(order(sub)) + 1; order = 0 if no proper sub-SO (→ EPM).

Public API
----------
compute_epms(rn, ercs, hier, syn_result, comp_result, *, counters=None) -> EPMResult
compute_espm(rn, ercs, hier, syn_result, comp_result, epm_result,
             *, max_order=10, verbose=False, counters=None) -> ESPMResult
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field

from .fundamental_graph import FundamentalGraph, DFSState


# ---------------------------------------------------------------------------
# Result data types
# ---------------------------------------------------------------------------

@dataclass
class EPMResult:
    """
    Output of compute_epms.

    External fields (species bitmasks, compatible with run_network.py)
    ------------------------------------------------------------------
    single_epm_indices : list[int]   — ERC list positions that are single-ERC EPMs
    single_epm_masks   : list[int]   — species bitmasks of those EPMs
    multi_epm_masks    : list[int]   — multi-ERC EPMs (order-0 from DFS traversal)
    all_epm_masks      : list[int]   — single + multi, sorted
    leaf_masks         : list[int]   — dead-end closures (req≠0, no progress)
    stats              : dict        — traversal counters

    Internal state (reused by compute_espm)
    ----------------------------------------
    _graph       : FundamentalGraph
    _visited_sp  : set[int]         — species masks already explored
    _sp_to_state : dict[int, DFSState]  — sp_mask → DFSState for each found SSM
    _so_order    : dict[int, int]   — sp_mask → SO order (0 = EPM)
    """
    single_epm_indices: list[int]
    single_epm_masks:   list[int]
    multi_epm_masks:    list[int]
    all_epm_masks:      list[int]
    leaf_masks:         list[int]
    stats:              dict

    _graph:        object = field(default=None,            repr=False)
    _visited_sp:   set    = field(default_factory=set,     repr=False)
    _sp_to_state:  dict   = field(default_factory=dict,    repr=False)
    _so_order:     dict   = field(default_factory=dict,    repr=False)

    def __len__(self) -> int:
        return len(self.all_epm_masks)


@dataclass
class ESPMResult:
    """
    Output of compute_espm.

    epm_masks           : list[int]              — order-0 SOs (EPMs)
    espm_by_order       : dict[int, list[int]]   — order ≥ 1 → [species masks]
    all_so_masks        : list[int]              — all SOs (EPMs + ESPMs)
    leaf_masks_by_order : dict[int, list[int]]   — order → [dead-end species masks]
    stats_by_order      : dict[int, dict]        — order → traversal counters
    """
    epm_masks:           list[int]
    espm_by_order:       dict[int, list[int]]
    all_so_masks:        list[int]
    leaf_masks_by_order: dict[int, list[int]]
    stats_by_order:      dict[int, dict]

    def max_order(self) -> int:
        return max(self.espm_by_order.keys(), default=0)

    def total_espm(self) -> int:
        return sum(len(v) for v in self.espm_by_order.values())


# ---------------------------------------------------------------------------
# Shared bit-iteration helper
# ---------------------------------------------------------------------------

def _bits(mask: int):
    """
    Yield the index (bit position) of every set bit in `mask`, ascending.

    Uses the two's-complement trick: m & (-m) isolates the lowest set bit.
    .bit_length()-1 converts it to an index; m &= m-1 clears it.

    Example: mask=0b1010 → yields 1, then 3.
    """
    m = mask
    while m:
        lsb = m & (-m)
        yield lsb.bit_length() - 1
        m &= m - 1


# ---------------------------------------------------------------------------
# Optional post-hoc connectivity check (NOT used in the hot path)
# ---------------------------------------------------------------------------

def _is_connected(supp_q, prod_q, n_rxn: int, n_species: int, X: int) -> bool:
    """
    Return True if species set X is connected through the network's reactions.

    NOT used during normal EPM/ESPM computation — kept for optional verification.
    ---------------------------------------------------------------------------
    The traversal algorithm (Mode-1 DFS) only ever combines two closures when
    there is a fundamental synergy or fundamental complementarity edge between
    them.  Both relation types structurally guarantee that the combined set is
    connected:

      • Complementarity (E_prod → E_cons via species s): s is produced by E_prod
        and required by a reaction in E_cons.  That reaction participates in both
        species groups, directly linking E_prod to E_cons.

      • Synergy (E_i + E_j → E_k): there exists a reaction whose support spans
        both E_i and E_j, directly linking species from both groups.

    Connectivity is therefore guaranteed by the combination rules themselves.
    This function exists only as a post-hoc sanity check for auditing.

    Input
    -----
    supp_q    : tuple[int] — quotiented support bitmask per reaction
    prod_q    : tuple[int] — quotiented product bitmask per reaction
    n_rxn     : int        — number of reactions
    n_species : int        — total species count
    X         : int        — species bitmask to test

    Output
    ------
    bool — True if X is a single connected component.
    """
    sp_list = [s for s in range(n_species) if (X >> s) & 1]
    if len(sp_list) <= 1:
        return True

    adj: dict[int, set[int]] = {s: set() for s in sp_list}
    for r in range(n_rxn):
        sq = supp_q[r]
        if not sq or (sq & X) != sq:
            continue
        pq_r = prod_q[r]
        rxn_sp = [s for s in sp_list if ((sq >> s) & 1) or ((pq_r >> s) & 1)]
        for a_i, a in enumerate(rxn_sp):
            for b in rxn_sp[a_i + 1:]:
                adj[a].add(b)
                adj[b].add(a)

    start = sp_list[0]
    visited = {start}
    queue = [start]
    while queue:
        curr = queue.pop()
        for nb in adj[curr]:
            if nb not in visited:
                visited.add(nb)
                queue.append(nb)
    return len(visited) == len(sp_list)


# ---------------------------------------------------------------------------
# Stage 1: single-ERC EPMs (minimal persistent ERCs)
# ---------------------------------------------------------------------------

def _single_erc_epms(ercs, hier) -> tuple[list[int], list[int]]:
    """
    Find all single-ERC EPMs: the ⊆-minimal persistent ERCs.

    Background
    ----------
    A persistent ERC (P-ERC) has req_mask == 0: it is self-sustaining.
    A single-ERC EPM is a P-ERC with no smaller P-ERC strictly inside it.
    Larger P-ERCs have a sub-SO, so they are ESPMs of order ≥ 1, not EPMs.

    Algorithm
    ---------
    Iterate ERCs smallest-first (they are already sorted by species_mask).
    When a P-ERC with no P-ERC below it is found → it is a single-ERC EPM.
    Immediately disqualify all its ancestors (supersets) via hier.ancestors[i].

    Parameters
    ----------
    ercs : list[ERCData]    — sorted by species_mask (output of compute_ercs)
    hier : HierarchyData    — hier.ancestors[i] = frozenset of strict supersets

    Returns
    -------
    (epm_idx, epm_masks)
      epm_idx   : list[int] — positions in ercs list
      epm_masks : list[int] — corresponding species bitmasks
    """
    disqualified: set[int] = set()
    epm_idx: list[int] = []

    for i in range(len(ercs)):
        if not ercs[i].is_persistent():
            continue
        if i in disqualified:
            continue
        epm_idx.append(i)
        disqualified.update(hier.ancestors[i])

    epm_masks = [ercs[i].species_mask for i in epm_idx]
    return epm_idx, epm_masks


# ---------------------------------------------------------------------------
# Order assignment (post-traversal, on species bitmasks)
# ---------------------------------------------------------------------------

def _assign_orders(all_sp_masks: list[int], *, verbose: bool = False) -> dict[int, int]:
    """
    Assign SO order by species-mask containment.

    Sort by popcount ascending so sub-SOs are processed before supersets.
    order(C) = max order of proper sub-SOs + 1; order = 0 if no sub-SO.

    Complexity: O(n_ssm²) — can be slow for thousands of SSMs.

    Parameters
    ----------
    all_sp_masks : list[int] — species bitmasks of all found SSMs
    verbose      : bool — print progress every 10% (useful when n_ssm is large)

    Returns
    -------
    dict[int, int] — species bitmask → order
    """
    n = len(all_sp_masks)
    if verbose:
        print(f"  [Stage 3] Assigning orders over {n} SSMs  (O(n²) = {n*n:,} checks)...")
        t0 = time.perf_counter()
        report_every = max(1, n // 10)   # print at every ~10%

    sorted_masks = sorted(all_sp_masks, key=lambda m: bin(m).count('1'))
    so_order: dict[int, int] = {}
    for idx, cl in enumerate(sorted_masks):
        max_sub = -1
        for sub_cl, sub_ord in so_order.items():
            if sub_cl != cl and (sub_cl & cl) == sub_cl and sub_ord > max_sub:
                max_sub = sub_ord
        so_order[cl] = max_sub + 1
        if verbose and (idx + 1) % report_every == 0:
            pct = (idx + 1) / n * 100
            elapsed = time.perf_counter() - t0
            print(f"    {elapsed:6.1f}s  {pct:5.1f}%  ({idx+1}/{n} SSMs assigned)")

    if verbose:
        elapsed = time.perf_counter() - t0
        print(f"    {elapsed:6.1f}s  100.0%  order assignment complete")

    return so_order


# ---------------------------------------------------------------------------
# Mode-1 DFS in ERC-index space
# ---------------------------------------------------------------------------

def _mode1_dfs(
    seed_states: list[DFSState],
    g: FundamentalGraph,
    visited_sp: set,
    stats: dict,
    *,
    verbose: bool = False,
) -> tuple[list[DFSState], list[DFSState]]:
    """
    DFS from each seed state following fundamental comp + syn edges until SSM.

    State is a DFSState (ERC-index set + incremental req/prod/sp).
    Deduplication is by species mask (state.sp): same species set = same state.

    No closure computation, no reaction scans.
    Connectivity is NOT checked (structurally guaranteed by the combination rules).

    Option A (complementarity): O(|req_species| × avg_comp_degree) per state.
    Option B (synergy):         O(|erc_set|     × avg_syn_degree)  per state.
      — Note: |erc_set| replaces the old O(n_ercs) full scan.

    Parameters
    ----------
    seed_states : list[DFSState] — starting states (one per non-P ERC)
    g           : FundamentalGraph
    visited_sp  : set[int] — species masks already explored (modified in-place)
    stats       : dict — counters (modified in-place)
    verbose     : bool — print progress and events as the DFS runs

    Returns
    -------
    (ssm_states, leaf_states)
      ssm_states  : new SSMs found (DFSState list)
      leaf_states : dead ends — req≠0, no extension possible (DFSState list)
    """
    ssm_states:  list[DFSState] = []
    leaf_states: list[DFSState] = []

    # ── Canonical-ordering dominance tracking (per-closure Pareto frontier)
    # visited_sp answers "has ANY representative of this closure been
    # explored" (unchanged contract, still used by callers). But different
    # candidate ERCs can produce DFSStates sharing the same .sp with
    # DIFFERENT min_ext (Rule 1 threshold) and min_seed=min(erc_set) (Rule 2
    # threshold) — measured on real networks: 17-25% of branching states
    # have such a collision among their own raw candidates. Whichever
    # representative got explored first previously "won" by pure processing
    # order, silently making every later, possibly more-permissive sibling
    # for the same .sp redundant-looking (`if new_state.sp not in
    # visited_sp` would reject it) even when it could reach genuinely more.
    # explored_sig tracks, per .sp, every (min_ext, min_seed) signature that
    # has already had its children generated; a new candidate is skipped
    # only if some already-explored signature for the same .sp dominates it
    # (lower-or-equal min_ext AND higher-or-equal min_seed — strictly at
    # least as permissive on both canonical-ordering rules). Because this
    # is one continuous stack (the whole DFS is a single call), no special
    # batching/eviction is needed the way compute_espm's per-round version
    # requires: pushing every not-yet-dominated candidate and re-checking
    # dominance at pop time is sufficient and correct.
    explored_sig: dict[int, list[tuple[int, int]]] = {}

    def _is_dominated(sp: int, min_ext: int, min_seed: int) -> bool:
        for (e_ext, e_seed) in explored_sig.get(sp, ()):
            if e_ext <= min_ext and e_seed >= min_seed:
                return True
        return False

    def _sig(state: DFSState) -> tuple[int, int]:
        return (state.min_ext, min(state.erc_set))

    stack: list[DFSState] = [s for s in seed_states if not _is_dominated(s.sp, *_sig(s))]

    # ── Verbose setup ─────────────────────────────────────────────────────
    # Progress metric: explored / (explored + stack_size).
    # This is the fraction of CURRENTLY KNOWN work done *so far*.
    # It can drop when new children are pushed (stack grows), and momentarily
    # hit 100% when the stack empties between batches of children being pushed.
    # To avoid false "done" signals:  capped at 99% inside the loop; only the
    # final post-loop summary prints 100%.
    if verbose:
        n_syn_edges  = sum(len(v) for v in g.syn_from.values())
        n_comp_pairs = sum(len(v) for v in g.comp_by_species.values())
        print(f"  [DFS] seeds={len(seed_states)}  n_ercs={g.n}"
              f"  syn_edges={n_syn_edges}  comp_pairs={n_comp_pairs}")
        print(f"  [DFS] Legend:")
        print(f"          explored  = DFS states popped and processed (each is a unique ERC-set)")
        print(f"          pending   = states waiting on the stack (search frontier)")
        print(f"          stack↑/↓  = how pending changed since last report (+growing / -shrinking)")
        print(f"          SSMs      = semi-self-maintaining states found (req==0; EPM candidates)")
        print(f"          dead-ends = states where req≠0 and no comp/syn extension was possible")
        print(f"          comp-ext  = extensions via complementarity (added a producer for a needed species)")
        print(f"          syn-ext   = extensions via synergy (added a partner ERC implied by a syn triple)")
        print(f"          states/s  = processing rate; slowing toward 0 means each state is expensive")
        print(f"          % done    = explored/(explored+pending); NOT monotone — can drop when stack grows")
        print(f"          branching = explored - SSMs - dead-ends (states that pushed children but aren't EPMs)")
        print(f"          pushes    = comp-ext + syn-ext (total children pushed; one state can push many)")
        print(f"          pruned    = extensions skipped by canonical ordering (ordering rule 1 + synergy rule 2)")
        t0            = time.perf_counter()
        t_last        = t0
        TICK_S        = 3.0        # print at least every N seconds
        PCT_STEP      = 5.0        # also print at each 5% milestone
        last_pct      = -PCT_STEP  # force first print immediately
        last_stack_sz = len(stack) # track stack trend
        MAX_EVENTS    = 8          # show individual SSM/leaf lines up to this many
        n_ev_ssm      = 0
        n_ev_leaf     = 0
        ev_suppressed = False
        progress_count = 0         # how many progress lines have been printed
        DIST_EVERY    = 20         # print SSM length distribution every N progress lines

        def _progress_line(explored, stack_sz, t_now, *, final=False):
            nonlocal last_stack_sz, progress_count
            total_known = explored + stack_sz
            pct = explored / max(total_known, 1) * 100
            if not final:
                pct = min(pct, 99.0)   # never show 100% mid-loop
            elapsed   = t_now - t0
            rate      = explored / max(elapsed, 1e-6)
            ssm_ct    = stats.get('ssm_found', 0)
            leaf_ct   = stats.get('leaves_found', 0)
            comp_ct   = stats.get('comp_extensions', 0)
            syn_ct    = stats.get('syn_extensions', 0)
            branching = explored - ssm_ct - leaf_ct   # states that pushed children
            delta     = stack_sz - last_stack_sz
            trend     = f"{'↑' if delta >= 0 else '↓'}{delta:+d}"
            last_stack_sz  = stack_sz
            progress_count += 1
            label = "DONE " if final else "     "
            print(f"  [DFS {label}{elapsed:7.1f}s]  {pct:5.1f}% of known work"
                  f"  |  explored={explored:8d}  pending={stack_sz:6d} (stack {trend})"
                  f"  |  SSMs={ssm_ct:7d}  dead-ends={leaf_ct:8d}  branching={branching:8d}"
                  f"  |  pushes={comp_ct + syn_ct:8d} (comp={comp_ct:7d} syn={syn_ct:7d})"
                  f"  |  pruned={stats.get('canonical_pruned', 0):8d}  |  {rate:6.0f} states/s")
            # every DIST_EVERY lines, print length distributions for all three state types
            if final or (progress_count % DIST_EVERY == 0):
                _print_len_dists(stats)
            return pct

        def _print_len_dists(stats):
            def _dist_str(d):
                return "  ".join(f"len={k}:{v:7d}" for k, v in sorted(d.items()))

            by_ssm  = stats.get('ssm_by_length', {})
            by_dead = stats.get('dead_by_length', {})
            by_bran = stats.get('branching_by_length', {})
            req_blen = stats.get('branching_req_by_length', {})  # {erc_len: {req_cnt: count}}

            if by_ssm:
                print(f"  [DFS len-dist]  SSMs       by #ERCs:  {_dist_str(by_ssm)}")
            if by_dead:
                print(f"  [DFS len-dist]  Dead-ends  by #ERCs:  {_dist_str(by_dead)}")
            if by_bran:
                print(f"  [DFS len-dist]  Branching  by #ERCs:  {_dist_str(by_bran)}")
            if req_blen:
                # For each erc_set length, show the req-count distribution compactly
                parts = []
                for elen in sorted(req_blen.keys()):
                    req_d = req_blen[elen]
                    inner = ",".join(f"req{r}:{c}" for r, c in sorted(req_d.items()))
                    parts.append(f"len={elen}:[{inner}]")
                print(f"  [DFS len-dist]  Branching  req/len:   " + "  ".join(parts))

    # ── Main DFS loop ─────────────────────────────────────────────────────
    while stack:
        state = stack.pop()

        state_sig = _sig(state)
        if _is_dominated(state.sp, *state_sig):
            continue
        explored_sig.setdefault(state.sp, []).append(state_sig)
        visited_sp.add(state.sp)
        explored = stats.get('states_explored', 0) + 1
        stats['states_explored'] = explored

        # ── Verbose: periodic progress line ───────────────────────────
        if verbose:
            t_now = time.perf_counter()
            total_known = explored + len(stack)
            cur_pct = min(explored / max(total_known, 1) * 100, 99.0)
            if cur_pct >= last_pct + PCT_STEP or (t_now - t_last) >= TICK_S:
                last_pct = _progress_line(explored, len(stack), t_now)
                t_last = t_now

        if state.is_ssm:
            ssm_states.append(state)
            n_ssm = stats.get('ssm_found', 0) + 1
            stats['ssm_found'] = n_ssm
            # track distribution of generator lengths (number of ERCs in the SSM)
            n_erc_key = len(state.erc_set)
            by_len = stats.setdefault('ssm_by_length', {})
            by_len[n_erc_key] = by_len.get(n_erc_key, 0) + 1
            # ── Verbose: SSM event ─────────────────────────────────
            if verbose:
                n_ev_ssm += 1
                if n_ev_ssm <= MAX_EVENTS:
                    n_sp   = bin(state.sp).count('1')
                    n_prod = bin(state.prod).count('1')
                    n_erc  = len(state.erc_set)
                    elapsed = time.perf_counter() - t0
                    print(f"  [DFS   {elapsed:7.2f}s] ✓ SSM #{n_ssm:4d}"
                          f"  — {n_erc} ERCs, {n_sp} species present, {n_prod} species produced"
                          f"  req=0 (self-sustaining)  erc_indices={sorted(state.erc_set)}")
                elif not ev_suppressed:
                    print(f"  [DFS] (individual SSM/leaf lines suppressed after {MAX_EVENTS}"
                          f"; totals tracked in progress lines)")
                    ev_suppressed = True
            continue

        made_progress = False

        # ── Canonical ordering setup ───────────────────────────────────
        # min_seed: smallest ERC index in this state.  Any explicit extension
        # new_idx must satisfy new_idx >= state.min_ext (rule 1).  After
        # synergy closure, if any newly-added non-persistent ERC has index <
        # min_seed, we're building an SSM whose canonical seed is that smaller
        # ERC — prune this branch (rule 2).  Persistent ERCs are exempt from
        # rule 2 because they are never seeds (pre-seeded in visited_sp).
        # (Reuses state_sig's min(erc_set) computed above at pop-time — no
        # need to recompute the same O(|erc_set|) min() a second time.)
        min_seed = state_sig[1]

        # ── Option A: complementarity ──────────────────────────────────
        # For each species s still required by this state, find minimal
        # producer ERCs and extend with each candidate. Candidates are
        # deduplicated by target ERC index first (the same producer can be
        # proposed for several different req species bits — extend_state is
        # a pure function of (state, candidate), so re-deriving an already-
        # seen candidate's result is pure waste; measured at 32-62% of raw
        # candidate edges on real networks).
        comp_cand: dict[int, DFSState] = {}
        for s_bit in _bits(state.req):
            for prod_idx in g.comp_by_species.get(s_bit, []):
                # Rule 1: canonical ordering — only extend with ERC ≥ min_ext
                if prod_idx < state.min_ext:
                    stats['canonical_pruned'] = stats.get('canonical_pruned', 0) + 1
                    continue
                if (g.species_mask[prod_idx] & state.sp) == g.species_mask[prod_idx]:
                    continue  # E_{prod_idx} species already covered by state
                if prod_idx not in comp_cand:
                    comp_cand[prod_idx] = g.extend_state(state, prod_idx)
        for new_state in comp_cand.values():
            # Rule 2: if synergy closure pulled in a non-persistent ERC
            # smaller than min_seed, the canonical path is from that ERC.
            newly_added = new_state.erc_set - state.erc_set
            if any(i < min_seed and not g.is_persistent[i] for i in newly_added):
                stats['canonical_pruned'] = stats.get('canonical_pruned', 0) + 1
                continue
            if not _is_dominated(new_state.sp, *_sig(new_state)):
                stack.append(new_state)
                made_progress = True
                stats['comp_extensions'] = stats.get('comp_extensions', 0) + 1

        # ── Option B: synergy ──────────────────────────────────────────
        # Iterate only ERCs currently IN the state (not all n_ercs).
        # For each E_i ∈ state, look for outward synergy partners E_j ∉ state.
        # Same per-target dedup as Option A: the same partner can be reached
        # via several different E_i already in the state.
        syn_cand: dict[int, DFSState] = {}
        for erc_i in state.erc_set:
            for (j, k) in g.syn_from.get(erc_i, []):
                if j in state.erc_set:
                    continue  # E_j already active
                # Rule 1: canonical ordering — only extend with ERC ≥ min_ext
                if j < state.min_ext:
                    stats['canonical_pruned'] = stats.get('canonical_pruned', 0) + 1
                    continue
                if not (g.prod_mask[k] & state.req):
                    continue  # E_k's production doesn't satisfy any req
                if j not in syn_cand:
                    syn_cand[j] = g.extend_state(state, j)
        for new_state in syn_cand.values():
            # Rule 2: synergy-closure pruning (same as Option A)
            newly_added = new_state.erc_set - state.erc_set
            if any(i < min_seed and not g.is_persistent[i] for i in newly_added):
                stats['canonical_pruned'] = stats.get('canonical_pruned', 0) + 1
                continue
            if not _is_dominated(new_state.sp, *_sig(new_state)):
                stack.append(new_state)
                made_progress = True
                stats['syn_extensions'] = stats.get('syn_extensions', 0) + 1

        if made_progress:
            # track branching state: erc_set length + how many species still required
            n_erc_key = len(state.erc_set)
            n_req_key = bin(state.req).count('1')
            by_bran = stats.setdefault('branching_by_length', {})
            by_bran[n_erc_key] = by_bran.get(n_erc_key, 0) + 1
            req_blen = stats.setdefault('branching_req_by_length', {})
            req_d = req_blen.setdefault(n_erc_key, {})
            req_d[n_req_key] = req_d.get(n_req_key, 0) + 1
        else:
            leaf_states.append(state)
            n_leaf = stats.get('leaves_found', 0) + 1
            stats['leaves_found'] = n_leaf
            # track dead-end length
            n_erc_key = len(state.erc_set)
            by_dead = stats.setdefault('dead_by_length', {})
            by_dead[n_erc_key] = by_dead.get(n_erc_key, 0) + 1
            # ── Verbose: leaf event ────────────────────────────────
            if verbose:
                n_ev_leaf += 1
                if n_ev_leaf <= MAX_EVENTS and not ev_suppressed:
                    n_req  = bin(state.req).count('1')
                    n_sp   = bin(state.sp).count('1')
                    n_erc  = len(state.erc_set)
                    elapsed = time.perf_counter() - t0
                    print(f"  [DFS   {elapsed:7.2f}s] ✗ DEAD-END #{n_leaf:4d}"
                          f"  — {n_erc} ERCs, {n_sp} species present,"
                          f" still needs {n_req} species externally"
                          f"  (no comp-producer or syn-partner can satisfy req)"
                          f"  erc_indices={sorted(state.erc_set)}")

    # ── Verbose: final summary ────────────────────────────────────────────
    if verbose:
        t_end    = time.perf_counter()
        explored = stats.get('states_explored', 0)
        _progress_line(explored, 0, t_end, final=True)   # only place that shows 100%
        ssm_ct    = stats.get('ssm_found', 0)
        leaf_ct   = stats.get('leaves_found', 0)
        comp_ct   = stats.get('comp_extensions', 0)
        syn_ct    = stats.get('syn_extensions', 0)
        branching = explored - ssm_ct - leaf_ct
        total_ext = comp_ct + syn_ct
        pruned_ct = stats.get('canonical_pruned', 0)
        print(f"  [DFS SUMMARY]"
              f"  total_time={t_end - t0:.2f}s"
              f"  |  {explored} states explored"
              f"    = {ssm_ct} SSMs + {leaf_ct} dead-ends + {branching} branching"
              f"  |  {total_ext} total child-pushes ({comp_ct} comp, {syn_ct} syn)"
              f"  |  {pruned_ct} pruned by canonical ordering")
        _print_len_dists(stats)

    return ssm_states, leaf_states


# ---------------------------------------------------------------------------
# Public API — EPM
# ---------------------------------------------------------------------------

def compute_epms(
    rn,  # noqa: ARG001 — kept for API compatibility; no reaction data needed post-refactor
    ercs,
    hier,
    syn_result=None,
    comp_result=None,
    *,
    counters=None,
    verbose: bool = False,
) -> EPMResult:
    """
    Compute EPMs (Def 29) via the fundamental adjacency graph traversal.

    Stage 1: single-ERC EPMs — minimal persistent ERCs (no DFS needed).
    Stage 2: build FundamentalGraph, run Mode-1 DFS from each non-P ERC.
    Stage 3: assign orders; collect multi-ERC EPMs (order-0 SSMs).

    Parameters
    ----------
    rn          : RNData          — compiled reaction network (supp_q, prod_q, …)
    ercs        : list[ERCData]   — all ERCs sorted by species_mask
    hier        : HierarchyData   — ERC containment hierarchy
    syn_result  : SynergyResult   — fundamental synergies (required for Stage 2)
    comp_result : CompResult      — fundamental complementarities (required for Stage 2)
    counters    : optional Counters

    Returns
    -------
    EPMResult  (carries internal state for compute_espm to reuse)
    """
    # ── Stage 1: single-ERC EPMs ──────────────────────────────────────────
    print("  [Stage 1] Finding single-ERC EPMs...")
    single_idx, single_masks = _single_erc_epms(ercs, hier)
    single_sp_set = set(single_masks)

    if syn_result is None or comp_result is None:
        result = EPMResult(
            single_epm_indices=single_idx,
            single_epm_masks=single_masks,
            multi_epm_masks=[],
            all_epm_masks=sorted(single_sp_set),
            leaf_masks=[],
            stats={'note': 'adjacency skipped (no syn/comp result)'},
        )
        if counters is not None:
            counters.inc('epm.n_single', len(single_idx))
        return result

    # ── Stage 2: build graph and run Mode-1 DFS ───────────────────────────
    print("  [Stage 2] Building fundamental graph...")
    g = FundamentalGraph(ercs, hier, syn_result, comp_result)

    visited_sp: set[int] = set()
    stats1: dict = {}

    # Pre-visit single-ERC EPMs — they are already found; skip them in DFS.
    for sp in single_masks:
        visited_sp.add(sp)

    # Seed: one initial state per non-persistent ERC.
    non_p_seeds = [
        g.make_seed_state(i)
        for i in range(len(ercs))
        if not ercs[i].is_persistent()
    ]

    print("  [Stage 2] Running Mode-1 DFS...")
    ssm_states, leaf_states = _mode1_dfs(non_p_seeds, g, visited_sp, stats1, verbose=verbose)

    # ── Stage 3: order assignment ─────────────────────────────────────────
    all_sp_masks = list(single_sp_set) + [s.sp for s in ssm_states]
    so_order = _assign_orders(all_sp_masks, verbose=verbose)

    for sp in single_masks:
        if sp not in so_order:
            so_order[sp] = 0

    multi_epm_masks = [
        s.sp for s in ssm_states
        if so_order.get(s.sp, -1) == 0 and s.sp not in single_sp_set
    ]
    all_epm_masks = sorted(single_sp_set | set(multi_epm_masks))

    # Map species mask → DFSState for all found SSMs (needed by compute_espm
    # for Mode-2 traversal: we need erc_set to find outward synergy edges).
    sp_to_state: dict[int, DFSState] = {}
    for s in ssm_states:
        sp_to_state[s.sp] = s
    # Single-ERC EPMs: build their trivial states.
    for i in single_idx:
        sp = ercs[i].species_mask
        if sp not in sp_to_state:
            sp_to_state[sp] = g.make_seed_state(i)

    result = EPMResult(
        single_epm_indices=single_idx,
        single_epm_masks=single_masks,
        multi_epm_masks=multi_epm_masks,
        all_epm_masks=all_epm_masks,
        leaf_masks=[s.sp for s in leaf_states],
        stats=stats1,
        _graph=g,
        _visited_sp=visited_sp,
        _sp_to_state=sp_to_state,
        _so_order=so_order,
    )

    if counters is not None:
        counters.inc('epm.n_single', len(single_idx))
        counters.inc('epm.n_multi',  len(multi_epm_masks))
        counters.inc('epm.n_total',  len(all_epm_masks))
        counters.inc('epm.n_leaves', len(leaf_states))

    return result


# ---------------------------------------------------------------------------
# Public API — ESPM
# ---------------------------------------------------------------------------

def compute_espm(
    rn,  # noqa: ARG001 — kept for API compatibility; no reaction data needed post-refactor
    ercs,
    hier,
    syn_result,
    comp_result,
    epm_result: EPMResult,
    *,
    max_order: int = 10,
    counters=None,
    verbose: bool = False,
) -> ESPMResult:
    """
    Compute ESPMs (Def 31) by BFS extension of known SOs via Mode-2 + Mode-1.

    For each known SO at order k-1, find outward fundamental synergy edges
    (Mode 2) and run Mode-1 DFS from each extension.  New SSMs are assigned
    orders by sub-SO containment.

    Mode-2 scan: O(|erc_set| × syn_degree) per SO — not O(n_ercs).

    Shared state (visited_sp, so_order, sp_to_state) is reused from
    compute_epms to avoid re-exploring already-processed states.

    Candidate dedup + canonical-ordering fix (see inline comments in the
    round loop below): raw Mode-2 candidates are deduplicated per-so_state
    by target ERC (many edges propose the identical ERC — a pure function
    of (so_state, candidate), so re-deriving it is pure waste), and then
    merged across so_states by resulting closure, keeping the Pareto
    frontier of (min_ext, min_seed)-non-dominated representatives rather
    than an arbitrary single winner. The latter was a real, measured bug:
    different candidate ERCs that converge on the same closure can carry
    different canonical-ordering thresholds, and letting an arbitrary
    processing order pick one silently dropped legitimate ESPMs on some
    networks (confirmed via oracle + soundness validation — this fix only
    ever recovers previously-missed, genuinely closed/SSM species sets, it
    does not change what counts as valid).

    Parameters
    ----------
    rn          : RNData
    ercs        : list[ERCData]
    hier        : HierarchyData
    syn_result  : SynergyResult
    comp_result : CompResult
    epm_result  : EPMResult    — output of compute_epms
    max_order   : int          — stop after this order (safety cap)
    counters    : optional Counters
    verbose     : bool         — print per-round progress

    Returns
    -------
    ESPMResult
    """
    # ── Reuse or rebuild shared state ────────────────────────────────────
    if epm_result._graph is not None:
        g             = epm_result._graph
        visited_sp    = epm_result._visited_sp
        sp_to_state   = epm_result._sp_to_state
        so_order      = epm_result._so_order
    else:
        g = FundamentalGraph(ercs, hier, syn_result, comp_result)
        visited_sp  = set()
        sp_to_state = {}
        so_order    = {sp: 0 for sp in epm_result.all_epm_masks}
        for sp in epm_result.all_epm_masks:
            # Reconstruct trivial single-ERC states (best effort)
            for i, e in enumerate(ercs):
                if e.species_mask == sp and e.is_persistent():
                    sp_to_state[sp] = g.make_seed_state(i)
                    break

    espm_by_order:       dict[int, list[int]] = {}
    leaf_masks_by_order: dict[int, list[int]] = {}
    stats_by_order:      dict[int, dict]      = {}

    # ── Carry over higher-order SOs found during compute_epms ────────────
    for sp, ord_k in so_order.items():
        if ord_k >= 1:
            espm_by_order.setdefault(ord_k, []).append(sp)
    for k in espm_by_order:
        espm_by_order[k].sort()

    if verbose and espm_by_order:
        ph1 = "  ".join(f"o{k}:{len(v)}" for k, v in sorted(espm_by_order.items()))
        print(f"  [Phase-1 carry-over]: {ph1}")

    # ── BFS by order ──────────────────────────────────────────────────────
    # current_layer: DFSState objects for all SOs at the previous order.
    # At order 1, we extend from EPMs (order 0).
    current_layer: list[DFSState] = [
        sp_to_state[sp]
        for sp in epm_result.all_epm_masks
        if sp in sp_to_state
    ]

    for order in range(1, max_order + 1):
        if not current_layer:
            break

        stats_k: dict = {'mode2_seeds': 0}

        # ── Mode-2: find outward growth seeds from each current SO ───────
        # For each SO state: iterate its ERC-set (not all n_ercs!) and try
        # the three fundamental extension moves that can grow an already-SSM
        # module (req == 0, so nothing is "needed" — growth is exploratory):
        #
        #   (a) Horizontal — synergy: E_j outside the SO forms a fundamental
        #       synergy with some E_i already in the SO.
        #   (b) Horizontal — complementarity (consumer side): E_j outside the
        #       SO fundamentally requires a species the SO already produces.
        #       (Mode-1's Option A only ever attaches *producers* of an open
        #       requirement; an SSM state has no open requirement, so the
        #       inverse — attaching *consumers* of what is already produced —
        #       is the move that discovers e.g. M1/M2 in the worked example.)
        #   (c) Vertical lift: replace/absorb E_i already in the SO with a
        #       direct hierarchy ancestor E_a (E_a ⊋ E_i).  E_a already
        #       incorporates every reaction of E_i (and hence of the SO's
        #       other constituents built on E_i), so this is just another
        #       extend_state() call — no separate "absorb" bookkeeping is
        #       needed.  This discovers e.g. M3 in the worked example, which
        #       is unreachable via (a)/(b) because E_a is dominated by E_i as
        #       a producer/consumer and therefore never appears as a
        #       *fundamental* complementarity partner itself.
        # Candidate collection is deduplicated in two layers before any
        # closure-chasing happens:
        #   1. Per-so_state, per-candidate-ERC dedup (`local_cand`): the same
        #      target ERC is frequently reachable via several different edges
        #      (e.g. as a synergy partner of two different members of the
        #      same SO) — extend_state() is a pure function of (so_state,
        #      candidate), so calling it more than once per pair is wasted
        #      work. Measured on real networks: 32-62% of all raw candidate
        #      edges are exactly this kind of duplicate.
        #   2. Cross-so_state / cross-candidate convergence: different
        #      candidate ERCs (or the same target reached from different SOs
        #      in this round) frequently converge on the IDENTICAL resulting
        #      closure anyway — measured at 42-65% of the already-deduplicated
        #      candidates. But the resulting DFSStates are not interchangeable:
        #      min_ext = new_idx + 1 depends on which specific ERC produced
        #      the extension, and min_seed = min(erc_set) can differ too, so
        #      two states with the SAME .sp can differ in which future
        #      canonical-ordering rules (Rule 1 on min_ext, Rule 2 on
        #      min_seed) apply to them going forward. Naively keeping just
        #      one arbitrary representative per closure — which is what the
        #      pre-fix code effectively did, since _mode1_dfs dedups its
        #      shared stack purely by .sp — silently let an arbitrary
        #      processing order decide which candidate "won", each of which
        #      can prune a genuinely different (and sometimes non-overlapping)
        #      part of the reachable search space. Confirmed on e_coli_core:
        #      up to 8 distinct min_ext values reachable for a single closure
        #      in one round.
        #   Fix: keep the Pareto frontier of non-dominated (min_ext, min_seed)
        #   representatives per closure (a state dominates another sharing
        #   its .sp iff its min_ext is <= AND its min_seed is >=  — lower
        #   min_ext is always more permissive for Rule 1, higher min_seed is
        #   always more permissive for Rule 2). Singleton frontiers (the
        #   overwhelming majority) go through the normal batched _mode1_dfs
        #   call. Genuine ties are rare (~2-3% of distinct closures per
        #   round, measured) and are processed sequentially, each given a
        #   real chance to expand by temporarily evicting its .sp from
        #   visited_sp for its own turn — anything an earlier tied rep
        #   already found stays correctly deduplicated; only genuinely new
        #   descendants reachable due to THIS rep's own permissiveness get
        #   added on top. This only ever adds legitimately-reachable SSMs
        #   that arbitrary ordering was silently dropping before — it cannot
        #   remove anything a stricter ordering choice would have found.
        best_by_sp: dict[int, list[DFSState]] = {}
        for so_state in current_layer:
            local_cand: dict[int, DFSState] = {}
            for i in so_state.erc_set:
                # (a) synergy partners
                for (j, _) in g.syn_from.get(i, []):
                    if j in so_state.erc_set:
                        continue  # E_j already in SO
                    if j not in local_cand:
                        local_cand[j] = g.extend_state(so_state, j)

                # (c) vertical lift: direct hierarchy ancestors of E_i
                for a in g.parents[i]:
                    if a in so_state.erc_set:
                        continue
                    if (g.species_mask[a] & so_state.sp) == g.species_mask[a]:
                        continue  # E_a already fully covered
                    if a not in local_cand:
                        local_cand[a] = g.extend_state(so_state, a)

            # (b) complementarity consumers of species already produced
            for s_bit in _bits(so_state.prod):
                for cons_idx in g.comp_consumers_by_species.get(s_bit, []):
                    if cons_idx in so_state.erc_set:
                        continue
                    if (g.species_mask[cons_idx] & so_state.sp) == g.species_mask[cons_idx]:
                        continue  # already fully covered
                    if cons_idx not in local_cand:
                        local_cand[cons_idx] = g.extend_state(so_state, cons_idx)

            for ext_state in local_cand.values():
                if ext_state.sp in visited_sp:
                    continue
                frontier = best_by_sp.setdefault(ext_state.sp, [])
                new_seed = min(ext_state.erc_set)
                dominated = False
                survivors = []
                for cand in frontier:
                    cand_seed = min(cand.erc_set)
                    if cand.min_ext <= ext_state.min_ext and cand_seed >= new_seed:
                        dominated = True
                        survivors.append(cand)
                    elif ext_state.min_ext <= cand.min_ext and new_seed >= cand_seed:
                        continue  # cand is dominated by ext_state — drop it
                    else:
                        survivors.append(cand)  # neither dominates — keep both
                if not dominated:
                    survivors.append(ext_state)
                best_by_sp[ext_state.sp] = survivors

        singleton_seeds: list[DFSState] = []
        tied_groups: list[list[DFSState]] = []
        for frontier in best_by_sp.values():
            if len(frontier) == 1:
                singleton_seeds.append(frontier[0])
            else:
                tied_groups.append(frontier)

        if not singleton_seeds and not tied_groups:
            if verbose:
                print(f"  [Round {order:>2}]: no Mode-2 seeds from"
                      f" {len(current_layer)} SOs — done.")
            break

        stats_k['mode2_seeds'] = len(singleton_seeds) + sum(len(t) for t in tied_groups)
        stats_k['pareto_ties'] = len(tied_groups)

        # ── Mode-1 DFS from each Mode-2 extension ────────────────────────
        new_ssm_states: list[DFSState] = []
        new_leaf_states: list[DFSState] = []
        if singleton_seeds:
            ssms, leaves = _mode1_dfs(singleton_seeds, g, visited_sp, stats_k, verbose=verbose)
            new_ssm_states.extend(ssms)
            new_leaf_states.extend(leaves)
        for tied in tied_groups:
            for rep in tied:
                if rep.sp in visited_sp:
                    visited_sp.discard(rep.sp)
                ssms, leaves = _mode1_dfs([rep], g, visited_sp, stats_k, verbose=False)
                new_ssm_states.extend(ssms)
                new_leaf_states.extend(leaves)

        # ── Update so_order with new SSMs ────────────────────────────────
        for s in sorted(new_ssm_states, key=lambda st: bin(st.sp).count('1')):
            if s.sp not in so_order:
                max_sub = -1
                for sub_sp, sub_ord in so_order.items():
                    if sub_sp != s.sp and (sub_sp & s.sp) == sub_sp and sub_ord > max_sub:
                        max_sub = sub_ord
                so_order[s.sp] = max_sub + 1
            sp_to_state[s.sp] = s  # store for next Mode-2 round

        if new_leaf_states:
            leaf_masks_by_order.setdefault(order, []).extend(
                s.sp for s in new_leaf_states
            )
        stats_by_order[order] = stats_k

        # ── Verbose round summary ─────────────────────────────────────────
        if verbose:
            n_seeds  = stats_k.get('mode2_seeds', 0)
            n_cls    = stats_k.get('states_explored', 0)
            n_leaves = len(new_leaf_states)
            if new_ssm_states:
                ord_counts: dict[int, int] = {}
                for s in new_ssm_states:
                    ok = so_order.get(s.sp, -1)
                    ord_counts[ok] = ord_counts.get(ok, 0) + 1
                ord_str = "  ".join(f"o{k}:{c}" for k, c in sorted(ord_counts.items()))
            else:
                ord_str = "none"
            print(f"  [Round {order:>2}]: {n_seeds:4d} seeds → {len(new_ssm_states):4d}"
                  f" new SOs  [{ord_str}]  {n_leaves} leaves  {n_cls} states")

        current_layer = new_ssm_states
        if not new_ssm_states:
            break

    # ── Final: rebuild espm_by_order from so_order ───────────────────────
    espm_by_order.clear()
    for sp, ord_k in so_order.items():
        if ord_k >= 1:
            espm_by_order.setdefault(ord_k, []).append(sp)
    for k in list(espm_by_order):
        espm_by_order[k].sort()

    if counters is not None:
        for k, masks in espm_by_order.items():
            counters.inc(f'espm.n_order{k}', len(masks))

    return ESPMResult(
        epm_masks=sorted(epm_result.all_epm_masks),
        espm_by_order=espm_by_order,
        all_so_masks=sorted(so_order.keys()),
        leaf_masks_by_order=leaf_masks_by_order,
        stats_by_order=stats_by_order,
    )


# ---------------------------------------------------------------------------
# Public API — latent join (on demand, NOT searched by compute_epms/compute_espm)
# ---------------------------------------------------------------------------

def latent_join(rn, X: int, Y: int) -> int | None:
    """
    Reconstruct the join of two persistent modules X, Y on demand.

    compute_epms / compute_espm deliberately never search for pure "latent"
    joins: unions of already-persistent modules that share no synergy and no
    complementarity (e.g. two persistent ERCs that merely overlap in one
    shared, never-required species).  Lemma req_containment guarantees such a
    join is automatically SSM whenever X and Y are, so there is nothing to
    "discover" by search — it can always be recomputed for O(1) closures on
    request instead of being enumerated ahead of time.  This is what keeps
    the ESPM search confined to the generated core (Sosgen) rather than all
    of Sos, which is exponentially larger for disjoint/loosely-connected
    modules (see companion paper, Latent–generated factorization).

    Parameters
    ----------
    rn : RNData
    X, Y : int — species bitmasks of two already-known persistent modules

    Returns
    -------
    int | None
        The species bitmask of clos(X ∪ Y) if that closure is itself a
        genuine persistent module (closed, SSM, connected, reactive);
        None otherwise (e.g. X and Y are disconnected, or the join
        activates a reaction that leaves an unmet requirement).
    """
    from .closure import build_inv_idx, closure_opt, is_ssm

    inv_idx = build_inv_idx(rn.supp_q, rn.n_species)
    joined = closure_opt(rn.supp_q, rn.prod_q, X | Y, inv_idx)

    if not is_ssm(rn.supp_q, rn.prod_q, joined):
        return None
    if not _is_connected(rn.supp_q, rn.prod_q, rn.n_reactions, rn.n_species, joined):
        return None
    return joined
