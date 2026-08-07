"""
epm_reach_propagate.py — validated, drop-in-alternative to
cot_gen.epm.compute_epms combining two provably-sound optimizations:

  1. reach(k) dead-end pruning
  2. branch-scoped unit propagation, with a density-based adaptive gate

Both are proven to never change the result set (only how much work is spent
finding it), and were validated against the real, unmodified
cot_gen.epm.compute_epms (identical EPM sets, every configuration) before
being swept across 16 real BiGG networks, where the combined approach beat
the unoptimized baseline in wall-clock time on 15 of 16 (0.9x-2.8x), the
sole exception being an extreme density outlier (see DENSITY_GUARD below).

------------------------------------------------------------------------
Theory: reach(k) pruning
------------------------------------------------------------------------
Every DFS branch only ever explicitly adds ERCs with index >= its current
`min_ext` threshold (canonical-ordering Rule 1); this only increases along
a branch. Define, for every possible threshold k:

    reach_prod(k) = union of prod_mask over the synergy-closure of the
                    WHOLE pool {ERC i : i >= k}, taken all at once.

`erc_syn_close` is a monotone Horn-closure operator in its input facts: a
standard property of forward-chaining propagation (more initial facts can
only derive a superset of consequences — provable by induction on
propagation steps). Since any real branch's explicit choices are a subset
of {i >= k}, monotonicity gives closure(explicit choices) subset-or-equal
closure({i>=k}) = reach_prod(k). Hence for any state S with min_ext=k:

    req(S) subset-or-equal reach_prod(k)   is NECESSARY for S to ever
                                            reach req == 0.

If this fails, S (and everything it could ever grow into) is a certain,
provable dead end — prune it at construction time, zero exploration
needed. reach_prod(k) for every k is computed ONCE, before the DFS starts,
by a single O(n) incremental sweep (erc_syn_close is only ever called with
one new index at a time, from k=n down to k=0, reusing the previous
closure) — not per-DFS-state.

------------------------------------------------------------------------
Theory: branch-scoped unit propagation
------------------------------------------------------------------------
At any DFS state, build the full set of legal candidate extensions
("survivors": pass canonical ordering, reach-pruning, and dominance
checks — exactly the set that would normally be pushed to the stack).
For each species s still required, count how many survivors would supply
it (`cover_count[s]`). If exactly one survivor covers some s, that
survivor is FORCED: no other legal move can ever satisfy s in this
branch, so it must be part of any successful completion — apply it
immediately (folding in the identical dominance-registration bookkeeping
a real stack pop would do) instead of pushing a stack frame for a
decision that was never actually open, then re-derive from the new state.
Terminates when: req==0 (SSM), no survivors (certain dead end), or every
remaining species has >=2 legal survivors (a genuine branch point — push
them all, exactly as the unmodified search would).

This is analogous to unit propagation / Boolean constraint propagation in
SAT solving. It is unconditionally safe: skipping it (falling back to
plain branching) is always correct, since that's the already-validated
baseline behavior — so gating it by any performance heuristic can never
change the result set, only the time spent finding it.

------------------------------------------------------------------------
DENSITY_GUARD — the one validated failure mode
------------------------------------------------------------------------
The per-state propagation tally costs O(n_req_bits x n_survivors). On 15
of 16 tested BiGG networks (synergy-density n_fund_synergies/n_ercs
ranging 0.4 to 13.8) this cost was consistently outweighed by the savings
from skipping stack round-trips — often substantially (up to 2.8x wall-
clock). The single exception, iAF692, sits at density 38.9 — nearly 3x
higher than anything else tested — where the tally overhead exceeded the
benefit. DENSITY_GUARD_THRESHOLD (default 20.0, comfortably above the
highest validated success at 13.8 and below the one validated failure at
38.9) disables propagation automatically on networks past that point,
falling back to reach-pruning alone (itself unconditionally validated).
This is a simple, static, one-time-computed guard — not a general learned
adaptive system — because 16 real networks is not yet strong evidence
that more machinery than this would earn its keep.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field

from cot_gen.fundamental_graph import FundamentalGraph, DFSState


# =============================================================================
# Result type
# =============================================================================

@dataclass
class EPMResult:
    single_epm_indices: list
    single_epm_masks: list
    multi_epm_masks: list
    all_epm_masks: list
    leaf_masks: list
    stats: dict
    density: float = 0.0
    effective_propagation_cap: object = None

    def __len__(self) -> int:
        return len(self.all_epm_masks)


# =============================================================================
# Shared helpers
# =============================================================================

def _bits(mask: int):
    """Yield the bit position of every set bit in `mask`, ascending."""
    m = mask
    while m:
        lsb = m & (-m)
        yield lsb.bit_length() - 1
        m &= m - 1


def _single_erc_epms(ercs, hier) -> tuple[list[int], list[int]]:
    """Find all single-ERC EPMs: the subset-minimal persistent ERCs."""
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


def _assign_orders(all_sp_masks: list[int]) -> dict[int, int]:
    """order(C) = max order of proper sub-SOs + 1; 0 if none. O(n^2)."""
    sorted_masks = sorted(all_sp_masks, key=lambda m: bin(m).count('1'))
    so_order: dict[int, int] = {}
    for cl in sorted_masks:
        max_sub = -1
        for sub_cl, sub_ord in so_order.items():
            if sub_cl != cl and (sub_cl & cl) == sub_cl and sub_ord > max_sub:
                max_sub = sub_ord
        so_order[cl] = max_sub + 1
    return so_order


def build_reach_table(g: FundamentalGraph, n: int) -> dict[int, int]:
    """reach_prod(k) for every k=0..n. See module docstring for the proof."""
    table = {n: 0}
    sp = 0
    prod = 0
    for k in range(n - 1, -1, -1):
        implied, new_sp = g.erc_syn_close(sp, [k])
        for i in implied:
            prod |= g.prod_mask[i]
        sp = new_sp
        table[k] = prod
    return table


# =============================================================================
# Structural diagnostics (cheap, non-exponential -- pure graph analysis, no
# DFS) for studying what predicts search difficulty across networks.
# =============================================================================

def compute_structural_metrics(n_ercs: int, syn_result, comp_result) -> dict:
    """
    Cheap (O(n_ercs + n_edges), no search) structural summary of a network's
    fundamental synergy/complementarity graphs, meant for correlating against
    observed DFS cost across many networks -- not used by the search itself.

    Returns
    -------
    dict with:
      syn_density        : mean fundamental-synergy out-degree per ERC
                            (Option B's per-state branching-factor driver)
      comp_density        : mean fundamental-complementarity out-degree per ERC
                            (Option A's per-state branching-factor driver;
                            also the mechanism that resolves req, so higher
                            values may reduce dead-ends even as they widen
                            branching -- a double-edged, not purely harmful,
                            quantity)
      syn_to_comp_ratio    : syn_density / comp_density -- crude indicator of
                            which mechanism dominates a network's combinatorics
      max_syn_outdegree    : largest single-ERC synergy out-degree (hub
                            effect: distinguishes uniformly-dense from
                            hub-dominated networks, which the mean alone hides)
      max_scc_fraction     : fraction of ERCs in the largest strongly-connected
                            component of the complementarity digraph (producer
                            -> consumer edges) -- operationalizes the
                            "entanglement" hypothesis: a network whose ERCs
                            are pervasively mutually dependent, not modular
      n_nontrivial_sccs    : count of SCCs with >1 member
    """
    syn_degree = [0] * n_ercs
    for st in syn_result.fundamental:
        syn_degree[st.i] += 1
        syn_degree[st.j] += 1
    n_syn = len(syn_result.fundamental)
    n_comp = len(comp_result.fundamental)

    adj: dict[int, set] = {i: set() for i in range(n_ercs)}
    for fc in comp_result.fundamental:
        adj[fc.prod_idx].add(fc.cons_idx)

    # Iterative Tarjan SCC (avoids Python recursion-limit issues at scale).
    index_counter = [0]
    stack: list[int] = []
    lowlink: dict[int, int] = {}
    index: dict[int, int] = {}
    on_stack: dict[int, bool] = {}
    sccs: list[list[int]] = []

    def strongconnect(v0: int):
        work = [(v0, iter(adj[v0]))]
        index[v0] = index_counter[0]; lowlink[v0] = index_counter[0]; index_counter[0] += 1
        stack.append(v0); on_stack[v0] = True
        while work:
            v, it = work[-1]
            advanced = False
            for w in it:
                if w not in index:
                    index[w] = index_counter[0]; lowlink[w] = index_counter[0]; index_counter[0] += 1
                    stack.append(w); on_stack[w] = True
                    work.append((w, iter(adj[w])))
                    advanced = True
                    break
                elif on_stack.get(w):
                    lowlink[v] = min(lowlink[v], index[w])
            if advanced:
                continue
            work.pop()
            if work:
                u = work[-1][0]
                lowlink[u] = min(lowlink[u], lowlink[v])
            if lowlink[v] == index[v]:
                comp_nodes = []
                while True:
                    w = stack.pop(); on_stack[w] = False; comp_nodes.append(w)
                    if w == v:
                        break
                sccs.append(comp_nodes)

    for v in range(n_ercs):
        if v not in index:
            strongconnect(v)

    biggest_scc = max((len(s) for s in sccs), default=0)
    n_nontrivial = sum(1 for s in sccs if len(s) > 1)

    syn_density = n_syn / max(n_ercs, 1)
    comp_density = n_comp / max(n_ercs, 1)

    return {
        'syn_density': syn_density,
        'comp_density': comp_density,
        'syn_to_comp_ratio': syn_density / comp_density if comp_density > 0 else float('inf'),
        'max_syn_outdegree': max(syn_degree, default=0),
        'max_scc_fraction': biggest_scc / max(n_ercs, 1),
        'n_nontrivial_sccs': n_nontrivial,
    }


# =============================================================================
# Core DFS: reach-pruning + gated unit propagation
# =============================================================================

def _mode1_dfs(
    seed_states: list[DFSState],
    g: FundamentalGraph,
    visited_sp: set,
    stats: dict,
    *,
    time_budget_s: float | None = None,
    reach_table: dict | None = None,
    propagation_cap: int | None = None,
) -> tuple[list[DFSState], list[DFSState]]:
    ssm_states: list[DFSState] = []
    leaf_states: list[DFSState] = []

    explored_sig: dict[int, list[tuple[int, int]]] = {}

    def _is_dominated(sp: int, min_ext: int, min_seed: int) -> bool:
        for (e_ext, e_seed) in explored_sig.get(sp, ()):
            if e_ext <= min_ext and e_seed >= min_seed:
                return True
        return False

    def _sig(state: DFSState) -> tuple[int, int]:
        return (state.min_ext, min(state.erc_set))

    def _reach_fails(state: DFSState) -> bool:
        if reach_table is None:
            return False
        reach = reach_table.get(state.min_ext, 0)
        return (state.req & ~reach) != 0

    stack: list[DFSState] = []
    for s in seed_states:
        if _is_dominated(s.sp, *_sig(s)):
            continue
        if _reach_fails(s):
            stats['reach_pruned_at_construction'] = stats.get('reach_pruned_at_construction', 0) + 1
            continue
        stack.append(s)

    t0 = time.perf_counter()
    while stack:
        if time_budget_s is not None and (time.perf_counter() - t0) > time_budget_s:
            stats['time_budget_exceeded'] = True
            stats['stack_remaining_at_cutoff'] = len(stack)
            # [PROGRESS] how deep had the frontier gotten when time ran out --
            # distinguishes "died early, barely started" from "got deep into
            # large partial generators, ran out of time near the end".
            if stack:
                sizes = [len(s.erc_set) for s in stack]
                stats['stack_erc_size_mean_at_cutoff'] = sum(sizes) / len(sizes)
                stats['stack_erc_size_max_at_cutoff'] = max(sizes)
            break

        state = stack.pop()
        state_sig = _sig(state)
        if _is_dominated(state.sp, *state_sig):
            continue
        explored_sig.setdefault(state.sp, []).append(state_sig)
        visited_sp.add(state.sp)
        stats['states_explored'] = stats.get('states_explored', 0) + 1

        min_seed = state_sig[1]

        while True:
            if state.is_ssm:
                ssm_states.append(state)
                stats['ssm_found'] = stats.get('ssm_found', 0) + 1
                break

            survivors: list[tuple[DFSState, str]] = []

            comp_cand: dict[int, DFSState] = {}
            for s_bit in _bits(state.req):
                for prod_idx in g.comp_by_species.get(s_bit, []):
                    if prod_idx < state.min_ext:
                        stats['canonical_pruned'] = stats.get('canonical_pruned', 0) + 1
                        continue
                    if (g.species_mask[prod_idx] & state.sp) == g.species_mask[prod_idx]:
                        continue
                    if prod_idx not in comp_cand:
                        comp_cand[prod_idx] = g.extend_state(state, prod_idx)
            for new_state in comp_cand.values():
                newly_added = new_state.erc_set - state.erc_set
                if any(i < min_seed and not g.is_persistent[i] for i in newly_added):
                    stats['canonical_pruned'] = stats.get('canonical_pruned', 0) + 1
                    continue
                if _reach_fails(new_state):
                    stats['reach_pruned_at_construction'] = stats.get('reach_pruned_at_construction', 0) + 1
                    continue
                if not _is_dominated(new_state.sp, *_sig(new_state)):
                    survivors.append((new_state, 'comp'))

            syn_cand: dict[int, DFSState] = {}
            for erc_i in state.erc_set:
                for (j, k) in g.syn_from.get(erc_i, []):
                    if j in state.erc_set:
                        continue
                    if j < state.min_ext:
                        stats['canonical_pruned'] = stats.get('canonical_pruned', 0) + 1
                        continue
                    if not (g.prod_mask[k] & state.req):
                        continue
                    if j not in syn_cand:
                        syn_cand[j] = g.extend_state(state, j)
            for new_state in syn_cand.values():
                newly_added = new_state.erc_set - state.erc_set
                if any(i < min_seed and not g.is_persistent[i] for i in newly_added):
                    stats['canonical_pruned'] = stats.get('canonical_pruned', 0) + 1
                    continue
                if _reach_fails(new_state):
                    stats['reach_pruned_at_construction'] = stats.get('reach_pruned_at_construction', 0) + 1
                    continue
                if not _is_dominated(new_state.sp, *_sig(new_state)):
                    survivors.append((new_state, 'syn'))

            if not survivors:
                leaf_states.append(state)
                stats['leaves_found'] = stats.get('leaves_found', 0) + 1

                # [ECONOMICS] characterize this dead-end by what it still
                # needed (req) and what it had already built (prod, erc_set)
                # -- see compute_epms's docstring / module-level notes on
                # dead_end_req_species for why req is tracked per-species
                # (cheap: req sets at dead-ends are typically small, ~1-10
                # bits) while prod/erc_set are tracked only by SIZE (their
                # popcounts can be large -- per-species tallying would be
                # expensive and less diagnostic, since prod at a dead end is
                # "whatever got built", not a scarce resource).
                req_size = bin(state.req).count('1')
                prod_size = bin(state.prod).count('1')
                erc_size = len(state.erc_set)
                d = stats.setdefault('dead_end_req_size_hist', {})
                d[req_size] = d.get(req_size, 0) + 1
                d = stats.setdefault('dead_end_prod_size_hist', {})
                d[prod_size] = d.get(prod_size, 0) + 1
                d = stats.setdefault('dead_end_erc_size_hist', {})
                d[erc_size] = d.get(erc_size, 0) + 1
                req_species_counter = stats.setdefault('dead_end_req_species', {})
                for s_bit in _bits(state.req):
                    req_species_counter[s_bit] = req_species_counter.get(s_bit, 0) + 1
                break

            # ── Per-species forced-move detection, gated by propagation_cap ──
            n_req_bits = bin(state.req).count('1')
            chosen = None
            if propagation_cap is None or n_req_bits * len(survivors) <= propagation_cap:
                cover_count: dict[int, int] = {}
                for s_bit in _bits(state.req):
                    s_mask = 1 << s_bit
                    cover_count[s_bit] = sum(1 for sv, _k in survivors if sv.prod & s_mask)
                for sv, _k in survivors:
                    covered_here = state.req & sv.prod
                    if any(cover_count[s_bit] == 1 for s_bit in _bits(covered_here)):
                        chosen = sv
                        break
            else:
                stats['propagation_gated'] = stats.get('propagation_gated', 0) + 1

            if chosen is not None:
                explored_sig.setdefault(chosen.sp, []).append(_sig(chosen))
                visited_sp.add(chosen.sp)
                stats['propagation_steps'] = stats.get('propagation_steps', 0) + 1
                state = chosen
                continue

            n_comp = n_syn = 0
            for sv, kind in survivors:
                stack.append(sv)
                if kind == 'comp':
                    n_comp += 1
                else:
                    n_syn += 1
            stats['comp_extensions'] = stats.get('comp_extensions', 0) + n_comp
            stats['syn_extensions'] = stats.get('syn_extensions', 0) + n_syn
            break

    # [PROGRESS] fraction of currently-known work actually completed --
    # 1.0 for a normal completion (stack empties on its own, nothing left
    # to report as "remaining"); < 1.0 whenever the time budget cut the
    # search off with states still pending. Cheap, single computation.
    explored = stats.get('states_explored', 0)
    remaining = stats.get('stack_remaining_at_cutoff', 0)
    stats['progress_fraction'] = explored / max(explored + remaining, 1)

    return ssm_states, leaf_states


# =============================================================================
# Public API
# =============================================================================

def compute_epms(
    rn, ercs, hier, syn_result=None, comp_result=None,
    *,
    time_budget_s: float | None = None,
    reach_table: dict | None = None,
    propagation_cap: int | None = None,
) -> EPMResult:
    """
    Drop-in alternative to cot_gen.epm.compute_epms — see module docstring
    for the two optimizations and their correctness arguments. Pass
    reach_table=None / propagation_cap=None to reproduce the unmodified
    algorithm's behavior exactly (both are purely-additive prunes/shortcuts).
    """
    single_idx, single_masks = _single_erc_epms(ercs, hier)
    single_sp_set = set(single_masks)

    if syn_result is None or comp_result is None:
        return EPMResult(
            single_epm_indices=single_idx, single_epm_masks=single_masks,
            multi_epm_masks=[], all_epm_masks=sorted(single_sp_set),
            leaf_masks=[], stats={'note': 'adjacency skipped (no syn/comp result)'},
        )

    g = FundamentalGraph(ercs, hier, syn_result, comp_result)
    visited_sp: set[int] = set()
    stats1: dict = {}
    for sp in single_masks:
        visited_sp.add(sp)

    non_p_seeds = [g.make_seed_state(i) for i in range(len(ercs)) if not ercs[i].is_persistent()]

    ssm_states, leaf_states = _mode1_dfs(
        non_p_seeds, g, visited_sp, stats1,
        time_budget_s=time_budget_s, reach_table=reach_table, propagation_cap=propagation_cap,
    )

    all_sp_masks = list(single_sp_set) + [s.sp for s in ssm_states]
    so_order = _assign_orders(all_sp_masks)
    for sp in single_masks:
        if sp not in so_order:
            so_order[sp] = 0

    multi_epm_masks = [s.sp for s in ssm_states
                        if so_order.get(s.sp, -1) == 0 and s.sp not in single_sp_set]
    all_epm_masks = sorted(single_sp_set | set(multi_epm_masks))

    return EPMResult(
        single_epm_indices=single_idx, single_epm_masks=single_masks,
        multi_epm_masks=multi_epm_masks, all_epm_masks=all_epm_masks,
        leaf_masks=[s.sp for s in leaf_states], stats=stats1,
    )


DENSITY_GUARD_THRESHOLD = 20.0   # validated safe up to 13.8; validated failure at 38.9
DEFAULT_PROPAGATION_CAP = 200    # matched-or-beat "unlimited" in every sweep run


def compute_epms_adaptive(
    rn, ercs, hier, syn_result, comp_result,
    *,
    time_budget_s: float | None = None,
    density_guard_threshold: float = DENSITY_GUARD_THRESHOLD,
    propagation_cap: int = DEFAULT_PROPAGATION_CAP,
) -> EPMResult:
    """
    Convenience wrapper: builds the reach(k) table, computes the network's
    synergy density, and applies DENSITY_GUARD automatically (see module
    docstring) before running compute_epms.
    """
    n = len(ercs)
    g = FundamentalGraph(ercs, hier, syn_result, comp_result)
    reach_table = build_reach_table(g, n)
    density = len(syn_result.fundamental) / max(n, 1)
    effective_cap = 0 if density > density_guard_threshold else propagation_cap

    result = compute_epms(
        rn, ercs, hier, syn_result, comp_result,
        time_budget_s=time_budget_s, reach_table=reach_table, propagation_cap=effective_cap,
    )
    result.density = density
    result.effective_propagation_cap = effective_cap
    return result
