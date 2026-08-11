"""
epm_risk_analysis.py — instruments Mode-1 DFS with a refined safe/mid/risky
classification for every synergy and complementarity extension move, plus
exploration-volume accounting (not just outcome counting) so contribution
share can be weighed against effort share.

Classification (as specified, Tomas Veloz)
--------------------------------------------
SYNERGY move (existing member `i`, new candidate `j`, triggering (i,j)->k):
  producible(x) := species_mask[x] (req_mask[x] | prod_mask[x]) subset-or-
                   equal state.prod  -- x's own need AND output are already
                   backed by production established so far, not merely
                   "present" via req.
  safe synergy      : producible(i) AND producible(j)   -- always true when
                       i, j both belong to an already self-maintaining SO
                       (proven this session: req=0 implies every member's
                       species_mask subset-or-equal that SO's prod).
  mid-safe synergy   : exactly one of {i, j} producible.
  risky synergy      : neither producible -- pure speculation on both sides.

COMPLEMENTARITY move (receiving state, newly-implied requirement add_req =
union of req_mask over everything erc_syn_close pulls in when candidate c
is added -- the "giver"):
  (SUPERSEDES an earlier version of this criterion -- total_req < parts_req,
  i.e. "any reduction vs. the naive disjoint sum" -- which Tomas Veloz
  identified as vacuous: the DFS's own selection rule guarantees the
  triggering species always cancels, so that test was structurally
  guaranteed to always read "safe" and carried no signal. Counterexample
  that exposed it: state {a,b}->{c,d} (req size 2) pulls in a producer of
  `a`, itself needing {x,y,z} (req size 3); combined req becomes {b,x,y,z}
  (size 4) -- strictly WORSE than the 2 we started with, even though it's
  still less than the naive disjoint sum of 2+3=5. "Better than pretending
  there's no overlap" is not the same question as "was this actually a
  good trade".)
  mutual        := add_req & state.prod != 0 -- the giver ALSO needs
                   something the receiver already produces (a genuine
                   two-way exchange, not one-directional charity).
  giver_size = popcount(add_req); receiver_size = popcount(state.req).
  safe complementarity : giver_size == 0 (a pure gift -- the giver has no
                          outstanding need of its own once combined, so it
                          can never be a bad trade), OR mutual, OR
                          giver_size < receiver_size (the giver is "less
                          needy" than what we were already carrying).
  risky complementarity : none of the above -- a one-directional gift from
                          a needier producer nets a LARGER outstanding
                          requirement than the one being resolved.

A branch is tagged with the WORST classification any move along its
ancestry has used: risky > mid-safe > safe (synergy and complementarity
share one branch-level scale: any risky synergy OR risky complementarity
makes a branch "risky"; mid-safe synergy with no risky move makes it
"mid"; otherwise "safe"). ALONGSIDE the category, each branch also carries
the MOVE TYPE ('syn' or 'comp') that is currently responsible for that
category -- updated whenever a move's own category is at least as bad as
everything inherited so far (so for an all-safe branch, type tracks the
most recent move; for a branch that escalated to mid/risky, type tracks
the most recent move AT that worst level). This is what lets every
downstream stat (exploration volume, contribution, success rate) be split
by "was this driven by a synergy step or a complementarity step", not just
by category -- without that, the outcome-level plots couldn't be traced
back to the move-classification plot, which already had this split.
Tagging is a side dict keyed by id(state), never touching DFSState's
fields/hashing -- purely additive, does not alter control flow, pruning,
or dominance logic (an unmodified copy of epm_reach_propagate.py's
validated _mode1_dfs).

Both OUTCOMES (SSM / dead end) and EXPLORATION VOLUME (every state actually
popped and processed, whether or not it's a terminal state) are tallied by
(category, move type) -- the volume side is what answers "what proportion
of the total search effort produced that contribution", which a pure
outcome/success-rate view can't answer on its own.

Validated by exact-match against epm_reach_propagate.compute_epms (same
SSM/leaf sp-mask sets, same states_explored) before trusting any risk
statistic -- see epm_risk_analysis_validate.py.
"""
from __future__ import annotations

import time
from dataclasses import dataclass

from pyCOT.analysis.organizations.fundamental_graph import FundamentalGraph, DFSState

RISK_ORDER = {'safe': 0, 'mid': 1, 'risky': 2}


def _bits(mask: int):
    m = mask
    while m:
        lsb = m & (-m)
        yield lsb.bit_length() - 1
        m &= m - 1


def _popcount(mask: int) -> int:
    return bin(mask).count('1')


def _assign_orders(all_sp_masks: list[int]) -> dict[int, int]:
    sorted_masks = sorted(all_sp_masks, key=lambda m: bin(m).count('1'))
    so_order: dict[int, int] = {}
    for cl in sorted_masks:
        max_sub = -1
        for sub_cl, sub_ord in so_order.items():
            if sub_cl != cl and (sub_cl & cl) == sub_cl and sub_ord > max_sub:
                max_sub = sub_ord
        so_order[cl] = max_sub + 1
    return so_order


def _single_erc_epms(ercs, hier):
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


def build_reach_table(g: FundamentalGraph, n: int) -> dict[int, int]:
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


def _producible(g: FundamentalGraph, x: int, prod: int) -> bool:
    return (g.species_mask[x] & ~prod) == 0


def _classify_synergy(g: FundamentalGraph, i: int, j: int, state_prod: int) -> str:
    pi = _producible(g, i, state_prod)
    pj = _producible(g, j, state_prod)
    if pi and pj:
        return 'safe'
    if pi or pj:
        return 'mid'
    return 'risky'


def _implied_add_req(g: FundamentalGraph, state: DFSState, new_idx: int) -> int:
    implied, _ = g.erc_syn_close(state.sp, [new_idx])
    add_req = 0
    for i in implied:
        add_req |= g.req_mask[i]
    return add_req


def _classify_complementarity(state: DFSState, add_req: int, trigger_bit: int) -> str:
    """
    trigger_bit: the specific species whose presence in state.req/state.prod
    caused this candidate to be proposed in the first place. It is EXCLUDED
    from the mutual-overlap check below: by construction of comp_by_species
    (Mode-1: producers of a species the receiver needs) and especially
    comp_consumers_by_species (Mode-2: consumers of a species the receiver
    already produces), that one species is guaranteed to sit in both add_req
    and state.prod regardless of anything else -- checking it "overlaps"
    would be tautological, not evidence of a genuine two-way exchange (this
    is exactly the same vacuity bug the original total_req<parts_req
    criterion had, just one level deeper -- caught empirically: Mode-2
    complementarity stayed 100% "safe" even after the first fix, whereas
    Mode-1 (where the giver's own req doesn't structurally have to contain
    the triggering species) showed real variance, exposing the asymmetry).
    """
    giver_size = _popcount(add_req)
    if giver_size == 0:
        return 'safe'                       # pure gift, no outstanding need of its own
    residual_overlap = add_req & state.prod & ~(1 << trigger_bit)
    if residual_overlap != 0:
        return 'safe'                       # mutual: giver ALSO needs something else we already have
    receiver_size = _popcount(state.req)
    return 'safe' if giver_size < receiver_size else 'risky'


@dataclass
class RiskEPMResult:
    all_epm_masks: list
    leaf_masks: list
    stats: dict


def _mode1_dfs_risk(
    seed_states: list[DFSState],
    g: FundamentalGraph,
    visited_sp: set,
    stats: dict,
    *,
    time_budget_s: float | None = None,
    reach_table: dict | None = None,
    propagation_cap: int | None = None,
    seed_info: dict[int, tuple[str, str | None]] | None = None,
) -> tuple[list[DFSState], list[DFSState]]:
    """
    seed_info: optional {id(seed_state): (category, move_type)}, used by the
    ESPM/Mode-2 caller to seed each starting state with the (category,
    move_type) of the Mode-2 move that produced it (instead of the default
    ('safe', None) used for Mode-1's own from-scratch, single-ERC seeds).
    """
    ssm_states: list[DFSState] = []
    leaf_states: list[DFSState] = []

    explored_sig: dict[int, list[tuple[int, int]]] = {}
    # id(state) -> (worst category on its ancestry, move type responsible for it)
    branch_info: dict[int, tuple[str, str | None]] = dict(seed_info) if seed_info else {}

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

    def _tag(parent: DFSState, child: DFSState, move_cat: str, move_type: str) -> None:
        inherited_cat, inherited_type = branch_info.get(id(parent), ('safe', None))
        if RISK_ORDER[move_cat] >= RISK_ORDER[inherited_cat]:
            branch_info[id(child)] = (move_cat, move_type)
        else:
            branch_info[id(child)] = (inherited_cat, inherited_type)
        key = f'{move_type}_moves_{move_cat}'
        stats[key] = stats.get(key, 0) + 1

    def _tally_states_explored(state: DFSState) -> None:
        cat, mtype = branch_info.get(id(state), ('safe', None))
        stats[f'states_explored_{cat}'] = stats.get(f'states_explored_{cat}', 0) + 1
        if mtype is not None:   # None only for never-extended seeds -- see module docstring
            key = f'states_explored_{mtype}_{cat}'
            stats[key] = stats.get(key, 0) + 1

    def _tally_outcome(state: DFSState, outcome: str) -> None:
        cat, mtype = branch_info.get(id(state), ('safe', None))
        stats[f'{outcome}_via_{cat}'] = stats.get(f'{outcome}_via_{cat}', 0) + 1
        if mtype is not None:
            key = f'{outcome}_via_{mtype}_{cat}'
            stats[key] = stats.get(key, 0) + 1

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
            break

        state = stack.pop()
        state_sig = _sig(state)
        if _is_dominated(state.sp, *state_sig):
            continue
        explored_sig.setdefault(state.sp, []).append(state_sig)
        visited_sp.add(state.sp)
        stats['states_explored'] = stats.get('states_explored', 0) + 1
        _tally_states_explored(state)

        min_seed = state_sig[1]

        while True:
            if state.is_ssm:
                ssm_states.append(state)
                stats['ssm_found'] = stats.get('ssm_found', 0) + 1
                _tally_outcome(state, 'ssm')
                break

            survivors: list[tuple[DFSState, str]] = []

            comp_cand: dict[int, tuple[DFSState, str]] = {}
            for s_bit in _bits(state.req):
                for prod_idx in g.comp_by_species.get(s_bit, []):
                    if prod_idx < state.min_ext:
                        stats['canonical_pruned'] = stats.get('canonical_pruned', 0) + 1
                        continue
                    if (g.species_mask[prod_idx] & state.sp) == g.species_mask[prod_idx]:
                        continue
                    if prod_idx not in comp_cand:
                        add_req = _implied_add_req(g, state, prod_idx)
                        new_state = g.extend_state(state, prod_idx)
                        cat = _classify_complementarity(state, add_req, s_bit)
                        comp_cand[prod_idx] = (new_state, cat)
            for new_state, cat in comp_cand.values():
                newly_added = new_state.erc_set - state.erc_set
                if any(i < min_seed and not g.is_persistent[i] for i in newly_added):
                    stats['canonical_pruned'] = stats.get('canonical_pruned', 0) + 1
                    continue
                if _reach_fails(new_state):
                    stats['reach_pruned_at_construction'] = stats.get('reach_pruned_at_construction', 0) + 1
                    continue
                if not _is_dominated(new_state.sp, *_sig(new_state)):
                    _tag(state, new_state, cat, 'comp')
                    survivors.append((new_state, 'comp'))

            syn_cand: dict[int, tuple[DFSState, str]] = {}
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
                        cat = _classify_synergy(g, erc_i, j, state.prod)
                        new_state = g.extend_state(state, j)
                        syn_cand[j] = (new_state, cat)
            for new_state, cat in syn_cand.values():
                newly_added = new_state.erc_set - state.erc_set
                if any(i < min_seed and not g.is_persistent[i] for i in newly_added):
                    stats['canonical_pruned'] = stats.get('canonical_pruned', 0) + 1
                    continue
                if _reach_fails(new_state):
                    stats['reach_pruned_at_construction'] = stats.get('reach_pruned_at_construction', 0) + 1
                    continue
                if not _is_dominated(new_state.sp, *_sig(new_state)):
                    _tag(state, new_state, cat, 'syn')
                    survivors.append((new_state, 'syn'))

            if not survivors:
                leaf_states.append(state)
                stats['leaves_found'] = stats.get('leaves_found', 0) + 1
                _tally_outcome(state, 'dead_end')
                break

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

            for sv, kind in survivors:
                stack.append(sv)
            break

    explored = stats.get('states_explored', 0)
    remaining = stats.get('stack_remaining_at_cutoff', 0)
    stats['progress_fraction'] = explored / max(explored + remaining, 1)

    return ssm_states, leaf_states


def compute_epms_risk(
    rn, ercs, hier, syn_result, comp_result,
    *,
    time_budget_s: float | None = None,
    reach_table: dict | None = None,
    propagation_cap: int | None = None,
) -> RiskEPMResult:
    single_idx, single_masks = _single_erc_epms(ercs, hier)
    single_sp_set = set(single_masks)

    g = FundamentalGraph(ercs, hier, syn_result, comp_result)
    visited_sp: set[int] = set()
    stats: dict = {}
    for sp in single_masks:
        visited_sp.add(sp)

    non_p_seeds = [g.make_seed_state(i) for i in range(len(ercs)) if not ercs[i].is_persistent()]

    ssm_states, leaf_states = _mode1_dfs_risk(
        non_p_seeds, g, visited_sp, stats,
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

    return RiskEPMResult(
        all_epm_masks=all_epm_masks,
        leaf_masks=[s.sp for s in leaf_states],
        stats=stats,
    )
