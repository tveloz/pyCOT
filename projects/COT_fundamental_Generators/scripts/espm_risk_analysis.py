"""
espm_risk_analysis.py — extends epm_risk_analysis.py's safe/mid/risky move
classification into Mode-2 (ESPM/BFS-by-order growth from already-known SOs).

Structural prediction to check (not assumed): Mode-2 only ever grows
outward from an ALREADY self-maintaining SO (so_state.req == 0, by
construction -- current_layer is seeded from confirmed EPMs and only ever
refilled with confirmed new SSM states). That means, for every Mode-2
SYNERGY move (i, j) where i is an existing member of so_state.erc_set: i is
*always* producible (species_mask[i] subset-or-equal so_state.prod, proven
this session for any member of a self-maintaining set). So Mode-2 synergy
should be structurally incapable of "risky" (which needs BOTH sides
unproducible) -- only "safe" (j also producible) or "mid" (j not) should
ever occur. Verified empirically here, not presumed.

COMPLEMENTARITY in Mode-2 is a DIFFERENT prediction than Mode-1, under the
corrected criterion (see epm_risk_analysis.py's module docstring for why
the original "any reduction vs. naive sum" test was vacuous and got
replaced with mutual-exchange / giver-vs-receiver-size). so_state.req is
always 0 here (already self-maintaining), so the receiver_size term is
always 0 -- "giver_size < receiver_size" can never fire. That means every
Mode-2 complementarity move is "safe" only if the consumer candidate is a
pure gift (giver_size == 0, i.e. its own requirement is already fully
covered by so_state.prod) or mutual (it also needs something so_state
already produces) -- otherwise "risky", since ANY nonzero, non-mutual
requirement is strictly worse than the SO's current zero debt. Unlike
Mode-1, this is NOT expected to come out overwhelmingly safe -- verify,
don't assume.

HIERARCHY-LIFT moves (via g.parents -- vertical lift to a covering ERC) are
a third Mode-2 move type the user's synergy/complementarity taxonomy does
not cover. Rather than invent a criterion, this module tracks them as a
separate, UNCLASSIFIED count (mode2_hierarchy_moves) and treats them as
category-neutral (they don't change a branch's inherited risk category) --
flagged explicitly rather than silently folded into either bucket.

Non-checkpointed (unlike espm_resumable.py): this is an analysis/statistics
tool, not the production resumable pipeline -- a time_budget_s safety net
is provided instead. Validated by exact-match (SO species-mask sets, by
order) against espm_resumable.compute_espm_resumable before trusting any
risk statistic -- see espm_risk_analysis_validate.py.
"""
from __future__ import annotations

import time
from dataclasses import dataclass

from cot_gen.fundamental_graph import FundamentalGraph, DFSState
from cot_gen.epm import EPMResult

from epm_risk_analysis import (
    _bits, _popcount, _producible, _classify_synergy, _classify_complementarity,
    _implied_add_req, _mode1_dfs_risk, RISK_ORDER,
)


@dataclass
class RiskESPMResult:
    espm_by_order: dict            # order -> list of species masks
    stats: dict


def compute_espm_risk(
    rn, ercs, hier, syn_result, comp_result, epm_result: EPMResult,
    *,
    max_order: int = 10,
    time_budget_s: float | None = None,
) -> RiskESPMResult:
    t0 = time.perf_counter()
    stats: dict = {}

    if epm_result._graph is not None:
        g           = epm_result._graph
        visited_sp  = set(epm_result._visited_sp)
        sp_to_state = dict(epm_result._sp_to_state)
        so_order    = dict(epm_result._so_order)
    else:
        g = FundamentalGraph(ercs, hier, syn_result, comp_result)
        visited_sp  = set()
        sp_to_state = {}
        so_order    = {sp: 0 for sp in epm_result.all_epm_masks}
        for sp in epm_result.all_epm_masks:
            for i, e in enumerate(ercs):
                if e.species_mask == sp and e.is_persistent():
                    sp_to_state[sp] = g.make_seed_state(i)
                    break

    current_layer = [sp_to_state[sp] for sp in epm_result.all_epm_masks if sp in sp_to_state]

    order = 1
    while order <= max_order:
        if time_budget_s is not None and (time.perf_counter() - t0) > time_budget_s:
            stats['time_budget_exceeded'] = True
            break
        if not current_layer:
            break

        # ── candidate collection, classified at the point of creation ──────
        # id(ext_state) -> (worst category so far, move type ('syn'/'comp') responsible for it)
        cand_info: dict[int, tuple[str, str | None]] = {}
        best_by_sp: dict[int, list[DFSState]] = {}
        for so_state in current_layer:
            local_cand: dict[int, DFSState] = {}
            local_cat: dict[int, str | None] = {}
            local_type: dict[int, str | None] = {}

            for i in so_state.erc_set:
                for (j, _) in g.syn_from.get(i, []):
                    if j in so_state.erc_set:
                        continue
                    if j not in local_cand:
                        cat = _classify_synergy(g, i, j, so_state.prod)
                        local_cand[j] = g.extend_state(so_state, j)
                        local_cat[j] = cat
                        local_type[j] = 'syn'
                        stats[f'mode2_syn_moves_{cat}'] = stats.get(f'mode2_syn_moves_{cat}', 0) + 1
                for a in g.parents[i]:
                    if a in so_state.erc_set:
                        continue
                    if (g.species_mask[a] & so_state.sp) == g.species_mask[a]:
                        continue
                    if a not in local_cand:
                        local_cand[a] = g.extend_state(so_state, a)
                        local_cat[a] = None   # unclassified -- see module docstring
                        local_type[a] = None
                        stats['mode2_hierarchy_moves'] = stats.get('mode2_hierarchy_moves', 0) + 1

            for s_bit in _bits(so_state.prod):
                for cons_idx in g.comp_consumers_by_species.get(s_bit, []):
                    if cons_idx in so_state.erc_set:
                        continue
                    if (g.species_mask[cons_idx] & so_state.sp) == g.species_mask[cons_idx]:
                        continue
                    if cons_idx not in local_cand:
                        add_req = _implied_add_req(g, so_state, cons_idx)
                        ext_state = g.extend_state(so_state, cons_idx)
                        cat = _classify_complementarity(so_state, add_req, s_bit)
                        local_cand[cons_idx] = ext_state
                        local_cat[cons_idx] = cat
                        local_type[cons_idx] = 'comp'
                        stats[f'mode2_comp_moves_{cat}'] = stats.get(f'mode2_comp_moves_{cat}', 0) + 1

            for idx, ext_state in local_cand.items():
                cat, mtype = local_cat[idx], local_type[idx]
                # inherit so_state's own (category, move type) -- ('safe', None)
                # default for order-0 EPM seeds, which carry no risk of their
                # own, only whatever this NEW move contributes
                parent_cat, parent_type = cand_info.get(id(so_state), ('safe', None))
                if cat is None:
                    worst_cat, worst_type = parent_cat, parent_type   # hierarchy: neutral, don't escalate
                elif RISK_ORDER[cat] >= RISK_ORDER[parent_cat]:
                    worst_cat, worst_type = cat, mtype
                else:
                    worst_cat, worst_type = parent_cat, parent_type

                if ext_state.sp in visited_sp:
                    continue
                frontier = best_by_sp.setdefault(ext_state.sp, [])
                new_seed = min(ext_state.erc_set)
                dominated = False
                survivors = []
                for c in frontier:
                    cand_seed = min(c.erc_set)
                    if c.min_ext <= ext_state.min_ext and cand_seed >= new_seed:
                        dominated = True
                        survivors.append(c)
                    elif ext_state.min_ext <= c.min_ext and new_seed >= cand_seed:
                        continue
                    else:
                        survivors.append(c)
                if not dominated:
                    survivors.append(ext_state)
                    cand_info[id(ext_state)] = (worst_cat, worst_type)
                best_by_sp[ext_state.sp] = survivors

        singleton_seeds: list[DFSState] = []
        tied_groups: list[list[DFSState]] = []
        for frontier in best_by_sp.values():
            if len(frontier) == 1:
                singleton_seeds.append(frontier[0])
            else:
                tied_groups.append(frontier)

        if not singleton_seeds and not tied_groups:
            break

        new_ssm_states: list[DFSState] = []
        stats_k: dict = {}
        if singleton_seeds:
            seed_info = {id(s): cand_info.get(id(s), ('safe', None)) for s in singleton_seeds}
            ssms, _leaves = _mode1_dfs_risk(singleton_seeds, g, visited_sp, stats_k,
                                             seed_info=seed_info)
            new_ssm_states.extend(ssms)
        for tied in tied_groups:
            for rep in tied:
                if rep.sp in visited_sp:
                    visited_sp.discard(rep.sp)
                seed_info = {id(rep): cand_info.get(id(rep), ('safe', None))}
                ssms, _leaves = _mode1_dfs_risk([rep], g, visited_sp, stats_k,
                                                 seed_info=seed_info)
                new_ssm_states.extend(ssms)

        for key, val in stats_k.items():
            if isinstance(val, (int, float)):
                stats[key] = stats.get(key, 0) + val

        for s in sorted(new_ssm_states, key=lambda st: bin(st.sp).count('1')):
            if s.sp not in so_order:
                max_sub = -1
                for sub_sp, sub_ord in so_order.items():
                    if sub_sp != s.sp and (sub_sp & s.sp) == sub_sp and sub_ord > max_sub:
                        max_sub = sub_ord
                so_order[s.sp] = max_sub + 1
            sp_to_state[s.sp] = s

        current_layer = new_ssm_states
        order += 1
        if not new_ssm_states:
            break

    espm_by_order: dict[int, list[int]] = {}
    for sp, ord_k in so_order.items():
        if ord_k >= 1:
            espm_by_order.setdefault(ord_k, []).append(sp)
    for k in espm_by_order:
        espm_by_order[k].sort()

    return RiskESPMResult(espm_by_order=espm_by_order, stats=stats)
