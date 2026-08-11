"""
deep_report.py — Deep structural analysis + presentation-quality reporting
for a single reaction network's full COT pipeline.

Built on top of the validated computation (cot_gen.epm.compute_epms /
compute_espm, cot_gen.max_semiorg), never reimplementing the search itself
-- this module only INSTRUMENTS the existing, oracle-validated traversal to
record extra bookkeeping (candidate-convergence stats, per-move-type
provenance) for reporting, and adds visualization on top.

Public API
----------
compute_hierarchy_stats(ercs, hier, syn, comp) -> HierarchyStats
compute_epms_instrumented(ercs, hier, syn, comp) -> (EPMResult, DegeneracyStats)
compute_espm_instrumented(ercs, hier, syn, comp, epm_result, max_order=10)
    -> (ESPMResult, dict[order -> MoveTypeCounts], DegeneracyStats)
build_so_lattice(ercs, epm_result, espm_result, so_order) -> SOLattice

plot_hierarchy_overview(...) -> html path
plot_epm_hierarchy(...) -> html path
plot_so_lattice(...) -> html path
plot_degeneracy(...) -> png path
plot_espm_composition(...) -> png path
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field

from pyCOT.analysis.organizations.epm import _bits, _mode1_dfs, _single_erc_epms, EPMResult, ESPMResult
from pyCOT.analysis.organizations.fundamental_graph import FundamentalGraph

# Shared palette -- consistent with cot_gen/explorer.py's relation colors.
COLOR_CONTAINMENT = "#7f8c8d"
COLOR_SYNERGY = "#e67e22"
COLOR_COMPLEMENTARITY = "#2980b9"
COLOR_VERTICAL_LIFT = "#8e44ad"
COLOR_EPM = "#27ae60"
COLOR_ESPM = "#f1c40f"
COLOR_NEUTRAL = "#d0d3d4"
COLOR_MAXSO = "#c0392b"

# Fixed categorical order for the four SO-discovery sources (never reordered).
# "carryover" = SOs of order >=1 that Mode-1's OWN DFS already reached
# directly (a larger SSM containing a smaller one as a proper subset,
# found by requirement-driven closure alone, with no Mode-2 growth step
# involved at all) -- these are real order-k SOs and must be counted, or
# per-order totals silently undercount relative to espm_by_order.
MOVE_TYPE_ORDER = ["synergy", "complementarity", "vertical_lift", "carryover"]
MOVE_TYPE_COLOR = {
    "synergy": COLOR_SYNERGY,
    "complementarity": COLOR_COMPLEMENTARITY,
    "vertical_lift": COLOR_VERTICAL_LIFT,
    "carryover": "#95a5a6",
}


# ---------------------------------------------------------------------------
# Hierarchy statistics
# ---------------------------------------------------------------------------

@dataclass
class HierarchyStats:
    n_ercs: int
    n_persistent: int
    sizes: list[int]                 # species count per ERC
    levels: list[int]                # containment level per ERC (0 = minimal)
    n_levels: int
    syn_degree: list[int]            # fundamental synergies touching each ERC
    comp_out_degree: list[int]       # as fundamental producer
    comp_in_degree: list[int]        # as fundamental consumer
    n_fundamental_syn: int
    n_fundamental_comp: int
    n_distinct_syn_pairs: int        # deduplicated (i,j) regardless of target count
    hasse_edges: int
    n_comparable_pairs: int
    n_incomparable_pairs: int


def _compute_containment_levels(hier, n: int) -> list[int]:
    """Level 0 = subset-minimal ERCs; level(i) = 1 + max(level(child)) for
    children in hier.children, processed in increasing descendant-count
    order (matches cot_gen/explorer.py's identical routine)."""
    order = sorted(range(n), key=lambda i: len(hier.descendants[i]))
    level = [0] * n
    children_of = [[] for _ in range(n)]
    for i in range(n):
        for p in hier.parents[i]:
            children_of[p].append(i)
    for i in order:
        if children_of[i]:
            level[i] = 1 + max(level[c] for c in children_of[i])
    return level


def compute_hierarchy_stats(ercs, hier, syn, comp) -> HierarchyStats:
    n = len(ercs)
    sizes = [e.size() for e in ercs]
    levels = _compute_containment_levels(hier, n)
    syn_degree = [0] * n
    syn_pairs = set()
    for st in syn.fundamental:
        syn_degree[st.i] += 1
        syn_degree[st.j] += 1
        syn_pairs.add((min(st.i, st.j), max(st.i, st.j)))
    comp_out = [0] * n
    comp_in = [0] * n
    for fc in comp.fundamental:
        comp_out[fc.prod_idx] += 1
        comp_in[fc.cons_idx] += 1
    n_pairs = n * (n - 1) // 2
    n_comp_pairs = sum(len(hier.ancestors[i]) for i in range(n))
    hasse_edges = sum(len(hier.parents[i]) for i in range(n))
    return HierarchyStats(
        n_ercs=n,
        n_persistent=sum(1 for e in ercs if e.is_persistent()),
        sizes=sizes, levels=levels, n_levels=(max(levels) + 1 if levels else 0),
        syn_degree=syn_degree, comp_out_degree=comp_out, comp_in_degree=comp_in,
        n_fundamental_syn=len(syn.fundamental), n_fundamental_comp=len(comp.fundamental),
        n_distinct_syn_pairs=len(syn_pairs),
        hasse_edges=hasse_edges, n_comparable_pairs=n_comp_pairs,
        n_incomparable_pairs=n_pairs - n_comp_pairs,
    )


# ---------------------------------------------------------------------------
# EPM computation instrumented with degeneracy tracking
# ---------------------------------------------------------------------------

@dataclass
class DegeneracyStats:
    state_records: list[tuple[int, int, float]] = field(default_factory=list)  # (N, distinct, top_share)
    edge_counter: dict[int, int] = field(default_factory=dict)   # target sp -> incoming edge count

    def summary(self) -> dict:
        if not self.state_records:
            return {"n_branching_states": 0}
        ratios = [d / n for (n, d, _) in self.state_records]
        counts = list(self.edge_counter.values())
        return {
            "n_branching_states": len(self.state_records),
            "mean_ratio": sum(ratios) / len(ratios),
            "median_ratio": sorted(ratios)[len(ratios) // 2],
            "n_distinct_targets": len(counts),
            "total_edges": sum(counts),
            "gini": _gini(counts),
            "top20_share": _top_k_share(counts, 0.2),
        }


def _gini(counts: list[int]) -> float:
    if not counts:
        return 0.0
    values = sorted(counts)
    n = len(values)
    total = sum(values)
    if total == 0:
        return 0.0
    weighted = sum((i + 1) * v for i, v in enumerate(values))
    return (2 * weighted) / (n * total) - (n + 1) / n


def _top_k_share(counts: list[int], frac: float) -> float:
    values = sorted(counts, reverse=True)
    total = sum(values)
    if total == 0:
        return 0.0
    k = max(1, int(round(len(values) * frac)))
    return sum(values[:k]) / total


def _instrumented_mode1_dfs(seed_states, g, visited_sp, deg: DegeneracyStats, *, lineage=None, stats=None):
    """Same logic as cot_gen.epm._mode1_dfs, plus per-state raw candidate
    bookkeeping for degeneracy reporting. Never changes which states get
    explored or which SSMs get found -- purely additive instrumentation.

    `lineage`, if given, is a dict[sp, set[str]] of provenance tags (e.g.
    which Mode-2 move type proposed the seed this state descends from).
    It is READ for the state being expanded and WRITTEN (propagated,
    union'd with anything already there) for every child produced -- this
    is necessary because a Mode-2 seed's own .sp is frequently NOT the
    final SSM's .sp (further Option A/B closure can add more ERCs before
    req reaches 0), so tagging only the seed and checking `final.sp in
    seed_moves` silently drops the tag for any SO that needed extra
    Mode-1 chasing beyond its immediate seed.
    """
    ssm_states, leaf_states = [], []
    stack = [s for s in seed_states if s.sp not in visited_sp]
    while stack:
        state = stack.pop()
        if state.sp in visited_sp:
            continue
        visited_sp.add(state.sp)
        if stats is not None:
            stats['states_explored'] = stats.get('states_explored', 0) + 1
        if state.is_ssm:
            ssm_states.append(state)
            if stats is not None:
                stats['ssm_found'] = stats.get('ssm_found', 0) + 1
                by_len = stats.setdefault('ssm_by_length', {})
                by_len[len(state.erc_set)] = by_len.get(len(state.erc_set), 0) + 1
            continue
        made_progress = False
        min_seed = min(state.erc_set)
        raw_targets: list[int] = []
        parent_tags = lineage.get(state.sp, set()) if lineage is not None else None

        for s_bit in _bits(state.req):
            for prod_idx in g.comp_by_species.get(s_bit, []):
                if prod_idx < state.min_ext:
                    if stats is not None:
                        stats['canonical_pruned'] = stats.get('canonical_pruned', 0) + 1
                    continue
                if (g.species_mask[prod_idx] & state.sp) == g.species_mask[prod_idx]:
                    continue
                new_state = g.extend_state(state, prod_idx)
                newly_added = new_state.erc_set - state.erc_set
                if any(i < min_seed and not g.is_persistent[i] for i in newly_added):
                    if stats is not None:
                        stats['canonical_pruned'] = stats.get('canonical_pruned', 0) + 1
                    continue
                raw_targets.append(new_state.sp)
                deg.edge_counter[new_state.sp] = deg.edge_counter.get(new_state.sp, 0) + 1
                if new_state.sp not in visited_sp:
                    stack.append(new_state)
                    made_progress = True
                    if stats is not None:
                        stats['comp_extensions'] = stats.get('comp_extensions', 0) + 1
                    if lineage is not None and parent_tags:
                        lineage.setdefault(new_state.sp, set()).update(parent_tags)

        for erc_i in state.erc_set:
            for (j, k) in g.syn_from.get(erc_i, []):
                if j in state.erc_set:
                    continue
                if j < state.min_ext:
                    if stats is not None:
                        stats['canonical_pruned'] = stats.get('canonical_pruned', 0) + 1
                    continue
                if not (g.prod_mask[k] & state.req):
                    continue
                new_state = g.extend_state(state, j)
                newly_added = new_state.erc_set - state.erc_set
                if any(i < min_seed and not g.is_persistent[i] for i in newly_added):
                    if stats is not None:
                        stats['canonical_pruned'] = stats.get('canonical_pruned', 0) + 1
                    continue
                raw_targets.append(new_state.sp)
                deg.edge_counter[new_state.sp] = deg.edge_counter.get(new_state.sp, 0) + 1
                if new_state.sp not in visited_sp:
                    stack.append(new_state)
                    made_progress = True
                    if stats is not None:
                        stats['syn_extensions'] = stats.get('syn_extensions', 0) + 1
                    if lineage is not None and parent_tags:
                        lineage.setdefault(new_state.sp, set()).update(parent_tags)

        if len(raw_targets) >= 2:
            n = len(raw_targets)
            distinct_targets = set(raw_targets)
            top_count = max(raw_targets.count(t) for t in distinct_targets)
            deg.state_records.append((n, len(distinct_targets), top_count / n))

        if not made_progress:
            leaf_states.append(state)
            if stats is not None:
                stats['leaves_found'] = stats.get('leaves_found', 0) + 1
                by_dead = stats.setdefault('dead_by_length', {})
                by_dead[len(state.erc_set)] = by_dead.get(len(state.erc_set), 0) + 1

    return ssm_states, leaf_states


def compute_epms_instrumented(ercs, hier, syn, comp) -> tuple[EPMResult, DegeneracyStats]:
    """Mirrors cot_gen.epm.compute_epms exactly (Stage 1 + Stage 2 + Stage 3
    order assignment), with degeneracy stats recorded during Mode-1 DFS."""
    from pyCOT.analysis.organizations.epm import _assign_orders

    single_idx, single_masks = _single_erc_epms(ercs, hier)
    single_sp_set = set(single_masks)
    g = FundamentalGraph(ercs, hier, syn, comp)
    visited_sp: set[int] = set(single_masks)
    deg = DegeneracyStats()
    stats: dict = {}

    non_p_seeds = [g.make_seed_state(i) for i in range(len(ercs)) if not ercs[i].is_persistent()]
    ssm_states, leaf_states = _instrumented_mode1_dfs(non_p_seeds, g, visited_sp, deg, stats=stats)

    all_sp_masks = list(single_sp_set) + [s.sp for s in ssm_states]
    so_order = _assign_orders(all_sp_masks, verbose=False)
    for sp in single_masks:
        so_order.setdefault(sp, 0)

    multi_epm_masks = [s.sp for s in ssm_states if so_order.get(s.sp, -1) == 0 and s.sp not in single_sp_set]
    all_epm_masks = sorted(single_sp_set | set(multi_epm_masks))

    sp_to_state = {s.sp: s for s in ssm_states}
    for i in single_idx:
        sp = ercs[i].species_mask
        sp_to_state.setdefault(sp, g.make_seed_state(i))

    result = EPMResult(
        single_epm_indices=single_idx, single_epm_masks=single_masks,
        multi_epm_masks=multi_epm_masks, all_epm_masks=all_epm_masks,
        leaf_masks=[s.sp for s in leaf_states], stats=stats,
        _graph=g, _visited_sp=visited_sp, _sp_to_state=sp_to_state, _so_order=so_order,
    )
    return result, deg


# ---------------------------------------------------------------------------
# ESPM computation instrumented with per-order move-type provenance
# ---------------------------------------------------------------------------

@dataclass
class MoveTypeCounts:
    synergy: int = 0
    complementarity: int = 0
    vertical_lift: int = 0
    carryover: int = 0
    n_new_so: int = 0

    def as_dict(self) -> dict:
        return {"synergy": self.synergy, "complementarity": self.complementarity,
                "vertical_lift": self.vertical_lift, "carryover": self.carryover}


def compute_espm_instrumented(ercs, hier, syn, comp, epm_result: EPMResult, *, max_order=10):
    """Mirrors cot_gen.epm.compute_espm's Mode-2/Mode-1 loop exactly (same
    candidate generation, same seeds reach Mode-1), but tags each surviving
    Mode-2 seed with which move type(s) proposed it, so each order's new SOs
    can be attributed to synergy / complementarity-consumer / vertical-lift
    for reporting."""
    g = epm_result._graph
    visited_sp = epm_result._visited_sp
    sp_to_state = epm_result._sp_to_state
    so_order = epm_result._so_order

    espm_by_order: dict[int, list[int]] = {}
    leaf_masks_by_order: dict[int, list[int]] = {}
    move_counts_by_order: dict[int, MoveTypeCounts] = {}
    deg = DegeneracyStats()

    for sp, ord_k in so_order.items():
        if ord_k >= 1:
            espm_by_order.setdefault(ord_k, []).append(sp)
    for k in espm_by_order:
        espm_by_order[k].sort()
        mc = MoveTypeCounts()
        mc.carryover = len(espm_by_order[k])
        mc.n_new_so = len(espm_by_order[k])
        move_counts_by_order[k] = mc

    current_layer = [sp_to_state[sp] for sp in epm_result.all_epm_masks if sp in sp_to_state]

    for order in range(1, max_order + 1):
        if not current_layer:
            break

        # seed_sp -> set of move types that proposed it this round
        seed_moves: dict[int, set[str]] = {}
        mode2_seeds = []
        raw_targets_by_state: list[list[int]] = []

        for so_state in current_layer:
            raw_targets = []
            for i in so_state.erc_set:
                for (j, _) in g.syn_from.get(i, []):
                    if j in so_state.erc_set:
                        continue
                    ext_state = g.extend_state(so_state, j)
                    raw_targets.append(ext_state.sp)
                    deg.edge_counter[ext_state.sp] = deg.edge_counter.get(ext_state.sp, 0) + 1
                    if ext_state.sp not in visited_sp:
                        mode2_seeds.append(ext_state)
                        seed_moves.setdefault(ext_state.sp, set()).add("synergy")
                for a in g.parents[i]:
                    if a in so_state.erc_set:
                        continue
                    if (g.species_mask[a] & so_state.sp) == g.species_mask[a]:
                        continue
                    ext_state = g.extend_state(so_state, a)
                    raw_targets.append(ext_state.sp)
                    deg.edge_counter[ext_state.sp] = deg.edge_counter.get(ext_state.sp, 0) + 1
                    if ext_state.sp not in visited_sp:
                        mode2_seeds.append(ext_state)
                        seed_moves.setdefault(ext_state.sp, set()).add("vertical_lift")
            for s_bit in _bits(so_state.prod):
                for cons_idx in g.comp_consumers_by_species.get(s_bit, []):
                    if cons_idx in so_state.erc_set:
                        continue
                    if (g.species_mask[cons_idx] & so_state.sp) == g.species_mask[cons_idx]:
                        continue
                    ext_state = g.extend_state(so_state, cons_idx)
                    raw_targets.append(ext_state.sp)
                    deg.edge_counter[ext_state.sp] = deg.edge_counter.get(ext_state.sp, 0) + 1
                    if ext_state.sp not in visited_sp:
                        mode2_seeds.append(ext_state)
                        seed_moves.setdefault(ext_state.sp, set()).add("complementarity")

            if len(raw_targets) >= 2:
                distinct = set(raw_targets)
                top_count = max(raw_targets.count(t) for t in distinct)
                deg.state_records.append((len(raw_targets), len(distinct), top_count / len(raw_targets)))

        if not mode2_seeds:
            break

        new_ssm_states, new_leaf_states = _instrumented_mode1_dfs(
            mode2_seeds, g, visited_sp, deg, lineage=seed_moves)

        # Bucket by the TRUE assigned so_order (max(sub-SO order) + 1), not
        # by the BFS round index -- they usually coincide but can diverge
        # (a state reached in round r can have sub-SOs implying a higher
        # order), exactly as the real compute_espm's own verbose reporting
        # already accounts for via ord_counts. Merge into (not overwrite)
        # whatever that order's MoveTypeCounts already holds from the
        # carry-over pre-population above.
        for s in sorted(new_ssm_states, key=lambda st: bin(st.sp).count('1')):
            if s.sp not in so_order:
                max_sub = -1
                for sub_sp, sub_ord in so_order.items():
                    if sub_sp != s.sp and (sub_sp & s.sp) == sub_sp and sub_ord > max_sub:
                        max_sub = sub_ord
                so_order[s.sp] = max_sub + 1
            sp_to_state[s.sp] = s
            mc = move_counts_by_order.setdefault(so_order[s.sp], MoveTypeCounts())
            # An SO can genuinely be reachable via more than one move type
            # at once (that convergence is itself real, interesting
            # information -- see the degeneracy analysis) but for a
            # readable stacked bar each SO must land in exactly one
            # category, so assign by a fixed priority: complementarity
            # (cheapest, most "fundamental" per the theoretical
            # discussion) first, then vertical lift, then synergy.
            moves = seed_moves.get(s.sp, set())
            if "complementarity" in moves:
                mc.complementarity += 1
            elif "vertical_lift" in moves:
                mc.vertical_lift += 1
            elif "synergy" in moves:
                mc.synergy += 1
            else:
                mc.carryover += 1
            mc.n_new_so += 1

        if new_leaf_states:
            leaf_masks_by_order.setdefault(order, []).extend(s.sp for s in new_leaf_states)

        current_layer = new_ssm_states
        if not new_ssm_states:
            break

    espm_by_order.clear()
    for sp, ord_k in so_order.items():
        if ord_k >= 1:
            espm_by_order.setdefault(ord_k, []).append(sp)
    for k in list(espm_by_order):
        espm_by_order[k].sort()

    result = ESPMResult(
        epm_masks=sorted(epm_result.all_epm_masks),
        espm_by_order=espm_by_order,
        all_so_masks=sorted(so_order.keys()),
        leaf_masks_by_order=leaf_masks_by_order,
        stats_by_order={},
    )
    return result, move_counts_by_order, deg


# ---------------------------------------------------------------------------
# SO lattice structure (containment graph among ALL discovered SOs)
# ---------------------------------------------------------------------------

@dataclass
class SOLattice:
    nodes: list[int]                        # species masks, one per SO
    order_of: dict[int, int]                 # sp -> order
    parents_of: dict[int, list[int]]         # sp -> immediate sub-SO(s) (Hasse edges, order k-1 only)


def build_so_lattice(all_so_masks: list[int], so_order: dict[int, int]) -> SOLattice:
    """Immediate-parent Hasse edges: for each SO of order k, its parents are
    the order-(k-1) SOs that are proper subsets of it. (Not exhaustive
    subset search across all pairs -- restricted to the adjacent order,
    which is what "immediate" sub-SO means for this order definition.)"""
    by_order: dict[int, list[int]] = {}
    for sp in all_so_masks:
        by_order.setdefault(so_order.get(sp, 0), []).append(sp)

    parents_of: dict[int, list[int]] = {sp: [] for sp in all_so_masks}
    for k in sorted(by_order.keys()):
        if k == 0:
            continue
        prev = by_order.get(k - 1, [])
        for sp in by_order[k]:
            parents_of[sp] = [p for p in prev if p != sp and (p & sp) == p]

    return SOLattice(nodes=list(all_so_masks), order_of=dict(so_order), parents_of=parents_of)
