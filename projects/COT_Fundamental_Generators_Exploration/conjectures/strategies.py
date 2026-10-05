"""
strategies.py — the paper's two §6.4 conjectures, factored into directly
comparable strategies over the SAME fundamental-graph search engine
(pyCOT.analysis.organizations.so_search).

Both strategies share every piece of machinery (FundamentalGraph, DFSState,
extend_state, Horn synergy-closure, canonical-ordering dedup) EXCEPT one
move: whether Mode-2 (the outward-growth step run once an SSM / semi-
organization has already been reached) is allowed to use hierarchy vertical
lift. That single flag, `use_vertical_lift` on
so_search.compute_so_hierarchy, IS the entire difference between the two
conjectures:

  Conjecture 1 — run_contained_bfs (use_vertical_lift=True)
    Extend via fundamental synergy and fundamental complementarity, and
    additionally allow upward containment through the ERC hierarchy — but
    only once already at a semi-organization, never mid-construction of one.

  Conjecture 2 — run_independent_seed (use_vertical_lift=False)
    Extend via fundamental synergy and fundamental complementarity ONLY.
    The hierarchy is never touched directly; every semi-organization must be
    reachable purely through synergy/complementarity chains, independently
    seeded per ERC (compute_elementary_sos already seeds Mode-1 DFS from
    every non-persistent ERC, and Mode-2 starts a fresh outward search from
    every elementary SO found — that per-ERC/per-elementary-SO independence
    is what the paper's second conjecture proposes instead of lifting
    through a single generator's construction).

Isolating the difference to one flag means any gap in what the two
strategies find, or in how much search effort they spend, is attributable
to that one conceptual choice — not to two different code paths of
differing quality.

Both `run_contained_bfs` and `run_independent_seed` return a `StrategyRun`
with directly comparable stats (states explored, extension counts by type,
wall-clock time), aggregated across every DFS/BFS round.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field

from pyCOT.analysis.organizations.so_search import (
    compute_elementary_sos,
    compute_so_hierarchy,
)

_STAT_KEYS = (
    'states_explored', 'comp_extensions', 'syn_extensions',
    'canonical_pruned', 'lift_candidates', 'ssm_found', 'leaves_found',
)


@dataclass
class StrategyRun:
    """Uniform result shape for both conjecture strategies."""
    strategy:            str
    all_so_masks:        list = field(default_factory=list)
    so_by_order:         dict = field(default_factory=dict)
    max_order_reached:   int = 0
    n_elementary:        int = 0
    wall_time_s:         float = 0.0

    states_explored:     int = 0
    comp_extensions:     int = 0
    syn_extensions:      int = 0
    lift_candidates:     int = 0   # always 0 for run_independent_seed
    canonical_pruned:    int = 0
    ssm_found:           int = 0
    leaves_found:        int = 0

    def n_so_total(self) -> int:
        return len(self.all_so_masks)


def _aggregate_stats(elementary_result, so_hierarchy_result) -> dict:
    agg = {k: 0 for k in _STAT_KEYS}
    for k in _STAT_KEYS:
        agg[k] += elementary_result.stats.get(k, 0)
    for round_stats in so_hierarchy_result.stats_by_order.values():
        for k in _STAT_KEYS:
            agg[k] += round_stats.get(k, 0)
    return agg


def _run(
    strategy_name: str,
    rn, ercs, hier, syn_result, comp_result,
    *,
    use_vertical_lift: bool,
    max_order: int = 10,
    verbose: bool = False,
) -> StrategyRun:
    t0 = time.perf_counter()
    elementary_result = compute_elementary_sos(
        rn, ercs, hier, syn_result, comp_result, verbose=verbose,
    )
    so_hierarchy_result = compute_so_hierarchy(
        rn, ercs, hier, syn_result, comp_result, elementary_result,
        max_order=max_order, use_vertical_lift=use_vertical_lift, verbose=verbose,
    )
    wall = time.perf_counter() - t0
    agg = _aggregate_stats(elementary_result, so_hierarchy_result)

    return StrategyRun(
        strategy=strategy_name,
        all_so_masks=list(so_hierarchy_result.all_so_masks),
        so_by_order={k: list(v) for k, v in so_hierarchy_result.so_by_order.items()},
        max_order_reached=so_hierarchy_result.max_order(),
        n_elementary=len(elementary_result.all_elementary_masks),
        wall_time_s=wall,
        **agg,
    )


def run_contained_bfs(rn, ercs, hier, syn_result, comp_result, *, max_order=10, verbose=False) -> StrategyRun:
    """Conjecture 1: horizontal moves + hierarchy vertical lift (Mode-2 only)."""
    return _run(
        "contained_bfs", rn, ercs, hier, syn_result, comp_result,
        use_vertical_lift=True, max_order=max_order, verbose=verbose,
    )


def run_independent_seed(rn, ercs, hier, syn_result, comp_result, *, max_order=10, verbose=False) -> StrategyRun:
    """Conjecture 2: horizontal moves only, independently seeded per ERC."""
    return _run(
        "independent_seed", rn, ercs, hier, syn_result, comp_result,
        use_vertical_lift=False, max_order=max_order, verbose=verbose,
    )
