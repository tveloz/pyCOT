"""
org_graph.py — Induced Hasse diagram over ACTUAL organizations (is_organization
== True only, not all semi-organizations), and maximally-divergent root-to-top
chain selection within it.

Containment edges are computed by DIRECT, exhaustive pairwise bitmask-subset
comparison of each organization's FULL species footprint (sp | E_mask |
F_mask) -- not by walking SOLattice.parents_of. That lattice-level structure
only records containment against the immediately preceding cot_gen "order"
(build_so_lattice's own docstring: "restricted to the adjacent order...
not exhaustive subset search across all pairs"), where "order" is a
structural/ERC-combination depth, not species-set size. Two organizations
can be genuine species-set subset/superset of one another with NO
intermediate semi-organization sitting at the adjacent order to link them,
in which case the old order-chain walk silently missed the edge, making a
node LOOK isolated/maximal in the resulting Hasse diagram when it was
actually a proper subset of several other organizations (confirmed
empirically on e_coli_core: an 11-species organization is a species-set
subset of all 12 other organizations found, yet the order-chain walk
reported zero ancestors for it). Species-set inclusion is the only notion
of containment this diagram is meant to depict, so it is computed directly
rather than borrowed from a different, coarser relation.
"""
from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class OrgGraph:
    nodes: list[int]                       # organization sp_masks
    parents_of: dict[int, list[int]]       # induced immediate-parent organizations
    children_of: dict[int, list[int]]      # reverse of parents_of
    rank_of: dict[int, int]                # containment-poset rank (see build_org_graph)


def build_org_graph(org_results: dict, so_lattice=None) -> OrgGraph:
    """org_results: dict[sp_mask -> DecompositionResult], pre-filtered to
    is_organization. `so_lattice` is accepted (and ignored) only to keep
    the existing call signature stable for callers that still pass it.

    rank_of is a CONTAINMENT-POSET rank, deliberately NOT cot_gen's own
    "order" (an ERC-combination-depth notion, unrelated to species-set
    size): 0 for an organization containing no other organization, else
    1 + max(rank of the organizations it contains. Two organizations can
    be genuinely comparable (one a species-subset of the other) while
    sharing the same cot_gen order -- e.g. the bare-food-only organization
    is a subset of every other order-0 organization on e_coli_core -- and
    plotting by cot_gen order then draws a same-level "containment" edge,
    which is geometrically nonsensical for a Hasse diagram (a proper
    subset must render strictly below its superset). Ranking by the
    containment structure itself rules this out by construction: a
    contained organization is always strictly smaller, so processing
    nodes in increasing size order lets rank_of[sp] be computed from
    already-known ranks in one pass, no fixpoint needed.
    """
    nodes = list(org_results.keys())
    full_mask = {sp: (sp | r.E_mask | r.F_mask) for sp, r in org_results.items()}

    parents_of: dict[int, list[int]] = {}
    for sp in nodes:
        m = full_mask[sp]
        contained = [
            other for other in nodes
            if other != sp and full_mask[other] != m and (full_mask[other] & m) == full_mask[other]
        ]
        # Maximal elements of `contained` under species-set inclusion (a is
        # dominated if some other contained organization properly contains it).
        maximal = [
            a for a in contained
            if not any(b != a and (full_mask[a] & full_mask[b]) == full_mask[a] for b in contained)
        ]
        parents_of[sp] = maximal

    children_of: dict[int, list[int]] = {sp: [] for sp in nodes}
    for sp, parents in parents_of.items():
        for p in parents:
            children_of[p].append(sp)

    rank_of: dict[int, int] = {}
    for sp in sorted(nodes, key=lambda s: bin(full_mask[s]).count('1')):
        subs = parents_of[sp]
        rank_of[sp] = 0 if not subs else 1 + max(rank_of[s] for s in subs)

    return OrgGraph(nodes=nodes, parents_of=parents_of, children_of=children_of, rank_of=rank_of)


def enumerate_root_to_leaf_chains(graph: OrgGraph, *, max_chains: int = 5000) -> list[list[int]]:
    """
    All maximal chains: paths from a root (no induced parents) to a leaf (no
    induced children) via induced-parent edges, returned top-down (index 0 =
    the leaf/top node, last index = the root/bottom node) -- convenient for
    the "shared prefix from the top" divergence scoring used by
    select_divergent_chains.
    """
    leaves = [sp for sp in graph.nodes if not graph.children_of.get(sp)]
    chains: list[list[int]] = []

    def dfs(node: int, path_top_down: list[int]):
        if len(chains) >= max_chains:
            return
        parents = graph.parents_of.get(node, [])
        path_top_down = path_top_down + [node]
        if not parents:
            chains.append(path_top_down)
            return
        for p in parents:
            dfs(p, path_top_down)

    for leaf in leaves:
        dfs(leaf, [])
    return chains


def _common_prefix_len(a: list[int], b: list[int]) -> int:
    n = 0
    for x, y in zip(a, b):
        if x != y:
            break
        n += 1
    return n


def select_divergent_chains(chains: list[list[int]], n: int = 5) -> list[list[int]]:
    """
    Greedily pick up to `n` chains (each top-down: [leaf, ..., root]) that
    diverge from each other as close to the top as possible: repeatedly add
    the remaining chain whose worst-case (max) shared-top-prefix length
    against the already-chosen set is smallest, breaking ties in favor of
    longer chains (more levels shown = more informative).
    """
    if not chains:
        return []
    remaining = list(chains)
    # Seed with the longest chain (most informative single chain).
    remaining.sort(key=len, reverse=True)
    chosen = [remaining.pop(0)]

    while remaining and len(chosen) < n:
        def score(c):
            worst = max(_common_prefix_len(c, s) for s in chosen)
            return (worst, -len(c))
        remaining.sort(key=score)
        chosen.append(remaining.pop(0))

    return chosen
