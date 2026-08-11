"""
erc_hierarchy_viz.py — ERC hierarchy + fundamental relations visualization.

Renders the generative structure of a reaction network — the object the
Persistent_Modules_Generator.py module docstring points to as worth
studying on its own, before any organization/EPM/ESPM computation — as an
interactive graph:

  - gray edges   : ERC containment (the Hasse diagram of the hierarchy)
  - orange nodes : fundamental synergies (diamond nodes between the pair
                    of ERCs that jointly unlock a third)
  - blue dashed  : fundamental complementarities (producer ERC -> consumer
                    ERC, labeled with the species being supplied)

Ported from projects/COT_Fundamental_Generators_Exploration/cot_gen/deep_report.py and
deep_report_viz.py (compute_hierarchy_stats + plot_hierarchy_overview),
trimmed to just the hierarchy/relations visualization (the degeneracy-
tracking and SO-lattice plotting in the original stay in that project,
which is a research/instrumentation tool rather than a library-visualization
concern).

Public API
----------
compute_hierarchy_stats(ercs, hier, syn_result, comp_result) -> HierarchyStats
plot_erc_hierarchy(ercs, hier, syn_result, comp_result, out_path, *,
                    stats=None, highlight_ercs=None, highlight_label="highlighted",
                    max_nodes=250, show_synergy=True, show_complementarity=True,
                    title="ERC hierarchy") -> str (path to the written HTML file)

Quick start
-----------
    from pyCOT.io.functions import read_txt
    from pyCOT.analysis.organizations import (
        build_rndata, compute_ercs, build_hierarchy,
        compute_synergies_basis_first, compute_complementarities,
    )
    from pyCOT.visualization.erc_hierarchy_viz import plot_erc_hierarchy

    rn = read_txt("data/Examples_tests/Centler2006_EcoliSugar/centler_glucose.txt")
    rn_data = build_rndata(rn, network_id="centler_glucose")
    ercs = compute_ercs(rn_data)
    hier = build_hierarchy(ercs)
    syn = compute_synergies_basis_first(ercs, hier)
    comp = compute_complementarities(ercs, hier, syn)
    plot_erc_hierarchy(ercs, hier, syn, comp, "erc_hierarchy.html")
"""
from __future__ import annotations

import os
from dataclasses import dataclass

# Shared palette -- consistent with cot_gen/explorer.py's relation colors,
# so plots from either tool read the same way.
COLOR_CONTAINMENT = "#7f8c8d"
COLOR_SYNERGY = "#e67e22"
COLOR_COMPLEMENTARITY = "#2980b9"
COLOR_HIGHLIGHT = "#27ae60"
COLOR_NEUTRAL = "#d0d3d4"


@dataclass
class HierarchyStats:
    """Per-ERC structural statistics, used both for the plot and for reports."""
    n_ercs: int
    n_persistent: int
    sizes: list          # species count per ERC
    levels: list          # containment level per ERC (0 = subset-minimal)
    n_levels: int
    syn_degree: list      # fundamental synergies touching each ERC
    comp_out_degree: list  # as fundamental producer
    comp_in_degree: list   # as fundamental consumer
    n_fundamental_syn: int
    n_fundamental_comp: int
    hasse_edges: int
    n_comparable_pairs: int
    n_incomparable_pairs: int


def _compute_containment_levels(hier, n: int) -> list:
    """Level 0 = subset-minimal ERCs; level(i) = 1 + max(level(child))."""
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


def compute_hierarchy_stats(ercs, hier, syn_result, comp_result) -> HierarchyStats:
    """Compute per-ERC level/degree statistics used by plot_erc_hierarchy."""
    n = len(ercs)
    sizes = [e.size() for e in ercs]
    levels = _compute_containment_levels(hier, n)
    syn_degree = [0] * n
    for st in syn_result.fundamental:
        syn_degree[st.i] += 1
        syn_degree[st.j] += 1
    comp_out = [0] * n
    comp_in = [0] * n
    for fc in comp_result.fundamental:
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
        n_fundamental_syn=len(syn_result.fundamental),
        n_fundamental_comp=len(comp_result.fundamental),
        hasse_edges=hasse_edges, n_comparable_pairs=n_comp_pairs,
        n_incomparable_pairs=n_pairs - n_comp_pairs,
    )


def _new_pyvis_network(height="800px", *, level_separation=140, node_spacing=110):
    from pyvis.network import Network
    net = Network(height=height, width="100%", directed=True, notebook=False, bgcolor="#ffffff")
    net.set_options(f"""
    {{
      "layout": {{"hierarchical": {{"enabled": true, "direction": "UD",
                                    "sortMethod": "hubsize",
                                    "levelSeparation": {level_separation}, "nodeSpacing": {node_spacing}}}}},
      "physics": {{"enabled": false}},
      "interaction": {{"hover": true, "navigationButtons": true, "keyboard": true}}
    }}
    """)
    return net


def _inject_title(html_path: str, title: str) -> None:
    with open(html_path, "r", encoding="utf-8") as f:
        content = f.read()
    banner = (f'<div style="font-family:sans-serif;padding:10px 16px;'
              f'background:#f4f4f4;border-bottom:1px solid #ddd;font-size:16px;">'
              f'<b>{title}</b></div>')
    content = content.replace("<body>", "<body>" + banner, 1)
    with open(html_path, "w", encoding="utf-8") as f:
        f.write(content)


def plot_erc_hierarchy(
    ercs, hier, syn_result, comp_result, out_path: str,
    *,
    stats: HierarchyStats | None = None,
    highlight_ercs=None,
    highlight_label: str = "highlighted",
    max_nodes: int = 250,
    show_synergy: bool = True,
    show_complementarity: bool = True,
    title: str = "ERC hierarchy",
) -> str:
    """
    Render the ERC containment hierarchy with fundamental synergies and
    complementarities overlaid, as an interactive pyvis HTML graph.

    Nodes are leveled by containment depth (subset-minimal ERCs at the
    bottom). Gray edges are containment (child -> parent). Orange diamond
    nodes mark fundamental synergies (edges to both contributing ERCs).
    Blue dashed edges are fundamental complementarities (producer -> consumer).

    Parameters
    ----------
    ercs, hier, syn_result, comp_result : outputs of
        pyCOT.analysis.organizations.compute_ercs / build_hierarchy /
        compute_synergies_basis_first / compute_complementarities
    out_path : where to write the HTML file (parent dirs created if needed)
    stats : precomputed HierarchyStats; computed automatically if omitted
    highlight_ercs : iterable of ERC indices to color distinctly (e.g. the
        members of one particular EPM/ESPM), so their position in the full
        hierarchy is immediately visible
    max_nodes : caps total nodes shown for readability on genome-scale
        networks -- keeps the highlighted set plus a level-balanced sample
    show_synergy, show_complementarity : toggle each relation overlay
    title : banner text written into the HTML page

    Returns
    -------
    str : out_path, for chaining
    """
    if stats is None:
        stats = compute_hierarchy_stats(ercs, hier, syn_result, comp_result)

    n = len(ercs)
    highlight = set(highlight_ercs or [])

    if n > max_nodes:
        by_level: dict = {}
        for i in range(n):
            by_level.setdefault(stats.levels[i], []).append(i)
        keep = set(highlight)
        budget = max_nodes - len(keep)
        levels_sorted = sorted(by_level.keys())
        per_level = max(1, budget // max(1, len(levels_sorted)))
        for lvl in levels_sorted:
            for i in by_level[lvl]:
                if len(keep) >= max_nodes:
                    break
                if i not in keep and len([x for x in keep if stats.levels[x] == lvl]) < per_level:
                    keep.add(i)
        node_set = keep
        capped_note = f" (showing {len(node_set)} of {n} ERCs, level-balanced sample)"
    else:
        node_set = set(range(n))
        capped_note = ""

    net = _new_pyvis_network()
    for i in node_set:
        is_hl = i in highlight
        color = COLOR_HIGHLIGHT if is_hl else COLOR_NEUTRAL
        size = 22 if is_hl else 14
        title_txt = (f"E{i}  |  {stats.sizes[i]} species  |  level {stats.levels[i]}<br>"
                     f"synergy degree: {stats.syn_degree[i]}<br>"
                     f"complementarity out/in: {stats.comp_out_degree[i]}/{stats.comp_in_degree[i]}")
        if is_hl:
            title_txt += f"<br><b>{highlight_label}</b>"
        net.add_node(i, label=f"E{i}", level=stats.levels[i], color=color, size=size,
                     title=title_txt, shape="dot")

    for i in node_set:
        for p in hier.parents[i]:
            if p in node_set:
                net.add_edge(i, p, color=COLOR_CONTAINMENT, width=1.2, arrows="to")

    if show_synergy:
        seen_pairs: dict = {}
        for st in syn_result.fundamental:
            if st.i in node_set and st.j in node_set:
                seen_pairs.setdefault((min(st.i, st.j), max(st.i, st.j)), []).append(st.k)
        for (i, j), targets in seen_pairs.items():
            d_id = f"syn_{i}_{j}"
            lvl = min(stats.levels[i], stats.levels[j])
            net.add_node(d_id, label="+", shape="diamond", color=COLOR_SYNERGY, size=10,
                         level=max(0, lvl - 1),
                         title=f"fundamental synergy E{i}+E{j} -> " + ", ".join(f"E{k}" for k in targets))
            net.add_edge(i, d_id, color=COLOR_SYNERGY, width=1, arrows="")
            net.add_edge(j, d_id, color=COLOR_SYNERGY, width=1, arrows="")

    if show_complementarity:
        for fc in comp_result.fundamental:
            if fc.prod_idx in node_set and fc.cons_idx in node_set:
                net.add_edge(fc.prod_idx, fc.cons_idx, color=COLOR_COMPLEMENTARITY, width=1,
                             dashes=True, arrows="to",
                             title=f"E{fc.prod_idx} supplies species {fc.species} to E{fc.cons_idx}")

    out_dir = os.path.dirname(out_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    net.write_html(out_path, open_browser=False, notebook=False)
    _inject_title(out_path, title + capped_note)
    return out_path
