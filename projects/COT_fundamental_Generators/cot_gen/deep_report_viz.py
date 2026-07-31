"""
deep_report_viz.py — Presentation-quality plots for cot_gen.deep_report.

Hierarchy-style graphs (pyvis, interactive HTML) reuse the same visual
language as cot_gen/explorer.py (gray containment, orange synergy diamonds,
blue dashed complementarity) for consistency across the tool suite.
Statistical plots (matplotlib, static PNG) follow a fixed categorical color
order and single-hue sequential encodings, never a rainbow palette.
"""
from __future__ import annotations

import os
import webbrowser

from .deep_report import (
    COLOR_CONTAINMENT, COLOR_SYNERGY, COLOR_COMPLEMENTARITY,
    COLOR_EPM, COLOR_NEUTRAL, COLOR_MAXSO,
    MOVE_TYPE_ORDER, MOVE_TYPE_COLOR,
    _gini, _top_k_share,
)


# ---------------------------------------------------------------------------
# Shared pyvis scaffold
# ---------------------------------------------------------------------------

def _new_network(height="800px", *, level_separation=140, node_spacing=110,
                  sort_method="hubsize", direction="DU"):
    """
    direction="DU" is kept as the default to preserve plot_hierarchy_overview's
    existing look (unchanged/not reported as wrong). plot_so_lattice instead
    passes direction="UD" (vis-network's well-documented default: level 0 at
    the TOP, increasing levels downward) with EXPLICITLY INVERTED level
    numbers (level = max_level - k) -- "DU"'s exact vis.js semantics turned
    out to be easy to get backwards in practice (that's what produced the
    reversed so_lattice levels this fixes), so the more reliable fix is to
    stop depending on it there rather than trust it silently elsewhere too.
    """
    from pyvis.network import Network
    net = Network(height=height, width="100%", directed=True, notebook=False, bgcolor="#ffffff")
    net.set_options(f"""
    {{
      "layout": {{"hierarchical": {{"enabled": true, "direction": "{direction}",
                                    "sortMethod": "{sort_method}",
                                    "levelSeparation": {level_separation}, "nodeSpacing": {node_spacing}}}}},
      "physics": {{"enabled": false}},
      "interaction": {{"hover": true, "navigationButtons": true, "keyboard": true}}
    }}
    """)
    return net


def _open(path):
    try:
        webbrowser.open("file://" + os.path.abspath(path))
    except Exception:
        pass


# ---------------------------------------------------------------------------
# ERC hierarchy overview (optionally EPM-highlighted)
# ---------------------------------------------------------------------------

def plot_hierarchy_overview(ercs, hier, syn, comp, stats, out_path, *,
                             highlight_ercs=None, highlight_label="EPM member",
                             max_nodes=250, show_synergy=True, show_complementarity=True,
                             title="ERC hierarchy"):
    """
    Whole-network ERC hierarchy: nodes leveled by containment depth, gray
    containment edges, orange synergy diamonds, blue dashed complementarity
    edges. If `highlight_ercs` is given, those nodes are colored distinctly
    (green) with a size bump, so their position within the full hierarchy
    is immediately visible.

    Caps total nodes at `max_nodes` for readability on genome-scale
    networks -- keeps the highlighted set (if any) plus a level-balanced
    sample of the rest.
    """
    n = len(ercs)
    highlight = set(highlight_ercs or [])

    if n > max_nodes:
        by_level: dict[int, list[int]] = {}
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
                if i not in keep and len(keep) < max_nodes:
                    if len([x for x in keep if stats.levels[x] == lvl]) < per_level:
                        keep.add(i)
        node_set = keep
        capped_note = f" (showing {len(node_set)} of {n} ERCs, level-balanced sample)"
    else:
        node_set = set(range(n))
        capped_note = ""

    net = _new_network()
    for i in node_set:
        is_hl = i in highlight
        color = COLOR_EPM if is_hl else COLOR_NEUTRAL
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
        seen_pairs: dict[tuple[int, int], list[int]] = {}
        for st in syn.fundamental:
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
        for fc in comp.fundamental:
            if fc.prod_idx in node_set and fc.cons_idx in node_set:
                net.add_edge(fc.prod_idx, fc.cons_idx, color=COLOR_COMPLEMENTARITY, width=1,
                             dashes=True, arrows="to",
                             title=f"E{fc.prod_idx} supplies species {fc.species} to E{fc.cons_idx}")

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    net.write_html(out_path, open_browser=False, notebook=False)
    _inject_title(out_path, title + capped_note)
    return out_path


def _inject_title(html_path, title):
    with open(html_path, "r", encoding="utf-8") as f:
        content = f.read()
    banner = (f'<div style="font-family:sans-serif;padding:10px 16px;'
              f'background:#f4f4f4;border-bottom:1px solid #ddd;font-size:16px;">'
              f'<b>{title}</b></div>')
    content = content.replace("<body>", "<body>" + banner, 1)
    with open(html_path, "w", encoding="utf-8") as f:
        f.write(content)


# ---------------------------------------------------------------------------
# SO lattice (Hasse diagram over discovered semi-organizations)
# ---------------------------------------------------------------------------

# Above this many nodes, vis-network's browser-side layout (hierarchical
# edge-crossing minimization, still O(n*e)-ish even with physics off) stops
# finishing in any reasonable time -- empirically confirmed hung/blank at
# 4244 nodes / 24267 edges, fine at 723 / 3641. Past this threshold we skip
# the interactive pyvis HTML entirely and render a static, physics-free
# layered PNG instead (matplotlib LineCollection scales to tens of
# thousands of edges in seconds).
SO_LATTICE_INTERACTIVE_MAX_NODES = 800


def _order_shade(k: int, max_order: int) -> str:
    """Sequential single-hue ramp (light -> dark), order 0 lightest."""
    t = k / max(max_order, 1)
    r0, g0, b0 = 0xfd, 0xeb, 0xd0  # light gold
    r1, g1, b1 = 0x9c, 0x64, 0x0c  # deep amber
    r = int(r0 + (r1 - r0) * t)
    g = int(g0 + (g1 - g0) * t)
    b = int(b0 + (b1 - b0) * t)
    return f"#{r:02x}{g:02x}{b:02x}"


def _adaptive_lattice_layout(n_shown: int) -> dict:
    """
    Node size and spacing that scale with how many nodes actually have to
    fit in the view, instead of one fixed formula for every network size.
    Area-preserving inverse-sqrt scaling: a lattice with 4x as many nodes
    gets base dots half the diameter, so total ink stays roughly constant
    and containment structure stays legible whether the network has a
    dozen SOs or several hundred.
    """
    scale = 1.0 / max(1.0, n_shown) ** 0.5
    return {
        "base_size": max(7, min(38, 260 * scale)),
        "level_separation": max(70, min(260, 4000 * scale)),
        "node_spacing": max(50, min(200, 3000 * scale)),
        "edge_width": max(0.5, min(2.2, 180 * scale / max(1, n_shown) ** 0.15)),
    }


def plot_so_lattice(so_lattice, rn, out_path, *, max_order_shown=None, title="Semi-organization lattice"):
    """
    One node per discovered SO (species mask), leveled by order (order 0 at
    the bottom, increasing upward). Edges connect each SO to its immediate
    order-(k-1) sub-SO(s), drawn as black dotted lines so containment stays
    visible against the order-shaded nodes. Node size scales with species
    count relative to the largest SO shown, on top of a BASE size that
    itself adapts to how many nodes are in the view (_adaptive_lattice_layout)
    -- a network with a handful of SOs gets large, clearly-separated dots;
    one with hundreds gets smaller dots and tighter spacing so it still fits,
    rather than one fixed size that is too small for small networks and too
    crowded for large ones.

    Always renders the static, physics-free PNG (matplotlib -- fast at any
    size, and directly viewable without a browser). Additionally renders the
    interactive pyvis HTML when the lattice is small enough for the browser
    to lay out (see SO_LATTICE_INTERACTIVE_MAX_NODES) -- above that, the
    browser-side layout hangs/blank-renders, so it's skipped entirely.

    Returns a list of the path(s) written (PNG always; HTML too if small).
    """
    orders_present = sorted(set(so_lattice.order_of.values()))
    if max_order_shown is not None:
        orders_present = [o for o in orders_present if o <= max_order_shown]
    max_order = max(orders_present) if orders_present else 0

    shown = [sp for sp in so_lattice.nodes if so_lattice.order_of.get(sp, 0) in orders_present]

    paths = [_plot_so_lattice_static(so_lattice, shown, max_order, out_path, title=title)]

    if len(shown) > SO_LATTICE_INTERACTIVE_MAX_NODES:
        return paths

    layout = _adaptive_lattice_layout(len(shown))
    base_size = layout["base_size"]
    max_sp_shown = max((bin(sp).count('1') for sp in shown), default=1)

    net = _new_network(height="850px", direction="UD", sort_method="directed",
                        level_separation=layout["level_separation"],
                        node_spacing=layout["node_spacing"])

    for sp in shown:
        k = so_lattice.order_of.get(sp, 0)
        n_sp = bin(sp).count('1')
        size_frac = n_sp / max_sp_shown if max_sp_shown else 0.0
        size = base_size * (1 + 0.7 * size_frac)
        # Invert: order 0 gets the HIGHEST level number, which direction="UD"
        # places at the bottom (level 0 = top, increasing levels go down).
        net.add_node(sp, label=f"o{k}", level=max_order - k, color=_order_shade(k, max_order),
                     size=size, shape="dot",
                     title=f"order {k}  |  {n_sp} species")

    shown_set = set(shown)
    for sp in shown:
        for parent_sp in so_lattice.parents_of.get(sp, []):
            if parent_sp in shown_set:
                net.add_edge(parent_sp, sp, color="#000000", width=layout["edge_width"],
                             dashes=True, arrows="to")

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    net.write_html(out_path, open_browser=False, notebook=False)
    note = f" ({len(shown)} SOs, orders 0-{max_order}, order 0 at bottom)"
    _inject_title(out_path, title + note)
    paths.append(out_path)
    return paths


def _plot_so_lattice_static(so_lattice, shown, max_order, out_path, *, title):
    """
    Physics-free fallback for large lattices: deterministic layered layout
    (y = order, x = arbitrary spread within order), edges drawn as a single
    low-alpha LineCollection. Scales to tens of thousands of edges in
    seconds, unlike browser-side force/hierarchical layout.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    shown_set = set(shown)
    by_level: dict[int, list[int]] = {}
    for sp in shown:
        by_level.setdefault(so_lattice.order_of.get(sp, 0), []).append(sp)

    pos: dict[int, tuple[float, float]] = {}
    for lvl, sps in by_level.items():
        sps_sorted = sorted(sps)
        count = len(sps_sorted)
        for i, sp in enumerate(sps_sorted):
            x = (i - count / 2) / max(count, 1) * max(count ** 0.5, 1) * 3.0
            pos[sp] = (x, lvl)

    segs = []
    for sp in shown:
        for parent_sp in so_lattice.parents_of.get(sp, []):
            if parent_sp in shown_set:
                segs.append([pos[parent_sp], pos[sp]])

    layout = _adaptive_lattice_layout(len(shown))
    max_sp_shown = max((bin(sp).count('1') for sp in shown), default=1)
    # Marker area (matplotlib `s`) scales as size^2, so square the same
    # adaptive base used by the interactive view to keep the two visually
    # consistent rather than picking an unrelated formula here.
    base_area = layout["base_size"] ** 2 / 3.5
    edge_alpha = max(0.08, min(0.6, 60.0 / max(1, len(segs)) ** 0.5))

    fig, ax = plt.subplots(figsize=(24, 14))
    ax.add_collection(LineCollection(segs, colors="#000000", linewidths=layout["edge_width"] * 0.6,
                                      linestyles="dotted", alpha=edge_alpha, zorder=1))

    xs = [pos[sp][0] for sp in shown]
    ys = [pos[sp][1] for sp in shown]
    colors = [_order_shade(so_lattice.order_of.get(sp, 0), max_order) for sp in shown]
    sizes = [base_area * (1 + 0.7 * (bin(sp).count('1') / max_sp_shown if max_sp_shown else 0.0))
             for sp in shown]
    ax.scatter(xs, ys, c=colors, s=sizes, zorder=2, edgecolors="none")

    ax.set_xlabel("(arbitrary horizontal spread within each order)")
    ax.set_ylabel("order (0 = EPM, minimal organizations, at the bottom)")
    ax.set_yticks(range(0, max_order + 1))
    ax.set_xticks([])
    ax.set_title(f"{title}\n({len(shown)} SOs, orders 0-{max_order}; {len(segs)} containment edges; "
                 f"static layered layout, no physics -- too large for interactive rendering)")

    png_path = os.path.splitext(out_path)[0] + ".png"
    os.makedirs(os.path.dirname(png_path), exist_ok=True)
    plt.tight_layout()
    plt.savefig(png_path, dpi=150)
    plt.close(fig)
    return png_path


# ---------------------------------------------------------------------------
# Statistical plots (matplotlib)
# ---------------------------------------------------------------------------

def plot_degeneracy(deg_stats, out_path, *, title="Generative degeneracy"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    state_records = deg_stats.state_records
    counts = list(deg_stats.edge_counter.values())
    if not state_records or not counts:
        return None

    ratios = [d / n for (n, d, _) in state_records]
    sorted_counts = sorted(counts, reverse=True)
    gini = _gini(counts)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    fig.suptitle(title, fontsize=13, fontweight="bold")

    axes[0].hist(ratios, bins=30, color=COLOR_COMPLEMENTARITY, edgecolor="white")
    axes[0].set_xlabel("distinct targets / candidates tried (1.0 = no convergence)")
    axes[0].set_ylabel("branching states")
    axes[0].set_title("Per-state convergence")

    ranks = list(range(1, len(sorted_counts) + 1))
    axes[1].loglog(ranks, sorted_counts, marker="o", markersize=3, linestyle="none", color=COLOR_MAXSO)
    axes[1].set_xlabel("rank of target closure (log)")
    axes[1].set_ylabel("incoming candidate edges (log)")
    axes[1].set_title("Rank-frequency (power-law check)")

    total = sum(sorted_counts)
    cum = 0.0
    xs, ys = [0.0], [0.0]
    for i, c in enumerate(sorted_counts, start=1):
        cum += c
        xs.append(i / len(sorted_counts))
        ys.append(cum / total if total else 0)
    axes[2].plot(xs, ys, color=COLOR_EPM, linewidth=2, label="observed")
    axes[2].plot([0, 1], [0, 1], "--", color="#95a5a6", linewidth=1, label="perfect equality")
    axes[2].set_xlabel("fraction of distinct target closures")
    axes[2].set_ylabel("cumulative fraction of edges")
    axes[2].set_title(f"Lorenz curve (Gini={gini:.3f})")
    axes[2].legend(frameon=False)

    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.tight_layout()
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    return out_path


def plot_espm_composition(move_counts_by_order, out_path, *, title="ESPM construction by order"):
    """Stacked bar: x=order, y=new SOs, stacked by which move type(s)
    contributed. Fixed categorical color order (synergy, complementarity,
    vertical_lift), never reordered."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    orders = sorted(move_counts_by_order.keys())
    if not orders:
        return None

    fig, ax = plt.subplots(figsize=(max(6, len(orders) * 1.1), 5))
    bottoms = [0] * len(orders)
    for move in MOVE_TYPE_ORDER:
        values = [getattr(move_counts_by_order[k], move) for k in orders]
        ax.bar([str(k) for k in orders], values, bottom=bottoms,
               label=move.replace("_", " "), color=MOVE_TYPE_COLOR[move],
               edgecolor="white", linewidth=1)
        bottoms = [b + v for b, v in zip(bottoms, values)]

    ax.set_xlabel("order (k)")
    ax.set_ylabel("new semi-organizations")
    ax.set_title(title)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False, title="reached via")

    fig.tight_layout()
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    return out_path


def plot_order_size_distribution(so_lattice, out_path, *, title="Semi-organization sizes by order"):
    """Box-plot-style: species count distribution per order. Sequential
    single-hue by order (magnitude encoding, light->dark)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    by_order: dict[int, list[int]] = {}
    for sp, k in so_lattice.order_of.items():
        by_order.setdefault(k, []).append(bin(sp).count('1'))
    orders = sorted(by_order.keys())
    if not orders:
        return None

    fig, ax = plt.subplots(figsize=(max(6, len(orders) * 0.9), 5))
    data = [by_order[k] for k in orders]
    bp = ax.boxplot(data, tick_labels=[str(k) for k in orders], patch_artist=True,
                     medianprops={"color": "#2c3e50"})
    max_k = max(orders)
    for i, box in enumerate(bp["boxes"]):
        t = orders[i] / max(max_k, 1)
        r = int(0xfd + (0x9c - 0xfd) * t)
        g = int(0xeb + (0x64 - 0xeb) * t)
        b = int(0xd0 + (0x0c - 0xd0) * t)
        box.set_facecolor(f"#{r:02x}{g:02x}{b:02x}")

    ax.set_xlabel("order (k)")
    ax.set_ylabel("species count")
    ax.set_title(title)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    return out_path
