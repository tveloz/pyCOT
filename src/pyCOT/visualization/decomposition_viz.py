"""
decomposition_viz.py — Presentation-quality plots for pyCOT.analysis.decomposition.

Views:

  plot_organization_hasse / plot_organization_chains
      How |E|, |F|, and fragile-circuit membership evolve across the
      EPM/ESPM hierarchy's organizations (Hasse diagram and top-5
      most-divergent root-to-leaf chains), each organization drawn as an
      E/F/circuit stacked bar with containment arrows tracking which
      group remains/expands vs. gets decomposed going up.

  plot_decomposition_evolution(results, so_lattice, ...)
      How |E|, |F|, and total fragile-circuit species evolve across the
      EPM/ESPM hierarchy's orders (stacked bar), plus the fraction of SOs
      that are full organizations per order.

Follows the same visual conventions as cot_gen/deep_report_viz.py: matplotlib
for statistics, fixed categorical color order, single-hue sequential ramps
for magnitude — never a rainbow.
"""
from __future__ import annotations

import os
import webbrowser


# ---------------------------------------------------------------------------
# Color scheme (fixed, never cycled arbitrarily)
# ---------------------------------------------------------------------------

COLOR_CATALYST = "#16a085"        # E   -- teal: inert, structural
COLOR_OVERPRODUCED = "#e67e22"    # F   -- orange: abundant/generative
COLOR_FAIL_BORDER = "#e74c3c"     # red border overlay: circuit fails self-maintenance
COLOR_OK_BORDER = "#2c3e50"

# Fixed-order categorical palette for individual fragile circuits (D1, D2, ...).
# Cycles (with a neutral fallback) only past this many distinct circuits.
CIRCUIT_PALETTE = [
    "#8e44ad", "#2980b9", "#27ae60", "#d35400", "#c0392b",
    "#16a085", "#2c3e50", "#f39c12", "#7f8c8d", "#1abc9c",
]
COLOR_CIRCUIT_OVERFLOW = "#34495e"


def _circuit_color(i: int) -> str:
    if i < len(CIRCUIT_PALETTE):
        return CIRCUIT_PALETTE[i]
    return COLOR_CIRCUIT_OVERFLOW


def _mask_to_indices(mask: int) -> list[int]:
    out = []
    m = mask
    while m:
        lsb = m & (-m)
        out.append(lsb.bit_length() - 1)
        m &= m - 1
    return out


def _total_mask(r) -> int:
    """E | F | every circuit's species -- equals X_full (sp_mask | E0_mask) by
    construction (structural sanity, validated in tests/validate_decomposition.py)."""
    m = r.E_mask | r.F_mask
    for c in r.circuits:
        m |= c.species_mask
    return m


def _total_n(r) -> int:
    return bin(_total_mask(r)).count('1')


def _node_segments(r):
    """[(label, color, count, species_mask, is_fail_border)] for one
    DecompositionResult, in fixed stacking order: E, F, then circuits
    largest-first (numbered C1, C2, ... by that rank)."""
    segs = []
    if r.E_mask:
        segs.append(("E", COLOR_CATALYST, bin(r.E_mask).count('1'), r.E_mask, False))
    if r.F_mask:
        segs.append(("F", COLOR_OVERPRODUCED, bin(r.F_mask).count('1'), r.F_mask, False))
    for i, c in enumerate(sorted(r.circuits, key=lambda c: -c.size())):
        segs.append((f"C{i+1}", _circuit_color(i), c.size(), c.species_mask, not c.is_self_maintaining))
    return segs


def _open(path):
    try:
        webbrowser.open("file://" + os.path.abspath(path))
    except Exception:
        pass


# ---------------------------------------------------------------------------
# Evolution across the hierarchy
# ---------------------------------------------------------------------------

def _draw_stacked_bar(ax, x_left, y_center, total_width, segments, height, *, label_min_frac=0.13):
    """
    Draw one HORIZONTAL stacked bar starting at x_left, vertically centered
    on y_center, total_width wide, proportioned left-to-right by
    segments=[(label,color,count,mask,is_fail)] (E, F, C1, C2, ... in that
    fixed order — see _node_segments). Horizontal stacking (rather than
    vertical) is deliberate: once every organization is a left-to-right
    bar, containment/lineage arrows between organizations read naturally
    as left-to-right connectors between blocks, instead of needing to
    dodge around a vertical bar's own stacking axis.

    Returns dict[label -> (x0, x1, color)] for ribbon-connector use by
    callers that track lineage across bars (plot_organization_chains).
    """
    import matplotlib.pyplot as plt

    total_count = sum(s[2] for s in segments) or 1
    spans = {}
    x = x_left
    y0 = y_center - height / 2
    for label, color, count, mask, is_fail in segments:
        w = total_width * (count / total_count)
        edge = COLOR_FAIL_BORDER if is_fail else "#ffffff"
        lw = 2.0 if is_fail else 0.6
        ax.add_patch(plt.Rectangle((x, y0), w, height, facecolor=color,
                                    edgecolor=edge, linewidth=lw, zorder=3))
        if w / max(total_width, 1e-9) >= label_min_frac:
            ax.annotate(f"{label}\n{count}", (x + w / 2, y_center), ha="center", va="center",
                        fontsize=10, color="white", fontweight="bold", zorder=4, clip_on=True)
        spans[label] = (x, x + w, color)
        x += w
    return spans


def _segment_containment_edges(child_segments, parent_segments):
    """
    Overlap at the DECOMPOSITION-GROUP level, not the whole-organization
    level: for each child segment (E, F, or one fragile circuit), find
    EVERY parent segment its species have a NON-NULL INTERSECTION with --
    not only a segment that fully contains it.

    Since a parent's segments partition X_full_parent (E/F/circuits are
    pairwise disjoint by construction), the child segment's species are
    partitioned across exactly the parent segments returned here: if it
    lands entirely in one, that is the only edge (the "remains/simply
    expands" case); if it gets redistributed across >=2 parent segments
    going up (part absorbed into F, part remaining its own circuit, etc.),
    one edge is drawn to EACH destination instead of none, so the diagram
    can trace partial absorption rather than silently dropping it. Two
    DIFFERENT child segments both overlapping the SAME parent segment is a
    merge, and both legitimately get an edge into it.

    Returns list[(child_label, parent_label)], possibly several parent
    labels per child label.
    """
    edges = []
    for clabel, ccolor, ccount, cmask, cfail in child_segments:
        if cmask == 0:
            continue
        for plabel, pcolor, pcount, pmask, pfail in parent_segments:
            if pmask and (cmask & pmask):
                edges.append((clabel, plabel))
    return edges


def _draw_segment_arrow(ax, xy_from, xy_to, color, *, lw=1.3, alpha=0.8):
    ax.annotate("", xy=xy_to, xytext=xy_from,
                arrowprops=dict(arrowstyle="-|>", color=color, lw=lw, alpha=alpha,
                                 shrinkA=0, shrinkB=0, mutation_scale=9),
                zorder=2)


def _org_legend_handles():
    import matplotlib.patches as mpatches
    return [
        mpatches.Patch(color=COLOR_CATALYST, label="E (catalysts)"),
        mpatches.Patch(color=COLOR_OVERPRODUCED, label="F (overproduced)"),
        mpatches.Patch(color=CIRCUIT_PALETTE[0], label="fragile circuits C1, C2, ... (size-ranked; red border = fails self-maintenance)"),
    ]


# ---------------------------------------------------------------------------
# Hasse diagram of ACTUAL organizations (is_organization == True only),
# each node rendered as an E/F/circuit stacked-bar glyph
# ---------------------------------------------------------------------------

def plot_organization_hasse(org_results, so_lattice, out_path, *,
                             title="Hasse diagram of organizations",
                             max_bar_width=1.6, bar_height=0.75, level_spacing=2.4):
    """
    org_results: dict[sp_mask -> DecompositionResult], PRE-FILTERED to
    is_organization == True (this is a diagram of organizations, not all
    semi-organizations -- see plot_decomposition_evolution / the SO lattice
    for the full picture).

    Containment edges skip non-organization intermediates: node A's parents
    here are the maximal organization ancestors of A (decomp.org_graph).

    Layout is the conventional Hasse-diagram orientation: level runs
    VERTICALLY (rank 0 -- organizations containing no other organization --
    at the bottom, increasing upward), with same-level nodes spread out
    horizontally to avoid overlap. Each node is still drawn as a HORIZONTAL
    E|F|C1|C2... stacked bar -- the node's own internal orientation is
    independent of the diagram's vertical axis, and conflating the two
    (level AND composition both reading left-to-right) is exactly what
    made the previous version unreadable.

    The vertical level is the CONTAINMENT-POSET rank (decomp.org_graph's
    rank_of), not cot_gen's own ERC-combination "order": the two need not
    agree (e.g. the bare-food-only organization can be a species-subset of
    several other same-order organizations), and plotting by cot_gen order
    can then draw a same-level "containment" edge, which renders as a
    geometrically nonsensical near-horizontal line for what is supposed to
    be a strict subset relation. Ranking by containment itself rules this
    out: a contained organization is always placed strictly below.

    Containment is drawn PER DECOMPOSITION GROUP (decomp.viz._segment_
    containment_edges), not per whole organization: an arrow runs from a
    lower node's E/F/circuit segment up to EVERY parent segment its species
    have non-null overlap with -- i.e. where that group's species end up,
    including when they get redistributed across >=2 parent segments going
    up (part absorbed into F, part remaining a smaller circuit, etc.), in
    which case one arrow is drawn to each destination. Two lower segments
    both overlapping the same upper segment (a merge) both get arrows into it.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if not org_results:
        return None

    from pyCOT.analysis.decomposition.org_graph import build_org_graph
    graph = build_org_graph(org_results, so_lattice)

    by_level: dict[int, list[int]] = {}
    for sp in org_results:
        by_level.setdefault(graph.rank_of[sp], []).append(sp)
    levels = sorted(by_level.keys())
    level_max = {k: max(_total_n(org_results[sp]) for sp in nodes) for k, nodes in by_level.items()}

    pos_x: dict[int, float] = {}
    for k in levels:
        nodes = sorted(by_level[k], key=lambda sp: -_total_n(org_results[sp]))
        n = len(nodes)
        spread = max(1.8, n * 1.3)
        for i, sp in enumerate(nodes):
            pos_x[sp] = (i - (n - 1) / 2.0) * (spread / max(n, 1))

    bar_spans: dict[int, tuple[float, float]] = {}
    node_y: dict[int, float] = {}
    for k in levels:
        y0 = k * level_spacing
        for sp in by_level[k]:
            w = max_bar_width * (_total_n(org_results[sp]) / level_max[k])
            cx = pos_x[sp]
            bar_spans[sp] = (cx - w / 2, cx + w / 2)
            node_y[sp] = y0

    fig, ax = plt.subplots(figsize=(max(11, max(len(v) for v in by_level.values()) * 2.2),
                                     max(7, len(levels) * 2.4)))

    segments_of = {sp: _node_segments(r) for sp, r in org_results.items()}
    spans_of: dict[int, dict[str, tuple[float, float, str]]] = {}
    for sp, r in org_results.items():
        x0, x1 = bar_spans[sp]
        spans_of[sp] = _draw_stacked_bar(ax, x0, node_y[sp], x1 - x0, segments_of[sp], bar_height)
        ax.annotate(f"{_total_n(r)} sp.", ((x0 + x1) / 2, node_y[sp] + bar_height / 2),
                    textcoords="offset points", xytext=(0, 4),
                    ha="center", va="bottom", fontsize=10, zorder=4)

    # decomp.org_graph's "parents_of" means the maximal organizations
    # CONTAINED WITHIN this one (by direct species-set inclusion), i.e. the
    # node(s) BELOW in this rank-ordered vertical layout -- so containment
    # is checked, and arrows drawn, from the smaller precursor UP to the
    # bigger node.
    for sp in org_results:
        for precursor in graph.parents_of.get(sp, []):
            below_spans, above_spans = spans_of[precursor], spans_of[sp]
            for below_label, above_label in _segment_containment_edges(segments_of[precursor], segments_of[sp]):
                bx0, bx1, bcolor = below_spans[below_label]
                ax0, ax1, _acolor = above_spans[above_label]
                _draw_segment_arrow(
                    ax,
                    ((bx0 + bx1) / 2, node_y[precursor] + bar_height / 2),
                    ((ax0 + ax1) / 2, node_y[sp] - bar_height / 2),
                    bcolor,
                )

    ax.set_yticks([k * level_spacing for k in levels])
    ax.set_yticklabels([str(k) for k in levels], fontsize=12)
    ax.set_ylabel("containment rank (0 = contains no other organization)", fontsize=13)
    ax.set_xticks([])
    ax.set_xlim(min(pos_x.values(), default=-1) - max_bar_width - 1.0,
                max(pos_x.values(), default=1) + max_bar_width + 1.0)
    ax.set_ylim(-level_spacing * 0.6, max(levels) * level_spacing + level_spacing * 0.6)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.legend(handles=_org_legend_handles(), loc="upper left", fontsize=11, frameon=False)
    ax.set_title(f"{title} ({len(org_results)} organizations)\n"
                 f"an arrow means that group remains/expands above -- "
                 f"no arrow means it was decomposed there", fontsize=13)

    fig.tight_layout()
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


# ---------------------------------------------------------------------------
# Same Hasse-diagram rendering as plot_organization_hasse, restricted to a
# handful of maximally-divergent root-to-top organization chains
# ---------------------------------------------------------------------------

def plot_organization_chains(org_results, so_lattice, out_path, *, n_chains=5,
                              title="Organization chains", **hasse_kwargs):
    """
    Selects up to `n_chains` root-to-top organization chains that diverge
    from each other as close to the top as possible (decomp.org_graph), then
    renders their union with EXACTLY plot_organization_hasse's layout
    (containment rank VERTICAL, bottom=organizations containing nothing
    else, same-rank nodes spread horizontally, each node an E/F/circuit
    stacked bar, per-segment overlap arrows) -- same diagram, just
    restricted to a few illustrative paths instead of every organization,
    so it doesn't get overwhelmed by the full lattice's branching on larger
    networks.
    """
    from pyCOT.analysis.decomposition.org_graph import build_org_graph, enumerate_root_to_leaf_chains, select_divergent_chains

    if not org_results:
        return None

    graph = build_org_graph(org_results, so_lattice)
    all_chains = enumerate_root_to_leaf_chains(graph)
    chosen = select_divergent_chains(all_chains, n_chains)
    if not chosen:
        return None

    keep = {sp for chain in chosen for sp in chain}
    filtered = {sp: org_results[sp] for sp in keep}
    chain_title = f"{title}\n({len(chosen)} maximally-divergent root-to-top chains)"
    return plot_organization_hasse(filtered, so_lattice, out_path, title=chain_title, **hasse_kwargs)


def plot_decomposition_evolution(results, so_lattice, out_path, *,
                                  title="Decomposition structure vs. hierarchy order"):
    """
    Left panel: stacked bar of avg |E|, avg |F|, avg total-circuit-species
    per order (fixed categorical colors, matching the E/F/circuit roles
    used throughout this module).
    Right panel: fraction of SOs that are full organizations per order
    (bar; a single sequential-magnitude series, not a second axis on the
    same panel as the left one).
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    by_order: dict[int, list] = {}
    for sp, r in results.items():
        by_order.setdefault(so_lattice.order_of.get(sp, 0), []).append(r)
    orders = sorted(by_order.keys())
    if not orders:
        return None

    avg_e, avg_f, avg_c, org_frac = [], [], [], []
    for k in orders:
        rs = by_order[k]
        n = len(rs)
        avg_e.append(sum(bin(r.E_mask).count('1') for r in rs) / n)
        avg_f.append(sum(bin(r.F_mask).count('1') for r in rs) / n)
        avg_c.append(sum(sum(c.size() for c in r.circuits) for r in rs) / n)
        org_frac.append(sum(1 for r in rs if r.is_organization) / n)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle(title, fontsize=13, fontweight="bold")

    xs = [str(k) for k in orders]
    ax0 = axes[0]
    ax0.bar(xs, avg_e, label="E (catalysts)", color=COLOR_CATALYST, edgecolor="white")
    ax0.bar(xs, avg_f, bottom=avg_e, label="F (overproduced)", color=COLOR_OVERPRODUCED, edgecolor="white")
    bottoms2 = [e + f for e, f in zip(avg_e, avg_f)]
    ax0.bar(xs, avg_c, bottom=bottoms2, label="fragile circuits", color=CIRCUIT_PALETTE[0], edgecolor="white")
    ax0.set_xlabel("order (k)")
    ax0.set_ylabel("avg species count")
    ax0.set_title("Composition")
    ax0.legend(frameon=False)
    ax0.spines["top"].set_visible(False)
    ax0.spines["right"].set_visible(False)

    ax1 = axes[1]
    max_k = max(orders)
    colors = []
    for k in orders:
        t = k / max(max_k, 1)
        r = int(0xfd + (0x9c - 0xfd) * t)
        g = int(0xeb + (0x64 - 0xeb) * t)
        b = int(0xd0 + (0x0c - 0xd0) * t)
        colors.append(f"#{r:02x}{g:02x}{b:02x}")
    ax1.bar(xs, org_frac, color=colors, edgecolor="white")
    ax1.set_ylim(0, 1.05)
    ax1.set_xlabel("order (k)")
    ax1.set_ylabel("fraction of SOs that are full organizations")
    ax1.set_title("Organization rate")
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    fig.tight_layout()
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    return out_path
