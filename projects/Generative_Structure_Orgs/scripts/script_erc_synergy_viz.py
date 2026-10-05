#!/usr/bin/env python3
"""
script_erc_synergy_viz.py
=========================
Plot the ERC hierarchy of a single reaction network and overlay synergy
arrows between the reactant ERCs and the synergy target ERC.

Synergy rendering:
  Each synergy (E1, E2) → T is drawn as a small "junction" diamond node
  positioned at the centroid of E1, E2 and T.  Three arrows are drawn:
    E1  → junction   (two converging lines)
    E2  → junction
    junction → T     (single outgoing line)
  This makes the merging topology visually explicit.

Junction and arrow colours by synergy classification (each synergy drawn
once at its highest level):
  pale blue  (#A8D8EA)  — basic synergy only
  orange     (#E67E22)  — maximal synergy  (not fundamental)
  green      (#27AE60)  — fundamental synergy

Node colour:
  green  — self-maintaining  (SM)
  orange — semi-self-maintaining  (SSM only)
  steel  — neither / not reactive

Node size ∝ closure size.

Correctness note
----------------
The library's get_maximal_synergies() has an inverted containment check.
Corrected versions are implemented here as _compute_basic() /
_filter_maximal_v2() / _filter_fundamental().
"""

import os
import sys
import time
from itertools import combinations
from collections import defaultdict

import numpy as np
from scipy.optimize import linprog
from matplotlib.patches import FancyArrowPatch
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import networkx as nx

# -- Path setup ----------------------------------------------------------------
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_PYCOT_ROOT = os.path.normpath(os.path.join(_SCRIPT_DIR, '..', '..', '..'))
sys.path.insert(0, os.path.join(_PYCOT_ROOT, 'src'))

from pyCOT.io.functions import read_txt
from pyCOT.analysis.ERC_Hierarchy import ERC, ERC_Hierarchy, species_list_to_names
from pyCOT.analysis.SORN_Generators import is_semi_self_maintaining

# -- Network file --------------------------------------------------------------
RN_FILE = os.path.join(_PYCOT_ROOT, 'data', 'biomodels',
                       'biomodels_interesting', 'bigg_iAF692.txt')
# Uncomment alternatives:
# RN_FILE = os.path.join(_PYCOT_ROOT, 'networks', 'testing', 'Farm.txt')
# RN_FILE = 'data\\Examples_tests\\testing\\ERC_synergy0.txt'
RN_FILE = 'data/Examples_tests/Fundamental_Synergy_Example.txt'

# -- Visual parameters ---------------------------------------------------------
# Target print size: ~7 lines of body text tall x 40% of the page width
# (roughly 1.2in x 2.6in at typical 10-11pt article settings) -- everything
# below is sized so it survives that shrink, not to look right at 100%
# on screen. Bigger source fonts/nodes/lines = more pixels-per-element
# once LaTeX scales the raster down to the small print footprint.
NODE_SIZE_BASE  = 900    # was 400
NODE_SIZE_SCALE = 220    # was 120
NODE_FONT_SIZE  = 20     # ERC label font (was 12)
COL_SM   = '#27AE60'
COL_SSM  = '#E67E22'
COL_DEF  = '#AED6F1'
SM_EPS   = 1e-6

# Synergy colours (one per classification level)
COL_SYN = {

    #'basic' is normal blue
    'basic':       "#1C5FB7",   # blue
    'maximal':     "#D76E0B",   # orange
    'fundamental': '#27AE60',   # green
}

# Junction diamond marker size (scatter units = pts²)
JUNC_SIZE = 180   # was 80

# Arrow line widths and opacities per level -- boosted across the board;
# containment edges (drawn separately below) get an even bigger boost
# since those were flagged as the least visible at print size.
ARROW_LW    = {'basic': 1.6, 'maximal': 2.0, 'fundamental': 2.4}
ARROW_ALPHA = {'basic': 0.80, 'maximal': 0.88, 'fundamental': 0.95}

# Containment-edge (ERC hierarchy) styling -- the "especially containment
# ones" arrows called out for visibility.
CONTAIN_EDGE_COLOR    = '#333333'   # was '#555555' -- darker, more contrast
CONTAIN_EDGE_WIDTH    = 2.6         # was 1.1
CONTAIN_EDGE_ALPHA    = 0.9         # was 0.65
CONTAIN_ARROWSIZE     = 30          # was 14

# Shrink (in display pts) to pull arrow endpoints back from node centres.
# Formula: ~sqrt(scatter_size / π) gives the marker radius in pts.
# Rescaled to match the larger NODE_SIZE_BASE/JUNC_SIZE above.
SHRINK_ERC  = 20   # was 12
SHRINK_JUNC =  8   # was 5

# Set False to hide non-maximal basic synergies (less visual clutter)
SHOW_BASIC = True

# Set False to hide any ERC whose closure is the empty set (e.g. E0 in a
# system where the empty set is itself closed) -- off by default, since
# that node carries no species and is usually not worth the space.
SHOW_EMPTY_ERC = False

# Vertical spacing between hierarchy levels, in data units (was a flat 2.0
# in the position formula below). Compressed so the whole figure comes out
# wide-and-short -- matching a 40%-width / 7-line-tall print target --
# instead of the taller, more square layout the original 2.0 spacing gives.
LEVEL_Y_SPACING = 1.2
NODE_X_SPACING  = 2.4   # widened slightly to compensate for bigger node size

# How far to push junctions off an ERC level line when the computed y lands
# on one, and the same-level junction offset just below -- both were tuned
# as fractions of the old flat 2.0 level spacing, rescaled to match
# LEVEL_Y_SPACING above.
PUSH_OFF_LEVEL = 0.40 * (LEVEL_Y_SPACING / 2.0)


# ===========================================================================
# Node classification helpers
# ===========================================================================

def _is_reactive(RN, species_set):
    if not species_set:
        return False
    return any(r.support_indices() for r in RN.get_reactions_from_species(species_set))


def _is_sm(RN, species_set):
    sub = RN.sub_reaction_network(species_set)
    S   = np.asarray(sub.stoichiometry_matrix(), dtype=float)
    _, n_rx = S.shape
    if n_rx == 0:
        return False
    res = linprog(np.zeros(n_rx),
                  A_ub=-S, b_ub=SM_EPS * S.sum(axis=1),
                  bounds=[(0, None)] * n_rx, method='highs')
    return res.status == 0


# ===========================================================================
# Corrected synergy detection
# ===========================================================================

def _compute_basic(erc1, erc2, hierarchy, RN):
    """
    Basic synergies (Def. 18): joint closure covers a minimal generator of T
    that neither erc1 nor erc2 alone covers.
    """
    if (erc1 in hierarchy.get_contain(erc2) or
            erc2 in hierarchy.get_contain(erc1)):
        return []

    cl1   = erc1.get_closure_names(RN)
    cl2   = erc2.get_closure_names(RN)
    joint = cl1 | cl2
    sub1  = {e.label for e in hierarchy.get_contain(erc1)}
    sub2  = {e.label for e in hierarchy.get_contain(erc2)}

    result = []
    for target in hierarchy.ercs:
        if target is erc1 or target is erc2:
            continue
        if target.label in sub1 or target.label in sub2:
            continue
        for gen in target.min_generators:
            gen_sp = set(species_list_to_names(gen))
            if gen_sp.issubset(joint) and not gen_sp.issubset(cl1) and not gen_sp.issubset(cl2):
                result.append((erc1, erc2, target))
                break
    return result


def _filter_maximal_v2(basics, RN):
    """
    Definition 19: (E1,E2)→T is maximal iff no other (E1,E2)→T' exists
    with cl(T) ⊊ cl(T') (T' strictly bigger).
    Comparison done directly on cached closure name sets.
    """
    if not basics:
        return []
    by_pair = defaultdict(list)
    for e1, e2, t in basics:
        by_pair[tuple(sorted([e1.label, e2.label]))].append((e1, e2, t))

    result = []
    for syns in by_pair.values():
        cls = {s[2].label: s[2].get_closure_names(RN) for s in syns}
        for e1, e2, T in syns:
            if not any(cls[T.label] < cls[T2.label]
                       for _, _, T2 in syns if T2.label != T.label):
                result.append((e1, e2, T))
    return result


def _filter_fundamental(maximals, hierarchy, RN):
    """
    Definition 20: (E1,E2)→T is fundamental iff no strictly smaller pair
    (E1'⊆E1, E2'⊆E2, at least one strict) has a maximal synergy to T.
    """
    if not maximals:
        return []

    def _desc(erc):
        if erc.label not in hierarchy.graph:
            return {erc.label}
        return {erc.label} | set(nx.descendants(hierarchy.graph, erc.label))

    desc  = {}
    l2erc = {e.label: e for e in hierarchy.ercs}
    for e1, e2, _ in maximals:
        for e in (e1, e2):
            if e.label not in desc:
                desc[e.label] = _desc(e)

    result = []
    for e1, e2, target in maximals:
        orig = tuple(sorted([e1.label, e2.label]))
        is_fund = True
        for l1 in desc[e1.label]:
            if not is_fund:
                break
            for l2 in desc[e2.label]:
                if l1 == l2 or tuple(sorted([l1, l2])) == orig:
                    continue
                s1, s2 = l2erc.get(l1), l2erc.get(l2)
                if s1 is None or s2 is None:
                    continue
                sub_b = _compute_basic(s1, s2, hierarchy, RN)
                sub_m = _filter_maximal_v2(sub_b, RN)
                if any(s[2].label == target.label for s in sub_m):
                    is_fund = False
                    break
        if is_fund:
            result.append((e1, e2, target))
    return result


def compute_all_synergies(ercs, hierarchy, RN):
    """Returns (all_basics, all_maximals, all_fundamentals) as lists of (e1,e2,target)."""
    all_b, all_m, all_f = [], [], []
    pairs = list(combinations(ercs, 2))
    print(f"  Checking {len(pairs)} ERC pairs...")
    for e1, e2 in pairs:
        b = _compute_basic(e1, e2, hierarchy, RN)
        if not b:
            continue
        all_b.extend(b)
        m = _filter_maximal_v2(b, RN)
        all_m.extend(m)
        f = _filter_fundamental(m, hierarchy, RN)
        all_f.extend(f)
    return all_b, all_m, all_f


# ===========================================================================
# Load, compute ERCs, build hierarchy
# ===========================================================================
print(f"Loading: {RN_FILE}")
RN = read_txt(RN_FILE)
print(f"  {len(RN.species())} species, {len(RN.reactions())} reactions")

print("Computing ERCs...")
t0   = time.time()
ercs = ERC.ERCs(RN)
print(f"  {len(ercs)} ERCs  ({time.time()-t0:.1f}s)")
for erc in ercs:
    print(f"  {erc.label}: closure={species_list_to_names(erc.get_closure(RN))}")

print("Building hierarchy...")
hierarchy = ERC_Hierarchy(RN, ercs)
G = hierarchy.graph.copy()   # copy: safe to mutate without touching hierarchy's own graph

if not SHOW_EMPTY_ERC:
    empty_labels = [erc.label for erc in ercs if len(erc.get_closure(RN)) == 0]
    if empty_labels:
        print(f"  SHOW_EMPTY_ERC=False: hiding empty-closure ERC(s) {empty_labels}")
        G.remove_nodes_from(empty_labels)
        # Downstream synergy-drawing code already skips any junction/edge
        # whose endpoint isn't in `pos` (built from G's nodes below), so
        # removing the node from G here is enough to also drop every
        # containment edge and synergy arrow that touched it.

# Classify nodes
node_sizes, node_colors = {}, {}
for erc in ercs:
    cl   = erc.get_closure(RN)
    n_sp = len(cl)
    node_sizes[erc.label] = NODE_SIZE_BASE + NODE_SIZE_SCALE * n_sp
    if not _is_reactive(RN, cl):
        node_colors[erc.label] = COL_DEF
    elif not is_semi_self_maintaining(RN, cl):
        node_colors[erc.label] = COL_DEF
    elif _is_sm(RN, cl):
        node_colors[erc.label] = COL_SM
    else:
        node_colors[erc.label] = COL_SSM

# Compute synergies
print("Computing synergies...")
t0 = time.time()
all_b, all_m, all_f = compute_all_synergies(ercs, hierarchy, RN)
print(f"  Basic={len(all_b)}  Maximal={len(all_m)}  Fundamental={len(all_f)}"
      f"  ({time.time()-t0:.1f}s)")

# Assign each synergy its highest classification level
# key: (sorted_pair_tuple, target_label)  →  'basic' | 'maximal' | 'fundamental'
syn_level = {}
for e1, e2, t in all_b:
    k = (tuple(sorted([e1.label, e2.label])), t.label)
    syn_level[k] = 'basic'
for e1, e2, t in all_m:
    k = (tuple(sorted([e1.label, e2.label])), t.label)
    syn_level[k] = 'maximal'
for e1, e2, t in all_f:
    k = (tuple(sorted([e1.label, e2.label])), t.label)
    syn_level[k] = 'fundamental'

# Recover (e1, e2) label order from the original lists for junction placement
# We need the two individual reactant labels (not just sorted pair).
syn_reactants = {}   # key → (e1_label, e2_label)
for e1, e2, t in all_b:
    k = (tuple(sorted([e1.label, e2.label])), t.label)
    syn_reactants[k] = (e1.label, e2.label)
for e1, e2, t in all_m:
    k = (tuple(sorted([e1.label, e2.label])), t.label)
    syn_reactants[k] = (e1.label, e2.label)
for e1, e2, t in all_f:
    k = (tuple(sorted([e1.label, e2.label])), t.label)
    syn_reactants[k] = (e1.label, e2.label)


# ===========================================================================
# Layout: level-based (same as script_erc_hierarchy.py)
# ===========================================================================
levels      = ERC.get_node_levels(G)
level_nodes = defaultdict(list)
for node, lvl in levels.items():
    level_nodes[lvl].append(node)

pos = {}
for lvl, nodes in level_nodes.items():
    nodes.sort(key=lambda n: len(nx.ancestors(G, n)), reverse=True)
    for i, node in enumerate(nodes):
        pos[node] = np.array([(i - (len(nodes) - 1) / 2) * NODE_X_SPACING, lvl * LEVEL_Y_SPACING])

# Junction positions: between the reactant midpoint and the target,
# always kept between levels (never on top of an ERC node).
#
# X: midpoint of E1 and E2 — keeps the junction visually between the reactants.
# Y: midpoint between max(y_E1, y_E2) and y_T — keeps the junction between
#    the highest reactant level and the target level.
# Push-off: if the computed y falls on an ERC level line, shift toward target.
# Jitter: when the same reactant pair has multiple synergies, spread junctions
#         horizontally so they don't stack.

level_ys = sorted(set(float(p[1]) for p in pos.values()))

# Proximity/offset thresholds below scale with LEVEL_Y_SPACING (all were
# originally tuned as fractions of the old flat 2.0 spacing) so junction
# placement stays proportionally correct now that levels are compressed.
_LEVEL_SCALE = LEVEL_Y_SPACING / 2.0
_LEVEL_CLEARANCE = 0.22 * _LEVEL_SCALE

def _clear_of_levels(y, level_ys, y_target, step=PUSH_OFF_LEVEL, max_iter=8):
    """Shift y toward y_target until it is not within _LEVEL_CLEARANCE of any ERC level."""
    direction = 1 if y_target >= y else -1
    for _ in range(max_iter):
        if not any(abs(y - ly) < _LEVEL_CLEARANCE for ly in level_ys):
            break
        y += direction * step
    return y

junc_pos = {}
pair_junc_count = defaultdict(int)   # sorted-pair → # junctions already placed

for key, (l1, l2) in syn_reactants.items():
    _, t_label = key
    if l1 not in pos or l2 not in pos or t_label not in pos:
        continue

    p1 = pos[l1]
    p2 = pos[l2]
    pt = pos[t_label]

    # X: centred between the two reactants (target x does not influence this)
    x_junc = (p1[0] + p2[0]) / 2.0

    # Y: midpoint between the higher reactant and the target
    y_hi  = max(p1[1], p2[1])
    y_T   = float(pt[1])
    if abs(y_T - y_hi) > 0.5 * _LEVEL_SCALE:
        y_junc = (y_hi + y_T) / 2.0
    else:
        # Reactants and target at the same level — push junction in target direction
        y_sign = 1 if y_T >= y_hi else -1
        y_junc = y_hi + y_sign * 0.7 * _LEVEL_SCALE

    # Ensure junction is not sitting on an ERC level line
    y_junc = _clear_of_levels(y_junc, level_ys, y_T)

    base = np.array([x_junc, y_junc])

    # Horizontal jitter so multiple synergies for the same reactant pair
    # don't pile up at the same x position
    pair_key = tuple(sorted([l1, l2]))
    count = pair_junc_count[pair_key]
    if count > 0:
        sign   = 1 if count % 2 == 1 else -1
        offset = ((count + 1) // 2) * 0.42 * (NODE_X_SPACING / 2.0)
        base   = base + np.array([sign * offset, 0.0])
    pair_junc_count[pair_key] += 1

    junc_pos[key] = base


# ===========================================================================
# Plot
# ===========================================================================
ordered_nodes = list(G.nodes())
sizes_list    = [node_sizes.get(n, NODE_SIZE_BASE) for n in ordered_nodes]
colors_list   = [node_colors.get(n, COL_DEF)        for n in ordered_nodes]

# Wide-and-short canvas to match the 40%-page-width / 7-line-tall print
# target (bbox_inches='tight' below will crop to actual content, but this
# sets the aspect ratio spacing/legends are laid out against).
fig, ax = plt.subplots(figsize=(13, 5.5))

# -- ERC hierarchy (nodes + containment edges) --------------------------------
scatter = nx.draw_networkx_nodes(G, pos, ax=ax,
                                 nodelist=ordered_nodes,
                                 node_color=colors_list,
                                 node_size=sizes_list,
                                 alpha=0.92)
scatter.set_zorder(3)
scatter.set_edgecolor('#222222')
scatter.set_linewidth(1.4)
nx.draw_networkx_edges(G, pos, ax=ax,
                       edge_color=CONTAIN_EDGE_COLOR, arrows=True,
                       arrowsize=CONTAIN_ARROWSIZE,
                       width=CONTAIN_EDGE_WIDTH, alpha=CONTAIN_EDGE_ALPHA,
                       node_size=sizes_list)
nx.draw_networkx_labels(G, pos, ax=ax,
                        font_size=NODE_FONT_SIZE, font_weight='bold')

# -- Helper: draw a FancyArrowPatch -------------------------------------------
def _farrow(src, dst, col, lw, alpha, shrinkA=SHRINK_ERC, shrinkB=SHRINK_ERC):
    """
    Draw an arrow from src to dst (both numpy arrays in data coordinates).
    shrinkA / shrinkB are in display pts; they pull the arrow endpoints back
    from the node centres so the arrowhead lands at the visual boundary.
    """
    patch = FancyArrowPatch(
        posA=tuple(src), posB=tuple(dst),
        arrowstyle='->', lw=lw, color=col, alpha=alpha,
        shrinkA=shrinkA, shrinkB=shrinkB,
        mutation_scale=22,    # controls arrowhead size (was 10)
        zorder=1,
    )
    ax.add_patch(patch)

# -- Draw synergies -----------------------------------------------------------
for key, stype in syn_level.items():
    if stype == 'basic' and not SHOW_BASIC:
        continue
    if key not in junc_pos:
        continue
    l1, l2 = syn_reactants[key]
    _, t_label = key
    if l1 not in pos or l2 not in pos or t_label not in pos:
        continue

    col   = COL_SYN[stype]
    lw    = ARROW_LW[stype]
    alpha = ARROW_ALPHA[stype]
    junc  = junc_pos[key]

    # Draw junction diamond marker
    ax.scatter(*junc, s=JUNC_SIZE, c=col, marker='D',
               edgecolors='white', linewidths=0.8, alpha=0.95, zorder=3)

    # E1 → junction  (shrinkB is small: junction node is tiny)
    _farrow(pos[l1], junc, col, lw, alpha, shrinkA=SHRINK_ERC, shrinkB=SHRINK_JUNC)
    # E2 → junction
    _farrow(pos[l2], junc, col, lw, alpha, shrinkA=SHRINK_ERC, shrinkB=SHRINK_JUNC)
    # junction → T  (shrinkA small, shrinkB matches ERC node size)
    _farrow(junc, pos[t_label], col, lw, alpha, shrinkA=SHRINK_JUNC, shrinkB=SHRINK_ERC)

# -- Hover tooltip (ERC nodes) ------------------------------------------------
pos_display = {n: p for n, p in pos.items()}   # reuse

annot = ax.annotate("", xy=(0, 0), xytext=(15, 15),
                    textcoords="offset points",
                    bbox=dict(boxstyle="round", fc="white", ec="0.5", alpha=0.9),
                    visible=False, zorder=10)
erc_by_label = {erc.label: erc for erc in ercs}

def _on_hover(event):
    if event.inaxes != ax:
        return
    cont, ind = scatter.contains(event)
    if cont:
        node = ordered_nodes[ind["ind"][0]]
        erc  = erc_by_label.get(node)
        cl   = species_list_to_names(erc.get_closure(RN)) if erc else []
        annot.xy = tuple(pos[node])
        annot.set_text(f"{node}\n{cl}")
        annot.set_visible(True)
    else:
        annot.set_visible(False)
    fig.canvas.draw_idle()

fig.canvas.mpl_connect("motion_notify_event", _on_hover)

# Pad the data limits a bit so the two corner-anchored legends that stay
# INSIDE the axes (closure size, synergy type) have clear space regardless
# of the specific network's hierarchy shape.
ax.margins(x=0.15, y=0.15)

# -- Legends ------------------------------------------------------------------
# 1) Maintenance class (node colour). Anchored OUTSIDE the axes (to the
# right, bottom-aligned) rather than an inside corner like the other two:
# an inside 'lower right' collided with whichever node ends up in that
# corner, which corner that is depends on the network's own hierarchy
# shape (and shifts further when SHOW_EMPTY_ERC changes which nodes
# exist) -- not something a fixed inside-corner placement can be made
# robust to, confirmed by hitting the collision twice while tuning
# margins. Outside the axes there is nothing to collide with.
# Uses fig.legend() rather than ax.legend()+add_artist(): confirmed by a
# standalone test that ax.legend()+add_artist() with an outside-axes
# bbox_to_anchor silently drops the legend from bbox_inches='tight''s
# saved extent (renders fine on screen, vanishes from the saved PNG) --
# fig.legend() does not have that problem since it isn't competing with
# the axes' own single "current legend" slot that add_artist works around.
maintenance_legend = fig.legend(
    handles=[
        mpatches.Patch(facecolor=COL_SM,  label='Self-maintaining  (SM)'),
        mpatches.Patch(facecolor=COL_SSM, label='Semi-self-maintaining  (SSM)'),
        mpatches.Patch(facecolor=COL_DEF, label='Neither / not reactive'),
    ],
    loc='lower left', bbox_to_anchor=(0.92, 0.02), bbox_transform=fig.transFigure,
    fontsize=12, framealpha=0.9,
    title='Maintenance class', title_fontsize=12,
)

# 2) Synergy type (arrow + junction colour)
syn_handles = [
    Line2D([0], [0], color=COL_SYN['fundamental'], lw=2.5,
           label=f'Fundamental synergy  ({len(all_f)})'),
    Line2D([0], [0], color=COL_SYN['maximal'],     lw=1.8,
           label=f'Maximal synergy  ({len(all_m)})'),
]
if SHOW_BASIC:
    syn_handles.append(
        Line2D([0], [0], color=COL_SYN['basic'], lw=1.0,
               label=f'Basic synergy  ({len(all_b)})'))
syn_legend = ax.legend(handles=syn_handles, loc='upper right', fontsize=12,
                       framealpha=0.9, title='Synergy type', title_fontsize=12)
ax.add_artist(syn_legend)

# 3) Closure size (node size)
all_counts = sorted({len(erc.get_closure(RN)) for erc in ercs if erc.label in G})
if len(all_counts) > 4:
    tick_idx   = [0, len(all_counts)//3, 2*len(all_counts)//3, -1]
    size_ticks = sorted({all_counts[i] for i in tick_idx})
else:
    size_ticks = all_counts
size_handles = [
    ax.scatter([], [], s=NODE_SIZE_BASE + NODE_SIZE_SCALE * n,
               color='#888888', alpha=0.85, label=f'{n} species')
    for n in size_ticks
]
ax.legend(handles=size_handles, loc='upper left', fontsize=12, framealpha=0.9,
          title='Closure size', title_fontsize=12,
          labelspacing=1.2, handletextpad=1.0)

# -- Title --------------------------------------------------------------------
net_name = os.path.splitext(os.path.basename(RN_FILE))[0]
ax.set_title(
    f"ERC Hierarchy + Synergies — {net_name}\n"
    f"Basic: {len(all_b)}   Maximal: {len(all_m)}   Fundamental: {len(all_f)}   "
    f"◆ = synergy junction",
    fontsize=12)
ax.axis('off')
plt.tight_layout()
_out_dir = os.path.join(_SCRIPT_DIR, '..', 'outputs', 'synergy_viz')
os.makedirs(_out_dir, exist_ok=True)
_out_path = os.path.join(_out_dir, f'synergy_viz_{net_name}.png')
plt.savefig(_out_path, dpi=150, bbox_inches='tight')
print(f"Saved: {os.path.abspath(_out_path)}")
plt.show()
