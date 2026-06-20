#!/usr/bin/env python3
"""
script_erc_complementarity_viz.py
===================================
Plot the ERC hierarchy of a single reaction network and overlay
complementarity relationships as intermediate junction nodes.

Complementarity rendering
--------------------------
Each directional supply  E_producer ——[s]——> E_consumer  is drawn as a
small square junction node placed at the midpoint between the two ERCs.
Two arrows are drawn:
    E_producer  →  junction   (producer side)
    junction    →  E_consumer (consumer side)
The junction is labelled with the supplied species name(s).

Junction and arrow colours by classification (each direction drawn once
at its highest applicable level):
    blue   (#3498DB) — complementary (also synergetic, lowest level)
    orange (#E67E22) — purely complementary (not synergetic)
    green  (#27AE60) — fundamental complementarity (E_prod ∈ minprod(s),
                        E_cons ∈ mincons(s))

Mutual complementarity (E1 supplies to E2 AND E2 supplies to E1) produces
two separate junctions; a small perpendicular offset separates them.

Node colour (maintenance class):
    green  — self-maintaining (SM)
    orange — semi-self-maintaining only (SSM)
    steel  — neither / not reactive

Node size ∝ closure size.

Correctness note
-----------------
The supply-based complementarity definition is implemented directly from
Section 5 of the paper:
    supl(E, E') = prod(R_E) ∩ req(E')   with  req(E) = supp(R_E) \\ prod(R_E)
Synergy detection uses the corrected _compute_basic() function from
script_erc_synergy_viz.py (library's get_maximal_synergies is inverted).
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
# Alternatives:
# RN_FILE = os.path.join(_PYCOT_ROOT, 'networks', 'testing', 'Farm.txt')
# RN_FILE = 'data\\Examples_tests\\testing\\ERC_synergy0.txt'

# -- Visual parameters ---------------------------------------------------------
NODE_SIZE_BASE  = 400
NODE_SIZE_SCALE = 120
COL_SM   = '#27AE60'
COL_SSM  = '#E67E22'
COL_DEF  = '#AED6F1'
SM_EPS   = 1e-6

# Complementarity level colours
COL_COMP = {
    'complementary': '#3498DB',   # blue — complementary (also synergetic)
    'pure':          '#E67E22',   # orange — purely complementary
    'fundamental':   '#27AE60',   # green — fundamental
}

# Junction square marker size (scatter units = pts²)
JUNC_SIZE = 90

# Arrow widths and opacities
ARROW_LW    = {'complementary': 0.7, 'pure': 0.9, 'fundamental': 1.1}
ARROW_ALPHA = {'complementary': 0.65, 'pure': 0.75, 'fundamental': 0.90}

# Shrink from node centres to arrowhead (display pts)
SHRINK_ERC  = 12
SHRINK_JUNC = 5

# How far to push junctions off an ERC level line when the midpoint lands on one.
# ERC levels sit at y = 0, 2, 4, ... so 0.8 places the junction safely between levels.
PUSH_OFF_LEVEL = 0.8

# Show only fundamental junctions by default to reduce clutter.
# Set to False to also draw non-fundamental complementary junctions.
SHOW_FUNDAMENTAL_ONLY = False


# =============================================================================
# Node classification helpers (identical to synergy_viz)
# =============================================================================

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
    rhs    = SM_EPS * S.sum(axis=1)
    res    = linprog(np.zeros(n_rx),
                     A_ub=-S, b_ub=rhs,
                     bounds=[(0, None)] * n_rx, method='highs')
    return res.status == 0


# =============================================================================
# Synergy helper (for "purely complementary" classification)
# =============================================================================

def _has_basic_synergy(erc1, erc2, hierarchy, RN):
    """Return True iff (erc1,erc2) has at least one basic synergy."""
    if (erc1 in hierarchy.get_contain(erc2) or
            erc2 in hierarchy.get_contain(erc1)):
        return False
    cl1   = erc1.get_closure_names(RN)
    cl2   = erc2.get_closure_names(RN)
    joint = cl1 | cl2
    sub1  = {e.label for e in hierarchy.get_contain(erc1)}
    sub2  = {e.label for e in hierarchy.get_contain(erc2)}
    for target in hierarchy.ercs:
        if target is erc1 or target is erc2:
            continue
        if target.label in sub1 or target.label in sub2:
            continue
        for gen in target.min_generators:
            gen_sp = set(species_list_to_names(gen))
            if (gen_sp.issubset(joint) and
                    not gen_sp.issubset(cl1) and
                    not gen_sp.issubset(cl2)):
                return True
    return False


# =============================================================================
# Supply and complementarity computation
# =============================================================================

def _is_incomparable(erc1, erc2, hierarchy):
    return (erc1 not in hierarchy.get_contain(erc2) and
            erc2 not in hierarchy.get_contain(erc1))


def _supply(erc_prod, erc_cons, RN):
    """supl(erc_prod, erc_cons) = prod(R_{erc_prod}) ∩ req(erc_cons)."""
    return erc_prod.get_produced_species(RN) & erc_cons.get_required_species(RN)


def _compute_minprod_mincons(ercs, hierarchy, RN):
    """
    Return (minprod, mincons) dicts mapping species_name → list of ERCs.
    minprod[s]: inclusion-minimal ERCs producing s.
    mincons[s]: inclusion-minimal ERCs requiring s.
    """
    producers = defaultdict(set)
    consumers = defaultdict(set)

    for erc in ercs:
        for s in erc.get_produced_species(RN):
            producers[s].add(erc.label)
        for s in erc.get_required_species(RN):
            consumers[s].add(erc.label)

    label_to_erc = {e.label: e for e in ercs}

    def _minimal(label_set):
        minimal = []
        for lbl in label_set:
            erc = label_to_erc[lbl]
            descendants = {e.label for e in hierarchy.get_contain(erc)}
            if not descendants.intersection(label_set):
                minimal.append(erc)
        return minimal

    minprod = {s: _minimal(lbls) for s, lbls in producers.items()}
    mincons = {s: _minimal(lbls) for s, lbls in consumers.items()}
    return minprod, mincons


def _classify_supply(e_prod, e_cons, supply_species, pair_is_syn, minprod, mincons):
    """
    Return the highest classification level for a directional supply
    (e_prod → e_cons via species in supply_species).

    Levels (highest to lowest): 'fundamental' > 'pure' > 'complementary'
    """
    # Fundamental: any supply species has e_prod ∈ minprod(s) AND e_cons ∈ mincons(s)
    for s in supply_species:
        if e_prod in minprod.get(s, []) and e_cons in mincons.get(s, []):
            return 'fundamental'
    # Purely complementary: not synergetic
    if not pair_is_syn:
        return 'pure'
    return 'complementary'


def compute_complementarity_entries(ercs, hierarchy, RN):
    """
    Return a list of dicts, one per directional supply (producer → consumer):
      {
        'producer':  ERC object,
        'consumer':  ERC object,
        'species':   set of supply species names,
        'level':     'fundamental' | 'pure' | 'complementary',
      }
    """
    minprod, mincons = _compute_minprod_mincons(ercs, hierarchy, RN)
    entries = []

    for e1, e2 in combinations(ercs, 2):
        if not _is_incomparable(e1, e2, hierarchy):
            continue

        s12 = _supply(e1, e2, RN)
        s21 = _supply(e2, e1, RN)

        if not s12 and not s21:
            continue

        pair_syn = _has_basic_synergy(e1, e2, hierarchy, RN)

        if s12:
            level = _classify_supply(e1, e2, s12, pair_syn, minprod, mincons)
            entries.append({'producer': e1, 'consumer': e2,
                            'species': s12, 'level': level})
        if s21:
            level = _classify_supply(e2, e1, s21, pair_syn, minprod, mincons)
            entries.append({'producer': e2, 'consumer': e1,
                            'species': s21, 'level': level})

    return entries


# =============================================================================
# Load, compute ERCs and hierarchy
# =============================================================================
print(f"Loading: {RN_FILE}")
RN = read_txt(RN_FILE)
print(f"  {len(RN.species())} species, {len(RN.reactions())} reactions")

print("Computing ERCs...")
t0   = time.time()
ercs = ERC.ERCs(RN)
print(f"  {len(ercs)} ERCs  ({time.time()-t0:.1f}s)")

print("Building hierarchy...")
hierarchy = ERC_Hierarchy(RN, ercs)
G = hierarchy.graph

# Classify nodes
node_sizes  = {}
node_colors = {}
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

# Compute complementarity relationships
print("Computing complementarities...")
t0 = time.time()
comp_entries = compute_complementarity_entries(ercs, hierarchy, RN)
t_comp = time.time() - t0

# Partition by level for statistics
n_fund = sum(1 for e in comp_entries if e['level'] == 'fundamental')
n_pure = sum(1 for e in comp_entries if e['level'] == 'pure')
n_comp = sum(1 for e in comp_entries if e['level'] == 'complementary')
print(f"  Fundamental={n_fund}  Purely complementary={n_pure}  "
      f"Complementary(+synergy)={n_comp}  ({t_comp:.1f}s)")

# =============================================================================
# Layout: level-based (same as synergy_viz and hierarchy scripts)
# =============================================================================
levels      = ERC.get_node_levels(G)
level_nodes = defaultdict(list)
for node, lvl in levels.items():
    level_nodes[lvl].append(node)

pos = {}
for lvl, nodes in level_nodes.items():
    nodes.sort(key=lambda n: len(nx.ancestors(G, n)), reverse=True)
    for i, node in enumerate(nodes):
        pos[node] = np.array([(i - (len(nodes) - 1) / 2) * 2.0, lvl * 2.0])

# =============================================================================
# Junction positions
#
# Each junction sits at the true midpoint between producer and consumer.
# ERC nodes live on discrete level lines:  y = 0, 2, 4, 6, ...
# A midpoint can land on a level line in two cases:
#   • same-level pairs  (y_prod == y_cons → midpoint at the same y)
#   • pairs two or more levels apart whose midpoint hits an intermediate level
# In both cases we push the junction perpendicular to the producer→consumer
# direction by PUSH_OFF_LEVEL units so it sits visually between two levels.
# For mutual pairs (E1→E2 and E2→E1 both present) an additional perpendicular
# jitter separates the two junctions.
# =============================================================================
level_ys   = sorted(set(float(p[1]) for p in pos.values()))
junc_pos   = {}
pair_count = defaultdict(int)

for idx, entry in enumerate(comp_entries):
    lp = entry['producer'].label
    lc = entry['consumer'].label
    if lp not in pos or lc not in pos:
        continue

    pp = pos[lp]
    pc = pos[lc]

    # True midpoint – places the junction halfway between the two ERC nodes
    mid = (pp + pc) * 0.5

    # Perpendicular unit vector (90° counter-clockwise rotation of prod→cons)
    direction = pc - pp
    dist = np.linalg.norm(direction)
    if dist > 1e-9:
        perp = np.array([-direction[1], direction[0]]) / dist
    else:
        perp = np.array([0.0, 1.0])

    # If the midpoint y coincides with an ERC level line, push it off so the
    # junction lands between two levels rather than on top of ERC nodes.
    mid_y   = float(mid[1])
    on_level = any(abs(mid_y - ly) < 0.25 for ly in level_ys)
    if on_level:
        # Alternate the push direction per pair so crowded regions spread out
        sign = 1 if hash(frozenset([lp, lc])) % 2 == 0 else -1
        mid  = mid + sign * PUSH_OFF_LEVEL * perp

    # Additional perpendicular jitter for mutual pairs (bidirectional supply)
    pair_key = frozenset([lp, lc])
    count    = pair_count[pair_key]
    if count > 0:
        sign = 1 if count % 2 == 1 else -1
        mid  = mid + sign * (count + 1) * 0.32 * perp
    pair_count[pair_key] += 1
    junc_pos[idx] = mid


# =============================================================================
# Plot
# =============================================================================
ordered_nodes = list(G.nodes())
sizes_list    = [node_sizes.get(n, NODE_SIZE_BASE) for n in ordered_nodes]
colors_list   = [node_colors.get(n, COL_DEF)       for n in ordered_nodes]

fig, ax = plt.subplots(figsize=(14, 9))

# -- ERC hierarchy nodes + containment edges ----------------------------------
scatter = nx.draw_networkx_nodes(G, pos, ax=ax,
                                 nodelist=ordered_nodes,
                                 node_color=colors_list,
                                 node_size=sizes_list,
                                 alpha=0.92)
scatter.set_zorder(3)
nx.draw_networkx_edges(G, pos, ax=ax,
                       edge_color='#555555', arrows=True, arrowsize=14,
                       width=1.1, alpha=0.65)
nx.draw_networkx_labels(G, pos, ax=ax,
                        font_size=12, font_weight='bold')


# -- FancyArrowPatch helper ---------------------------------------------------
def _farrow(src, dst, col, lw, alpha, style='->', shrinkA=SHRINK_ERC, shrinkB=SHRINK_ERC):
    patch = FancyArrowPatch(
        posA=tuple(src), posB=tuple(dst),
        arrowstyle=style,
        lw=lw, color=col, alpha=alpha,
        shrinkA=shrinkA, shrinkB=shrinkB,
        mutation_scale=10,
        zorder=1,
    )
    ax.add_patch(patch)


# -- Draw complementarity junctions and arrows --------------------------------
level_order = {'fundamental': 3, 'pure': 2, 'complementary': 1}

for idx, entry in enumerate(comp_entries):
    level = entry['level']

    if SHOW_FUNDAMENTAL_ONLY and level != 'fundamental':
        continue

    if idx not in junc_pos:
        continue

    lp = entry['producer'].label
    lc = entry['consumer'].label
    if lp not in pos or lc not in pos:
        continue

    col   = COL_COMP[level]
    lw    = ARROW_LW[level]
    alpha = ARROW_ALPHA[level]
    jpos  = junc_pos[idx]

    # For non-fundamental levels use dashed arrow style
    arrow_style = '->' if level == 'fundamental' else '->'

    # Junction square marker
    ax.scatter(*jpos, s=JUNC_SIZE,
               c=col, marker='s',
               edgecolors='white', linewidths=0.8,
               alpha=0.95, zorder=4)

    # Label junction with supply species (abbreviated if long)
    species_label = ', '.join(sorted(entry['species']))
    if len(species_label) > 20:
        species_label = species_label[:18] + '…'
    ax.text(jpos[0], jpos[1] + 0.18, species_label,
            ha='center', va='bottom', fontsize=6,
            color=col, zorder=5, alpha=min(alpha + 0.1, 1.0))

    # Producer → junction
    _farrow(pos[lp], jpos, col, lw, alpha, shrinkA=SHRINK_ERC, shrinkB=SHRINK_JUNC)
    # Junction → consumer
    _farrow(jpos, pos[lc], col, lw, alpha, shrinkA=SHRINK_JUNC, shrinkB=SHRINK_ERC)


# -- Hover tooltip (ERC nodes) ------------------------------------------------
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
        req  = sorted(erc.get_required_species(RN)) if erc else []
        prod = sorted(erc.get_produced_species(RN)) if erc else []
        annot.xy = tuple(pos[node])
        annot.set_text(
            f"{node}\n"
            f"closure: {cl}\n"
            f"req: {req}\n"
            f"prod: {prod}"
        )
        annot.set_visible(True)
    else:
        annot.set_visible(False)
    fig.canvas.draw_idle()


fig.canvas.mpl_connect("motion_notify_event", _on_hover)

# -- Legends ------------------------------------------------------------------

# 1) Maintenance class
maintenance_legend = ax.legend(
    handles=[
        mpatches.Patch(facecolor=COL_SM,  label='Self-maintaining  (SM)'),
        mpatches.Patch(facecolor=COL_SSM, label='Semi-self-maintaining  (SSM)'),
        mpatches.Patch(facecolor=COL_DEF, label='Neither / not reactive'),
    ],
    loc='upper right', fontsize=10, framealpha=0.9,
    title='Maintenance class', title_fontsize=10,
)
ax.add_artist(maintenance_legend)

# 2) Complementarity type
shown_fund = sum(1 for e in comp_entries if e['level'] == 'fundamental')
shown_pure = sum(1 for e in comp_entries if e['level'] == 'pure')
shown_comp = sum(1 for e in comp_entries if e['level'] == 'complementary')

comp_handles = [
    Line2D([0], [0], color=COL_COMP['fundamental'], lw=2.5,
           label=f'Fundamental  ({shown_fund} directions)'),
]
if not SHOW_FUNDAMENTAL_ONLY:
    comp_handles += [
        Line2D([0], [0], color=COL_COMP['pure'], lw=1.8,
               label=f'Purely complementary  ({shown_pure})'),
        Line2D([0], [0], color=COL_COMP['complementary'], lw=1.2,
               label=f'Complementary (+synergy)  ({shown_comp})'),
    ]

# Square marker proxy for junction
junc_handle = ax.scatter([], [], s=JUNC_SIZE,
                         color='#888888', marker='s', alpha=0.9,
                         label='■ = supply junction  (prod→cons)')
comp_handles.append(junc_handle)

comp_legend = ax.legend(handles=comp_handles,
                        loc='upper left', fontsize=10,
                        framealpha=0.9,
                        title='Complementarity type', title_fontsize=10)
ax.add_artist(comp_legend)

# 3) Closure size
all_counts = sorted({len(erc.get_closure(RN)) for erc in ercs})
if len(all_counts) > 4:
    tick_idx   = [0, len(all_counts) // 3, 2 * len(all_counts) // 3, -1]
    size_ticks = sorted({all_counts[i] for i in tick_idx})
else:
    size_ticks = all_counts

size_handles = [
    ax.scatter([], [], s=NODE_SIZE_BASE + NODE_SIZE_SCALE * n,
               color='#888888', alpha=0.85, label=f'{n} species')
    for n in size_ticks
]
ax.legend(handles=size_handles,
          loc='lower left', fontsize=10, framealpha=0.9,
          title='Closure size', title_fontsize=10,
          labelspacing=1.2, handletextpad=1.0)

# -- Title --------------------------------------------------------------------
net_name = os.path.splitext(os.path.basename(RN_FILE))[0]
n_comp_pairs = len({frozenset([e['producer'].label, e['consumer'].label])
                    for e in comp_entries})
n_fund_pairs = len({frozenset([e['producer'].label, e['consumer'].label])
                    for e in comp_entries if e['level'] == 'fundamental'})
n_pure_pairs = len({frozenset([e['producer'].label, e['consumer'].label])
                    for e in comp_entries if e['level'] == 'pure'})

ax.set_title(
    f"ERC Hierarchy + Complementarity — {net_name}\n"
    f"Complementary pairs: {n_comp_pairs}   "
    f"Fundamental edges: {n_fund_pairs}   "
    f"Purely complementary: {n_pure_pairs}   "
    f"■ = supply junction",
    fontsize=11)
ax.axis('off')
plt.tight_layout()
_out_dir = os.path.join(_SCRIPT_DIR, '..', 'outputs', 'complementarity_viz')
os.makedirs(_out_dir, exist_ok=True)
_out_path = os.path.join(_out_dir, f'complementarity_viz_{net_name}.png')
plt.savefig(_out_path, dpi=150, bbox_inches='tight')
print(f"Saved: {os.path.abspath(_out_path)}")
plt.show()
