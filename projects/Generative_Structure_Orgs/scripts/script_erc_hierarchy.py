#!/usr/bin/env python3
"""
script_erc_hierarchy.py
=======================
Compute the ERCs of a reaction network and plot the ERC hierarchy graph.

Node size  ∝  number of species in the ERC closure.
Node color:
  green  (#27AE60)  — self-maintaining  (SM)
  orange (#E67E22)  — semi-self-maintaining but not SM  (SSM only)
  steel  (#AED6F1)  — neither

Definitions
-----------
SSM : every species consumed by the enabled reactions is also produced.
      Checked combinatorially (no LP needed) — cheap.

SM  : there exists a flux vector v > 0 over ALL enabled reactions
      such that S·v ≥ 0 (COT definition: every species has non-negative
      net production when all reactions are active).
      LP check with lower bound ε on each flux — only run if SSM=True.

Note on the ε lower bound
--------------------------
The library's minimize_sv() uses v ≥ 0, which allows individual reactions
to have zero flux.  This gives a false positive for sets where the only
non-negative-production flux has some reactions silent (e.g. a wasteful
2-cycle that must be excluded to keep the rest balanced).  The COT
definition (Dittrich & di Fenizio 2007) requires v > 0 — all enabled
reactions must fire.  The local _is_sm() below enforces v ≥ SM_EPS.

Default network: data/biomodels/biomodels_interesting/BIOMD0000000237_manyOrgs.txt
"""

import os
import sys
import time
from collections import defaultdict

import numpy as np
from scipy.optimize import linprog
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
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

# Uncomment to use the Farm network instead:
#RN_FILE = os.path.join(_PYCOT_ROOT, 'networks', 'testing', 'Farm.txt')

RN_FILE = '../../networks/testing/LBCA.txt'

# -- Visual parameters ---------------------------------------------------------
NODE_SIZE_BASE  = 400   # minimum node area (matplotlib scatter units)
NODE_SIZE_SCALE = 120   # added per species in the closure
COL_SM  = '#27AE60'     # green  — self-maintaining
COL_SSM = '#E67E22'     # orange — semi-self-maintaining (not SM)
COL_DEF = '#AED6F1'     # steel blue — neither
fontsize_default = 14
# Minimum flux for each reaction in the SM LP (enforces v > 0 per COT def.)
SM_EPS = 1e-6

# ===========================================================================
# Reactivity check  (cheapest — runs before SSM and SM)
# ===========================================================================

def _is_reactive(RN, species_set):
    """
    A set is reactive iff it enables at least one non-inflow reaction.
    Inflow reactions (empty support) are excluded: they fire regardless of
    the species set and do not make it 'reactive' in the COT sense.
    The empty set is never reactive.
    """
    if not species_set:
        return False
    enabled = RN.get_reactions_from_species(species_set)
    return any(r.support_indices() for r in enabled)


# ===========================================================================
# Local SM check  (v ≥ SM_EPS — all enabled reactions must fire)
# ===========================================================================

def _is_sm(RN, species_set):
    """
    COT self-maintenance: ∃ v > 0 over all reactions enabled by species_set
    such that S·v ≥ 0.

    Uses scipy linprog with bounds (SM_EPS, None) so every reaction must
    carry at least SM_EPS flux.  This prevents the false-positive that
    arises when the only non-negative-production solution sets some
    reactions to zero (effectively ignoring part of the species set).

    Returns True iff the LP is feasible.
    """
    sub_RN = RN.sub_reaction_network(species_set)
    S      = np.asarray(sub_RN.stoichiometry_matrix(), dtype=float)
    n_sp, n_rx = S.shape

    if n_rx == 0:
        return False

    # Feasibility LP:  find v ≥ SM_EPS  s.t.  S·v ≥ 0
    # Rewrite as:  -S·v ≤ 0  with  v ≥ SM_EPS
    # Shift:  let w = v - SM_EPS  (w ≥ 0).
    # Then  S·(w + SM_EPS·1) ≥ 0  →  S·w ≥ -SM_EPS · S·1
    # i.e.  -S·w ≤ SM_EPS · S·1   (component-wise)
    rhs    = SM_EPS * S.sum(axis=1)          # S · (SM_EPS * ones)
    A_ub   = -S
    b_ub   = rhs
    bounds = [(0, None)] * n_rx
    # Objective: minimise 0 (pure feasibility)
    c      = np.zeros(n_rx)

    res = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method='highs')
    return res.status == 0   # 0 = optimal (feasible)


# ===========================================================================
# Load network
# ===========================================================================
print(f"Loading network: {RN_FILE}")
RN = read_txt(RN_FILE)
print(f"  Species:   {len(RN.species())}")
print(f"  Reactions: {len(RN.reactions())}")

# ===========================================================================
# Compute ERCs
# ===========================================================================
print("\nComputing ERCs...")
t0 = time.time()
ercs = ERC.ERCs(RN)
print(f"  Found {len(ercs)} ERCs  ({time.time() - t0:.1f}s)")

for erc in ercs:
    min_gens = species_list_to_names(erc.min_generators)
    cl       = species_list_to_names(erc.get_closure(RN))
    print(f"\n  ERC {erc.label}")
    print(f"    Min generators : {min_gens}")
    print(f"    Closure        : {cl}")

# ===========================================================================
# Build hierarchy and classify nodes
# SSM first (cheap, combinatorial); SM only if SSM passes (LP, expensive).
# ===========================================================================
print("\nBuilding ERC hierarchy...")
hierarchy = ERC_Hierarchy(RN, ercs)
G = hierarchy.graph

print("\nClassifying nodes (SSM → SM)...")
node_sizes  = {}
node_colors = {}

for erc in ercs:
    label   = erc.label
    closure = erc.get_closure(RN)
    n_sp    = len(closure)

    node_sizes[label] = NODE_SIZE_BASE + NODE_SIZE_SCALE * n_sp

    # Step 1: reactive — cheapest check (empty set or no enabled reactions)
    if not _is_reactive(RN, closure):
        node_colors[label] = COL_DEF
        tag = 'not reactive'
    else:
        # Step 2: SSM — cheap combinatorial check
        is_ssm = is_semi_self_maintaining(RN, closure)

        if not is_ssm:
            node_colors[label] = COL_DEF
            tag = 'neither'
        else:
            # Step 3: SM — LP check only for SSM=True nodes
            is_sm = _is_sm(RN, closure)
            if is_sm:
                node_colors[label] = COL_SM
                tag = 'self-maintaining'
            else:
                node_colors[label] = COL_SSM
                tag = 'semi-self-maintaining only'

    print(f"  ERC {label:>4}  |  {n_sp:2d} species  |  {tag}")

# ===========================================================================
# Layout: level-based (same algorithm as ERC_Hierarchy.plot_hierarchy)
# ===========================================================================
levels      = ERC.get_node_levels(G)
level_nodes = defaultdict(list)
for node, lvl in levels.items():
    level_nodes[lvl].append(node)

pos = {}
for lvl, nodes in level_nodes.items():
    nodes.sort(key=lambda n: len(nx.ancestors(G, n)), reverse=True)
    for i, node in enumerate(nodes):
        pos[node] = ((i - (len(nodes) - 1) / 2) * 2.0, lvl * 2.0)

# ===========================================================================
# Plot
# ===========================================================================
ordered_nodes = list(G.nodes())
sizes_list    = [node_sizes.get(n, NODE_SIZE_BASE) for n in ordered_nodes]
colors_list   = [node_colors.get(n, COL_DEF)        for n in ordered_nodes]

fig, ax = plt.subplots(figsize=(12, 8))

scatter = nx.draw_networkx_nodes(G, pos, ax=ax,
                                 nodelist=ordered_nodes,
                                 node_color=colors_list,
                                 node_size=sizes_list,
                                 alpha=0.90)

nx.draw_networkx_edges(G, pos, ax=ax,
                       edge_color='#555555',
                       arrows=True, arrowsize=18)

nx.draw_networkx_labels(G, pos, ax=ax,
                        font_size=10, font_weight='bold')

# -- Hover annotation ---------------------------------------------------------
annot = ax.annotate("", xy=(0, 0), xytext=(15, 15),
                    textcoords="offset points",
                    bbox=dict(boxstyle="round", fc="white", ec="0.5", alpha=0.9),
                    visible=False)

erc_by_label = {erc.label: erc for erc in ercs}

def _on_hover(event):
    if event.inaxes != ax:
        return
    cont, ind = scatter.contains(event)
    if cont:
        node = ordered_nodes[ind["ind"][0]]
        erc  = erc_by_label.get(node)
        cl   = species_list_to_names(erc.get_closure(RN)) if erc else []
        annot.xy = pos[node]
        annot.set_text(f"{node}\n{cl}")
        annot.set_visible(True)
    else:
        annot.set_visible(False)
    fig.canvas.draw_idle()

fig.canvas.mpl_connect("motion_notify_event", _on_hover)

# -- Colour legend ------------------------------------------------------------
from matplotlib.patches import Patch
colour_legend = ax.legend(
    handles=[
        Patch(facecolor=COL_SM,  label='Self-maintaining  (SM)'),
        Patch(facecolor=COL_SSM, label='Semi-self-maintaining  (SSM only)'),
        Patch(facecolor=COL_DEF, label='Neither / not reactive'),
    ],
    loc='upper right', fontsize=fontsize_default, framealpha=0.9,
    title='Maintenance class', title_fontsize=fontsize_default,
)
ax.add_artist(colour_legend)   # keep it when the size legend is added

# -- Size legend --------------------------------------------------------------
# Pick 3–4 representative species counts spanning the actual range.
all_counts = sorted({len(erc.get_closure(RN)) for erc in ercs})
if len(all_counts) <= 4:
    size_ticks = all_counts
else:
    lo, hi = all_counts[0], all_counts[-1]
    mid1   = all_counts[len(all_counts) // 3]
    mid2   = all_counts[2 * len(all_counts) // 3]
    size_ticks = sorted({lo, mid1, mid2, hi})

size_handles = [
    ax.scatter([], [], s=NODE_SIZE_BASE + NODE_SIZE_SCALE * n,
               color='#888888', alpha=0.85, label=f'{n} species')
    for n in size_ticks
]
ax.legend(handles=size_handles,
          loc='upper left', fontsize=fontsize_default, framealpha=0.9,
          title='Closure size', title_fontsize=fontsize_default,
          labelspacing=1.2, handletextpad=1.0)

net_name = os.path.splitext(os.path.basename(RN_FILE))[0]
ax.set_title(f"ERC Hierarchy — {net_name}\n"
             f"Node size ∝ closure size  ·  colour = maintenance class",
             fontsize=fontsize_default)
ax.axis('off')
plt.tight_layout()
plt.show()
