#!/usr/bin/env python3
"""
script_erc_hierarchy_tree_viz.py
=================================
Scatter plot: ERC hierarchy tree depth vs. branching.

Unit of observation: ONE TREE (connected component of the ERC hierarchy forest).
A single reaction network is a FOREST that may contain several trees; each tree
with >= 3 ERC nodes contributes one point to the plot.

  X    = tree depth     (longest root-to-leaf path within the tree)
  Y    = mean branching (mean out-degree of non-leaf nodes = avg children
                         per internal ERC node)
  Size = number of ERC nodes in the tree

Groups: BioModels (red circles), BiGG (dark blue stars).
Horizontal jitter added because depth is an integer.

Outputs:
  outputs/forest_stats/per_tree_stats.csv
  visualizations/forest_stats/erc_tree_depth_branching.png
"""

import os
import sys
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_PYCOT_ROOT = os.path.normpath(os.path.join(_SCRIPT_DIR, '..', '..', '..'))
sys.path.insert(0, os.path.join(_PYCOT_ROOT, 'src'))
sys.path.insert(0, _SCRIPT_DIR)

from pyCOT.io.functions import read_txt
from pyCOT.analysis.ERC_Hierarchy import ERC, ERC_Hierarchy
from utils_ercs import load_ercs

# ── Configuration ─────────────────────────────────────────────────────────────
_BIOMD = os.path.join(_PYCOT_ROOT, 'data', 'biomodels')

SCAN_FOLDERS = {
    os.path.join(_BIOMD, 'BioMD_metabolic'):       'BioMD_metabolic',
    os.path.join(_BIOMD, 'BioMD_cell_cycle'):       'BioMD_cell_cycle',
    os.path.join(_BIOMD, 'BioMD_circadian'):        'BioMD_circadian',
    os.path.join(_BIOMD, 'BioMD_signaling'):        'BioMD_signaling',
    os.path.join(_BIOMD, 'BioMD_gene_regulation'):  'BioMD_gene_regulation',
    os.path.join(_BIOMD, 'BioMD_apoptosis'):        'BioMD_apoptosis',
    os.path.join(_BIOMD, 'BioMD_immune'):           'BioMD_immune',
    os.path.join(_BIOMD, 'BioMD_other'):            'BioMD_other',
    os.path.join(_BIOMD, 'BiGG'):                   'BiGG',
    os.path.join(_BIOMD, 'Other'):                  'Other',
}

MAX_REACTIONS = 1000  # skip networks with more than this many reactions (pre-ERC filter)
MAX_ERCS      = 900   # skip networks with more than this many ERCs (post-ERC filter)
MIN_NODES     = 3     # trees with fewer ERC nodes are excluded
MIN_ERCS      = 4     # networks with fewer ERCs are skipped entirely

# ── Font and figure size (adjust here to rescale the whole plot) ──────────────
FONT_SCALE = 2.0          # multiply all font sizes by this factor
_FBASE_LAB = 12           # base axis-label fontsize  (pt, before scaling)
_FBASE_TIT = 11           # base title fontsize
_FBASE_LEG = 9            # base legend fontsize
_FBASE_TCK = 10           # base tick-label fontsize
FIG_W, FIG_H = 14, 10    # figure size in inches (increase if text overlaps)

OUT_DIR = os.path.normpath(os.path.join(_SCRIPT_DIR, '..', 'outputs', 'forest_stats'))
VIZ_DIR = os.path.normpath(os.path.join(_SCRIPT_DIR, '..', 'visualizations', 'forest_stats'))
CSV_OUT = os.path.join(OUT_DIR, 'per_tree_stats.csv')
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(VIZ_DIR, exist_ok=True)


# ── Helpers ───────────────────────────────────────────────────────────────────

def assign_group(dataset):
    if str(dataset).startswith('BioMD_'):
        return 'BioModels'
    if dataset == 'BiGG':
        return 'BiGG'
    return 'Other'


def _node_depths(G):
    """Longest path from any root (in-degree 0) down to each node. Root = 0."""
    depth = {n: 0 for n in G.nodes()}
    for node in nx.topological_sort(G):
        for child in G.successors(node):
            if depth[node] + 1 > depth[child]:
                depth[child] = depth[node] + 1
    return depth


def per_tree_stats(G_sub):
    """
    Stats for one connected component (tree) of the ERC hierarchy DAG.

    Returns (depth, branching, n_nodes):
      depth     — longest path from root ERC to leaf ERC
      branching — mean out-degree of non-leaf nodes (nodes that have children)
      n_nodes   — number of ERC nodes in this tree
    """
    n_nodes = G_sub.number_of_nodes()
    depths = _node_depths(G_sub)
    depth = max(depths.values()) if depths else 0
    out_degrees = [G_sub.out_degree(n) for n in G_sub.nodes() if G_sub.out_degree(n) > 0]
    branching = float(np.mean(out_degrees)) if out_degrees else 0.0
    return depth, branching, n_nodes


# ── Compute per-tree statistics ───────────────────────────────────────────────

records = []

for folder, dataset_label in SCAN_FOLDERS.items():
    if not os.path.isdir(folder):
        continue
    txt_files = sorted(f for f in os.listdir(folder) if f.endswith('.txt'))
    print(f"\n=== {dataset_label}: {len(txt_files)} files ===")

    for fname in txt_files:
        txt_path = os.path.join(folder, fname)
        try:
            RN = read_txt(txt_path)
            n_rxn = len(RN.reactions())
            if n_rxn > MAX_REACTIONS:
                print(f"  SKIP {fname}: {n_rxn} reactions > {MAX_REACTIONS}")
                continue
            ercs, _ = load_ercs(txt_path, RN, ERC)
            # Exclude E_∅ (empty-closure ERC)
            ercs = [e for e in ercs if len(e.get_closure_names(RN)) > 0]
            if len(ercs) < MIN_ERCS:
                continue
            if len(ercs) > MAX_ERCS:
                print(f"  SKIP {fname}: {len(ercs)} ERCs > {MAX_ERCS}")
                continue

            hierarchy = ERC_Hierarchy(RN, ercs)
            G = hierarchy.graph
            G_undir = G.to_undirected()

            n_qualifying = 0
            for comp in nx.connected_components(G_undir):
                if len(comp) < MIN_NODES:
                    continue
                sub = G.subgraph(comp)
                depth, branching, n_nodes = per_tree_stats(sub)
                if depth == 0:
                    # Linear chain of 1 level — no internal branching structure
                    continue
                records.append({
                    'file':      fname,
                    'dataset':   dataset_label,
                    'group':     assign_group(dataset_label),
                    'n_nodes':   n_nodes,
                    'depth':     depth,
                    'branching': round(branching, 4),
                })
                n_qualifying += 1

            print(f"  {fname}: {len(ercs)} ERCs → {n_qualifying} qualifying trees")

        except Exception as exc:
            print(f"  ERR {fname}: {exc}")

df = pd.DataFrame(records)
df.to_csv(CSV_OUT, index=False)
print(f"\nSaved {len(df)} tree records to {CSV_OUT}")
print()
print(df[['n_nodes', 'depth', 'branching']].describe().round(3).to_string())
print()
print(f"Depth distribution:\n{df['depth'].value_counts().sort_index().to_string()}")


# ── Plot ──────────────────────────────────────────────────────────────────────

GROUP_STYLE = {
    'BioModels': {'color': '#E74C3C', 'marker': 'o', 'label': 'BioModels'},
    'BiGG':      {'color': '#2C3E50', 'marker': '*', 'label': 'BiGG'},
    'Other':     {'color': '#7F8C8D', 'marker': 's', 'label': 'Other'},
}

# Size: scale logarithmically with n_nodes
n_min = df['n_nodes'].min()
n_max = df['n_nodes'].max()
log_min = np.log(n_min)
log_max = np.log(max(n_max, n_min + 1))

def node_size(n):
    """Map n_nodes -> scatter marker area (pts²). Range: 30–500."""
    t = (np.log(np.asarray(n, dtype=float)) - log_min) / max(log_max - log_min, 1)
    return 30 + 470 * t

np.random.seed(42)

_FL  = int(round(_FBASE_LAB * FONT_SCALE))
_FT  = int(round(_FBASE_TIT * FONT_SCALE))
_FLG = int(round(_FBASE_LEG * FONT_SCALE))
_FCK = int(round(_FBASE_TCK * FONT_SCALE))

fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))

for grp, sty in GROUP_STYLE.items():
    sub = df[df['group'] == grp]
    if len(sub) == 0:
        continue
    jitter = np.random.normal(0, 0.18, size=len(sub))
    x = sub['depth'].values + jitter
    y = sub['branching'].values
    s = node_size(sub['n_nodes'].values)

    ax.scatter(x, y,
               s=s, c=sty['color'], marker=sty['marker'],
               alpha=0.55, edgecolors='white', linewidths=0.3,
               zorder=3, label=sty['label'])

# ── Depth tick marks at integer positions ────────────────────────────────────
depth_vals = sorted(df['depth'].unique())
ax.set_xticks(depth_vals)
ax.set_xticklabels([str(d) for d in depth_vals], fontsize=_FCK)
ax.tick_params(axis='y', labelsize=_FCK)

ax.set_xlabel('Tree depth  (longest root-to-leaf path within the tree)', fontsize=_FL)
ax.set_ylabel('Mean branching  (avg children per internal ERC node)', fontsize=_FL)
ax.set_title(
    f'ERC hierarchy: depth vs. branching per tree  '
    f'({len(df)} trees from {df["file"].nunique()} networks,  '
    f'$\\geq {MIN_NODES}$ nodes)',
    fontsize=_FT, fontweight='bold',
)

# ── Size legend ───────────────────────────────────────────────────────────────
size_ticks = [3, 10, 50, 200]
size_ticks = [v for v in size_ticks if v <= n_max]
size_handles = [
    ax.scatter([], [], s=node_size(v), c='#888888', alpha=0.7,
               edgecolors='white', lw=0.3, label=f'{v} ERCs')
    for v in size_ticks
]

# ── Group legend ──────────────────────────────────────────────────────────────
group_handles = [
    Line2D([0], [0],
           marker=sty['marker'], color='w',
           markerfacecolor=sty['color'], markeredgecolor='white',
           markersize=9, label=sty['label'])
    for grp, sty in GROUP_STYLE.items()
    if grp in df['group'].values
]

leg1 = ax.legend(handles=group_handles,
                 title='Database', title_fontsize=_FLG,
                 fontsize=_FLG, loc='upper right', framealpha=0.90)
ax.add_artist(leg1)
ax.legend(handles=size_handles,
          title='Tree size (ERCs)', title_fontsize=_FLG,
          fontsize=_FLG, loc='upper left', framealpha=0.90)

ax.grid(True, alpha=0.20, which='both')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
out_path = os.path.join(VIZ_DIR, 'erc_tree_depth_branching.png')
plt.savefig(out_path, dpi=150, bbox_inches='tight')
print(f"\nSaved: {os.path.abspath(out_path)}")
