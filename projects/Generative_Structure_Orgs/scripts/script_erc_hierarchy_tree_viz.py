#!/usr/bin/env python3
"""
script_erc_hierarchy_tree_viz.py
=================================
Scatter plot: ERC hierarchy tree depth vs. network ERC count.

Unit of observation: ONE TREE (connected component of the ERC hierarchy forest).
A single reaction network is a FOREST that may contain several trees; each tree
with >= 3 ERC nodes contributes one point to the plot.

  X    = network ERC count  (n_ercs, log scale)
  Y    = tree depth         (longest root-to-leaf path, log scale)
  Size = mean branching     (mean out-degree of non-leaf ERC nodes;
                             b=1 → tiny, b=3 → large)

Groups: BioModels (red circles), BiGG (dark blue stars).

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
                    'n_ercs':    len(ercs),   # network-level total ERCs
                    'depth':     depth,
                    'branching': round(branching, 4),
                })
                n_qualifying += 1

            print(f"  {fname}: {len(ercs)} ERCs -> {n_qualifying} qualifying trees")

        except Exception as exc:
            print(f"  ERR {fname}: {exc}")

df = pd.DataFrame(records)
df.to_csv(CSV_OUT, index=False)
print(f"\nSaved {len(df)} tree records to {CSV_OUT}")
print()
print(df[['n_ercs', 'n_nodes', 'depth', 'branching']].describe().round(3).to_string())
print()
print(f"Depth distribution:\n{df['depth'].value_counts().sort_index().to_string()}")
print()
# Quick correlation table
import scipy.stats as _sc
for _c1, _c2 in [('depth', 'branching'), ('n_ercs', 'depth'), ('n_ercs', 'branching')]:
    _r, _p = _sc.pearsonr(df[_c1], df[_c2])
    print(f"  r({_c1}, {_c2}) = {_r:+.3f}  p={_p:.2e}")


# ── Plot: X = n_ercs (log), Y = depth (log), size = branching ────────────────

GROUP_STYLE = {
    'BioModels': {'color': '#E74C3C', 'marker': 'o', 'label': 'BioModels', 'lw': 0.3},
    'BiGG':      {'color': '#2C3E50', 'marker': '*', 'label': 'BiGG',      'lw': 1.5},
    'Other':     {'color': '#7F8C8D', 'marker': 's', 'label': 'Other',     'lw': 0.3},
}


def branching_size(b):
    """Marker area from mean branching.  b=1->15, b=1.5->57, b=2->135, b=3->354."""
    b = np.asarray(b, dtype=float)
    return np.maximum(15.0, 15.0 + 120.0 * (b - 1.0) ** 1.5)


# ── Per-tree depth scaling fit: depth ~ A * n_nodes^gamma ────────────────────
from scipy import stats as _stats
_mask_fit = (df['n_nodes'] > 0) & (df['depth'] > 0)
_lx_f = np.log10(df.loc[_mask_fit, 'n_nodes'].astype(float))
_ly_f = np.log10(df.loc[_mask_fit, 'depth'].astype(float))
_sl, _ic, _r_f, _p_f, _se_f = _stats.linregress(_lx_f, _ly_f)
_A_fit  = 10 ** _ic
_gamma  = _sl
_ci_gam = 1.96 * _se_f
print(f"\nDepth scaling: depth ~ {_A_fit:.3f} * n_nodes^{_gamma:.3f}"
      f"  R2={_r_f**2:.3f}  p={_p_f:.2e}  CI: {_gamma:.3f}+/-{_ci_gam:.3f}")

# ── Binned medians of depth by n_nodes (log-spaced, per tree) ─────────────────
_edges = np.logspace(np.log10(df['n_nodes'].min()),
                     np.log10(df['n_nodes'].max() * 1.001), 8)
_bin_rows = []
for _lo, _hi in zip(_edges[:-1], _edges[1:]):
    _idx = (df['n_nodes'] >= _lo) & (df['n_nodes'] < _hi)
    if _idx.sum() < 3:
        continue
    _sd = df.loc[_idx, 'depth']
    _bin_rows.append({
        'x':   float(np.sqrt(_lo * _hi)),
        'med': float(_sd.median()),
        'q25': float(_sd.quantile(0.25)),
        'q75': float(_sd.quantile(0.75)),
        'n':   int(_idx.sum()),
    })

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
    # Multiplicative jitter on x (log scale); small additive jitter on y (linear scale)
    _xj = sub['n_nodes'].values.astype(float) * np.exp(np.random.normal(0, 0.05, len(sub)))
    _yj = sub['depth'].values.astype(float) + np.random.normal(0, 0.12, len(sub))
    _s  = branching_size(sub['branching'].values)
    ax.scatter(_xj, _yj,
               s=_s, c=sty['color'], marker=sty['marker'],
               alpha=0.50, edgecolors='white', linewidths=sty['lw'],
               zorder=3)

# ── Fit curve ─────────────────────────────────────────────────────────────────
_n_range = np.logspace(np.log10(df['n_nodes'].min()), np.log10(df['n_nodes'].max()), 300)
_d_range = _A_fit * _n_range ** _gamma
_lbl_fit = (f'Fit: $d = {_A_fit:.2f}\\cdot |\\mathcal{{E}}_T|^{{{_gamma:.2f}}}$'
            f'  ($R^2 = {_r_f**2:.2f}$)')
_fit_h, = ax.plot(_n_range, _d_range, color='#111111', lw=2.5, ls='--',
                   zorder=6, label=_lbl_fit)

# ── Binned medians +/- IQR ────────────────────────────────────────────────────
from matplotlib.lines import Line2D as _L2D
for _row in _bin_rows:
    _fill = _row['n'] >= 10
    ax.errorbar(_row['x'], _row['med'],
                yerr=[[max(_row['med'] - _row['q25'], 0)],
                      [max(_row['q75'] - _row['med'], 0)]],
                fmt='D', color='#555555',
                ms=int(7 * FONT_SCALE * 0.6),
                lw=1.8 * FONT_SCALE * 0.5,
                capsize=4 * FONT_SCALE * 0.5,
                elinewidth=1.5 * FONT_SCALE * 0.5,
                mfc='#555555' if _fill else 'white',
                mec='#555555', mew=1.5,
                zorder=7)
_bm_h = _L2D([0], [0], marker='D', color='#555555', ms=8, lw=0,
             mfc='#555555', mec='#555555',
             label='Median $\\pm$ IQR per bin\n(hollow: $n < 10$ networks)')

# ── Axes ──────────────────────────────────────────────────────────────────────
ax.set_xscale('log')
ax.set_yscale('linear')
ax.set_xlabel('Number of ERCs in tree  ($|\\mathcal{E}_T|$)', fontsize=_FL)
ax.set_ylabel('Tree depth  (longest root-to-leaf path)', fontsize=_FL)
ax.set_title(
    f'ERC hierarchy: depth scales with tree size  '
    f'({len(df)} trees, {df["file"].nunique()} networks;  '
    f'point size $\\propto$ branching)',
    fontsize=_FT, fontweight='bold',
)
ax.tick_params(labelsize=_FCK)

# ── Group + fit legend (upper left) ──────────────────────────────────────────
_grp_handles = [
    Line2D([0], [0], marker=sty['marker'], color='w',
           markerfacecolor=sty['color'], markeredgecolor=sty['color'],
           markersize=9, label=sty['label'])
    for grp, sty in GROUP_STYLE.items()
    if grp in df['group'].values
]
leg1 = ax.legend(handles=_grp_handles + [_fit_h, _bm_h],
                 title='Database / Fit', title_fontsize=_FLG,
                 fontsize=_FLG, loc='upper left', framealpha=0.90)
ax.add_artist(leg1)

# ── Branching size legend (lower right) ──────────────────────────────────────
_b_ticks = [1.0, 1.5, 2.0, 3.0]
_sz_handles = [
    ax.scatter([], [], s=branching_size(b), c='#888888', alpha=0.7,
               edgecolors='white', lw=0.3, label=f'$b = {b:.1f}$')
    for b in _b_ticks
]
ax.legend(handles=_sz_handles,
          title='Mean branching  ($b$)', title_fontsize=_FLG,
          fontsize=_FLG, loc='lower right', framealpha=0.90)

ax.grid(True, alpha=0.20, which='both')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
out_path = os.path.join(VIZ_DIR, 'erc_tree_depth_branching.png')
plt.savefig(out_path, dpi=150, bbox_inches='tight')
print(f"\nSaved: {os.path.abspath(out_path)}")
