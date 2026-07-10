#!/usr/bin/env python3
"""
04_plot_erc_hierarchy_tree.py
==============================
Scatter plot: ERC hierarchy tree depth vs. ERCs per tree.

Unit of observation: one tree (connected component of the ERC hierarchy forest).
X = number of ERCs in the tree  (|E_T|, log scale)
Y = tree depth (longest root-to-leaf path, linear scale)
Size = mean branching (mean out-degree of non-leaf nodes; b=1→tiny, b=3→large)

Reads:  outputs/per_tree_stats.csv   (produced by 01_compute_erc_hierarchy.py)
Output: visualizations/erc_tree_depth_branching.png
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats as _stats

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _SCRIPT_DIR)

from config import (
    FONT_SCALE, FBASE_LAB, FBASE_TIT, FBASE_LEG, FBASE_TCK,
    FIG_W, FIG_H, OUT_DIR, VIZ_DIR, GROUP_STYLE,
)

os.makedirs(VIZ_DIR, exist_ok=True)

CSV_FILE = os.path.join(OUT_DIR, 'per_tree_stats.csv')
OUT_PNG  = os.path.join(VIZ_DIR, 'erc_tree_depth_branching.png')

# ── Load data ──────────────────────────────────────────────────────────────────
df = pd.read_csv(CSV_FILE)
print(f"Loaded {len(df)} tree records from {len(df['file'].unique())} networks")
print()
print(df[['n_ercs', 'n_nodes', 'depth', 'branching']].describe().round(3).to_string())
print()
print(f"Depth distribution:\n{df['depth'].value_counts().sort_index().to_string()}")
print()

import scipy.stats as _sc
for _c1, _c2 in [('depth', 'branching'), ('n_ercs', 'depth'), ('n_ercs', 'branching')]:
    _r, _p = _sc.pearsonr(df[_c1], df[_c2])
    print(f"  r({_c1}, {_c2}) = {_r:+.3f}  p={_p:.2e}")

# ── Per-tree depth scaling fit: depth ~ A * n_nodes^gamma ─────────────────────
_mask_fit = (df['n_nodes'] > 0) & (df['depth'] > 0)
_lx_f = np.log10(df.loc[_mask_fit, 'n_nodes'].astype(float))
_ly_f = np.log10(df.loc[_mask_fit, 'depth'].astype(float))
_sl, _ic, _r_f, _p_f, _se_f = _stats.linregress(_lx_f, _ly_f)
_A_fit  = 10 ** _ic
_gamma  = _sl
_ci_gam = 1.96 * _se_f
print(f"\nDepth scaling: depth ~ {_A_fit:.3f} * n_nodes^{_gamma:.3f}"
      f"  R2={_r_f**2:.3f}  p={_p_f:.2e}  CI: {_gamma:.3f}+/-{_ci_gam:.3f}")

# ── Binned medians (log-spaced over n_nodes) ───────────────────────────────────
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

# ── Marker area from mean branching ───────────────────────────────────────────
def branching_size(b):
    """b=1→15, b=1.5→57, b=2→135, b=3→354"""
    b = np.asarray(b, dtype=float)
    return np.maximum(15.0, 15.0 + 120.0 * (b - 1.0) ** 1.5)

# ── Scaled font sizes ──────────────────────────────────────────────────────────
_FL  = int(round(FBASE_LAB * FONT_SCALE))
_FT  = int(round(FBASE_TIT * FONT_SCALE))
_FLG = int(round(FBASE_LEG * FONT_SCALE))
_FCK = int(round(FBASE_TCK * FONT_SCALE))

# ── Figure ────────────────────────────────────────────────────────────────────
np.random.seed(42)
fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))

for grp, sty in GROUP_STYLE.items():
    sub = df[df['group'] == grp]
    if len(sub) == 0:
        continue
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

# ── Binned medians ± IQR ──────────────────────────────────────────────────────
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
_bm_h = Line2D([0], [0], marker='D', color='#555555', ms=8, lw=0,
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

# ── Group + fit legend (upper left) ───────────────────────────────────────────
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

# ── Branching size legend (lower right) ───────────────────────────────────────
_b_ticks   = [1.0, 1.5, 2.0, 3.0]
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
plt.savefig(OUT_PNG, dpi=150, bbox_inches='tight')
print(f"\nSaved: {os.path.abspath(OUT_PNG)}")
