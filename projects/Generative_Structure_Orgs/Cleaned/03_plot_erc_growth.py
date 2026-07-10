#!/usr/bin/env python3
"""
03_plot_erc_growth.py
======================
Plot ERC count vs. reaction count (log-log) for all networks.

Reads:  outputs/forest_stats.csv   (produced by 01_compute_erc_hierarchy.py)
Output: visualizations/erc_growth_tendency.png

Groups: BioModels (red circles), BiGG (dark-blue stars).
Each group: scatter + per-group power-law fit line + binned medians ± IQR.
Overall fit across all networks: dashed black line.
Reference line: ERCs = reactions (diagonal).
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _SCRIPT_DIR)

from config import (
    FONT_SCALE, FBASE_LAB, FBASE_TIT, FBASE_LEG, FBASE_TCK,
    FIG_W, FIG_H, OUT_DIR, VIZ_DIR, assign_group,
)

import os
os.makedirs(VIZ_DIR, exist_ok=True)

CSV_FILE = os.path.join(OUT_DIR, 'forest_stats.csv')
OUT_PNG  = os.path.join(VIZ_DIR, 'erc_growth_tendency.png')

# ── Groups ─────────────────────────────────────────────────────────────────────
# (name, color, marker, scatter_size, linewidth)
GROUPS = [
    ('BioModels', '#E74C3C', 'o',  80, 0),
    ('BiGG',      '#2C3E50', '*', 500, 1.5),
]

# ── Load & filter ──────────────────────────────────────────────────────────────
df_all = pd.read_csv(CSV_FILE)
df = df_all[(df_all['n_ercs'] >= 4) & (df_all['n_reactions'] > 0)].copy()
df['group'] = df['dataset'].map(assign_group)
print(f"Loaded {len(df)} networks (n_ercs >= 4)")
print(df.groupby('group')['n_ercs'].count().to_string())
print()

# ── Power-law fit (OLS on log-log) ────────────────────────────────────────────
def powerlaw_fit(x_arr, y_arr):
    x_arr = np.asarray(x_arr, dtype=float)
    y_arr = np.asarray(y_arr, dtype=float)
    mask  = (y_arr > 0) & (x_arr > 0)
    if mask.sum() < 5:
        return None
    slope, intercept, r, _, se = stats.linregress(np.log10(x_arr[mask]),
                                                   np.log10(y_arr[mask]))
    return {'alpha': slope, 'logA': intercept, 'R2': r**2,
            'n': int(mask.sum()), 'se_alpha': se}

# ── Binned medians (log-spaced over x) ────────────────────────────────────────
def binned_medians(x_arr, y_arr, n_bins=7):
    x_arr = np.asarray(x_arr, dtype=float)
    y_arr = np.asarray(y_arr, dtype=float)
    mask  = (x_arr > 0) & (y_arr > 0)
    x, y  = x_arr[mask], y_arr[mask]
    if len(x) < 5:
        return []
    edges = np.logspace(np.log10(x.min()), np.log10(x.max() * 1.001), n_bins + 1)
    rows  = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        idx = (x >= lo) & (x < hi)
        if idx.sum() < 3:
            continue
        rows.append({
            'x':   float(np.median(x[idx])),
            'med': float(np.median(y[idx])),
            'q25': float(np.percentile(y[idx], 25)),
            'q75': float(np.percentile(y[idx], 75)),
            'n':   int(idx.sum()),
        })
    return rows

# ── Scaled font sizes ──────────────────────────────────────────────────────────
_FL  = int(round(FBASE_LAB * FONT_SCALE))
_FT  = int(round(FBASE_TIT * FONT_SCALE))
_FLG = int(round(FBASE_LEG * FONT_SCALE))
_FCK = int(round(FBASE_TCK * FONT_SCALE))

# ── Figure ────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
handles = []
labels  = []

print("Power-law fits  n_ercs ~ A · n_reactions^alpha  (OLS on log-log):")

for grp_name, color, marker, ms, lw in GROUPS:
    sub = df[df['group'] == grp_name]
    if len(sub) == 0:
        continue

    x = sub['n_reactions'].values.astype(float)
    y = sub['n_ercs'].values.astype(float)

    ax.scatter(x, y, color=color, marker=marker, s=ms,
               alpha=0.45, zorder=2, linewidths=lw)

    fit = powerlaw_fit(x, y)
    if fit is not None:
        A   = 10 ** fit['logA']
        x_f = np.logspace(np.log10(max(x.min(), 0.5)), np.log10(x.max()), 250)
        ax.plot(x_f, A * x_f ** fit['alpha'], color=color, lw=2.2, zorder=3, alpha=0.92)
        ci = 1.96 * fit['se_alpha']
        print(f"  {grp_name:<12}  alpha={fit['alpha']:.3f} +/- {ci:.3f}"
              f"  R2={fit['R2']:.3f}  n={fit['n']}")
        lbl = f"{grp_name}  ($\\alpha = {fit['alpha']:.2f}$,  n={fit['n']})"
    else:
        print(f"  {grp_name:<12}  -- sparse --")
        lbl = f"{grp_name}  (n={len(sub)}, sparse)"

    bm = binned_medians(x, y, n_bins=6)
    for row in bm:
        ax.errorbar(row['x'], row['med'],
                    yerr=[[max(row['med'] - row['q25'], 0)],
                          [max(row['q75'] - row['med'], 0)]],
                    fmt=marker, color=color, ms=8, lw=1.6,
                    capsize=4, zorder=5, elinewidth=1.4,
                    mfc='white', mec=color, mew=1.5)

    handles.append(Line2D([0], [0], color=color, lw=2.2,
                           marker=marker, ms=6, mfc=color, mec=color))
    labels.append(lbl)

# ── Overall fit ───────────────────────────────────────────────────────────────
x_all = df['n_reactions'].values.astype(float)
y_all = df['n_ercs'].values.astype(float)
fit_all = powerlaw_fit(x_all, y_all)
if fit_all is not None:
    A_all = 10 ** fit_all['logA']
    x_fa  = np.logspace(np.log10(max(x_all[x_all > 0].min(), 0.5)),
                        np.log10(x_all.max()), 250)
    line_all, = ax.plot(x_fa, A_all * x_fa ** fit_all['alpha'],
                        color='#111111', lw=2.5, ls='--', zorder=6, alpha=0.88)
    handles.append(line_all)
    labels.append(f"All  ($\\alpha = {fit_all['alpha']:.2f}$,  n={fit_all['n']})")
    ci_all = 1.96 * fit_all['se_alpha']
    print(f"  {'All':<12}  alpha={fit_all['alpha']:.3f} +/- {ci_all:.3f}"
          f"  R2={fit_all['R2']:.3f}  n={fit_all['n']}")

bm_proxy = Line2D([0], [0], marker='o', color='#555555', ms=7, lw=0,
                  mfc='white', mec='#555555', mew=1.4)
handles.append(bm_proxy)
labels.append('Median ± IQR (binned, hollow)')

xr = np.array([max(x_all[x_all > 0].min(), 0.5), x_all.max()])
ref_h, = ax.plot(xr, xr, ':', color='#aaaaaa', lw=1.0, alpha=0.6, zorder=0)
handles.append(ref_h)
labels.append('ERCs = reactions')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Number of reactions  ($n_r$)', fontsize=_FL)
ax.set_ylabel('Number of ERCs', fontsize=_FL)
ax.set_title(
    f'ERC count growth vs reaction network size  ({len(df)} networks)\n'
    r'log–log axes;  lines: power-law fit  $y = A \cdot n_r^{\,\alpha}$',
    fontsize=_FT, fontweight='bold')
ax.tick_params(labelsize=_FCK)
ax.legend(handles=handles, labels=labels,
          fontsize=_FLG, loc='upper left', framealpha=0.88,
          borderpad=0.7, labelspacing=0.55)
ax.grid(True, alpha=0.20, which='both')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig(OUT_PNG, dpi=150, bbox_inches='tight')
print(f"\nSaved: {os.path.abspath(OUT_PNG)}")
