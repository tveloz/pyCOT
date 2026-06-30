#!/usr/bin/env python3
"""
script_erc_forest_viz.py
========================
Growth-tendency visualisation for ERC hierarchy forest statistics.

Single figure: n_ercs vs n_reactions on log-log axes.
  Groups: BioModels (all BioMD_* sub-types combined), BiGG, All (aggregated).
  Each group: scatter + power-law fit line y = A · n^α + binned medians ± IQR.
  Overall fit across all networks shown as a dashed black line.

Output: visualizations/forest_stats/erc_growth_tendency.png
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats

# ── Font and figure size (adjust here to rescale all text) ────────────────────
FONT_SCALE = 2.0          # multiply all font sizes by this factor
_FBASE_LAB = 12           # base axis-label fontsize
_FBASE_TIT = 12           # base title fontsize
_FBASE_LEG = 9.5          # base legend fontsize
_FBASE_TCK = 10           # base tick-label fontsize
FIG_W, FIG_H = 14, 10    # figure size in inches (increase if text overlaps)

# ── Paths ──────────────────────────────────────────────────────────────────────
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

CSV_FILE = os.path.normpath(os.path.join(
    _SCRIPT_DIR, '..', 'outputs', 'forest_stats', 'forest_stats.csv'))

OUT_DIR = os.path.normpath(os.path.join(
    _SCRIPT_DIR, '..', 'visualizations', 'forest_stats'))
os.makedirs(OUT_DIR, exist_ok=True)

# ── Groups ─────────────────────────────────────────────────────────────────────
# (name, color, marker, scatter_size)
GROUPS = [
    ('BioModels', '#E74C3C', 'o', 80,  0),    # (name, color, marker, s, linewidth)
    ('BiGG',      '#2C3E50', '*', 500, 1.5),  # stars need large s and linewidth to be visible
]

def assign_group(dataset):
    if str(dataset).startswith('BioMD_'):
        return 'BioModels'
    if dataset == 'BiGG':
        return 'BiGG'
    return 'Other'

# ── Load & filter ──────────────────────────────────────────────────────────────
df_all = pd.read_csv(CSV_FILE)
df = df_all[(df_all['n_ercs'] >= 4) & (df_all['n_reactions'] > 0)].copy()
df['group'] = df['dataset'].map(assign_group)

print(f"Loaded {len(df)} networks (n_ercs >= 4)")
print(df.groupby('group')['n_ercs'].count().to_string())
print()

# ── Power-law fit (OLS on log-log) ────────────────────────────────────────────
def powerlaw_fit(x_arr, y_arr):
    """OLS  y = A · x^α  on log-log (y > 0, x > 0 only)."""
    x_arr = np.asarray(x_arr, dtype=float)
    y_arr = np.asarray(y_arr, dtype=float)
    mask = (y_arr > 0) & (x_arr > 0)
    if mask.sum() < 5:
        return None
    lx = np.log10(x_arr[mask])
    ly = np.log10(y_arr[mask])
    slope, intercept, r, _, se = stats.linregress(lx, ly)
    return {'alpha': slope, 'logA': intercept, 'R2': r**2,
            'n': int(mask.sum()), 'se_alpha': se}

# ── Binned medians (log-spaced over x) ────────────────────────────────────────
def binned_medians(x_arr, y_arr, n_bins=7):
    """Returns list of dicts {x, med, q25, q75, n} using log-spaced bins."""
    x_arr = np.asarray(x_arr, dtype=float)
    y_arr = np.asarray(y_arr, dtype=float)
    mask = (x_arr > 0) & (y_arr > 0)
    x, y = x_arr[mask], y_arr[mask]
    if len(x) < 5:
        return []
    edges = np.logspace(np.log10(x.min()), np.log10(x.max() * 1.001), n_bins + 1)
    rows = []
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

# ── Figure ────────────────────────────────────────────────────────────────────
_FL  = int(round(_FBASE_LAB * FONT_SCALE))
_FT  = int(round(_FBASE_TIT * FONT_SCALE))
_FLG = int(round(_FBASE_LEG * FONT_SCALE))
_FCK = int(round(_FBASE_TCK * FONT_SCALE))

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

    # Scatter
    ax.scatter(x, y, color=color, marker=marker, s=ms,
               alpha=0.45, zorder=2, linewidths=lw)

    # Power-law fit line
    fit = powerlaw_fit(x, y)
    if fit is not None:
        A   = 10 ** fit['logA']
        x_f = np.logspace(np.log10(max(x.min(), 0.5)), np.log10(x.max()), 250)
        ax.plot(x_f, A * x_f ** fit['alpha'],
                color=color, lw=2.2, zorder=3, alpha=0.92)
        ci = 1.96 * fit['se_alpha']
        print(f"  {grp_name:<12}  alpha={fit['alpha']:.3f} +/- {ci:.3f}"
              f"  R2={fit['R2']:.3f}  n={fit['n']}")
        lbl = f"{grp_name}  ($\\alpha = {fit['alpha']:.2f}$,  n={fit['n']})"
    else:
        print(f"  {grp_name:<12}  -- sparse --")
        lbl = f"{grp_name}  (n={len(sub)}, sparse)"

    # Binned medians ± IQR
    bm = binned_medians(x, y, n_bins=6)
    for row in bm:
        yerr_lo = max(row['med'] - row['q25'], 0)
        yerr_hi = max(row['q75'] - row['med'], 0)
        ax.errorbar(row['x'], row['med'],
                    yerr=[[yerr_lo], [yerr_hi]],
                    fmt=marker, color=color, ms=8, lw=1.6,
                    capsize=4, zorder=5, elinewidth=1.4,
                    mfc='white', mec=color, mew=1.5)

    proxy = Line2D([0], [0], color=color, lw=2.2,
                   marker=marker, ms=6, mfc=color, mec=color)
    handles.append(proxy)
    labels.append(lbl)

# ── Overall (all groups) power-law fit ────────────────────────────────────────
x_all = df['n_reactions'].values.astype(float)
y_all = df['n_ercs'].values.astype(float)
fit_all = powerlaw_fit(x_all, y_all)

if fit_all is not None:
    A_all = 10 ** fit_all['logA']
    x_f_all = np.logspace(
        np.log10(max(x_all[x_all > 0].min(), 0.5)),
        np.log10(x_all.max()), 250)
    line_all, = ax.plot(x_f_all, A_all * x_f_all ** fit_all['alpha'],
                        color='#111111', lw=2.5, ls='--', zorder=6, alpha=0.88)
    handles.append(line_all)
    labels.append(
        f"All  ($\\alpha = {fit_all['alpha']:.2f}$,  n={fit_all['n']})")
    ci_all = 1.96 * fit_all['se_alpha']
    print(f"  {'All':<12}  alpha={fit_all['alpha']:.3f} +/- {ci_all:.3f}"
          f"  R2={fit_all['R2']:.3f}  n={fit_all['n']}")

# ── Binned medians legend proxy ────────────────────────────────────────────────
bm_proxy = Line2D([0], [0], marker='o', color='#555555', ms=7, lw=0,
                  mfc='white', mec='#555555', mew=1.4)
handles.append(bm_proxy)
labels.append('Median ± IQR (binned, hollow)')

# ── Reference line  n_ercs = n_reactions ──────────────────────────────────────
xr = np.array([max(x_all[x_all > 0].min(), 0.5), x_all.max()])
ref_h, = ax.plot(xr, xr, ':', color='#aaaaaa', lw=1.0, alpha=0.6, zorder=0)
handles.append(ref_h)
labels.append('ERCs = reactions')

# ── Axes & legend ─────────────────────────────────────────────────────────────
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
out_path = os.path.join(OUT_DIR, 'erc_growth_tendency.png')
plt.savefig(out_path, dpi=150, bbox_inches='tight')
print(f"\nSaved: {os.path.abspath(out_path)}")
