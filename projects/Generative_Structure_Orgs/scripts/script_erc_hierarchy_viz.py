#!/usr/bin/env python3
"""
script_erc_hierarchy_viz.py
===========================
Visualises how ERC hierarchy depth and branching grow with the number of ERCs.

Two panels (log-log axes):
  Left  — mean tree depth  (per-tree max depth averaged across trees)
  Right — mean branching   (mean out-degree among non-leaf ERC nodes,
                            i.e. average number of children per internal node)

Groups: BioModels (all BioMD_* combined, red), BiGG (dark), All (black dashed).
Power-law fits y = A * n^alpha by OLS on log-log data (y > 0 only).
Binned medians +/- IQR shown as hollow markers.

Output: visualizations/forest_stats/erc_hierarchy_structure.png
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
_FBASE_TIT = 11           # base title fontsize
_FBASE_LEG = 9            # base legend fontsize
_FBASE_SUP = 13           # base suptitle fontsize
FIG_W, FIG_H = 22, 10    # figure size in inches (increase if text overlaps)

# ── Paths ──────────────────────────────────────────────────────────────────────
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

CSV_FILE = os.path.normpath(os.path.join(
    _SCRIPT_DIR, '..', 'outputs', 'forest_stats', 'forest_stats.csv'))

OUT_DIR = os.path.normpath(os.path.join(
    _SCRIPT_DIR, '..', 'visualizations', 'forest_stats'))
os.makedirs(OUT_DIR, exist_ok=True)

# ── Groups ─────────────────────────────────────────────────────────────────────
GROUPS = [
    ('BioModels', '#E74C3C', 'o', 14),
    ('BiGG',      '#2C3E50', '*', 35),
]

def assign_group(dataset):
    if str(dataset).startswith('BioMD_'):
        return 'BioModels'
    if dataset == 'BiGG':
        return 'BiGG'
    return 'Other'

# ── Load & filter ──────────────────────────────────────────────────────────────
df_all = pd.read_csv(CSV_FILE)
df = df_all[df_all['n_ercs'] >= 4].copy()
df['group'] = df['dataset'].map(assign_group)

print(f"Loaded {len(df)} networks (n_ercs >= 4)")
print(df.groupby('group')['n_ercs'].count().to_string())

# ── Power-law fit (OLS on log-log, y > 0 only) ───────────────────────────────
def powerlaw_fit(x_arr, y_arr):
    x_arr = np.asarray(x_arr, dtype=float)
    y_arr = np.asarray(y_arr, dtype=float)
    mask = (y_arr > 0) & (x_arr > 0) & np.isfinite(y_arr) & np.isfinite(x_arr)
    if mask.sum() < 5:
        return None
    lx = np.log10(x_arr[mask])
    ly = np.log10(y_arr[mask])
    slope, intercept, r, _, se = stats.linregress(lx, ly)
    return {'alpha': slope, 'logA': intercept, 'R2': r**2,
            'n': int(mask.sum()), 'se_alpha': se}

# ── Binned medians (log-spaced bins over x) ───────────────────────────────────
def binned_medians(x_arr, y_arr, n_bins=6):
    x_arr = np.asarray(x_arr, dtype=float)
    y_arr = np.asarray(y_arr, dtype=float)
    mask = (x_arr > 0) & (y_arr > 0) & np.isfinite(x_arr) & np.isfinite(y_arr)
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
        })
    return rows


def draw_panel(ax, x_col, y_col, title, xlabel, ylabel, tag=''):
    handles, labels = [], []
    print(f"\n[{tag}]  {y_col} ~ A * {x_col}^alpha  (OLS log-log)")

    for grp_name, color, marker, ms in GROUPS:
        sub = df[df['group'] == grp_name].dropna(subset=[x_col, y_col])
        if len(sub) == 0:
            continue

        x = sub[x_col].values.astype(float)
        y = sub[y_col].values.astype(float)

        # Scatter (all points, including y=0, for context)
        ax.scatter(x, y, color=color, marker=marker, s=ms,
                   alpha=0.30, zorder=2, linewidths=0)

        fit = powerlaw_fit(x, y)
        if fit is not None:
            A = 10 ** fit['logA']
            x_pos = x[x > 0]
            x_f = np.logspace(np.log10(x_pos.min()), np.log10(x_pos.max()), 300)
            ax.plot(x_f, A * x_f ** fit['alpha'],
                    color=color, lw=2.2, zorder=3, alpha=0.92)
            ci = 1.96 * fit['se_alpha']
            print(f"  {grp_name:<12}  alpha={fit['alpha']:.3f} +/- {ci:.3f}"
                  f"  R2={fit['R2']:.3f}  n={fit['n']}")
            lbl = f"{grp_name}  ($\\alpha={fit['alpha']:.2f}$, n={fit['n']})"
        else:
            print(f"  {grp_name:<12}  -- sparse --")
            lbl = f"{grp_name}  (sparse)"

        # Binned medians +/- IQR as hollow markers
        bm = binned_medians(x, y)
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

    # Overall (all groups) fit
    x_all = df[x_col].values.astype(float)
    y_all = df.dropna(subset=[y_col])[y_col].values.astype(float)
    x_all2 = df.dropna(subset=[y_col])[x_col].values.astype(float)
    fit_all = powerlaw_fit(x_all2, y_all)
    if fit_all is not None:
        A_all = 10 ** fit_all['logA']
        xp = x_all2[x_all2 > 0]
        x_f = np.logspace(np.log10(xp.min()), np.log10(xp.max()), 300)
        line_all, = ax.plot(x_f, A_all * x_f ** fit_all['alpha'],
                            color='#111111', lw=2.5, ls='--', zorder=6, alpha=0.88)
        ci = 1.96 * fit_all['se_alpha']
        print(f"  {'All':<12}  alpha={fit_all['alpha']:.3f} +/- {ci:.3f}"
              f"  R2={fit_all['R2']:.3f}  n={fit_all['n']}")
        handles.append(line_all)
        labels.append(f"All  ($\\alpha={fit_all['alpha']:.2f}$, n={fit_all['n']})")

    # Legend proxy for binned medians
    bm_proxy = Line2D([0], [0], marker='o', color='#555555', ms=7, lw=0,
                      mfc='white', mec='#555555', mew=1.4)
    handles.append(bm_proxy)
    labels.append('Median ± IQR (binned, hollow)')

    _FL  = int(round(_FBASE_LAB * FONT_SCALE))
    _FT  = int(round(_FBASE_TIT * FONT_SCALE))
    _FLG = int(round(_FBASE_LEG * FONT_SCALE))
    _FCK = int(round(10 * FONT_SCALE))
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(xlabel, fontsize=_FL)
    ax.set_ylabel(ylabel, fontsize=_FL)
    ax.set_title(title, fontsize=_FT, fontweight='bold')
    ax.tick_params(labelsize=_FCK)
    ax.legend(handles=handles, labels=labels,
              fontsize=_FLG, loc='upper left', framealpha=0.88,
              borderpad=0.7, labelspacing=0.55)
    ax.grid(True, alpha=0.20, which='both')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


# ── Figure ────────────────────────────────────────────────────────────────────
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(FIG_W, FIG_H))

draw_panel(
    ax1,
    x_col='n_ercs',
    y_col='mean_tree_depth',
    title='ERC hierarchy depth vs ERC count\n'
          r'log$-$log;  fit $y = A \cdot n^{\,\alpha}$',
    xlabel='Number of ERCs  ($n$)',
    ylabel='Mean tree depth',
    tag='Depth',
)

draw_panel(
    ax2,
    x_col='n_ercs',
    y_col='mean_out_degree',
    title='ERC hierarchy branching vs ERC count\n'
          r'log$-$log;  fit $y = A \cdot n^{\,\alpha}$',
    xlabel='Number of ERCs  ($n$)',
    ylabel='Mean out-degree of non-leaf nodes',
    tag='Branching',
)

fig.suptitle(
    f'ERC hierarchy structure: depth and branching  ({len(df)} networks, '
    r'$n_{\mathrm{ERC}} \geq 4$)',
    fontsize=int(round(_FBASE_SUP * FONT_SCALE)), fontweight='bold',
)
plt.tight_layout()

out_path = os.path.join(OUT_DIR, 'erc_hierarchy_structure.png')
plt.savefig(out_path, dpi=150, bbox_inches='tight')
print(f"\nSaved: {os.path.abspath(out_path)}")
