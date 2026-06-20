#!/usr/bin/env python3
"""
script_erc_growth_tendency.py
=====================================================================
Analyses how synergy and complementarity pair counts grow with the
number of ERCs, using:
  1. Sampling-quality assessment (networks per ERC bin)
  2. Power-law fit: count ~ A * n^alpha (OLS on log-log, y > 0 only)
  3. Numpy-based LOWESS smooth (no external dependency)
  4. Binned medians ± IQR  (solid = n>=MIN_N_RELIABLE, hollow = sparse)

Output:
  outputs/growth_tendency/growth_tendency.png
  (console: power-law exponents + sampling report)
"""

import os, sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import minimize

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

SYN_CSV  = os.path.join(_SCRIPT_DIR, '..', 'outputs', 'synergy_stats',
                        'synergy_stats.csv')
COMP_CSV = os.path.join(_SCRIPT_DIR, '..', 'outputs', 'complementarity_stats',
                        'complementarity_stats.csv')
OUT_DIR  = os.path.join(_SCRIPT_DIR, '..', 'outputs', 'growth_tendency')
os.makedirs(OUT_DIR, exist_ok=True)

MIN_N_RELIABLE = 15   # networks per bin for a reliable median

# ── Source classification (determines scatter marker) ─────────────────────────
# BioModels files: BIOMD*.txt  → circle 'o'
# BiGG files:      bigg_*.txt  → star   '*'
# Other:           anything else → square 's'
SOURCE_MARKER = {'biomodels': 'o', 'bigg': '*', 'other': 's'}
SOURCE_SIZE   = {'biomodels': 10,  'bigg': 22,  'other': 12}
SOURCE_LABEL  = {'biomodels': 'BioModels', 'bigg': 'BiGG', 'other': 'Other'}


def _get_source(fname):
    fn = str(fname).lower()
    if fn.startswith('bigg_'):
        return 'bigg'
    if fn.startswith('biomd'):
        return 'biomodels'
    return 'other'


# ── Load & filter ─────────────────────────────────────────────────────────────
s  = pd.read_csv(SYN_CSV)
c  = pd.read_csv(COMP_CSV)

# At least 4 ERCs so C(n,2) >= 10
sn = s[s['n_ercs'] >= 2].copy()
cn = c[c['n_ercs'] >= 2].copy()

# ── ERC bins ──────────────────────────────────────────────────────────────────
BIN_EDGES  = [2, 5, 10, 15, 25, 50, 10_000]
BIN_LABELS = ['2–4', '5–9', '10–14', '15–24', '25–49', '50+']
BIN_COLORS = ['#3498DB', '#27AE60', '#E67E22', '#8E44AD', '#E74C3C', '#2C3E50']

def assign_bins(df):
    df = df.copy()
    df['erc_bin'] = pd.cut(df['n_ercs'], bins=BIN_EDGES,
                           labels=BIN_LABELS, right=False)
    return df

sn = assign_bins(sn)
cn = assign_bins(cn)

# ── Joint classification (both / syn_only / comp_only / neither) ─────────────
# Merge on 'file' — both CSVs cover the same BioModels set
_merged = sn[['file', 'erc_bin', 'n_basic_pairs']].merge(
    cn[['file', 'n_complementary_pairs']],
    on='file', how='left'
)
_merged['n_complementary_pairs'] = _merged['n_complementary_pairs'].fillna(0)
_merged['has_syn']  = _merged['n_basic_pairs'] > 0
_merged['has_comp'] = _merged['n_complementary_pairs'] > 0

# Category order: bottom → top in stacked bar
CAT_ORDER  = ['both', 'comp_only', 'syn_only', 'neither']
CAT_COLORS = {'both': '#8E44AD', 'comp_only': '#3498DB',
              'syn_only': '#E74C3C', 'neither': '#BDC3C7'}
CAT_LABELS = {'both': 'Both (syn & comp)', 'comp_only': 'Complementarity only',
              'syn_only': 'Synergy only',    'neither': 'Neither'}

def _cat(row):
    if row['has_syn'] and row['has_comp']:  return 'both'
    if row['has_comp']:                     return 'comp_only'
    if row['has_syn']:                      return 'syn_only'
    return 'neither'

_merged['cat'] = _merged.apply(_cat, axis=1)

# Pre-compute per-bin counts for all four categories
_bin_cat_counts = {
    lbl: {cat: int((_merged[_merged['erc_bin'] == lbl]['cat'] == cat).sum())
          for cat in CAT_ORDER}
    for lbl in BIN_LABELS
}

# ── LOWESS (numpy-only, tricube weights, log-log or log-linear) ───────────────
def _lowess(y, x, frac=0.45, n_grid=80):
    """
    Locally-weighted linear regression on the given (x, y) data.
    Returns (x_grid, y_smooth) evaluated on a uniform grid over x.
    """
    n = len(x)
    k = max(int(n * frac), 5)
    x_grid = np.linspace(x.min(), x.max(), n_grid)
    y_grid  = np.empty(n_grid)
    for i, xi in enumerate(x_grid):
        dists = np.abs(x - xi)
        idx   = np.argsort(dists)[:k]
        dk    = dists[idx[-1]] + 1e-12
        w     = np.clip((1 - (dists[idx] / dk)**3)**3, 0, None)
        xk, yk = x[idx], y[idx]
        sw = w.sum()
        xm = np.dot(w, xk) / sw
        ym = np.dot(w, yk) / sw
        sxx = np.dot(w, (xk - xm)**2)
        sxy = np.dot(w, (xk - xm) * (yk - ym))
        if sxx < 1e-12:
            y_grid[i] = ym
        else:
            b = sxy / sxx
            y_grid[i] = ym + b * (xi - xm)
    return x_grid, y_grid


def lowess_loglog(x_arr, y_arr, frac=0.45):
    """LOWESS on log10(x) vs log10(y), returns (x_orig_scale, y_orig_scale)."""
    mask = (y_arr > 0) & (x_arr > 0)
    if mask.sum() < 10:
        return None, None
    lx = np.log10(x_arr[mask].astype(float))
    ly = np.log10(y_arr[mask].astype(float))
    gx, gy = _lowess(ly, lx, frac=frac)
    return 10**gx, 10**gy


def lowess_log_lin(x_arr, y_arr, frac=0.45):
    """LOWESS on log10(x) vs y (linear), returns (x_orig_scale, y_smooth)."""
    valid = np.isfinite(y_arr) & (x_arr > 0) & ~np.isnan(y_arr)
    if valid.sum() < 10:
        return None, None
    lx = np.log10(x_arr[valid].astype(float))
    ly = y_arr[valid].astype(float)
    gx, gy = _lowess(ly, lx, frac=frac)
    return 10**gx, gy

# ── Power-law fits ────────────────────────────────────────────────────────────
def powerlaw_fit(x_arr, y_arr):
    """OLS on log10-log10, y > 0 only (biased — zero networks excluded)."""
    mask = (y_arr > 0) & (x_arr > 0)
    if mask.sum() < 10:
        return None
    lx = np.log10(x_arr[mask].astype(float))
    ly = np.log10(y_arr[mask].astype(float))
    slope, intercept, r, _, se = stats.linregress(lx, ly)
    return {'alpha': slope, 'logA': intercept, 'R2': r**2,
            'n': int(mask.sum()), 'se_alpha': se}


def poisson_fit(x_arr, y_arr):
    """
    Poisson GLM: E[y] = A * n^alpha — MLE on ALL observations (zeros included).

    This is the methodologically correct estimator for count-data power laws:
    it avoids the upward-bias from conditioning on y > 0, and correctly weights
    zero-count small networks against high-count large networks.

    Returns dict with alpha, A, logA (base-10 log of A), ok; or None.
    """
    x = np.log(x_arr.astype(float))   # natural log for optimisation
    y = y_arr.astype(float)

    def nll(params):
        b0, b1 = params
        log_mu = np.clip(b0 + b1 * x, -700, 700)
        mu = np.exp(log_mu)
        # 0 * log(mu) = 0 by convention even when mu→0
        terms = np.where(y > 0, y * log_mu, 0.0)
        return float(np.sum(mu - terms))

    mask = y > 0
    if mask.sum() < 10:
        return None
    # Initialise from OLS on non-zeros
    b1_0, b0_0 = np.polyfit(x[mask], np.log(y[mask]), 1)
    res = minimize(nll, [b0_0, b1_0], method='L-BFGS-B',
                   options={'maxiter': 5000, 'ftol': 1e-14, 'gtol': 1e-8})
    b0, b1 = res.x
    A = float(np.exp(b0))
    return {'alpha': b1, 'A': A, 'logA': b0 / np.log(10), 'ok': res.success}

# ── Binned medians ─────────────────────────────────────────────────────────────
def binned_medians(df, col):
    rows = []
    for lbl in BIN_LABELS:
        sub = df.loc[df['erc_bin'] == lbl, col].dropna()
        ercs = df.loc[df['erc_bin'] == lbl, 'n_ercs']
        n = len(sub)
        if n == 0:
            continue
        rows.append({
            'bin': lbl, 'n': n,
            'x_med': float(ercs.median()),
            'med': float(sub.median()),
            'q25': float(sub.quantile(0.25)),
            'q75': float(sub.quantile(0.75)),
            'reliable': n >= MIN_N_RELIABLE,
        })
    return pd.DataFrame(rows)

# ── Sampling quality report ────────────────────────────────────────────────────
print("=" * 60)
print("SAMPLING QUALITY  (n_ercs >= 2)")
print("=" * 60)
print(f"{'Bin':>7}  {'N':>5}  {'N_syn>0':>8}  {'N_comp>0':>9}  {'Reliable?':>10}")
for lbl in BIN_LABELS:
    ns  = (sn['erc_bin'] == lbl).sum()
    nsb = int((sn[sn['erc_bin'] == lbl]['n_basic_pairs'] > 0).sum())
    ncb = int((cn[cn['erc_bin'] == lbl]['n_complementary_pairs'] > 0).sum())
    flag = 'YES' if ns >= MIN_N_RELIABLE else 'sparse'
    print(f"  {lbl:>6}  {ns:5d}  {nsb:8d}  {ncb:9d}  {flag:>10}")
print()

_CASES = [
    ('syn_basic',  sn, 'n_basic_pairs',            'Synergy (basic/maximal)'),
    ('syn_fund',   sn, 'n_fundamental_pairs',       'Synergy (fundamental)'),
    ('cmp_comp',   cn, 'n_complementary_pairs',     'Complementarity (all)'),
    ('cmp_pure',   cn, 'n_pure_complementary_pairs','Complementarity (pure)'),
    ('cmp_fund',   cn, 'n_fundamental_edges',       'Complementarity (fund)'),
]

# 3-synergy columns — only present after script_erc_synergy_stats has been re-run
_HAS_3SYN = 'n_basic_3syn_triples' in sn.columns
if _HAS_3SYN:
    sn3 = sn[sn['n_basic_3syn_triples'] >= 0].copy()
    _CASES += [
        ('syn3_basic', sn3, 'n_basic_3syn_triples',       'Ternary synergy (basic)'),
        ('syn3_fund',  sn3, 'n_fundamental_3syn_triples', 'Ternary synergy (fund)'),
    ]
else:
    sn3 = None

# Poisson GLM fits (primary — includes y=0 networks, no selection bias)
print("POWER-LAW FITS   Poisson GLM  y ~ A·n^alpha  (ALL networks, zeros included)")
print("=" * 75)
print("  Note: alpha>2 means synergy DENSITY increases with n in BioModels")
print("        (curated database is biased toward complex, integrated networks)")
print("=" * 75)
pfits = {}        # tag  → poisson fit dict
pfits_by_col = {} # col  → poisson fit dict  (for figure drawing)
for tag, df, col, label in _CASES:
    f = poisson_fit(df['n_ercs'].values, df[col].values)
    pfits[tag] = f
    pfits_by_col[col] = f
    if f:
        sign = '+' if f['alpha'] >= 2 else ''
        print(f"  {label:<38} alpha={f['alpha']:+.2f}  "
              f"(2 {'+' if f['alpha']>=2 else ''}{f['alpha']-2:+.2f})  "
              f"converged={f['ok']}")
    else:
        print(f"  {label:<38} -- insufficient data --")
print()

# OLS-on-nonzeros (legacy reference — excluded zeros cause upward bias)
print("REFERENCE  OLS on log-log  (y>0 only — biased upward for small-n types)")
print("=" * 75)
fits = {}
for tag, df, col, label in _CASES:
    f = powerlaw_fit(df['n_ercs'].values, df[col].values)
    fits[tag] = f
    if f:
        ci = 1.96 * f['se_alpha']
        print(f"  {label:<38} alpha={f['alpha']:.2f} +/- {ci:.2f}  "
              f"R2={f['R2']:.2f}  n={f['n']}")
    else:
        print(f"  {label:<38} -- insufficient data --")
print()
print("  Quadratic growth (alpha=2.00) <=> constant density = same fraction of C(n,2)")
print()

# ── FIGURE ───────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(17, 10))
gs  = fig.add_gridspec(2, 3, hspace=0.38, wspace=0.30,
                        left=0.06, right=0.97, top=0.91, bottom=0.08)

ax_samp = fig.add_subplot(gs[0, 0])
ax_syn  = fig.add_subplot(gs[0, 1])
ax_cmp  = fig.add_subplot(gs[0, 2])
ax_leg  = fig.add_subplot(gs[1, 0])
ax_sratio = fig.add_subplot(gs[1, 1])
ax_cratio = fig.add_subplot(gs[1, 2])
ax_leg.axis('off')

# ── Panel A: stacked sampling quality ────────────────────────────────────────
x_pos  = np.arange(len(BIN_LABELS))
bottoms = np.zeros(len(BIN_LABELS))
for cat in CAT_ORDER:
    heights = np.array([_bin_cat_counts[lbl][cat] for lbl in BIN_LABELS],
                       dtype=float)
    ax_samp.bar(x_pos, heights, bottom=bottoms,
                color=CAT_COLORS[cat], alpha=0.88,
                edgecolor='#333333', lw=0.5,
                label=CAT_LABELS[cat])
    # annotate each non-zero segment with its count (centred inside segment)
    for xi, h, bot in zip(x_pos, heights, bottoms):
        if h >= 8:   # only label segments tall enough to read
            ax_samp.text(xi, bot + h / 2, str(int(h)),
                         ha='center', va='center',
                         fontsize=8, fontweight='bold', color='white')
    bottoms += heights

# Total-count label above each bar
bin_ns = [int((sn['erc_bin'] == lbl).sum()) for lbl in BIN_LABELS]
for xi, n in zip(x_pos, bin_ns):
    ax_samp.text(xi, n + 2, str(n),
                 ha='center', va='bottom', fontsize=9, fontweight='bold')

ax_samp.axhline(MIN_N_RELIABLE, color='darkred', ls='--', lw=1.2,
                label=f'min reliable ({MIN_N_RELIABLE})', zorder=5)
ax_samp.set_xticks(x_pos)
ax_samp.set_xticklabels(BIN_LABELS, fontsize=9)
ax_samp.set_xlabel('Number of ERCs (bin)', fontsize=11)
ax_samp.set_ylabel('Networks in bin', fontsize=11)
ax_samp.set_title('Sampling quality & synergy/complementarity coverage',
                  fontsize=10.5, fontweight='bold')
ax_samp.legend(fontsize=8, loc='upper right', framealpha=0.88)
ax_samp.grid(True, axis='y', alpha=0.25)
ax_samp.set_xlim(-0.6, len(BIN_LABELS) - 0.4)

# ── Helper: draw absolute-count growth panel ──────────────────────────────────
def _scatter_by_source(ax, df, n_arr, y_arr, color, zorder=1):
    """Scatter points using per-source markers: o=BioModels, *=BiGG, s=Other."""
    sources = df['file'].map(_get_source) if 'file' in df.columns else None
    mask = y_arr > 0
    if sources is None:
        ax.scatter(n_arr[mask], y_arr[mask], c=color, s=10, alpha=0.22, zorder=zorder)
        return
    for src, mkr in SOURCE_MARKER.items():
        src_mask = mask & (sources == src)
        if src_mask.any():
            ax.scatter(n_arr[src_mask], y_arr[src_mask],
                       c=color, marker=mkr, s=SOURCE_SIZE[src],
                       alpha=0.28, zorder=zorder)


def draw_abs_panel(ax, df, series, title, show_c3_max=False):
    """
    series: list of (col, color, label, marker)
    Scatter markers encode data source: o=BioModels, *=BiGG, s=Other.
    Uses Poisson GLM exponents (pfits_by_col) for legend labels and fit lines.
    show_c3_max: also draw C(n,3) dashed line for ternary reference.
    """
    n_arr = df['n_ercs'].values.astype(float)

    # Theoretical max C(n,2)
    x_th = np.logspace(np.log10(2), np.log10(n_arr.max()), 250)
    ax.plot(x_th, x_th * (x_th - 1) / 2,
            'k--', lw=1.2, alpha=0.45, zorder=0, label='$\\binom{n}{2}$ (max)')
    if show_c3_max:
        ax.plot(x_th, x_th * (x_th - 1) * (x_th - 2) / 6,
                'k:', lw=1.0, alpha=0.35, zorder=0, label='$\\binom{n}{3}$ (max)')

    for col, color, label, marker in series:
        y_raw = df[col].values.astype(float)
        # exclude -1 sentinel (3-synergy not computed for this network)
        valid  = y_raw >= 0
        n_arr2 = n_arr[valid]
        y_arr  = y_raw[valid]
        df2_   = df[valid].copy()

        fit = pfits_by_col.get(col) or powerlaw_fit(n_arr2, y_arr)

        # Scatter by source
        _scatter_by_source(ax, df2_, n_arr2, y_arr, color, zorder=1)

        # LOWESS
        xs, ys = lowess_loglog(n_arr2, y_arr, frac=0.45)
        if xs is not None:
            ax.plot(xs, ys, color=color, lw=2.8, zorder=3)

        # Power-law fit line
        if fit and n_arr2.max() > 0:
            x_f = np.logspace(np.log10(2), np.log10(n_arr2.max()), 150)
            y_f = fit['A'] * x_f**fit['alpha']
            ax.plot(x_f, y_f, color=color, lw=1.1, ls=':', alpha=0.9, zorder=2)
            lbl_txt = f"{label}  ($\\alpha$={fit['alpha']:.2f})"
        else:
            lbl_txt = label

        # Binned medians
        df3 = df2_.copy()
        df3['_y'] = np.where(y_arr > 0, y_arr, np.nan)
        df3['n_ercs'] = n_arr2
        bm = binned_medians(df3, '_y')
        for _, row in bm.iterrows():
            if np.isnan(row['med']) or row['med'] <= 0:
                continue
            kw = dict(fmt=marker, color=color, ms=8, lw=1.8,
                      capsize=4, zorder=5, elinewidth=1.5)
            yerr_lo = max(row['med'] - row['q25'], 0)
            yerr_hi = max(row['q75'] - row['med'], 0)
            ax.errorbar(row['x_med'], row['med'],
                        yerr=[[yerr_lo], [yerr_hi]],
                        mfc=color if row['reliable'] else 'white', **kw)

        ax.scatter([], [], c=color, marker=marker, s=60, label=lbl_txt)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Number of ERCs  ($n$)', fontsize=11)
    ax.set_ylabel('Count', fontsize=11)
    ax.set_title(title, fontsize=11, fontweight='bold')
    ax.legend(fontsize=8.5, loc='upper left', framealpha=0.85)
    ax.grid(True, alpha=0.18, which='both')

_syn_series = [
    ('n_basic_pairs',       '#E74C3C', 'Binary basic / maximal', 'o'),
    ('n_fundamental_pairs', '#27AE60', 'Binary fundamental',     's'),
]
if _HAS_3SYN and sn3 is not None and len(sn3) > 0:
    _syn_series += [
        ('n_basic_3syn_triples',       '#9B59B6', 'Ternary basic',       'o'),
        ('n_fundamental_3syn_triples', '#1ABC9C', 'Ternary fundamental', 's'),
    ]
draw_abs_panel(ax_syn, sn, _syn_series,
               title='Synergy: count growth vs ERCs  (log–log)',
               show_c3_max=_HAS_3SYN)

draw_abs_panel(ax_cmp, cn, [
    ('n_complementary_pairs',       '#E74C3C', 'Complementary',        'o'),
    ('n_pure_complementary_pairs',  '#E67E22', 'Purely complementary', '^'),
    ('n_fundamental_edges',         '#27AE60', 'Fundamental',          's'),
], title='Complementarity: count growth vs ERCs  (log–log)')

# ── Helper: draw ratio panel (count / C(n,2) vs n_ERCs) ───────────────────────
def draw_ratio_panel(ax, df, series, title):
    """
    series: list of (ratio_col, color, label, marker)
    ratio_col in [0,1]; -1 sentinel means not computed (excluded).
    Scatter markers encode data source: o=BioModels, *=BiGG, s=Other.
    """
    n_arr = df['n_ercs'].values.astype(float)

    for ratio_col, color, label, marker in series:
        r_raw = df[ratio_col].values.astype(float)
        # Exclude -1 sentinels (3-synergy not computed)
        valid  = (r_raw >= 0) & np.isfinite(r_raw) & (n_arr > 0)
        n_v    = n_arr[valid]
        r_arr  = r_raw[valid]
        df_v   = df[valid].copy()

        # Scatter by source
        sources = df_v['file'].map(_get_source) if 'file' in df_v.columns else None
        if sources is None:
            ax.scatter(n_v, r_arr, c=color, s=10, alpha=0.22, zorder=1)
        else:
            for src, mkr in SOURCE_MARKER.items():
                src_mask = sources == src
                if src_mask.any():
                    ax.scatter(n_v[src_mask.values], r_arr[src_mask.values],
                               c=color, marker=mkr, s=SOURCE_SIZE[src],
                               alpha=0.28, zorder=1)

        # LOWESS on (log n, ratio)
        xs, ys = lowess_log_lin(n_v, r_arr, frac=0.50)
        if xs is not None:
            ax.plot(xs, ys, color=color, lw=2.8, zorder=3, label=label)

        # Binned medians
        df_v2 = df_v.copy()
        df_v2['n_ercs'] = n_v
        bm = binned_medians(df_v2, ratio_col)
        for _, row in bm.iterrows():
            kw = dict(fmt=marker, color=color, ms=8, lw=1.8,
                      capsize=4, zorder=5, elinewidth=1.5)
            yerr_lo = max(row['med'] - row['q25'], 0)
            yerr_hi = max(row['q75'] - row['med'], 0)
            ax.errorbar(row['x_med'], row['med'],
                        yerr=[[yerr_lo], [yerr_hi]],
                        mfc=color if row['reliable'] else 'white', **kw)

    ax.set_xscale('log')
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel('Number of ERCs  ($n$)', fontsize=11)
    ax.set_ylabel('Pairs / $\\binom{n}{2}$', fontsize=11)
    ax.set_title(title, fontsize=11, fontweight='bold')
    ax.legend(fontsize=8.5, loc='upper right', framealpha=0.85)
    ax.grid(True, alpha=0.2)

_sratio_series = [
    ('ratio_basic',       '#E74C3C', 'Binary basic / maximal', 'o'),
    ('ratio_fundamental', '#27AE60', 'Binary fundamental',     's'),
]
if _HAS_3SYN and sn3 is not None and len(sn3) > 0:
    _sratio_series += [
        ('ratio_3syn_basic',       '#9B59B6', 'Ternary basic / C(n,3)',       'o'),
        ('ratio_3syn_fundamental', '#1ABC9C', 'Ternary fundamental / C(n,3)', 's'),
    ]
draw_ratio_panel(ax_sratio, sn, _sratio_series,
                 title='Synergy fraction vs ERCs')

draw_ratio_panel(ax_cratio, cn, [
    ('ratio_complementary', '#E74C3C', 'Complementary',        'o'),
    ('ratio_pure',          '#E67E22', 'Purely complementary', '^'),
    ('ratio_fundamental',   '#27AE60', 'Fundamental',          's'),
], title='Complementarity fraction of $\\binom{n}{2}$  vs ERCs')

# ── Legend panel ──────────────────────────────────────────────────────────────
from matplotlib.lines import Line2D

legend_items = [
    Line2D([0],[0], color='gray', lw=2.8,
           label='LOWESS smooth'),
    Line2D([0],[0], color='gray', lw=1.1, ls=':',
           label='Power-law fit  ($y \\sim n^{\\alpha}$,  Poisson GLM)'),
    Line2D([0],[0], color='k', lw=1.2, ls='--', alpha=0.5,
           label='$\\binom{n}{2}$ (theoretical max)'),
    Line2D([0],[0], marker='o', color='gray', ms=8, lw=0,
           mfc='gray', label=f'Bin median ± IQR  (n $\\geq$ {MIN_N_RELIABLE})'),
    Line2D([0],[0], marker='o', color='gray', ms=8, lw=0,
           mfc='white', label=f'Bin median ± IQR  (n < {MIN_N_RELIABLE}, sparse)'),
    # Data source markers
    Line2D([0],[0], marker='o', color='gray', ms=7, lw=0, alpha=0.7,
           label='BioModels network'),
    Line2D([0],[0], marker='*', color='gray', ms=10, lw=0, alpha=0.7,
           label='BiGG network'),
    Line2D([0],[0], marker='s', color='gray', ms=7, lw=0, alpha=0.7,
           label='Other network'),
]
ax_leg.legend(handles=legend_items, fontsize=10.5, loc='center',
              title='Symbol guide', title_fontsize=11, framealpha=0.9,
              borderpad=1.0, labelspacing=0.9)

fig.suptitle(
    f'Growth tendency of synergy and complementarity with ERC count\n'
    f'({len(sn)} non-trivial networks from BioModels;  '
    f'solid markers = reliable bins  $n\\geq{MIN_N_RELIABLE}$)',
    fontsize=12, fontweight='bold'
)

out_path = os.path.join(OUT_DIR, 'growth_tendency.png')
plt.savefig(out_path, dpi=150, bbox_inches='tight')
print(f"Saved: {os.path.abspath(out_path)}")
plt.show()
