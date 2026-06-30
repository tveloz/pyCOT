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
import matplotlib.ticker
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

# ── Font size (adjust here to rescale all text in the figure) ─────────────────
FONT_SCALE = 1.0   # multiply all font sizes by this factor

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
sn = s[s['n_ercs'] >= 4].copy()
cn = c[c['n_ercs'] >= 4].copy()

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
# Layout: if ternary data exist, split synergy into a 4-column figure;
# otherwise keep 3 columns.  Binary and ternary always get separate panels.
_has_ternary = _HAS_3SYN and sn3 is not None and len(sn3) > 0

if _has_ternary:
    # 4-col layout; figsize scaled up to accommodate larger text
    _b = [16, 13, 10, 18, 10]   # base: FLAB, FTIT, FLEG, FSUP, FTCK
    _FLAB, _FTIT, _FLEG, _FSUP, _FTCK = [int(round(v * FONT_SCALE)) for v in _b]
    fig = plt.figure(figsize=(22, 14))
    gs  = fig.add_gridspec(2, 4, hspace=0.65, wspace=0.45,
                            left=0.05, right=0.97, top=0.89, bottom=0.08)
    ax_samp    = fig.add_subplot(gs[0, 0])
    ax_syn_bin = fig.add_subplot(gs[0, 1])   # binary synergy absolute
    ax_syn_ter = fig.add_subplot(gs[0, 2])   # ternary synergy absolute
    ax_cmp     = fig.add_subplot(gs[0, 3])   # complementarity absolute
    ax_leg     = fig.add_subplot(gs[1, 0])
    ax_sratio  = fig.add_subplot(gs[1, 1])   # binary synergy ratio
    ax_sratio3 = fig.add_subplot(gs[1, 2])   # ternary synergy ratio
    ax_cratio  = fig.add_subplot(gs[1, 3])   # complementarity ratio
else:
    # 3-col layout; figsize scaled up to accommodate larger text
    _b = [12, 10, 8, 15, 9]     # base: FLAB, FTIT, FLEG, FSUP, FTCK
    _FLAB, _FTIT, _FLEG, _FSUP, _FTCK = [int(round(v * FONT_SCALE)) for v in _b]
    fig = plt.figure(figsize=(16, 12))
    gs  = fig.add_gridspec(2, 3, hspace=0.62, wspace=0.45,
                            left=0.06, right=0.97, top=0.89, bottom=0.08)
    ax_samp    = fig.add_subplot(gs[0, 0])
    ax_syn_bin = fig.add_subplot(gs[0, 1])
    ax_cmp     = fig.add_subplot(gs[0, 2])
    ax_leg     = fig.add_subplot(gs[1, 0])
    ax_sratio  = fig.add_subplot(gs[1, 1])
    ax_cratio  = fig.add_subplot(gs[1, 2])
    ax_syn_ter = None
    ax_sratio3 = None

ax_leg.axis('off')

# ── Panel A: stacked sampling quality ────────────────────────────────────────
x_pos   = np.arange(len(BIN_LABELS))
bottoms = np.zeros(len(BIN_LABELS))
for cat in CAT_ORDER:
    heights = np.array([_bin_cat_counts[lbl][cat] for lbl in BIN_LABELS],
                       dtype=float)
    ax_samp.bar(x_pos, heights, bottom=bottoms,
                color=CAT_COLORS[cat], alpha=0.88,
                edgecolor='#333333', lw=0.5,
                label=CAT_LABELS[cat])
    for xi, h, bot in zip(x_pos, heights, bottoms):
        if h >= 8:
            ax_samp.text(xi, bot + h / 2, str(int(h)),
                         ha='center', va='center',
                         fontsize=_FTCK - 2, fontweight='bold', color='white')
    bottoms += heights

bin_ns = [int((sn['erc_bin'] == lbl).sum()) for lbl in BIN_LABELS]
for xi, n in zip(x_pos, bin_ns):
    ax_samp.text(xi, n + 2, str(n),
                 ha='center', va='bottom', fontsize=_FTCK, fontweight='bold')

ax_samp.axhline(MIN_N_RELIABLE, color='darkred', ls='--', lw=1.2,
                label=f'min reliable ({MIN_N_RELIABLE})', zorder=5)
ax_samp.set_xticks(x_pos)
ax_samp.set_xticklabels(BIN_LABELS, fontsize=_FTCK)
ax_samp.set_xlabel('Number of ERCs (bin)', fontsize=_FLAB)
ax_samp.set_ylabel('Networks in bin', fontsize=_FLAB)
ax_samp.set_title('Sampling quality\n& synergy/complementarity coverage',
                  fontsize=_FTIT, fontweight='bold')
ax_samp.legend(fontsize=_FLEG, loc='upper right', framealpha=0.88)
ax_samp.tick_params(labelsize=_FTCK)
ax_samp.grid(True, axis='y', alpha=0.25)
ax_samp.set_xlim(-0.6, len(BIN_LABELS) - 0.4)

# ── Helpers ───────────────────────────────────────────────────────────────────
def _scatter_by_source(ax, df, n_arr, y_arr, color, zorder=1):
    """Scatter y>0 points using per-source markers."""
    sources = df['file'].map(_get_source) if 'file' in df.columns else None
    mask = y_arr > 0
    if sources is None:
        ax.scatter(n_arr[mask], y_arr[mask], c=color, s=10, alpha=0.20, zorder=zorder)
        return
    for src, mkr in SOURCE_MARKER.items():
        src_mask = mask & (sources == src)
        if src_mask.any():
            ax.scatter(n_arr[src_mask], y_arr[src_mask],
                       c=color, marker=mkr, s=SOURCE_SIZE[src],
                       alpha=0.25, zorder=zorder)


def _fit_line(fit):
    """Return the amplitude A from either a Poisson GLM or an OLS fit dict."""
    if fit is None:
        return None
    return fit.get('A') or 10 ** fit['logA']


# ── Absolute-count panel ──────────────────────────────────────────────────────
def draw_abs_panel(ax, df, series, title, max_line='n2'):
    """
    series : list of (col, color, label, marker)
    max_line: 'n2' = C(n,2)  |  'n3' = C(n,3)  |  None = skip

    Main visual: power-law fit line  y = A · nᵅ  (straight line on log–log).
    The exponent α is printed in each legend entry.
    Binned medians ± IQR serve as empirical anchors.
    """
    n_arr = df['n_ercs'].values.astype(float)
    x_max = max(n_arr.max(), 5)
    x_th  = np.logspace(np.log10(2), np.log10(x_max), 300)

    if max_line == 'n2':
        ax.plot(x_th, x_th * (x_th - 1) / 2,
                color='#555555', lw=1.0, ls='--', alpha=0.35, zorder=0,
                label='$\\binom{n}{2}$ max')
    elif max_line == 'n3':
        ax.plot(x_th, x_th * (x_th - 1) * (x_th - 2) / 6,
                color='#555555', lw=1.0, ls='--', alpha=0.35, zorder=0,
                label='$\\binom{n}{3}$ max')

    for col, color, label, marker in series:
        y_raw  = df[col].values.astype(float)
        valid  = y_raw >= 0                   # exclude -1 sentinels
        n_v    = n_arr[valid]
        y_arr  = y_raw[valid]
        df_v   = df[valid].copy()

        # Power-law fit: Poisson GLM preferred, OLS fallback
        fit = pfits_by_col.get(col) or powerlaw_fit(n_v, y_arr)
        A   = _fit_line(fit)

        # Scatter (y > 0 only)
        _scatter_by_source(ax, df_v, n_v, y_arr, color, zorder=1)

        # Power-law line  y = A · nᵅ  — a straight line on log–log axes
        if fit is not None and A is not None:
            x_f = np.logspace(np.log10(max(n_v.min(), 1.5)), np.log10(n_v.max()), 250)
            ax.plot(x_f, A * x_f ** fit['alpha'],
                    color=color, lw=2.2, zorder=3, alpha=0.90)
            lbl_txt = f"{label}  ($\\alpha = {fit['alpha']:.2f}$)"
        else:
            lbl_txt = label

        # Binned medians ± IQR  (empirical check against the fit line)
        df3 = df_v.copy()
        df3['_y']     = np.where(y_arr > 0, y_arr, np.nan)
        df3['n_ercs'] = n_v
        bm = binned_medians(df3, '_y')
        for _, row in bm.iterrows():
            if np.isnan(row['med']) or row['med'] <= 0:
                continue
            yerr_lo = max(row['med'] - row['q25'], 0)
            yerr_hi = max(row['q75'] - row['med'], 0)
            ax.errorbar(row['x_med'], row['med'],
                        yerr=[[yerr_lo], [yerr_hi]],
                        fmt=marker, color=color, ms=8, lw=1.8,
                        capsize=4, zorder=5, elinewidth=1.5,
                        mfc=color if row['reliable'] else 'white')

        # Legend proxy (line + marker together)
        ax.plot([], [], color=color, lw=2.2, marker=marker, ms=8, label=lbl_txt)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Number of ERCs  ($n$)', fontsize=_FLAB)
    ax.set_ylabel('Count', fontsize=_FLAB)
    ax.set_title(title, fontsize=_FTIT, fontweight='bold')
    ax.legend(fontsize=_FLEG, loc='upper left', framealpha=0.82,
              borderpad=0.6, labelspacing=0.5)
    ax.tick_params(labelsize=_FTCK)
    ax.grid(True, alpha=0.18, which='both')


# ── Ratio panel (log–log y-axis) ──────────────────────────────────────────────
def draw_ratio_panel(ax, df, series, title, denom_label='$\\binom{n}{2}$'):
    """
    Log–log ratio panel.  Main visual: power-law fit line on ratio data.

    ratio ~ B · nᵝ  (straight line on log–log).
    For binary counts:   β = α_count − 2  (since C(n,2) ~ n²/2).
    For ternary counts:  β = α_count − 3  (since C(n,3) ~ n³/6).
    So β < 0 means the fraction shrinks with network size (fundamental become rarer).

    Networks with ratio = 0 are excluded from the fit and scatter (counted in label).
    """
    n_arr  = df['n_ercs'].values.astype(float)
    y_mins = []

    for ratio_col, color, label, marker in series:
        r_raw  = df[ratio_col].values.astype(float)
        valid  = (r_raw > 0) & np.isfinite(r_raw) & (n_arr > 0)
        n_zero = int((r_raw >= 0).sum()) - int(valid.sum())
        n_v    = n_arr[valid]
        r_arr  = r_raw[valid]
        df_v   = df[valid].copy()

        if len(r_arr) < 5:
            continue

        y_mins.append(r_arr.min())

        # Scatter (by source)
        sources = df_v['file'].map(_get_source) if 'file' in df_v.columns else None
        if sources is None:
            ax.scatter(n_v, r_arr, c=color, s=10, alpha=0.20, zorder=1)
        else:
            for src, mkr in SOURCE_MARKER.items():
                sm = sources == src
                if sm.any():
                    ax.scatter(n_v[sm.values], r_arr[sm.values],
                               c=color, marker=mkr, s=SOURCE_SIZE[src],
                               alpha=0.25, zorder=1)

        # OLS power-law fit directly on ratio data (ratio ~ B · nᵝ)
        rfit = powerlaw_fit(n_v, r_arr)
        zero_note = f'  ({n_zero} zeros)' if n_zero > 0 else ''
        if rfit is not None:
            B   = 10 ** rfit['logA']
            x_f = np.logspace(np.log10(max(n_v.min(), 1.5)), np.log10(n_v.max()), 250)
            ax.plot(x_f, B * x_f ** rfit['alpha'],
                    color=color, lw=2.2, zorder=3, alpha=0.90)
            lbl = f"{label}  ($\\beta = {rfit['alpha']:.2f}$){zero_note}"
        else:
            lbl = f"{label}{zero_note}"

        # Legend proxy
        ax.plot([], [], color=color, lw=2.2, marker=marker, ms=8, label=lbl)

        # Binned medians ± IQR
        df_v2 = df_v.copy()
        df_v2['n_ercs'] = n_v
        bm = binned_medians(df_v2, ratio_col)
        for _, row in bm.iterrows():
            if row['med'] <= 0 or np.isnan(row['med']):
                continue
            y_lo_clamp = max(row['q25'], 1e-6)
            yerr_lo = max(row['med'] - y_lo_clamp, 0)
            yerr_hi = max(row['q75'] - row['med'], 0)
            ax.errorbar(row['x_med'], row['med'],
                        yerr=[[yerr_lo], [yerr_hi]],
                        fmt=marker, color=color, ms=8, lw=1.8,
                        capsize=4, zorder=5, elinewidth=1.5,
                        mfc=color if row['reliable'] else 'white')

    # Reference: ratio = 1 (all pairs have the property)
    ax.axhline(1.0, color='#555555', ls='--', lw=0.9, alpha=0.4, zorder=0)

    ax.set_xscale('log')
    ax.set_yscale('log')
    y_lo = max(min(y_mins) * 0.15, 1e-5) if y_mins else 1e-4
    ax.set_ylim(y_lo, 2.5)
    ax.set_xlabel('Number of ERCs  ($n$)', fontsize=_FLAB)
    ax.set_ylabel(f'Fraction of {denom_label}', fontsize=_FLAB)
    ax.set_title(title, fontsize=_FTIT, fontweight='bold')
    ax.legend(fontsize=_FLEG, loc='upper right', framealpha=0.82,
              borderpad=0.6, labelspacing=0.5)
    ax.tick_params(labelsize=_FTCK)
    ax.grid(True, alpha=0.20, which='both')
    ax.yaxis.set_minor_locator(matplotlib.ticker.LogLocator(
        base=10, subs=np.arange(2, 10) * 0.1, numticks=20))
    ax.grid(True, which='minor', axis='y', alpha=0.08)


# ── Draw panels ───────────────────────────────────────────────────────────────
# Binary synergy counts
draw_abs_panel(ax_syn_bin, sn,
               [('n_basic_pairs',       '#E74C3C', 'Basic / maximal', 'o'),
                ('n_fundamental_pairs', '#27AE60', 'Fundamental',     's')],
               title='Binary synergy: count vs ERCs',
               max_line='n2')

# Complementarity counts
draw_abs_panel(ax_cmp, cn,
               [('n_complementary_pairs',      '#E74C3C', 'All complementary',  'o'),
                ('n_pure_complementary_pairs', '#E67E22', 'Pure complementary', '^'),
                ('n_fundamental_edges',        '#27AE60', 'Fundamental',        's')],
               title='Complementarity: count vs ERCs',
               max_line='n2')

# Ternary synergy counts (separate panel, separate scale)
if _has_ternary and ax_syn_ter is not None:
    draw_abs_panel(ax_syn_ter, sn3,
                   [('n_basic_3syn_triples',       '#9B59B6', 'Ternary basic',       'o'),
                    ('n_fundamental_3syn_triples', '#1ABC9C', 'Ternary fundamental', 's')],
                   title='Ternary synergy: count vs ERCs',
                   max_line='n3')

# Binary synergy ratio (log–log)
draw_ratio_panel(ax_sratio, sn,
                 [('ratio_basic',       '#E74C3C', 'Basic / maximal', 'o'),
                  ('ratio_fundamental', '#27AE60', 'Fundamental',     's')],
                 title='Binary synergy fraction vs ERCs',
                 denom_label='$\\binom{n}{2}$')

# Complementarity ratio (log–log)
draw_ratio_panel(ax_cratio, cn,
                 [('ratio_complementary', '#E74C3C', 'All complementary',  'o'),
                  ('ratio_pure',          '#E67E22', 'Pure complementary', '^'),
                  ('ratio_fundamental',   '#27AE60', 'Fundamental',        's')],
                 title='Complementarity fraction vs ERCs',
                 denom_label='$\\binom{n}{2}$')

# Ternary synergy ratio (log–log)
if _has_ternary and ax_sratio3 is not None:
    draw_ratio_panel(ax_sratio3, sn3,
                     [('ratio_3syn_basic',       '#9B59B6', 'Ternary basic',       'o'),
                      ('ratio_3syn_fundamental', '#1ABC9C', 'Ternary fundamental', 's')],
                     title='Ternary synergy fraction vs ERCs',
                     denom_label='$\\binom{n}{3}$')

# ── Legend panel ──────────────────────────────────────────────────────────────
from matplotlib.lines import Line2D

legend_items = [
    Line2D([0], [0], color='gray', lw=2.2,
           label='Power-law fit  $y = A \\cdot n^{\\alpha}$\n'
                 '(straight line on log–log axes)'),
    Line2D([0], [0], color='#555555', lw=1.0, ls='--', alpha=0.5,
           label='Theoretical max  $\\binom{n}{k}$'),
    Line2D([0], [0], marker='o', color='gray', ms=8, lw=0,
           mfc='gray',  label=f'Bin median ± IQR  ($n \\geq {MIN_N_RELIABLE}$, reliable)'),
    Line2D([0], [0], marker='o', color='gray', ms=8, lw=0,
           mfc='white', label=f'Bin median ± IQR  ($n < {MIN_N_RELIABLE}$, sparse)'),
    Line2D([0], [0], marker='o', color='gray', ms=7,  lw=0, alpha=0.7,
           label='BioModels network'),
    Line2D([0], [0], marker='*', color='gray', ms=10, lw=0, alpha=0.7,
           label='BiGG network'),
    Line2D([0], [0], marker='s', color='gray', ms=7,  lw=0, alpha=0.7,
           label='Other network'),
    Line2D([0], [0], color='none',
           label='Count panels: $\\alpha$ (Poisson GLM, incl. zeros)\n'
                 'Ratio panels:  $\\beta = \\alpha - k$ '
                 '($k=2$ binary, $k=3$ ternary);\n'
                 '$\\beta < 0$ means fraction shrinks with $n$'),
]
ax_leg.legend(handles=legend_items, fontsize=_FLEG, loc='center',
              title='Symbol guide', title_fontsize=_FLEG + 1,
              framealpha=0.92, borderpad=1.0, labelspacing=0.85)

fig.suptitle(
    f'Growth tendency of synergy and complementarity with ERC count\n'
    f'({len(sn)} non-trivial networks;  '
    f'lines are power-law fits $y \\propto n^{{\\alpha}}$;  '
    f'ratio panels use log–log axes)',
    fontsize=_FSUP, fontweight='bold'
)

out_path = os.path.join(OUT_DIR, 'growth_tendency.png')
plt.savefig(out_path, dpi=200, bbox_inches='tight')
print(f"Saved: {os.path.abspath(out_path)}")
plt.show()
