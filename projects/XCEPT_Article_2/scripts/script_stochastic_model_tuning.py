#!/usr/bin/env python3
"""
script_stochastic_model_tuning.py  —  Stochastic APS Model Tuning
==================================================================

PURPOSE
-------
Run the stochastic budget-allocation simulation from the mean initial
conditions of each regime (R1–R6) and inspect trajectories visually.
Use this to tune kinetic parameters (k, K values) and see how different
climate conditions shift the dynamics before running Monte Carlo analysis.

WORKFLOW
--------
1. Load DEFAULT_KINETICS from aps_lib_stochastic (or override via --set).
2. Run N stochastic replicates from mean IC of every selected regime,
   under every selected climate condition.
3. Plot:
   (a) Species time series  — mean ± std band across replicates
   (b) Flux time series     — mean reaction rates per group
   (c) Regime trajectory    — classified regime at each time step (heatmap)
   (d) Budget utilisation   — idle fraction per species over time
4. Print summary table: classified regime at t=5/15/T_max per (regime, climate).
5. Optionally save plots to visualizations/stochastic_tuning/.

CLI EXAMPLES
------------
  # Quick check — all regimes, neutral climate, 10 replicates:
  python script_stochastic_model_tuning.py

  # Three climate conditions, 30 replicates, save figures:
  python script_stochastic_model_tuning.py --climate 0.2 0.5 0.8 --runs 30 --save

  # Only regimes 4-6, dry climate, longer horizon:
  python script_stochastic_model_tuning.py --regime 4 5 6 --climate 0.2 --T 40

  # Override two kinetic parameters:
  python script_stochastic_model_tuning.py --set r7.k=0.3 --set r26.k=5.0

  # Faster time step for finer resolution:
  python script_stochastic_model_tuning.py --dt 0.02 --T 20 --runs 5
"""

import os
import sys
import argparse
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
from typing import Dict, List, Optional

warnings.filterwarnings('ignore', category=RuntimeWarning)

# ── library import ────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)

from aps_lib_stochastic import (
    SPECIES, REACTION_DEFS,
    REGIME_NAMES, REGIME_COLORS, REGIME_IC_RANGES,
    DEFAULT_KINETICS, _CI_THRESHOLDS,
    regime_mean_ic, classify_regime, conflict_index, make_simulator,
)

VIS_DIR = os.path.join(SCRIPT_DIR, 'visualizations', 'stochastic_tuning')
os.makedirs(VIS_DIR, exist_ok=True)

# ─────────────────────────────────────────────────────────────────────────────
# REACTION GROUPS  (for grouped flux plots)
# ─────────────────────────────────────────────────────────────────────────────
REACTION_GROUPS = {
    'Land':       ['r1', 'r2', 'r3', 'r4', 'r5', 'r6'],
    'Economy':    ['r7', 'r8', 'r9', 'r10'],
    'Migration':  ['r12', 'r13', 'r14', 'r15', 'r16'],
    'Resilience': ['r17', 'r18', 'r19', 'r20'],
    'Trust':      ['r21', 'r22', 'r23', 'r24', 'r25'],
    'Violence':   ['r26', 'r27', 'r28', 'r29', 'r30'],
    'ArmedGroup': ['r31', 'r32', 'r33', 'r34'],
}

SPECIES_GROUPS = {
    'Population-RL':  ['SR_RL', 'WR_RL', 'AG_RL'],
    'Population-SL':  ['SR_SL', 'WR_SL', 'AG_SL'],
    'Land':           ['RL', 'SL'],
    'Economy':        ['E_peace', 'E_conflict'],
    'Soc/Gov':        ['T', 'V', 'Gov'],
}

SP_COLORS = {
    'SR_RL': '#27ae60', 'SR_SL': '#82e0aa',
    'WR_RL': '#e67e22', 'WR_SL': '#f0b27a',
    'AG_RL': '#c0392b', 'AG_SL': '#f1948a',
    'RL':    '#2980b9', 'SL':    '#85c1e9',
    'E_peace': '#1abc9c', 'E_conflict': '#e74c3c',
    'T':    '#9b59b6', 'V': '#e74c3c', 'Gov': '#2980b9',
}

# ─────────────────────────────────────────────────────────────────────────────
# SIMULATION HELPERS
# ─────────────────────────────────────────────────────────────────────────────

def run_replicates(regime: int, climate: float, kinetics: dict,
                   T_max: float, dt: float, n_runs: int,
                   base_seed: int = 0):
    """
    Run n_runs stochastic replicates from mean IC of a regime.

    Returns
    -------
    traj_list : list of DataFrames  (one per run, species columns)
    flux_list : list of DataFrames  (one per run, reaction columns)
    """
    ic = regime_mean_ic(regime)
    traj_list, flux_list = [], []
    for i in range(n_runs):
        sim = make_simulator(kinetics=kinetics, dt=dt, climate=climate,
                             rng_seed=base_seed + i)
        traj, flux = sim.simulate(ic, T_max)
        traj_list.append(traj)
        flux_list.append(flux)
    return traj_list, flux_list


def stack_mean_std(df_list: list):
    """
    Given a list of DataFrames with identical index and columns,
    return (mean_df, std_df).
    """
    arr = np.stack([df.values for df in df_list], axis=0)  # (runs, T, cols)
    mean_df = pd.DataFrame(arr.mean(axis=0),
                           index=df_list[0].index,
                           columns=df_list[0].columns)
    std_df  = pd.DataFrame(arr.std(axis=0),
                           index=df_list[0].index,
                           columns=df_list[0].columns)
    return mean_df, std_df


def classify_trajectory(traj_df: pd.DataFrame) -> pd.Series:
    """Return classified regime (1-6) at each time step."""
    return pd.Series(
        [classify_regime(dict(row)) for _, row in traj_df.iterrows()],
        index=traj_df.index, name='regime',
    )


# ─────────────────────────────────────────────────────────────────────────────
# CONFLICT INDICATORS
# ─────────────────────────────────────────────────────────────────────────────
#
# All indicators are ratios in [0, 1].  ε prevents division by zero.
# They are computed from the mean trajectory DataFrame row-by-row.
#
# Population groups used:
#   civilian_pop = SR_RL + SR_SL + WR_RL + WR_SL
#   ag_pop       = AG_RL + AG_SL
#   total_pop    = civilian_pop + ag_pop
#
# Indicator catalogue:
#   ag_power        AG_RL / (AG_RL + SR_RL)
#                   Who controls productive territory. 0 = civilians, 1 = AG.
#
#   poverty         WR_SL / (SR_RL + WR_SL)
#                   Worst-off civilians (poor on bad land) vs best-off
#                   (resilient on good land). Inequality / displacement proxy.
#
#   ag_infiltration AG_SL / (WR_SL + AG_SL)
#                   AG scouts among vulnerable communities — recruitment
#                   pressure & early-warning of escalation.
#
#   militarisation  (AG_RL + AG_SL) / total_pop
#                   Overall weight of the armed group in the population.
#
#   displacement    (SR_SL + WR_SL) / civilian_pop
#                   Fraction of civilians on scarce land — livelihood collapse.
#
#   resilience_ratio (SR_RL + SR_SL) / civilian_pop
#                   Fraction of civilians who can sustain themselves.
#                   Long-run stability predictor.
#
#   war_economy     E_conflict / (E_peace + E_conflict)
#                   Structural lock-in: once war economy dominates, peace
#                   dividends shrink and disarmament becomes harder.
#
#   social_cohesion T / (T + V)
#                   Trust vs violence balance. High = functional social fabric.

EPS = 1e-9

INDICATORS = {
    'ag_power':         ('AG Power\n(AG_RL / AG_RL+SR_RL)',        '#c0392b'),
    'poverty':          ('Poverty\n(WR_SL / SR_RL+WR_SL)',         '#e67e22'),
    'ag_infiltration':  ('AG Infiltration\n(AG_SL / WR_SL+AG_SL)', '#8e44ad'),
    'militarisation':   ('Militarisation\n(AG_pop / total_pop)',    '#2c3e50'),
    'displacement':     ('Displacement\n(civilians on SL)',         '#d35400'),
    'resilience_ratio': ('Resilience Ratio\n(SR_pop / civilian)',   '#27ae60'),
    'war_economy':      ('War Economy Share\n(Ec / Ep+Ec)',         '#e74c3c'),
    'social_cohesion':  ('Social Cohesion\n(T / T+V)',              '#2980b9'),
}


def compute_indicators(traj: pd.DataFrame) -> pd.DataFrame:
    """
    Compute all conflict indicators from a species time-series DataFrame.

    Parameters
    ----------
    traj : DataFrame with species as columns (output of simulate or mean_df)

    Returns
    -------
    DataFrame with indicator names as columns, same index as traj.
    """
    def col(name):
        return traj[name].values if name in traj.columns else np.zeros(len(traj))

    SR_RL = col('SR_RL');  SR_SL = col('SR_SL')
    WR_RL = col('WR_RL');  WR_SL = col('WR_SL')
    AG_RL = col('AG_RL');  AG_SL = col('AG_SL')
    Ep    = col('E_peace');  Ec   = col('E_conflict')
    T     = col('T');        V    = col('V')

    civilian_pop = SR_RL + SR_SL + WR_RL + WR_SL
    ag_pop       = AG_RL + AG_SL
    total_pop    = civilian_pop + ag_pop

    war_economy     = Ec / (Ep + Ec + EPS)
    poverty         = WR_SL / (SR_RL + WR_SL + EPS)
    social_cohesion = T / (T + V + EPS)

    data = {
        'ag_power':         AG_RL           / (AG_RL + SR_RL + EPS),
        'poverty':          poverty,
        'ag_infiltration':  AG_SL           / (WR_SL + AG_SL + EPS),
        'militarisation':   ag_pop          / (total_pop + EPS),
        'displacement':     (SR_SL + WR_SL) / (civilian_pop + EPS),
        'resilience_ratio': (SR_RL + SR_SL) / (civilian_pop + EPS),
        'war_economy':      war_economy,
        'social_cohesion':  social_cohesion,
        # Conflict Index: primary composite classifier (range [0, 3])
        'conflict_index':   war_economy + poverty + (1.0 - social_cohesion),
    }
    return pd.DataFrame(data, index=traj.index)


# ─────────────────────────────────────────────────────────────────────────────
# PLOT E: Conflict indicators
# One panel per indicator; all regimes as coloured lines (mean ± std band).
# ─────────────────────────────────────────────────────────────────────────────

def plot_indicators(results: dict, climate: float, regimes: list,
                    save: bool, suffix: str = ''):
    """
    Compute and plot all conflict indicators plus the Conflict Index (CI).

    results : {regime: (traj_mean, traj_std)}

    Layout: 3 rows × 3 cols.
      Panels 1-8  : individual indicators (ratio [0,1], one colour per regime)
      Panel 9     : Conflict Index (CI ∈ [0,3]), with regime-boundary bands and
                    5 horizontal threshold lines.  This is the primary summary panel.
    """
    ratio_inds = list(INDICATORS.keys())  # 8 ratio indicators
    n_cols = 3
    n_rows = 3  # 3×3 = 9 panels

    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(5.5 * n_cols, 3.8 * n_rows),
                             sharex=True)
    axes = np.array(axes).flatten()
    fig.suptitle(f'Conflict indicators — climate={climate:.2f}', fontsize=13)

    # Pre-compute indicators for each regime
    ind_data = {}
    for regime in regimes:
        mean_df, std_df = results[regime]
        ind_mean  = compute_indicators(mean_df)
        ind_upper = compute_indicators(mean_df + std_df.clip(lower=0))
        ind_lower = compute_indicators((mean_df - std_df).clip(lower=0))
        ind_data[regime] = (ind_mean, ind_lower, ind_upper)

    # ── Panels 0-7: ratio indicators ─────────────────────────────────────────
    for ax_i, ind_name in enumerate(ratio_inds):
        ax = axes[ax_i]
        ind_label, _ = INDICATORS[ind_name]

        for regime in regimes:
            ind_mean, ind_lower, ind_upper = ind_data[regime]
            col  = REGIME_COLORS[regime]
            time = ind_mean.index
            mu   = ind_mean[ind_name]
            lo   = ind_lower[ind_name]
            hi   = ind_upper[ind_name]
            ax.plot(time, mu, color=col, lw=1.5,
                    label=f'R{regime}' if ax_i == 0 else None)
            ax.fill_between(time, lo.clip(0, 1), hi.clip(0, 1),
                            color=col, alpha=0.12)

        ax.axhline(0.5, color='#cccccc', lw=0.8, ls='--')
        ax.set_ylim(0, 1)
        ax.set_title(ind_label, fontsize=8.5)
        ax.set_ylabel('ratio [0–1]', fontsize=7)
        ax.set_xlabel('Time (years)', fontsize=7)
        ax.grid(True, lw=0.35, alpha=0.45)
        ax.tick_params(labelsize=7)

    axes[0].legend(fontsize=7, ncol=2, loc='upper right')

    # ── Panel 8: Conflict Index ───────────────────────────────────────────────
    ax = axes[8]

    # Shade regime bands (background stripes)
    band_edges = [0.0] + list(_CI_THRESHOLDS) + [3.0]
    for r_i in range(6):
        lo_b, hi_b = band_edges[r_i], band_edges[r_i + 1]
        ax.axhspan(lo_b, hi_b,
                   color=REGIME_COLORS[r_i + 1], alpha=0.08)
        ax.text(0.01, (lo_b + hi_b) / 2,
                f'R{r_i + 1}', transform=ax.get_yaxis_transform(),
                fontsize=7, color=REGIME_COLORS[r_i + 1],
                va='center', fontweight='bold')

    # Threshold lines
    for thresh in _CI_THRESHOLDS:
        ax.axhline(thresh, color='#999999', lw=0.8, ls='--')

    for regime in regimes:
        ind_mean, ind_lower, ind_upper = ind_data[regime]
        col  = REGIME_COLORS[regime]
        time = ind_mean.index
        mu   = ind_mean['conflict_index']
        lo   = ind_lower['conflict_index']
        hi   = ind_upper['conflict_index']
        ax.plot(time, mu, color=col, lw=2.0,
                label=f'R{regime}')
        ax.fill_between(time, lo.clip(0, 3), hi.clip(0, 3),
                        color=col, alpha=0.12)

    ax.set_ylim(0, 3)
    ax.set_title('Conflict Index  (war_economy + poverty + 1–cohesion)',
                 fontsize=8.5)
    ax.set_ylabel('CI [0–3]', fontsize=7)
    ax.set_xlabel('Time (years)', fontsize=7)
    ax.legend(fontsize=7, ncol=2, loc='upper left')
    ax.grid(True, lw=0.35, alpha=0.45)
    ax.tick_params(labelsize=7)

    plt.tight_layout()
    if save:
        fname = os.path.join(VIS_DIR,
                             f'indicators_climate{climate:.2f}{suffix}.png')
        fig.savefig(fname, dpi=130, bbox_inches='tight')
        print(f'  Saved: {fname}')
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# PLOT A: Species time series  (one panel per regime)
# ─────────────────────────────────────────────────────────────────────────────

def plot_species_ts(results: dict, climate: float, regimes: list,
                    save: bool, suffix: str = '',
                    plot_species: list = None, tag: str = ''):
    """
    results: {regime: (traj_mean, traj_std)}
    One 6-panel figure, each panel shows key species groups.

    plot_species : list of species names to plot (default: all key species)
    tag          : short label appended to filename to distinguish multiple calls
    """
    if plot_species is None:
        plot_species = ['SR_RL', 'WR_RL', 'SR_SL', 'WR_SL', 'AG_RL', 'AG_SL',
                        'E_peace', 'T', 'V', 'Gov']

    n_reg  = len(regimes)
    n_cols = min(3, n_reg)
    n_rows = int(np.ceil(n_reg / n_cols))

    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(6 * n_cols, 4 * n_rows),
                             sharex=True)
    axes = np.array(axes).flatten()
    title_vars = ', '.join(plot_species[:4]) + ('...' if len(plot_species) > 4 else '')
    fig.suptitle(f'Species dynamics [{title_vars}] — climate={climate:.2f}', fontsize=11)

    for ax_i, regime in enumerate(regimes):
        ax = axes[ax_i]
        mean_df, std_df = results[regime]
        time = mean_df.index

        for sp in plot_species:
            if sp not in mean_df.columns:
                continue
            col = SP_COLORS.get(sp, '#555555')
            mu  = mean_df[sp]
            sd  = std_df[sp]
            ax.plot(time, mu, color=col, lw=1.5, label=sp)
            ax.fill_between(time, (mu - sd).clip(0),
                            mu + sd, color=col, alpha=0.15)

        ax.set_title(REGIME_NAMES[regime], fontsize=9,
                     color=REGIME_COLORS[regime])
        ax.set_ylabel('Amount')
        ax.set_xlabel('Time (years)')
        ax.grid(True, lw=0.4, alpha=0.5)
        if ax_i == 0:
            ax.legend(fontsize=6, ncol=2, loc='upper right')

    for ax_i in range(len(regimes), len(axes)):
        axes[ax_i].set_visible(False)

    plt.tight_layout()
    if save:
        fname = os.path.join(VIS_DIR, f'species_ts{tag}_climate{climate:.2f}{suffix}.png')
        fig.savefig(fname, dpi=120, bbox_inches='tight')
        print(f'  Saved: {fname}')
    return fig

# ─────────────────────────────────────────────────────────────────────────────
# PLOT B: Flux time series  (one panel per regime, grouped reactions)
# ─────────────────────────────────────────────────────────────────────────────

def plot_flux_ts(flux_results: dict, climate: float, regimes: list,
                 save: bool, suffix: str = ''):
    """
    flux_results: {regime: (flux_mean, flux_std)}
    """
    n_reg  = len(regimes)
    n_cols = min(3, n_reg)
    n_rows = int(np.ceil(n_reg / n_cols))

    group_colors = {
        'Land': '#3498db', 'Economy': '#1abc9c', 'Migration': '#e67e22',
        'Resilience': '#27ae60', 'Trust': '#9b59b6',
        'Violence': '#e74c3c', 'ArmedGroup': '#2c3e50',
    }

    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(6 * n_cols, 4 * n_rows),
                             sharex=True)
    axes = np.array(axes).flatten()
    fig.suptitle(f'Total flux per group — climate={climate:.2f}', fontsize=13)

    for ax_i, regime in enumerate(regimes):
        ax = axes[ax_i]
        flux_mean, _ = flux_results[regime]
        time = flux_mean.index

        for grp, rxns in REACTION_GROUPS.items():
            present = [r for r in rxns if r in flux_mean.columns]
            if not present:
                continue
            total = flux_mean[present].sum(axis=1)
            ax.plot(time, total, color=group_colors[grp], lw=1.5, label=grp)

        ax.set_title(REGIME_NAMES[regime], fontsize=9,
                     color=REGIME_COLORS[regime])
        ax.set_ylabel('Events / year')
        ax.set_xlabel('Time (years)')
        ax.grid(True, lw=0.4, alpha=0.5)
        if ax_i == 0:
            ax.legend(fontsize=6, ncol=1, loc='upper right')

    for ax_i in range(len(regimes), len(axes)):
        axes[ax_i].set_visible(False)

    plt.tight_layout()
    if save:
        fname = os.path.join(VIS_DIR, f'flux_ts_climate{climate:.2f}{suffix}.png')
        fig.savefig(fname, dpi=120, bbox_inches='tight')
        print(f'  Saved: {fname}')
    return fig

# ─────────────────────────────────────────────────────────────────────────────
# PLOT C: Regime trajectory heatmap
# ─────────────────────────────────────────────────────────────────────────────

def plot_regime_heatmap(traj_results: dict, climate: float, regimes: list,
                        save: bool, suffix: str = ''):
    """
    Heatmap: rows = starting regime, columns = time, colour = classified regime.
    """
    traj_classes = {}
    for regime in regimes:
        mean_df, _ = traj_results[regime]
        traj_classes[regime] = classify_trajectory(mean_df)

    all_times  = traj_classes[regimes[0]].index
    cmap_vals  = np.array([[traj_classes[r].values] for r in regimes])[:, 0, :]

    fig, ax = plt.subplots(figsize=(10, 2 + 0.6 * len(regimes)))
    im = ax.imshow(cmap_vals, aspect='auto', origin='upper',
                   extent=[all_times[0], all_times[-1], len(regimes) + 0.5, 0.5],
                   vmin=1, vmax=6,
                   cmap='RdYlGn_r')

    ax.set_yticks(range(1, len(regimes) + 1))
    ax.set_yticklabels([REGIME_NAMES[r] for r in regimes], fontsize=8)
    ax.set_xlabel('Time (years)')
    ax.set_title(f'Classified regime over time — climate={climate:.2f}', fontsize=11)

    cbar = fig.colorbar(im, ax=ax, ticks=range(1, 7), fraction=0.03)
    cbar.ax.set_yticklabels([f'R{i}' for i in range(1, 7)], fontsize=7)

    plt.tight_layout()
    if save:
        fname = os.path.join(VIS_DIR, f'regime_traj_climate{climate:.2f}{suffix}.png')
        fig.savefig(fname, dpi=120, bbox_inches='tight')
        print(f'  Saved: {fname}')
    return fig

# ─────────────────────────────────────────────────────────────────────────────
# PLOT D: Budget utilisation (idle fraction per key species)
# ─────────────────────────────────────────────────────────────────────────────

def plot_budget_utilisation(traj_results: dict, flux_results: dict,
                             climate: float, regimes: list,
                             kinetics: dict, dt: float,
                             save: bool, suffix: str = ''):
    """
    For a subset of species, show fraction of budget idle vs. allocated
    across time, for each regime (mean trajectory).
    Key species: SR_RL, WR_RL, WR_SL, E_peace, Gov, V.
    """
    track_species = ['AG_RL', 'E_conflict', 'V', 'T', 'E_peace', 'Gov']

    # Pre-build: for each tracked species, which reactions consume it?
    consuming: Dict[str, List[str]] = {sp: [] for sp in track_species}
    for rxn_id, rxn in REACTION_DEFS.items():
        for sp in rxn.get('consumed', {}):
            if sp in consuming:
                consuming[sp].append(rxn_id)

    n_sp   = len(track_species)
    n_reg  = len(regimes)
    fig, axes = plt.subplots(n_sp, n_reg,
                             figsize=(3.5 * n_reg, 2.5 * n_sp),
                             sharex=True, sharey=True)
    if n_reg == 1:
        axes = axes[:, np.newaxis]
    if n_sp == 1:
        axes = axes[np.newaxis, :]

    fig.suptitle(f'Budget utilisation (idle fraction) — climate={climate:.2f}',
                 fontsize=12)

    for col, regime in enumerate(regimes):
        mean_df, _      = traj_results[regime]
        flux_mean, _    = flux_results[regime]
        time            = mean_df.index[1:]  # flux starts at dt

        for row, sp in enumerate(track_species):
            ax = axes[row, col]
            rxns = consuming[sp]
            if not rxns or sp not in mean_df.columns:
                ax.text(0.5, 0.5, 'n/a', transform=ax.transAxes,
                        ha='center', va='center', fontsize=8)
                continue

            # Compute allocated fraction: Σ_r stoich * flux_r * dt / max(amount, ε)
            amounts = mean_df[sp].values[1:]  # align with flux times
            allocated = np.zeros_like(amounts)
            for rxn_id in rxns:
                stoich = REACTION_DEFS[rxn_id]['consumed'][sp]
                if rxn_id in flux_mean.columns:
                    allocated += stoich * flux_mean[rxn_id].values * dt

            idle_frac = np.clip(1.0 - allocated / np.maximum(amounts, 1e-9), 0, 1)

            ax.fill_between(time, 0, idle_frac, alpha=0.6,
                            color=SP_COLORS.get(sp, '#999999'), label='idle')
            ax.fill_between(time, idle_frac, 1, alpha=0.4,
                            color='#555555', label='allocated')
            ax.set_ylim(0, 1)
            ax.set_xlim(time[0], time[-1])

            if col == 0:
                ax.set_ylabel(sp, fontsize=8)
            if row == 0:
                ax.set_title(f'R{regime}', fontsize=9,
                             color=REGIME_COLORS[regime])
            ax.grid(True, lw=0.3, alpha=0.4)

    # Shared legend
    legend_patches = [
        Patch(color='#27ae60', alpha=0.6, label='idle'),
        Patch(color='#555555', alpha=0.4, label='allocated'),
    ]
    fig.legend(handles=legend_patches, loc='lower right', fontsize=8)

    plt.tight_layout(rect=[0, 0.02, 1, 1])
    if save:
        fname = os.path.join(VIS_DIR,
                             f'budget_util_climate{climate:.2f}{suffix}.png')
        fig.savefig(fname, dpi=120, bbox_inches='tight')
        print(f'  Saved: {fname}')
    return fig

# ─────────────────────────────────────────────────────────────────────────────
# SUMMARY TABLE
# ─────────────────────────────────────────────────────────────────────────────

def print_summary(all_results: dict, T_max: float, regimes: list,
                  climates: list):
    """
    Print: for each (regime, climate), classified regime at t=5, 15, T_max.
    all_results: {climate: {regime: (traj_mean, traj_std)}}
    """
    checkpoints = [cp for cp in [5.0, 15.0, T_max] if cp <= T_max]
    header_cols = ' | '.join([f't={cp:.0f}' for cp in checkpoints])
    print(f'\n{"="*72}')
    print(f'  Summary: classified regime at checkpoints')
    print(f'  {"Regime":<28}  {"Climate":>7}  | {header_cols}')
    print(f'  {"-"*68}')

    for climate in climates:
        for regime in regimes:
            mean_df, _ = all_results[climate][regime]
            parts = []
            for cp in checkpoints:
                idx = mean_df.index.get_indexer([cp], method='nearest')[0]
                state = dict(mean_df.iloc[idx])
                r = classify_regime(state)
                parts.append(f'  R{r} ({REGIME_NAMES[r][:12]:<12})')
            label = REGIME_NAMES[regime][:26]
            print(f'  {label:<28}  {climate:>7.2f}  |' + ' |'.join(parts))
        print(f'  {"-"*68}')
    print(f'{"="*72}\n')

# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description='Stochastic APS model tuning script')
    p.add_argument('--regime', type=int, nargs='+', default=list(range(1, 7)),
                   help='Which regimes to simulate (default: 1-6)')
    p.add_argument('--climate', type=float, nargs='+', default=[0.1,0.4, 0.7],
                   help='Climate values in [0,1] (default: 0.05)')
    p.add_argument('--dt', type=float, nargs='+', default=0.05,
                   help='Time step in years (default: 0.05)')
    p.add_argument('--T', type=float, default=15.0,
                   help='Simulation horizon in years (default: 30)')
    p.add_argument('--runs', type=int, default=20,
                   help='Stochastic replicates per (regime, climate) (default: 10)')
    p.add_argument('--seed', type=int, default=42,
                   help='Base random seed (default: 42)')
    p.add_argument('--set', type=str, nargs='+', default=[],
                   metavar='rxn.param=value',
                   help='Override kinetic params, e.g. --set r7.k=0.3 r26.K.WR_SL=15')
    p.add_argument('--save', action='store_true', default=True,
                   help='Save plots to visualizations/stochastic_tuning/')
    p.add_argument('--no-show', action='store_true',
                   help='Do not call plt.show() (useful in headless envs)')
    return p.parse_args()


def apply_overrides(kinetics: dict, overrides: List[str]) -> dict:
    """
    Parse strings like 'r7.k=0.3' or 'r26.K.WR_SL=15.0' and update kinetics.
    """
    import copy
    kin = copy.deepcopy(kinetics)
    for spec in overrides:
        try:
            key, val = spec.split('=')
            val = float(val)
            parts = key.split('.')
            rxn_id = parts[0]
            if rxn_id not in kin:
                kin[rxn_id] = {}
            if len(parts) == 2:           # e.g. r7.k
                kin[rxn_id][parts[1]] = val
            elif len(parts) == 3:         # e.g. r26.K.WR_SL
                kin[rxn_id].setdefault(parts[1], {})[parts[2]] = val
            else:
                print(f'  [!] Unrecognised override format: {spec}')
        except Exception as e:
            print(f'  [!] Could not parse override "{spec}": {e}')
    return kin

# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main():
    args    = parse_args()
    regimes = args.regime
    climates= args.climate

    # Apply user overrides
    kinetics = apply_overrides(dict(DEFAULT_KINETICS), args.set)
    if args.set:
        print(f'  Applied overrides: {args.set}')

    print(f'\n{"="*60}')
    print(f'  Stochastic APS Tuning')
    print(f'  Regimes : {regimes}')
    print(f'  Climates: {climates}')
    print(f'  Runs    : {args.runs}  |  dt={args.dt}  |  T={args.T}')
    print(f'{"="*60}')

    all_results      = {}   # {climate: {regime: (traj_mean, traj_std)}}
    all_flux_results = {}   # {climate: {regime: (flux_mean, flux_std)}}

    for climate in climates:
        all_results[climate]      = {}
        all_flux_results[climate] = {}
        print(f'\n  Climate = {climate:.2f}')

        for regime in regimes:
            print(f'    R{regime} ({REGIME_NAMES[regime]}) ...', end=' ', flush=True)
            trajs, fluxes = run_replicates(
                regime, climate, kinetics,
                T_max=args.T, dt=args.dt, n_runs=args.runs,
                base_seed=args.seed,
            )
            traj_mean, traj_std = stack_mean_std(trajs)
            flux_mean, flux_std = stack_mean_std(fluxes)
            all_results[climate][regime]      = (traj_mean, traj_std)
            all_flux_results[climate][regime] = (flux_mean, flux_std)

            # Quick status: starting vs ending classified regime
            r_start = classify_regime(dict(traj_mean.iloc[0]))
            r_end   = classify_regime(dict(traj_mean.iloc[-1]))
            print(f'R{r_start} -> R{r_end}')

    # ── Plots ────────────────────────────────────────────────────────────────
    figs = []
    for climate in climates:
        sfx = f'_c{climate:.2f}'
        figs.append(
            plot_species_ts(all_results[climate], climate, regimes,
                            args.save, sfx,
                            plot_species=['SR_RL', 'WR_RL', 'SR_SL', 'WR_SL', 'AG_RL', 'AG_SL'],
                            tag='_population'))
        figs.append(
            plot_species_ts(all_results[climate], climate, regimes,
                            args.save, sfx,
                            plot_species=['E_conflict', 'E_peace', 'T', 'V', 'Gov', 'AG_RL'],
                            tag='_social'))
        # figs.append(
        #     plot_flux_ts(all_flux_results[climate], climate, regimes,
        #                  args.save, sfx))
        figs.append(
            plot_regime_heatmap(all_results[climate], climate, regimes,
                                args.save, sfx))
        figs.append(
            plot_budget_utilisation(
                all_results[climate], all_flux_results[climate],
                climate, regimes, kinetics, args.dt, args.save, sfx))
        figs.append(
            plot_indicators(all_results[climate], climate, regimes,
                            args.save, sfx))

    # ── Summary table ────────────────────────────────────────────────────────
    print_summary(all_results, args.T, regimes, climates)

    if not args.no_show:
        plt.show()

    print('Done.')


if __name__ == '__main__':
    main()
