#!/usr/bin/env python3
"""
script_strategy_game.py  —  APS Strategic Interaction Model (Phase 1)
======================================================================

Models a two-player game between the Government (Gov) and Armed Group (AG).
Each player can boost the firing rate of their strategic reactions by a
factor `boost` above the natural kinetic baseline.

PLAYERS AND STRATEGIES
-----------------------
Government  : {None, Development, Defense, Balance}
  None        — no strategic boost; pure baseline kinetics (aps_lib defaults)
  Development — boosts r6, r10, r17, r18, r21, r22  (economy, trust, resilience)
  Defense     — boosts r20, r28  (DDR & counter-insurgency)
  Balance     — boosts all Gov strategic reactions equally

Armed Group : {None, Coerce, Extract, Balance}
  None    — no strategic boost; pure baseline kinetics
  Coerce  — boosts r16, r27, r29, r31  (weaken opponent capacity)
  Extract — boosts r14, r15, r32, r33, r34  (capture resources/territory)
  Balance — boosts all AG strategic reactions equally

IMPLEMENTATION
--------------
Strategy boost = multiply k_r by `boost` for each reaction in the chosen
strategy.  Natural reactions (r1–r5, r7–r9, r12–r13, r19, r23–r26, r30)
are NEVER boosted regardless of strategy.

PHASE 1 (this script)
  Both players choose one fixed strategy and hold it for all N periods.
  Every combination (4×4 = 16) is run and compared.

OUTPUTS (per regime × climate)
  1. Trajectory plots  — chosen metric over time, 4×4 panel grid
  2. Payoff matrices   — 4×4 heatmaps (one per outcome metric, Δ vs None×None)
  3. Reaction activity — stacked bar: firing events/yr per reaction group
  4. Regime evolution  — heatmap: regime at end of each period, all 16 combos
  5. Console tables    — 4×4 payoff table per metric

CLI EXAMPLES
------------
  # All 16 combos, default regime R3, neutral climate, 10 periods of 1 year
  python script_strategy_game.py

  # Regimes R3 and R5, two climates, stronger boost
  python script_strategy_game.py --regime 3 5 --climate 0.2 0.8 --boost 2.0

  # Only Gov=Development vs AG=Extract, longer game
  python script_strategy_game.py --gov None Development --ag None Extract --periods 20

  # Show V (violence) trajectories instead of E_peace
  python script_strategy_game.py --metric V --regime 4
"""

import os
import sys
import copy
import argparse
import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from typing import Dict, List, Tuple

warnings.filterwarnings('ignore', category=RuntimeWarning)

# ── library import ────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)

from aps_lib_stochastic import (
    REACTION_DEFS, DEFAULT_KINETICS,
    REGIME_NAMES, REGIME_COLORS,
    regime_mean_ic, classify_regime, make_simulator,
    _CI_THRESHOLDS, conflict_index,
)

VIS_DIR    = os.path.join(SCRIPT_DIR, 'visualizations', 'strategy_game')
REPORT_DIR = os.path.join(SCRIPT_DIR, 'reports', 'adaptive')
os.makedirs(VIS_DIR,    exist_ok=True)
os.makedirs(REPORT_DIR, exist_ok=True)

# Short labels for rule names used in matrix plot axes
RULE_SHORT = {
    'reactive':          'Rct',
    'development_first': 'Dev1',
    'defense_first':     'Def1',
    'adaptive_balance':  'ABal',
    'opportunistic':     'Opp',
    'territorial':       'Ter',
    'economic':          'Eco',
    'pressure':          'Prs',
}

# Ordered lists (same order as dicts, guaranteed Python 3.7+)
GOV_RULES_LIST = list({'reactive': None, 'development_first': None,
                        'defense_first': None, 'adaptive_balance': None})
AG_RULES_LIST  = list({'opportunistic': None, 'territorial': None,
                        'economic': None, 'pressure': None})

# ─────────────────────────────────────────────────────────────────────────────
# STRATEGY DEFINITIONS
# ─────────────────────────────────────────────────────────────────────────────
# Natural reactions (never boosted): r1-r5, r7-r9, r12-r13, r19, r23-r26, r30
STRATEGY_REACTIONS: Dict[str, Dict[str, List[str]]] = {
    'Gov': {
        'None':        [],
        'Development': ['r6', 'r10', 'r17', 'r18', 'r21', 'r22'],
        'Defense':     ['r20', 'r28'],
        'Balance':     ['r6', 'r10', 'r17', 'r18', 'r21', 'r22', 'r20', 'r28'],
    },
    'AG': {
        'None':    [],
        'Coerce':  ['r16', 'r27', 'r29', 'r31'],
        'Extract': ['r14', 'r15', 'r32', 'r33', 'r34'],
        'Balance': ['r14', 'r15', 'r16', 'r27', 'r29', 'r31', 'r32', 'r33', 'r34'],
    },
}

GOV_STRATEGIES = ['None', 'Development', 'Defense', 'Balance']
AG_STRATEGIES  = ['None', 'Coerce', 'Extract', 'Balance']

REACTION_GROUPS: Dict[str, List[str]] = {
    'Land':       ['r1', 'r2', 'r3', 'r4', 'r5', 'r6'],
    'Economy':    ['r7', 'r8', 'r9', 'r10'],
    'Migration':  ['r12', 'r13', 'r14', 'r15', 'r16'],
    'Resilience': ['r17', 'r18', 'r19', 'r20'],
    'Trust':      ['r21', 'r22', 'r23', 'r24', 'r25'],
    'Violence':   ['r26', 'r27', 'r28', 'r29', 'r30'],
    'ArmedGroup': ['r31', 'r32', 'r33', 'r34'],
}

GROUP_COLORS = {
    'Land': '#3498db', 'Economy': '#1abc9c', 'Migration': '#e67e22',
    'Resilience': '#27ae60', 'Trust': '#9b59b6',
    'Violence': '#e74c3c', 'ArmedGroup': '#2c3e50',
}

# Outcome metrics: name -> (display label, sign)
# sign = +1 means higher is better for Gov; -1 means lower is better for Gov
METRICS: Dict[str, Tuple[str, int]] = {
    'T':             ('Trust',          +1),
    'V':             ('Violence',       -1),
    'Gov':           ('Governance',     +1),
    'E_peace':       ('Peace Econ.',    +1),
    'E_conflict':    ('Conflict Econ.', -1),
    'conflict_index':('Conflict Index', -1),
}

# Short strategy labels for dense displays
SHORT = {'None': 'Nil', 'Development': 'Dev', 'Defense': 'Def',
         'Balance': 'Bal', 'Coerce': 'Crc', 'Extract': 'Ext'}

# Strategy colours for adaptive timeline visualisation
STRATEGY_COLORS = {
    'Gov': {
        'None':        '#aaaaaa',
        'Development': '#27ae60',
        'Defense':     '#2980b9',
        'Balance':     '#8e44ad',
    },
    'AG': {
        'None':    '#aaaaaa',
        'Coerce':  '#e74c3c',
        'Extract': '#e67e22',
        'Balance': '#c0392b',
    },
}

# ─────────────────────────────────────────────────────────────────────────────
# INDICATOR FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────

INDICATOR_NAMES = [
    'war_economy', 'poverty', 'social_cohesion', 'ag_power',
    'militarisation', 'ag_infiltration', 'displacement', 'resilience_ratio',
]

INDICATOR_LABELS = {
    'war_economy':     'War Economy',
    'poverty':         'Poverty',
    'social_cohesion': 'Social Cohesion',
    'ag_power':        'AG Power',
    'militarisation':  'Militarisation',
    'ag_infiltration': 'AG Infiltration',
    'displacement':    'Displacement',
    'resilience_ratio':'Resilience Ratio',
}


def indicators_from_state(state: dict) -> dict:
    """Compute all 8 conflict indicators + CI from a state dict."""
    eps = 1e-9
    SR_RL = state.get('SR_RL', 0.0);  WR_RL = state.get('WR_RL', 0.0)
    SR_SL = state.get('SR_SL', 0.0);  WR_SL = state.get('WR_SL', 0.0)
    AG_RL = state.get('AG_RL', 0.0);  AG_SL = state.get('AG_SL', 0.0)
    Ep    = state.get('E_peace', 0.0); Ec    = state.get('E_conflict', 0.0)
    T_    = state.get('T', 0.0);       V_    = state.get('V', 0.0)
    Gv    = state.get('Gov', 0.0)
    total_pop = SR_RL + WR_RL + SR_SL + WR_SL + AG_RL + AG_SL + eps
    total_civ = SR_RL + WR_RL + SR_SL + WR_SL + eps
    return {
        'war_economy':      Ec  / (Ep + Ec + eps),
        'poverty':          WR_SL / (SR_RL + WR_SL + eps),
        'social_cohesion':  T_  / (T_ + V_ + eps),
        'ag_power':         (AG_RL + AG_SL) / total_pop,
        'militarisation':   V_  / (T_ + V_ + Gv + eps),
        'ag_infiltration':  AG_RL / (AG_RL + AG_SL + eps),
        'displacement':     (SR_SL + WR_SL) / total_civ,
        'resilience_ratio': SR_RL / (SR_RL + WR_RL + eps),
        'conflict_index':   conflict_index(state),
    }


def _ci_from_traj(traj: pd.DataFrame) -> np.ndarray:
    """Vectorised CI computation over a trajectory DataFrame."""
    eps = 1e-9
    def col(name):
        return traj[name].values if name in traj.columns else np.zeros(len(traj))
    Ep = col('E_peace'); Ec = col('E_conflict')
    SR = col('SR_RL');   WS = col('WR_SL')
    T_ = col('T');       V_ = col('V')
    war_economy     = Ec / (Ep + Ec + eps)
    poverty         = WS / (SR + WS  + eps)
    social_cohesion = T_ / (T_ + V_  + eps)
    return war_economy + poverty + (1.0 - social_cohesion)


def _ind_from_traj(traj: pd.DataFrame, key: str) -> np.ndarray:
    """Compute a single indicator as an array over a trajectory DataFrame."""
    eps = 1e-9
    def col(name):
        return traj[name].values if name in traj.columns else np.zeros(len(traj))
    SR_RL = col('SR_RL'); WR_RL = col('WR_RL')
    SR_SL = col('SR_SL'); WR_SL = col('WR_SL')
    AG_RL = col('AG_RL'); AG_SL = col('AG_SL')
    Ep = col('E_peace'); Ec = col('E_conflict')
    T_ = col('T'); V_ = col('V'); Gv = col('Gov')
    total_pop = SR_RL + WR_RL + SR_SL + WR_SL + AG_RL + AG_SL + eps
    total_civ = SR_RL + WR_RL + SR_SL + WR_SL + eps
    if key == 'war_economy':      return Ec  / (Ep + Ec + eps)
    if key == 'poverty':          return WR_SL / (SR_RL + WR_SL + eps)
    if key == 'social_cohesion':  return T_ / (T_ + V_ + eps)
    if key == 'ag_power':         return (AG_RL + AG_SL) / total_pop
    if key == 'militarisation':   return V_ / (T_ + V_ + Gv + eps)
    if key == 'ag_infiltration':  return AG_RL / (AG_RL + AG_SL + eps)
    if key == 'displacement':     return (SR_SL + WR_SL) / total_civ
    if key == 'resilience_ratio': return SR_RL / (SR_RL + WR_RL + eps)
    return np.zeros(len(traj))


# ─────────────────────────────────────────────────────────────────────────────
# ADAPTIVE STRATEGY RULES
# ─────────────────────────────────────────────────────────────────────────────
# Each rule function: (indicators: dict, current_strat: str) -> next_strat: str
# Thresholds refer to _CI_THRESHOLDS = (0.45, 0.94, 1.39, 1.84, 2.31)
#   index: [0]=R1/R2 boundary, [1]=R2/R3, [2]=R3/R4, [3]=R4/R5, [4]=R5/R6

def _gov_reactive(ind, cur):
    """Switch based on current CI level: low->Develop, mid->Balance, high->Defend."""
    ci = ind['conflict_index']
    if ci > _CI_THRESHOLDS[3]: return 'Defense'       # R5-R6 zone
    if ci > _CI_THRESHOLDS[1]: return 'Balance'       # R3-R4 zone
    return 'Development'                               # R1-R2 zone

def _gov_development_first(ind, cur):
    """Always invest in development; pivot to defense only in extreme conflict."""
    return 'Defense' if ind['conflict_index'] > _CI_THRESHOLDS[4] else 'Development'

def _gov_defense_first(ind, cur):
    """Maintain security focus; relax only in near-peace conditions."""
    return 'Development' if ind['conflict_index'] < _CI_THRESHOLDS[0] else 'Defense'

def _gov_adaptive_balance(ind, cur):
    """Context-aware: rebuild trust->Dev; AG threat->Def; otherwise Balance."""
    if ind['social_cohesion'] < 0.35: return 'Development'
    if ind['ag_power']        > 0.25: return 'Defense'
    return 'Balance'

GOV_RULES = {
    'reactive':          _gov_reactive,
    'development_first': _gov_development_first,
    'defense_first':     _gov_defense_first,
    'adaptive_balance':  _gov_adaptive_balance,
}

def _ag_opportunistic(ind, cur):
    """Exploit poverty for recruitment; shift to extraction when war economy is high."""
    if ind['poverty']     > 0.55: return 'Coerce'
    if ind['war_economy'] > 0.45: return 'Extract'
    return 'None'

def _ag_territorial(ind, cur):
    """Prioritise territorial control (Coerce) regardless of conditions."""
    return 'Coerce'

def _ag_economic(ind, cur):
    """Prioritise resource extraction regardless of conditions."""
    return 'Extract'

def _ag_pressure(ind, cur):
    """Coerce in low-conflict (soft targets); extract in high-conflict (war economy)."""
    return 'Coerce' if ind['conflict_index'] < _CI_THRESHOLDS[2] else 'Extract'

AG_RULES = {
    'opportunistic': _ag_opportunistic,
    'territorial':   _ag_territorial,
    'economic':      _ag_economic,
    'pressure':      _ag_pressure,
}


def capacity_check(state: dict, player: str, threshold: float, rng) -> bool:
    """
    Stochastic implementation gate: p = min(1, level / threshold).

    Gov  -> level = state['Gov']
    AG   -> level = state['AG_RL'] + state['AG_SL']

    threshold = 0  (default): always True — full capacity, deterministic.
    threshold > 0: probabilistic — low Gov/AG stock -> low chance to act.
    """
    if threshold <= 0:
        return True
    level = (state.get('Gov', 0.0) if player == 'Gov'
             else state.get('AG_RL', 0.0) + state.get('AG_SL', 0.0))
    return bool(rng.random() < min(1.0, level / threshold))


def apply_rule(player: str, indicators: dict,
               current_strat: str, rule_name: str) -> str:
    """Look up and apply a named rule; hold current strategy if rule is unknown."""
    rules = GOV_RULES if player == 'Gov' else AG_RULES
    fn = rules.get(rule_name)
    return fn(indicators, current_strat) if fn is not None else current_strat


# ─────────────────────────────────────────────────────────────────────────────
# CORE FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────

def apply_strategy_boost(kinetics: dict,
                          gov_strat: str, ag_strat: str,
                          boost: float) -> dict:
    """
    Return a deep copy of kinetics with k values multiplied by `boost`
    for every reaction belonging to the chosen strategies.
    Natural reactions (empty strategy lists) are always left unchanged.
    """
    kin = copy.deepcopy(kinetics)
    boosted = set(STRATEGY_REACTIONS['Gov'][gov_strat] +
                  STRATEGY_REACTIONS['AG'][ag_strat])
    for rxn_id in boosted:
        if rxn_id in kin:
            kin[rxn_id]['k'] = kin[rxn_id]['k'] * boost
    return kin


def compute_metrics(state: dict) -> dict:
    """Extract scalar outcome metrics and classified regime from a state dict."""
    return {
        'T':             state.get('T', 0.0),
        'V':             state.get('V', 0.0),
        'Gov':           state.get('Gov', 0.0),
        'E_peace':       state.get('E_peace', 0.0),
        'E_conflict':    state.get('E_conflict', 0.0),
        'conflict_index': conflict_index(state),
        'regime':        classify_regime(state),
    }


def run_period(ic: dict, climate: float, kinetics: dict,
               gov_strat: str, ag_strat: str, boost: float,
               T_period: float, dt: float, n_runs: int,
               base_seed: int) -> Tuple[dict, pd.DataFrame, pd.DataFrame]:
    """
    Run one game period with strategy boost applied.

    Returns
    -------
    mean_end_state : dict  (species -> mean final amount)
    traj_mean      : DataFrame  (time series, species columns)
    flux_mean      : DataFrame  (time series, reaction columns)
    """
    boosted_kin = apply_strategy_boost(kinetics, gov_strat, ag_strat, boost)
    trajs, fluxes = [], []
    for i in range(n_runs):
        sim = make_simulator(kinetics=boosted_kin, dt=dt,
                             climate=climate, rng_seed=base_seed + i)
        traj, flux = sim.simulate(ic, T_period)
        trajs.append(traj)
        fluxes.append(flux)

    arr_t = np.stack([df.values for df in trajs], axis=0)
    arr_f = np.stack([df.values for df in fluxes], axis=0)
    traj_mean = pd.DataFrame(arr_t.mean(0), index=trajs[0].index,
                             columns=trajs[0].columns)
    flux_mean = pd.DataFrame(arr_f.mean(0), index=fluxes[0].index,
                             columns=fluxes[0].columns)
    return dict(traj_mean.iloc[-1]), traj_mean, flux_mean


def run_game(start_regime: int, climate: float, kinetics: dict,
             gov_strat: str, ag_strat: str, boost: float,
             n_periods: int, T_period: float,
             dt: float, n_runs: int, base_seed: int) -> List[dict]:
    """
    Chain n_periods game periods.  End state of period p becomes IC for p+1.

    Returns
    -------
    period_log : list of dicts, one per period:
        { period, gov_strat, ag_strat, t_offset, T_period, metrics, traj, flux }
    """
    ic = regime_mean_ic(start_regime)
    log = []
    for p in range(n_periods):
        seed_p = base_seed + p * n_runs
        end_state, traj, flux = run_period(
            ic, climate, kinetics, gov_strat, ag_strat, boost,
            T_period, dt, n_runs, seed_p)
        log.append({
            'period':    p + 1,
            'gov_strat': gov_strat,
            'ag_strat':  ag_strat,
            't_offset':  p * T_period,
            'T_period':  T_period,
            'metrics':   compute_metrics(end_state),
            'traj':      traj,
            'flux':      flux,
        })
        ic = end_state
    return log


def run_game_adaptive(
        start_regime: int, climate: float, kinetics: dict,
        gov_init: str, ag_init: str, boost: float,
        n_periods: int, T_period: float, dt: float, n_runs: int, base_seed: int,
        gov_rule: str, ag_rule: str,
        gov_capacity_threshold: float, ag_capacity_threshold: float,
        gov_period_mult: int, ag_period_mult: int,
) -> List[dict]:
    """
    Chain n_periods game periods with adaptive, rule-based strategy updates.

    Strategy update logic (applied at end of period p):
      1. Apply named rule to current indicators -> intended next strategy.
      2. Capacity check: Bernoulli gate with p = min(1, level/threshold).
         Gov  level = state['Gov'];  AG level = AG_RL + AG_SL.
         threshold=0 -> always succeeds (deterministic / full capacity).
      3. If check fails -> effective strategy for next period is 'None'.

    Period multipliers control re-evaluation frequency:
      gov_period_mult=2 -> Gov re-evaluates every 2 periods.
      ag_period_mult=1  -> AG re-evaluates every period (default).

    Returns
    -------
    log : list of dicts, one per period.
          Same schema as run_game plus 'indicators' (end-of-period indicator dict).
    """
    rng = np.random.default_rng(base_seed)
    ic      = regime_mean_ic(start_regime)
    gov_eff = gov_init
    ag_eff  = ag_init
    log     = []

    for p in range(n_periods):
        seed_p = base_seed + p * n_runs
        end_state, traj, flux = run_period(
            ic, climate, kinetics,
            gov_eff, ag_eff, boost,
            T_period, dt, n_runs, seed_p)
        ind = indicators_from_state(end_state)
        log.append({
            'period':     p + 1,
            'gov_strat':  gov_eff,
            'ag_strat':   ag_eff,
            't_offset':   p * T_period,
            'T_period':   T_period,
            'metrics':    compute_metrics(end_state),
            'indicators': ind,
            'traj':       traj,
            'flux':       flux,
        })
        ic = end_state

        # Gov re-evaluation
        if (p + 1) % gov_period_mult == 0:
            nxt     = apply_rule('Gov', ind, gov_eff, gov_rule)
            gov_eff = nxt if capacity_check(end_state, 'Gov',
                                            gov_capacity_threshold, rng) else 'None'

        # AG re-evaluation
        if (p + 1) % ag_period_mult == 0:
            nxt    = apply_rule('AG', ind, ag_eff, ag_rule)
            ag_eff = nxt if capacity_check(end_state, 'AG',
                                           ag_capacity_threshold, rng) else 'None'

    return log


# ─────────────────────────────────────────────────────────────────────────────
# PLOT 1: Trajectory grid
# 4×4 panels (Gov strategy × AG strategy); one metric per figure.
# Each panel: full chained time series with period boundary markers.
# Background tinted by final classified regime.
# ─────────────────────────────────────────────────────────────────────────────

def plot_strategy_trajectories(all_logs: dict, metric: str,
                                climate: float, regime: int,
                                gov_strats: list, ag_strats: list,
                                T_period: float, save: bool):
    """4×4 panel grid showing one outcome metric over the full game horizon."""
    n_gov = len(gov_strats)
    n_ag  = len(ag_strats)
    label, sign = METRICS[metric]

    fig, axes = plt.subplots(n_gov, n_ag,
                             figsize=(3.8 * n_ag, 2.8 * n_gov),
                             sharex=True, sharey=True)
    axes = np.array(axes).reshape(n_gov, n_ag)
    fig.suptitle(
        f'{label} — game trajectory  |  R{regime}  climate={climate:.2f}',
        fontsize=11)

    # Collect global y-range for consistent scale
    all_vals = []
    for gs in gov_strats:
        for ag in ag_strats:
            for entry in all_logs.get((gs, ag), []):
                if metric == 'conflict_index':
                    all_vals.extend(_ci_from_traj(entry['traj']).tolist())
                elif metric in entry['traj'].columns:
                    all_vals.extend(entry['traj'][metric].values.tolist())
    if metric == 'conflict_index':
        y_max = 3.05
    else:
        y_max = max(all_vals) * 1.05 if all_vals else 1.0

    # Pre-compute CI zone label midpoints for annotation
    _ci_bounds = [0.0] + list(_CI_THRESHOLDS) + [3.0]
    _ci_zone_mids = [(_ci_bounds[i] + _ci_bounds[i+1]) / 2 for i in range(6)]

    for r, gs in enumerate(gov_strats):
        for c, ag in enumerate(ag_strats):
            ax = axes[r, c]
            log = all_logs.get((gs, ag), [])

            # Stitch periods into one continuous time series
            t_all, v_all = [], []
            for entry in log:
                t_local = entry['traj'].index.values + entry['t_offset']
                if metric == 'conflict_index':
                    v_local = _ci_from_traj(entry['traj'])
                else:
                    v_local = (entry['traj'][metric].values
                               if metric in entry['traj'].columns
                               else np.zeros(len(t_local)))
                t_all.extend(t_local.tolist())
                v_all.extend(v_local.tolist())

            if t_all:
                ax.plot(t_all, v_all, color='#2c3e50', lw=1.2)

            # CI threshold lines + zone labels
            if metric == 'conflict_index':
                for thr in _CI_THRESHOLDS:
                    ax.axhline(thr, color='#aaaaaa', lw=0.5, ls=':', zorder=1)
                if c == n_ag - 1:   # rightmost column only
                    for zi, zmid in enumerate(_ci_zone_mids):
                        ax.annotate(f'R{zi+1}',
                                    xy=(1.02, zmid),
                                    xycoords=('axes fraction', 'data'),
                                    fontsize=5.5, va='center', ha='left',
                                    color='#777777')

            # Period boundary lines
            for p in range(1, len(log)):
                ax.axvline(p * T_period, color='#bbbbbb', lw=0.7, ls='--')

            # Background tint = final regime colour
            if log:
                r_end = log[-1]['metrics']['regime']
                ax.set_facecolor(REGIME_COLORS[r_end] + '25')

            ax.set_ylim(0, y_max)
            ax.tick_params(labelsize=6)
            ax.grid(True, lw=0.3, alpha=0.35)

            # Row / column labels
            if r == 0:
                ax.set_title(f'AG: {SHORT[ag]}', fontsize=8, pad=3)
            if c == 0:
                ax.set_ylabel(f'{SHORT[gs]}\n{label}', fontsize=7)
            if r == n_gov - 1:
                ax.set_xlabel('years', fontsize=6)

    plt.tight_layout()
    if save:
        fname = os.path.join(VIS_DIR,
                             f'traj_{metric}_R{regime}_c{climate:.2f}.png')
        fig.savefig(fname, dpi=130, bbox_inches='tight')
        print(f'  Saved: {fname}')
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# PLOT 2: Payoff matrices
# One heatmap per outcome metric; cell value = final amount Δ vs None×None.
# ─────────────────────────────────────────────────────────────────────────────

def plot_payoff_matrix(all_logs: dict, climate: float, regime: int,
                        gov_strats: list, ag_strats: list, save: bool):
    """
    5-panel figure: one 4×4 payoff heatmap per outcome metric.
    Colour encodes deviation from the baseline (None×None); absolute final
    value is annotated in each cell.
    """
    n_metrics = len(METRICS)
    ncols = 3
    nrows = (n_metrics + ncols - 1) // ncols
    fig, axes_grid = plt.subplots(nrows, ncols,
                                  figsize=(4.5 * ncols, 4.8 * nrows))
    axes_flat = np.array(axes_grid).flatten()
    fig.suptitle(
        f'Payoff matrix — end of game  |  R{regime}  climate={climate:.2f}',
        fontsize=12)

    baseline_log = all_logs.get(('None', 'None'), [])
    baseline = baseline_log[-1]['metrics'] if baseline_log else {}

    for ax_i, (metric, (label, sign)) in enumerate(METRICS.items()):
        ax = axes_flat[ax_i]
        n_gov = len(gov_strats)
        n_ag  = len(ag_strats)
        abs_data   = np.zeros((n_gov, n_ag))
        delta_data = np.zeros((n_gov, n_ag))
        ref = baseline.get(metric, 0.0)

        for r, gs in enumerate(gov_strats):
            for c, ag in enumerate(ag_strats):
                log = all_logs.get((gs, ag), [])
                val = log[-1]['metrics'].get(metric, 0.0) if log else 0.0
                abs_data[r, c]   = val
                delta_data[r, c] = (val - ref) * sign  # positive = good for Gov

        vmax = max(abs(delta_data).max(), 1e-9)
        # For CI: invert colormap (lower CI = green = better for Gov)
        cmap = 'RdYlGn'
        im = ax.imshow(delta_data, cmap=cmap,
                       vmin=-vmax, vmax=vmax, aspect='auto')

        ax.set_xticks(range(n_ag))
        ax.set_xticklabels([SHORT[ag] for ag in ag_strats],
                           rotation=30, ha='right', fontsize=8)
        ax.set_yticks(range(n_gov))
        ax.set_yticklabels([SHORT[gs] for gs in gov_strats], fontsize=8)
        ax.set_xlabel('AG strategy', fontsize=8)
        if ax_i % ncols == 0:
            ax.set_ylabel('Gov strategy', fontsize=8)
        ax.set_title(label, fontsize=9)

        # Annotate with absolute final value
        for r in range(n_gov):
            for c in range(n_ag):
                ax.text(c, r, f'{abs_data[r, c]:.2f}',
                        ha='center', va='center', fontsize=7.5,
                        color='black',
                        fontweight='bold' if abs(delta_data[r, c]) == vmax else 'normal')

        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04,
                     label='delta vs None×None\n(green = better for Gov)')

    # Hide spare axes
    for ax_i in range(n_metrics, nrows * ncols):
        axes_flat[ax_i].set_visible(False)

    plt.tight_layout()
    if save:
        fname = os.path.join(VIS_DIR,
                             f'payoff_matrix_R{regime}_c{climate:.2f}.png')
        fig.savefig(fname, dpi=130, bbox_inches='tight')
        print(f'  Saved: {fname}')
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# PLOT 2b: Combined CI payoff matrix — all regimes × climates, shared color scale
# ─────────────────────────────────────────────────────────────────────────────

def plot_combined_ci_matrix(
        results_by_rc: dict,   # {(regime, climate): all_logs}
        regimes: list,
        climates: list,
        gov_strats: list,
        ag_strats: list,
        save: bool):
    """
    Single figure with one 4×4 CI payoff heatmap per (regime, climate) combination.

    Layout: rows = regimes, columns = climates.
    Color encodes the ABSOLUTE CI value on a shared scale across all panels
    (RdYlGn_r: red = high CI = bad, green = low CI = good for peace).
    Cell annotations show the same absolute CI value, so numbers and colors
    are directly consistent.  A single colorbar is placed to the right,
    in its own axes, so it never overlaps a matrix panel.
    """
    n_reg = len(regimes)
    n_cli = len(climates)
    n_gov = len(gov_strats)
    n_ag  = len(ag_strats)

    # ── Collect absolute CI matrices and find global min/max ──────────────────
    abs_grids = {}
    for regime in regimes:
        for climate in climates:
            all_logs = results_by_rc.get((regime, climate), {})
            abs_m = np.zeros((n_gov, n_ag))
            for r, gs in enumerate(gov_strats):
                for c, ag in enumerate(ag_strats):
                    log = all_logs.get((gs, ag), [])
                    abs_m[r, c] = (log[-1]['metrics'].get('conflict_index', 0.0)
                                   if log else 0.0)
            abs_grids[(regime, climate)] = abs_m

    all_vals    = np.concatenate([m.ravel() for m in abs_grids.values()])
    global_vmin = float(all_vals.min())
    global_vmax = float(all_vals.max())
    # Round to one decimal for a clean colorbar
    global_vmin = max(0.0, global_vmin - 0.05)
    global_vmax = min(3.0, global_vmax + 0.05)

    # ── Build figure with manual axes to avoid colorbar overlap ───────────────
    cell_w, cell_h = 3.2, 3.2
    cbar_w  = 0.5          # inches reserved for colorbar strip
    left_m  = 0.7          # left margin (ylabel space)
    right_m = cbar_w + 0.4 # right margin (colorbar + gap)
    top_m   = 0.55          # suptitle
    bot_m   = 0.45          # xlabel

    fig_w = left_m + cell_w * n_cli + right_m
    fig_h = top_m  + cell_h * n_reg + bot_m

    fig = plt.figure(figsize=(fig_w, fig_h))
    fig.suptitle(
        'Conflict Index — payoff matrices (shared color scale)\n'
        'Cell value = final CI.   Green = low CI (better).   Red = high CI (worse).',
        fontsize=10, y=1.0 - 0.01 / fig_h,
    )

    # Compute normalized positions for matrix axes
    ax_left   = left_m  / fig_w
    ax_bottom = bot_m   / fig_h
    ax_width  = cell_w  / fig_w
    ax_height = cell_h  / fig_h
    h_gap     = 0.06 / fig_w   # horizontal gap between panels
    v_gap     = 0.45 / fig_h   # vertical gap between panels

    axes = np.empty((n_reg, n_cli), dtype=object)
    im_ref = None

    for ri, regime in enumerate(regimes):
        # rows drawn top-to-bottom: ri=0 is the topmost row
        y0 = ax_bottom + (n_reg - 1 - ri) * (ax_height + v_gap)
        for ci_col, climate in enumerate(climates):
            x0  = ax_left + ci_col * (ax_width + h_gap)
            ax  = fig.add_axes([x0, y0, ax_width, ax_height])
            axes[ri, ci_col] = ax

            abs_m = abs_grids[(regime, climate)]
            im = ax.imshow(abs_m, cmap='RdYlGn_r',
                           vmin=global_vmin, vmax=global_vmax,
                           aspect='auto')
            im_ref = im

            # Tick labels
            ax.set_xticks(range(n_ag))
            ax.set_yticks(range(n_gov))
            ax.set_xticklabels([SHORT[ag] for ag in ag_strats],
                               rotation=30, ha='right', fontsize=7)
            ax.set_yticklabels([SHORT[gs] for gs in gov_strats], fontsize=7)

            # Axis labels only on edges
            if ci_col == 0:
                ax.set_ylabel(f'R{regime}\nGov strat.', fontsize=8)
            else:
                ax.set_yticklabels([])
            if ri == n_reg - 1:
                ax.set_xlabel('AG strategy', fontsize=8)
            else:
                ax.set_xticklabels([])

            ax.set_title(f'R{regime}  climate={climate:.2f}', fontsize=8.5, pad=3)

            # Annotate each cell with absolute CI value
            for r in range(n_gov):
                for c in range(n_ag):
                    val = abs_m[r, c]
                    # Choose white text when background is very dark
                    norm_val = (val - global_vmin) / max(global_vmax - global_vmin, 1e-9)
                    txt_color = 'white' if norm_val > 0.85 else 'black'
                    ax.text(c, r, f'{val:.2f}',
                            ha='center', va='center', fontsize=7.5,
                            color=txt_color, fontweight='bold')

    # ── Single shared colorbar in its own axes ────────────────────────────────
    cbar_left   = ax_left + n_cli * (ax_width + h_gap) + 0.02
    cbar_bottom = ax_bottom
    cbar_height = n_reg * ax_height + (n_reg - 1) * v_gap
    cbar_ax     = fig.add_axes([cbar_left, cbar_bottom,
                                cbar_w / fig_w, cbar_height])
    cbar = fig.colorbar(im_ref, cax=cbar_ax)
    cbar.set_label('Conflict Index\n(0 = peace,  3 = fragmentation)', fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    # Add regime boundary ticks
    for thr in _CI_THRESHOLDS:
        if global_vmin <= thr <= global_vmax:
            cbar.ax.axhline((thr - global_vmin) / (global_vmax - global_vmin),
                            color='white', lw=1.0, ls='--')

    if save:
        fname = os.path.join(VIS_DIR, 'combined_ci_matrix.png')
        fig.savefig(fname, dpi=150, bbox_inches='tight')
        print(f'  Saved: {fname}')
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# PLOT 3: Reaction group activity profile  (replaces budget utilisation)
#
# Design rationale
# ----------------
# Budget utilisation (idle-fraction fill plots) was hard to read because it
# mixed scale effects with allocation effects.  Here we instead show the
# *total volume of events* per reaction group accumulated over the game,
# expressed in events / year.  This directly answers "which processes are
# activated by each strategy choice?" without the per-species decomposition.
#
# Layout: one row per Gov strategy (4 rows), x-axis = AG strategies (4 bars),
# stacked bars coloured by reaction group.  A second thin bar for each combo
# overlays the Natural-only baseline (Gov=None, AG=None) as an outline,
# making it easy to see what strategy adds on top of the background.
# ─────────────────────────────────────────────────────────────────────────────

def _group_activity(log: List[dict]) -> Dict[str, float]:
    """
    Total flux (events/year) per reaction group, averaged over all periods.
    """
    total_T = sum(e['T_period'] for e in log)
    if total_T == 0:
        return {grp: 0.0 for grp in REACTION_GROUPS}
    activity = {}
    for grp, rxns in REACTION_GROUPS.items():
        cumul = 0.0
        for entry in log:
            flux = entry['flux']
            cols = [r for r in rxns if r in flux.columns]
            if cols:
                cumul += float(flux[cols].values.sum())
        activity[grp] = cumul / total_T
    return activity


def plot_reaction_activity(all_logs: dict, climate: float, regime: int,
                            gov_strats: list, ag_strats: list, save: bool):
    """
    Reaction group activity profile: stacked bar chart.

    Rows    = Gov strategy
    x-axis  = AG strategy
    y-axis  = mean events / year accumulated across all periods
    Stacked segments = reaction groups
    Grey outline bar = None×None baseline (natural background activity)
    """
    groups =  ['Economy',
    'Migration',
    'Resilience',   
    'Trust',
    'Violence',
    'ArmedGroup']
    n_gov  = len(gov_strats)
    n_ag   = len(ag_strats)
    x      = np.arange(n_ag)
    width  = 0.55

    # Baseline (None × None) activity for reference outline
    baseline_log = all_logs.get(('None', 'None'), [])
    baseline_act = _group_activity(baseline_log) if baseline_log else {}
    baseline_total = [sum(baseline_act.get(g, 0) for g in groups)] * n_ag

    fig, axes = plt.subplots(n_gov, 1,
                             figsize=(max(6, 2 * n_ag), 2.6 * n_gov),
                             sharex=False)
    if n_gov == 1:
        axes = [axes]
    fig.suptitle(
        f'Reaction group activity (events/yr)  |  R{regime}  climate={climate:.2f}',
        fontsize=11)

    # Collect global y max for consistent y-axes
    all_totals = []
    for gs in gov_strats:
        for ag in ag_strats:
            act = _group_activity(all_logs.get((gs, ag), []))
            all_totals.append(sum(act.values()))
    y_max = max(all_totals) * 1.12 if all_totals else 1.0

    for row, gs in enumerate(gov_strats):
        ax = axes[row]
        bottoms = np.zeros(n_ag)

        for grp in groups:
            heights = []
            for ag in ag_strats:
                act = _group_activity(all_logs.get((gs, ag), []))
                heights.append(act.get(grp, 0.0))
            ax.bar(x, heights, width, bottom=bottoms,
                   label=grp, color=GROUP_COLORS[grp],
                   edgecolor='white', linewidth=0.4)
            bottoms += np.array(heights)

        # Overlay baseline total as a step-outline (None×None reference)
        for xi, bt in enumerate(baseline_total):
            ax.plot([xi - width / 2, xi + width / 2], [bt, bt],
                    color='black', lw=1.4, ls='--')

        ax.set_ylim(0, y_max)
        ax.set_ylabel(f'Gov: {gs}\nevents/yr', fontsize=8)
        ax.set_xticks(x)
        ax.set_xticklabels([f'AG: {SHORT[ag]}' for ag in ag_strats], fontsize=9)
        ax.grid(axis='y', lw=0.4, alpha=0.4)
        ax.tick_params(labelsize=7)
        if row == 0:
            ax.legend(loc='upper right', fontsize=7, ncol=len(groups))

    # Dashed line legend entry
    axes[0].plot([], [], color='black', lw=1.4, ls='--', label='None×None baseline')
    axes[0].legend(loc='upper right', fontsize=7,
                   ncol=max(1, (len(groups) + 2) // 2))

    axes[-1].set_xlabel('AG strategy', fontsize=9)
    plt.tight_layout()
    if save:
        fname = os.path.join(VIS_DIR,
                             f'reaction_activity_R{regime}_c{climate:.2f}.png')
        fig.savefig(fname, dpi=130, bbox_inches='tight')
        print(f'  Saved: {fname}')
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# PLOT 4: Regime evolution heatmap
# Rows = all 16 strategy combos; columns = period; colour = classified regime.
# ─────────────────────────────────────────────────────────────────────────────

def plot_regime_evolution(all_logs: dict, climate: float, regime: int,
                           gov_strats: list, ag_strats: list, save: bool):
    """
    Heatmap: classified regime at the end of each period for every combo.
    """
    combos = [(gs, ag) for gs in gov_strats for ag in ag_strats]
    labels = [f'{SHORT[gs]}/{SHORT[ag]}' for gs, ag in combos]
    n_periods = max((len(all_logs.get(c, [])) for c in combos), default=1)

    data = np.full((len(combos), n_periods), np.nan)
    for i, combo in enumerate(combos):
        for p, entry in enumerate(all_logs.get(combo, [])):
            data[i, p] = entry['metrics']['regime']

    fig_h = max(5, len(combos) * 0.38)
    fig_w = max(7, n_periods * 0.7)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    im = ax.imshow(data, aspect='auto', origin='upper',
                   vmin=1, vmax=6, cmap='RdYlGn_r',
                   extent=[0.5, n_periods + 0.5,
                           len(combos) + 0.5, 0.5])

    ax.set_yticks(range(1, len(combos) + 1))
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xticks(range(1, n_periods + 1))
    ax.set_xticklabels([f'P{p}' for p in range(1, n_periods + 1)], fontsize=8)
    ax.set_xlabel('Period', fontsize=9)
    ax.set_title(
        f'Regime evolution by strategy combo  |  starting R{regime}  climate={climate:.2f}',
        fontsize=10)

    # Annotate regime number in each cell
    for i in range(len(combos)):
        for p in range(n_periods):
            if not np.isnan(data[i, p]):
                ax.text(p + 1, i + 1, f'R{int(data[i, p])}',
                        ha='center', va='center', fontsize=6.5, color='white',
                        fontweight='bold')

    cbar = fig.colorbar(im, ax=ax, ticks=range(1, 7), fraction=0.02, pad=0.01)
    cbar.ax.set_yticklabels(
        [f'R{i}: {REGIME_NAMES[i][:16]}' for i in range(1, 7)], fontsize=6.5)

    # Horizontal dividers between Gov strategy blocks
    n_ag = len(ag_strats)
    for i in range(n_ag, len(combos), n_ag):
        ax.axhline(i + 0.5, color='white', lw=1.5)

    plt.tight_layout()
    if save:
        fname = os.path.join(VIS_DIR,
                             f'regime_evolution_R{regime}_c{climate:.2f}.png')
        fig.savefig(fname, dpi=130, bbox_inches='tight')
        print(f'  Saved: {fname}')
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# PLOT 5: Adaptive game — CI trajectory + strategy timeline + indicator panels
# ─────────────────────────────────────────────────────────────────────────────

def plot_adaptive_game(log: List[dict], climate: float, regime: int,
                       gov_rule: str, ag_rule: str, save: bool):
    """
    Adaptive game visualisation for a single (gov_rule x ag_rule) run.

    Layout (3 rows):
      Row 0 (tall) : Conflict Index trajectory over time.
                     Two thin coloured strips at the top of the panel encode
                     Gov (upper strip) and AG (lower strip) active strategy
                     per period.  Period backgrounds are tinted by end-regime.
      Row 1        : war_economy, poverty, social_cohesion, ag_power
      Row 2        : militarisation, ag_infiltration, displacement, resilience_ratio

    How to read the strategy strips
    --------------------------------
    Look at the top ~14 % of the CI panel.  The upper band (Gov) and lower
    band (AG) change colour whenever the actor switches strategy.  Colour key:
      Gov: grey=None  green=Development  blue=Defense  purple=Balance
      AG:  grey=None  red=Coerce         orange=Extract dark-red=Balance
    A grey strip means the actor failed the capacity check that period.
    """
    import matplotlib.patches as mpatches

    fig = plt.figure(figsize=(14, 10))
    gs_outer = fig.add_gridspec(
        3, 1, height_ratios=[3, 2, 2], hspace=0.55)

    # ── Row 0: CI + strategy timeline ────────────────────────────────────────
    ax_ci = fig.add_subplot(gs_outer[0])

    t_all, ci_all = [], []
    for entry in log:
        t_local = entry['traj'].index.values + entry['t_offset']
        t_all.extend(t_local.tolist())
        ci_all.extend(_ci_from_traj(entry['traj']).tolist())

    ax_ci.plot(t_all, ci_all, color='#2c3e50', lw=1.4, zorder=5)

    # CI threshold lines + regime zone labels
    _ci_bounds = [0.0] + list(_CI_THRESHOLDS) + [3.0]
    for thr in _CI_THRESHOLDS:
        ax_ci.axhline(thr, color='#888888', lw=0.6, ls=':', alpha=0.7, zorder=2)
    for zi in range(6):
        zmid = (_ci_bounds[zi] + _ci_bounds[zi + 1]) / 2
        ax_ci.annotate(f'R{zi+1}', xy=(1.01, zmid),
                       xycoords=('axes fraction', 'data'),
                       fontsize=7, va='center', ha='left', color='#555555')

    # Per-period decorations
    for entry in log:
        t0    = entry['t_offset']
        t1    = t0 + entry['T_period']
        r_end = entry['metrics']['regime']
        # Faint regime background (lower 86 % of axis height)
        ax_ci.axvspan(t0, t1, ymin=0.0, ymax=0.86,
                      color=REGIME_COLORS[r_end], alpha=0.09, zorder=1)
        # Gov strategy strip (upper 8 %)
        ax_ci.axvspan(t0, t1, ymin=0.92, ymax=1.0,
                      color=STRATEGY_COLORS['Gov'][entry['gov_strat']],
                      alpha=0.85, zorder=3)
        # AG strategy strip (next 6 %)
        ax_ci.axvspan(t0, t1, ymin=0.86, ymax=0.92,
                      color=STRATEGY_COLORS['AG'][entry['ag_strat']],
                      alpha=0.85, zorder=3)
        # Period boundary line
        if entry['period'] > 1:
            ax_ci.axvline(t0, color='#bbbbbb', lw=0.6, ls='--', zorder=4)

    ax_ci.set_ylim(0, 3.05)
    ax_ci.set_ylabel('Conflict Index', fontsize=9)
    ax_ci.set_title(
        f'Adaptive Game  |  R{regime}  climate={climate:.2f}  '
        f'Gov-rule: {gov_rule}  AG-rule: {ag_rule}',
        fontsize=10)
    ax_ci.tick_params(labelsize=7)
    ax_ci.grid(True, axis='y', lw=0.3, alpha=0.35)

    # Strip role labels
    ax_ci.text(0.002, 0.96, 'Gov', transform=ax_ci.transAxes,
               fontsize=6, va='center', color='white', fontweight='bold')
    ax_ci.text(0.002, 0.89, 'AG',  transform=ax_ci.transAxes,
               fontsize=6, va='center', color='white', fontweight='bold')

    # Strategy colour legend
    gov_patches = [mpatches.Patch(color=c, label=f'G:{s}')
                   for s, c in STRATEGY_COLORS['Gov'].items()]
    ag_patches  = [mpatches.Patch(color=c, label=f'A:{s}')
                   for s, c in STRATEGY_COLORS['AG'].items()]
    ax_ci.legend(handles=gov_patches + ag_patches,
                 loc='upper left', fontsize=6.5, ncol=4,
                 bbox_to_anchor=(0.0, 0.86), framealpha=0.75)

    # ── Rows 1-2: indicator panels ───────────────────────────────────────────
    for row_i in range(2):
        gs_inner = gs_outer[1 + row_i].subgridspec(1, 4, wspace=0.38)
        for col_i in range(4):
            ind_key = INDICATOR_NAMES[row_i * 4 + col_i]
            ax_i = fig.add_subplot(gs_inner[col_i])
            for entry in log:
                t_local = entry['traj'].index.values + entry['t_offset']
                vals    = _ind_from_traj(entry['traj'], ind_key)
                ax_i.plot(t_local, vals, color='#2c3e50', lw=0.9)
                # Regime-tinted background per period
                ax_i.axvspan(t_local[0], t_local[-1],
                             color=REGIME_COLORS[entry['metrics']['regime']],
                             alpha=0.08)
            ax_i.set_ylim(0, 1.05)
            ax_i.set_ylabel(INDICATOR_LABELS[ind_key], fontsize=7)
            ax_i.tick_params(labelsize=6)
            ax_i.grid(True, lw=0.3, alpha=0.35)
            if row_i == 1:
                ax_i.set_xlabel('years', fontsize=6)

    if save:
        fname = os.path.join(
            VIS_DIR,
            f'adaptive_R{regime}_c{climate:.2f}_{gov_rule}_{ag_rule}.png')
        fig.savefig(fname, dpi=130, bbox_inches='tight')
        print(f'  Saved: {fname}')
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# PLOT 9: Rule marginal performance
# Mean final CI per rule averaged over all opposing rules, regimes, climates.
# ─────────────────────────────────────────────────────────────────────────────

def plot_rule_marginal_performance(results_by_rule_rc, gov_rules, ag_rules,
                                   save=True):
    """
    Two-panel bar chart: left = Gov rules, right = AG rules.
    Each bar = mean final CI averaged over all opposing rules, regimes, climates.
    Bars sorted ascending by mean CI; error bars = ±1 std.
    Saved as adaptive_rule_marginal_performance.png.
    """
    gov_cis = {gr: [] for gr in gov_rules}
    ag_cis  = {ar: [] for ar in ag_rules}
    for (_regime, _climate), results in results_by_rule_rc.items():
        for gr in gov_rules:
            for ar in ag_rules:
                log = results.get((gr, ar), [])
                if log:
                    ci = log[-1]['indicators']['conflict_index']
                    gov_cis[gr].append(ci)
                    ag_cis[ar].append(ci)

    fig, (ax_gov, ax_ag) = plt.subplots(1, 2, figsize=(10, 4.5))
    fig.suptitle(
        'Rule marginal performance — mean final CI across all combinations\n'
        '(lower CI = better for peace)',
        fontsize=11)

    for ax, rules, cis_dict, title in [
        (ax_gov, gov_rules, gov_cis, 'Government rules'),
        (ax_ag,  ag_rules,  ag_cis,  'Armed Group rules'),
    ]:
        shorts = [RULE_SHORT.get(r, r) for r in rules]
        means  = np.array([np.mean(cis_dict[r]) if cis_dict[r] else 0.0
                           for r in rules])
        stds   = np.array([np.std(cis_dict[r])  if cis_dict[r] else 0.0
                           for r in rules])
        order  = np.argsort(means)
        s_shorts = [shorts[i] for i in order]
        s_means  = means[order]
        s_stds   = stds[order]
        n = len(rules)
        colors = plt.cm.RdYlGn_r(np.linspace(0.2, 0.8, n))
        bars = ax.bar(range(n), s_means, color=colors, edgecolor='white', lw=0.5,
                      yerr=s_stds, capsize=5,
                      error_kw={'elinewidth': 1.2, 'ecolor': '#555555'})
        for bar, m, s in zip(bars, s_means, s_stds):
            ax.text(bar.get_x() + bar.get_width() / 2, m + s + 0.02,
                    f'{m:.2f}', ha='center', va='bottom', fontsize=8.5)
        ax.set_xticks(range(n))
        ax.set_xticklabels(s_shorts, fontsize=9)
        ax.set_ylabel('Mean final CI', fontsize=9)
        ax.set_title(title, fontsize=10)
        y_top = min(3.0, float((s_means + s_stds).max()) + 0.35) if n else 3.0
        ax.set_ylim(0, y_top)
        for thr in _CI_THRESHOLDS:
            if thr <= y_top:
                ax.axhline(thr, color='#cccccc', lw=0.6, ls=':', zorder=1)
        ax.grid(axis='y', lw=0.3, alpha=0.4)
        ax.tick_params(labelsize=8)

    plt.tight_layout()
    if save:
        fname = os.path.join(VIS_DIR, 'adaptive_rule_marginal_performance.png')
        fig.savefig(fname, dpi=130, bbox_inches='tight')
        print(f'  Saved: {fname}')
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# PLOT 10: Best vs worst trajectories
# CI over time for top-n and bottom-n rule pairs per (regime × climate).
# ─────────────────────────────────────────────────────────────────────────────

def plot_best_worst_trajectories(results_by_rule_rc, regimes, climates,
                                  gov_rules, ag_rules, n_top=3, save=True):
    """
    One figure per regime; columns = climates.
    Top-n pairs (lowest final CI) drawn in cool colors.
    Bottom-n pairs (highest final CI) drawn in warm colors.
    Saved as adaptive_best_worst_traj_R{regime}.png.
    """
    COOL_COLORS = ['#1abc9c', '#2980b9', '#8e44ad']
    WARM_COLORS = ['#e74c3c', '#e67e22', '#f39c12']
    _ci_bounds = [0.0] + list(_CI_THRESHOLDS) + [3.0]

    for regime in regimes:
        climates_here = [c for c in climates if (regime, c) in results_by_rule_rc]
        if not climates_here:
            continue
        n_cli = len(climates_here)
        fig, axes = plt.subplots(1, n_cli, figsize=(5.5 * n_cli, 5),
                                 sharey=True)
        if n_cli == 1:
            axes = [axes]
        fig.suptitle(
            f'Best vs worst rule pairs — CI trajectory  |  R{regime}\n'
            f'Cool colors = best (low CI)   Warm colors = worst (high CI)',
            fontsize=10)

        for ci_col, climate in enumerate(climates_here):
            ax = axes[ci_col]
            results = results_by_rule_rc[(regime, climate)]
            ranked = sorted(
                [(gr, ar, results[(gr, ar)][-1]['indicators']['conflict_index'])
                 for gr in gov_rules for ar in ag_rules
                 if results.get((gr, ar))],
                key=lambda x: x[2])
            best  = ranked[:n_top]
            worst = ranked[-n_top:]

            for thr in _CI_THRESHOLDS:
                ax.axhline(thr, color='#dddddd', lw=0.5, ls=':', zorder=1)
            for zi in range(6):
                zmid = (_ci_bounds[zi] + _ci_bounds[zi + 1]) / 2
                ax.annotate(f'R{zi+1}', xy=(1.01, zmid),
                            xycoords=('axes fraction', 'data'),
                            fontsize=6.5, va='center', ha='left', color='#888888')

            for pairs, colors in [(best, COOL_COLORS), (worst, WARM_COLORS)]:
                for (gr, ar, ci_end), color in zip(pairs, colors):
                    log = results[(gr, ar)]
                    t_all, ci_all = [], []
                    for entry in log:
                        t_local = entry['traj'].index.values + entry['t_offset']
                        t_all.extend(t_local.tolist())
                        ci_all.extend(_ci_from_traj(entry['traj']).tolist())
                    lbl = (f'{RULE_SHORT.get(gr, gr)}/{RULE_SHORT.get(ar, ar)}'
                           f'  CI={ci_end:.2f}')
                    ax.plot(t_all, ci_all, color=color, lw=1.8, label=lbl, zorder=5)
                    for entry in log[1:]:
                        ax.axvline(entry['t_offset'],
                                   color='#dddddd', lw=0.4, ls='--', zorder=1)

            ax.set_ylim(0, 3.05)
            ax.set_title(f'climate={climate:.2f}', fontsize=9)
            ax.set_xlabel('years', fontsize=8)
            if ci_col == 0:
                ax.set_ylabel('Conflict Index', fontsize=9)
            ax.legend(fontsize=7, loc='upper left')
            ax.tick_params(labelsize=7)
            ax.grid(True, axis='y', lw=0.3, alpha=0.35)

        plt.tight_layout()
        if save:
            fname = os.path.join(VIS_DIR, f'adaptive_best_worst_traj_R{regime}.png')
            fig.savefig(fname, dpi=130, bbox_inches='tight')
            print(f'  Saved: {fname}')
        plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# PLOT 11: Strategy usage per rule pair
# Stacked bar of Gov strategy fraction, pairs sorted by final CI.
# ─────────────────────────────────────────────────────────────────────────────

def plot_strategy_usage(results_by_rule_rc, regimes, climates,
                         gov_rules, ag_rules, save=True):
    """
    Per regime: stacked bar of Gov strategy fraction per rule pair.
    Pairs sorted by ascending final CI (best→worst left→right).
    One figure per regime with one column per climate.
    Saved as adaptive_strategy_usage_R{regime}.png.
    """
    GOV_STRAT_ORDER  = ['Development', 'Balance', 'Defense', 'None']
    GOV_STRAT_COLORS = {
        'Development': '#27ae60',
        'Balance':     '#8e44ad',
        'Defense':     '#2980b9',
        'None':        '#aaaaaa',
    }

    for regime in regimes:
        climates_here = [c for c in climates if (regime, c) in results_by_rule_rc]
        if not climates_here:
            continue
        n_cli = len(climates_here)
        fig, axes = plt.subplots(1, n_cli, figsize=(max(7, 6 * n_cli), 4.5),
                                 sharey=True)
        if n_cli == 1:
            axes = [axes]
        fig.suptitle(
            f'Gov strategy usage per rule pair  |  R{regime}\n'
            f'Pairs sorted by final CI (best=lowest, left to right)'
            f'   — annotated CI value above each bar',
            fontsize=10)

        for ci_col, climate in enumerate(climates_here):
            ax = axes[ci_col]
            results = results_by_rule_rc[(regime, climate)]
            ranked = sorted(
                [(gr, ar, results[(gr, ar)][-1]['indicators']['conflict_index'])
                 for gr in gov_rules for ar in ag_rules
                 if results.get((gr, ar))],
                key=lambda x: x[2])
            pair_labels = [f'{RULE_SHORT.get(gr, gr)}/{RULE_SHORT.get(ar, ar)}'
                           for gr, ar, _ in ranked]
            n_pairs = len(ranked)
            x = np.arange(n_pairs)
            bottoms = np.zeros(n_pairs)

            for strat in GOV_STRAT_ORDER:
                fracs = []
                for gr, ar, _ in ranked:
                    log = results[(gr, ar)]
                    n_periods = max(len(log), 1)
                    fracs.append(
                        sum(1 for e in log if e['gov_strat'] == strat) / n_periods)
                fracs = np.array(fracs)
                ax.bar(x, fracs, bottom=bottoms,
                       color=GOV_STRAT_COLORS[strat], label=strat,
                       edgecolor='white', lw=0.4)
                bottoms += fracs

            # Final CI annotation
            for i, (_gr, _ar, ci_end) in enumerate(ranked):
                ax.text(i, 1.02, f'{ci_end:.2f}',
                        ha='center', va='bottom', fontsize=6.5, color='#333333')

            ax.set_xticks(x)
            ax.set_xticklabels(pair_labels, rotation=45, ha='right', fontsize=7)
            ax.set_title(f'climate={climate:.2f}', fontsize=9)
            if ci_col == 0:
                ax.set_ylabel('Gov strategy fraction', fontsize=9)
                ax.legend(loc='upper left', fontsize=7.5)
            ax.set_ylim(0, 1.12)
            ax.grid(axis='y', lw=0.3, alpha=0.4)
            ax.tick_params(labelsize=7)

        plt.tight_layout()
        if save:
            fname = os.path.join(VIS_DIR,
                                 f'adaptive_strategy_usage_R{regime}.png')
            fig.savefig(fname, dpi=130, bbox_inches='tight')
            print(f'  Saved: {fname}')
        plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# PLOT 12: Capacity robustness
# Final CI vs gov-capacity threshold, bars grouped by climate (per regime).
# ─────────────────────────────────────────────────────────────────────────────

def plot_capacity_robustness(logs_by_thresh_climate, thresholds, climates,
                              gov_rule, ag_rule, regime, save=True):
    """
    Bar chart: x = capacity threshold, bars grouped by climate.
    Annotated with final CI value above each bar.
    Saved as capacity_robustness_R{regime}_{gov_rule}_{ag_rule}.png.
    """
    CLIMATE_COLORS = ['#2980b9', '#27ae60', '#e74c3c', '#8e44ad']
    n_cli = len(climates)
    n_thr = len(thresholds)
    bar_w = 0.7 / max(n_cli, 1)
    x = np.arange(n_thr)

    fig, ax = plt.subplots(figsize=(max(6, 2.5 * n_thr), 4.5))
    fig.suptitle(
        f'Capacity robustness  |  R{regime}  '
        f'Gov: {gov_rule}   AG: {ag_rule}\n'
        f'Final CI by capacity threshold and climate  (lower = better for peace)',
        fontsize=10)

    for ci_i, climate in enumerate(climates):
        color = CLIMATE_COLORS[ci_i % len(CLIMATE_COLORS)]
        ci_vals = []
        for thresh in thresholds:
            log = logs_by_thresh_climate.get((climate, thresh), [])
            ci_vals.append(
                log[-1]['indicators']['conflict_index'] if log else np.nan)
        offset = (ci_i - (n_cli - 1) / 2.0) * bar_w
        bars = ax.bar(x + offset, ci_vals, bar_w * 0.88,
                      color=color, label=f'climate={climate:.2f}',
                      edgecolor='white', lw=0.5)
        for bar, val in zip(bars, ci_vals):
            if not np.isnan(val):
                ax.text(bar.get_x() + bar.get_width() / 2, val + 0.02,
                        f'{val:.2f}', ha='center', va='bottom', fontsize=8)

    for thr in _CI_THRESHOLDS:
        ax.axhline(thr, color='#dddddd', lw=0.6, ls=':', zorder=1)

    thresh_labels = ['det.\n(thresh=0)' if t == 0 else f'thresh={t:.0f}'
                     for t in thresholds]
    ax.set_xticks(x)
    ax.set_xticklabels(thresh_labels, fontsize=9)
    ax.set_ylabel('Final Conflict Index', fontsize=9)
    ax.set_ylim(0, 3.05)
    ax.legend(fontsize=8, loc='upper left')
    ax.grid(axis='y', lw=0.3, alpha=0.4)
    ax.tick_params(labelsize=8)

    plt.tight_layout()
    if save:
        fname = os.path.join(
            VIS_DIR,
            f'capacity_robustness_R{regime}_{gov_rule}_{ag_rule}.png')
        fig.savefig(fname, dpi=130, bbox_inches='tight')
        print(f'  Saved: {fname}')
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# REPORT WRITERS  (adaptive mode)
# ─────────────────────────────────────────────────────────────────────────────

def save_adaptive_detail_report(log: List[dict], gov_rule: str, ag_rule: str,
                                 regime: int, climate: float,
                                 gov_cap: float, ag_cap: float,
                                 gov_pm: int, ag_pm: int) -> pd.DataFrame:
    """
    Per-period CSV for a single adaptive run.

    Columns: regime_start, climate, gov_rule, ag_rule, gov_cap_thresh,
             ag_cap_thresh, gov_period_mult, ag_period_mult,
             period, gov_strat, ag_strat, CI, conflict_regime,
             war_economy, poverty, social_cohesion, ag_power,
             militarisation, ag_infiltration, displacement, resilience_ratio.
    """
    rows = []
    for entry in log:
        ind = entry.get('indicators', {})
        row = {
            'regime_start':    regime,
            'climate':         climate,
            'gov_rule':        gov_rule,
            'ag_rule':         ag_rule,
            'gov_cap_thresh':  gov_cap,
            'ag_cap_thresh':   ag_cap,
            'gov_period_mult': gov_pm,
            'ag_period_mult':  ag_pm,
            'period':          entry['period'],
            'gov_strat':       entry['gov_strat'],
            'ag_strat':        entry['ag_strat'],
            'CI':              ind.get('conflict_index',
                                      entry['metrics'].get('conflict_index', float('nan'))),
            'conflict_regime': entry['metrics']['regime'],
        }
        for k in INDICATOR_NAMES:
            row[k] = ind.get(k, float('nan'))
        rows.append(row)
    df = pd.DataFrame(rows)
    fname = os.path.join(
        REPORT_DIR,
        f'detail_R{regime}_c{climate:.2f}_{gov_rule}_{ag_rule}'
        f'_gcap{gov_cap:.0f}_gpm{gov_pm}.csv')
    df.to_csv(fname, index=False)
    print(f'  Report: {os.path.basename(fname)}')
    return df


def save_adaptive_grid_report(results_by_rule: dict, gov_rules: list, ag_rules: list,
                               regime: int, climate: float,
                               gov_cap: float, ag_cap: float,
                               gov_pm: int, ag_pm: int) -> pd.DataFrame:
    """
    Summary CSV for a full (gov_rule × ag_rule) grid run.

    One row per rule pair. Columns: gov_rule, ag_rule, final_CI, final_regime,
    min_CI, max_CI, mean_CI, gov_strategy_switches, ag_strategy_switches,
    gov_none_periods (periods where Gov defaulted to None due to capacity).
    Sorted ascending by final_CI (best pair first).
    """
    rows = []
    for gr in gov_rules:
        for ar in ag_rules:
            log = results_by_rule.get((gr, ar), [])
            if not log:
                continue
            ci_vals   = [e['indicators']['conflict_index'] for e in log]
            last      = log[-1]
            gov_sw    = sum(1 for i in range(1, len(log))
                            if log[i]['gov_strat'] != log[i-1]['gov_strat'])
            ag_sw     = sum(1 for i in range(1, len(log))
                            if log[i]['ag_strat']  != log[i-1]['ag_strat'])
            gov_none  = sum(1 for e in log if e['gov_strat'] == 'None')
            rows.append({
                'regime_start':          regime,
                'climate':               climate,
                'gov_rule':              gr,
                'ag_rule':               ar,
                'gov_cap_thresh':        gov_cap,
                'ag_cap_thresh':         ag_cap,
                'gov_period_mult':       gov_pm,
                'ag_period_mult':        ag_pm,
                'final_CI':              last['indicators']['conflict_index'],
                'final_regime':          last['metrics']['regime'],
                'min_CI':                min(ci_vals),
                'max_CI':                max(ci_vals),
                'mean_CI':               sum(ci_vals) / len(ci_vals),
                'gov_strategy_switches': gov_sw,
                'ag_strategy_switches':  ag_sw,
                'gov_none_periods':      gov_none,
            })
    df = pd.DataFrame(rows).sort_values('final_CI').reset_index(drop=True)
    fname = os.path.join(
        REPORT_DIR,
        f'grid_R{regime}_c{climate:.2f}_gcap{gov_cap:.0f}_gpm{gov_pm}.csv')
    df.to_csv(fname, index=False)
    print(f'  Report: {os.path.basename(fname)}')
    return df


def save_capacity_sweep_report(logs_by_thresh: dict, thresholds: list,
                                gov_rule: str, ag_rule: str,
                                regime: int, climate: float) -> pd.DataFrame:
    """
    Summary CSV for a capacity-sweep run.  One row per (threshold × period).
    """
    rows = []
    for thresh in thresholds:
        log = logs_by_thresh.get(thresh, [])
        for entry in log:
            ind = entry.get('indicators', {})
            rows.append({
                'regime_start':   regime,
                'climate':        climate,
                'gov_rule':       gov_rule,
                'ag_rule':        ag_rule,
                'gov_cap_thresh': thresh,
                'period':         entry['period'],
                'gov_strat':      entry['gov_strat'],
                'ag_strat':       entry['ag_strat'],
                'CI':             ind.get('conflict_index', float('nan')),
                'conflict_regime':entry['metrics']['regime'],
            })
    df = pd.DataFrame(rows)
    fname = os.path.join(
        REPORT_DIR,
        f'capacity_sweep_R{regime}_c{climate:.2f}_{gov_rule}_{ag_rule}.csv')
    df.to_csv(fname, index=False)
    print(f'  Report: {os.path.basename(fname)}')
    return df


def save_freq_sweep_report(logs_by_mult: dict, mults: list,
                            gov_rule: str, ag_rule: str,
                            regime: int, climate: float) -> pd.DataFrame:
    """
    Summary CSV for a frequency-sweep run.  One row per (mult × period).
    """
    rows = []
    for mult in mults:
        log = logs_by_mult.get(mult, [])
        for entry in log:
            ind = entry.get('indicators', {})
            rows.append({
                'regime_start':    regime,
                'climate':         climate,
                'gov_rule':        gov_rule,
                'ag_rule':         ag_rule,
                'gov_period_mult': mult,
                'period':          entry['period'],
                'gov_strat':       entry['gov_strat'],
                'ag_strat':        entry['ag_strat'],
                'CI':              ind.get('conflict_index', float('nan')),
                'conflict_regime': entry['metrics']['regime'],
            })
    df = pd.DataFrame(rows)
    fname = os.path.join(
        REPORT_DIR,
        f'freq_sweep_R{regime}_c{climate:.2f}_{gov_rule}_{ag_rule}.csv')
    df.to_csv(fname, index=False)
    print(f'  Report: {os.path.basename(fname)}')
    return df


# ─────────────────────────────────────────────────────────────────────────────
# PLOT 6: Adaptive rule matrix
# 4×4 heatmap of final CI by (gov_rule × ag_rule).
# Multiple panels for different (regime, climate) combinations.
# ─────────────────────────────────────────────────────────────────────────────

def plot_adaptive_rule_matrix(
        results_by_rule_rc: dict,   # {(regime, climate): {(gov_rule, ag_rule): log}}
        regimes: list,
        climates: list,
        gov_rules: list,
        ag_rules: list,
        save: bool):
    """
    One 4×4 heatmap per (regime, climate) combination.
    Rows = Gov rules, cols = AG rules.
    Color + annotation = final CI (absolute, shared scale across all panels).
    Layout mirrors plot_combined_ci_matrix: rows=regimes, cols=climates.
    """
    n_reg = len(regimes)
    n_cli = len(climates)
    n_gov = len(gov_rules)
    n_ag  = len(ag_rules)

    # Collect CI grids
    ci_grids = {}
    for regime in regimes:
        for climate in climates:
            results = results_by_rule_rc.get((regime, climate), {})
            m = np.zeros((n_gov, n_ag))
            for ri, gr in enumerate(gov_rules):
                for ci_i, ar in enumerate(ag_rules):
                    log = results.get((gr, ar), [])
                    if log:
                        m[ri, ci_i] = log[-1]['indicators']['conflict_index']
            ci_grids[(regime, climate)] = m

    all_vals    = np.concatenate([m.ravel() for m in ci_grids.values()])
    global_vmin = max(0.0, float(all_vals.min()) - 0.05)
    global_vmax = min(3.0, float(all_vals.max()) + 0.05)

    cell_w, cell_h = 3.6, 3.6
    cbar_w  = 0.5
    left_m  = 0.85
    right_m = cbar_w + 0.5
    top_m   = 0.75
    bot_m   = 0.55

    fig_w = left_m + cell_w * n_cli + right_m
    fig_h = top_m  + cell_h * n_reg + bot_m

    fig = plt.figure(figsize=(fig_w, fig_h))
    fig.suptitle(
        'Adaptive rule matrix — final Conflict Index (shared scale)\n'
        'Rows = Gov rules   |   Cols = AG rules   |   '
        'Green = low CI (better)   Red = high CI (worse)',
        fontsize=10,
    )

    ax_left   = left_m / fig_w
    ax_bottom = bot_m  / fig_h
    ax_width  = cell_w / fig_w
    ax_height = cell_h / fig_h
    h_gap     = 0.08 / fig_w
    v_gap     = 0.55 / fig_h

    im_ref = None
    for ri, regime in enumerate(regimes):
        y0 = ax_bottom + (n_reg - 1 - ri) * (ax_height + v_gap)
        for ci_col, climate in enumerate(climates):
            x0 = ax_left + ci_col * (ax_width + h_gap)
            ax = fig.add_axes([x0, y0, ax_width, ax_height])
            m  = ci_grids[(regime, climate)]

            im = ax.imshow(m, cmap='RdYlGn_r',
                           vmin=global_vmin, vmax=global_vmax, aspect='auto')
            im_ref = im

            ax.set_xticks(range(n_ag))
            ax.set_yticks(range(n_gov))
            ax.set_xticklabels([RULE_SHORT.get(ar, ar) for ar in ag_rules],
                               rotation=30, ha='right', fontsize=7)
            ax.set_yticklabels([RULE_SHORT.get(gr, gr) for gr in gov_rules], fontsize=7)

            if ci_col == 0:
                ax.set_ylabel(f'R{regime}\nGov rule', fontsize=8)
            else:
                ax.set_yticklabels([])
            if ri == n_reg - 1:
                ax.set_xlabel('AG rule', fontsize=8)
            else:
                ax.set_xticklabels([])

            ax.set_title(f'R{regime}  climate={climate:.2f}', fontsize=8.5, pad=3)

            for r in range(n_gov):
                for c in range(n_ag):
                    val = m[r, c]
                    norm_val = (val - global_vmin) / max(global_vmax - global_vmin, 1e-9)
                    txt_color = 'white' if norm_val > 0.82 else 'black'
                    ax.text(c, r, f'{val:.2f}',
                            ha='center', va='center', fontsize=7.5,
                            color=txt_color, fontweight='bold')

    # Shared colorbar
    cbar_left   = ax_left + n_cli * (ax_width + h_gap) + 0.02
    cbar_bottom = ax_bottom
    cbar_height = n_reg * ax_height + (n_reg - 1) * v_gap
    cbar_ax     = fig.add_axes([cbar_left, cbar_bottom, cbar_w / fig_w, cbar_height])
    cbar = fig.colorbar(im_ref, cax=cbar_ax)
    cbar.set_label('Conflict Index\n(0 = peace,  3 = fragmentation)', fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    for thr in _CI_THRESHOLDS:
        if global_vmin <= thr <= global_vmax:
            norm = (thr - global_vmin) / max(global_vmax - global_vmin, 1e-9)
            cbar.ax.axhline(norm, color='white', lw=0.8, ls='--')

    if save:
        regs_str = ''.join(str(r) for r in regimes)
        cli_str  = '_'.join(f'{c:.2f}' for c in climates)
        fname = os.path.join(VIS_DIR,
                             f'adaptive_rule_matrix_R{regs_str}_c{cli_str}.png')
        fig.savefig(fname, dpi=150, bbox_inches='tight')
        print(f'  Saved: {fname}')
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# PLOT 7: Capacity constraint comparison
# ─────────────────────────────────────────────────────────────────────────────

def plot_capacity_comparison(logs_by_thresh: dict, thresholds: list,
                              gov_rule: str, ag_rule: str,
                              regime: int, climate: float, save: bool):
    """
    Top panel  : CI trajectories overlaid (one line per threshold).
    Lower strips: Gov strategy executed per period for each threshold.
                  Grey = None (capacity check failed).
    """
    THRESH_COLORS = ['#27ae60', '#e67e22', '#c0392b']   # green / orange / red
    n_t = len(thresholds)

    fig, axes = plt.subplots(
        n_t + 1, 1, figsize=(12, 3.5 + 1.6 * n_t),
        gridspec_kw={'height_ratios': [3] + [1] * n_t},
        sharex=False)
    ax_ci = axes[0]

    _ci_bounds = [0.0] + list(_CI_THRESHOLDS) + [3.0]
    for thr in _CI_THRESHOLDS:
        ax_ci.axhline(thr, color='#cccccc', lw=0.5, ls=':', zorder=1)
    for zi in range(6):
        zmid = (_ci_bounds[zi] + _ci_bounds[zi + 1]) / 2
        ax_ci.annotate(f'R{zi+1}', xy=(1.005, zmid),
                       xycoords=('axes fraction', 'data'),
                       fontsize=6.5, va='center', ha='left', color='#777777')

    t_max = 0
    for ti, (thresh, color) in enumerate(zip(thresholds, THRESH_COLORS)):
        log = logs_by_thresh.get(thresh, [])
        if not log:
            continue
        t_all, ci_all = [], []
        for entry in log:
            t_local = entry['traj'].index.values + entry['t_offset']
            t_all.extend(t_local.tolist())
            ci_all.extend(_ci_from_traj(entry['traj']).tolist())
        t_max = max(t_max, t_all[-1])
        lbl = 'deterministic (thresh=0)' if thresh == 0 else f'thresh={thresh:.0f}'
        ax_ci.plot(t_all, ci_all, color=color, lw=1.6, label=lbl, zorder=5)
        for entry in log:
            if entry['period'] > 1:
                ax_ci.axvline(entry['t_offset'], color='#dddddd', lw=0.5, ls='--')

        # Strategy strip
        ax_s = axes[ti + 1]
        for entry in log:
            t0 = entry['t_offset']
            t1 = t0 + entry['T_period']
            ax_s.axvspan(t0, t1, color=STRATEGY_COLORS['Gov'][entry['gov_strat']],
                         alpha=0.9)
            if entry['period'] > 1:
                ax_s.axvline(t0, color='white', lw=0.5)
        ax_s.set_yticks([])
        lbl_short = 'det.' if thresh == 0 else f'thr={thresh:.0f}'
        ax_s.set_ylabel(lbl_short, fontsize=7, rotation=0, ha='right', va='center')
        ax_s.set_xlim(0, t_max)
        ax_s.tick_params(labelsize=6)
        if ti < n_t - 1:
            ax_s.set_xticklabels([])

    ax_ci.set_ylim(0, 3.05)
    ax_ci.set_xlim(0, t_max)
    ax_ci.set_ylabel('Conflict Index', fontsize=9)
    ax_ci.set_title(
        f'Capacity constraint experiment  |  R{regime}  climate={climate:.2f}\n'
        f'Gov rule: {gov_rule}   AG rule: {ag_rule}   '
        f'(Strategy strips below — grey = capacity failure)',
        fontsize=9)
    ax_ci.legend(loc='upper right', fontsize=8)
    ax_ci.tick_params(labelsize=7)
    ax_ci.grid(True, axis='y', lw=0.3, alpha=0.35)
    ax_ci.set_xticklabels([])
    axes[-1].set_xlabel('years', fontsize=8)

    plt.tight_layout()
    if save:
        fname = os.path.join(
            VIS_DIR,
            f'capacity_comparison_R{regime}_c{climate:.2f}_{gov_rule}_{ag_rule}.png')
        fig.savefig(fname, dpi=130, bbox_inches='tight')
        print(f'  Saved: {fname}')
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# PLOT 8: Re-evaluation frequency comparison
# ─────────────────────────────────────────────────────────────────────────────

def plot_freq_comparison(logs_by_mult: dict, mults: list,
                          gov_rule: str, ag_rule: str,
                          regime: int, climate: float, save: bool):
    """
    Top panel  : CI trajectories overlaid (one line per gov_period_mult).
    Lower strips: Gov strategy per period; black tick = re-evaluation point.
    """
    MULT_COLORS = ['#2980b9', '#8e44ad', '#c0392b']
    n_m = len(mults)

    fig, axes = plt.subplots(
        n_m + 1, 1, figsize=(12, 3.5 + 1.6 * n_m),
        gridspec_kw={'height_ratios': [3] + [1] * n_m},
        sharex=False)
    ax_ci = axes[0]

    _ci_bounds = [0.0] + list(_CI_THRESHOLDS) + [3.0]
    for thr in _CI_THRESHOLDS:
        ax_ci.axhline(thr, color='#cccccc', lw=0.5, ls=':', zorder=1)
    for zi in range(6):
        zmid = (_ci_bounds[zi] + _ci_bounds[zi + 1]) / 2
        ax_ci.annotate(f'R{zi+1}', xy=(1.005, zmid),
                       xycoords=('axes fraction', 'data'),
                       fontsize=6.5, va='center', ha='left', color='#777777')

    t_max = 0
    for mi, (mult, color) in enumerate(zip(mults, MULT_COLORS)):
        log = logs_by_mult.get(mult, [])
        if not log:
            continue
        t_all, ci_all = [], []
        for entry in log:
            t_local = entry['traj'].index.values + entry['t_offset']
            t_all.extend(t_local.tolist())
            ci_all.extend(_ci_from_traj(entry['traj']).tolist())
        t_max = max(t_max, t_all[-1])
        ax_ci.plot(t_all, ci_all, color=color, lw=1.6,
                   label=f'gov re-eval every {mult} period(s)', zorder=5)
        for entry in log:
            if entry['period'] > 1:
                ax_ci.axvline(entry['t_offset'], color='#dddddd', lw=0.5, ls='--')

        # Strategy strip
        ax_s = axes[mi + 1]
        for entry in log:
            t0 = entry['t_offset']
            t1 = t0 + entry['T_period']
            ax_s.axvspan(t0, t1, color=STRATEGY_COLORS['Gov'][entry['gov_strat']],
                         alpha=0.9)
            if entry['period'] > 1:
                ax_s.axvline(t0, color='white', lw=0.5)
            # Black tick at re-evaluation point
            if entry['period'] % mult == 0:
                ax_s.axvline(t1, color='black', lw=1.2, zorder=5)
        ax_s.set_yticks([])
        ax_s.set_ylabel(f'mult={mult}', fontsize=7, rotation=0, ha='right', va='center')
        ax_s.set_xlim(0, t_max)
        ax_s.tick_params(labelsize=6)
        if mi < n_m - 1:
            ax_s.set_xticklabels([])

    ax_ci.set_ylim(0, 3.05)
    ax_ci.set_xlim(0, t_max)
    ax_ci.set_ylabel('Conflict Index', fontsize=9)
    ax_ci.set_title(
        f'Re-evaluation frequency experiment  |  R{regime}  climate={climate:.2f}\n'
        f'Gov rule: {gov_rule}   AG rule: {ag_rule}   '
        f'(Strategy strips below — black tick = re-evaluation)',
        fontsize=9)
    ax_ci.legend(loc='upper right', fontsize=8)
    ax_ci.tick_params(labelsize=7)
    ax_ci.grid(True, axis='y', lw=0.3, alpha=0.35)
    ax_ci.set_xticklabels([])
    axes[-1].set_xlabel('years', fontsize=8)

    plt.tight_layout()
    if save:
        fname = os.path.join(
            VIS_DIR,
            f'freq_comparison_R{regime}_c{climate:.2f}_{gov_rule}_{ag_rule}.png')
        fig.savefig(fname, dpi=130, bbox_inches='tight')
        print(f'  Saved: {fname}')
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# CONSOLE OUTPUT
# ─────────────────────────────────────────────────────────────────────────────

def print_payoff_table(all_logs: dict, gov_strats: list, ag_strats: list):
    """Print 4×4 outcome tables for each metric."""
    col_w = 12
    sep = '  ' + '-' * (16 + col_w * len(ag_strats))
    for metric, (label, sign) in METRICS.items():
        direction = 'higher=Gov wins' if sign > 0 else 'lower=Gov wins'
        print(f'\n  {label}  ({direction})')
        print(sep)
        header = f'  {"Gov \\ AG":<14}' + ''.join(f'{SHORT[ag]:>{col_w}}' for ag in ag_strats)
        print(header)
        print(sep)
        for gs in gov_strats:
            row = f'  {SHORT[gs]:<14}'
            for ag in ag_strats:
                log = all_logs.get((gs, ag), [])
                val = (log[-1]['metrics'].get(metric, float('nan'))
                       if log else float('nan'))
                row += f'{val:>{col_w}.2f}'
            print(row)
        print(sep)


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description='APS Strategy Game  (fixed or adaptive modes)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
MODES
  fixed     Run all strategy combos (default 4x4=16) and compare payoff matrices.
  adaptive  One run per (regime x climate): each actor picks strategies via named
            rules evaluated at the end of every period (or every N periods).

ADAPTIVE RULES — Gov
  reactive          Dev<R3, Balance R3-R4, Defense>R4 based on CI
  development_first Always Development; Defense only in extreme conflict (R6)
  defense_first     Always Defense; Development only in near-peace (R1)
  adaptive_balance  Social cohesion low->Dev; AG threat->Def; else Balance

ADAPTIVE RULES — AG
  opportunistic  Exploit poverty->Coerce; war economy->Extract; else None
  territorial    Always Coerce
  economic       Always Extract
  pressure       Coerce in low conflict, Extract in high conflict

CAPACITY THRESHOLDS
  --gov-capacity T   p(Gov implements) = min(1, Gov_stock/T). T=0 -> always.
  --ag-capacity  T   p(AG  implements) = min(1, (AG_RL+AG_SL)/T). T=0 -> always.

PERIOD MULTIPLIERS (adaptive only)
  --gov-period-mult N   Gov re-evaluates strategy every N periods (default 1).
  --ag-period-mult  N   AG  re-evaluates strategy every N periods (default 1).

CLI EXAMPLES
  # Fixed mode (default): all 16 combos, regime R3, 15 periods
  python script_strategy_game.py --regime 3

  # Adaptive mode: reactive Gov vs opportunistic AG, starting R4
  python script_strategy_game.py --mode adaptive --regime 4 --gov-rule reactive --ag-rule opportunistic

  # Adaptive with stochastic capacity (Gov needs Gov>=5 to act):
  python script_strategy_game.py --mode adaptive --gov-capacity 5

  # Gov re-evaluates every 3 periods, AG every period:
  python script_strategy_game.py --mode adaptive --gov-period-mult 3
""")

    # ── Shared args ───────────────────────────────────────────────────────────
    p.add_argument('--regime', type=int, nargs='+', default=[1,2,3],
                   help='Starting regime(s)')
    p.add_argument('--climate', type=float, nargs='+', default=[0.1,0.4,0.7],
                   help='Climate value(s) in [0,1]')
    p.add_argument('--boost', type=float, default=5,
                   help='Strategy rate multiplier (default: 5)')
    p.add_argument('--periods', type=int, default=15,
                   help='Number of game periods (default: 15)')
    p.add_argument('--T-period', type=float, default=1.0, dest='T_period',
                   help='Years per period (default: 1.0)')
    p.add_argument('--runs', type=int, default=10,
                   help='Stochastic replicates per period (default: 10)')
    p.add_argument('--dt', type=float, default=0.05,
                   help='Time step in years (default: 0.05)')
    p.add_argument('--seed', type=int, default=42, help='Base random seed')
    p.add_argument('--save', action='store_true', default=True,
                   help='Save all plots to visualizations/strategy_game/')
    p.add_argument('--no-show', action='store_true',
                   help='Suppress plt.show() (headless mode)')

    # ── Fixed-mode args ───────────────────────────────────────────────────────
    p.add_argument('--gov', type=str, nargs='+', default=GOV_STRATEGIES,
                   choices=GOV_STRATEGIES, metavar='STRAT',
                   help='[fixed] Gov strategies to run (default: all 4)')
    p.add_argument('--ag', type=str, nargs='+', default=AG_STRATEGIES,
                   choices=AG_STRATEGIES, metavar='STRAT',
                   help='[fixed] AG strategies to run (default: all 4)')
    p.add_argument('--metric', type=str, default='conflict_index',
                   choices=list(METRICS.keys()),
                   help='[fixed] Metric shown in trajectory grid (default: conflict_index)')

    # ── Adaptive-mode args ────────────────────────────────────────────────────
    p.add_argument('--mode', choices=['fixed', 'adaptive'], default='fixed',
                   help='Execution mode (default: fixed)')
    p.add_argument('--gov-rule', default='reactive',
                   choices=list(GOV_RULES.keys()),
                   help='[adaptive] Gov decision rule (default: reactive)')
    p.add_argument('--ag-rule', default='opportunistic',
                   choices=list(AG_RULES.keys()),
                   help='[adaptive] AG  decision rule (default: opportunistic)')
    p.add_argument('--gov-init', default='Development',
                   choices=GOV_STRATEGIES,
                   help='[adaptive] Gov initial strategy (default: Development)')
    p.add_argument('--ag-init', default='Coerce',
                   choices=AG_STRATEGIES,
                   help='[adaptive] AG  initial strategy (default: Coerce)')
    p.add_argument('--gov-capacity', type=float, default=0.0, dest='gov_capacity',
                   help='[adaptive] Gov capacity threshold (0=deterministic)')
    p.add_argument('--ag-capacity',  type=float, default=0.0, dest='ag_capacity',
                   help='[adaptive] AG  capacity threshold (0=deterministic)')
    p.add_argument('--gov-period-mult', type=int, default=1, dest='gov_period_mult',
                   help='[adaptive] Gov re-evaluates every N periods (default: 1)')
    p.add_argument('--ag-period-mult',  type=int, default=1, dest='ag_period_mult',
                   help='[adaptive] AG  re-evaluates every N periods (default: 1)')

    # ── Adaptive batch-mode flags ─────────────────────────────────────────────
    p.add_argument('--all-rules', action='store_true', default=False,
                   dest='all_rules',
                   help='[adaptive] Run all 4x4 rule combinations; '
                        'generates rule matrix + per-pair plots + reports')
    p.add_argument('--capacity-sweep', action='store_true', default=False,
                   dest='capacity_sweep',
                   help='[adaptive] Sweep gov-capacity thresholds (see --cap-thresholds)')
    p.add_argument('--cap-thresholds', type=float, nargs='+', default=[0, 5, 15],
                   dest='cap_thresholds',
                   help='[adaptive] Gov capacity thresholds for sweep (default: 0 5 15)')
    p.add_argument('--freq-sweep', action='store_true', default=False,
                   dest='freq_sweep',
                   help='[adaptive] Sweep gov-period-mult values (see --freq-mults)')
    p.add_argument('--freq-mults', type=int, nargs='+', default=[1, 3, 6],
                   dest='freq_mults',
                   help='[adaptive] Gov period multipliers for sweep (default: 1 3 6)')

    return p.parse_args()


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    print(f'\n{"="*68}')
    print(f'  APS Strategy Game  ({args.mode.upper()} mode)')
    print(f'  Starting regimes : {args.regime}')
    print(f'  Climates         : {args.climate}')
    print(f'  Boost            : {args.boost}x  |  Periods: {args.periods}'
          f'  |  T/period: {args.T_period} yr')
    print(f'  Runs/period      : {args.runs}  |  dt: {args.dt}')
    if args.mode == 'adaptive':
        print(f'  Gov rule/init    : {args.gov_rule} / {args.gov_init}'
              f'  (every {args.gov_period_mult} period(s)'
              f', capacity threshold={args.gov_capacity})')
        print(f'  AG  rule/init    : {args.ag_rule} / {args.ag_init}'
              f'  (every {args.ag_period_mult} period(s)'
              f', capacity threshold={args.ag_capacity})')
    else:
        print(f'  Gov strategies   : {args.gov}')
        print(f'  AG  strategies   : {args.ag}')
    print(f'{"="*68}\n')

    # ── FIXED MODE ────────────────────────────────────────────────────────────
    if args.mode == 'fixed':
        gov_strats = args.gov
        ag_strats  = args.ag

        results_by_rc: Dict[Tuple[int, float], dict] = {}   # for combined CI plot

        for regime in args.regime:
            for climate in args.climate:
                print(f'  -- R{regime}  climate={climate:.2f} ' + '-' * 36)

                all_logs: Dict[Tuple[str, str], List[dict]] = {}
                n_combos = len(gov_strats) * len(ag_strats)
                done = 0
                for gs in gov_strats:
                    for ag in ag_strats:
                        done += 1
                        print(f'    [{done:2d}/{n_combos}]'
                              f'  Gov={gs:<12}  AG={ag:<8} ...',
                              end=' ', flush=True)
                        log = run_game(
                            start_regime=regime,
                            climate=climate,
                            kinetics=dict(DEFAULT_KINETICS),
                            gov_strat=gs,
                            ag_strat=ag,
                            boost=args.boost,
                            n_periods=args.periods,
                            T_period=args.T_period,
                            dt=args.dt,
                            n_runs=args.runs,
                            base_seed=args.seed,
                        )
                        all_logs[(gs, ag)] = log
                        r_end = log[-1]['metrics']['regime']
                        ep    = log[-1]['metrics']['E_peace']
                        ec    = log[-1]['metrics']['E_conflict']
                        print(f'-> R{r_end}  Ep={ep:.1f}  Ec={ec:.1f}')

                results_by_rc[(regime, climate)] = all_logs

                print_payoff_table(all_logs, gov_strats, ag_strats)
                print(f'\n  Plotting...')
                plot_strategy_trajectories(
                    all_logs, args.metric, climate, regime,
                    gov_strats, ag_strats, args.T_period, args.save)
                plot_payoff_matrix(
                    all_logs, climate, regime,
                    gov_strats, ag_strats, args.save)
                plot_reaction_activity(
                    all_logs, climate, regime,
                    gov_strats, ag_strats, args.save)
                plot_regime_evolution(
                    all_logs, climate, regime,
                    gov_strats, ag_strats, args.save)
                plt.close('all')
                print()

        # ── Combined CI matrix across all regimes × climates ──────────────────
        print('  Plotting combined CI matrix...')
        plot_combined_ci_matrix(
            results_by_rc, args.regime, args.climate,
            gov_strats, ag_strats, args.save)
        plt.close('all')

    # ── ADAPTIVE MODE ─────────────────────────────────────────────────────────
    else:
        def _run_adaptive(regime, climate, gov_rule, ag_rule,
                          gov_cap, ag_cap, gov_pm, ag_pm):
            """Helper: run one adaptive game and return log."""
            return run_game_adaptive(
                start_regime=regime, climate=climate,
                kinetics=dict(DEFAULT_KINETICS),
                gov_init=args.gov_init, ag_init=args.ag_init,
                boost=args.boost, n_periods=args.periods,
                T_period=args.T_period, dt=args.dt,
                n_runs=args.runs, base_seed=args.seed,
                gov_rule=gov_rule, ag_rule=ag_rule,
                gov_capacity_threshold=gov_cap,
                ag_capacity_threshold=ag_cap,
                gov_period_mult=gov_pm,
                ag_period_mult=ag_pm,
            )

        def _print_log_summary(log):
            print(f'  {"P":<4} {"Gov":<14} {"AG":<10} {"CI":>5} {"Regime":>8}')
            print('  ' + '-' * 48)
            for entry in log:
                ci = entry['indicators']['conflict_index']
                r  = entry['metrics']['regime']
                print(f'  P{entry["period"]:<3} {entry["gov_strat"]:<14}'
                      f' {entry["ag_strat"]:<10} {ci:>5.2f}  R{r}')

        # ── Batch A: all-rules grid ────────────────────────────────────────
        if args.all_rules:
            results_by_rule_rc = {}
            for regime in args.regime:
                for climate in args.climate:
                    print(f'\n  -- ALL-RULES GRID  R{regime}  climate={climate:.2f} '
                          + '-' * 20)
                    results_by_rule = {}
                    n_combos = len(GOV_RULES_LIST) * len(AG_RULES_LIST)
                    done = 0
                    for gr in GOV_RULES_LIST:
                        for ar in AG_RULES_LIST:
                            done += 1
                            print(f'    [{done:2d}/{n_combos}]'
                                  f'  Gov={gr:<20} AG={ar:<14} ...',
                                  end=' ', flush=True)
                            log = _run_adaptive(regime, climate, gr, ar,
                                                args.gov_capacity, args.ag_capacity,
                                                args.gov_period_mult, args.ag_period_mult)
                            results_by_rule[(gr, ar)] = log
                            ci_end = log[-1]['indicators']['conflict_index']
                            r_end  = log[-1]['metrics']['regime']
                            print(f'-> CI={ci_end:.2f}  R{r_end}')
                            # Detail report
                            save_adaptive_detail_report(
                                log, gr, ar, regime, climate,
                                args.gov_capacity, args.ag_capacity,
                                args.gov_period_mult, args.ag_period_mult)

                    results_by_rule_rc[(regime, climate)] = results_by_rule
                    # Grid summary report
                    df_grid = save_adaptive_grid_report(
                        results_by_rule, GOV_RULES_LIST, AG_RULES_LIST,
                        regime, climate, args.gov_capacity, args.ag_capacity,
                        args.gov_period_mult, args.ag_period_mult)
                    print(f'\n  Best rule pair: '
                          f'Gov={df_grid.iloc[0]["gov_rule"]}  '
                          f'AG={df_grid.iloc[0]["ag_rule"]}  '
                          f'CI={df_grid.iloc[0]["final_CI"]:.3f}')
                    print(f'  Worst rule pair: '
                          f'Gov={df_grid.iloc[-1]["gov_rule"]}  '
                          f'AG={df_grid.iloc[-1]["ag_rule"]}  '
                          f'CI={df_grid.iloc[-1]["final_CI"]:.3f}')
                    # Individual detailed plot only for best and worst pairs
                    print('  Plotting best/worst individual runs...')
                    for row_idx in [0, -1]:
                        gr_pair = df_grid.iloc[row_idx]['gov_rule']
                        ar_pair = df_grid.iloc[row_idx]['ag_rule']
                        plot_adaptive_game(results_by_rule[(gr_pair, ar_pair)],
                                           climate, regime,
                                           gr_pair, ar_pair, args.save)
                        plt.close('all')

                # Per-regime aggregate plots (across all climates)
                print(f'\n  Aggregate plots for R{regime}...')
                regime_rc = {k: v for k, v in results_by_rule_rc.items()
                             if k[0] == regime}
                plot_best_worst_trajectories(
                    regime_rc, [regime], args.climate,
                    GOV_RULES_LIST, AG_RULES_LIST, save=args.save)
                plot_strategy_usage(
                    regime_rc, [regime], args.climate,
                    GOV_RULES_LIST, AG_RULES_LIST, save=args.save)
                plt.close('all')

            # Combined rule matrix + marginal performance (all regimes × climates)
            print('\n  Plotting combined adaptive rule matrix...')
            plot_adaptive_rule_matrix(
                results_by_rule_rc, args.regime, args.climate,
                GOV_RULES_LIST, AG_RULES_LIST, args.save)
            print('  Plotting rule marginal performance...')
            plot_rule_marginal_performance(
                results_by_rule_rc, GOV_RULES_LIST, AG_RULES_LIST, save=args.save)
            plt.close('all')

        # ── Batch B: capacity sweep ────────────────────────────────────────
        elif args.capacity_sweep:
            thresholds = args.cap_thresholds
            for regime in args.regime:
                logs_all_climates = {}   # {(climate, thresh): log} for robustness plot
                for climate in args.climate:
                    print(f'\n  -- CAPACITY SWEEP  R{regime}  climate={climate:.2f}'
                          f'  Gov:{args.gov_rule}  AG:{args.ag_rule} ' + '-' * 20)
                    logs_by_thresh = {}
                    for thresh in thresholds:
                        print(f'    gov-capacity={thresh:.0f} ...', end=' ', flush=True)
                        log = _run_adaptive(regime, climate,
                                            args.gov_rule, args.ag_rule,
                                            thresh, args.ag_capacity,
                                            args.gov_period_mult, args.ag_period_mult)
                        logs_by_thresh[thresh] = log
                        logs_all_climates[(climate, thresh)] = log
                        ci_end = log[-1]['indicators']['conflict_index']
                        print(f'-> CI={ci_end:.2f}  R{log[-1]["metrics"]["regime"]}')
                        save_adaptive_detail_report(
                            log, args.gov_rule, args.ag_rule, regime, climate,
                            thresh, args.ag_capacity,
                            args.gov_period_mult, args.ag_period_mult)
                    save_capacity_sweep_report(logs_by_thresh, thresholds,
                                               args.gov_rule, args.ag_rule,
                                               regime, climate)
                    print('  Plotting capacity comparison...')
                    plot_capacity_comparison(logs_by_thresh, thresholds,
                                             args.gov_rule, args.ag_rule,
                                             regime, climate, args.save)
                    plt.close('all')
                    print()
                # Robustness plot across all climates for this regime
                print(f'  Plotting capacity robustness for R{regime}...')
                plot_capacity_robustness(
                    logs_all_climates, thresholds, args.climate,
                    args.gov_rule, args.ag_rule, regime, save=args.save)
                plt.close('all')

        # ── Batch C: frequency sweep ───────────────────────────────────────
        elif args.freq_sweep:
            mults = args.freq_mults
            for regime in args.regime:
                for climate in args.climate:
                    print(f'\n  -- FREQ SWEEP  R{regime}  climate={climate:.2f}'
                          f'  Gov:{args.gov_rule}  AG:{args.ag_rule} ' + '-' * 20)
                    logs_by_mult = {}
                    for mult in mults:
                        print(f'    gov-period-mult={mult} ...', end=' ', flush=True)
                        log = _run_adaptive(regime, climate,
                                            args.gov_rule, args.ag_rule,
                                            args.gov_capacity, args.ag_capacity,
                                            mult, args.ag_period_mult)
                        logs_by_mult[mult] = log
                        ci_end = log[-1]['indicators']['conflict_index']
                        print(f'-> CI={ci_end:.2f}  R{log[-1]["metrics"]["regime"]}')
                        save_adaptive_detail_report(
                            log, args.gov_rule, args.ag_rule, regime, climate,
                            args.gov_capacity, args.ag_capacity,
                            mult, args.ag_period_mult)
                    save_freq_sweep_report(logs_by_mult, mults,
                                           args.gov_rule, args.ag_rule,
                                           regime, climate)
                    print('  Plotting frequency comparison...')
                    plot_freq_comparison(logs_by_mult, mults,
                                         args.gov_rule, args.ag_rule,
                                         regime, climate, args.save)
                    plt.close('all')
                    print()

        # ── Single run (original behaviour) ───────────────────────────────
        else:
            for regime in args.regime:
                for climate in args.climate:
                    print(f'  -- R{regime}  climate={climate:.2f}'
                          f'  Gov:{args.gov_rule}  AG:{args.ag_rule} ' + '-' * 20)
                    log = _run_adaptive(regime, climate,
                                        args.gov_rule, args.ag_rule,
                                        args.gov_capacity, args.ag_capacity,
                                        args.gov_period_mult, args.ag_period_mult)
                    _print_log_summary(log)
                    save_adaptive_detail_report(
                        log, args.gov_rule, args.ag_rule, regime, climate,
                        args.gov_capacity, args.ag_capacity,
                        args.gov_period_mult, args.ag_period_mult)
                    print(f'\n  Plotting...')
                    plot_adaptive_game(log, climate, regime,
                                       args.gov_rule, args.ag_rule, args.save)
                    plt.close('all')
                    print()

    if not args.no_show:
        plt.show()

    print('Done.')


if __name__ == '__main__':
    main()
