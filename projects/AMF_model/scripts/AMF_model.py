#!/usr/bin/env python3
"""
Publication figures and tables for the AMF trade-balance reaction-network model
(network AMF_RN_MAK_v6.txt, 32 reactions, 12 species; Chemical Organisation Theory, pyCOT).

SELF-CONTAINED REPLICATION SCRIPT
---------------------------------
This single script reproduces every figure and table of the article from scratch. It needs
only the network file and the pyCOT library: the calibrated rate constants, the N and P supply
of each scenario, the initial conditions and the run settings are written below in fixed form,
so no other script, parameter file or previous output is read. The network is read with
pyCOT's read_txt, every trajectory is integrated with pyCOT's simulation (mass-action kinetics,
LSODA) and the organisations are computed with pyCOT's brute_force_organizations. All runs are
serial and every random draw comes from a fixed seed, so two executions give identical results.

Parameter sets (section 1)
--------------------------
PARAMS_BASE  base calibration (CMA-ES, Johnson 2009 constraints) with luxury P supply
             k20 = 100 in scenarios III and IV. Used in every figure.
PARAMS_ALT   alternative calibration forcing fungal N delivery >= 5 % in scenarios I and II.
             Used only in Fig. 6 (N-P trade-off).

Computations (run in this order by main())
------------------------------------------
1. Monte Carlo of the base calibration: n runs per scenario, every initial condition and rate
   constant perturbed uniformly by +/- pct %, each run paired with a non-mycorrhizal control
   (same perturbed values, symbiont species set to zero).
2. N x P supply phase diagram with the base parameters (horizon PHASE_TMULT x t_end).
3. Partial rank correlation coefficients (PRCC) of the rate constants on the Monte Carlo runs.
4. Monte Carlo of the alternative calibration (same seed, n, pct and horizon).
5. Supply sweeps of the base calibration (P supply k20 in III and IV; N supply k19 in II and IV),
   with the same perturbations at every level (paired comparison).

Figures (<outdir>/figuras, PDF vector + PNG)
--------------------------------------------
Fig. 1  Model structure: (a) plant-fungus-soil scheme of the reaction network; (b-d) the three
        nested organisations reached from different initial conditions with identical parameters.
Fig. 2  Monte Carlo dynamics of every species group in every scenario (base calibration).
Fig. 3  Reproduction of Johnson's (2010) trade-balance scenarios: (a) mycorrhizal growth response
        (MGR), (b) C supply, C demand and P benefit in the 2 x 2 N-P supply design, (c) MGR vs
        tissue N:P of the plant grown without the fungus (thresholds of Koerselman & Meuleman
        1996; the model N:P is in model units, not a mass ratio).
Fig. 4  N x P supply phase diagram: MGR in coexistence and regions of fungal exclusion and host
        collapse.
Fig. 5  Global sensitivity (PRCC) of the MGR and of fungal biomass.
Fig. 6  N-P trade-off: base calibration vs calibration forcing fungal N delivery.
Fig. 7  Sensitivity of the scenarios to P and N supply.
Fig. 8  Robustness: outcome of every Monte Carlo run and success vs commensalism tolerance.

Tables (<outdir>/tablas)
------------------------
1_Table_scenario_summary.xlsx   medians [IQR] and outcome frequencies per scenario
2_Table_success_tolerance.xlsx  success rate per scenario for several commensalism tolerances (Fig. 8)
3_Table_organizations.xlsx      organisations of the network (pyCOT) and destination of each IC
4_reproducibility.json          settings, parameters, file hashes and package versions
5_Table_MGR_harvest_times.xlsx  MGR and outcome class at experimental harvest times
6_Table_PRCC_full.xlsx          all PRCC values (Fig. 5 shows the most influential ones)
7_Table_MC_runs.xlsx            metrics and outcome of every Monte Carlo run (base and alternative)
                                and the phase-diagram grid
8_Table_supply_sweeps.xlsx      summary of the supply sweeps (Fig. 7)

Usage
-----
Place this file in pyCOT/projects/AMF_model/scripts (the network is read from
pyCOT/data/Ecological_models/AMF_RN_MAK_v6.txt) and run

    python figuras_articulo.py

The run settings (horizon, perturbation, number of runs, seed, sweeps, figures to skip, output
folder) are at the beginning of main(). Model-specific content (species groups, scheme layout,
reaction labels) is declared in the CONFIGURATION blocks and can be adapted to another pyCOT
network without touching the analysis or plotting code.
"""
import os
import sys
import json
import time
import hashlib
import platform
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.colors import TwoSlopeNorm, ListedColormap
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Patch
from matplotlib.lines import Line2D
from scipy import stats

# Add pyCOT to the path (this script lives in pyCOT/projects/AMF_model/scripts)
_script_dir = os.path.dirname(os.path.abspath(__file__))
_pycot_root = os.path.normpath(os.path.join(_script_dir, '..', '..', '..'))
sys.path.insert(0, os.path.join(_pycot_root, 'src'))

from pyCOT.io.functions import read_txt
from pyCOT.simulations.ode import simulation

NETWORK_FILE = os.path.join(_pycot_root, 'data', 'Ecological_models', 'AMF_RN_MAK_v6.txt')

# ======================================================================
# 1. PARAMETERS AND INITIAL CONDITIONS (fixed)
# ======================================================================
# Base calibration: CMA-ES in log10(k) against the trade-balance scenarios of Johnson (2009, 2010).
# Free constants: r9-r18, r26-r30; the others keep their literature-based values. One shared vector
# for the four scenarios; only the external N (r19) and P (r20) supplies change between scenarios.
# Calibrated luxury P supply was 73.75; the article uses k20 = 100 in III and IV (base case A100).
PARAMS_BASE = {
    # --- 1. C fixation and plant dynamics ---
    'r1': 0.8585840472, 'r2': 5.759302053e-05, 'r3': 0.003651665857, 'r4': 0.226801881,
    # --- 2. root uptake and plant recycling ---
    'r5': 0.3433837821, 'r6': 0.02972854989, 'r7': 2.764655722e-05, 'r8': 0.01631043586,
    # --- 3. hyphal uptake and mycelium dynamics ---
    'r9': 0.25597763317872, 'r10': 0.8177143367007289, 'r11': 42.327448817787676,
    'r12': 0.08583546417955697, 'r27': 0.30018070705448724,
    # --- 4. symbiotic exchange, regulation and intraradical fungal dynamics ---
    'r13': 0.015444527176148151, 'r14': 0.00011331018989759174, 'r15': 8.2974829290558,
    'r16': 28.050684097685146, 'r17': 0.003068367276602748, 'r18': 0.04347158619645481,
    # --- 5. external soil supply (overwritten per scenario, see SUPPLY_BASE) ---
    'r19': 1.0, 'r20': 1.0,
    # --- 6. self-limitation reactions 2X => X ---
    'r21': 0.1742528886, 'r22': 0.00284359432, 'r23': 0.1153093823, 'r24': 0.103839988,
    'r25': 10.61036832, 'r26': 1.6822986542229676, 'r28': 2.1980386831059064,
    'r29': 0.00454957491299379, 'r30': 77.65965886049912, 'r31': 0.0005002545511, 'r32': 0.08727425082,
}
# External N (r19) and P (r20) supply per scenario: 2 x 2 limitation/luxury design
SUPPLY_BASE = {
    'I':   {'r19': 0.7300480596, 'r20': 3.115516772},
    'II':  {'r19': 3.199996005,  'r20': 3.115516772},
    'III': {'r19': 0.7300480596, 'r20': 100.0},
    'IV':  {'r19': 3.199996005,  'r20': 100.0},
}

# Alternative calibration (Fig. 6): same objective plus fungal N delivery in [5, 40] % in I and II,
# started from the base calibration and with the same supply (k20 = 100 in III and IV).
PARAMS_ALT = {
    'r1': 0.8585840472, 'r2': 5.759302053e-05, 'r3': 0.003651665857, 'r4': 0.226801881,
    'r5': 0.3433837821, 'r6': 0.02972854989, 'r7': 2.764655722e-05, 'r8': 0.01631043586,
    'r9': 1.1493407485554916, 'r10': 0.8868310514285082, 'r11': 28.631711846791077,
    'r12': 0.0502535583213668, 'r27': 2.9371420252168083,
    'r13': 0.03414834387017321, 'r14': 7.598496209442223e-05, 'r15': 2.904417187494706,
    'r16': 23.25677345222023, 'r17': 0.0050773869739926735, 'r18': 0.042210653351482094,
    'r19': 1.0, 'r20': 1.0,
    'r21': 0.1742528886, 'r22': 0.00284359432, 'r23': 0.1153093823, 'r24': 0.103839988,
    'r25': 10.61036832, 'r26': 1.164762262156337, 'r28': 10.029287740614276,
    'r29': 0.0034389909632820577, 'r30': 59.434668923716, 'r31': 0.0005002545511, 'r32': 0.08727425082,
}
SUPPLY_ALT = {s: dict(v) for s, v in SUPPLY_BASE.items()}

CALIBRATIONS = {
    'base': {'label': 'base calibration', 'params': PARAMS_BASE, 'supply': SUPPLY_BASE},
    'alt':  {'label': 'forced N trade (≥ 5 %)', 'params': PARAMS_ALT, 'supply': SUPPLY_ALT},
}

# Seedling initial condition with a fungal inoculum above the invasion threshold. Soil N and P
# start at SOIL_IC_FACTOR times their abiotic equilibrium sqrt(k_supply / k_self-limitation).
x0_dict_base = {'Pact': 0.4, 'Plim': 0.1, 'Cp': 0.5, 'Np': 0.1, 'Pp': 0.05,
                'F': 0.1, 'M': 0.1, 'Cf': 0.1, 'Nf': 0.1, 'Pf': 0.1}
SOIL_IC_FACTOR = 0.1

# Active parameter set (filled by set_calibration; every analysis function reads these two)
params_dict_base = {}
SCENARIO_SUPPLY = {}


def set_calibration(name):
    """Make one of CALIBRATIONS the active parameter set."""
    cal = CALIBRATIONS[name]
    params_dict_base.clear()
    params_dict_base.update({r: [float(v)] for r, v in cal['params'].items()})
    SCENARIO_SUPPLY.clear()
    SCENARIO_SUPPLY.update({s: {r: [float(v)] for r, v in d.items()} for s, d in cal['supply'].items()})


set_calibration('base')

# ======================================================================
# 2. MODEL CONFIGURATION (species groups, criteria; adapt if the network changes)
# ======================================================================
PLANT_BIOMASS   = ['Pact', 'Plim']                    # total plant biomass (MGR)
PLANT_VIABLE    = ['Pact']                            # biomass able to take up nutrients and grow (classifier)
PLANT_SPECIES   = ['Pact', 'Plim', 'Cp', 'Np', 'Pp']  # plant-soil organisation (without soil)
SYMBIONT        = ['F', 'M', 'Cf', 'Nf', 'Pf']        # set to zero in the non-mycorrhizal control
FUNGAL_BIOMASS  = ['F', 'M']
PLANT_N, PLANT_P = 'Np', 'Pp'
FLUX_GROUPS = {
    'C_fix':       ['r1'],            # C fixed by the plant
    'C_to_fungus': ['r13', 'r15'],    # plant C transferred to the fungus
    'P_direct':    ['r6'],            # P taken up by the root
    'P_fungal':    ['r15'],           # P delivered by the fungus
    'N_direct':    ['r5'],            # N taken up by the root
    'N_fungal':    ['r13', 'r14'],    # N delivered by the fungus
    'N_fungal_uptake': ['r9'],        # N taken up from the soil by the mycelium
    'P_fungal_uptake': ['r10'],       # P taken up from the soil by the mycelium
}
FUNGAL_N         = 'Nf'
N_BENEFIT_LIT    = (0.05, 0.40)   # fraction of plant N supplied by the fungus (reference band)
MGR_COMM_TOL     = 0.05           # |MGR| <= tol: commensalism
P_BENEFIT_HIGH   = 0.50           # mycorrhizal pathway dominant for plant P
SUCCESS_CRITERIA = {
    'I':   [('MGR > tol (mutualism)',            lambda r, m: r['MGR'] > MGR_COMM_TOL),
            ('P_benefit >= 0.5',                 lambda r, m: r['P_fungal_frac'] >= P_BENEFIT_HIGH),
            ('MGR < MGR(II): C-limited',         lambda r, m: r['MGR'] < m['II']['MGR'])],
    'II':  [('MGR > tol (mutualism)',            lambda r, m: r['MGR'] > MGR_COMM_TOL),
            ('P_benefit >= 0.5',                 lambda r, m: r['P_fungal_frac'] >= P_BENEFIT_HIGH),
            ('MGR > MGR(I): strong mutualism',   lambda r, m: r['MGR'] > m['I']['MGR'])],
    'III': [('|MGR| <= tol (commensalism)',      lambda r, m: abs(r['MGR']) <= MGR_COMM_TOL),
            ('P_benefit < 0.5',                  lambda r, m: r['P_fungal_frac'] < P_BENEFIT_HIGH),
            ('C sink < C sink(II) and (IV)',     lambda r, m: r['C_to_fungus'] < min(m['II']['C_to_fungus'],
                                                                                     m['IV']['C_to_fungus']))],
    'IV':  [('MGR < -tol (parasitism)',          lambda r, m: r['MGR'] < -MGR_COMM_TOL),
            ('P_benefit < 0.5',                  lambda r, m: r['P_fungal_frac'] < P_BENEFIT_HIGH),
            ('C sink > C sink(II)',              lambda r, m: r['C_to_fungus'] > m['II']['C_to_fungus'])],
}
C_DEMAND_LIT     = (0.04, 0.20)   # fraction of photosynthate allocated to AMF (literature)
ORG_DEMO_SCENARIO = 'II'          # scenario of the organisation demonstration (Fig. 1)
INOCULUM_SUBTHRESHOLD = 1e-3      # inoculum below the invasion threshold (Fig. 1)
ZERO_TOL         = 1e-6           # absent species
PRESENCE_TOL     = 0.05           # minimum biomass to count as present
DECAY_TOL        = 0.20           # relative drop in the final window that indicates ongoing extinction
CONV_TOL         = 0.02           # maximum relative change in the final window for "converged"
STEADY_FRAC      = 0.2            # final fraction of the horizon used for steady-state averages
PHASE_N_GRID     = 21             # points per axis of the phase diagram (Fig. 4)
PHASE_SPAN       = 1.0            # log10 extension beyond the scenario supply levels
PHASE_TMULT      = 4              # phase-diagram horizon = PHASE_TMULT x t_end
SUPPLY_REACTIONS = {'N': 'r19', 'P': 'r20'}
SOIL_CONTROL     = {'N': 'r31', 'P': 'r32'}
ORG_CLASSIFY_T   = 20000          # extra days used to classify the destination (Fig. 1)
PRCC_OUTPUTS     = {'MGR': 'O3', 'Fungal_biomass': 'all', 'C_demand': 'O3'}  # output: subset of runs
PRCC_MIN_EXTRA   = 10             # PRCC computed only if n >= n_parameters + PRCC_MIN_EXTRA
SIM_TOL          = {'rtol': 1e-6, 'atol': 1e-8}    # Monte Carlo and sweeps
SIM_TOL_PHASE    = {'rtol': 1e-7, 'atol': 1e-10}   # phase diagram

SCENARIOS = ['I', 'II', 'III', 'IV']
STATE_ORDER = ['éxito', 'fracaso', 'exclusión del hongo', 'colapso del hospedero']   # internal keys
ORG_TO_STATE = {'O2': 'exclusión del hongo', 'O1': 'colapso del hospedero'}


# ======================================================================
# 3. SIMULATION AND CLASSIFICATION (pyCOT)
# ======================================================================
def network_orders(rn):
    species_order = [s.name for s in rn.species()]
    reaction_order = [r.name() if callable(r.name) else r.name for r in rn.reactions()]
    return species_order, reaction_order


def simulate(rn, x0, p, t_span, n_steps, rtol=1e-6, atol=1e-8, allow_fail=False):
    """One pyCOT mass-action simulation; returns (time series, fluxes) or (None, None) if it fails
    and allow_fail is set."""
    try:
        return simulation(rn, x0=x0, rate=['mak'] * len(p), spec_vector=p, t_span=t_span,
                          n_steps=n_steps, method='LSODA', rtol=rtol, atol=atol, verbose=False)
    except Exception as e:
        if allow_fail:
            return None, None
        raise RuntimeError(f'Simulation failed: {e}') from e


class Progress:
    """Prints the progress of a serial loop every 10 %."""
    def __init__(self, total, desc):
        self.total, self.desc, self.done, self.t0 = total, desc, 0, time.perf_counter()
        self.step = max(1, total // 10)

    def tick(self):
        self.done += 1
        if self.done % self.step == 0 or self.done == self.total:
            el = time.perf_counter() - self.t0
            print(f'   [{self.desc}] {self.done}/{self.total} ({100 * self.done / self.total:.0f} %)  '
                  f'{el:.0f} s  ETA {el / self.done * (self.total - self.done):.0f} s', flush=True)


def parse_network_stoich(path):
    """{reaction: (reactants, products)} read from the pyCOT .txt network file."""
    def side(s):
        out = {}
        for t in s.split('+'):
            t = t.strip()
            if not t:
                continue
            k = 0
            while k < len(t) and (t[k].isdigit() or t[k] == '.'):
                k += 1
            c = float(t[:k]) if k else 1.0
            sp = t[k:].strip()
            out[sp] = out.get(sp, 0.0) + c
        return out
    st = {}
    for line in open(path, encoding='utf-8'):
        line = line.split('#')[0].strip().rstrip(';').strip()
        if not line or ':' not in line or '=>' not in line:
            continue
        name, eq = line.split(':', 1)
        lhs, rhs = eq.split('=>')
        st[name.strip()] = (side(lhs), side(rhs))
    return st


def get_scenario_config(scenario, species_order, reaction_order):
    """Shared rate constants + N/P supply of the scenario; seedling initial condition."""
    local_p = {r: list(v) for r, v in params_dict_base.items()}
    local_p.update({r: list(v) for r, v in SCENARIO_SUPPLY[scenario].items()})
    local_x0 = dict(x0_dict_base)
    local_x0['Ns'] = SOIL_IC_FACTOR * (local_p[SUPPLY_REACTIONS['N']][0] / local_p[SOIL_CONTROL['N']][0]) ** 0.5
    local_x0['Ps'] = SOIL_IC_FACTOR * (local_p[SUPPLY_REACTIONS['P']][0] / local_p[SOIL_CONTROL['P']][0]) ** 0.5
    x0_list = [local_x0.get(s, 0.0) for s in species_order]
    p_list = [list(local_p[r]) for r in reaction_order]
    return x0_list, p_list, ['mak'] * len(reaction_order)


def get_perturbed_value(val, pct_variation):
    """Uniform perturbation of +/- pct_variation %."""
    f = pct_variation / 100.0
    return max(0.0, np.random.uniform(val * (1 - f), val * (1 + f)))


def _window_index(t, t_span, steady_frac):
    return np.searchsorted(t, t_span[1] - steady_frac * (t_span[1] - t_span[0]))


def _ss_mean(df, col, t_span, steady_frac=STEADY_FRAC):
    """Mean of a column over the final window."""
    if col not in df.columns:
        return np.nan
    n = _window_index(df['Time'].values, t_span, steady_frac)
    return float(np.mean(df[col].values[n:]))


def _ss_flux(flux_df, group, t_span, steady_frac=STEADY_FRAC):
    return float(np.nansum([_ss_mean(flux_df, r, t_span, steady_frac) for r in FLUX_GROUPS[group]]))


def make_nm_x0(x0_list, species_order):
    """Initial condition of the non-mycorrhizal control: symbiont set to zero."""
    x0_nm = list(x0_list)
    for s in SYMBIONT:
        if s in species_order:
            x0_nm[species_order.index(s)] = 0.0
    return x0_nm


def _present(ts, species, t_span, steady_frac=STEADY_FRAC):
    """Present if the sum exceeds PRESENCE_TOL at the end and does not drop by more than
    DECAY_TOL over the final window (slow extinction)."""
    n = _window_index(ts['Time'].values, t_span, steady_frac)
    y = sum(ts[q].values for q in species if q in ts.columns)
    if y[-1] < PRESENCE_TOL:
        return False
    return (y[n] - y[-1]) / (y[n] + 1e-12) < DECAY_TOL


def classify_run(ts, t_span, steady_frac=STEADY_FRAC):
    """'O3' viable plant and fungus; 'O2' viable plant only; 'O1' no viable plant. Plant viability
    is evaluated on PLANT_VIABLE (Pact): residual Plim cannot be reactivated without Pact."""
    plant = _present(ts, PLANT_VIABLE, t_span, steady_frac)
    fungus = _present(ts, FUNGAL_BIOMASS, t_span, steady_frac)
    return 'O3' if (plant and fungus) else ('O2' if plant else 'O1')


def _converged(ts, species, t_span, steady_frac=STEADY_FRAC):
    n = _window_index(ts['Time'].values, t_span, steady_frac)
    y = sum(ts[q].values for q in species if q in ts.columns)
    if y[-1] < PRESENCE_TOL:
        return True
    return abs(y[n] - y[-1]) / (abs(y[-1]) + 1e-12) < CONV_TOL


def compute_johnson_metrics(all_ts_data, all_flux_data, all_ts_nm, t_span, steady_frac=STEADY_FRAC):
    """One row per (scenario, run): MGR, fluxes, organisation and convergence."""
    eps = 1e-12
    rows = []
    for scen in SCENARIOS:
        for i, (ts, fl, ts_nm) in enumerate(zip(all_ts_data[scen], all_flux_data[scen], all_ts_nm[scen])):
            B_am = sum(_ss_mean(ts, s, t_span, steady_frac) for s in PLANT_BIOMASS)
            B_nm = sum(_ss_mean(ts_nm, s, t_span, steady_frac) for s in PLANT_BIOMASS)
            c_fix, c_f = _ss_flux(fl, 'C_fix', t_span, steady_frac), _ss_flux(fl, 'C_to_fungus', t_span, steady_frac)
            p_d, p_f = _ss_flux(fl, 'P_direct', t_span, steady_frac), _ss_flux(fl, 'P_fungal', t_span, steady_frac)
            n_d, n_f = _ss_flux(fl, 'N_direct', t_span, steady_frac), _ss_flux(fl, 'N_fungal', t_span, steady_frac)
            org_am, org_nm = classify_run(ts, t_span, steady_frac), classify_run(ts_nm, t_span, steady_frac)
            c_dem = c_f / (c_fix + eps)
            rows.append({
                'Scenario': scen, 'run': i + 1,
                'B_AM': B_am, 'B_NM': B_nm,
                'MGR': float(np.log((B_am + eps) / (B_nm + eps))),
                'C_fix': c_fix,
                'C_to_fungus': c_f,                      # absolute fungal C sink
                'C_demand': c_dem,                       # fraction of fixed C
                'C_demand_en_rango': bool(C_DEMAND_LIT[0] <= c_dem <= C_DEMAND_LIT[1]) if org_am == 'O3' else np.nan,
                'P_fungal_frac': p_f / (p_d + p_f + eps),
                'N_fungal_frac': n_f / (n_d + n_f + eps),
                'N_direct': n_d, 'N_fungal': n_f,
                'N_fungal_uptake': _ss_flux(fl, 'N_fungal_uptake', t_span, steady_frac),
                'NP_plant_NM': _ss_mean(ts_nm, PLANT_N, t_span, steady_frac)
                               / (_ss_mean(ts_nm, PLANT_P, t_span, steady_frac) + eps),
                'Fungal_biomass': sum(_ss_mean(ts, s, t_span, steady_frac) for s in FUNGAL_BIOMASS),
                'Org_AM': org_am, 'Org_NM': org_nm,
                'colapso_hospedero': bool(org_am == 'O1' and org_nm != 'O1'),
                'convergido': bool(_converged(ts, PLANT_BIOMASS, t_span, steady_frac)
                                   and _converged(ts, FUNGAL_BIOMASS, t_span, steady_frac)
                                   and _converged(ts_nm, PLANT_BIOMASS, t_span, steady_frac)),
            })
    return pd.DataFrame(rows)


def evaluate_success(metrics):
    """Outcome of every run: success / failure (O3) / fungal exclusion (O2) / host collapse (O1).
    Between-scenario criteria use the medians of the O3 runs of each scenario."""
    o3 = metrics[metrics['Org_AM'] == 'O3']
    ref = o3.groupby('Scenario')[['MGR', 'C_demand', 'C_to_fungus', 'P_fungal_frac']].median().to_dict('index')
    for sc in SCENARIOS:
        ref.setdefault(sc, {'MGR': np.nan, 'C_demand': np.nan, 'C_to_fungus': np.nan, 'P_fungal_frac': np.nan})
    rows = []
    for _, r in metrics.iterrows():
        row = {'Scenario': r['Scenario'], 'run': r['run'], 'Org_AM': r['Org_AM']}
        crit = SUCCESS_CRITERIA[r['Scenario']]
        if r['Org_AM'] != 'O3':
            row['estado'] = ORG_TO_STATE[r['Org_AM']]
            for name, _f in crit:
                row[name] = np.nan
        else:
            oks = []
            for name, f in crit:
                v = bool(f(r, ref)); row[name] = v; oks.append(v)
            row['estado'] = 'éxito' if all(oks) else 'fracaso'
        rows.append(row)
    return pd.DataFrame(rows)


def compute_organizations(rn):
    """Organisations of the network with pyCOT (list of species sets, smallest first) or None."""
    try:
        from pyCOT.analysis.Persistent_Modules_Generator import brute_force_organizations
        res = brute_force_organizations(rn, verbose=False)
        orgs = []
        for o in res.get('organizations', []):
            sp = o.get('species', o) if isinstance(o, dict) else getattr(o, 'species', o)
            orgs.append(frozenset(x if isinstance(x, str) else x.name for x in sp))
        return sorted(orgs, key=len)
    except Exception as e:
        print(f'  Warning: pyCOT organisations could not be computed ({e})')
        return None


def _classify_support(final_state, orgs, decaying=frozenset()):
    supp = frozenset(s for s, v in final_state.items() if v > ZERO_TOL and s not in decaying)
    if not orgs:
        return None, supp
    for i, O in enumerate(orgs):
        if supp == O:
            return i, supp
    return int(np.argmin([len(supp ^ O) for O in orgs])), supp


# ======================================================================
# 4. ANALYSES (Monte Carlo, phase diagram, PRCC, supply sweeps)
# ======================================================================
def run_monte_carlo(rn, n_sims, pct, seed, t_end, desc='Monte Carlo'):
    """Paired AM / NM Monte Carlo of the active calibration. Every initial condition and rate
    constant is perturbed uniformly by +/- pct % (global seed; scenarios in order; for each run
    first the initial condition, then the rate constants)."""
    species_order, reaction_order = network_orders(rn)
    t_span, n_steps = (0, t_end), int(t_end) + 1
    np.random.seed(seed)
    cond = {sc: {'x0': [], 'p': []} for sc in SCENARIOS}
    for scen in SCENARIOS:
        x0_base, p_base, _ = get_scenario_config(scen, species_order, reaction_order)
        for _k in range(n_sims):
            cond[scen]['x0'].append([get_perturbed_value(v, pct) for v in x0_base])
            cond[scen]['p'].append([[get_perturbed_value(v, pct) for v in p] for p in p_base])
    ts_am = {sc: [] for sc in SCENARIOS}; fl_am = {sc: [] for sc in SCENARIOS}
    ts_nm = {sc: [] for sc in SCENARIOS}
    prog = Progress(2 * n_sims * len(SCENARIOS), desc)
    for scen in SCENARIOS:
        for x0, p in zip(cond[scen]['x0'], cond[scen]['p']):
            ts, fl = simulate(rn, x0, p, t_span, n_steps, **SIM_TOL); prog.tick()
            tn, _ = simulate(rn, make_nm_x0(x0, species_order), p, t_span, n_steps, **SIM_TOL); prog.tick()
            ts_am[scen].append(ts); fl_am[scen].append(fl); ts_nm[scen].append(tn)
    metrics = compute_johnson_metrics(ts_am, fl_am, ts_nm, t_span)
    success = evaluate_success(metrics)
    metrics = metrics.merge(success[['Scenario', 'run', 'estado']], on=['Scenario', 'run'], how='left')
    return {'metrics': metrics, 'success': success, 'cond': cond, 'ts_am': ts_am, 'ts_nm': ts_nm,
            'fl_am': fl_am, 't_end': t_end, 'phase': None, 'prcc': None}


def compute_phase_diagram(rn, t_end, n_grid=None):
    """MGR and final organisation on a log10(k_N) x log10(k_P) grid, unperturbed active parameters."""
    species_order, reaction_order = network_orders(rn)
    n_grid = n_grid or PHASE_N_GRID
    kN = sorted({v[SUPPLY_REACTIONS['N']][0] for v in SCENARIO_SUPPLY.values()})
    kP = sorted({v[SUPPLY_REACTIONS['P']][0] for v in SCENARIO_SUPPLY.values()})
    gN = np.logspace(np.log10(kN[0]) - PHASE_SPAN, np.log10(kN[-1]) + PHASE_SPAN, n_grid)
    gP = np.logspace(np.log10(kP[0]) - PHASE_SPAN, np.log10(kP[-1]) + PHASE_SPAN, n_grid)
    T = PHASE_TMULT * t_end
    ridx = {r: i for i, r in enumerate(reaction_order)}
    rows = []
    prog = Progress(2 * n_grid * n_grid, 'Phase diagram')
    for i, kp in enumerate(gP):
        for j, kn in enumerate(gN):
            p = [list(params_dict_base[r]) for r in reaction_order]
            p[ridx[SUPPLY_REACTIONS['N']]] = [kn]; p[ridx[SUPPLY_REACTIONS['P']]] = [kp]
            x0 = dict(x0_dict_base)
            x0['Ns'] = SOIL_IC_FACTOR * (kn / p[ridx[SOIL_CONTROL['N']]][0]) ** 0.5
            x0['Ps'] = SOIL_IC_FACTOR * (kp / p[ridx[SOIL_CONTROL['P']]][0]) ** 0.5
            x0l = [x0.get(sp, 0.0) for sp in species_order]
            ts, _ = simulate(rn, x0l, p, (0, T), 401, allow_fail=True, **SIM_TOL_PHASE); prog.tick()
            tn, _ = simulate(rn, make_nm_x0(x0l, species_order), p, (0, T), 401, allow_fail=True,
                             **SIM_TOL_PHASE); prog.tick()
            row = {'k_N': kn, 'k_P': kp, 'MGR': np.nan, 'Org_AM': None, 'Org_NM': None}
            if ts is not None and tn is not None:
                b = float(sum(ts[q].iloc[-1] for q in PLANT_BIOMASS))
                bn = float(sum(tn[q].iloc[-1] for q in PLANT_BIOMASS))
                row.update(MGR=np.log((b + 1e-9) / (bn + 1e-9)), Org_AM=classify_run(ts, (0, T)),
                           Org_NM=classify_run(tn, (0, T)))
            rows.append(row)
    return pd.DataFrame(rows)


def prcc(X, y):
    """Partial rank correlation coefficients (Marino et al. 2008) and p-values."""
    n, p = X.shape
    R = np.column_stack([stats.rankdata(X[:, j]) for j in range(p)])
    ry = stats.rankdata(y)
    coef, pval = np.full(p, np.nan), np.full(p, np.nan)
    for j in range(p):
        Z = np.column_stack([np.ones(n), np.delete(R, j, axis=1)])
        rx = R[:, j] - Z @ np.linalg.lstsq(Z, R[:, j], rcond=None)[0]
        rr = ry - Z @ np.linalg.lstsq(Z, ry, rcond=None)[0]
        if rx.std() == 0 or rr.std() == 0:
            continue
        c = float(np.corrcoef(rx, rr)[0, 1])
        df = n - 2 - (p - 1)
        if df > 0 and abs(c) < 1:
            t = c * np.sqrt(df / (1 - c ** 2))
            pval[j] = 2 * stats.t.sf(abs(t), df)
        coef[j] = c
    return coef, pval


def compute_prcc(metrics, cond, reaction_order, outputs=None):
    """PRCC of every perturbed rate constant on the outputs of PRCC_OUTPUTS, per scenario."""
    outputs = outputs or PRCC_OUTPUTS
    rows = []
    for scen in SCENARIOS:
        P = np.array([[v[0] for v in p] for p in cond[scen]['p']])
        keep = [j for j in range(P.shape[1]) if np.std(P[:, j]) > 0]
        names = [reaction_order[j] for j in keep]
        ms = metrics[metrics['Scenario'] == scen].reset_index(drop=True)
        for out, subset in outputs.items():
            mask = (ms['Org_AM'] == 'O3').values if subset == 'O3' else np.ones(len(ms), bool)
            y = ms.loc[mask, out].values
            X = P[mask][:, keep]
            if len(y) < len(keep) + PRCC_MIN_EXTRA:
                print(f'  Warning PRCC: scenario {scen}, output {out}: n = {len(y)} < {len(keep) + PRCC_MIN_EXTRA} '
                      f'(parameters + {PRCC_MIN_EXTRA}); skipped. Increase n_sims.')
                continue
            c, pv = prcc(X, y)
            for nm, ci, pi in zip(names, c, pv):
                rows.append({'Escenario': scen, 'Salida': out, 'Parámetro': nm, 'PRCC': ci, 'p': pi, 'n': len(y)})
    return pd.DataFrame(rows)


def run_supply_sweep(rn, reaction, scenarios, levels, mode, n_sims, pct, seed, t_end):
    """Monte Carlo at several values of one rate constant applied in `scenarios` (mode 'values':
    absolute values; 'factors': factors of each scenario's base value). The same perturbations are
    used at every level (paired comparison); the soil initial condition follows the supply, as in the
    Monte Carlo. Scenarios not swept are simulated once and repeated at every level."""
    species_order, reaction_order = network_orders(rn)
    R = reaction
    base = {sc: (SCENARIO_SUPPLY[sc][R][0] if R in SCENARIO_SUPPLY[sc] else params_dict_base[R][0])
            for sc in SCENARIOS}
    val_of = (lambda sc, lv: base[sc] * lv) if mode == 'factors' else (lambda sc, lv: lv)
    T, n_steps, f = (0, t_end), int(t_end) + 1, pct / 100
    np.random.seed(seed)
    draws = {}
    for sc in SCENARIOS:
        x0b, pb, _ = get_scenario_config(sc, species_order, reaction_order)
        draws[sc] = [(np.random.uniform(1 - f, 1 + f, len(x0b)), np.random.uniform(1 - f, 1 + f, len(pb)))
                     for _ in range(n_sims)]

    def run_level(sc, value):
        saved = SCENARIO_SUPPLY[sc].get(R)
        SCENARIO_SUPPLY[sc][R] = [value]
        x0b, pb, _ = get_scenario_config(sc, species_order, reaction_order)
        if saved is None:
            SCENARIO_SUPPLY[sc].pop(R, None)
        else:
            SCENARIO_SUPPLY[sc][R] = saved
        out = []
        for fx, fp in draws[sc]:
            x0 = list(np.asarray(x0b) * fx); p = [[v * fp[j] for v in pv] for j, pv in enumerate(pb)]
            a = simulate(rn, x0, p, T, n_steps, allow_fail=True, **SIM_TOL); prog.tick()
            b = simulate(rn, make_nm_x0(x0, species_order), p, T, n_steps, allow_fail=True, **SIM_TOL); prog.tick()
            out.append((a, b))
        return out

    n_jobs = 2 * n_sims * (len(scenarios) * len(levels) + len(SCENARIOS) - len(scenarios))
    prog = Progress(n_jobs, f'Sweep {R}')
    fixed = {sc: run_level(sc, base[sc]) for sc in SCENARIOS if sc not in scenarios}
    rows = []
    for lv in levels:
        ts_am = {sc: [] for sc in SCENARIOS}; fl_am = {sc: [] for sc in SCENARIOS}
        ts_nm = {sc: [] for sc in SCENARIOS}
        for sc in SCENARIOS:
            for a, b in (run_level(sc, val_of(sc, lv)) if sc in scenarios else fixed[sc]):
                if a[0] is None or b[0] is None:
                    continue
                ts_am[sc].append(a[0]); fl_am[sc].append(a[1]); ts_nm[sc].append(b[0])
        met = compute_johnson_metrics(ts_am, fl_am, ts_nm, T)
        suc = evaluate_success(met)
        st = suc.groupby('Scenario')['estado'].value_counts(normalize=True).unstack(fill_value=0) \
                .reindex(index=SCENARIOS, columns=STATE_ORDER, fill_value=0)
        o3 = met[met['Org_AM'] == 'O3']
        med = o3.groupby('Scenario')[['MGR', 'C_demand', 'C_to_fungus', 'P_fungal_frac', 'N_fungal_frac']]
        q, Q, md = (med.quantile(0.25).reindex(SCENARIOS), med.quantile(0.75).reindex(SCENARIOS),
                    med.median().reindex(SCENARIOS))
        for sc in SCENARIOS:
            rows.append({'nivel': lv, 'Escenario': sc, f'valor_{R}': val_of(sc, lv) if sc in scenarios else base[sc],
                         '% éxito': 100 * st.loc[sc, 'éxito'],
                         '% exclusión hongo': 100 * st.loc[sc, 'exclusión del hongo'],
                         '% colapso hospedero': 100 * st.loc[sc, 'colapso del hospedero'],
                         'MGR_med': md.loc[sc, 'MGR'], 'MGR_q25': q.loc[sc, 'MGR'], 'MGR_q75': Q.loc[sc, 'MGR'],
                         'C_demand_med': md.loc[sc, 'C_demand'], 'C_to_fungus_med': md.loc[sc, 'C_to_fungus'],
                         'P_benefit_med': md.loc[sc, 'P_fungal_frac'], 'N_benefit_med': md.loc[sc, 'N_fungal_frac'],
                         'n_runs': len(met[met['Scenario'] == sc])})
    return {'summary': pd.DataFrame(rows), 'reaction': R, 'mode': mode}


def compute_species_fate(stoich, metrics, all_flux_data, species, t_span, steady_frac=STEADY_FRAC):
    """Steady-state fate of a species (median of the O3 runs): for every reaction with net consumption
    of the species, consumed flux = rate x net coefficient, and its share of the total consumption."""
    net = {r: re_.get(species, 0) - pr.get(species, 0) for r, (re_, pr) in stoich.items()}
    cons = [r for r, c in net.items() if c > 0]
    rows = []
    for scen in SCENARIOS:
        runs = metrics[(metrics['Scenario'] == scen) & (metrics['Org_AM'] == 'O3')]['run'].values - 1
        if len(runs) == 0 or not cons:
            continue
        rates = pd.DataFrame([{r: _ss_mean(all_flux_data[scen][i], r, t_span, steady_frac) for r in cons}
                              for i in runs]).median()
        tot = float(sum(rates[r] * net[r] for r in cons)) + 1e-12
        for r in cons:
            f = float(rates[r] * net[r])
            rows.append({'Escenario': scen, 'Especie': species, 'Reacción': r, 'Coef. neto consumido': net[r],
                         'Velocidad (mediana)': float(rates[r]), 'Flujo consumido': f, 'Fracción': f / tot})
    return pd.DataFrame(rows)


# ======================================================================
# 5. FIGURES AND TABLES
# ======================================================================
# ======================================================================
# CONFIGURATION (model-specific; everything below it is generic)
# ======================================================================
SCEN = SCENARIOS                                        # ['I', 'II', 'III', 'IV']
# Okabe-Ito based categorical palette, validated for colour-vision deficiency
# (adjacent-pair CVD dE >= 8.6, normal-vision dE >= 18.7 against a white background).
SC_COLOR = {'I': '#0072B2', 'II': '#009E73', 'III': '#E69F00', 'IV': '#B2182B'}
SC_NAME = {'I': 'I  C-limited mutualism', 'II': 'II  Strong mutualism',
           'III': 'III  Commensalism', 'IV': 'IV  Parasitism'}
SC_SHORT = {'I': 'I\n(N$-$, P$-$)', 'II': 'II\n(N+, P$-$)', 'III': 'III\n(N$-$, P+)', 'IV': 'IV\n(N+, P+)'}
JOHNSON_GRID = {'III': (0, 0), 'IV': (0, 1), 'I': (1, 0), 'II': (1, 1)}   # rows: P luxury, P limitation

STATE_EN = {'éxito': 'Success', 'fracaso': 'Failure', 'exclusión del hongo': 'Fungal exclusion (O$_2$)',
            'colapso del hospedero': 'Host collapse (O$_1$)'}
STATE_COLOR = {'éxito': '#4D4D4D', 'fracaso': '#BDBDBD', 'exclusión del hongo': '#9ECAE1',
               'colapso del hospedero': '#67000D'}
STATE_HATCH = {'éxito': '', 'fracaso': '', 'exclusión del hongo': '//', 'colapso del hospedero': 'xx'}

# Short English labels for the reactions of AMF_RN_MAK_v6 (unknown ids fall back to an automatic
# label built from the stoichiometry, so the script works with other networks).
REACTION_LABELS = {
    'r1': 'C fixation', 'r2': 'plant reactivation', 'r3': 'plant nutrient limitation',
    'r4': 'plant growth', 'r5': 'root N uptake', 'r6': 'root P uptake',
    'r7': 'active-plant recycling', 'r8': 'limited-plant recycling',
    'r9': 'hyphal N uptake', 'r10': 'hyphal P uptake', 'r11': 'extraradical mycelium growth',
    'r12': 'mycelium turnover', 'r13': 'C-for-N exchange', 'r14': 'N transfer to plant',
    'r15': 'C-for-P exchange', 'r16': 'intraradical fungal growth',
    'r17': 'P-feedback on colonisation', 'r18': 'fungal mortality',
    'r19': 'external N supply', 'r20': 'external P supply',
}
SPECIES_EN = {'Pact': 'active plant', 'Plim': 'limited plant', 'Cp': 'plant C', 'Np': 'plant N',
              'Pp': 'plant P', 'F': 'intraradical fungus', 'M': 'extraradical mycelium',
              'Cf': 'fungal C', 'Nf': 'fungal N', 'Pf': 'fungal P', 'Ns': 'soil N', 'Ps': 'soil P'}

# Scheme of Fig. 1a: compartments (axes-fraction boxes) and flows (start, end, label, reactions,
# curvature). Reaction ids are checked against the network file.
# 'lines': species legend of the box; 'notes': reactions whose stoichiometry is printed in the box
# (read from the network file, so the scheme always matches the model).
SCHEME_BOXES = {
    'Plant':  {'xy': (0.02, 0.50), 'wh': (0.34, 0.42), 'color': '#E5F5E0',
               'species': ['Pact', 'Plim', 'Cp', 'Np', 'Pp'],
               'lines': ['Pact, Plim: active, nutrient-limited plant', 'Cp, Np, Pp: plant C, N, P'],
               'notes': ['r1', 'r4']},
    'Fungus': {'xy': (0.64, 0.50), 'wh': (0.34, 0.42), 'color': '#FEE6CE',
               'species': ['F', 'M', 'Cf', 'Nf', 'Pf'],
               'lines': ['F, M: intraradical, extraradical fungus', 'Cf, Nf, Pf: fungal C, N, P'],
               'notes': ['r11', 'r16']},
    'Soil':   {'xy': (0.02, 0.02), 'wh': (0.96, 0.20), 'color': '#EFEDF5', 'species': ['Ns', 'Ps'],
               'lines': ['Ns, Ps: available soil N, P'], 'notes': []},
}
# (start, end, label, reactions, arc curvature[, label position]) in axes fractions
SCHEME_FLOWS = [
    ((0.36, 0.80), (0.64, 0.80), 'C', ['r13', 'r15'], -0.25, (0.50, 0.985)),
    ((0.64, 0.64), (0.36, 0.64), 'N, P', ['r13', 'r14', 'r15'], -0.25, (0.50, 0.475)),
    ((0.13, 0.22), (0.13, 0.50), 'N, P uptake', ['r5', 'r6'], 0.0),
    ((0.25, 0.50), (0.25, 0.22), 'recycling', ['r7', 'r8'], 0.0),
    ((0.75, 0.22), (0.75, 0.50), 'N, P uptake', ['r9', 'r10'], 0.0),
    ((0.87, 0.50), (0.87, 0.22), 'turnover', ['r12', 'r17', 'r18'], 0.0),
    ((0.19, 1.04), (0.19, 0.92), 'CO$_2$', ['r1'], 0.0),
    ((0.50, -0.09), (0.50, 0.02), 'external N, P supply', ['r19', 'r20'], 0.0),
]
PLANT_BIOMASS_SP = PLANT_BIOMASS            # ['Pact', 'Plim']
FUNGAL_BIOMASS_SP = FUNGAL_BIOMASS          # ['F', 'M']
SOIL_SP = ['Ns', 'Ps']
TOLERANCES = [0.025, 0.05, 0.075, 0.10]        # commensalism bands evaluated in Fig. 8 / Table 2
PRCC_MIN_ABS = 0.2                             # parameters shown in Fig. 5: |PRCC| >= this and p < 0.05
PRCC_TOP = 12                                  # at most this many parameters in Fig. 5 (full table in tablas/)
PRCC_ALWAYS = FLUX_GROUPS.get('N_fungal', [])   # always shown (fungal N delivery), even if PRCC ~ 0
PRCC_EXTRA = ['r20']                           # always shown in Fig. 5 (no dagger), e.g. external P supply
NP_THRESHOLDS = (14.0, 16.0)                   # tissue N:P below/above: N/P limitation (Koerselman & Meuleman 1996)
HARVEST_DAYS = [30, 60, 90]                    # experimental harvest times reported in Table 5 (days)
MGR_T_SKIP = 30                                # MGR(t) y-axis scaled to t >= this (days); transient clipped
# Fig. 2 rows: (label, species); every species of a row gets its own colour (MC_SPECIES_COLORS, in order)
MC_ROWS = [('Plant biomass', ['Pact', 'Plim']), ('Fungal biomass', ['F', 'M']), ('Carbon', ['Cp', 'Cf']),
           ('Nitrogen', ['Np', 'Nf']), ('Phosphorus', ['Pp', 'Pf']), ('Soil nutrients', ['Ns', 'Ps'])]
MC_SPECIES_COLORS = ['#08519C', '#B35806', '#542788']   # plant / fungus / third species of a row
MC_MGR_ROW = False                             # add the MGR(t) row to Fig. 2
FULL_W, ONE_W = 190.0, 90.0                    # Elsevier full-page and single-column widths (mm)

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 7, 'axes.labelsize': 7.5, 'axes.titlesize': 7.5, 'xtick.labelsize': 6.5,
    'ytick.labelsize': 6.5, 'legend.fontsize': 6.5, 'axes.linewidth': 0.6, 'lines.linewidth': 1.0,
    'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'axes.spines.top': False,
    'axes.spines.right': False, 'pdf.fonttype': 42, 'ps.fonttype': 42, 'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.02, 'axes.grid': False,
})


def mm(w, h):
    return (w / 25.4, h / 25.4)


def panel(ax, letter, x=-0.12, y=1.04):
    ax.text(x, y, f'({letter})', transform=ax.transAxes, fontsize=8.5, fontweight='bold', va='bottom', ha='left')


def reaction_label(r, stoich=None):
    if r in REACTION_LABELS:
        return f'{r} {REACTION_LABELS[r]}'
    if stoich and r in stoich:
        re_, pr = stoich[r]
        if is_self_limitation(stoich[r]):
            s = list(re_)[0]
            return f'{r} self-limitation {SPECIES_EN.get(s, s)}'
        return f'{r} {equation(r, stoich)}'
    return r



# ======================================================================
# OUTPUT FILES
# ======================================================================
class Outputs:
    """Numbered output files in <outdir>/figuras and <outdir>/tablas."""
    def __init__(self, root, formats, dpi):
        self.fig_dir, self.tab_dir = os.path.join(root, 'figuras'), os.path.join(root, 'tablas')
        os.makedirs(self.fig_dir, exist_ok=True); os.makedirs(self.tab_dir, exist_ok=True)
        self.formats, self.dpi, self.written = formats, dpi, []

    def figure(self, fig, n, name):
        for ext in self.formats:
            p = os.path.join(self.fig_dir, f'{n}_{name}.{ext}')
            fig.savefig(p, dpi=self.dpi if ext != 'pdf' else None)
            self.written.append(p)
        plt.close(fig)
        print(f'  Figure {n} -> {os.path.join(self.fig_dir, f"{n}_{name}")}.{{{",".join(self.formats)}}}')

    def table(self, n, name, ext='xlsx'):
        p = os.path.join(self.tab_dir, f'{n}_{name}.{ext}')
        self.written.append(p)
        return p


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def scenario_median_iqr(df, col):
    g = df.groupby('Scenario')[col]
    return g.median().reindex(SCEN), g.quantile(0.25).reindex(SCEN), g.quantile(0.75).reindex(SCEN)


# ======================================================================
# FIGURE 1  model structure + organisations
# ======================================================================
def simulate_organisations(network, t_end):
    """Three initial conditions, identical parameters (demonstration scenario ORG_DEMO_SCENARIO):
    bare soil, plant at its fungus-free equilibrium + sub-threshold inoculum, plant + inoculum."""
    rn = read_txt(network)
    sp = [s.name for s in rn.species()]
    rx = [r.name() if callable(r.name) else r.name for r in rn.reactions()]
    x0, p, rate = get_scenario_config(ORG_DEMO_SCENARIO, sp, rx)
    idx = {s: i for i, s in enumerate(sp)}
    sim = lambda x, T, n: simulation(rn, x0=x, rate=rate, spec_vector=p, t_span=(0, T), n_steps=n,
                                        method='LSODA', rtol=1e-8, atol=1e-12, verbose=False)[0]
    x_bare = list(x0)
    for s in PLANT_SPECIES + SYMBIONT:
        x_bare[idx[s]] = 0.0
    ts_nm = sim(make_nm_x0(x0, sp), t_end, 501)
    x_sub = [float(ts_nm[s].iloc[-1]) for s in sp]
    for s in FUNGAL_BIOMASS_SP:
        x_sub[idx[s]] = INOCULUM_SUBTHRESHOLD
    ics = [('Bare soil', x_bare), (f'Plant at equilibrium +\ninoculum {INOCULUM_SUBTHRESHOLD:g}', x_sub),
           ('Seedling + inoculum', list(x0))]
    orgs = compute_organizations(rn)
    res, rows = [], []
    for name, x in ics:
        ts = sim(x, t_end, 501)
        ts_long = sim([float(ts[s].iloc[-1]) for s in sp], ORG_CLASSIFY_T, 201)
        dest = classify_run(ts_long, (0, ORG_CLASSIFY_T))
        final = {s: float(ts_long[s].iloc[-1]) for s in sp}
        k, supp = _classify_support(final, orgs) if orgs else (None, frozenset(s for s, v in final.items() if v > ZERO_TOL))
        res.append((name, ts, dest))
        rows.append({'Initial condition': name.replace('\n', ' '), 'Destination (classifier)': dest,
                     'Matched organisation (pyCOT)': (f'#{k + 1} of {len(orgs)}' if k is not None else 'n/a'),
                     'Species present': ', '.join(sorted(supp))})
    org_rows = []
    if orgs:
        for i, O in enumerate(orgs):
            org_rows.append({'Organisation': f'#{i + 1}', 'Size': len(O), 'Species': ', '.join(sorted(O)),
                             'Contains previous': bool(i == 0 or orgs[i - 1] <= O)})
    return res, pd.DataFrame(rows), pd.DataFrame(org_rows), sp


def is_self_limitation(st):
    re_, pr = st
    return len(re_) == 1 and list(re_.values())[0] == 2 and pr == {list(re_)[0]: 1}


def equation(r, stoich):
    side = lambda d: ' + '.join(f"{'' if c == 1 else (int(c) if float(c).is_integer() else c)}{'' if c == 1 else ' '}{sp}"
                                for sp, c in d.items()) or '∅'
    re_, pr = stoich[r]
    return f'{side(re_)} → {side(pr)}'


def draw_scheme(ax, stoich):
    ax.set_xlim(0, 1); ax.set_ylim(-0.16, 1.08); ax.axis('off')
    for name, b in SCHEME_BOXES.items():
        (x, y), (w, h) = b['xy'], b['wh']
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.005,rounding_size=0.02',
                                    fc=b['color'], ec='#636363', lw=0.6))
        ax.text(x + w / 2, y + h - 0.02, name, ha='center', va='top', fontsize=8, fontweight='bold')
        lines = b.get('lines') or [f'{s}: {SPECIES_EN.get(s, s)}' for s in b['species']]
        sep = '\n' if h > 0.3 else '      '
        ax.text(x + w / 2, y + h - 0.105, sep.join(lines), ha='center', va='top', fontsize=6, linespacing=1.2)
        notes = [f'{r}: {equation(r, stoich)}' for r in b.get('notes', []) if r in stoich]
        if notes:
            ax.text(x + w / 2, y + 0.012, '\n'.join(notes), ha='center', va='bottom', fontsize=5.6,
                    linespacing=1.2, style='italic', color='#404040')
    ctrl = [r for r in stoich if is_self_limitation(stoich[r])]
    if ctrl:
        ax.text(0.5, -0.135, f'Every species X also has a self-limitation reaction 2X → X '
                f'({ctrl[0]}-{ctrl[-1]}).', ha='center', va='top', fontsize=5.8, color='#525252')
    missing = []
    for flow in SCHEME_FLOWS:
        (x0, y0), (x1, y1), lab, rxs, rad = flow[:5]
        label_xy = flow[5] if len(flow) > 5 else None
        missing += [r for r in rxs if r not in stoich]
        ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle='-|>', mutation_scale=8, lw=1.0,
                                     color='#252525', connectionstyle=f'arc3,rad={rad}', shrinkA=1, shrinkB=1))
        xm, ym = label_xy if label_xy else ((x0 + x1) / 2, (y0 + y1) / 2)
        horiz = abs(y1 - y0) < 1e-6
        one_line = horiz or abs(y1 - y0) < 0.15
        ax.text(xm + (0 if horiz else 0.012), ym, f'{lab} ({", ".join(rxs)})' if one_line else f'{lab}\n({", ".join(rxs)})',
                fontsize=5.8, ha='center' if horiz else 'left', va='center', color='#252525')
    if missing:
        print(f'  Warning: scheme reactions not found in the network: {sorted(set(missing))}')


def fig1(network, t_end, stoich, out):
    res, ic_table, org_table, sp = simulate_organisations(network, t_end)
    fig = plt.figure(figsize=mm(FULL_W, 135))
    gs = GridSpec(2, 3, height_ratios=[1.1, 1], hspace=0.38, wspace=0.28, figure=fig)
    ax = fig.add_subplot(gs[0, :]); draw_scheme(ax, stoich); panel(ax, 'a', x=0.0, y=0.98)
    org_lbl = {'O1': 'O$_1$ abiotic soil', 'O2': 'O$_2$ plant-soil', 'O3': 'O$_3$ full symbiosis'}
    axes = [fig.add_subplot(gs[1, j]) for j in range(3)]
    for j, (ax, (name, ts, dest)) in enumerate(zip(axes, res)):
        t = ts['Time'].values
        series = [('Plant biomass', PLANT_BIOMASS_SP, '#1B7837', '-'),
                  ('Fungal biomass', FUNGAL_BIOMASS_SP, '#B35806', '-'),
                  ('Soil N', ['Ns'], '#542788', '--'), ('Soil P', ['Ps'], '#8073AC', ':')]
        for lab, cols, c, ls in series:
            y = sum(ts[s].values for s in cols if s in ts.columns)
            ax.plot(t, y, color=c, ls=ls, lw=1.1, label=lab)
        ax.set_yscale('symlog', linthresh=1e-3); ax.set_ylim(0, None)
        ax.set_xlabel('Time (days)')
        ax.set_title(f'{name}\n→ {org_lbl[dest]}', fontsize=6.8)
        panel(ax, 'bcd'[j], x=-0.22)
    axes[0].set_ylabel('Concentration (symlog)')
    axes[2].legend(loc='center right', frameon=False)
    # fig.text(0.5, 0.005, f'Panels (b)-(d): identical parameters (scenario {ORG_DEMO_SCENARIO}); '
    #          'only the initial condition differs.', ha='center', fontsize=6, color='#525252')
    out.figure(fig, 1, 'Fig1_model_structure')
    with pd.ExcelWriter(out.table(3, 'Table_organizations'), engine='openpyxl') as w:
        ic_table.to_excel(w, sheet_name='Initial_conditions', index=False)
        (org_table if not org_table.empty else pd.DataFrame({'note': ['pyCOT organisation analysis unavailable']})) \
            .to_excel(w, sheet_name='Organizations_pyCOT', index=False)


# ======================================================================
# FIGURE 2  Monte Carlo dynamics
# ======================================================================
def fig_mc_dynamics(base, n_sims, pct, out, n):
    """Trajectories of the base Monte Carlo (the same runs summarised in every other figure)."""
    m = base['metrics']
    runs = {s: {'AM': base['ts_am'][s], 'NM': base['ts_nm'][s]} for s in SCEN}
    pairs = {s: {k: {'AM': a, 'NM': b} for k, (a, b) in enumerate(zip(base['ts_am'][s], base['ts_nm'][s]))}
             for s in SCEN}
    # MGR(t) = ln(B_AM(t) / B_NM(t)) for every paired run
    mgr_t, t_grid = {}, None
    for s in SCEN:
        M = []
        for k in sorted(pairs[s]):
            d = pairs[s][k]
            if 'AM' in d and 'NM' in d:
                bam = sum(d['AM'][q].values for q in PLANT_BIOMASS_SP)
                bnm = sum(d['NM'][q].values for q in PLANT_BIOMASS_SP)
                M.append(np.log((bam + 1e-12) / (bnm + 1e-12)))
                t_grid = d['AM']['Time'].values
        mgr_t[s] = np.vstack(M) if M else None
    # ---- figure: one row per species group (MC_ROWS), one column per scenario, AM runs ----
    n_rows = len(MC_ROWS) + (1 if MC_MGR_ROW else 0)
    fig, axes = plt.subplots(n_rows, len(SCEN), figsize=mm(FULL_W, 34 * n_rows + 15), sharex=True, squeeze=False)
    row_max = [0.0] * len(MC_ROWS)
    for j, s in enumerate(SCEN):
        orgs = m[m['Scenario'] == s]['Org_AM'].value_counts(normalize=True)
        am_runs = runs[s]['AM']
        for i, (lab, species) in enumerate(MC_ROWS):
            ax = axes[i, j]
            for k, sp in enumerate(species):
                Y = [ts[sp].values for ts in am_runs if sp in ts.columns]
                if not Y:
                    continue
                t = am_runs[0]['Time'].values
                c = MC_SPECIES_COLORS[k % len(MC_SPECIES_COLORS)]
                Ym = np.vstack(Y)
                row_max[i] = max(row_max[i], float(np.nanpercentile(Ym, 99.5)))
                for y in Y:
                    ax.plot(t, y, color=c, alpha=0.05, lw=0.4, rasterized=True)
                ax.plot(t, np.median(Ym, axis=0), color=c, lw=1.2,
                        label=f'{sp} ({SPECIES_EN.get(sp, sp)})')
            # for h in HARVEST_DAYS:
            #     ax.axvline(h, color='#969696', lw=0.4, ls=':', zorder=0)
            if i == 0:
                ax.set_title(f'{SC_NAME[s]}') #\nO$_3$ {100 * orgs.get("O3", 0):.0f} %, O$_2$ {100 * orgs.get("O2", 0):.0f} %, '
                            #  f'O$_1$ {100 * orgs.get("O1", 0):.0f} %', fontsize=6.5)
            if j == 0:
                ax.set_ylabel(lab)
            # if j == len(SCEN) - 1:
            #     ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1.0), frameon=False, fontsize=5.6,
            #               handlelength=1.4, borderaxespad=0)
            if j == 0:
                ax.legend(loc='upper left', frameon=True, framealpha=0.85, edgecolor='none',
                          fontsize=5.4, handlelength=1.2, borderaxespad=0.3, borderpad=0.3,
                          labelspacing=0.25)            
        panel(axes[0, j], 'abcd'[j], x=-0.2, y=1.22)
        if MC_MGR_ROW:
            ax = axes[len(MC_ROWS), j]
            if mgr_t[s] is not None:
                q1, med, q3 = np.percentile(mgr_t[s], [25, 50, 75], axis=0)
                ax.fill_between(t_grid, q1, q3, color=SC_COLOR[s], alpha=0.25, lw=0)
                ax.plot(t_grid, med, color=SC_COLOR[s], lw=1.2)
            ax.axhspan(-MGR_COMM_TOL, MGR_COMM_TOL, color='#D9D9D9', lw=0, zorder=0)
            ax.axhline(0, color='#525252', lw=0.5)
            # for h in HARVEST_DAYS:
            #     ax.axvline(h, color='#969696', lw=0.4, ls=':', zorder=0)
            if j == 0:
                ax.set_ylabel('MGR(t)')
        axes[-1, j].set_xlabel('Time (days)')
    for i in range(len(MC_ROWS)):        # same y scale within a row so scenarios are comparable
        for j in range(len(SCEN)):
            axes[i, j].set_ylim(0, 1.05 * row_max[i] if row_max[i] > 0 else 1)
            if j:
                axes[i, j].tick_params(labelleft=False)
    if MC_MGR_ROW:
        # y range from t >= MGR_T_SKIP: the initial transient would otherwise compress the values
        late = t_grid >= MGR_T_SKIP
        lo_ = min(np.nanmin(np.percentile(mgr_t[s][:, late], 25, axis=0)) for s in SCEN if mgr_t[s] is not None)
        hi_ = max(np.nanmax(np.percentile(mgr_t[s][:, late], 75, axis=0)) for s in SCEN if mgr_t[s] is not None)
        lo_, hi_ = min(lo_, -MGR_COMM_TOL), max(hi_, MGR_COMM_TOL)
        pad = 0.1 * (hi_ - lo_)
        for j in range(len(SCEN)):
            axes[-1, j].set_ylim(lo_ - pad, hi_ + pad)
            if j:
                axes[-1, j].tick_params(labelleft=False)
    handles = [Line2D([0], [0], color='#525252', lw=1.2, label='Median'),
               Line2D([0], [0], color='#969696', lw=0.6, alpha=0.6, label=f'Individual runs (n = {n_sims}, ±{pct:g} %)')]
    if MC_MGR_ROW:
        handles += [Patch(fc='#BDBDBD', label='MGR(t): median and IQR'),
                    Patch(fc='#D9D9D9', label=f'Commensalism band ±{MGR_COMM_TOL:g}')]
    fig.legend(handles=handles, loc='lower center', ncol=3, frameon=False, bbox_to_anchor=(0.45, 0.0))
    plt.tight_layout(rect=[0, 0.03, 1, 1], h_pad=0.6, w_pad=0.6)
    out.figure(fig, n, f'Fig{n}_Monte_Carlo_dynamics')
    # table: MGR at experimental harvest times and at the end of the horizon
    rows_t = []
    tol = MGR_COMM_TOL
    for s in SCEN:
        if mgr_t[s] is None:
            continue
        for h in HARVEST_DAYS + [float(t_grid[-1])]:
            v = mgr_t[s][:, int(np.argmin(np.abs(t_grid - h)))]
            rows_t.append({'Scenario': SC_NAME[s], 'Time (days)': h, 'MGR median': np.median(v),
                           'MGR Q1': np.percentile(v, 25), 'MGR Q3': np.percentile(v, 75),
                           'Mutualism, MGR > tol (%)': 100 * np.mean(v > tol),
                           'Commensalism, |MGR| <= tol (%)': 100 * np.mean(np.abs(v) <= tol),
                           'Parasitism, MGR < -tol (%)': 100 * np.mean(v < -tol)})
    pd.DataFrame(rows_t).round(4).to_excel(out.table(5, 'Table_MGR_harvest_times'), index=False)


# ======================================================================
# FIGURE 3  Johnson validation
# ======================================================================
def strip_box(ax, data, positions, colors, width=0.55):
    bp = ax.boxplot(data, positions=positions, widths=width, patch_artist=True, showfliers=False,
                    medianprops={'color': 'black', 'lw': 1.2}, whiskerprops={'lw': 0.6}, capprops={'lw': 0.6},
                    boxprops={'lw': 0.6})
    for patch, c in zip(bp['boxes'], colors):
        patch.set_facecolor(c); patch.set_alpha(0.25)
    rng = np.random.default_rng(1)
    for x, d, c in zip(positions, data, colors):
        ax.scatter(x + rng.uniform(-0.18, 0.18, len(d)), d, s=3, color=c, alpha=0.6, lw=0, zorder=3)


def fig2(base, out, n=3):
    m = base['metrics']; o3 = m[m['Org_AM'] == 'O3'].copy()
    tol = MGR_COMM_TOL
    fig = plt.figure(figsize=mm(FULL_W, 170))
    gs = GridSpec(3, 2, width_ratios=[1.25, 1], wspace=0.32, hspace=0.55, figure=fig)
    x = np.arange(len(SCEN)); cols = [SC_COLOR[s] for s in SCEN]

    # (a) MGR
    ax = fig.add_subplot(gs[0, 0])
    tol = MGR_COMM_TOL
    ax.axhspan(-tol, tol, color='#D9D9D9', lw=0, zorder=0)
    ax.axhline(0, color='#525252', lw=0.5)
    strip_box(ax, [o3.loc[o3['Scenario'] == s, 'MGR'].values for s in SCEN], x, cols)

    # for xi, s in zip(x, SCEN):
    #     ms = o3[o3['Scenario'] == s]  # CAMBIO: usar o3, no m
    #     ax.text(
    #         xi, 1.0,
    #         f"{100 * (ms['estado'] == 'éxito').mean():.0f} %",
    #         transform=ax.get_xaxis_transform(),
    #         ha='center', va='bottom', fontsize=6.5
    #     )

    # ax.text(
    #     1.02, 1.0, 'success',
    #     transform=ax.transAxes,
    #     ha='left', va='bottom',
    #     fontsize=6, color='#525252'
    # )  # CAMBIO: separar un poco del borde

    ax.set_xticks(x)
    ax.set_xticklabels([SC_SHORT[s] for s in SCEN])
    ax.set_ylabel('MGR = ln(B$_{AM}$ / B$_{NM}$)')
    ax.legend(
        handles=[Patch(fc='#D9D9D9', label=f'Commensalism band |MGR| ≤ {tol:g}')],
        loc='center right', bbox_to_anchor=(1.0, 1.02), frameon=False
    )
    panel(ax, 'a') 

    # (c) independent test: MGR vs tissue N:P of the plant grown without the fungus (Johnson 2010, Fig. 2)
    ax = fig.add_subplot(gs[1, 0])  # CAMBIO: ocupa filas 1 y 2 para eliminar el hueco

    if 'NP_plant_NM' in o3.columns and o3['NP_plant_NM'].notna().any():
        lo_np, hi_np = NP_THRESHOLDS

        # CAMBIO: filtrar > 0 para eje log y Spearman consistente
        dnp = o3[o3['NP_plant_NM'].notna() & (o3['NP_plant_NM'] > 0)].copy()

        if len(dnp):
            xs = dnp['NP_plant_NM']
            xmin = min(xs.min(), lo_np) * 0.7
            xmax = max(xs.max(), hi_np) * 1.4

            ax.axvspan(xmin, lo_np, color='#FDE0DD', lw=0, zorder=0)
            ax.axvspan(hi_np, xmax, color='#DEEBF7', lw=0, zorder=0)
            ax.axhspan(-tol, tol, color='#D9D9D9', lw=0, zorder=0)
            ax.axhline(0, color='#525252', lw=0.5)

            for s in SCEN:
                d = dnp[dnp['Scenario'] == s]
                ax.scatter(
                    d['NP_plant_NM'], d['MGR'],
                    s=4, color=SC_COLOR[s], alpha=0.6, lw=0, label=s
                )

            ax.set_xscale('log')
            ax.set_xlim(xmin, xmax)

            from scipy.stats import spearmanr
            rho, pval = spearmanr(dnp['NP_plant_NM'], dnp['MGR'])

            ax.text(
                0.72, 1.04,
                f'Spearman ρ = {rho:.2f} (p {"< 0.001" if pval < 0.001 else f"= {pval:.3f}"})',
                transform=ax.transAxes, va='bottom', ha='right', fontsize=6
            )

            ax.text(
                lo_np, 0.97, f'N-limited  < {lo_np:g} ',
                transform=ax.get_xaxis_transform(),
                ha='right', va='top', fontsize=5.0, color='#99000D'
            )
            ax.text(
                hi_np, 0.97, f' > {hi_np:g}  P-limited',
                transform=ax.get_xaxis_transform(),
                ha='left', va='top', fontsize=5.0, color='#08519C'
            )

            ax.set_xlabel('Tissue N:P of the plant without fungus')
            ax.set_ylabel('MGR')

            ax.legend(
                loc='upper left', bbox_to_anchor=(0.0, 1.0), frameon=False, markerscale=2.5, ncol=1,
                title='Scenario', title_fontsize=6,
                handletextpad=0.1, labelspacing=0.3, borderaxespad=0.3
            )

            o3_np = dnp.groupby('Scenario')['NP_plant_NM'].median().reindex(SCEN)
            print(
                '  Fig. 3c, median N:P of the NM plant: ' +
                ', '.join(f'{s} {v:.1f}' for s, v in o3_np.items())
            )
        else:
            ax.text(
                0.5, 0.5, 'No hay NP_plant_NM > 0',
                ha='center', va='center', transform=ax.transAxes
            )
    else:
        ax.text(
            0.5, 0.5, 'NP_plant_NM not found in the metrics table',
            ha='center', va='center', transform=ax.transAxes
        )

    panel(ax, 'c')

    # (b) Johnson 2x2
    sub = GridSpecFromSubplotSpec(2, 2, subplot_spec=gs[0:2, 1], hspace=0.45, wspace=0.18)

    o3['Cs'] = o3['C_fix'] / o3['C_fix'].max()
    o3['Cd'] = o3['C_to_fungus'] / o3['C_to_fungus'].max()
 
    keys = [('Cs', 'C$_{supply}$'), ('Cd', 'C$_{demand}$'), ('P_fungal_frac', 'P$_{benefit}$')]
    bar_c = ['#1B7837', '#B35806', '#2166AC']

    first = None
    for s, (r, c) in JOHNSON_GRID.items():
        ax = fig.add_subplot(sub[r, c], sharey=first)
        if first is None:
            first = ax

        d = o3[o3['Scenario'] == s]
        med = [d[k].median() for k, _ in keys]
        q1 = [d[k].quantile(0.25) for k, _ in keys]
        q3 = [d[k].quantile(0.75) for k, _ in keys]

        ax.bar(
            range(3), med, color=bar_c, width=0.68,
            yerr=[np.subtract(med, q1), np.subtract(q3, med)],
            error_kw={'lw': 0.6, 'capsize': 1.5}
        )

        for i, (v, top) in enumerate(zip(med, q3)):
            ax.text(i, top + 0.03, f'{v:.2f}', ha='center', va='bottom', fontsize=5.6)

        ax.set_xticks(range(3))
        ax.set_xticklabels([l for _, l in keys], fontsize=5.8)
        ax.set_ylim(0, 1.12)
        ax.set_title(SC_NAME[s], fontsize=5.5, fontweight='bold')

        if c == 1:
            plt.setp(ax.get_yticklabels(), visible=False)

        if r == 0 and c == 0:
            panel(ax, 'b', x=-0.32)

    # first.set_ylabel('Relative value')  # CAMBIO: etiqueta del eje y de (b)

    fig.text(0.755, 0.355, 'N supply: limitation  →  luxury', ha='center', fontsize=6.5)
    fig.text(0.545, 0.66, 'P supply: limitation  →  luxury',
             va='center', rotation=90, fontsize=6.5)

    # fig.text(
    #     0.785, 0.29,
    #     'Bars: median over O$_3$ runs; whiskers: IQR.\n'
    #     'C supply and C demand are normalised by\n'
    #     'their maximum over all O$_3$ runs (relative\n'
    #     'values, as in Johnson\'s conceptual model).',
    #     ha='center', va='top', fontsize=6, color='#525252'
    # )

    out.figure(fig, n, f'Fig{n}_Johnson_validation')


# ======================================================================
# FIGURE 4  phase diagram
# ======================================================================
def fig3(base, out, n=4):
    ph = base['phase']
    if ph is None or ph.empty:
        print('  Phase diagram skipped: not computed.')
        return
    gN, gP = np.sort(ph['k_N'].unique()), np.sort(ph['k_P'].unique())
    piv = lambda c: ph.pivot(index='k_P', columns='k_N', values=c).reindex(index=gP, columns=gN).values
    M, ORG, NM = piv('MGR').astype(float), piv('Org_AM'), piv('Org_NM')
    M = np.where(ORG == 'O3', M, np.nan)
    tol = MGR_COMM_TOL
    vmax = max(float(np.nanpercentile(np.abs(M), 95)) if np.isfinite(M).any() else 1.0, 2 * tol)
    fig, ax = plt.subplots(figsize=mm(140, 112))
    pc = ax.pcolormesh(gN, gP, np.clip(M, -vmax, vmax), cmap='RdBu', norm=TwoSlopeNorm(0, -vmax, vmax),
                       shading='nearest', rasterized=True)
    induced = (ORG == 'O1') & (NM != 'O1')
    cat = np.where(ORG == 'O2', 0.0, np.where(induced, 1.0, np.where(ORG == 'O1', 2.0, np.nan))).astype(float)
    cat_c = ['#D9D9D9', '#67000D', '#525252']
    ax.pcolormesh(gN, gP, cat, cmap=ListedColormap(cat_c), vmin=0, vmax=2, shading='nearest', rasterized=True)
    XX, YY = np.meshgrid(gN, gP)
    if np.isfinite(M).sum() > 4:
        ax.contour(XX, YY, M, levels=[0.0], colors='black', linewidths=1.2)
        ax.contour(XX, YY, M, levels=[-tol, tol], colors='black', linewidths=0.6, linestyles='--')
    for s in SCEN:
        kn = SCENARIO_SUPPLY[s][SUPPLY_REACTIONS['N']][0]; kp = SCENARIO_SUPPLY[s][SUPPLY_REACTIONS['P']][0]
        ax.plot(kn, kp, 'o', ms=6, mfc=SC_COLOR[s], mec='white', mew=1.0, zorder=5)
        ax.annotate(s, (kn, kp), xytext=(5, 4), textcoords='offset points', fontsize=7.5, fontweight='bold',
                    color='black', zorder=6)
    ax.set_xscale('log'); ax.set_yscale('log')
    rN, rP = SUPPLY_REACTIONS['N'], SUPPLY_REACTIONS['P']
    ax.set_xlabel(f'External N supply, k$_{{{rN[1:]}}}$'); ax.set_ylabel(f'External P supply, k$_{{{rP[1:]}}}$')
    cb = fig.colorbar(pc, ax=ax, extend='both', fraction=0.045, pad=0.02)
    cb.set_label('MGR in coexistence (O$_3$)'); cb.outline.set_linewidth(0.5)
    handles = [Line2D([0], [0], color='black', lw=1.2, label='MGR = 0'),
               Line2D([0], [0], color='black', lw=0.6, ls='--', label=f'MGR = ±{tol:g}'),
               Patch(fc=cat_c[0], label='Fungal exclusion (O$_2$)'),
               Patch(fc=cat_c[1], label='Fungus-induced host collapse (O$_1$)'),
               Patch(fc=cat_c[2], label='Plant not viable without fungus (O$_1$)')]
    ax.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, -0.13), ncol=3, frameon=False)
    out.figure(fig, n, f'Fig{n}_phase_diagram')


# ======================================================================
# FIGURE 5  PRCC
# ======================================================================
def fig4(base, stoich, out, n=5, outputs=('MGR', 'Fungal_biomass')):
    if not base.get('prcc'):
        print('  PRCC figure skipped: PRCC not computed.')
        return
    df = list(base['prcc'].values())[0]
    outputs = [o for o in outputs if o in set(df['Salida'])]
    if not outputs:
        print('  PRCC figure skipped: PRCC table has none of the requested outputs.'); return
    titles = {'MGR': 'MGR (O$_3$ runs)', 'Fungal_biomass': 'Fungal biomass (all runs)', 'C_demand': 'C demand (O$_3$ runs)'}
    sig = df[(df['PRCC'].abs() >= PRCC_MIN_ABS) & (df['p'] < 0.05) & df['Salida'].isin(outputs)]
    ranking = sig.groupby('Parámetro')['PRCC'].apply(lambda v: v.abs().max()).sort_values(ascending=False)
    keep = list(ranking.index[:PRCC_TOP])
    forced = [r for r in PRCC_ALWAYS if r in set(df['Parámetro']) and r not in keep]
    extra = [r for r in PRCC_EXTRA if r in set(df['Parámetro']) and r not in keep + forced]
    keep += forced + extra
    if not keep:
        print('  PRCC figure skipped: no parameter passes the PRCC threshold.'); return
    keep = sorted(keep, key=lambda r: int(r[1:]) if r[1:].isdigit() else 999)
    n_total = df['Parámetro'].nunique()
    fig, axes = plt.subplots(1, len(outputs), figsize=mm(FULL_W, 20 + 5.2 * len(keep)), squeeze=False)
    for j, (ax, o) in enumerate(zip(axes[0], outputs)):
        d = df[df['Salida'] == o]
        P = d.pivot(index='Parámetro', columns='Escenario', values='PRCC').reindex(index=keep, columns=SCEN)
        Q = d.pivot(index='Parámetro', columns='Escenario', values='p').reindex(index=keep, columns=SCEN)
        im = ax.imshow(P.values, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
        for i in range(P.shape[0]):
            for k in range(P.shape[1]):
                v, pv = P.values[i, k], Q.values[i, k]
                if np.isfinite(v):
                    star = '**' if pv < 0.001 else ('*' if pv < 0.05 else '')
                    ax.text(k, i, f'{v:.2f}{star}', ha='center', va='center', fontsize=5.6,
                            color='white' if abs(v) > 0.6 else 'black')
        for k, s in enumerate(SCEN):
            if not np.isfinite(P[s].values).any():
                ax.text(k, (len(keep) - 1) / 2, 'n too small', rotation=90, ha='center', va='center', fontsize=5.6,
                        color='#737373')
        ax.set_xticks(range(len(SCEN))); ax.set_xticklabels(SCEN)
        ax.set_yticks(range(len(keep)))
        ax.set_yticklabels([reaction_label(r, stoich) + (' †' if r in forced else '') for r in keep] if j == 0 else [])
        ax.tick_params(length=0); [s.set_visible(False) for s in ax.spines.values()]
        ax.set_title(titles.get(o, o)); panel(ax, 'ab'[j], x=-0.02 if j else -0.75)
    cb = fig.colorbar(im, ax=axes[0].tolist(), fraction=0.03, pad=0.02); cb.set_label('PRCC'); cb.outline.set_linewidth(0.5)
    # fig.text(0.0, -0.01, f'Shown: the {min(PRCC_TOP, len(ranking))} of {n_total} parameters with the largest |PRCC| '
    #          f'(≥ {PRCC_MIN_ABS}, p < 0.05 in at least one scenario-output pair)'
    #          + ('; † shown for reference (fungal N delivery)' if forced else '')
    #          + '. * p < 0.05; ** p < 0.001. Full table: tablas/6_Table_PRCC_full.xlsx.', fontsize=5.8, color='#525252')
    df.to_excel(out.table(6, 'Table_PRCC_full'), index=False)
    out.figure(fig, n, f'Fig{n}_PRCC')


# ======================================================================
# FIGURE 6  N-P trade-off
# ======================================================================
def state_bars(ax, metrics, x, width, label_top=None):
    for i, s in enumerate(SCEN):
        ms = metrics[metrics['Scenario'] == s]
        fr = ms['estado'].value_counts(normalize=True)
        bottom = 0.0
        for st in STATE_ORDER:
            v = float(fr.get(st, 0.0))
            if v > 0:
                ax.bar(x[i], v, width, bottom=bottom, color=STATE_COLOR[st], hatch=STATE_HATCH[st],
                       edgecolor='white', lw=0.5)
            bottom += v
        if label_top:
            ax.text(x[i], 1.02, label_top, ha='center', va='bottom', fontsize=5.4, rotation=90)


def fig5(base, alt, alt_label, base_label, stoich, out, n=6):
    if alt is None:
        print('  Trade-off figure skipped: alternative calibration not run.'); return
    fig = plt.figure(figsize=mm(FULL_W, 70))
    gs = GridSpec(1, 3, width_ratios=[1.1, 1, 1], wspace=0.42, figure=fig)
    # (a) outcomes
    ax = fig.add_subplot(gs[0, 0]); x = np.arange(len(SCEN)); w = 0.36
    state_bars(ax, base['metrics'], x - w / 2 - 0.02, w)
    state_bars(ax, alt['metrics'], x + w / 2 + 0.02, w)
    for i in range(len(SCEN)):
        ax.text(x[i] - w / 2 - 0.02, -0.04, 'B', ha='center', va='top', fontsize=5.8, transform=ax.get_xaxis_transform())
        ax.text(x[i] + w / 2 + 0.02, -0.04, 'N', ha='center', va='top', fontsize=5.8, transform=ax.get_xaxis_transform())
    ax.set_xticks(x); ax.set_xticklabels(SCEN); ax.tick_params(axis='x', pad=9)
    ax.set_ylim(0, 1); ax.set_ylabel('Fraction of Monte Carlo runs')
    ax.legend(handles=[Patch(fc=STATE_COLOR[s], hatch=STATE_HATCH[s], ec='white', label=STATE_EN[s]) for s in STATE_ORDER],
              loc='upper center', bbox_to_anchor=(0.5, -0.17), ncol=2, frameon=False, fontsize=5.8)
    ax.set_title(f'B: {base_label}   N: {alt_label}', fontsize=6.2)
    panel(ax, 'a', x=-0.28)
    # (b) N benefit vs P benefit
    ax = fig.add_subplot(gs[0, 1])
    ax.axvspan(*N_BENEFIT_LIT, color='#F0F0F0', lw=0, zorder=0)
    ax.axhline(P_BENEFIT_HIGH, color='#969696', lw=0.6, ls=':')
    for s in SCEN:
        pts = []
        for res, mk in ((base, 'o'), (alt, 's')):
            o3 = res['metrics'][(res['metrics']['Org_AM'] == 'O3') & (res['metrics']['Scenario'] == s)]
            if o3.empty:
                continue
            xm, ym = o3['N_fungal_frac'].median(), o3['P_fungal_frac'].median()
            ax.plot(xm, ym, mk, ms=5, mfc=SC_COLOR[s] if mk == 'o' else 'white', mec=SC_COLOR[s], mew=1.2, zorder=3)
            pts.append((xm, ym))
        if len(pts) == 2:
            ax.annotate('', xy=pts[1], xytext=pts[0], arrowprops={'arrowstyle': '-|>', 'color': SC_COLOR[s], 'lw': 0.8,
                                                                  'shrinkA': 4, 'shrinkB': 4})
            # ax.annotate(s, pts[0], xytext=(-9, 2), textcoords='offset points', fontsize=6.5, color=SC_COLOR[s],
            #             fontweight='bold')
            ax.annotate(
                s, pts[0],
                xytext=(0, 7), textcoords='offset points',
                ha='center', va='bottom',
                fontsize=6.5, color=SC_COLOR[s], fontweight='bold',
                zorder=4
            )            
    ax.set_xscale('symlog', linthresh=1e-3); ax.set_xlim(0, 0.5)
    ax.set_xlabel('N benefit (N via fungus / plant N uptake)'); ax.set_ylabel('P benefit (P via fungus / plant P uptake)')
    ax.legend(handles=[Line2D([0], [0], marker='o', ls='', mfc='#525252', mec='#525252', label=base_label),
                       Line2D([0], [0], marker='s', ls='', mfc='white', mec='#525252', label=alt_label)],
              loc='upper center', bbox_to_anchor=(0.5, -0.17), ncol=1, frameon=False, fontsize=5.8)
    ax.text(N_BENEFIT_LIT[0], 0.02, ' N reference\n band', transform=ax.get_xaxis_transform(), fontsize=5.5,
            va='bottom', color='#737373')
    panel(ax, 'b', x=-0.25)
    # (c) fate of fungal N
    ax = fig.add_subplot(gs[0, 2])
    nt = alt.get('ntrade')
    fate = next((v for k, v in (nt or {}).items() if k.startswith('Destino')), None)
    if fate is not None and not fate.empty:
        piv = fate.pivot(index='Reacción', columns='Escenario', values='Fracción').reindex(columns=SCEN).fillna(0)
        main = piv[piv.max(axis=1) >= 0.03].copy()
        if len(main) < len(piv):
            main.loc['other'] = piv[piv.max(axis=1) < 0.03].sum()
        main = main.loc[main.mean(axis=1).sort_values(ascending=False).index]
        trade = set(FLUX_GROUPS.get('N_fungal', []))
        grey = ['#252525', '#636363', '#969696', '#BDBDBD', '#D9D9D9']
        bottom = np.zeros(len(SCEN)); gi = 0
        for r, vals in main.iterrows():
            c = '#B35806' if r in trade else grey[min(gi, len(grey) - 1)]
            gi += r not in trade
            ax.bar(np.arange(len(SCEN)), vals.values, 0.62, bottom=bottom, color=c, edgecolor='white', lw=0.5,
                   label=reaction_label(r, stoich) if r != 'other' else 'other')
            bottom += vals.values
        ax.set_ylim(0, 1); ax.set_xticks(range(len(SCEN))); ax.set_xticklabels(SCEN)
        ax.set_ylabel(f'Fate of fungal N ({FUNGAL_N}), fraction')
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12), ncol=1, frameon=False, fontsize=5.6)
        ax.set_title(alt_label, fontsize=6.2)
    else:
        ax.text(0.5, 0.5, 'Fate of fungal N not computed', ha='center', va='center', transform=ax.transAxes)
    panel(ax, 'c', x=-0.25)
    out.figure(fig, n, f'Fig{n}_N_P_tradeoff')


# ======================================================================
# FIGURE 7  supply sweeps
# ======================================================================
def sweep_panels(ax_m, ax_s, sw, base_level, xlabel):
    rs, R = sw['summary'], sw['reaction']
    swept = [s for s in SCEN if rs[rs['Escenario'] == s][f'valor_{R}'].nunique() > 1]
    for s in SCEN:
        d = rs[rs['Escenario'] == s].sort_values('nivel')
        ls = '-' if s in swept else ':'
        ax_m.plot(d['nivel'], d['MGR_med'], ls, marker='o', ms=3, color=SC_COLOR[s], label=s)
        ax_m.fill_between(d['nivel'], d['MGR_q25'], d['MGR_q75'], color=SC_COLOR[s], alpha=0.15, lw=0)
        ax_s.plot(d['nivel'], d['% éxito'], ls, marker='o', ms=3, color=SC_COLOR[s], label=s)
    ax_m.axhspan(-MGR_COMM_TOL, MGR_COMM_TOL, color='#D9D9D9', lw=0, zorder=0)
    ax_m.axhline(0, color='#525252', lw=0.5)
    lv = np.sort(rs['nivel'].unique())
    for a in (ax_m, ax_s):
        if lv.max() / max(lv.min(), 1e-12) > 4:
            a.set_xscale('log')
        a.set_xticks(lv); a.set_xticklabels([f'{v:g}' for v in lv]); a.minorticks_off()
        if base_level is not None:
            a.axvline(base_level, color='#525252', lw=0.6, ls='--')
        a.set_xlabel(xlabel)
    ax_m.set_ylabel('MGR (median, IQR)'); ax_s.set_ylabel('Successful runs (%)'); ax_s.set_ylim(-3, 103)
    return swept


def fig6(sweeps, out, n=7):
    sweeps = [s for s in sweeps if s[0] is not None]
    if not sweeps:
        print('  Supply figure skipped: no supply sweep run.'); return
    fig, axes = plt.subplots(2, len(sweeps), figsize=mm(FULL_W if len(sweeps) > 1 else ONE_W + 20, 110), squeeze=False)
    letters = iter('abcd')
    for j, (sw, base_level, xlabel) in enumerate(sweeps):
        swept = sweep_panels(axes[0, j], axes[1, j], sw, base_level, xlabel)
        axes[0, j].set_title(f'{sw["reaction"]} varied in {", ".join(swept)} (solid lines)', fontsize=6.5)
    for j in range(len(sweeps)):
        panel(axes[0, j], next(letters), x=-0.2)
    for j in range(len(sweeps)):
        panel(axes[1, j], next(letters), x=-0.2)
    plt.tight_layout(rect=[0, 0.06, 1, 1])
    fig.legend(handles=[Line2D([0], [0], color=SC_COLOR[s], marker='o', ms=3, label=SC_NAME[s]) for s in SCEN],
               loc='lower center', ncol=4, frameon=False, bbox_to_anchor=(0.5, 0.0))
    # fig.text(0.5, -0.025, 'Solid: scenarios whose supply is varied; dotted: unchanged scenarios. '
    #          'Dashed vertical line: base case. Grey band: commensalism tolerance.', ha='center', fontsize=6, color='#525252')
    out.figure(fig, n, f'Fig{n}_supply_sensitivity')


# ======================================================================
# FIGURE 8  robustness + Tables 1-2
# ======================================================================
def success_vs_tolerance(metrics):
    global MGR_COMM_TOL       # evaluate_success() reads the module-level tolerance
    saved = MGR_COMM_TOL; rows = []
    try:
        for tol in TOLERANCES:
            MGR_COMM_TOL = tol
            st = evaluate_success(metrics.drop(columns=['estado'], errors='ignore'))
            for s in SCEN:
                ss = st[st['Scenario'] == s]
                rows.append({'Scenario': s, 'Tolerance': tol, 'Success (%)': 100 * (ss['estado'] == 'éxito').mean()})
    finally:
        MGR_COMM_TOL = saved
    return pd.DataFrame(rows)


def fig7(base, out, n=8):
    m = base['metrics']
    tt = success_vs_tolerance(m)
    fig, axes = plt.subplots(1, 2, figsize=mm(FULL_W, 62), gridspec_kw={'width_ratios': [1, 1.1], 'wspace': 0.35})
    ax = axes[0]; x = np.arange(len(SCEN))
    state_bars(ax, m, x, 0.6)
    for i, s in enumerate(SCEN):
        v = 100 * (m[m['Scenario'] == s]['estado'] == 'éxito').mean()
        ax.text(x[i], 1.01, f'{v:.0f} %', ha='center', va='bottom', fontsize=6.3)
    ax.set_xticks(x); ax.set_xticklabels([SC_SHORT[s] for s in SCEN]); ax.set_ylim(0, 1.08)
    ax.set_ylabel(f'Fraction of runs (n = {int(m.groupby("Scenario").size().min())} per scenario)')
    ax.legend(handles=[Patch(fc=STATE_COLOR[s], hatch=STATE_HATCH[s], ec='white', label=STATE_EN[s]) for s in STATE_ORDER],
              loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2, frameon=False, fontsize=6)
    panel(ax, 'a', x=-0.2)
    ax = axes[1]
    for s in SCEN:
        d = tt[tt['Scenario'] == s]
        ax.plot(d['Tolerance'], d['Success (%)'], marker='o', ms=3.5, color=SC_COLOR[s], label=SC_NAME[s])
    ax.axvline(MGR_COMM_TOL, color='#525252', lw=0.6, ls='--')
    ax.set_xlabel('Commensalism tolerance |MGR| ≤ tol'); ax.set_ylabel('Successful runs (%)'); ax.set_ylim(-3, 103)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2, frameon=False)
    panel(ax, 'b', x=-0.18)
    out.figure(fig, n, f'Fig{n}_robustness')
    return tt


def table_summary(base, out, tt):
    m = base['metrics']; o3 = m[m['Org_AM'] == 'O3']
    rows = []
    fmt = lambda s, c: (f"{o3[o3['Scenario'] == s][c].median():.3g} "
                        f"[{o3[o3['Scenario'] == s][c].quantile(0.25):.3g}, {o3[o3['Scenario'] == s][c].quantile(0.75):.3g}]")
    for s in SCEN:
        ms = m[m['Scenario'] == s]; fr = ms['estado'].value_counts(normalize=True)
        row = {'Scenario': SC_NAME[s], 'Runs': len(ms)}
        row.update({f'{STATE_EN[k].split(" (")[0]} (%)': round(100 * float(fr.get(k, 0)), 1) for k in STATE_ORDER})
        for c, lab in (('MGR', 'MGR'), ('C_demand', 'C demand (fraction of fixed C)'), ('C_to_fungus', 'C sink (absolute)'),
                       ('P_fungal_frac', 'P benefit'), ('N_fungal_frac', 'N benefit')):
            row[f'{lab}, median [IQR]'] = fmt(s, c)
        rows.append(row)
    with pd.ExcelWriter(out.table(1, 'Table_scenario_summary'), engine='openpyxl') as w:
        pd.DataFrame(rows).to_excel(w, sheet_name='Summary', index=False)
        pd.DataFrame({'Note': ['Medians and interquartile ranges over Monte Carlo runs ending in coexistence (O3); '
                               'outcome percentages over all runs. Success follows the criteria of evaluate_success() '
                               'in figuras_articulo.py.']}).to_excel(w, sheet_name='Notes', index=False)
    tt.pivot(index='Scenario', columns='Tolerance', values='Success (%)').reindex(SCEN).round(1) \
      .to_excel(out.table(2, 'Table_success_tolerance'))



# ======================================================================
# MAIN
# ======================================================================
def main():
    start_time = time.perf_counter()

    # ---- RUN SETTINGS (fixed; the article uses these values) ----
    t_end = 100            # horizon of the Monte Carlo, the sweeps and Fig. 1 (days)
    pct_var = 10           # Monte Carlo perturbation of every initial condition and rate constant (+/- %)
    n_sims = 100           # Monte Carlo runs per scenario (each paired with its non-mycorrhizal control)
    random_seed = 123      # seed of every random draw
    SWEEPS = [             # supply sweeps of Fig. 7 (base calibration)
        {'nutrient': 'P', 'reaction': 'r20', 'scenarios': ['III', 'IV'], 'mode': 'values',
         'levels': [50, 73.75, 100, 150, 200, 300]},
        {'nutrient': 'N', 'reaction': 'r19', 'scenarios': ['II', 'IV'], 'mode': 'factors',
         'levels': [0.5, 0.75, 1, 1.25, 1.5, 2, 3]},
    ]
    skip = []              # figures to skip, e.g. [7]; their computations are skipped as well
    formats, dpi = ['pdf', 'png'], 600
    save_dir = os.path.join(_pycot_root, 'projects', 'AMF_model', 'outputs_figures_article')

    # ---- network ----
    try:
        rn = read_txt(NETWORK_FILE)
    except FileNotFoundError:
        sys.exit(f'Network file not found: {NETWORK_FILE}')
    species_order, reaction_order = network_orders(rn)
    stoich = parse_network_stoich(NETWORK_FILE)
    missing = [r for r in reaction_order if r not in PARAMS_BASE or r not in PARAMS_ALT]
    if missing:
        sys.exit(f'Reactions without a rate constant in section 1: {missing}')
    print(f'Network: {NETWORK_FILE}\nSpecies: {species_order}\nReactions: {reaction_order}')
    print(f'Settings: t_end = {t_end} d, +/- {pct_var} %, n = {n_sims} runs per scenario, seed = {random_seed}')
    out = Outputs(save_dir, formats, dpi)

    # ---- 1. Monte Carlo, phase diagram and PRCC of the base calibration ----
    print('\n' + '=' * 70 + '\n 1. BASE CALIBRATION: MONTE CARLO, PHASE DIAGRAM, PRCC\n' + '=' * 70)
    set_calibration('base')
    base = run_monte_carlo(rn, n_sims, pct_var, random_seed, t_end, desc='Monte Carlo (base)')
    print('\nOutcome of the runs per scenario (base calibration):')
    print(base['metrics'].groupby('Scenario')['estado'].value_counts(normalize=True).unstack(fill_value=0)
          .reindex(index=SCENARIOS, columns=STATE_ORDER, fill_value=0).round(2).rename(columns=STATE_EN))
    print(f"Non-converged runs: {int((~base['metrics']['convergido']).sum())} of {len(base['metrics'])}")
    if 4 not in skip:
        base['phase'] = compute_phase_diagram(rn, t_end)
    if 5 not in skip:
        base['prcc'] = {'PRCC': compute_prcc(base['metrics'], base['cond'], reaction_order)}

    # ---- 2. Monte Carlo of the alternative calibration (Fig. 6) ----
    alt = None
    if 6 not in skip:
        print('\n' + '=' * 70 + '\n 2. ALTERNATIVE CALIBRATION (forced N trade): MONTE CARLO\n' + '=' * 70)
        set_calibration('alt')
        alt = run_monte_carlo(rn, n_sims, pct_var, random_seed, t_end, desc='Monte Carlo (alt)')
        alt['ntrade'] = {f'Destino_{FUNGAL_N}': compute_species_fate(stoich, alt['metrics'], alt['fl_am'],
                                                                     FUNGAL_N, (0, t_end))}
        set_calibration('base')

    # ---- 3. Supply sweeps of the base calibration (Fig. 7) ----
    sweeps = {}
    if 7 not in skip:
        print('\n' + '=' * 70 + '\n 3. SUPPLY SWEEPS (base calibration)\n' + '=' * 70)
        for sw in SWEEPS:
            sweeps[sw['nutrient']] = run_supply_sweep(rn, sw['reaction'], sw['scenarios'], sw['levels'], sw['mode'],
                                                      n_sims, pct_var, random_seed, t_end)

    def sweep_axis(sw, nutrient):
        """Base-case level and x label of a sweep (absolute values or factors of the base supply)."""
        if sw is None:
            return None, ''
        R = sw['reaction']
        if sw['mode'] == 'values':
            s = next(s for s in SCEN if R in SCENARIO_SUPPLY[s] and
                     sw['summary'][sw['summary']['Escenario'] == s][f'valor_{R}'].nunique() > 1)
            return SCENARIO_SUPPLY[s][R][0], f'External {nutrient} supply, k$_{{{R[1:]}}}$'
        return 1.0, f'{nutrient} supply (factor of base k$_{{{R[1:]}}}$)'

    # ---- 4. Figures and tables ----
    print('\n' + '=' * 70 + '\n 4. FIGURES AND TABLES\n' + '=' * 70)
    if 1 not in skip:
        fig1(NETWORK_FILE, t_end, stoich, out)
    if 2 not in skip:
        fig_mc_dynamics(base, n_sims, pct_var, out, 2)
    if 3 not in skip:
        fig2(base, out, 3)
    if 4 not in skip:
        fig3(base, out, 4)
    if 5 not in skip:
        fig4(base, stoich, out, 5)
    if 6 not in skip:
        fig5(base, alt, CALIBRATIONS['alt']['label'], CALIBRATIONS['base']['label'], stoich, out, 6)
    if 7 not in skip:
        fig6([(sweeps.get(k), *sweep_axis(sweeps.get(k), k)) for k in ('P', 'N')], out, 7)
    tt = fig7(base, out, 8) if 8 not in skip else success_vs_tolerance(base['metrics'])
    table_summary(base, out, tt)

    # Table 7: every Monte Carlo run (base and alternative) and the phase-diagram grid
    with pd.ExcelWriter(out.table(7, 'Table_MC_runs'), engine='openpyxl') as w:
        for name, res in (('base', base), ('alt', alt)):
            if res is None:
                continue
            m = res['metrics'].copy()
            m['Outcome'] = m['estado'].map(lambda v: STATE_EN[v].split(' (')[0])
            m.drop(columns=['estado']).to_excel(w, sheet_name=f'Runs_{name}', index=False)
            res['success'].to_excel(w, sheet_name=f'Criteria_{name}', index=False)
        if base['phase'] is not None:
            base['phase'].to_excel(w, sheet_name='Phase_diagram', index=False)
        if alt is not None:
            alt['ntrade'][f'Destino_{FUNGAL_N}'].to_excel(w, sheet_name=f'Fate_{FUNGAL_N}_alt', index=False)
    # Table 8: supply sweeps
    if sweeps:
        with pd.ExcelWriter(out.table(8, 'Table_supply_sweeps'), engine='openpyxl') as w:
            for k, sw in sweeps.items():
                sw['summary'].to_excel(w, sheet_name=f'{k}_supply_{sw["reaction"]}', index=False)
            pd.DataFrame(SWEEPS).to_excel(w, sheet_name='Settings', index=False)

    # Table 4: reproducibility record
    import scipy
    versions = {'python': platform.python_version(), 'numpy': np.__version__, 'pandas': pd.__version__,
                'scipy': scipy.__version__, 'matplotlib': matplotlib.__version__}
    try:
        import pyCOT
        versions['pyCOT'] = getattr(pyCOT, '__version__', 'unknown (local source)')
    except Exception:
        pass
    repro = {'script': os.path.basename(__file__), 'script_sha256': sha256(os.path.abspath(__file__)),
             'date': time.strftime('%Y-%m-%d %H:%M'), 'versions': versions, 'platform': platform.platform(),
             'network': {'path': NETWORK_FILE, 'sha256': sha256(NETWORK_FILE)},
             'settings': {'t_end': t_end, 'pct': pct_var, 'n_sims': n_sims, 'seed': random_seed,
                          'steady_frac': STEADY_FRAC, 'sim_tolerances': SIM_TOL,
                          'phase': {'n_grid': PHASE_N_GRID, 'span': PHASE_SPAN, 'tmult': PHASE_TMULT,
                                    'tolerances': SIM_TOL_PHASE},
                          'sweeps': SWEEPS, 'harvest_days': HARVEST_DAYS, 'skip': skip},
             'initial_condition': {**x0_dict_base, 'soil': f'{SOIL_IC_FACTOR} x sqrt(k_supply / k_self-limitation)'},
             'calibrations': {k: {'label': v['label'], 'params': v['params'], 'supply': v['supply']}
                              for k, v in CALIBRATIONS.items()},
             'commensalism_tolerance': MGR_COMM_TOL, 'C_demand_literature': list(C_DEMAND_LIT),
             'outputs': out.written}
    p = out.table(4, 'reproducibility', ext='json')
    with open(p, 'w', encoding='utf-8') as fh:
        json.dump(repro, fh, indent=2, ensure_ascii=False, default=str)
    print(f'Tables and reproducibility record in {out.tab_dir}')
    el = time.perf_counter() - start_time
    print(f'\nTOTAL TIME: {int(el // 3600):02d} h {int(el % 3600 // 60):02d} min {el % 60:04.1f} s')


if __name__ == '__main__':
    main()