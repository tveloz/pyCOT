# ========================================
# AMF – TRADE BALANCE MODEL
# ========================================
# Species notation:
# Pact = Plant_active    Plim = Plant_limited
# Cp   = C_plant         Np   = N_plant        Pp = P_plant
# F    = Fungus          M    = Mycelium
# Cf   = C_fungus        Nf   = N_fungus       Pf = P_fungus
# Ns   = N_soil          Ps   = P_soil
# ========================================

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os
import sys
from matplotlib.colors import TwoSlopeNorm

# Add pyCOT to path - locate relative to this script
_script_dir = os.path.dirname(os.path.abspath(__file__))
_pycot_root = os.path.normpath(os.path.join(_script_dir, '..', '..', '..'))
sys.path.insert(0, os.path.join(_pycot_root, 'src'))

from pyCOT.io.functions import read_txt
from pyCOT.simulations.ode import simulation

VIS_DIR = 'projects/AMF/outputs'
os.makedirs(VIS_DIR, exist_ok=True)

# ========================================
# SCENARIO ORDER (columns left to right)
# ========================================
SCENARIOS = ['I', 'II', 'III', 'IV']
SCENARIO_LABELS = {
    'I':   'I – C-limited mutualism\n(N↓P↓)',
    'II':  'II – Strong mutualism\n(N↑P↓)',
    'III': 'III – Commensalism\n(N↓P↑)',
    'IV':  'IV – Parasitism\n(N↑P↑)',
}
COLORS = {
    'I':   '#4A90D9',
    'II':  '#27AE60',
    'III': '#F39C12',
    'IV':  '#E74C3C',
}

# ========================================
# PLOTTING FUNCTIONS
# ========================================

def plot_amf_dynamics_combined(ts_dict, save_path=None):
    """
    Figure 2: 4 rows x 4 columns.
    Rows    : Biomass | Carbon Pools | Nitrogen | Phosphorus
    Columns : Scenario I | II | III | IV
    """
    groups = {
        "Biomass (Plant & Fungus)": ['Pact', 'Plim', 'F', 'M'],
        "Carbon Pools":             ['Cp', 'Cf'],
        "Nitrogen":                 ['Np', 'Nf', 'Ns'],
        "Phosphorus":               ['Pp', 'Pf', 'Ps'],
    }

    nrows = len(groups)    # 4
    ncols = len(SCENARIOS) # 4

    fig, axes = plt.subplots(nrows, ncols, figsize=(16, 14), sharey=False)
    row_labels = list(groups.keys())

    for col_idx, scen in enumerate(SCENARIOS):
        time_series = ts_dict[scen]
        time        = time_series['Time'].values
        all_species = [c for c in time_series.columns if c != 'Time']

        for row_idx, (group, group_vars) in enumerate(groups.items()):
            ax = axes[row_idx, col_idx]

            for var in [v for v in group_vars if v in all_species]:
                lw = 2.5 if var in ['Pact', 'F', 'Cp'] else 1.8
                ax.plot(time, time_series[var], linewidth=lw, label=var)

            if row_idx == 0:
                ax.set_title(
                    SCENARIO_LABELS[scen],
                    fontweight='bold', fontsize=11,
                    color=COLORS[scen], pad=6
                )

            if col_idx == 0:
                ax.set_ylabel(row_labels[row_idx], fontsize=10, fontweight='bold')
            else:
                ax.set_ylabel("")

            if row_idx == nrows - 1:
                ax.set_xlabel("Time (days)", fontsize=9)
            else:
                ax.set_xlabel("")

            ax.legend(fontsize=8, loc='best')
            ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {save_path}")

    plt.close()


def plot_flux_dynamics_combined(flux_dict, save_path=None):
    """
    Figure 3: 6 rows x 4 columns.
    Rows    : Plant Core & Growth | Mortality | Mycelium Uptake |
              Plant Direct Uptake | Symbiotic Exchange | Inputs
    Columns : Scenario I | II | III | IV
    """
    groups = {
        "Plant Core & Growth": ['r1', 'r2', 'r3', 'r4', 'r11', 'r16'],
        "Mortality":           ['r7', 'r8', 'r12', 'r17'],
        "Mycelium Uptake":     ['r9', 'r10'],
        "Plant Direct Uptake": ['r5', 'r6'],
        "Symbiotic Exchange":  ['r14', 'r15', 'r13'],
        "Inputs":              ['r18', 'r19'],
    }

    nrows = len(groups)    # 6
    ncols = len(SCENARIOS) # 4

    fig, axes = plt.subplots(nrows, ncols, figsize=(16, 20), sharey=False)
    row_labels = list(groups.keys())

    for col_idx, scen in enumerate(SCENARIOS):
        flux_df       = flux_dict[scen]
        time          = flux_df['Time'].values
        all_reactions = [c for c in flux_df.columns if c != 'Time']

        for row_idx, (group, rxns) in enumerate(groups.items()):
            ax = axes[row_idx, col_idx]
            plotted_rxns = [r for r in rxns if r in all_reactions]

            for rxn in plotted_rxns:
                ls    = '--' if rxn in ['r4', 'r11', 'r13', 'r16'] else '-'
                label = 'v' + rxn[1:]
                ax.plot(time, flux_df[rxn], linewidth=1.8, linestyle=ls, label=label)

            if row_idx == 0:
                ax.set_title(
                    SCENARIO_LABELS[scen],
                    fontweight='bold', fontsize=11,
                    color=COLORS[scen], pad=6
                )

            if col_idx == 0:
                ax.set_ylabel(row_labels[row_idx], fontsize=10, fontweight='bold')
            else:
                ax.set_ylabel("")

            if row_idx == nrows - 1:
                ax.set_xlabel("Time (days)", fontsize=9)
            else:
                ax.set_xlabel("")

            if group == "Inputs":
                y_values = flux_df[plotted_rxns].values
                ymax     = max(0.51, 1.01 * np.max(y_values))
                ax.set_ylim(-0.02, ymax)

            ax.legend(fontsize=8, loc='best')
            ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {save_path}")

    plt.close()


# ========================================
# BIOLOGICAL BASE PARAMETERS  (calibrated – Scenario IV)
# ========================================

_SPECIES_ORDER = [
    'Pact', 'Plim', 'Cp', 'Np', 'Pp',
    'F', 'M', 'Cf', 'Nf', 'Pf',
    'Ns', 'Ps'
]


def create_rate_list():
    return [
        'mmk', 'mak', 'mak', 'mak', 'mak', 'mak', 'mak', 'mak',
        'mmk', 'mmk', 'mak', 'mak', 'mak', 'mmk', 'mmk', 'mak',
        'mak', 'mak', 'mak'
    ]

# ------------------------------------------------------------------
# x0_dict: condiciones iniciales calibradas (Escenario IV).
# Los valores de Ns y Ps son sobreescritos por get_scenario_config()
# según el escenario activo.
# ------------------------------------------------------------------
# x0_dict = {
#     'Pact': 5.5,   # Calibrado
#     'Plim': 1.5,   # Calibrado
#     'Cp':   45.0,   # Calibrado
#     'Np':   2.0,   # Calibrado
#     'Pp':   0.4,   # Calibrado
#     'F':    20.0,   # Calibrado
#     'M':    0.5,   # Calibrado
#     'Cf':   40.0,    # Calibrado
#     'Nf':   5.0,   # Calibrado
#     'Pf':   0.5,   # Calibrado
#     'Ns':   30.0,   # Calibrado (sobreescrito en escenarios I, II, III)
#     'Ps':   25     # Calibrado (sobreescrito en escenarios I, II, III)
# }


# def create_base_params():
#     return [
#         [3.5, 10.0],  # r1(mmk):  Photosynthesis      (Vmax, Km)
#         [1e-05],               # r2(mak):  Activation
#         [0.015],                # r3(mak):  Senescence
#         [0.0001],               # r4(mak):  Plant growth
#         [0.0005],               # r5(mak):  Direct N uptake by roots
#         [0.005],               # r6(mak):  Direct P uptake by roots
#         [5e-08],               # r7(mak):  Pact recycling
#         [0.005],               # r8(mak):  Plim recycling
#         [1.5, 30.0],  # r9(mmk):  Mycelial N uptake    (Vmax, Km)
#         [1.2, 10.0],  # r10(mmk): Mycelial P uptake    (Vmax, Km)
#         [5e-07],               # r11(mak): Mycelium growth
#         [0.01],               # r12(mak): Mycelium mortality
#         [0.0007],               # r13(mak): Fungal C uptake
#         [0.8, 20.0],  # r14(mmk): F to plant N transfer (Vmax, Km)
#         [0.5, 20.0],  # r15(mmk): F to plant P transfer (Vmax, Km)
#         [6e-07],               # r16(mak): Fungal growth
#         [0.001],               # r17(mak): Fungal mortality
#         [0.0],                       # r18(mak): N input ← sobreescrito por escenario
#         [0.0],                       # r19(mak): P input ← sobreescrito por escenario
#     ]


# # ========================================
# # SCENARIO CONFIGURATION
# # ========================================

# def get_scenario_config(scenario: str, species_order=None):

#     if species_order is None:
#         species_order = _SPECIES_ORDER

#     params  = create_base_params()
#     IDX_R18 = 17
#     IDX_R19 = 18

#     if scenario == 'I':
#         x0_dict['Ns'] = 10
#         x0_dict['Ps'] = 5
#         params[IDX_R18][0] = 0.001
#         params[IDX_R19][0] = 0.01
#         params_name   = "I – N↓ P↓"
#         scenario_name = "I – C-limited mutualism"

#     elif scenario == 'II':
#         x0_dict['Ns'] = 45
#         x0_dict['Ps'] = 5
#         params[IDX_R18][0] = 0.5
#         params[IDX_R19][0] = 0.1
#         params_name   = "II – N↑ P↓"
#         scenario_name = "II – Strong mutualism"

#     elif scenario == 'III':
#         x0_dict['Ns'] = 10
#         x0_dict['Ps'] = 45
#         params[IDX_R18][0] = 0.001
#         params[IDX_R19][0] = 2.0
#         params_name   = "III – N↓ P↑"
#         scenario_name = "III – Commensalism"

#     elif scenario == 'IV':
#         # Valores óptimos del escenario IV (Calibrados)
#         x0_dict['Ns'] = 90   
#         x0_dict['Ps'] = 90 
#         params[IDX_R18][0] = 0.01  
#         params[IDX_R19][0] = 4.5  
#         params_name   = "IV – N↑ P↑"
#         scenario_name = "IV – Parasitism"

#     else:
#         raise ValueError("Invalid scenario")

#     x0 = [x0_dict[s] for s in species_order]

#     return params, x0, params_name, scenario_name

x0_dict = {
    'Pact': 6.010169486,   # Experimental
    'Plim': 1.537650907,   # Experimental
    'Cp':   40.95843139,   # TR: Plant Trait Database
    'Np':   2.104856174,   # Experimental
    'Pp':   0.419372747,   # Experimental
    'F':    16.77984151,   # Experimental
    'M':    0.557905064,   # Experimental
    'Cf':   29.02, # 37.42861371,   # Zhang et al. 2025
    'Nf':   2.72,  # 5.379405917,   # Zhang et al. 2025
    'Pf':   2.27,  # 0.45137232,    # Zhang et al. 2025
    'Ns':   91.64412942,   # Vidal (sobreescrito en escenarios I, II, III, IV)
    'Ps':   78.40399121    # Vidal (sobreescrito en escenarios I, II, III, IV)
}


def create_base_params():
    return [
        [3.524262147, 11.04293728],  # r1(mmk):  Photosynthesis      (Vmax, Km)
        [1.09521e-05],               # r2(mak):  Activation
        [0.016776613],               # r3(mak):  Senescence
        [0.000103496],               # r4(mak):  Plant growth
        [0.000493452],               # r5(mak):  Direct N uptake by roots
        [0.004409373],               # r6(mak):  Direct P uptake by roots
        [3.95453e-08],               # r7(mak):  Pact recycling
        [0.005492019],               # r8(mak):  Plim recycling
        [1.616387089, 29.45262544],  # r9(mmk):  Mycelial N uptake    (Vmax, Km)
        [1.183213002, 9.616567305],  # r10(mmk): Mycelial P uptake    (Vmax, Km)
        [4.69544e-07],               # r11(mak): Mycelium growth
        [0.008382089],               # r12(mak): Mycelium mortality
        [0.000744158],               # r13(mak): Fungal C uptake
        [0.951736387, 24.35077111],  # r14(mmk): F to plant N transfer (Vmax, Km)
        [0.600544796, 17.59860306],  # r15(mmk): F to plant P transfer (Vmax, Km)
        [6.29787e-07],               # r16(mak): Fungal growth
        [0.001030202],               # r17(mak): Fungal mortality
        [0.0],                       # r18(mak): N input ← sobreescrito por escenario
        [0.0],                       # r19(mak): P input ← sobreescrito por escenario
    ]


# ========================================
# SCENARIO CONFIGURATION
# ========================================

def get_scenario_config(scenario: str, species_order=None):

    if species_order is None:
        species_order = _SPECIES_ORDER

    params  = create_base_params()
    IDX_R18 = 17
    IDX_R19 = 18

    if scenario == 'I':
        x0_dict['Ns'] = 10
        x0_dict['Ps'] = 5
        params[IDX_R18][0] = 0.001
        params[IDX_R19][0] = 0.01
        params_name   = "I – N↓ P↓"
        scenario_name = "I – C-limited mutualism"

    elif scenario == 'II':
        x0_dict['Ns'] = 45
        x0_dict['Ps'] = 5
        params[IDX_R18][0] = 0.5
        params[IDX_R19][0] = 0.1
        params_name   = "II – N↑ P↓"
        scenario_name = "II – Strong mutualism"

    elif scenario == 'III':
        x0_dict['Ns'] = 10
        x0_dict['Ps'] = 45
        params[IDX_R18][0] = 0.001
        params[IDX_R19][0] = 2.0
        params_name   = "III – N↓ P↑"
        scenario_name = "III – Commensalism"

    elif scenario == 'IV':
        # Valores óptimos del escenario IV (Calibrados de la lista)
        x0_dict['Ns'] = 92.64412942   
        x0_dict['Ps'] = 78.40399121 
        params[IDX_R18][0] = 0.01  
        params[IDX_R19][0] = 4.5 
        params_name   = "IV – N↑ P↑"
        scenario_name = "IV – Parasitism"

    else:
        raise ValueError("Invalid scenario")

    x0 = [x0_dict[s] for s in species_order]

    return params, x0, params_name, scenario_name

# ========================================
# INDICATORS
# ========================================

def calculate_amf_indicators(flux_vector, t_span, steady_frac=0.2):
    """Calculate mycorrhizal function indicators at steady state."""
    t_ss_init = t_span[1] - steady_frac * (t_span[1] - t_span[0])
    time_arr  = flux_vector['Time'].values
    n_crit    = np.searchsorted(time_arr, t_ss_init)

    def ss(rxn):
        return float(flux_vector[rxn].iloc[n_crit:].mean())

    v1  = ss('r1')
    v6  = ss('r6')
    v13 = ss('r13')
    v15 = ss('r15')

    C_demand  = v13 / (v1  + 1e-10)
    P_benefit = v15 / (v6  + v15 + 1e-10)

    ind = {'C_supply': v1, 'C_demand': C_demand, 'P_benefit': P_benefit}
    return ind, n_crit


# ========================================
# PLOT TRADE BALANCE MODEL
# ========================================

def plot_trade_balance_model(rn, save_path=None):

    rate_list     = create_rate_list()
    n             = 2
    t_span        = (0, 600)
    n_steps       = 601
    species_order = [s.name for s in rn.species()]
    results       = {}

    ind_colors = {
        'C_supply':  '#2ECC71',
        'C_demand':  '#E67E22',
        'P_benefit': '#3498DB'
    }

    labels = {
        'I':   'I – C-limited mutualism (N↓P↓)',
        'II':  'II – Strong mutualism (N↑P↓)',
        'III': 'III – Commensalism (N↓P↑)',
        'IV':  'IV – Parasitism (N↑P↑)',
    }

    print("=" * 65)
    print("  AMF TRADE BALANCE MODEL (Biological Parameters)")
    print("=" * 65)

    ts_dict   = {}
    flux_dict = {}

    for scen in SCENARIOS:

        params, x0, params_name, scenario_name = get_scenario_config(scen, species_order)
        print(f"\n[{scen}] {scenario_name}")

        time_series, flux_vector = simulation(
            rn, x0=x0, rate=rate_list, spec_vector=params,
            t_span=t_span, n_steps=n_steps,
            method='LSODA', rtol=1e-8, atol=1e-10
        )

        ts_dict[scen]   = time_series
        flux_dict[scen] = flux_vector

        ind, _ = calculate_amf_indicators(flux_vector, t_span)
        results[scen] = ind

    # ── Figure 4: Indicator temporal trajectories (2x2) ─────
    layout_ind  = {'III': (0, 0), 'IV': (0, 1), 'I': (1, 0), 'II': (1, 1)}
    fig1, axes1 = plt.subplots(2, 2, figsize=(14, 10))

    for scen, (row, col) in layout_ind.items():
        ax  = axes1[row][col]
        _, _, _, scenario_name = get_scenario_config(scen, species_order)
        t   = flux_dict[scen]['Time'].values
        v1  = flux_dict[scen]['r1'].values
        v13 = flux_dict[scen]['r13'].values
        v15 = flux_dict[scen]['r15'].values
        v6  = flux_dict[scen]['r6'].values
        pb  = v15 / (v15 + v6 + 1e-10)

        ax.plot(t, v1,  color=ind_colors['C_supply'],  lw=2.2, label='C supply')
        ax.plot(t, v13, color=ind_colors['C_demand'],  lw=2.2, label='C demand')
        ax.plot(t, pb,  color=ind_colors['P_benefit'], lw=2.2, label='P benefit')

        ax.set_title(scenario_name, fontweight='bold', fontsize=12, color=COLORS[scen])
        ax.set_xlabel("Time (days)", fontsize=10)
        ax.set_ylabel("Indicator values", fontsize=10)
        ax.set_ylim(0, max(1, max(np.max(v1), np.max(v13), np.max(pb)) * 1.05))
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.25)

    fig1.text(0.53, 0.02,
              r"Limitation $\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad$ Luxury",
              ha='center', fontsize=14, color='#555')
    fig1.text(0.53, 0.001, r"N supply",
              ha='center', fontsize=18, fontweight='bold', color='#555')
    fig1.text(0.001, 0.53, r"P supply",
              va='center', rotation='vertical', fontsize=18, fontweight='bold', color='#555')
    fig1.text(0.02, 0.52,
              r"Limitation $\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad$ Luxury",
              va='center', rotation='vertical', fontsize=14, color='#555')

    plt.tight_layout(rect=[0.04, 0.04, 1, 1])

    if save_path:
        p1 = save_path + '_indicators_ts.png'
        os.makedirs(os.path.dirname(p1), exist_ok=True)
        fig1.savefig(p1, dpi=300, bbox_inches='tight')
        print(f"Saved: {p1}")

    plt.close(fig1)

    # ── Figure 2: combined 4x4 state variables ───────────────
    fig2_path = os.path.join(
        os.path.dirname(save_path), 'figure2_state_variables.png'
    ) if save_path else None

    plot_amf_dynamics_combined(ts_dict, save_path=fig2_path)

    # ── Figure 3: combined 6x4 reaction fluxes ───────────────
    fig3_path = os.path.join(
        os.path.dirname(save_path), 'figure3_reaction_fluxes.png'
    ) if save_path else None

    plot_flux_dynamics_combined(flux_dict, save_path=fig3_path)

    # ── Figure 4: Phase map 2x2 ──────────────────────────────
    layout      = {'III': (0, 0), 'IV': (0, 1), 'I': (1, 0), 'II': (1, 1)}
    fig4, axes4 = plt.subplots(2, 2, figsize=(14, 10))
    bar_keys    = ['C_supply', 'C_demand', 'P_benefit']
    bar_xlabels = ['C supply', 'C demand', 'P benefit']
    bar_cols    = [ind_colors[k] for k in bar_keys]

    for scen, (row, col) in layout.items():

        ax   = axes4[row][col]
        vals = [results[scen][k] for k in bar_keys]
        bars = ax.bar(
            bar_xlabels, vals,
            color=bar_cols, alpha=0.85,
            edgecolor='white', linewidth=1.2, width=0.55
        )

        for bar, val in zip(bars, vals):
            ypos = val + 0.003 if val >= 0 else val - 0.006
            ax.text(
                bar.get_x() + bar.get_width() / 2, ypos,
                f"{val:.4f}",
                ha='center',
                va='bottom' if val >= 0 else 'top',
                fontsize=9, fontweight='bold'
            )

        ax.axhline(0, color='black', lw=1.0, ls='--', alpha=0.6)
        ax.set_ylabel("Steady State Value", fontsize=9)
        ax.set_title(labels[scen], fontweight='bold', fontsize=11, color=COLORS[scen])
        ax.grid(True, axis='y', alpha=0.25)
        ax.set_ylim(0, max(1, max(vals) * 1.05))

    fig4.text(0.53, 0.02,
              r"Limitation $\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad$ Luxury",
              ha='center', fontsize=14, color='#555')
    fig4.text(0.53, 0.001, r"N supply",
              ha='center', fontsize=18, fontweight='bold', color='#555')
    fig4.text(0.001, 0.53, r"P supply",
              va='center', rotation='vertical', fontsize=18, fontweight='bold', color='#555')
    fig4.text(0.02, 0.52,
              r"Limitation $\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad$ Luxury",
              va='center', rotation='vertical', fontsize=14, color='#555')

    plt.tight_layout(rect=[0.04, 0.04, 1, 1])

    if save_path:
        p4 = save_path + '_phase_map.png'
        fig4.savefig(p4, dpi=300, bbox_inches='tight')
        print(f"Saved: {p4}")

    plt.show()

    return results


# ========================================
# PLOT MYCORRHIZAL FUNCTION
# ========================================

def plot_mycorrhizal_function(rn, save_path=None):

    rate_list     = create_rate_list()
    n             = 2
    t_span        = (0, 600)
    n_steps       = 601
    species_order = [s.name for s in rn.species()]
    steady_frac   = 0.2

    col_v1  = '#2ECC71'
    col_v13 = '#E67E22'
    col_pr  = '#3498DB'

    print("\n" + "=" * 65)
    print("  AMF MYCORRHIZAL FUNCTION – MODULATION BY N AND P")
    print("=" * 65)

    mf_results = {}

    for scen in SCENARIOS:

        params, x0, _, _ = get_scenario_config(scen, species_order)

        _, flux_vector = simulation(
            rn, x0=x0, rate=rate_list, spec_vector=params,
            t_span=t_span, n_steps=n_steps,
            method='LSODA', rtol=1e-8, atol=1e-10
        )

        t_ss = t_span[1] - steady_frac * (t_span[1] - t_span[0])
        n_c  = np.searchsorted(flux_vector['Time'].values, t_ss)
        ss   = lambda r: float(flux_vector[r].iloc[n_c:].mean())

        v1  = ss('r1')
        v13 = ss('r13')
        v15 = ss('r15')
        v6  = ss('r6')

        C_demand = v13 / (v1 + 1e-10)
        pb       = v15 / (v15 + v6 + 1e-10)

        mf_results[scen] = {
            'v1': v1, 'v13': v13, 'C_demand': C_demand,
            'v15': v15, 'v6': v6, 'P_benefit': pb
        }

    # Heatmap
    layout = {'I': (1, 0), 'III': (0, 0), 'II': (1, 1), 'IV': (0, 1)}
    metrics    = ['v1', 'C_demand', 'P_benefit']
    met_labels = [
        'C supply\n(v1 = r1)\nPhotosynthesis',
        'C demand\n(v13/v1 = r13/r1)\nFungal C uptake',
        'P exchange\n(v15/(v15+v6))\nMycorrhizal P benefit'
    ]
    met_colors = [col_v1, col_v13, col_pr]
    matrices   = {m: np.zeros((2, 2)) for m in metrics}

    for scen, (row, col) in layout.items():
        for m in metrics:
            matrices[m][row, col] = mf_results[scen][m]

    global_vmin = min(matrices[m].min() for m in metrics)
    global_vmax = max(matrices[m].max() for m in metrics)
    global_vc   = 0.0 if global_vmin < 0 < global_vmax else (global_vmin + global_vmax) / 2

    fig2, axes2 = plt.subplots(1, 3, figsize=(17, 6))
    scen_in_cell = {v: k for k, v in layout.items()}

    im_last = None

    for ax, m, mlabel, mcol in zip(axes2, metrics, met_labels, met_colors):
        mat = matrices[m]

        try:
            norm   = TwoSlopeNorm(vmin=global_vmin - 1e-6, vcenter=global_vc, vmax=global_vmax + 1e-6)
            cmap   = 'RdYlGn' if m != 'C_demand' else 'RdYlGn_r'
            im     = ax.imshow(mat, cmap=cmap, norm=norm, aspect='auto')
        except Exception:
            im     = ax.imshow(mat, cmap='YlOrRd', aspect='auto')

        im_last = im

        for r in range(2):
            for c in range(2):
                val = mat[r, c]
                ax.text(c, r, f"{val:.4f}",
                        ha='center', va='center', fontsize=13, fontweight='bold',
                        color='white' if abs(val) > 0.4 * max(abs(global_vmin), abs(global_vmax)) else 'black')
                ax.text(c, r - 0.38, f"Scen. {scen_in_cell.get((r,c),'?')}",
                        ha='center', va='center', fontsize=8, color='white', style='italic')

        ax.set_xticks([0, 1])
        ax.set_xticklabels(['N limitation', 'N luxury'])
        ax.set_yticks([0, 1])
        ax.set_yticklabels(['P luxury', 'P limitation'])
        ax.set_title(mlabel, fontweight='bold', fontsize=11, color=mcol, pad=10)

        for (row, col), scen in scen_in_cell.items():
            ax.add_patch(plt.Rectangle(
                (col - 0.5, row - 0.5), 1, 1,
                linewidth=3, edgecolor=COLORS[scen], facecolor='none', zorder=5
            ))

    fig2.subplots_adjust(right=0.88)
    cbar_ax = fig2.add_axes([0.90, 0.15, 0.02, 0.7])
    fig2.colorbar(im_last, cax=cbar_ax, label='Indicator value')

    titles = {
        'I':   'I – C-lim. (N↓P↓)',
        'II':  'II – Mutualism (N↑P↓)',
        'III': 'III – Commens. (N↓P↑)',
        'IV':  'IV – Parasitism (N↑P↑)',
    }

    legend_patches = [
        mpatches.Patch(color=COLORS[e], label=titles[e])
        for e in SCENARIOS
    ]

    fig2.legend(
        handles=legend_patches, loc='lower center', ncol=4,
        fontsize=9, bbox_to_anchor=(0.5, -0.06),
        title='Scenarios', title_fontsize=9
    )

    fig2.suptitle(
        "AMF Mycorrhizal Function – N x P Phase Map\n"
        "C supply, C demand, and P exchange ratio (steady state)",
        fontsize=13, fontweight='bold'
    )

    if save_path:
        p2 = save_path + '_mycorrhizal_function_heatmap.png'
        fig2.savefig(p2, dpi=300, bbox_inches='tight')
        print(f"Saved: {p2}")

    plt.close(fig2)

    # Summary table and validation (Johnson 2009)
    print("\n" + "=" * 80)
    print("  EXPECTED RANGE EVALUATION (Johnson 2009, Fig. 1)")
    print("=" * 80)

    ranges = {
        'I':   {'v1': [0.10, 0.40], 'ratio': [0.05, 0.20], 'pb': [0.70, 1]},
        'II':  {'v1': [0.70, 2.50], 'ratio': [0.20, 0.70], 'pb': [0.70, 1]},
        'III': {'v1': [0.10, 0.40], 'ratio': [0.05, 0.20], 'pb': [0.05, 0.50]},
        'IV':  {'v1': [0.20, 0.30], 'ratio': [0.30, 0.50], 'pb': [0.05, 0.20]}
    }

    scen_titles = {
        'I':   'I - C-limited mutualism (N↓ P↓)',
        'II':  'II - Strong mutualism (N↑ P↓)',
        'III': 'III - Commensalism (N↓ P↑)',
        'IV':  'IV - Parasitism (N↑ P↑)',
    }

    def evaluate(val, r_min, r_max, is_ratio=False):
        if r_min <= val <= r_max:
            return "correct"
        elif val < r_min:
            if val < 0.01 and not is_ratio:
                return "near zero"
            if is_ratio:
                return "low ratio"
            return "low"
        else:
            if val > 0.95 and r_max <= 0.95:
                return "spurious/high"
            return "too high"

    for scen in SCENARIOS:

        print(f"\n{scen_titles[scen]}")
        print(f"{'Indicator':<25} | {'Expected':<15} | {'Obtained':<10} | {'Status'}")
        print("-" * 75)

        r  = mf_results[scen]
        v1 = r['v1']
        cd = r['C_demand']
        pb = r['P_benefit']

        rng_v1  = ranges[scen]['v1']
        rng_rat = ranges[scen]['ratio']
        rng_pb  = ranges[scen]['pb']

        print(f"{'C supply (v1)':<25} | {rng_v1[0]:.2f} - {rng_v1[1]:.2f}     | {v1:<10.4f} | {evaluate(v1, rng_v1[0], rng_v1[1])}")
        print(f"{'C demand (v13/v1)':<25} | {rng_rat[0]:.2f} - {rng_rat[1]:.2f}     | {cd:<10.4f} | {evaluate(cd, rng_rat[0], rng_rat[1])}")
        print(f"{'P benefit (v15/(v15+v6))':<25} | {rng_pb[0]:.2f} - {rng_pb[1]:.2f}     | {pb:<10.4f} | {evaluate(pb, rng_pb[0], rng_pb[1])}")

    print("\n" + "=" * 80)

    return mf_results


# ========================================
# MAIN
# ========================================

if __name__ == '__main__':

    FILE_PATH = os.path.join(_pycot_root, 'data', 'Ecological_models', 'AMF.txt')
    rn = read_txt(FILE_PATH)

    print(f"Species:    {[s.name for s in rn.species()]}")
    print(f"Reactions:  {[r.name() for r in rn.reactions()]}")

    save_base = os.path.join(VIS_DIR, '3_trade_balance', 'trade_balance_v9')
    os.makedirs(os.path.dirname(save_base), exist_ok=True)

    results    = plot_trade_balance_model(rn, save_path=save_base)
    mf_results = plot_mycorrhizal_function(rn, save_path=save_base)

    print("\nSimulation completed successfully.")