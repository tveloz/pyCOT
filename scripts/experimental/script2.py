# Script 2: Simulations of Reaction Networks with pyCOT (UPDATED FOR NEW STRUCTURE)

# ========================================
# 1. LIBRARY LOADING AND CONFIGURATION
# ========================================
import os
import sys
import pandas as pd

# Add the root directory to the PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import pyCOT modules - UPDATED IMPORTS
from pyCOT.io.functions import read_txt
from pyCOT.visualization.plot_dynamics import plot_series_ode
from pyCOT.simulations.ode import simulation
from pyCOT.kinetics.deterministic_advanced import * # Updated: explicit import
from pyCOT.kinetics.deterministic_basic import *
from pyCOT.core.semantic_partition import define_semantic_categories
import matplotlib.pyplot as plt
# ========================================
# 2. CREATING THE REACTION_NETWORK OBJECT
# ========================================
# Alternative examples:
# file_path = 'Txt/BZ_cycle.txt'  
# file_path = 'networks/testing/Lotka_Volterra.txt'  
# file_path = 'networks/Riverland_model/Scenario1_baseline_only_reactions.txt'
# file_path = 'Txt/2019fig1.txt'
# file_path = 'Txt/2019fig2.txt'
# file_path = 'Txt/non_connected_example.txt' 
# file_path = 'Txt/PassiveUncomforableIndignated_problemsolution.txt'
# file_path = 'Txt/Farm.txt' 
# file_path = 'Txt/SEIR.txt' 
# file_path = 'Txt/2010Veloz_Ex_4.txt'
print(os.path.dirname(os.path.abspath(__file__)))

file_path = 'networks/Conflict_Theory/Resource_Scarcity_Toy_Model2.txt'

rn = read_txt(file_path)

additional_laws = {
    'saturated': rate_saturated,
    'threshold_memory': rate_threshold_memory,
    'cosine': rate_cosine
}

# ========================================
# DEFINE SEMANTIC CATEGORIES
# ========================================
species_list = ['SR', 'R', 'E', 'WR', 'DT', 'T', 'V']
category_dict = {
    'peace': ['SR', 'R', 'E', 'T'],
    'conflict': ['DT', 'V', 'WR']
}
semantic_partition = define_semantic_categories(species_list, category_dict)

# ========================================
# CUSTOM PLOTTING FUNCTION
# ========================================
def plot_dynamics_separated(time_series, semantic_partition, species_list, title="Dynamics", save_path=None):
    """
    Plot time series dynamics with peace and conflict categories separated into two subplots
    Both subplots share the same y-axis scale for easy comparison
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Get species indices for each category
    peace_species = semantic_partition.category_indices['peace']
    conflict_species = semantic_partition.category_indices['conflict']

    time = time_series['Time'].values

    # Calculate global min/max across all species for shared scale
    all_species_indices = peace_species + conflict_species
    all_values = []
    for idx in all_species_indices:
        species_name = species_list[idx]
        all_values.extend(time_series[species_name].values)

    y_min = min(all_values)
    y_max = max(all_values)
    y_margin = (y_max - y_min) * 0.05  # 5% margin
    y_lim = (y_min - y_margin, y_max + y_margin)

    # Peace species
    ax = axes[0]
    for idx in peace_species:
        species_name = species_list[idx]
        ax.plot(time, time_series[species_name], linewidth=2, label=species_name)

    ax.set_xlabel('Time', fontsize=12)
    ax.set_ylabel('Concentration', fontsize=12)
    ax.set_title('Peace Variables (SR, R, E, T)', fontsize=13, fontweight='bold', color='green')
    ax.set_ylim(y_lim)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    # Conflict species
    ax = axes[1]
    for idx in conflict_species:
        species_name = species_list[idx]
        ax.plot(time, time_series[species_name], linewidth=2, label=species_name)

    ax.set_xlabel('Time', fontsize=12)
    ax.set_ylabel('Concentration', fontsize=12)
    ax.set_title('Conflict Variables (WR, DT, V)', fontsize=13, fontweight='bold', color='red')
    ax.set_ylim(y_lim)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    # Overall title
    plt.suptitle(title, fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()

    if save_path:
        os.makedirs(os.path.dirname(save_path) if os.path.dirname(save_path) else '.', exist_ok=True)
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {save_path}")

    plt.show()
    return fig, axes

# ========================================
# CONFIGURE KINETICS
# ========================================

rate_list = [
    'saturated',           # r1:  SR + R => E + SR (production saturates)
    'threshold_memory',           # r2:  E + WR => SR (economic strengthening, LOW threshold)
    'threshold_memory',           # r3:  E + DT => WR (economic re-engagement, HIGH threshold)
    'threshold_memory',           # r4:  T + WR => SR (trust strengthening, LOW threshold)
    'threshold_memory',           # r5:  T + DT => WR (trust re-engagement, HIGH threshold)
    'mak',                 # r6:  SR => WR (natural degradation)
    'threshold_memory',    # r7:  WR + V => DT (violence-driven detachment)
    'saturated',           # r8:  SR + 2E => T + SR + 2E (trust generation needs prosperity)
    'mak',                 # r9:  V + T => (violence-trust annihilation)
    'mak',                 # r10: 2T => (trust decay)
    'mak',                 # r11: 2WR => 2WR + V (weak tensions - catalytic)
    'mak',                 # r12: DT + WR => DT + WR + V (detached-weak conflict)
    'mak',                 # r13: 2DT => 2DT + V (detached frustration)
    'mak',                 # r14: 2V => (violence decay)
    'cosine',              # r15: => R (seasonal resources)
    'mak',                 # r16: 2R => (resource depletion)
    'mak'                  # r17: 2E => (economic decay)
]

# ==========================================
# SHARED PARAMETER DEFINITIONS
# ==========================================

# Production & Saturation
Vmax_production = 0.5      # Maximum production rate
Km_production = 0.5         # Half-saturation for production
# Recovery Thresholds (KEY ASYMMETRY)
# Economic pathway: SAME FOR BOTH CASEE
E_threshold = 0.1    # LOW - easier to strengthen weak
# Trust pathway: ASYMMETRIC
T_threshold_weak_strong = 1.5    # LOW - easier with trust
T_threshold_detached_weak = 1.5  # HIGH - detached need more trust
# Base rates for recovery
k_recovery = 0.1               # Base recovery rate (same for all pathways)
# Degradation
k_degradation = 0.1            # Natural SR->WR rate
# Violence-driven detachment
detachment_threshold = 1.5      # WR*V threshold for detachment
k_detachment_base=5          # Base detachment rate
# Trust generation
Vmax_trust = 0.1              # Max trust generation rate
Km_trust_economy = 0.2         # Needs substantial economy (2E ~ 4)
# Trust-Violence annihilation
k_trust_destruction = 0.01     # Violence destroys trust
# Trust and violence decay
k_trust_decay = 0.01           # Natural trust erosion
k_violence_decay = 0.01        # Violence dissipation
# Violence generation
k_violence_weak = 0.01         # Weak-weak violence
k_violence_detached_weak = 0.05  # Detached-weak violence
k_violence_detached = 0.1     # Detached-detached violence
# Resources
R_amplitude = 0.5             # Seasonal variation amplitude
R_frequency = 1           # ~12 month period (2π/12)
k_resource_depletion = 0.02    # Resource consumption rate
# Economy
k_economic_decay = 0.02        # Economic output decay rate


# ==========================================
# SPEC_VECTOR CONSTRUCTION
# ==========================================
spec_vector = [
    # r1: SR + R => E + SR (saturated)
    [Vmax_production, Km_production],
    # r2: E + WR => SR (threshold - EASY)
    [E_threshold, k_recovery],
    # r3: E + DT => WR (threshold - HARD)
    [E_threshold, k_recovery],
    # r4: T + WR => SR (threshold - EASY)
    [T_threshold_weak_strong, k_recovery], 
    # r5: T + DT => WR (threshold - HARD)
    [T_threshold_detached_weak, k_recovery],
    # r6: SR => WR (mak)
    [k_degradation],
    # r7: WR + V => DT (threshold_memory)
    [detachment_threshold, k_detachment_base],
    # r8: SR + 2E => T + SR + 2E (saturated)
    [Vmax_trust, Km_trust_economy],    
    # r9: V + T => (mak)
    [k_trust_destruction],
    # r10: 2T => (mak)
    [k_trust_decay],
    # r11: 2WR => 2WR + V (mak)
    [k_violence_weak],
    # r12: DT + WR => DT + WR + V (mak)
    [k_violence_detached_weak],
    # r13: 2DT => 2DT + V (mak)
    [k_violence_detached],
    # r14: 2V => (mak)
    [k_violence_decay],
    # r15: => R (cosine)
    [R_amplitude, R_frequency],
    # r16: 2R => (mak)
    [k_resource_depletion],
    # r17: 2E => (mak)
    [k_economic_decay]
]

# ========================================
# SCENARIO 1: STABLE REGIME
# ========================================
print("\n" + "=" * 80)
print("SCENARIO 1: Stable Productive Regime")
print("=" * 80)

# SCENARIO 1: STABLE PRODUCTIVE REGIME
#[SR,R,E,WR,DT,T,V]
SR0=1.0
WR0=0.5
DT0=1.5
R0=0
E0=0
T0=0
V0=0
x0_stable = [SR0, R0, E0, WR0, DT0, T0, V0]
print(f"  Total Population        = {x0_stable[0] + x0_stable[3] + x0_stable[4]}")

ts_stable, fv_stable = simulation(
    rn,
    rate=rate_list,
    spec_vector=spec_vector,
    x0=x0_stable,
    t_span=(0, 200),
    n_steps=400,
    additional_laws=additional_laws
)


# # PARAMETRIZED SIMULATION
# time_series, flux_vector = simulation(
#     rn, 
# #    x0=x0, 
# #    spec_vector=spec_vector,
# #    rate=rate_list,
#     t_span=(0, 50), 
#     n_steps=200 
# )

# Plot with separated categories
print("\nPlotting dynamics with separated peace/conflict categories...")
plot_dynamics_separated(
    ts_stable,
    semantic_partition,
    species_list,
    title="Stable Regime Dynamics",
    save_path="visualizations/plot_series_ode/stable_regime_separated.png"
)
print("\nFinal state:")
print(ts_stable.tail(1))

# # Extract last state and continue simulation with intervention
# last_state = time_series[['G', 'R', 'V', 'N', 'P', 'F']].iloc[-1].values.tolist()
# x0 = last_state.copy()
# x0[4] = 1.1  # Increase peacekeeping (P)
# spec_vector[8][0] = 0.5  # Adjust peacekeeping rate

# print(f'New x0 after intervention 1: {x0}')

# time_series2, flux_vector2 = simulation(
#     rn, 
#     x0=x0, 
#     t_span=(50, 100), 
#     n_steps=200)

# # Second intervention: add funding
# last_state = time_series2[['G', 'R', 'V', 'N', 'P', 'F']].iloc[-1].values.tolist()
# x0 = last_state.copy()
# spec_vector[9][0] = 0.2  # Activate funding
# x0[5] = 1  # Add funding (F)

# time_series3, flux_vector3 = simulation(
#     rn, 
#     x0=x0, 
#     spec_vector=spec_vector,
#     rate=rate_list,
#     t_span=(100, 150), 
#     n_steps=200
# )

# # Combine time series and plot
# combined_df = pd.concat([time_series, time_series2, time_series3], ignore_index=True)

# color_mapping = {
#     'V': 'red',      # Violence
#     'G': 'blue',     # Grievances
#     'R': 'green',    # Resources
#     'N': 'orange',   # Narratives
#     'P': 'purple',   # Peacekeeping
#     'F': 'brown'     # Funding
# }

# fig, ax = plot_series_ode(combined_df, color_dict=color_mapping)

# ##################################################################################
# # Example 2: ODE simulation with specific parameters (Lotka-Volterra)
# ##################################################################################
# file_path = 'networks/testing/Lotka_Volterra.txt'
# rn = read_txt(file_path)
#
# x0 = [2, 3]  # Initial concentrations [Prey, Predator]
# rate_list = 'mak'
# spec_vector = [[0.5], [0.8], [0.1]]  # [birth, predation, death]
#
# time_series, flux_vector = simulation(
#     rn, 
#     rate=rate_list, 
#     spec_vector=spec_vector, 
#     x0=x0, 
#     t_span=(0, 100), 
#     n_steps=500
# )
#
# print("Lotka-Volterra Time Series:")
# print(time_series)
# plot_series_ode(time_series)

# ##################################################################################
# # Example 3: Mixed kinetics (MAK + MMK)
# ##################################################################################
# x0 = [0, 1, 0]
# rate_list = ['mak', 'mak', 'mmk', 'mmk', 'mak']
# spec_vector = [
#     [0.7],           # MAK: k
#     [0.5],           # MAK: k
#     [1.0, 0.3],      # MMK: [Vmax, Km]
#     [1.0, 0.4],      # MMK: [Vmax, Km]
#     [1.0]            # MAK: k
# ]
#
# time_series, flux_vector = simulation(
#     rn, 
#     rate=rate_list, 
#     spec_vector=spec_vector, 
#     x0=x0, 
#     t_span=(0, 50), 
#     n_steps=500
# )
#
# print("Mixed Kinetics Time Series:")
# print(time_series)
# plot_series_ode(time_series)

# ############################################################################
# # DEFINING CUSTOM KINETICS (UPDATED APPROACH)
# ############################################################################
# from pyCOT.kinetics.deterministic_basic import rate_mak
#
# # Custom kinetic: Ping-pong mechanism
# def rate_ping_pong(substrates, concentrations, species_idx, spec_vector):
#     """
#     Ping-pong enzyme mechanism.
#     Rate = Vmax * A * B / (KmA * B + KmB * A + A * B)
#     """
#     Vmax, KmA, KmB = spec_vector
#     if len(substrates) < 2:
#         return 0
#     
#     substrateA = substrates[0][0]
#     substrateB = substrates[1][0]
#     A = concentrations[species_idx[substrateA]]
#     B = concentrations[species_idx[substrateB]]
#     
#     return Vmax * A * B / (KmA * B + KmB * A + A * B)
#
# rate_ping_pong.expression = lambda substrates, reaction: (
#     "0 (ping-pong requires two substrates)"
#     if len(substrates) < 2 else
#     f"(Vmax_{reaction} * [{substrates[0][0]}] * [{substrates[1][0]}]) / "
#     f"(Km_{substrates[0][0]} * [{substrates[1][0]}] + "
#     f"Km_{substrates[1][0]} * [{substrates[0][0]}] + "
#     f"[{substrates[0][0]}] * [{substrates[1][0]}])"
# )
#
# # Custom kinetic: Threshold-based rate
# def rate_thresholds(reactants, concentrations, species_idx, spec_vector):
#     """
#     Rate with min/max thresholds.
#     """
#     k, threshold_min, threshold_max = spec_vector
#     
#     # Compute base MAK rate
#     rate_mak_value = rate_mak(reactants, concentrations, species_idx, [k])
#     
#     # Apply thresholds
#     if rate_mak_value < threshold_min:
#         return 0.0
#     elif rate_mak_value > threshold_max:
#         return threshold_max
#     else:
#         return rate_mak_value
#
# rate_thresholds.expression = lambda reactants, reaction: (
#     f"max(0, min(threshold_max_{reaction}, k_{reaction} * " + 
#     " * ".join(f"[{r}]^{c}" if c != 1 else f"[{r}]" for r, c in reactants) + 
#     "))"
# )
#
# # ##################################################################################
# # # Example 4: Simulation with custom kinetics
# # ##################################################################################
# x0 = [80, 50, 60]
# rate_list = ['mmk', 'hill', 'mak', 'ping_pong', 'threshold']
# spec_vector = [
#     [1.0, 0.3],      # MMK: [Vmax, Km]
#     [1.0, 2, 2],     # Hill: [Vmax, K, n]
#     [0.7],           # MAK: [k]
#     [1.0, 0.4, 0.6], # Ping-pong: [Vmax, KmA, KmB]
#     [0.1, 0.6, 8]    # Threshold: [k, min, max]
# ]
#
# # Register custom kinetics
# additional_laws = {
#     'ping_pong': rate_ping_pong, 
#     'threshold': rate_thresholds
# }
#
# time_series, flux_vector = simulation(
#     rn, 
#     rate=rate_list, 
#     spec_vector=spec_vector, 
#     x0=x0, 
#     t_span=(0, 50), 
#     n_steps=100, 
#     additional_laws=additional_laws
# )
#
# print("Custom Kinetics Time Series:")
# print(time_series)
# plot_series_ode(time_series)