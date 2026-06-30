# Script 2: Simulations of Reaction Networks with pyCOT (UPDATED FOR NEW STRUCTURE)

# ========================================
# 1. LIBRARY LOADING AND CONFIGURATION
# ========================================
import os
import sys
import pandas as pd

root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(root_dir, 'src'))  # Add src/ to PYTHONPATH for pyCOT package

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
file_path = 'networks/testing/Eigen_simple.txt'
file_path = 'networks/testing/autopoietic_ext.txt'
file_path = os.path.join(root_dir, file_path)
rn = read_txt(file_path)

additional_laws = {
    'saturated': rate_saturated,
    'threshold_memory': rate_threshold_memory,
    'cosine': rate_cosine
}

# ========================================
# DEFINE SEMANTIC CATEGORIES
# ========================================
# Species: l (source), a (autocatalyst), b (intermediate), c (catalyst), p (parasite)
species_list = ['l', 'a', 'b', 'c', 'p']
category_dict = {
    'base_network': ['l', 'a', 'b','c','p'],
    #'extended parasitic': ['p']
}
semantic_partition = define_semantic_categories(species_list, category_dict)

# ========================================
# CUSTOM PLOTTING FUNCTION
# ========================================
def plot_dynamics_single(time_series, species_list, title="Dynamics", save_path=None):
    """
    Plot time series dynamics for all species in a single plot.
    """
    fig, ax = plt.subplots(1, 1, figsize=(12, 6))

    time = time_series['Time'].values

    # Plot each species
    for species_name in species_list:
        if species_name in time_series.columns:
            ax.plot(time, time_series[species_name], linewidth=2.5, label=species_name, marker='o', markersize=3, markevery=20)

    ax.set_xlabel('Time', fontsize=12)
    ax.set_ylabel('Concentration', fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.legend(fontsize=11, loc='best')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if save_path:
        os.makedirs(os.path.dirname(save_path) if os.path.dirname(save_path) else '.', exist_ok=True)
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {save_path}")

    plt.show()
    return fig, ax

# ========================================
# CONFIGURE KINETICS
# ========================================
# 11 reactions: r1-r8 (base + decay), r9-r10 (parasitic), r11 (rescue)
rate_list = ['mak', 'mak', 'mak', 'mak', 'mak', 'mak', 'mak', 'mak', 'mak', 'mak', 'mak','mak']

# rate_list = [
#     'saturated',           # r1:  SR + R => E + SR (production saturates)
#     'threshold_memory',           # r2:  E + WR => SR (economic strengthening, LOW threshold)
#     'threshold_memory',           # r3:  E + DT => WR (economic re-engagement, HIGH threshold)
#     'threshold_memory',           # r4:  T + WR => SR (trust strengthening, LOW threshold)
#     'threshold_memory',           # r5:  T + DT => WR (trust re-engagement, HIGH threshold)
#     'mak',                 # r6:  SR => WR (natural degradation)
#     'threshold_memory',    # r7:  WR + V => DT (violence-driven detachment)
#     'saturated',           # r8:  SR + 2E => T + SR + 2E (trust generation needs prosperity)
#     'mak',                 # r9:  V + T => (violence-trust annihilation)
#     'mak',                 # r10: 2T => (trust decay)
#     'mak',                 # r11: 2WR => 2WR + V (weak tensions - catalytic)
#     'mak',                 # r12: DT + WR => DT + WR + V (detached-weak conflict)
#     'mak',                 # r13: 2DT => 2DT + V (detached frustration)
#     'mak',                 # r14: 2V => (violence decay)
#     'cosine',              # r15: => R (seasonal resources)
#     'mak',                 # r16: 2R => (resource depletion)
#     'mak'                  # r17: 2E => (economic decay)
# ]

# ==========================================
# SHARED PARAMETER DEFINITIONS
# ==========================================

# Production & Saturation
#Vmax_production = 0.5      # Maximum production rate
#Km_production = 0.5         # Half-saturation for production
# Recovery Thresholds (KEY ASYMMETRY)
# Economic pathway: SAME FOR BOTH CASEE
#E_threshold = 0.1    # LOW - easier to strengthen weak
# Trust pathway: ASYMMETRIC
#T_threshold_weak_strong = 1.5    # LOW - easier with trust
#T_threshold_detached_weak = 1.5  # HIGH - detached need more trust
# Base rates for recovery
#k_recovery = 0.1               # Base recovery rate (same for all pathways)
# Degradation
#k_degradation = 0.1            # Natural SR->WR rate
# Violence-driven detachment
#detachment_threshold = 1.5      # WR*V threshold for detachment
#k_detachment_base=5          # Base detachment rate
# Trust generation
#Vmax_trust = 0.1              # Max trust generation rate
#Km_trust_economy = 0.2         # Needs substantial economy (2E ~ 4)
# Trust-Violence annihilation
#k_trust_destruction = 0.01     # Violence destroys trust
# Trust and violence decay
#k_trust_decay = 0.01           # Natural trust erosion
#k_violence_decay = 0.01        # Violence dissipation
# Violence generation
#k_violence_weak = 0.01         # Weak-weak violence
#k_violence_detached_weak = 0.05  # Detached-weak violence
#k_violence_detached = 0.1     # Detached-detached violence
# Resources
#R_amplitude = 0.5             # Seasonal variation amplitude
#R_frequency = 1           # ~12 month period (2π/12)
#k_resource_depletion = 0.02    # Resource consumption rate
# Economy
#_economic_decay = 0.02        # Economic output decay rate
# Rate constants for 11 reactions
k_production_l = 1   # r1: =>l
k_autocatalysis = 1   # r2: l+a=>a+b
k_synthesis_a = 1     # r3: l+b=>a
k_catalysis = 1       # r4: a+c=>2c
k_decay_l = 1         # r5: l=>
k_decay_a = 0.2        # r6: a=>
k_decay_b = 0.2        # r7: b=>
k_decay_c = 0.2        # r8: c=>
k_parasitism = 2     # r9: a+p=>2p
k_decay_p = 0.2        # r10: p=>
k_rescue = 2          # r11: c+p=>a+p
k_give = 0.5
l00 = 1
a00 = 2
b00 = 3
c00 = 1
p00 = 1

tmaxs = 80
nsteps =160
# spec_vector for MAK kinetics: each reaction gets a [k] list
spec_vector = [
    [k_production_l],   # r1
    [k_autocatalysis],  # r2
    [k_synthesis_a],    # r3
    [k_catalysis],      # r4
    [k_decay_l],        # r5
    [k_decay_a],        # r6
    [k_decay_b],        # r7
    [k_decay_c],        # r8
    [k_parasitism],     # r9
    [k_decay_p],        # r10
    [k_rescue],         # r11  
    [k_give]       # r12
]

# ==========================================
# SPEC_VECTOR CONSTRUCTION
# ==========================================
# spec_vector = [
#     # r1: SR + R => E + SR (saturated)
#     [Vmax_production, Km_production],
#     # r2: E + WR => SR (threshold - EASY)
#     [E_threshold, k_recovery],
#     # r3: E + DT => WR (threshold - HARD)
#     [E_threshold, k_recovery],
#     # r4: T + WR => SR (threshold - EASY)
#     [T_threshold_weak_strong, k_recovery], 
#     # r5: T + DT => WR (threshold - HARD)
#     [T_threshold_detached_weak, k_recovery],
#     # r6: SR => WR (mak)
#     [k_degradation],
#     # r7: WR + V => DT (threshold_memory)
#     [detachment_threshold, k_detachment_base],
#     # r8: SR + 2E => T + SR + 2E (saturated)
#     [Vmax_trust, Km_trust_economy],    
#     # r9: V + T => (mak)
#     [k_trust_destruction],
#     # r10: 2T => (mak)
#     [k_trust_decay],
#     # r11: 2WR => 2WR + V (mak)
#     [k_violence_weak],
#     # r12: DT + WR => DT + WR + V (mak)
#     [k_violence_detached_weak],
#     # r13: 2DT => 2DT + V (mak)
#     [k_violence_detached],
#     # r14: 2V => (mak)
#     [k_violence_decay],
#     # r15: => R (cosine)
#     [R_amplitude, R_frequency],
#     # r16: 2R => (mak)
#     [k_resource_depletion],
#     # r17: 2E => (mak)
#     [k_economic_decay]
# ]

# ========================================
# SCENARIO 1: STABLE REGIME
# ========================================
print("\n" + "=" * 80)
print("SCENARIO 1: Stable Autopoietic Regime")
print("=" * 80)

# Species order: [l, a, b, c, p]
l0_s1 = l00      # source metabolite
a0_s1 = a00      # autocatalyst (self-replicating)
b0_s1 = b00      # intermediate
c0_s1 = 0.0      # catalyst (membrane-like)
p0_s1 = 0.0      # parasite (low initial amount)
x0_s1 = [l0_s1, a0_s1, b0_s1, c0_s1, p0_s1]
title_s1 = f"Scenario 1: Stable Regime (l={l0_s1}, a={a0_s1}, b={b0_s1}, c={c0_s1}, p={p0_s1})"

ts_s1, fv_s1 = simulation(
    rn,
    rate=rate_list,
    spec_vector=spec_vector,
    x0=x0_s1,
    t_span=(0, tmaxs),
    n_steps=nsteps
)

print(f"\nInitial state S1: {x0_s1}")
print(f"Final state S1:\n{ts_s1.tail(1)}")

# ========================================
# SCENARIO 2: PERTURBED WITH PARASITE
# ========================================
print("\n" + "=" * 80)
print("SCENARIO 2: Perturbed with High Parasite")
print("=" * 80)

# Species order: [l, a, b, c, p]
l0_s2 = l00      # source metabolite
a0_s2 = a00      # autocatalyst (self-replicating)
b0_s2 = b00      # intermediate
c0_s2 = 0      # catalyst (membrane-like)
p0_s2 = p00      # parasite (high initial amount)
x0_s2 = [l0_s2, a0_s2, b0_s2, c0_s2, p0_s2]
title_s2 = f"Scenario 2: High Parasite (l={l0_s2}, a={a0_s2}, b={b0_s2}, c={c0_s2}, p={p0_s2})"

ts_s2, fv_s2 = simulation(
    rn,
    rate=rate_list,
    spec_vector=spec_vector,
    x0=x0_s2,
    t_span=(0, tmaxs),
    n_steps=nsteps
)

print(f"\nInitial state S2: {x0_s2}")
print(f"Final state S2:\n{ts_s2.tail(1)}")

# ========================================
# SCENARIO 3: RESCUE WITH CATALYST
# ========================================
print("\n" + "=" * 80)
print("SCENARIO 3: Rescue with High Catalyst")
print("=" * 80)

# Species order: [l, a, b, c, p]
l0_s3 = l00      # source metabolite
a0_s3 = a00      # autocatalyst (self-replicating)
b0_s3 = b00      # intermediate
c0_s3 = c00      # catalyst (high - should rescue from parasite)
p0_s3 = p00      # parasite (high initial amount)
x0_s3 = [l0_s3, a0_s3, b0_s3, c0_s3, p0_s3]
title_s3 = f"Scenario 3: Rescue with Catalyst (l={l0_s3}, a={a0_s3}, b={b0_s3}, c={c0_s3}, p={p0_s3})"

ts_s3, fv_s3 = simulation(
    rn,
    rate=rate_list,
    spec_vector=spec_vector,
    x0=x0_s3,
    t_span=(0, tmaxs),
    n_steps=nsteps
)

print(f"\nInitial state S3: {x0_s3}")
print(f"Final state S3:\n{ts_s3.tail(1)}")

# ========================================
# PLOT ALL THREE SCENARIOS
# ========================================
print("\n" + "=" * 80)
print("PLOTTING RESULTS")
print("=" * 80)

print("\nPlotting Scenario 1...")
plot_dynamics_single(
    ts_s1,
    species_list,
    title=title_s1,
    save_path="visualizations/plot_series_ode/scenario_1_stable_regime.png"
)
plt.close()

print("\nPlotting Scenario 2...")
plot_dynamics_single(
    ts_s2,
    species_list,
    title=title_s2,
    save_path="visualizations/plot_series_ode/scenario_2_high_parasite.png"
)
plt.close()

print("\nPlotting Scenario 3...")
plot_dynamics_single(
    ts_s3,
    species_list,
    title=title_s3,
    save_path="visualizations/plot_series_ode/scenario_3_rescue_catalyst.png"
)
plt.close()

print("\n" + "=" * 80)
print("All plots saved successfully!")
print("=" * 80)
# # PARAMETRIZED SIMULATION
# time_series, flux_vector = simulation(
#     rn, 
# #    x0=x0, 
# #    spec_vector=spec_vector,
# #    rate=rate_list,
#     t_span=(0, 50), 
#     n_steps=200 
# )



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