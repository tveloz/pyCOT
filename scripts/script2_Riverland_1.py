# Script: Reaction Network Simulation with Extended Reactions
# Using three kinetic types: cosine, saturated, and mass action

# ========================================
# 1. LIBRARY LOADING AND CONFIGURATION
# ========================================
import os
import sys
import numpy as np

# Add the root directory to the PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import pyCOT modules
from pyCOT.io.functions import read_txt 
from pyCOT.plot_dynamics import plot_series_ode
from pyCOT.simulations import *
from pyCOT.plot_dynamics import *
from pyCOT.abstractions import abstraction_ordinary

# ========================================
# 2. CREATE EXTENDED REACTION NETWORK
# ======================================== 
network_file = 'Txt/Scenario1_baseline_only_reactions.txt'
rn = read_txt(network_file)

S = rn.stoichiometry_matrix()

# Print the reactions in the network
print("Loaded extended reaction network:")
print("R1: => W_V                    (seasonal water inflow)")
print("R2: W_V + V => V + V_W        (water processing by V)")
print("R3: W_V + V_n => V_n + V_W    (water processing with V needs)")
print("R4: W_V => W_R                (water transfer)")
print("R5: W_R + R => R + R_W        (water processing by R)")
print("R6: W_R + R_n => R_n + R_W    (water processing with R needs)")
print("R7: V_W + V_n => V            (satisfaction of V)")
print("R8: V => V_n                  (V need generation)")
print("R9: R_n + R_W => R            (satisfaction of R)")
print("R10: R => R_n                 (R need generation)")
print("R11: V_W + V_W =>             (V water degradation)")
print("R12: R_W + R_W =>             (R water degradation)")

# ========================================
# 3. SIMULATION PARAMETERS WITH SYMMETRY
# ========================================
# Initial Conditions # Species order: [W_V, V, V_W, V_n, W_R, R, R_W, R_n]
x0 = [0, 5, 0, 0, 0, 10, 0, 0]  # Start with populations but no water

# Define kinetic types for each reaction
rate_list = [
    'cosine',     # R1: => W_V (seasonal inflow)
    'saturated',  # R2: W_V + V => V + V_W 
    'saturated',  # R3: W_V + V_n => V_n + V_W
    'mak',        # R4: W_V => W_R (mass action)
    'saturated',  # R5: W_R + R => R + R_W
    'saturated',  # R6: W_R + R_n => R_n + R_W
    'saturated',  # R7: V_W + V_n => V
    'mak',        # R8: V => V_n (mass action)
    'saturated',  # R9: R_n + R_W => R
    'mak',        # R10: R => R_n (mass action)
    'mak',        # R11: V_W + V_W => (mass action degradation)
    'mak',        # R12: R_W + R_W => (mass action degradation)
    'mak',        # R13: 2W_R => (mass action degradation)
    'mak',        # R14: V_n => (mass action degradation)
    'mak',        # R15: R_n => (mass action degradation
    'mak',        # R16: V => (mass action degradation)
    'mak'         # R17: R => (mass action degradation)
]

# Define parameters with SYMMETRY between R and V reactions
# For testing oscillations, set most rates to zero except water inflow
water_processing_params = [2, 0.5]     # [Vmax, Km] for R2, R3, R5, R6 - SET TO ZERO FOR TESTING
satisfaction_params = [0.4, 1.0]         # [Vmax, Km] for R7, R9 - SET TO ZERO FOR TESTING
need_generation_rate = [0.2]             # [k] for R8, R10 - SET TO ZERO FOR TESTING
water_degradation_rate = [0.4]          # [k] for R11, R12, R13 - SET TO ZERO FOR TESTING
death_rate=[0]
birth_rate=[0]

# Specifying the parameter vector for each reaction
spec_vector = [
    [20, 1.0],               # R1:  [A, w] - amplitude=5.0, frequency=1.0 (clear oscillation)
    water_processing_params, # R2:  [Vmax, Km] - water processing by V
    water_processing_params, # R3:  [Vmax, Km] - need-based water processing by V
    [10],                    # R4:  [k] - water transfer rate - SET TO ZERO FOR TESTING
    water_processing_params, # R5:  [Vmax, Km] - water processing by R (SAME as R2)
    water_processing_params, # R6:  [Vmax, Km] - need-based water processing by R (SAME as R3)
    satisfaction_params,     # R7:  [Vmax, Km] - satisfaction of V
    need_generation_rate,    # R8:  [k] - need generation by V
    satisfaction_params,     # R9:  [Vmax, Km] - satisfaction of R (SAME as R7)
    need_generation_rate,    # R10: [k] - need generation by R (SAME as R8)
    water_degradation_rate,  # R11: [k] - V water degradation
    water_degradation_rate,  # R12: [k] - R water degradation (SAME as R11)
    water_degradation_rate,  # R13: [k] - 2W_R water degradation (SAME as R12)
    death_rate,              # R14: [k] - V_n death (SAME as R12)
    death_rate,              # R16: [k] - V_n death (SAME as R12)
    birth_rate,              # R15: [k] - R_n death (SAME as R12)
    birth_rate               # R17: [k] - R_n death (SAME as R12)
]

# Dictionary of additional kinetic laws
additional_laws = {
    'cosine': rate_cosine,
    'saturated': rate_saturated
}

# ========================================
# 4. RUN SIMULATION
# ========================================  
time_series, flux_vector = simulation(
    rn, 
    rate=rate_list, 
    spec_vector=spec_vector, 
    x0=x0, 
    t_span=(0, 200), 
    n_steps=1000+1,
    additional_laws=additional_laws
)

print("time_series0=\n",time_series)
print("flux_vector0=\n",flux_vector)
plot_series_ode(time_series, filename="Riverland_scenario1_extended_network_time_series.png", show_fig=True)
plot_series_ode(flux_vector, filename="Riverland_scenario1_extended_network_flux_vector.png", show_fig=True)

# ========================================
# 5. PLOTS WITH COGNITIVE DOMAIN INTERVALS
# ========================================
# Combined plot: Concentrations, Flows, Histogram
fig, axes = plt.subplots(1, 3, figsize=(18, 4))
plot_series_with_domain_intervals(time_series, flux_vector, S, title="Serie de Tiempo de Concentraciones", save_figure=False, ax=axes[0])
plot_flux_with_domain_intervals(flux_vector, S, title="Serie de Tiempo de Flujos", save_figure=False, ax=axes[1])
plot_process_types_histogram(flux_vector, S, title="Histograma de Tipos de Procesos", save_figure=False, ax=axes[2])
plt.tight_layout()
plt.savefig("./visualizations/process_classification/combined_plots.png", dpi=150, bbox_inches="tight")
plt.show()

# ========================================
# 6. DISPLAY RESULTS
# ========================================
print("\nSimulation completed successfully!")
print("\nTime Series Shape:", time_series.shape)
print("Species order: W_V, V, V_W, V_n, W_R, R, R_W, R_n")

print("\nFinal concentrations:")
final_concentrations = time_series.iloc[-1]
species_names = ['W_V', 'V', 'V_W', 'V_n', 'W_R', 'R', 'R_W', 'R_n']
for i, species in enumerate(species_names):
    print(f"{species}: {final_concentrations.iloc[i]:.4f}")

# Check oscillation in W_V (first 100 time points)
print(f"\nW_V oscillation check (first 10 values): {time_series.iloc[:10, 0].values}")

# ========================================
# 7. VISUALIZATION OF THE ABSTRACTION 
# ========================================  
abstract_time_series = abstraction_ordinary(time_series, threshold=0.05)

plot_abstraction_size(abstract_time_series)

plot_abstraction_sets(abstract_time_series)

plot_abstraction_graph_movie_html(abstract_time_series, filename="abstraction_graph_movie_html.html", interval=400, title="Abstraction Graph - Time")

# ========================================
# 8. PRINT RESULTS
# ======================================== 
print("\nFinal concentrations:")
final_concentrations = time_series.iloc[-1]
species_names = ['W_V']
for i, species in enumerate(species_names):
    print(f"{species}: {final_concentrations.iloc[i]:.4f}")

# Check oscillation in W_V (first 100 time points)
print(f"\nW_V oscillation check (first 10 values): {time_series.iloc[:3].values}")