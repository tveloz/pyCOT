# Script 2: Simulations of Reaction Networks with pyCOT

# ========================================
# 1. LIBRARY LOADING AND CONFIGURATION
# ========================================
# Import necessary libraries and modules
import os
import sys

# Add the root directory to the PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import pyCOT modules
from pyCOT.io.functions import read_txt 
from pyCOT.plot_dynamics import *
from pyCOT.simulations import *

# ========================================
# 2. CREATING THE REACTION_NETWORK OBJECT
# ========================================
# Alternative examples:
file_path = 'Txt/BZ_cycle.txt'  
file_path = 'networks/testing/Lotka_Volterra.txt'  

#file_path = 'networks/Riverland_model/Scenario1_baseline_only_reactions.txt'  # Autopoietic network
# file_path = 'Txt/2019fig1.txt'
# file_path = 'Txt/2019fig2.txt'
# file_path = 'Txt/non_connected_example.txt' 
# file_path = 'Txt/PassiveUncomforableIndignated_problemsolution.txt'
# file_path = 'Txt/Farm.txt' 
# file_path = 'Txt/SEIR.txt' 
# file_path = 'Txt/2010Veloz_Ex_4.txt'
file_path = 'networks/Conflict_Theory/conflict_toy_model0.txt'

rn = read_txt(file_path)  # Create ReactionNetwork object from file

# ========================================
# 3. SIMULATION OF THE REACTION NETWORK
# ========================================
# ###################################################################################
# # Example 1: ODE simulation with one kinetic equation for all reactions
# ###################################################################################

xg=2
xv=0.5
xr=5
xn=1
xp=0
xf=0
x0 = [xg,xr,xv,xn,xp,xf]  # Initial concentrations
#rate_list = 'mak'  # Kinetics for all reactions
spec_vector = [
    [0.2],#,r1:  G + R => V + R #{eq:mobilization}\\
    [1.0], #r2:  V => 2G #{eq:grievance_generation}\\
    [0.5],#r3:  V=> N + V#{eq:narrative_production}\\
    [0.1],#r4:  N + R => V + N #{eq:narrative_mobilization}\\
    [0.1],#r5:  G + N => R + N #{eq:resource_recruitment}
    [0.05], #G=>
    [0.1], #N=>
    [0.04],#V=>
    [0.5], #P+V=>P
    [0]] # F=>F+R
# Parameters for each reaction
rate_list = "mak"  # Kinetics for all reactions ("mak": mass action kinetics, "mmk": Michaelis-Menten kinetics, "hill": Hill kinetics)
#PARAMETRIZED SIMULATION

time_series, flux_vector = simulation(rn, x0=x0, spec_vector=spec_vector,rate=rate_list,t_span=(0, 50), n_steps=200)
#NON-PARAMETRIZED SIMULATION
#time_series, flux_vector = simulation(rn, rate=rate_list)

print("ODE Time Series:")

plot_series_ode(time_series)
print(time_series[-1:])


# Extract just the values you need
last_state = time_series[['G', 'R', 'V', 'N', 'P', 'F']].iloc[-1].values.tolist()

# Create new x0 with modifications
x0 = last_state.copy()  # Make a copy to avoid modifying the original

# Now modify F and P (assuming order is [G, R, V, N, P, F])
# Based on your column order in the selection
x0[4] = 1.1  # P is 4th in ['G','R','V','N','P','F']
spec_vector[8][0]=0.5
print(f'new x0 = {x0}')

time_series2, flux_vector2 = simulation(rn, x0=x0, spec_vector=spec_vector,rate=rate_list,t_span=(50, 100), n_steps=200)
# 1 — create new timescale (for example window 5)
 # Extract just the values you need
last_state = time_series2[['G', 'R', 'V', 'N', 'P', 'F']].iloc[-1].values.tolist()

# Create new x0 with modifications
x0 = last_state.copy()  # Make a copy to avoid modifying the original
spec_vector[9][0]=0.2
x0[5] = 1  # F is 5th in ['G','R','V','N','P','F']

time_series3, flux_vector3 = simulation(rn, x0=x0, spec_vector=spec_vector,rate=rate_list,t_span=(100, 150), n_steps=200)
combined_df = pd.concat([time_series, time_series2,time_series3], ignore_index=True)


color_mapping = {'V': 'red', 'G': 'blue', 'R': 'green', 'N': 'orange', 'P': 'purple', 'F': 'brown'}
fig, ax = plot_series_ode(combined_df, color_dict=color_mapping)

#process_rescaled = generate_rescaled_process_series(flux_vector, window_size=3)
# 2 — classify rescaled process vectors
#classifications = classify_process_series(process_rescaled, rn.stoichiometry_matrix())

# 3 — plot classification grid panel
#plot_process_classes_over_time(process_rescaled["Time"], classifications)
# ###################################################################################
# # Example 2: ODE simulation with one kinetic equation for all reactions and specific parameters
# ###################################################################################
# x0 = [2, 3]  # Initial concentrations
# rate_list = 'mak'  # Kinetics for all reactions
# spec_vector = [[0.5], [0.8], [0.1]]  # Parameters for each reaction

# # # Run simulation
# time_series, flux_vector = simulation(
#     rn, rate=rate_list, spec_vector=spec_vector, 
#     x0=x0, t_span=(0, 100), n_steps=500
# )

# print("ODE Time Series:")
#print(time_series)
#plot_series_ode(time_series, filename="time_series_plot.png")
# ##################################################################################
# # Example 3: ODE simulation using a mix of 'mak' and 'mmk' kinetics
# ##################################################################################
# x0 = [0, 1, 0]  # Initial concentrations
# rate_list = ['mak', 'mak', 'mmk', 'mmk', 'mak']
# spec_vector = [[0.7], [0.5], [1.0, 0.3], [1.0, 0.4], [1.0]]

# time_series, flux_vector = simulation(
#     rn, rate=rate_list, spec_vector=spec_vector, 
#     x0=x0, t_span=(0, 50), n_steps=500
# )

# print("ODE Time Series:")
# print(time_series)
# plot_series_ode(time_series)

# ############################################################################
# # DEFINING NEW KINETICS
# ############################################################################
# ###############################
# # Función nueva: rate_ping_pong
# def rate_ping_pong(substrates, concentrations, species_idx, spec_vector):
#     """
#     Calculate the rate of reaction using the ping-pong mechanism. The ping-pong mechanism is a type of enzyme kinetics where the enzyme alternates between two states. It assumes two substrates and a single product.
#     Parameters:
#     - substrates: List of tuples (species, stoichiometric coefficient) for the substrates. Example: substrates=[('A', 1), ('B', 1)]
#     - concentrations: List of current concentrations of the species. Example: concentrations = [0.1, 0.2, 0.3]
#     - species_idx: Dictionary mapping species names to their indices in the concentrations list. Example: {'A': 0, 'B': 1, 'C': 2}
#     - spec_vector: List of parameters for the reaction, including. Example: spec_vector = [Vmax, KmA, KmB]
#     Returns:
#     - rate: The reaction rate equation for the ping-pong mechanism. Example: rate = Vmax * A * B / (KmA * B + KmB * A + A * B)
#     """
#     Vmax, KmA, KmB = spec_vector # Extract the parameters from the spec_vector
#     if len(substrates) < 2:      # Check if there are at least two substrates
#         return 0  # Return 0 if not enough substrates are provided
    
#     substrateA = substrates[0][0]               # Assuming the first substrate is A. 
#     substrateB = substrates[1][0]               # Assuming the second substrate is B
#     A = concentrations[species_idx[substrateA]] # Concentration of substrate A
#     B = concentrations[species_idx[substrateB]] # Concentration of substrate B
#     return Vmax * A * B / (KmA * B + KmB * A + A * B) # The reaction rate equation for ping-pong mechanism
# rate_ping_pong.expression = lambda substrates, reaction: (
#     "0 (ping-pong requires two substrates)"
#     if len(substrates) < 2 else
#     f"(Vmax_{reaction} * [{substrates[0][0]}] * [{substrates[1][0]}]) / (Km_{substrates[0][0]} * [{substrates[1][0]}] + Km_{substrates[1][0]} * [{substrates[0][0]}] + [{substrates[0][0]}] * [{substrates[1][0]}])"
# ) # Expression of the equation for the ping-pong rate. 

# ################################
# # Función nueva: rate_thresholds 
# def rate_thresholds(reactants, concentrations, species_idx, spec_vector):
#     k, threshold_min, threshold_max= spec_vector 
    
#     rate_mak_value = rate_mak(reactants, concentrations, species_idx, [k]) 
    
#     if rate_mak_value < threshold_min:
#         return 0.0  
#     elif rate_mak_value > threshold_max:
#         return threshold_max
#     else:
#         return rate_mak_value
# rate_thresholds.expression = lambda reactants, reaction: (
#     f"0 if (rate < threshold_min_{reaction}); "
#     f"v{reaction[1]} = threshold_max_{reaction} if (rate > threshold_max_{reaction}); "
#     f"v{reaction[1]} = k_{reaction} * " + " * ".join(
#         f"[{reactant}]^{coef}" if coef != 1 else f"[{reactant}]"
#         for reactant, coef in reactants
#     ) + f" if (threshold_min_{reaction} <= rate <= threshold_max_{reaction})"
# )
# # ##################################################################################
# # # Example 4: Simulation with defined parameters for additional kinetics
# # ##################################################################################
# x0 = [80, 50, 60]
# rate_list = ['mmk', 'hill', 'mak', 'ping_pong', 'threshold']
# spec_vector = [
#     [1.0, 0.3],      # Parameters for 'mmk' (Vmax, Km)
#     [1.0, 2, 2],     # Parameters for 'hill' (Vmax, n, K)
#     [0.7],           # Parameter  for 'mak' (k)
#     [1.0, 0.4, 0.6], # Parameters for 'ping_pong' (Vmax, KmA, KmB)
#     [0.1, 0.6, 8]  # Parameters for 'threshold' (k, threshold_min, threshold_max)
# ]

# additional_laws = {'ping_pong': rate_ping_pong, 'threshold': rate_thresholds}

# time_series, flux_vector = simulation(rn, rate=rate_list, spec_vector=spec_vector, x0=x0, t_span=(0, 50), n_steps=100, additional_laws=additional_laws)


# print("ODE Time Series with Defined Parameters:")
# print(time_series)
# plot_series_ode(time_series) 

# ##################################################################################
# # Example 5: Simulation with additional kinetics and random parameters
# ##################################################################################
# rate_list = ['mmk', 'hill', 'mak', 'ping_pong', 'threshold']

# additional_laws = {'ping_pong': rate_ping_pong, 'threshold': rate_thresholds}

# time_series, flux_vector = simulation(rn, rate=rate_list, additional_laws=additional_laws)

# print("ODE Time Series with Custom Kinetics:")
# print(time_series)
# plot_series_ode(time_series)