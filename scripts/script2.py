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
from pyCOT.plot_dynamics import plot_series_ode
from pyCOT.simulations import *

# ========================================
# 2. CREATING THE REACTION_NETWORK OBJECT
# ========================================
file_path = 'Txt/autopoietic.txt'  
# Alternative examples:
# file_path = 'Txt/2019fig1.txt'
# file_path = 'Txt/2019fig2.txt'
# file_path = 'Txt/non_connected_example.txt' 
# file_path = 'Txt/PassiveUncomforableIndignated_problemsolution.txt'
# file_path = 'Txt/Farm.txt' 

rn = read_txt(file_path)  # Create ReactionNetwork object from file

# ========================================
# 3. SIMULATION OF THE REACTION NETWORK
# ========================================
###################################################################################
# Example 1: ODE simulation using only 'mak' kinetics
###################################################################################
# x0 = [0, 1, 0]  # Initial concentrations
# rate_list = 'mak'  # Kinetics for all reactions
# spec_vector = [[0.7], [0.5], [1.0], [1.0], [1.0]]  # Parameters for each reaction

# # Run simulation
# time_series, flux_vector = simulation(
#     rn, rate=rate_list, spec_vector=spec_vector, 
#     x0=x0, t_span=(0, 20), n_steps=500
# )

# print("ODE Time Series:")
# print(time_series)

# plot_series_ode(time_series)

###################################################################################
# Example 2: ODE simulation using a mix of 'mak' and 'mmk' kinetics
###################################################################################
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

###################################################################################
# Example 3: Simulation with additional kinetics and random parameters
###################################################################################
# x0 = [5, 7, 10]  # Initial concentrations
# rate_list = ['mak', 'mak', 'mmk', 'ping_pong', 'hill']

# additional_laws = {'ping_pong': rate_ping_pong}
# time_series, flux_vector = simulation(
#     rn, rate=rate_list, x0=x0, t_span=(0, 50), n_steps=500, 
#     additional_laws=additional_laws
# )

# print("ODE Time Series with Custom Kinetics:")
# print(time_series)
# plot_series_ode(time_series)

###################################################################################
# Example 4: Simulation with defined parameters for additional kinetics
###################################################################################
# x0 = [0.1, 0.5, 0.2]
# rate_list = ['mak', 'mak', 'mmk', 'ping_pong', 'hill']
# spec_vector = [
#     [0.5],             # mak
#     [0.3],             # mak
#     [1.2, 0.4],        # mmk
#     [0.9, 0.5, 0.1],   # ping_pong
#     [0.1, 0.9, 0.5]    # hill
# ]

# additional_laws = {'ping_pong': rate_ping_pong}
# time_series, flux_vector = simulation(
#     rn, rate=rate_list, spec_vector=spec_vector, 
#     x0=x0, t_span=(0, 50), n_steps=500, 
#     additional_laws=additional_laws
# )

# print("ODE Time Series with Defined Parameters:")
# print(time_series)
# plot_series_ode(time_series)

###################################################################################
# Example 5: Discrete stochastic simulation
###################################################################################
x0 = [5, 7, 10]  # Initial concentrations
ts, vs = simulate_discrete_random(rn, S=rn.stoichiometry_matrix(), x=x0, n_iter=11)

print("\nTime Series (Discrete Random Simulation):")
print(ts)
plot_series_ode(ts)

print("\nFlux Vector (Discrete Random Simulation):")
print(vs)
plot_series_ode(vs)