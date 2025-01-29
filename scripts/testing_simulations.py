import os
import sys
# Add the root directory to the PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pyCOT.simulations import *
from pyCOT.reaction_network import *
from pyCOT.closure_structure import * 
from pyCOT.file_manipulation import *   
from pyCOT.plot_dynamics import * 
from pyCOT.abstractions import abstraction_ordinary  
from pyCOT.rn_types import StoichiometryMatrix

# Load the reaction network
file_path = 'Txt/autopoietic.txt' 
# file_path = 'Txt/2019fig1.txt'
# # file_path = 'Txt/2019fig2.txt' # No funcionó
# file_path = 'Txt/non_connected_example.txt' 
# file_path = 'Txt/PassiveUncomforableIndignated_problemsolution.txt'
# file_path = 'Txt/Farm.txt'
# # file_path = 'Txt/RN_IN_04.txt' # No grafica la EDO, o depende de la escala 

# Generate the pyCOT object from the file containing the txt network
testRN = load_pyCOT_from_file(file_path)

#####################################################################################
# # Example 1.1: Time series of ODE with default parameters and initial conditions
###################################################################################### 
# # Solve EDO 
# time_series = simulate_ode_mak(testRN)
# print("Time series of ODE:")
# print(time_series)

# # Plots the time series of ODE
# plot_series_ode(time_series)

###################################################################################### 
# # Example 1.2: Time series of ODE with parameters and initial conditions with different seeds
###################################################################################### 
# # Creación automática de los valores de las constantes de velocidad
# # k0=generate_random_vector(len(testRN.RnStr),seed=3).round(2) 
# k0=generate_random_vector(len(testRN.RnStr)).round(4) 

# # Creación automática del diccionario
# k= {testRN.RnStr[i]: float(k0[i]) for i in range(len(testRN.RnStr))}

# # Mostrar el diccionario resultante
# print("Specific reaction rate constants:") 
# print(k)

# # Time span
# t_span = (0, 200)
 
# # Creación automática de los valores iniciales de las concentraciones
# # y0=generate_random_vector(len(testRN.SpStr),seed=3).round(1)
# x0=generate_random_vector(len(testRN.SpStr)).round(2)
# print("Initial concentrations for species:")
# print(x0)
 
# # Solve EDO
# time_series = simulate_ode_mak(testRN, k, x0, t_span, n_steps=500) 
# print("Time series of ODE:")
# print(time_series)

# # Plots the time series of ODE
# plot_series_ode(time_series)

######################################################################################
# # Example 1.3: PassiveUncomforableIndignated_problemsolution
######################################################################################
# # Load the reaction network
# file_path = 'Txt/PassiveUncomforableIndignated_problemsolution.txt' 
# testRN = load_pyCOT_from_file(file_path)
 
# # Define specific reaction rate constants
# k = {
#     testRN.RnStr[0]: 0.1,
#     testRN.RnStr[1]: 0.1,
#     testRN.RnStr[2]: 0.2,
#     testRN.RnStr[3]: 0.3,
#     testRN.RnStr[4]: 0.3, 
#     testRN.RnStr[5]: 0.3,
#     testRN.RnStr[6]: 0.3, 
#     testRN.RnStr[7]: 0.3, 
#     testRN.RnStr[8]: 0.15,
#     testRN.RnStr[9]: 0.15,                    
# }

# # Time span
# t_span = (0, 300)

# # Define initial concentrations for species
# y0 = [
#     0.3,  # Initial concentration for species 1: D
#     0.4,  # Initial concentration for species 2: U
#     0.3,  # Initial concentration for species 3: S
#     0.1,  # Initial concentration for species 4: s
#     0.9,  # Initial concentration for species 5: p
# ]

# ######################################################################################
# # Solve EDO
# time_series = simulate_ode_mak(testRN, k, t_span, x0,n_steps=1000)
# print(time_series)

# # Plots the time series of ODE
# plot_series_ode(time_series, xlabel="Time", ylabel="Concentration", title="Time Series of Concentrations")

########################################################################################
########################################################################################
########################################################################################
# # # Example 3: 
# # file_path = 'Txt/autopoietic.txt'
# # file_path = 'Txt/2019fig1.txt'
# # file_path = 'Txt/2019fig2.txt'
# # file_path = 'Txt/non_connected_example.txt'
# file_path = 'Txt/Farm.txt' 

# # # Generate the pyCOT object from the file containing the txt network
# testRN=load_pyCOT_from_file(file_path)  
# print(testRN.SpStr)
 
# # # Generates an initial state vector based on the number of species in the reaction network
# x_inicial = generate_random_vector(len(testRN.SpStr))  
# print(x_inicial)

# # # Calculates the stoichiometric matrix associated with the reaction network
# matrix_data=universal_stoichiometric_matrix(testRN) 
# print(matrix_data)

# S=StoichiometryMatrix(matrix_data, species=testRN.SpStr, reactions=testRN.RnStr) 
# print(S)

# # # Simulates the time series for 20 iterations, using the stoichiometric matrix and the initial vector of states
# time_series, reaction_rate_maks = simulate_discrete_random(S, x_inicial, testRN.SpStr, n_iter=300)
# print(time_series)
# print(reaction_rate_maks)

# # # Define a fixed threshold or per-species thresholds
# threshold = 0.05  # Fixed threshold
# # threshold = {"A": 0.4, "B": 0.5, "C": 0.4}  # Per-species thresholds

# # # Calculate abstractions
# # abstract_time_series = abstraction_ordinary(time_series)
# abstract_time_series = abstraction_ordinary(time_series, threshold)
# print(abstract_time_series)

# Plots the time series of ODE
# plot_series_ode(time_series, xlabel="Time", ylabel="Concentration", title="Time Series of Concentrations")
# plot_abstraction_size(abstract_time_series, xlabel="Time", ylabel="Number of Species", title="Number of species per abstraction over time", marker='o', label="Abstraction Size")
# plot_abstraction_sets(abstract_time_series)

# # # Graphs
# # Static graph
# plot_static_abstraction_graph(abstract_time_series)  

# # Movie
# plot_abstraction_graph_movie(abstract_time_series, interval=1000)

# # Histograms
# plot_join_concentration_histofram(time_series, bins=30, alpha=0.7, color='skyblue', edgecolor='black', xlabel="Concentration", ylabel="Frequency", title="Histogram of Species Concentrations")
# plot_reaction_rate_mak_histogram(reaction_rate_maks, bins=10, color='skyblue', edgecolor='black', alpha=0.7, xlabel="Reaction Rate", ylabel="Frequency", title="Histogram of Reaction Rates")
# plot_species_histograms(time_series, species_names=testRN.SpStr, bins=10, alpha=0.7, color="skyblue", edgecolor="black", xlabel="Concentration", ylabel="Frequency", title_plot="Histogram of")
# plot_combined_species_histogram(time_series, species_names=testRN.SpStr, bins=10, alpha=0.7, edgecolor="black", xlabel="Concentration", ylabel="Frequency", title_plot="Combined Histogram")

########################################################################################
########################################################################################
########################################################################################
# # # Example 4: Michaelis-Menten Kinetics
# Define specific reaction rate constants
k0 = {
    testRN.RnStr[0]: 0.65,
    testRN.RnStr[1]: 1.0,
    testRN.RnStr[2]: 2.0,
    testRN.RnStr[3]: 1.5,
    testRN.RnStr[4]: 2.5,    
}
 
# Define specific reaction rate constants
Vmax_dict = {
    testRN.SpStr[0]: 1.0,
    testRN.SpStr[1]: 1.2,
    testRN.SpStr[2]: 0.8,    
}

# Michaelis constants
Km_dict = {
    testRN.SpStr[0]: 0.5,
    testRN.SpStr[1]: 0.7,
    testRN.SpStr[2]: 0.6,    
}

# Initial conditions
x0 = [1.0, 2.0, 3.5]  # Initial concentrations of l, s1, s2

# Time range for the simulation
t_span = (0, 20)

n_steps = 200
#####################################################################################
# Solve the ODEs
time_series, flux_vector  = simulate_ode_mmk(testRN) # Default values for all Reaction Network
# time_series, flux_vector  = simulate_ode_mmk(testRN, x0) 
# time_series, flux_vector  = simulate_ode_mmk(testRN, x0, t_span, n_steps, k0)
# time_series, flux_vector  = simulate_ode_mmk(testRN, x0, t_span, n_steps, k0, Vmax_dict, Km_dict)

#####################################################################################
# Plot the time series of concentration and flux for the Michaelis-Menten system
plot_series_ode(time_series, title="Michaelis-Menten MMK: ODE")
# plot_series_ode(flux_vector, title="Michaelis-Menten MMK: ODE")

########################################################################################
########################################################################################
########################################################################################
# # # Example 5: Metapolulation
# # Initial configuration
# grid_shape=(1, 5)

# modules = simulate_odes_metapop_mak(testRN, grid_shape=grid_shape) 

# plot_series_MP(modules,testRN)