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

#####################################################################################
# # Example 1.1: autopoietic 
######################################################################################
# # Load the reaction network
# file_path = 'networks/testing/autopoietic.txt' 
# testRN = load_pyCOT_from_file(file_path)
 
# # Define specific reaction rate constants
# k = {
#     testRN.RnStr[0]: 0.5,
#     testRN.RnStr[1]: 1.0,
#     testRN.RnStr[2]: 2.0,
#     testRN.RnStr[3]: 1.5,
#     testRN.RnStr[4]: 2.5,    
# } 

# # Time span
# t_span = (0, 40)

# # Define initial concentrations for species
# y0 = [
#     0.8,  # Initial concentration for species 1: l
#     0.3,  # Initial concentration for species 2: s1
#     0.5,  # Initial concentration for species 3: s2
# ]
 
# # Solve EDO
# time_series = solve_ode(testRN, k, t_span, y0,n_steps=500)
# print(time_series)

# # Plots the time series of ODE
# plot_series_ode(time_series)

###################################################################################### 


 


# #####################################################################################
# # Example 1.2: PassiveUncomforableIndignated_problemsolution
######################################################################################
# # Load the reaction network
# file_path = 'networks/testing/PassiveUncomforableIndignated_problemsolution.txt' 
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
# time_series = solve_ode(testRN, k, t_span, y0,n_steps=1000)
# print(time_series)

# # Plots the time series of ODE
# plot_series_ode(time_series, xlabel="Time", ylabel="Concentration", title="Time Series of Concentrations")









########################################################################################
########################################################################################
########################################################################################
# # # Example 3: 
# file_path = 'Txt/autopoietic.txt'
# # file_path = 'Txt/2019fig1.txt'
# # file_path = 'Txt/2019fig2.txt'
# # file_path = 'Txt/non_connected_example.txt'
file_path = 'networks/testing/Farm.txt' 

# # Generate the pyCOT object from the file containing the txt network
testRN=load_pyCOT_from_file(file_path)  
print(testRN.SpStr)
 
# # Generates an initial state vector based on the number of species in the reaction network
x_inicial = generate_state_vector(len(testRN.SpStr))  
print(x_inicial)

# # Calculates the stoichiometric matrix associated with the reaction network
matrix_data=stoichiometric_matrix(testRN) 
print(matrix_data)
S=StoichiometryMatrix(matrix_data, species=testRN.SpStr, reactions=testRN.RnStr) 
print(S)

# # Simulates the time series for 20 iterations, using the stoichiometric matrix and the initial vector of states
time_series, reaction_rates = iterate_state_vector(S, x_inicial, testRN.SpStr, n_iter=20)
print(time_series)
print(reaction_rates)

# # Define a fixed threshold or per-species thresholds
threshold = 0.5  # Fixed threshold
# threshold = {"A": 0.4, "B": 0.5, "C": 0.4}  # Per-species thresholds

# # Calculate abstractions
# abstract_time_series = abstraction_ordinary(time_series)
abstract_time_series = abstraction_ordinary(time_series, threshold)
print(abstract_time_series)

# Plots the time series of ODE
plot_series_ode(time_series, xlabel="Time", ylabel="Concentration", title="Time Series of Concentrations")
plot_abstraction_size(abstract_time_series, xlabel="Time", ylabel="Number of Species", title="Number of species per abstraction over time", marker='o', label="Abstraction Size")
plot_abstraction_sets(abstract_time_series)

# # Graphs
# Static graph
plot_static_abstraction_graph(abstract_time_series)  

# Movie
plot_abstraction_graph_movie(abstract_time_series, interval=1000)

# Histograms
plot_join_concentration_histofram(time_series, bins=30, alpha=0.7, color='skyblue', edgecolor='black', xlabel="Concentration", ylabel="Frequency", title="Histogram of Species Concentrations")
plot_reaction_rate_histogram(reaction_rates, bins=10, color='skyblue', edgecolor='black', alpha=0.7, xlabel="Reaction Rate", ylabel="Frequency", title="Histogram of Reaction Rates")
plot_species_histograms(time_series, species_names=testRN.SpStr, bins=10, alpha=0.7, color="skyblue", edgecolor="black", xlabel="Concentration", ylabel="Frequency", title_plot="Histogram of")
plot_combined_species_histogram(time_series, species_names=testRN.SpStr, bins=10, alpha=0.7, edgecolor="black", xlabel="Concentration", ylabel="Frequency", title_plot="Combined Histogram")