##############################################################################
# This script uses the pyCOT library to analyse and visualise reaction networks 
# and their associated dynamics.
##############################################################################
import os   # Imports the 'os' module to interact with the operating system
import sys  # Imports the 'sys' module to interact with the Python interpreter
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # Adds the project's root directory to PYTHONPATH to allow importing modules from that location

# # Imports from the pyCOT library 
from pyCOT.reaction_network import *      # Provides tools for creating and manipulating reaction networks.
from pyCOT.closure_structure import *     # Contains functions to compute and analyze closure structures in networks.
from pyCOT.file_manipulation import *     # Includes utilities to load and save data or models for pyCOT.

from pyCOT.rn_visualize import *          # Functions for visualizing reaction networks as graphs or diagrams.
from pyCOT.hierarchy_visualize import *   # Tools for creating hierarchical visualizations of network structures.
from pyCOT.simulations import *           # General simulation tools for dynamic systems modeled in pyCOT.
from pyCOT.abstractions import *          # Provides methods to compute abstractions, such as threshold-based state simplifications.
from pyCOT.plot_dynamics import *         # Functions for plotting dynamic behaviors, such as time series or abstracted data.
from pyCOT.rn_types import StoichiometryMatrix

# # Load the path of the .txt file of the reaction network 
file_path = 'Txt/autopoietic.txt' 

# # Generate the pyCOT object from the .txt file containing the reaction network. 
testRN=load_pyCOT_from_file(file_path)  
print(testRN.SpStr)

##############################################################################
# # Examples
##############################################################################
# Generates an initial state vector based on the number of species in the reaction network
x_inicial = generate_state_vector(len(testRN.SpStr))  
print(x_inicial)

# Calculates the stoichiometric matrix associated with the reaction network
matrix_data=stoichiometric_matrix(testRN) 
print(matrix_data)
S=StoichiometryMatrix(matrix_data, species=testRN.SpStr, reactions=testRN.RnStr) 
print(S)

# Simulates the time series for 20 iterations, using the stoichiometric matrix and the initial vector of states
time_series, flux_raction = iterate_state_vector(S, x_inicial, testRN.SpStr, n_iter=20)
print(time_series)
print(flux_raction)

##############################################################################
# Calculate abstractions
##############################################################################
# # Define a fixed threshold or per-species thresholds
# threshold = 0.1  # Fixed threshold
threshold = {"l": 0.4, "s1": 0.5, "s2": 0.4}  # Per-species thresholds

# # Generates an abstraction with a threshold
# abstract_time_series = abstraction_ordinary(time_series)          # Generates an abstraction with a default threshold of 0.5
abstract_time_series = abstraction_ordinary(time_series, threshold) # Generates an abstraction with a defined threshold for each species.
print(abstract_time_series)

# print(abstract_time_series["Abstraction"])
# hierarchy_visualize(abstract_time_series["Abstraction"].apply(tuple).drop_duplicates())


###########################################################################
# PLOT
###########################################################################
# 1) Generates a plot of the time series
plot_series_ode(time_series)  

# 2) Creates a graph of the size of the set of abstractions over time
plot_abstraction_size(abstract_time_series)

# 3) Plot the abstractions sorted by size
plot_abstraction_sets(abstract_time_series)

###########################################################################
# Graphs
###########################################################################
# # 1) Static network 
# plot_static_abstraction_graph(abstract_time_series)

# # 2) Movie 
# plot_abstraction_graph_movie(abstract_time_series, interval=1000)

###########################################################################
# Solve ODE
###########################################################################
# Define specific reaction rate constants
k = {
    testRN.RnStr[0]: 0.5,
    testRN.RnStr[1]: 1.0,
    testRN.RnStr[2]: 2.0,
    testRN.RnStr[3]: 1.5,
    testRN.RnStr[4]: 2.5,    
} 

# Time span
t_span = (0, 40)

# Define initial concentrations for species
y0 = [
    0.8,  # Initial concentration for species 1: l
    0.3,  # Initial concentration for species 2: s1
    0.5,  # Initial concentration for species 3: s2
]
 
# Solve ODE
time_series = solve_ode(testRN, k, t_span, y0,n_steps=500)
print(time_series)

# Plots the time series of ODE
plot_series_ode(time_series)

###########################################################################
# Histograms
###########################################################################
# # 1) Generates a histogram of all species concentrations in the NR over time
# plot_join_concentration_histofram(time_series, bins=10)

# # Print reaction rates v(x,k)
# reaction_rates = flux_raction["Flux_0"] 
# print(reaction_rates)
# # 2) Generates histogram of reaction rates v(x,k)
# plot_reaction_rate_histogram(reaction_rates, bins=5)

# # 3) Generates individual histograms for each species
# plot_species_histograms(time_series, species_names=testRN.SpStr)

# # 4) Generates a combined histogram of the concentrations of all species.
# plot_combined_species_histogram(time_series, species_names=testRN.SpStr)