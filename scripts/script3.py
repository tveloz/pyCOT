# # # Script 3: Hierarchies Analysis
# This script focuses on the analysis of hierarchies within a reaction network (RN).
# The script is structured to:  
#   (a) load reaction network
#   (b) calculate semi-organisations
#   (c) plotting semi-organisations
#   (d) generate time series of abstractions for random dynamics and for mak dynamics
#   (e) plot static dynamics graphs in both cases
#   (f) plot dynamic movie in both cases

# # Note: To comment a line of code, use Ctrl + C, and to uncomment it, use Ctrl + K.

#####################################################################
# Imports from the pyCOT library
#####################################################################
# Import Python standard library modules
import os  # Imports the 'os' module to interact with the operating system
import sys # Imports the 'sys' module to interact with the Python interpreter

# Add the root directory to the PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))) # Adds the project's root directory to PYTHONPATH to allow importing modules from that location
sys.stdout.reconfigure(encoding='utf-8')                                     # Set the standard output encoding to UTF-8 (please comment out if using Windows)

# Imports from the pyCOT library 
from pyCOT.file_manipulation import load_pyCOT_from_file # Import only the load_pyCOT_from function of the file_manipulation module to load RN  
from pyCOT.rn_visualize import *                         # Imports all functions, classes and variables defined in the rn_visualize module

from pyCOT.simulations import *           # General simulation tools for dynamic systems modeled in pyCOT.
from pyCOT.abstractions import *          # Provides methods to compute abstractions, such as threshold-based state simplifications.
from pyCOT.plot_dynamics import *         # Functions for plotting dynamic behaviors, such as time series or abstracted data.
from pyCOT.rn_types import StoichiometryMatrix
from pyCOT.abstractions import abstraction_ordinary 
from pyCOT.closure_structure import *     # Contains functions to compute and analyze closure structures in networks.

#####################################################################
# (a) Load a reaction network
#####################################################################
# # File path
#file_path = '../Txt/autopoietic.txt' 
# file_path = 'Txt/2019fig1.txt'
# file_path = 'Txt/2019fig2.txt'
# file_path = '../Txt/non_connected_example.txt' 
# file_path = '../Txt/PassiveUncomforableIndignated_problemsolution.txt'
file_path = '../Txt/Farm.txt' 

# # Loads the RN from the specified file
testRN = load_pyCOT_from_file(file_path)           # Creates an object called testRN, which represents the RN

#####################################################################
# (b) calculate semi-organisations 
#####################################################################
# # Calculate the semi-organisations
input_data = reactive_semi_orgs(testRN)   
print("Semi-organisations:", input_data) 

#####################################################################
# (c) plotting semi-organisations 
#####################################################################
# # Hierarchy of semi-organisations 
# hierarchy_get_visualization_html(input_data)  # Generate the HTML of hierarchy visualisation
# hierarchy_visualize_html(input_data)          # Open the HTML hierarchy visualisation
"""
The hierarchy visualisation shows the semi-organisations of the reaction network (RN) as a directed graph.
The nodes represent the semi-organisations, and the edges represent the relationships between them.
By default in hierarchy_visualize_html():
    Nodes are represented by circular nodes colored 'cyan' and 
    Edges are represented by arrows colored 'gray'.
The node size by default is of size 20, and the edge width is of size 2.

"""

# # Hierarchy of semi-organisations with lst_color_subsets (Farm.txt)  
# lst_color_subsets = [('red', [['water'], 
#                               ['chickens', 'eggs', 'infr', 'water', 'grain', 'straw', 'fertilizer'], 
#                               ['chickens', 'eggs', 'infr', 'grass', 'grain', 'straw', 'water', 'fertilizer'], 
#                               ['chickens', 'money', 'infr', 'eggs', 'grass', 'water', 'grain', 'straw', 'fertilizer', 'farmer']])] # List of tuples with color and subsets of species
# hierarchy_visualize_html(input_data, node_color="blue", lst_color_subsets=lst_color_subsets) # Open the HTML hierarchy visualisation with the specified color for the nodes and list of color subsets

#####################################################################
# (d) generate time series of abstractions for random dynamics and for mak dynamics 
#####################################################################
# # Generates an initial state vector based on the number of species in the reaction network
x_inicial_random = generate_random_vector(len(testRN.SpStr))  
print(x_inicial_random)

# Calculates the stoichiometric matrix associated with the reaction network
matrix_data = universal_stoichiometric_matrix(testRN) # Calculates the universal stoichiometric matrix associated with the reaction network
print(matrix_data)
S = StoichiometryMatrix(matrix_data, species=testRN.SpStr, reactions=testRN.RnStr) # Creates a StoichiometryMatrix object with the required set of species and reactions
print(S)

############################################################
# # Simulates the time series of abstractions for random dynamics
############################################################ 
random_time_series, random_flux_series = simulate_discrete_random(testRN, S, x_inicial_random, n_iter=100)

# # Define a fixed threshold or per-species thresholds
threshold_random = 0.1  # Fixed threshold for random dynamics

# # Generates an abstraction for random dynamics 
# random_abstract_time_series = abstraction_ordinary(random_time_series)          # Generates an abstraction with a default threshold of 0.5
random_abstract_time_series = abstraction_ordinary(random_time_series, threshold_random) # Generates an abstraction with a defined threshold for each species.
print("Random_abstract_time_series:")
print(random_abstract_time_series)

############################################################
# # Simulates the time series of abstractions for mak dynamics
############################################################
mak_time_series, mak_flux_series = simulate_ode_mak(testRN) 

# # Define a fixed threshold or per-species thresholds
threshold_mak = 0.1  # Fixed threshold for mak dynamics

# # Generates an abstraction for mak dynamics 
# abstract_time_series = abstraction_ordinary(time_series)          # Generates an abstraction with a default threshold of 0.5
mak_abstract_time_series = abstraction_ordinary(mak_time_series, threshold_mak) # Generates an abstraction with a defined threshold for each species.
print("MAK_abstract_time_series:")
print(mak_abstract_time_series)

#####################################################################
# (e) plot static dynamics graphs in both cases 
##################################################################### 
# # Static graph with nodes of semi-organisations and abstractions 
plot_static_abstraction_graph_hierarchy_html(random_abstract_time_series)
plot_static_abstraction_graph_hierarchy_html(mak_abstract_time_series)            
 
##################################################################### 
# (f) plot dynamic movie in both cases
#####################################################################
# # # Movie  
plot_abstraction_graph_movie_html(random_abstract_time_series) # Open the HTML file with the animation of the abstraction graph
plot_abstraction_graph_movie_html(mak_abstract_time_series)    # Open the HTML file with the animation of the abstraction graph

# # # Movie with nodes of semi-organisations and abstractions 
film_semiorganizations_abstractions_html(random_abstract_time_series,input_data) # Open the HTML file with the animation of the abstraction graph
film_semiorganizations_abstractions_html(mak_abstract_time_series,input_data)    # Open the HTML file with the animation of the abstraction graph

# # Plot Hasse diagram 
plot_hasse_diagram(testRN.SpStr) # Only can be plotted for RN with less than 10 species, because the number of nodes would be 2^n.