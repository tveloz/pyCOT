# # # Script 1: Use of the pyCOT Library
# This script demonstrates the fundamental functionalities of the pyCOT library 
# for analyzing and visualizing reaction networks (RN). 
# The script is structured to: 
#    (a) load a reaction network, 
#    (b) plot the network to visualize its structure, 
#    (c) select a specific set of species from the network, and 
#    (d) generate a color-coded plot that highlights the selected species within the reaction network. 
# These steps provide a streamlined workflow for exploring and presenting reaction dynamics efficiently.

# # Note: To comment a line of code, use Ctrl + C, and to uncomment it, use Ctrl + K.

#####################################################################
# Imports from the pyCOT library
#####################################################################
# Import Python standard library modules
import os  # Imports the 'os' module to interact with the operating system
import sys # Imports the 'sys' module to interact with the Python interpreter

# Add the root directory to the PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))) # Adds the project's root directory to PYTHONPATH to allow importing modules from that location
sys.stdout.reconfigure(encoding='utf-8')                                     # Set the standard output encoding to UTF-8 (pls comment out if using Windows)

# Imports from the pyCOT library 
from pyCOT.file_manipulation import load_pyCOT_from_file # Import only the load_pyCOT_from function of the file_manipulation module to load RN  
from pyCOT.rn_visualize import *                         # Imports all functions, classes and variables defined in the rn_visualize module

#####################################################################
# (a) Load reaction network
#####################################################################
# # File path
#file_path = '../Txt/autopoietic.txt' 
#file_path = '../Txt/2019fig1.txt'
#file_path = '../Txt/2019fig2.txt'
# file_path = '../Txt/non_connected_example.txt' 
# file_path = '../Txt/PassiveUncomforableIndignated_problemsolution.txt'
file_path = '../Txt/Farm.txt' 

# # Loads the RN from the specified file
testRN = load_pyCOT_from_file(file_path)           # Creates an object called testRN, which represents the RN

print("Species:",testRN.SpStr)                   # List of species in the reaction network
print("Binary vector of species:",testRN.SpBt)   # Binary vector of species
print("Reactions:",testRN.RnStr)                 # List of reactions in the reaction network
print("Binary vector of reactions:",testRN.RnBt) # Binary vector of reactions
print("Reactant Coefficients:",testRN.RnMsupp)   # Stoichiometric matrix of reactants
print("Product Coefficients:",testRN.RnMprod)    # Stoichiometric matrix of products
print(dir(testRN))                               # List of all attributes and methods of the object testRN

#####################################################################
# (b) Plot reaction network
#####################################################################
# # Reaction network visualisations:
#rn_get_visualization(testRN) # Generates the html of the NR visualisation without showing it directly
rn_visualize_html(testRN)    # Open the html of the NR visualisation
"""
By default in rn_visualize_html(): 
  Species are represented by circular nodes colored 'cyan' and 
  Reactions are represented by rectangular nodes colored 'lightgray'.
Arrows represent interactions within the reaction network (RN):
  Red arrows: represent reactions that consume species.
  Green arrows: represent reactions that produce species.
In addition, the numbers above the arrows correspond to "stoichiometric coefficients". 
If no explicit coefficient is indicated, its value is assumed to be 1.
"""

# # Reaction network visualisations with options
# rn_visualize_html(testRN, node_size=30)                              # Visualise the reaction network with a specified node size (30 in this case)
# rn_visualize_html(testRN, node_size=25, shape_species_node='circle') # Visualise the reaction network with species nodes in a circular shape. Common options include 'dot', 'circle', 'box', 'ellipse', etc. Default is 'dot'.
#rn_visualize_html(testRN, physics_enabled=True)                      # Enables physical effects in the visualisation, so that moving a node moves the entire network
#rn_visualize_html(testRN, curvature=True)                            # Shows the network with curved arrows; ‘True’ applies a counter-clockwise curving
#rn_visualize_html(testRN, filename='my_rn_visualize_html.html')      # Saves the network visualisation as an HTML file with the specified name

#####################################################################
# (c) Select set of species
#####################################################################
# # Species and reaction set lists with specific colours (Farm.txt) 
# lst_color_specs = [("blue", ["water"]), ("yellow", ["eggs","chickens","infr"])]                      # Species set lists with specific colours
# lst_color_reacs = [("purple", ["R1"]), ("orange", ["R2","R3","R4"]), ("green", ["R12","R13","R14"])] # Reaction set lists with specific colours

#####################################################################
# (d) Plot reaction network and set of species with different colour
#####################################################################
# # Visualize the network with lst_color_spcs and lst_color_reacs applied  
# rn_visualize_html(testRN, lst_color_spcs=lst_color_specs)                                   # Visualises the reaction network with personalised colours for some species sets
# rn_visualize_html(testRN, lst_color_reacs=lst_color_reacs)                                  # Visualises the reaction network with personalised colours for some reaction sets
#rn_visualize_html(testRN, lst_color_spcs=lst_color_specs, lst_color_reacs=lst_color_reacs)  # Visualises the reaction network with personalised colours for both some sets of species and reactions

# # Reaction network visualisations with options 
#rn_visualize_html(testRN, global_species_color='yellow', global_reaction_color='orange')  # Visualiza la red de reacciones asignando el color amarillo a las especies y el color naranja a las reacciones
#rn_visualize_html(testRN, global_species_color='blue', global_reaction_color='orange',global_input_edge_color='blue', global_output_edge_color='gray')
#rn_visualize_html(testRN, lst_color_spcs=lst_color_specs, lst_color_reacs=lst_color_reacs, global_input_edge_color='blue', global_output_edge_color='gray')

