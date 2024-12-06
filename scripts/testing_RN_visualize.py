######################################################################################
import os
import sys
# Add the root directory to the PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pyCOT.rn_visualize import *
from pyCOT.reaction_network import *
from pyCOT.closure_structure import * 
from pyCOT.file_manipulation import *  

#####################################################################
# # File path
file_path = 'Txt/autopoietic.txt' 
# file_path = 'Txt/2019fig1.txt'
# file_path = 'Txt/2019fig2.txt'
# file_path = 'Txt/non_connected_example.txt' 
# file_path = 'Txt/PassiveUncomforableIndignated_problemsolution.txt'
# file_path = 'Txt/Farm.txt'

# # Loads the RN from the specified file
testRN=load_pyCOT_from_file(file_path) # Creates an object (testRN) representing the RN.

# print("Species:",testRN.SpStr)
# print("Binary vector of species:",testRN.SpBt)
# print("Reactions:",testRN.RnStr) 
# print("Binary vector of reactions:",testRN.RnBt)
# print("Reactant Coefficients:",testRN.RnMsupp) 
# print("Product Coefficients:",testRN.RnMprod) 
# print(dir(testRN)) # List of all attributes and methods of the object testRN

#####################################################################
# # Example 1: Reaction network visualisations
#####################################################################
# rn_get_visualization(testRN) #Generates the visualisation of the reaction network without showing it directly
rn_visualize(testRN) #Shows the reaction network visualisation file

# # Reaction network visualisations with options
# rn_visualize(testRN, physics_enabled=True) 
# rn_visualize(testRN, curvature=True) #None=False (Without curvature), True=curvedCCW (Counter Clockwise) 
# rn_visualize(testRN, filename='rn_visualize.html') 

# # Lists of species sets and reactions with colors
# S=[["water"],["eggs","chickens","infr"]]
# lst_color_specs = [("blue", S[0]), ("yellow", S[1])]
# lst_color_reacs = [("purple", ["R1"]), ("orange", ["R2","R3","R4"]), ("green", ["R12","R13","R14"])]

# # Visualize the network with lst_color_spcs and lst_color_reacs applied 
# rn_visualize(testRN, lst_color_spcs=[("yellow", ["l","s1"])], filename="rn_visualize.html")
# rn_visualize(testRN, lst_color_spcs=lst_color_specs, filename="rn_visualize.html")
# rn_visualize(testRN, lst_color_reacs=lst_color_reacs, filename="rn_visualize.html")
# rn_visualize(testRN, lst_color_spcs=lst_color_specs, lst_color_reacs=lst_color_reacs, filename="rn_visualize.html")
# rn_visualize(testRN, global_species_color='yellow', global_reaction_color='orange', filename="rn_visualize.html")
# rn_visualize(testRN, global_species_color='blue', global_reaction_color='orange',global_input_edge_color='blue', global_output_edge_color='gray', filename="rn_visualize.html")
# rn_visualize(testRN, lst_color_spcs=lst_color_specs, lst_color_reacs=lst_color_reacs, global_input_edge_color='blue', global_output_edge_color='gray', filename="rn_visualize.html")

##################################################################
# # Example 2: Hierarchy of semi-organisations
################################################################## 
# # Vector of Semi-Organisations calculated at brute force
# input_data = reactive_semi_orgs(testRN)   
# print("Semi-organisations:", input_data)
# hierarchy_visualize(input_data)
 
##################################################################
# # Example 3: Hierarchy of semi-organisations without duplicates
################################################################## 
# # # Convert each subset to a set and then to a list to remove duplicates
# unique_subsets = []
# for sublist in input_data:
#     if set(sublist) not in [set(x) for x in unique_subsets]:
#         unique_subsets.append(sublist)
# print("Number of Semi-organisations without duplicates=", len(unique_subsets))
# # Run hierarchy_visualize 
# hierarchy_visualize(unique_subsets)

##################################################################
# # Example 4: Hierarchy
##################################################################
# # input_data1 = [[], ['s1', 's2'], ['s1', 's3'], ['s1', 's2', 's3']]
# # input_data1 = [['a'], ['a', 'b'], ['a', 'c', 'd'], ['a', 'b', 'x'], ['a', 'c', 'd', 'y', 'z'], ['a', 'b', 'c', 'd', 'x', 'y', 'z']]
# # # Ejecutar hierarchy_visualize 
# # hierarchy_visualize(input_data1) 
# #####################################################################################
