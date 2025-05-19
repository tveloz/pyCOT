#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 15:28:44 2023

@author: tveloz
"""

import os
import sys
# Añadir el directorio raíz al PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# main.py (or another script)
from pyCOT.rn_visualize import rn_get_string
from pyCOT.reaction_network import *
from pyCOT.closure_structure import * 
import networkx as nx
from pyCOT.file_manipulation import *
import matplotlib.pyplot as plt
import time
import os

# Get the current working directory
#current_directory = os.getcwd()
#print("Current Directory:", current_directory)
#Create an instance of the HelloWorld class
# l1=["s1","s2"]
# l2=["s3","s2"]
# l=[l1,l2,["s3","s4"],["s1","s2"]]
# lp=power_list(l)
# for x in lp:
#     print(x)


# print(['s1','s2']==['s2','s1'])
# print(set(['s1','s2'])==set(['s2','s1']))

# SpStr= ['a', 'b', 'c', 'd']  # Default: ['a', 'b', 'c', 'd']
# SpBt=bt([True,True,True,True])  # Default: [a, b, c, d]
# RnStr= ['r1', 'r2']  # Updated: ['r1', 'r2']
# RnBt= bt([True, True])  # Default: [r1, r2]
# RnVecS = np.array([[1, 0, 0, 0], [0, 1, 0, 0]])  # Updated for reactions
# RnVecP= np.array([[0, 2, 0, 0], [0, 1, 1, 0]])  # Updated for reactions


from pyCOT.io.functions import read_txt
import numpy as np

# 1. Build a reaction network object from the file
file_path = "networks/testing/autopoietic.txt"
rn = read_txt(file_path)

# 2. Print the reactions as in the txt file
print(rn_get_string(rn))

# 3. Select a set of species (for example, s1 and l)
species = ["l","s2"]
selected_species =rn.generated_closure(species)
print(f"\nSelected species: {selected_species}")

# 4. Build the new reaction network object associated to the reactions from such species
# Get reactions activated by the selected species
activated_reactions = rn.get_reactions_from_species(selected_species)

# Create a new reaction network containing only the subset of species and their reactions
sub_rn = type(rn)()  # Create a new instance of the same class

# First add all involved species
all_species = set()
for species_name in selected_species:
    all_species.add(species_name)

# Then add species that appear in the activated reactions
for reaction in activated_reactions:
    # Add species from support
    for edge in reaction.support_edges():
        species = rn.get_species_by_index(edge.source_index)
        all_species.add(species.name)
    # Add species from products
    for edge in reaction.products_edges():
        species = rn.get_species_by_index(edge.target_index)
        all_species.add(species.name)

# Add all species to the new reaction network
for species_name in all_species:
    sub_rn.add_species(species_name)

# Add all activated reactions to the sub-network
for reaction in activated_reactions:
    # Get reaction details
    reaction_name = reaction.name()
    
    # Get support species and coefficients
    support_species = []
    support_coefficients = []
    for edge in reaction.support_edges():
        species = rn.get_species_by_index(edge.source_index).name
        support_species.append(species)
        support_coefficients.append(edge.coefficient)
    
    # Get product species and coefficients
    product_species = []
    product_coefficients = []
    for edge in reaction.products_edges():
        species = rn.get_species_by_index(edge.target_index).name
        product_species.append(species)
        product_coefficients.append(edge.coefficient)
    
    # Add the reaction to the sub-network
    sub_rn.add_reaction(
        reaction_name, 
        support_species, 
        product_species, 
        support_coefficients, 
        product_coefficients
    )

# Print the sub-network
print(rn_get_string(sub_rn))


# 5. Print the stoichiometric matrix of the loaded RN
print("\nStoichiometric Matrix of Original Reaction Network:")
stoich_matrix = rn.stoichiometry_matrix()
print("Species:", stoich_matrix.species)
print("Reactions:", stoich_matrix.reactions)
print("Matrix:")
print(stoich_matrix)

# 6. Print the stoichiometric matrix of the smaller RN for the selected set of species
print("\nStoichiometric Matrix of Sub-network:")
sub_stoich_matrix = sub_rn.stoichiometry_matrix()
print("Species:", sub_stoich_matrix.species)
print("Reactions:", sub_stoich_matrix.reactions)
print("Matrix:")
print(sub_stoich_matrix)
# for sorgs in reactive_sorgs:
#    print(sorgs)

#print(testRN.get_reactions_from_bt(bt('0100')))
# print("specs")
# print(testRN.SpStr)
# print("specsBt")
# print(testRN.SpBt)
# print("reacsBt")
# print(testRN.RnBt)
# print("Reac")
# print(testRN.RnStr)
# print("SuppVec")
# print(testRN.RnVecS)
# print("ProdVec")
# print(testRN.RnVecP)
#Example usage:

#file_path = '../../networks/testing/RedPDoSR00.txt'
#file_path = '../networks/biomodels_interesting/BIOMD0000000652_manyOrgs.xml'
#file_path = '../networks/biomodels_interesting/Biomodel_724.txt'
# file_path = '../networks/testing/ERCs_test.txt'
# file_path = '../networks/testing/Synergy_test.txt'
#file_path = '../networks/testing/Farm.txt'
# #file_path = '../networks/testing/RedPDoSR01.txt'
# #file_path = '../../networks/testing/MSORN_test1.txt'
# #file_path = '../../networks/testing/MSORN_test2.txt'
#file_path = '../networks/biomodels_interesting/central_ecoli.txt'
# file_path = '../networks/biomodels_all/BIOMD0000000086_url.xml'
#file_path = '../../networks/testing/MSORN_test2.txt'
#msorn= load_pyCOT_from_file(file_path)#load_pyCOT_from_Sbml(file_path)
#testRN2 = load_pyCOT_from_Sbml(file_path)
#testRN2 = load_pyCOT_from_file(file_path)
#print(str(testRN2.SpStr))
#print(str(testRN2.get_reactions_consuming_species("d")))
#print(str(msorn.get_connected_species_to_species("E4")))



