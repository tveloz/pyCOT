#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 15:28:44 2023

@author: tveloz
"""

# main.py (or another script)

from pyCOT.reaction_network import *
from pyCOT.closure_structure import *
import networkx as nx
from pyCOT.file_manipulation import *
import matplotlib.pyplot as plt
import time

# Create an instance of the HelloWorld class
# SpStr= ['a', 'b', 'c', 'd']  # Default: ['a', 'b', 'c', 'd']
# SpBt=bt([True,True,True,True])  # Default: [a, b, c, d]
# RnStr= ['r1', 'r2']  # Updated: ['r1', 'r2']
# RnBt= bt([True, True])  # Default: [r1, r2]
# RnVecS = np.array([[1, 0, 0, 0], [0, 1, 0, 0]])  # Updated for reactions
# RnVecP= np.array([[0, 2, 0, 0], [0, 1, 1, 0]])  # Updated for reactions

# testRN=pyCOT(SpStr,SpBt,RnStr,RnBt,RnVecS,RnVecP)
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
file_path = '../../networks/testing/in_out.txt'
file_path = '../../networks/testing/RedPDoSR00.txt'
#file_path = '../../networks/biomodels_interesting/BIOMD0000000652_manyOrgs.xml'
#file_path = '../../networks/biomodels_interesting/Biomodel_724.txt'
#file_path = '../../networks/testing/ERCs_test.txt'
#file_path = '../../networks/testing/Synergy_test.txt'
#file_path = '../../networks/testing/Farm.txt'
#file_path = '../../networks/testing/MSORN_test1.txt'
#file_path = '../../networks/testing/MSORN_test2.txt'
#file_path = '../../networks/RandomAlife/RN_Ns_40_Norg_24_id_316.txt'
#file_path = '../../networks/biomodels_interesting/central_ecoli.txt'
#file_path = '../../networks/biomodels_all/BIOMD0000000086_url.xml'
#file_path = '../../networks/biomodels_interesting/BIOMD0000000237_manyOrgs.xml'
testRN2 = load_pyCOT_from_file(file_path)
#testRN2 = load_pyCOT_from_Sbml(file_path)

# reaction_graph = build_graph(testRN2)
# visualize_graph(reaction_graph)
# print(testRN2.SpStr)
print("specs")
print(testRN2.SpStr)
print("specsBt")
print(testRN2.SpBt)
print("reacsBt")
print(testRN2.RnBt)
print("Reac")
print(testRN2.RnStr)
print("SuppVec")
print(testRN2.RnVecS)
print("ProdVec")
print(testRN2.RnVecP)
#reaction_graph = build_graph(testRN2)
#visualize_graph(testRN2)

testRN2.is_connected(["Od"])

#Example visualization of a reaction network
# reaction_network = [
#     ['A', 'B'],
#     ['B', 'C'],
#     ['A', 'D', 'E'],
#     ['B', 'D', 'F']
# ]

# reaction_network=gen[2]
# #print(reaction_network)
# # Build and visualize the graph
# reaction_graph = build_graph(reaction_network)
# visualize_graph(reaction_graph)

# Example usage
 # Node representing the species you're interested in

#CBL=Closed_by_levels(testRN2)
#print(str(CBL))

# Example usage:

# Accessing attributes of the instance
# print("Species Bitarray:", testRN2.SpBt)
# print("Species Vector:", RN.SpVec)
# print("Species Names:", RN.SpStr)
# print("Reaction Bitarray (Require Support):", RN.RnBtS)
# print("Reaction Bitarray (Produce Support):", RN.RnBtP)
# print("Reaction Vector (Require Support):", RN.RnVecS)
# print("Reaction Vector (Produce Support):", RN.RnVecP)
# print("Reaction Names:", RN.RnStr)
# mgen=minimal_generators(testRN2)
# print("the list of minimal generators")

# for g in mgen:
#     print(g)
# print(len(gen))
# print(len(mgen))
#print(len(mgen))
#print("specs")
#print(testRN2.SpStr)
# print("specsBt")
# print(testRN2.SpBt)
# print("reacsBt")
# print(testRN2.RnBt)
# print("Reac")
# print(testRN2.RnStr)
# print("SuppVec")
# print(testRN2.RnVecS)
# print("ProdVec")
# print(testRN2.RnVecP[1])
# print("bt_from_species")
# print(testRN2.get_bt_from_species([]))
# print("species_from_bt")
# print(testRN2.get_species_from_bt(bt('0001')))

# print("bt_from_reactions")
# print(testRN2.get_bt_from_reactions(['R1']))

# print("reactions_from_bt")
# print(testRN2.get_reactions_from_bt(bt('00000000')))
# print("support_bt_from_reaction")
#print(testRN2.get_supp_bt_from_reaction('R1'))
# print("products_bt_from_reaction")
# print(testRN2.get_prod_bt_from_reaction('R8'))
# print(str(testRN2.SpStr[1:3]))
# print("reactions_from_species")
# print(testRN2.get_reactions_from_species([]))

# print("get prod from reactions")
# print(testRN2.get_prod_from_reactions([]))
# print("get supp from reactions")
# print(testRN2.get_supp_from_reactions([]))
# file_path = '../../networks/autopoietic_ext2.txt'
# testRN2 = load_pyCOT_from_file(file_path)
# print("bt_from_reactions ")
# print(testRN2.get_bt_from_reactions([]))
# print("bt_from_species")
# print(testRN2.get_bt_from_species([]))


# print(testRN2.RnStr)

# print("get inflow")
# print(testRN2.get_inflow())

# print("get outflow")
# print(testRN2.get_outflow())


# print("get_supp_from_reactions")
# print(testRN2.get_supp_from_reactions(testRN2.RnStr[4]))

# print("get_prod_from_reactions "+str(testRN2.RnStr[0:2]))
# print(testRN2.get_prod_from_reactions(['R1']))
# print("get_reactions_from_species")
# print(testRN2.get_reactions_from_species([]))

# print("get get_prod_from_species")
# print(testRN2.get_prod_from_species(['l','s2']))
# print("get_reactions_consuming_species")
# print(testRN2.get_reactions_consuming_species(['l','s1']))      
# print("get_reactions_producing_species")
# print(testRN2.get_reactions_producing_species(['l','s1']))        

#print("get connected species to species")
#print(testRN2.get_connected_species_to_species(['l']))

# # print("get directly connected species to species")
# # print(testRN2.get_immediately_connected_species_to_species(['s2','x']))
# print("get forward connected species to species")
# print(testRN2.get_forward_connected_species_to_species('x'))
# print("get reactions consuming")
# print(testRN2.get_reactions_consuming_species('s1'))

# print("get reactions producing")
# print(testRN2.get_reactions_producing_species('s1'))

# print("is_closed")
# print(testRN2.is_closed(['s1','s2']))
# print("is_semi_self_maintaining")
# print(testRN2.is_semi_self_maintaining(['d','s1']))

