#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 15:28:44 2023

@author: tveloz
"""

# main.py (or another script)

from pyCOT_constructor import *
import networkx as nx
from File_manipulation import *

# Create an instance of the HelloWorld class
SpStr= ['a', 'b', 'c', 'd']  # Default: ['a', 'b', 'c', 'd']
SpBt=bt([True,True,True,True])  # Default: [a, b, c, d]
RnStr= ['r1', 'r2']  # Updated: ['r1', 'r2']
RnBt= bt([True, True])  # Default: [r1, r2]
RnVecS = np.array([[1, 0, 0, 0], [0, 1, 0, 0]])  # Updated for reactions
RnVecP= np.array([[0, 2, 0, 0], [0, 1, 1, 0]])  # Updated for reactions

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
# Example usage:
file_path = '../../networks/autopoietic_ext.txt'
testRN2 = load_pyCOT_from_file(file_path)
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
print(testRN2.get_bt_from_species(['d','l']))
print("species_from_bt")
# print(testRN2.get_species_from_bt(bt('1100')))

# print("bt_from_reactions")
# print(testRN2.get_bt_from_reactions(['R8','R4','R2','R9']))
# print("reactions_from_bt")
# print(testRN2.get_reactions_from_bt(bt('00110011')))
# print("support_bt_from_reaction")
# print(testRN2.get_supp_bt_from_reaction('R1'),1)
# print("products_bt_from_reaction")
# print(testRN2.get_prod_bt_from_reaction('R4'))
# print("reactions_from_species")
# print(testRN2.get_reactions_from_species(['l']))

# print("get prod from reactions")
# print(testRN2.get_prod_from_reactions(['R3','R2']))
# print("get supp from reactions")
# print(testRN2.get_supp_from_reactions(['R3','R1','R7']))
file_path = '../../networks/autopoietic_ext2.txt'
testRN2 = load_pyCOT_from_file(file_path)
# print(testRN2.RnStr)

# print("get inflow")
# print(testRN2.get_inflow())

# print("get outflow")
# print(testRN2.get_outflow([]))
print("get connected species to species")
print(testRN2.get_connected_species_to_species([]))

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

# print("is_simple_connected_directly")
# print(testRN2.is_simple_connected_directly('d','s2'))