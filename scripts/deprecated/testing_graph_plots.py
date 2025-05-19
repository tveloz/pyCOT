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
#file_path = '../../networks/biomodels_interesting/BIOMD0000000652_manyOrgs.xml'
#file_path = '../../networks/biomodels_interesting/Biomodel_724.txt'
#file_path = '../../networks/testing/ERCs_test.txt'
#file_path = '../../networks/testing/Synergy_test.txt'
file_path = '../../networks/testing/Farm.txt'
#file_path = '../../networks/testing/MSORN_test1.txt'
#file_path = '../../networks/testing/MSORN_test2.txt'
#file_path = '../../networks/RandomAlife/RN_Ns_40_Norg_24_id_316.txt'
#file_path = '../../networks/biomodels_interesting/central_ecoli.txt'
#file_path = '../../networks/biomodels_all/BIOMD0000000086_url.xml'
#file_path = '../../networks/biomodels_interesting/BIOMD0000000237_manyOrgs.xml'
testRN2 = load_pyCOT_from_file(file_path)
#testRN2 = load_pyCOT_from_Sbml(file_path)

G=testRN2.pyCOT_to_Graph()
species_nodes = [node for node, data in G.nodes(data=True) if data['type'] == 'species']
reaction_nodes = [node for node, data in G.nodes(data=True) if data['type'] == 'reaction']
pos = nx.bipartite_layout(G, species_nodes)
nx.draw(G, pos, with_labels=True, node_size=200, node_color='lightblue', font_size=10)
plt.show()
species= [node for node, attr in G.nodes(data=True) if attr.get('bipartite') == 1]
print(species[0])
cl_pth=find_closed_sets(G,species[0])
print(cl_pth)