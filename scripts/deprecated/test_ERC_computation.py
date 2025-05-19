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

import time

import networkx as nx
from pyCOT.file_manipulation import *
import matplotlib.pyplot as plt

from pyCOT.reaction_network import *
from pyCOT.closure_structure import *
from networkx.drawing.nx_agraph import graphviz_layout

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
file_path = 'Txt/autopoietic.txt' 
# file_path = '../../networks/testing/in_out.txt'
# file_path = '../../networks/testing/RedPDoSR00.txt'
#file_path = '../../networks/biomodels_interesting/BIOMD0000000652_manyOrgs.xml'
# file_path = '../../networks/biomodels_interesting/Biomodel_724.txt'
#file_path = '../../networks/testing/ERCs_test.txt'
#file_path = '../../networks/testing/Synergy_test.txt'
# file_path = 'networks/Navarino/RN_IN_05.txt'
#file_path = '../../networks/testing/MSORN_test1.txt'
#file_path = '../../networks/testing/MSORN_test2.txt'
#file_path = 'networks/RandomAlife/RN_Ns_40_Norg_20_id_396.txt'
# file_path = 'networks/biomodels_interesting/central_ecoli.txt'
# file_path = 'networks/biomodels_all_txt/BIOMD0000001044.txt'
#file_path = 'networks/biomodels_all_txt/BIOMD0000001002.txt'
#file_path = 'networks/GPT_Generated/EColi.txt'
#file_path = '../../networks/biomodels_interesting/BIOMD0000000237_manyOrgs.xml'
testRN2 = load_pyCOT_from_file(file_path)

# reaction_graph = build_graph(testRN2)
# visualize_graph(reaction_graph)
print("specs")
print(testRN2.SpStr)

# print(testRN2.SpStr)
# print("specsBt")
# print(testRN2.SpBt)
# print("reacsBt")
# print(testRN2.RnBt)
# print("Reac")
# print(testRN2.RnStr)
# print("MSupp")
# print(testRN2.RnMsupp)
# print("MProd")
# print(testRN2.RnMprod)
#reaction_graph = build_graph(testRN2)
#visualize_graph(testRN2)

erc=ERCs(testRN2)
print("ERCs created")
print("the list of "+str(len(erc))+" ERCs")
for g in erc:
      print(g)
print("#Reactions: "+str(len(testRN2.RnStr))+", #ERCs: "+str(len(erc)))

con=get_containments(testRN2,erc)
print("check containment")
for g in con:
      print(g)
print("the list of "+str(len(con))+" containments")

dc=get_direct_containments(testRN2,erc,con)
print("check direct containment")
for g in dc:
      print(g)
print("the list of "+str(len(dc))+" direct containments")
print("Number of species ="+str(len(testRN2.SpStr)))
print("Number of reactions ="+str(len(testRN2.RnStr)))
print("Number of ERCs ="+str(len(erc)))


#Plotting basics of ERC structure
erc_sizes = [len(entry[1]) for entry in erc]
reaction_counts = [len(entry[2]) for entry in erc]
erc_names = [entry[3] for entry in erc]
# Step 2: Sort the data
sorted_erc_sizes = sorted(erc_sizes)
sorted_reaction_counts = sorted(reaction_counts)

# Step 3: Plot the sorted reaction counts and ERC sizes
plt.figure(figsize=(15, 5))

# # Plot 1: Reaction counts sorted by size
# plt.subplot(1, 3, 1)
# plt.plot(sorted_reaction_counts, marker='o')
# plt.title("Reaction Counts (Sorted)")
# plt.xlabel("Index")
# plt.ylabel("Number of Reactions")

# # Plot 2: ERC sizes sorted by size
# plt.subplot(1, 3, 2)
# plt.plot(sorted_erc_sizes, marker='o')
# plt.title("ERC Sizes (Sorted)")
# plt.xlabel("Index")
# plt.ylabel("Size of ERC (Number of Species)")

# Plot 3: ERC size vs. Reaction count scatter plot
plt.subplot(1, 3, 1)
plt.scatter(erc_sizes, reaction_counts, alpha=0.7)
plt.title("ERC Size vs Number of Reactions")
plt.xlabel("Size of ERC (Number of Species)")
plt.ylabel("Number of Reactions")

# Plot 4: Histograms of both ERC sizes and reaction counts
plt.subplot(1, 3, 2)
plt.hist(erc_sizes, alpha=0.5, label='ERC Sizes')
plt.hist(reaction_counts, alpha=0.5, label='Reaction Counts')
plt.legend(loc='upper right')
plt.title("Histograms of ERC Sizes and Reaction Counts")
plt.xlabel("Size")
plt.ylabel("Frequency")

plt.show()

