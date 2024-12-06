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
file_path = '../../networks/biomodels_interesting/Biomodel_724.txt'
#file_path = '../../networks/testing/ERCs_test.txt'
#file_path = '../../networks/testing/Synergy_test.txt'
file_path = 'networks/testing/Farm.txt'
#file_path = '../../networks/testing/MSORN_test1.txt'
#file_path = '../../networks/testing/MSORN_test2.txt'
#file_path = 'networks/RandomAlife/RN_Ns_40_Norg_20_id_396.txt'
# file_path = 'networks/biomodels_interesting/central_ecoli.txt'
# file_path = 'networks/biomodels_all_txt/BIOMD0000001044.txt'
# file_path = 'networks/biomodels_all_txt/BIOMD0000001003.txt'
# file_path = 'networks/GPT_Generated/EColi.txt   '
#file_path = '../../networks/biomodels_interesting/BIOMD0000000237_manyOrgs.xml'
testRN2 = load_pyCOT_from_file(file_path)

# reaction_graph = build_graph(testRN2)
# visualize_graph(reaction_graph)
print(testRN2.SpStr)
print("specs")
print(testRN2.SpStr)
print("specsBt")
print(testRN2.SpBt)
print("reacsBt")
print(testRN2.RnBt)
print("Reac")
print(testRN2.RnStr)
print("MSupp")
print(testRN2.RnMsupp)
print("MProd")
print(testRN2.RnMprod)
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



# Example list of containment relationships (your ERCs)
containment_list = dc

# Step 1: Initialize a directed graph (forest)
G = nx.DiGraph()

# Step 2: Add edges based on containment list
for containment in containment_list:
    G.add_edge(containment[0], containment[1])

# Step 3: Visualize the forest
# Using a layout that is good for hierarchical structures
pos = nx.spring_layout(G, seed=42)  # You can also try nx.kamada_kawai_layout or nx.spectral_layout

# Draw the graph with labels
plt.figure(figsize=(8, 6))
nx.draw(G, pos, with_labels=True, node_size=2000, node_color="skyblue", font_size=10, font_weight="bold", arrows=True)
plt.title("ERC Containment Hierarchy (Forest)")
plt.show()


# # print("check synergies")
# syn=get_synergies(testRN2,erc)
# for g in syn:
#       print(g)
# print("the list of "+str(len(syn))+" synergies")

# d_syn=get_direct_synergies(testRN2,erc,dc,syn)

# for g in d_syn:
#       print(g)
# print("the list of "+str(len(d_syn))+" direct synergies")


# sorn=second_order_network(testRN2,erc,dc,d_syn)
# print("The list of second order network species")
# for i in range(len(sorn[0][0])):
#       print(str(sorn[0][0][i])+" has lenght "+str(len(sorn[0][1][i])))
            
# print("the list of "+str(len(sorn[1]))+" second order reaction network")
# for i in range(len(sorn[1])):    
#     print("r"+str(i)+": "+str(sorn[1][i]))

# msorn=MSORN(sorn)
# print(str(msorn.get_connected_species_to_species("E4")))
# start_time = time.time()
# cl_pth=Closed_paths_backward(msorn)
# end_time = time.time()
# execution_time = end_time - start_time
# print("Execution time computing closed paths", execution_time, "seconds")
# print("number of closed paths ="+str(len(cl_pth)))
# # for g in cl_pth:
# #     print(g)    
# gen=Closed_paths_to_Sp(cl_pth,sorn)
# print("************ANALYSIS RESULTS******************")

# print("We have "+str(len(gen))+" closed sets")
# gen_original=gen.copy()
# gen=remove_duplicates(gen)
# print("After eliminating duplicates we have "+str(len(gen))+" closed sets")
# for g in gen:
#       print(testRN2.get_reactions_from_species(g))
# SSM=0
# sizeSSM=0
# sizenSSM=0
# for g in gen:
#     if testRN2.is_semi_self_maintaining(g):
#         print("Closed set of size of size "+str(len(g))+ " is SSM")
#         SSM=SSM+1
#         sizeSSM=sizeSSM+len(g)
#     else:
#         print("Closed set of size "+str(len(g))+ " is not SSM")
#         sizenSSM=sizenSSM+len(g)
# nSSM=len(gen)-SSM
# print("Number of species ="+str(len(testRN2.SpStr)))
# print("Number of reactions ="+str(len(testRN2.RnStr)))
# print("Number of ERCs ="+str(len(msorn.SpStr)))
# print("SSM = "+str(SSM)+", not SSM ="+str(nSSM))
# print("relative size SSM ="+str(sizeSSM)+"/"+str(SSM))
# print("relative size nSSM ="+str(sizenSSM)+"/"+str(nSSM))



