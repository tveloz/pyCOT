#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 15:28:44 2023

@author: tveloz
"""

import sys
import os

# Añadir el directorio raíz al PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pyCOT.reaction_network import ReactionNetwork
from pyCOT.closure_structure import *
from  pyCOT.reactive_features import *
import networkx as nx
from pyCOT.file_manipulation import *
import matplotlib.pyplot as plt
import time
from collections import defaultdict
from collections import Counter # para el gráfico de cantidad de básicos vs repeticiones
from pathlib import Path
import pickle
import glob
from pathlib import Path

#from Closure_Str_RN_Computing import Closure_Str_RN_Computing

# This script aims at performing the analysis made in Closure_Str_RN_Computing.py in a systematic way. It reads a folder of
# RN in .txt format. The script can be executed in multiple runs as it stores everything in .pkl files, one for each
# .txt RN and performs a check to identify where should start to continue computing the missing parts.

#folder_path = 'networks/Navarino'
#txt_files = glob.glob(os.path.join(folder_path, '*.txt'))
input_path = Path("networks/Navarino/RN_IN_02_PyCOT.txt")
#input_path = Path('networks/biomodels_all_txt/BIOMD000000151.txt')
input_path = Path("networks/biomodels_interesting/bigg_iAF692.txt")
#input_path = Path("networks/biomodels_all_txt/BIOMD0000000183.txt")
#input_path = Path('networks/testing/autopoietic_ext.txt')
                  
pkl_file = input_path.with_name(input_path.stem + '.pkl')   
# Open the pickle file in read-binary mode
objects=[]
print("checking if "+str(pkl_file)+" exists")

if pkl_file.exists():
    print(str(pkl_file)+" found")
    # Load RN
    RN = load_pyCOT_from_file(input_path)
    if len(RN.RnStr)==0:
        print("zero reaction system")
    else:
        # Append basic information about the RN
        with open(pkl_file, 'rb') as file:
            while True:
                try:
                    # Attempt to load the next object
                    obj = pickle.load(file)
                    objects.append(obj)
                except EOFError:
                    # End of file reached, break out of the loop
                    break
        print("objeto de largo :"+ str(len(objects))+" retrieved")
        objects.append({"Species":[RN.SpStr]})
        objects.append({"Reactions":[RN.RnStr]})
        # for obj in objects:
        #     print(obj.keys())
else:
     print("file not found")

############ Identify variables retrieved #####################
ERCs=objects[0]["ERCs"]
ERCs_ct=objects[1]['ERCs exec time']
Containments=objects[2]['Containments']
Containments_ct=objects[3]['Containments exec time']
DirectContainments=objects[4]['Direct Containments']
DirectContainments_ct=objects[5]['Direct Containments exec time']
Synergies=objects[6]['Synergies']
Synergies_ct=objects[7]['Synergies exec time']
DirectSynergies=objects[8]['Direct Synergies']
DirectSynergies_ct=objects[9]['Direct Synergies exec time']
Nspecs=len(objects[10]['Species'][0])
Nreacs=len(objects[11]['Reactions'][0])
print("Species: "+str(len(RN.SpStr))+", Reactions: "+str(len(RN.RnStr)))
print(RN.SpStr)
############ ERCs structure in terms of species and reactions ############
ERCStr=[]
# store number of minimal generators per ERC
MinGensize=[len(erc[0]) for erc in ERCs]
# number species per ERC/number of species
ERCStr.append(MinGensize)
#proportion of species per ERC
PropSizeERC=[round(len(erc[1])/Nspecs,3) for erc in ERCs]
ERCStr.append(PropSizeERC)
# proportion of reactions in the equivalence class
PropReacsERC=[round(len(erc[2])/Nreacs,3) for erc in ERCs]
ERCStr.append(PropReacsERC)
# proportion of ERCs in which each species is part of
ERCclosures=[erc[1] for erc in ERCs] 
Specs_prop= []
for sp in objects[10]['Species'][0]:
    Specs_prop.append([sp,sum(sublist.count(sp) for sublist in ERCclosures)])
ERCStr.append(Specs_prop)
#print("ERCs")
# for el in ERCStr:
#     print(el)

########## Chain characterization in terms of levels #############
def build_chains(set_labels, containment_pairs):
    # Step 1: Build the adjacency list (graph)
    containment_graph = defaultdict(list)
    reverse_graph = defaultdict(list)
    
    for a, b in containment_pairs:
        containment_graph[a].append(b)
        reverse_graph[b].append(a)

    # Step 2: Identify root nodes (nodes that are not children of any other node)
    root_nodes = set(set_labels) - set(reverse_graph.keys())

    # Step 3: Function to recursively find chains via DFS
    def find_chains(node, current_chain, all_chains):
        current_chain.append(node)
        # If node has no children, it is the end of a chain
        if node not in containment_graph or not containment_graph[node]:
            all_chains.append(current_chain[:])  # Add a copy of the current chain
        else:
            # Recursively explore each child node
            for child in containment_graph[node]:
                find_chains(child, current_chain, all_chains)
        # Backtrack
        current_chain.pop()

    # Step 4: Find all chains starting from each root node
    all_chains = []
    for root in root_nodes:
        find_chains(root, [], all_chains)

    return all_chains

# Example usage

set_labels = [sublist[3] for sublist in ERCs]
containment_pairs = DirectContainments

chains = build_chains(set_labels, containment_pairs)
# print("ERCs")
# print(set_labels)
# print("Direct Containments")
# print(containment_pairs)
# print("Direct Containments resorted as chains")
# print(chains)

# Step 4: Calculate chain lengths
chain_lengths = [len(chain) for chain in chains]
#print("chain lenghts:")
#print(chain_lengths)
############## Level analysis ##################
# plt.hist(chain_lengths, bins=5, edgecolor='black')

############## Level analysis ##################


plt.hist(chain_lengths, bins=5, edgecolor='black')


# Add titles and labels
plt.title('Histogram of Chain Lengths')
plt.xlabel('Length')
plt.ylabel('Frequency')

# Show the plot
plt.show()
# Add titles and labels

# Step 1: Find the length of the longest chain
max_chain_length = max(len(chain) for chain in chains)

# Step 2: Initialize an empty list of lists to store chains by their level
chains_per_level = [[] for _ in range(max_chain_length)]

# Step 3: Sort the chains into the appropriate level based on their size
for chain in chains:
    chain_length = len(chain)
    # Store the chain in the corresponding level (subtract 1 to match 0-indexing)
    chains_per_level[chain_length - 1].append(chain)

# chains_per_level[0] will have chains of size 1, chains_per_level[1] will have chains of size 2, and so on
#for level in chains_per_level:
    #print(level)
# print("$$$$$$$$$$$$$$$$$$$$$$$$$$")
lengths_per_level=[]
closed_per_level=[]
reqs_per_level=[]
for level in range(0,max_chain_length):
    lengths_per_level.append([])
    closed_per_level.append([])
    reqs_per_level.append([])

    #print("level "+str(level))
    if len(chains_per_level[level])==0:
        lengths_per_level[level].append(0)
    else:
        for erc_indexes in chains_per_level[level]:
            head_erc=ERCs[erc_indexes[0]][1]
            #print(head_erc)
            closed_per_level[level].append(head_erc)
            lengths_per_level[level].append(len(head_erc))
#print("lengths per level")     
#print(lengths_per_level)

# Step 2: Calculate the averages and standard deviations
averages = [np.mean(lst) for lst in lengths_per_level]  # Mean (average)
std_devs = [np.std(lst) for lst in lengths_per_level]   # Standard deviation

# Step 3: Plot the averages with error bars for standard deviation
plt.errorbar(range(len(averages)), averages, yerr=std_devs, fmt='o', capsize=5, linestyle='-', color='b', label='Average with Std Dev')

# Step 4: Add labels and title
plt.title('Averages of Lists with Standard Deviation')
plt.xlabel('List index')
plt.ylabel('Average value')

# Step 5: Display the plot
plt.legend()
plt.show()
ERCs_indexes=[el[3] for el in ERCs]
chains_of_ERCs=[]
for erc_index in ERCs_indexes:
    chains_of_erc= [x for x in chains if erc_index in x]
    chains_of_ERCs.append([erc_index,chains_of_erc])
for i in range(len(chains_of_ERCs)):
    print(chains_of_ERCs[i])
print('###############')
ERCs_level=[]
#ERCs_exp_lenght=[]
for ch in chains_of_ERCs:
    erc=ch[0]
    chs=ch[1]    
    lev=max(len(x)-x.index(erc) for x in chs)
    ERCs_level.append([erc,lev])
print(ERCs_level)

plt.hist([x[1] for x in ERCs_level], bins=5, edgecolor='black')
# Add titles and labels
plt.title('Histogram of ERC Levels')
plt.xlabel('Level')
plt.ylabel('Frequency')
plt.show()