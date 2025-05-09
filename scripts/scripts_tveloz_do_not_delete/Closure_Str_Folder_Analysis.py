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
import pandas as pd
from collections import defaultdict
from collections import Counter # para el gráfico de cantidad de básicos vs repeticiones
from pathlib import Path
import pickle
import glob
from pathlib import Path
from Closure_Str_RN_Computing import Closure_Str_RN_Computing 


# Function to process each .pkl file and extract relevant data
def process_pkl_file(pkl_file):
    with open(pkl_file, 'rb') as file:
        objects = []
        while True:
            try:
                # Attempt to load the next object
                obj = pickle.load(file)
                objects.append(obj)
            except EOFError:
                break

    # Extract necessary data
    #RN = objects[0]['RN']  # Assuming 'RN' is the reaction network object
    ERCs = objects[0]['ERCs']
     # Loop through each sublist and add elements from the third sublist to the set
    unique_species = set()
    unique_reactions = set()
    for sublist in ERCs:
        unique_species.update(sublist[1])
        unique_reactions.update(sublist[2])

    # Convert the set back to a list (optional, if you need a list)
    Species = list(unique_species)
    Reactions = list(unique_reactions)
    if len(Reactions)>3:
        Containments = objects[2]['Containments']
        DirectContainments = objects[4]['Direct Containments']
        Synergies = objects[6]['Synergies']
        DirectSynergies = objects[8]['Direct Synergies']
        
        # Chain construction
        set_labels = [sublist[3] for sublist in ERCs]
        containment_pairs = DirectContainments
        chains = build_chains(set_labels, containment_pairs)

        # Number of chains per level
        max_chain_length = max(len(chain) for chain in chains)
        chains_per_level = [[] for _ in range(max_chain_length)]
        for chain in chains:
            chains_per_level[len(chain) - 1].append(chain)

        # Prepare data for DataFrame
        data = {
            'RN_id': pkl_file.stem,
            'Specs': len(Species),
            'Reacs': len(Reactions),
            'ERCs': ERCs,
            '#Containments': len(Containments),
            'DirectContainments': DirectContainments,
            '#Synergies': len(Synergies),
            'DirectSynergies': DirectSynergies,
            'chains_per_level': chains_per_level
        }
        
        return data
    else:
        data = {
            'RN_id': pkl_file.stem,
            'Specs': len(Species),
            'Reacs': len(Reactions),
            'ERCs': [],
            '#Containments': 0,
            'DirectContainments': [],
            '#Synergies': 0,
            'DirectSynergies': [],
            'chains_per_level': []
        }
        return data
    # Main function to process all .pkl files in a folder
def process_folder(folder_path):
    all_data = []
    for pkl_file in Path(folder_path).glob("*.pkl"):
        data = process_pkl_file(pkl_file)
        all_data.append(data)
    
    # Create DataFrame from all processed data
    df = pd.DataFrame(all_data)
    return df

# Example usage

folder_path = 'networks/biomodels_all_txt'
#txt_files = glob.glob(os.path.join(folder_path, '*.txt'))
df = process_folder(folder_path)

# Save the DataFrame to a CSV or do further analysis
df.to_csv("reaction_networks_analysis.csv", index=False)
print(df)   
# Handle cases where ERCs could be None or empty list
empty_erc_count = df[df['ERCs'].apply(lambda x: x is None or len(x) == 0)].shape[0]

# Print the result
print(f"Number of elements with an empty or None ERCs dictionary: {empty_erc_count}")
# Define cutoff values
min_specs = 0  # Set your minimum number of species
max_specs = 100  # Set your maximum number of species
min_reacs = 3    # Set your minimum number of reactions
max_reacs = 100  # Set your maximum number of reactions

# Filter the DataFrame based on the cutoff values
filtered_df = df[
    (df['Specs'] >= min_specs) & 
    (df['Specs'] <= max_specs) & 
    (df['Reacs'] >= min_reacs) & 
    (df['Reacs'] <= max_reacs)
]

# Plot Specs vs Reacs using the filtered DataFrame
# plt.figure(figsize=(8, 6))
# plt.scatter(filtered_df['Specs'], filtered_df['Reacs'], c='blue', label='Filtered Specs vs Reacs')

# # Add labels and title
# plt.xlabel('Number of Species (Specs)')
# plt.ylabel('Number of Reactions (Reacs)')
# plt.title(f'Scatter Plot of Specs vs Reacs\n(Cutoff: {min_specs} ≤ Specs ≤ {max_specs}, {min_reacs} ≤ Reacs ≤ {max_reacs})')

# Show the plot with a legend
# plt.legend()
# plt.grid(True)
# plt.show()

#Calculating number of ERCs per species

# If 'ERCs' contains lists or sets, calculate the number of ERCs per row
filtered_df['ERCs_count'] = filtered_df['ERCs'].apply(lambda x: len(x))
filtered_df['Direct_Containments_count'] = filtered_df['DirectContainments'].apply(lambda x: len(x))
filtered_df['DirectSynergies_count'] = filtered_df['DirectSynergies'].apply(lambda x: len(x))
#filtered_df['chain_lengths'] = filtered_df['chains_per_level'].apply(lambda x: [len(chain) for chain in x] if x else [])
# Now plot Specs (number of species) vs ERCs_count (number of ERCs)
plt.figure(figsize=(8, 6))
plt.scatter(filtered_df['Reacs'],filtered_df['DirectSynergies_count'],  c='blue', label='ERCs vs Synergies')
# plt.scatter(filtered_df['Reacs'],filtered_df['Direct_Containments_count'] , c='green', label='Direct Containments vs Containments')
#plt.scatter(filtered_df['ERCs_count'],filtered_df['chain_lengths'] , c='blue', label='Direct Containments vs Containments')
# plt.scatter(filtered_df['Reacs'],filtered_df['ERCs_count'] , c='yellow', label='Direct Containments vs Containments')
# Add labels and titlex
# plt.xlabel('Number of Species (Specs)')
# plt.ylabel('Number of ERCs')
# plt.title('Scatter Plot of Number of Species vs Number of ERCs')
plt.scatter(filtered_df['Direct_Containments_count'],filtered_df['Direct_Synergies_count'] , c='green', label='plot')
# plt.scatter(filtered_df['ERCs_count'],filtered_df['Direct_Containments_count'] , c='blue', label='Direct Containments vs Containments')
# plt.scatter(filtered_df['Reacs'],filtered_df['ERCs_count'] , c='yellow', label='Direct Containments vs Containments')
# Show the plot with a legend
plt.legend()
plt.grid(True)
plt.show()