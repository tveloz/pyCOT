#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 15:28:44 2023

@author: tveloz
"""

import sys
import os

# Add root directory to PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pyCOT.rn_rustworkx import ReactionNetwork
from pyCOT.ERC_Hierarchy import (
    ERC, species_list_to_names, generators
)
from pyCOT.io.functions import read_txt
import networkx as nx
import matplotlib.pyplot as plt
from collections import defaultdict
from collections import Counter # para el gráfico de cantidad de básicos vs repeticiones

file_path = 'networks/testing/autopoietic.txt'
file_path = 'networks/Navarino/RN_IN_05.txt'
#file_path = 'networks/testing/Farm.txt'
#file_path = 'networks/testing/autopoietic_ext.txt'
#file_path = 'networks/testing/Synergy_test.txt'
#file_path = 'networks/RandomAlife/RN_Ns_40_Norg_20_id_396.txt'
#file_path = 'networks/RandomAlife/RN_Ns_25_Norg_17_id_174.txt'
#file_path = 'networks/biomodels_interesting/central_ecoli.txt'


# Assuming load_pyCOT_from_file and reac_analysis are defined elsewhere in your code

RN = read_txt(file_path)

print("Reactions in the network:")
print("------------------------")


reactions = [reaction.name() for reaction in RN.reactions()]
print(reactions)

print("Species in the network:")
print("------------------------")
species = [specie.name for specie in RN.species()]
print(species)

gen = generators(RN)
print("Generators in the network:")
for g in gen:
    print(species_list_to_names(g))
ercs = ERC.ERCs(RN)  # Use ERC.ERCs instead of ERCs
print("ERCs in the network:")
for e in ercs:
    print(f"\nERC {e.label}:")
    print(f"Minimal generators: {species_list_to_names(e.min_generators)}")
    print(f"Closure: {species_list_to_names(e.get_closure(RN))}")
    print(f"All generators: {species_list_to_names(e.all_generators)}")
    print("------------------------")
 

# Create ERCs
ercs = ercs  # No need to recreate ERCs, they're already proper objects

# Test partially intersecting reactions for first ERC
# Test partially intersecting reactions for all ERCs
for test_erc in ercs:
    print("\nTesting partially intersecting reactions for ERC:", test_erc.label)
    print("Closure:", species_list_to_names(test_erc.get_closure(RN)), "min_generators:", species_list_to_names(test_erc.min_generators))
    print("------------------------")

    partial_reactions = RN.get_reactions_partially_intersecting_support(test_erc.get_closure(RN))
    if partial_reactions:
        print("\nPartially intersecting reactions:")
        for reaction in partial_reactions:
            support_species = species_list_to_names(RN.get_supp_from_reactions(reaction))
            intersection = set(support_species) & set(species_list_to_names(test_erc.get_closure(RN)))
            print(f"\n{reaction.name()}:")
            print(f"  Full support: {support_species}")
            print(f"  Intersection with closure: {list(intersection)}")
    else:
        print("No partially intersecting reactions found.")
    print("------------------------")

# Get original and modified graphs
original_graph = ERC.build_hierarchy_graph(ercs, RN)
modified_graph = ERC.get_modified_hierarchy(test_erc, ercs, RN)

print("\nOriginal hierarchy edges:")
print(original_graph.edges())
print("\nModified hierarchy edges:")
print(modified_graph.edges())
print("\nRemoved nodes:")
removed = set(original_graph.nodes()) - set(modified_graph.nodes())
print(removed)

# Create hierarchy graph first
hierarchy_graph = ERC.build_hierarchy_graph(ercs, RN)

# Test get_all_synergies between pairs of ERCs
print("\nSynergy Analysis Report:")
print("=======================")

# Statistics collectors
total_pairs = 0  # Total possible base-target pairs
total_synergetic_pairs = 0  # Pairs that have at least one synergy
synergies_per_base = defaultdict(int)  # Count of synergetic ERCs found for each base
synergies_per_target = defaultdict(int)  # Count of synergetic ERCs found for each target
coverage_ratios = []
total_generators_covered = 0

for base_erc in ercs:
    for target_erc in ercs:
        if base_erc != target_erc:
            total_pairs += 1
            synergy_details = ERC.get_all_synergies(base_erc, target_erc, hierarchy_graph, RN)
            if synergy_details:
                total_synergetic_pairs += 1
                # Count each synergetic ERC found
                synergies_per_base[base_erc.label] += 1
                synergies_per_target[target_erc.label] += 1
                
                # Collect coverage information
                total_target_generators = len(target_erc.min_generators)
                unique_covered_gens = set()
                for _, (covered_gens, _) in synergy_details.items():
                    # Convert generators to hashable tuples of species names
                    gen_tuples = {tuple(species_list_to_names(gen)) for gen in covered_gens}
                    unique_covered_gens.update(gen_tuples)
                coverage_ratio = len(unique_covered_gens) / total_target_generators
                coverage_ratios.append(coverage_ratio)
                total_generators_covered += len(unique_covered_gens)

# Calculate statistics
possible_targets = len(ercs) - 1  # Each ERC can target all others
avg_synergies_per_base = sum(synergies_per_base.values()) / len(ercs)
max_synergies_base = max(synergies_per_base.items(), key=lambda x: x[1]) if synergies_per_base else ('None', 0)
synergy_probability = total_synergetic_pairs / total_pairs if total_pairs > 0 else 0
avg_coverage = sum(coverage_ratios) / len(coverage_ratios) if coverage_ratios else 0

# Print report
print(f"\nOverall Statistics:")
print(f"- Total possible ERC pairs: {total_pairs}")
print(f"- Pairs with synergies: {total_synergetic_pairs}")
print(f"- Probability of finding synergy: {synergy_probability:.2%}")
print(f"- Average synergetic pairs per base ERC: {avg_synergies_per_base:.2f} out of {possible_targets} possible")
print(f"- Most synergetic base ERC: {max_synergies_base[0]} with {max_synergies_base[1]} synergetic pairs")
print(f"- Average generator coverage when synergy exists: {avg_coverage:.2%}")
# print("\nDetailed Base ERC Analysis:")
# print("-------------------------")
# for erc_label, count in sorted(synergies_per_base.items(), key=lambda x: x[1], reverse=True):
#     percentage = count / possible_targets
#     print(f"ERC {erc_label}: {count} synergetic pairs ({percentage:.2%} of possible targets)")

