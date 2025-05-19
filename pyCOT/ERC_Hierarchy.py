#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 29 18:55:14 2023

@author: tveloz
"""
#sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import pandas as pd
import re
from bitarray import bitarray as bt
from bitarray import frozenbitarray as fbt
from scipy.optimize import linprog
import random as rm
import matplotlib.pyplot as plt
from pyvis.network import Network
import networkx as nx
import itertools
from itertools import combinations
from itertools import chain
from collections import defaultdict, deque
from collections import Counter # para el gráfico de cantidad de básicos vs repeticiones
from pyCOT.io.functions import *
from pyCOT.rn_visualize import *
from pyCOT.rn_rustworkx import *
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout

# Auxiliary functions
def species_list_to_names(species_list):
    """Convert list of Species objects or list of lists of Species objects to list of species names"""
    if not species_list:
        return []
    if isinstance(species_list[0], list):
        return [sorted([sp.name for sp in sublist]) for sublist in species_list]
    return sorted([sp.name for sp in species_list])

def get_sorted_species_names(species_set, RN):
    """Get sorted species names from a closure"""
    return sorted([sp.name for sp in closure(RN, species_set)])

def generators(RN):
    """Get list of generators from reactions' support"""
    gen = []
    for reaction in RN.reactions():
        support = RN.get_supp_from_reactions(reaction)
        if support:  # Only add non-empty supports
            gen.append(support)
    return gen

def closure(RN, X):
    """Get closure of species set X"""
    temp = X
    CL = list(set(temp).union(set(RN.get_prod_from_species(temp))))
    while RN.get_reactions_from_species(CL) != RN.get_reactions_from_species(temp):
        temp = CL
        CL = list(set(temp).union(set(RN.get_prod_from_species(temp))))
    return sorted(CL, key=lambda sp: sp.name)

def closures(RN, ListX):
    """Get closures for a list of species sets"""
    return [closure(RN, X) for X in ListX]

class ERC:
    def __init__(self, min_generators, label, all_generators=None):
        """Initialize an ERC object."""
        self.label = label
        self.min_generators = min_generators
        self.all_generators = all_generators

    def get_closure(self, RN):
        return closure(RN, self.min_generators[0])
    
    def get_reacs(self, RN):
        specs = self.get_closure(RN)
        reacs = RN.get_reactions_from_species(specs)
        return reacs

    @staticmethod
    def build_hierarchy_graph(ercs, RN):
        """Build the minimal containment graph between ERCs"""
        graph = nx.DiGraph()
        for erc in ercs:
            graph.add_node(erc.label, erc=erc)
            
        # Get all containment relations
        containments = []
        for erc1 in ercs:
            closure1 = set(species_list_to_names(erc1.get_closure(RN)))
            for erc2 in ercs:
                if erc1 != erc2:
                    closure2 = set(species_list_to_names(erc2.get_closure(RN)))
                    if closure2.issubset(closure1) and closure1 != closure2:
                        containments.append((erc1.label, erc2.label))
        
        # Keep only direct containments
        for start, end in containments:
            is_direct = True
            for erc in ercs:
                mid = erc.label
                if (start, mid) in containments and (mid, end) in containments:
                    is_direct = False
                    break
            if is_direct:
                graph.add_edge(start, end)
                
        return graph

    @staticmethod
    def get_node_levels(graph):
        """Calculate the true level of each node based on longest containment chain"""
        levels = {}
        for node in graph.nodes():
            leaf_nodes = [n for n in graph.nodes() if graph.out_degree(n) == 0]
            all_paths = []
            for leaf in leaf_nodes:
                paths = list(nx.all_simple_paths(graph, node, leaf))
                all_paths.extend(paths)
            max_path_length = max([len(path)-1 for path in all_paths]) if all_paths else 0
            levels[node] = max_path_length
        return levels

    @staticmethod
    def plot_hierarchy(ercs, RN, graph=None, figsize=(10,10)):
        """Plot the ERC hierarchy"""
        if graph is None:
            graph = ERC.build_hierarchy_graph(ercs, RN)
            
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)
        
        levels = ERC.get_node_levels(graph)
        pos = {}
        level_nodes = defaultdict(list)
        
        for node, level in levels.items():
            level_nodes[level].append(node)
        
        # Position nodes by level
        for level, nodes in level_nodes.items():
            n_nodes = len(nodes)
            nodes.sort(key=lambda n: len(nx.ancestors(graph, n)), reverse=True)
            for i, node in enumerate(nodes):
                x = (i - (n_nodes-1)/2) * 2.0
                y = level * 2.0
                pos[node] = (x, y)

        # Draw graph components
        nodes = nx.draw_networkx_nodes(graph, pos,
                                     node_color='lightblue',
                                     node_size=2000)
        nx.draw_networkx_edges(graph, pos,
                             edge_color='gray',
                             arrows=True,
                             arrowsize=20)
        nx.draw_networkx_labels(graph, pos,
                              font_size=12,
                              font_weight='bold')
        
        # Add hover functionality
        annot = ax.annotate("", 
                           xy=(0,0), xytext=(20,20),
                           textcoords="offset points",
                           bbox=dict(boxstyle="round", fc="w", ec="0.5", alpha=0.9),
                           ha='center',
                           visible=False)

        def update_annot(ind):
            node = list(graph.nodes())[ind["ind"][0]]
            erc = next((e for e in ercs if e.label == node), None)
            if RN and erc:
                closure_text = f"Closure: {species_list_to_names(erc.get_closure(RN))}"
            else:
                closure_text = "Hover data not available (RN not provided)"
            pos_node = pos[node]
            annot.xy = pos_node
            annot.set_text(closure_text)

        def hover(event):
            if event.inaxes == ax:
                cont, ind = nodes.contains(event)
                annot.set_visible(cont)
                if cont:
                    update_annot(ind)
                    fig.canvas.draw_idle()

        fig.canvas.mpl_connect("motion_notify_event", hover)
        
        plt.title("ERC Hierarchy")
        plt.axis('off')
        plt.show()
        return graph

    @staticmethod
    def get_all_synergies(base_erc, target_erc, erc_hierarchy, RN):
        """
        Find synergetic ERCs and their coverage details with base_erc.
        
        Returns
        -------
        dict
            Dictionary mapping synergetic ERCs to tuples of:
            - list of covered generators
            - coverage ratio (covered/total generators)
        """
        base_closure = set(species_list_to_names(base_erc.get_closure(RN)))
        synergy_details = {}
        total_generators = len(target_erc.min_generators)
        
        other_labels = [node for node in erc_hierarchy.nodes() 
                       if node not in [base_erc.label, target_erc.label]]
        
        # For each potential synergetic ERC
        for other_label in other_labels:
            other_erc = erc_hierarchy.nodes[other_label]['erc']
            other_closure = set(species_list_to_names(other_erc.get_closure(RN)))
            covered_gens = []
            
            for gen in target_erc.min_generators:
                gen_species = set(species_list_to_names(gen))
                combined_closure = base_closure.union(other_closure)
                
                # Check for true synergy with this generator
                if (gen_species.issubset(combined_closure) and
                    not gen_species.issubset(base_closure) and
                    not gen_species.issubset(other_closure) and
                    len(gen_species & base_closure - other_closure) > 0 and
                    len(gen_species & other_closure - base_closure) > 0):
                    covered_gens.append(gen)
            
            # If we found any synergies with this ERC, record the details
            if covered_gens:
                coverage_ratio = len(covered_gens) / total_generators
                synergy_details[other_erc] = (covered_gens, coverage_ratio)
                
        return synergy_details

    @staticmethod
    def get_modified_hierarchy(base_erc, ercs, RN):
        """Remove ERCs that contain or are contained by base_erc's closure"""
        graph = ERC.build_hierarchy_graph(ercs, RN)
        base_closure = set(species_list_to_names(base_erc.get_closure(RN)))
        
        nodes_to_remove = []
        for erc in ercs:
            if erc == base_erc:
                continue
            erc_closure = set(species_list_to_names(erc.get_closure(RN)))
            if erc_closure.issubset(base_closure) or base_closure.issubset(erc_closure):
                nodes_to_remove.append(erc.label)
                
        for node in nodes_to_remove:
            graph.remove_node(node)
            
        return graph

    @staticmethod
    def get_top_level_ercs(graph, ercs):
        """Get ERCs that have no incoming edges (not contained by any other)"""
        top_labels = [node for node in graph.nodes() if graph.in_degree(node) == 0]
        return [erc for erc in ercs if erc.label in top_labels]

    @staticmethod
    def ERCs(RN):
        """Generate ERCs from reaction network"""
        List_gen = generators(RN)
        remaining_gens = List_gen.copy()
        used_closures = set()
        ERC_list = []
        counter = 0
        
        while remaining_gens:
            gen_species = remaining_gens[0]
            closure_species = tuple(sorted([sp.name for sp in closure(RN, gen_species)]))
            
            if closure_species not in used_closures:
                all_gens = []
                i = 0
                while i < len(remaining_gens):
                    gen = remaining_gens[i]
                    gen_closure = tuple(sorted([sp.name for sp in closure(RN, gen)]))
                    if gen_closure == closure_species:
                        all_gens.append(gen)
                        remaining_gens.pop(i)
                    else:
                        i += 1
                
                min_gens = []
                for gen1 in all_gens:
                    is_minimal = True
                    gen1_set = set(sp.name for sp in gen1)
                    for gen2 in all_gens:
                        if gen1 != gen2:
                            gen2_set = set(sp.name for sp in gen2)
                            if gen2_set.issubset(gen1_set):
                                is_minimal = False
                                break
                    if is_minimal:
                        min_gens.append(gen1)
                
                erc = ERC(min_generators=min_gens,
                         label=f"E{counter}",
                         all_generators=all_gens)
                ERC_list.append(erc)
                used_closures.add(closure_species)
                counter += 1
            else:
                remaining_gens.pop(0)
        
        return ERC_list

__all__ = [
    'ERC',
    'generators',
    'species_list_to_names',
    'closure',
    'closures'
]



