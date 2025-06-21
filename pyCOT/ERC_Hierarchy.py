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
import plotly.graph_objects as go
from plotly.subplots import make_subplots

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
        #print(f"Processing reaction: {reaction.name()}"+ " number of generatirs: "+str(len(gen)))
        support = RN.get_supp_from_reactions(reaction)
        if support:  # Only add non-empty supports
            gen.append(support)
    return gen

def closure(RN, X):
    """Get closure of species set X"""
    temp = X
    CL = list(set(temp).union(set(RN.get_prod_from_species(temp))))
    reacstemp=[reaction.name() for reaction in RN.get_reactions_from_species(temp)]
    reacsCL=[reaction.name() for reaction in RN.get_reactions_from_species(CL)]
    while set(reacstemp) != set(reacsCL):
        #print("is "+str(reacsCL) + " == "+ str(reacstemp))
        temp = CL
        CL = list(set(temp).union(set(RN.get_prod_from_species(temp))))
        reacstemp=[reaction.name() for reaction in RN.get_reactions_from_species(temp)]
        reacsCL=[reaction.name() for reaction in RN.get_reactions_from_species(CL)]
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
    def find_minimal_generators(generators):
        minimal = []
        gen_sets = [(gen, set(sp.name for sp in gen)) for gen in generators]

        for i, (gen1, set1) in enumerate(gen_sets):
            is_minimal = True
            for j, (gen2, set2) in enumerate(gen_sets):
                if i != j and set2.issubset(set1) and set2 != set1:
                    is_minimal = False
                    break
            if is_minimal:
                minimal.append(gen1)
        return minimal
    
    @staticmethod
    def ERCs(RN):
        # 1. Compute all closures once (O(n))
        generators_list = generators(RN)
        closure_cache = {id(gen): closure(RN, gen) for gen in generators_list}
        
        # 2. Group by closure signature (O(n))
        closure_groups = defaultdict(list)
        for gen in generators_list:
            sig = tuple(sorted(sp.name for sp in closure_cache[id(gen)]))
            closure_groups[sig].append(gen)
        
        # 3. Find minimal generators for each group (O(Σk²))
        ERC_list = []
        for counter, (closure_sig, gens_group) in enumerate(closure_groups.items()):
            min_gens = ERC.find_minimal_generators(gens_group)
            
            erc = ERC(min_generators=min_gens,
                    label=f"E{counter}",
                    all_generators=gens_group)
            ERC_list.append(erc)
        
        return ERC_list
    
    @staticmethod
    def build_hierarchy_graph(ercs, RN):
        """Build the minimal containment graph between ERCs"""
        #print("Building hierarchy graph...")
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
    def build_hierarchy_graph2(ercs, RN):
        """Build the minimal containment graph between ERCs - optimized version"""
        graph = nx.DiGraph()
        
        # Add all nodes
        for erc in ercs:
            graph.add_node(erc.label, erc=erc)
        
        if len(ercs) <= 1:
            return graph
        
        # 1. Compute all closures once and cache as sets
        print(f"Computing closures for {len(ercs)} ERCs...")
        closure_cache = {}
        for erc in ercs:
            closure_cache[erc.label] = set(species_list_to_names(erc.get_closure(RN)))
        
        # 2. Build all containment relations efficiently
        print("Finding containment relations...")
        containments = []
        erc_labels = [erc.label for erc in ercs]
        
        for i, erc1 in enumerate(ercs):
            closure1 = closure_cache[erc1.label]
            for j, erc2 in enumerate(ercs):
                if i != j:
                    closure2 = closure_cache[erc2.label]
                    if closure2.issubset(closure1) and closure1 != closure2:
                        containments.append((erc1.label, erc2.label))
        
        print(f"Found {len(containments)} containment relations")
        
        # 3. Use NetworkX transitive reduction for optimal direct containment detection
        try:
            # Build temporary graph with all containments
            temp_graph = nx.DiGraph()
            temp_graph.add_nodes_from(erc_labels)
            temp_graph.add_edges_from(containments)
            
            # Get transitive reduction (removes indirect edges)
            direct_graph = nx.transitive_reduction(temp_graph)
            
            # Add direct edges to final graph
            for edge in direct_graph.edges():
                graph.add_edge(edge[0], edge[1])
                
            print(f"Reduced to {len(direct_graph.edges())} direct containments")
            
        except AttributeError:
            # Fallback: manual direct containment detection (optimized)
            print("Using fallback method for direct containment detection...")
            containment_set = set(containments)
            
            for start, end in containments:
                is_direct = True
                # Check if there's any intermediate ERC
                for mid_label in erc_labels:
                    if mid_label != start and mid_label != end:
                        if (start, mid_label) in containment_set and (mid_label, end) in containment_set:
                            is_direct = False
                            break
                if is_direct:
                    graph.add_edge(start, end)
                    
            print(f"Reduced to {len(graph.edges())} direct containments")
        
        return graph

    @staticmethod  
    def build_hierarchy_graph3(ercs, RN):
        """Build the minimal containment graph between ERCs - ultra-optimized version"""
        graph = nx.DiGraph()
        
        # Add all nodes
        for erc in ercs:
            graph.add_node(erc.label, erc=erc)
        
        if len(ercs) <= 1:
            return graph
        
        # 1. Compute all closures once and cache as sets
        closure_cache = {}
        for erc in ercs:
            closure_cache[erc.label] = set(species_list_to_names(erc.get_closure(RN)))
        
        # 2. Sort ERCs by closure size for more efficient processing
        ercs_by_size = sorted(ercs, key=lambda erc: len(closure_cache[erc.label]), reverse=True)
        
        # 3. Build containment relations with early termination
        containments = []
        for i, erc1 in enumerate(ercs_by_size):
            closure1 = closure_cache[erc1.label]
            # Only check ERCs with smaller or equal closure sizes
            for j in range(i+1, len(ercs_by_size)):
                erc2 = ercs_by_size[j]
                closure2 = closure_cache[erc2.label]
                
                # Since we sorted by size, if closure2 is not a subset, 
                # no subsequent ERCs will be either
                if len(closure2) > len(closure1):
                    break
                    
                if closure2.issubset(closure1) and closure1 != closure2:
                    containments.append((erc1.label, erc2.label))
        
        # 4. Use NetworkX transitive reduction for optimal performance
        if containments:
            try:
                # Build temporary graph with all containments
                temp_graph = nx.DiGraph()
                temp_graph.add_nodes_from([erc.label for erc in ercs])
                temp_graph.add_edges_from(containments)
                
                # Get transitive reduction
                direct_graph = nx.transitive_reduction(temp_graph)
                
                # Add direct edges to final graph
                graph.add_edges_from(direct_graph.edges())
                
            except (AttributeError, ImportError):
                # Fallback: optimized manual detection
                containment_set = set(containments)
                erc_labels = [erc.label for erc in ercs]
                
                for start, end in containments:
                    is_direct = True
                    for mid_label in erc_labels:
                        if (mid_label != start and mid_label != end and 
                            (start, mid_label) in containment_set and 
                            (mid_label, end) in containment_set):
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
    def plot_hierarchy(ercs, RN, graph=None, figsize=(10,10), title="RC Hierarchy"):
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
        
        plt.title(title)
        plt.axis('off')
        plt.show()
        return graph

    
__all__ = [
    'ERC',
    'generators',
    'species_list_to_names',
    'closure',
    'closures'
]



