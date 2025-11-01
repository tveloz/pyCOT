#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 29 18:55:14 2023
ERC_Hierarchy - Optimized with Closure Caching
@author: tveloz
"""
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

def closure(RN, species_set):
    """
    Compute the closure of a set of species in a reaction network.
    
    Parameters
    ----------
    RN : ReactionNetwork
        The reaction network object
    species_set : list or set
        Set of species to compute closure for
        
    Returns
    -------
    list
        Closure of the species set
    """
    # Verify RN is a ReactionNetwork object
    if not hasattr(RN, 'get_prod_from_species'):
        raise TypeError("RN must be a ReactionNetwork object")
        
    # Convert input to list if needed
    temp = list(species_set)
    prev_len = -1
    
    while prev_len != len(temp):
        prev_len = len(temp)
        # Get products from current species set
        CL = list(set(temp).union(set(RN.get_prod_from_species(temp))))
        temp = CL
        
    return temp


def generators(RN):
    """
    Get generators from reaction network with their closures.
    
    Parameters
    ----------
    RN : ReactionNetwork
        The reaction network object
        
    Returns
    -------
    list of tuples
        List of (support, closure) pairs where each support is paired with its computed closure
    """
    # Verify RN is a ReactionNetwork object
    if not hasattr(RN, 'get_prod_from_species'):
        raise TypeError("RN must be a ReactionNetwork object")
        
    # Initialize with empty set and its closure
    empty_support = []
    empty_closure = closure(RN, empty_support)
    gen_pairs = [(empty_support, empty_closure)]
    
    # Track seen supports to avoid duplicates
    seen_support_signatures = set()
    
    # Add the empty support signature
    empty_support_sig = tuple(sorted([sp.name for sp in empty_support])) if empty_support else ()
    seen_support_signatures.add(empty_support_sig)
    
    for reaction in RN.reactions():
        support = RN.get_supp_from_reactions(reaction)
        
        if support:  # Only add non-empty supports
            # Create a hashable signature for the support
            support_signature = tuple(sorted([sp.name for sp in support]))
            
            # Only add if we haven't seen this support before
            if support_signature not in seen_support_signatures:
                # Compute closure for this support
                support_closure = closure(RN, support)
                
                # Store the pair
                gen_pairs.append((support, support_closure))
                seen_support_signatures.add(support_signature)
    
    return gen_pairs







class ERC:
    def __init__(self, min_generators, label, all_generators=None):
        """Initialize an ERC object."""
        self.label = label
        self.min_generators = min_generators
        self.all_generators = all_generators
        # Simple caching without RN reference tracking
        self._closure = None
        self._closure_names = None
        self._reactions = None
        self._required_species = None
        self._produced_species = None

    def get_closure(self, RN):
        """Get closure with simple caching"""
        if self._closure is None:
            self._closure = closure(RN, self.min_generators[0])
            self._closure_names = set(species_list_to_names(self._closure))
        return self._closure
    
    def get_closure_names(self, RN):
        """Get closure as set of species names with caching"""
        if self._closure_names is None:
            self.get_closure(RN)  # This will populate _closure_names
        return self._closure_names
    
    def get_reacs(self, RN):
        """Get reactions with caching"""
        if self._reactions is None:
            specs = self.get_closure(RN)
            self._reactions = RN.get_reactions_from_species(specs)
        return self._reactions

    def get_required_species(self, RN):
        """Get required species (consumed but not produced) with caching"""
        if self._required_species is None:
            reactions = self.get_reacs(RN)
            
            consumed = set()
            produced = set()
            
            for reaction in reactions:
                for edge in reaction.edges:
                    if edge.type == "reactant":
                        consumed.add(edge.species_name)
                    elif edge.type == "product":
                        produced.add(edge.species_name)
            
            self._required_species = consumed - produced
        return self._required_species

    def get_produced_species(self, RN):
        """Get produced species with caching"""
        if self._produced_species is None:
            reactions = self.get_reacs(RN)
            
            produced = set()
            
            for reaction in reactions:
                for edge in reaction.edges:
                    if edge.type == "product":
                        produced.add(edge.species_name)
            
            self._produced_species = produced
        return self._produced_species
    
    def clear_cache(self):
        """Clear all cached computations"""
        self._closure = None
        self._closure_names = None
        self._reactions = None
        self._required_species = None
        self._produced_species = None

    @staticmethod
    def find_minimal_generators(generator_pairs):
        """
        Find minimal generators from a list of (support, closure) pairs.
        All pairs should have the same closure.
        
        Parameters
        ----------
        generator_pairs : list of tuples
            List of (support, closure) pairs
            
        Returns
        -------
        list
            List of minimal supports (without their closures)
        """
        minimal = []
        
        # Extract supports and create sets of species names for comparison
        supports_with_sets = [(support, set(sp.name for sp in support)) 
                            for support, _ in generator_pairs]
        
        # Check each support to see if it's minimal
        for i, (support_i, set_i) in enumerate(supports_with_sets):
            is_minimal = True
            
            # Compare with all other supports
            for j, (support_j, set_j) in enumerate(supports_with_sets):
                if i != j:
                    # If another support is a proper subset, this one is not minimal
                    if set_j.issubset(set_i) and set_j != set_i:
                        is_minimal = False
                        break
            
            if is_minimal:
                minimal.append(support_i)
        
        return minimal
    
    @staticmethod
    def ERCs(RN):
        """
        Create ERCs with explicit support-closure pairing.
        
        Parameters
        ----------
        RN : ReactionNetwork
            The reaction network object
            
        Returns
        -------
        list
            List of ERC objects
        """
        print("Computing ERCs with explicit support-closure pairs...")
        
        # 1. Get all (support, closure) pairs
        generator_pairs = generators(RN)
        
        # 2. Group by closure signature
        closure_groups = defaultdict(list)
        for support, support_closure in generator_pairs:
            # Create signature from closure
            closure_sig = tuple(sorted(sp.name for sp in support_closure))
            # Store the pair in the appropriate group
            closure_groups[closure_sig].append((support, support_closure))
        
        # 3. Find minimal generators for each group
        ERC_list = []
        for counter, (closure_sig, pairs_group) in enumerate(closure_groups.items()):
            # Find minimal supports within this group
            min_gens = ERC.find_minimal_generators(pairs_group)
            
            # Extract all supports (for all_generators field)
            all_gens = [support for support, _ in pairs_group]
            
            # Create ERC object
            erc = ERC(min_generators=min_gens,
                    label=f"E{counter}",
                    all_generators=all_gens)
            
            # Pre-compute and cache closure (using the first pair's closure since all are identical)
            erc._closure = pairs_group[0][1]
            erc._closure_names = set(species_list_to_names(erc._closure))
            erc._RN_reference = id(RN)
            
            ERC_list.append(erc)
        
        print(f"Created {len(ERC_list)} ERCs")
        print(f"  Total generator pairs found: {len(generator_pairs)}")
        print(f"  Total minimal generators found: {sum(len(erc.min_generators) for erc in ERC_list)}")
        
        return ERC_list
    
    @staticmethod  
    def build_hierarchy_graph(ercs, RN):
        """Build the minimal containment graph between ERCs - ultra-optimized version"""
        graph = nx.DiGraph()
        
        # Add all nodes
        for erc in ercs:
            graph.add_node(erc.label, erc=erc)
        
        if len(ercs) <= 1:
            return graph
        
        print(f"Building hierarchy graph for {len(ercs)} ERCs...")
        
        # 1. Sort ERCs by closure size for more efficient processing (using cached closures)
        ercs_by_size = sorted(ercs, key=lambda erc: len(erc.get_closure_names(RN)), reverse=True)
        
        # 2. Build containment relations with early termination
        containments = []
        for i, erc1 in enumerate(ercs_by_size):
            closure1 = erc1.get_closure_names(RN)
            # Only check ERCs with smaller or equal closure sizes
            for j in range(i+1, len(ercs_by_size)):
                erc2 = ercs_by_size[j]
                closure2 = erc2.get_closure_names(RN)
                
                # Since we sorted by size, if closure2 is not a subset, 
                # no subsequent ERCs will be either
                if len(closure2) > len(closure1):
                    break
                    
                if closure2.issubset(closure1) and closure1 != closure2:
                    containments.append((erc1.label, erc2.label))
        
        print(f"Found {len(containments)} containment relations")
        
        # 3. Use NetworkX transitive reduction for optimal performance
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
                
                print(f"Reduced to {len(direct_graph.edges())} direct containments")
                
            except (AttributeError, ImportError):
                # Fallback: optimized manual detection
                print("Using fallback method...")
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
                        
                print(f"Reduced to {len(graph.edges())} direct containments")
        
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

    


class ERC_Hierarchy:
    """Helper class to manage ERC hierarchies and containment relationships"""
    
    def __init__(self, RN, ercs):
        """
        Initialize ERC hierarchy.
        
        Parameters
        ----------
        RN_or_ercs : ReactionNetwork or list
            Either a ReactionNetwork object or a list of pre-computed ERCs
        RN : ReactionNetwork, optional
            Required when passing a list of ERCs without RN references
        """
        self.RN = RN
        if ercs is None:
            ercs= ERC.ERCs(RN)
        self.ercs = ercs
        self.graph = None
        self.build_hierarchy_graph()

    def build_hierarchy_graph(self):
        """Build hierarchy graph using ERCs."""
        if not hasattr(self, 'RN'):
            raise AttributeError("RN not initialized")
        self.graph = ERC.build_hierarchy_graph(self.ercs, self.RN)
        self._build_containment_caches()
        return self.graph
    
    def _build_containment_caches(self):
        """Build caches for containment relationships"""
        if self.graph is None:
            return
            
        self._containment_cache = {}
        self._contained_cache = {}
        
        for erc in self.ercs:
            # Get all ERCs that this ERC contains (descendants)
            descendants = set()
            if erc.label in self.graph:
                descendants = set(nx.descendants(self.graph, erc.label))
            self._containment_cache[erc] = [e for e in self.ercs if e.label in descendants]
            
            # Get all ERCs that contain this ERC (ancestors)
            ancestors = set()
            if erc.label in self.graph:
                ancestors = set(nx.ancestors(self.graph, erc.label))
            self._contained_cache[erc] = [e for e in self.ercs if e.label in ancestors]
    
    def get_contain(self, erc):
        """Get all ERCs contained by the given ERC (cached)"""
        if erc not in self._containment_cache:
            self._build_containment_caches()
        return self._containment_cache.get(erc, [])
    
    def get_contained(self, erc):
        """Get all ERCs that contain the given ERC (cached)"""
        if erc not in self._contained_cache:
            self._build_containment_caches()
        return self._contained_cache.get(erc, [])
    def get_erc_by_label(self, label):
        """Get ERC object by its label"""
        for erc in self.ercs:
            if erc.label == label:
                return erc
        raise ValueError(f"ERC with label '{label}' not found")
    
    def get_erc_from_reaction(self, RN, hierarchy, reaction):
        """
        Find the unique ERC that is GENERATED by this reaction.
        
        This implements the theoretical definition: ERC(r) = closure(supp(r))
        
        The key difference from the old version:
        - OLD: Returns first ERC that contains the reaction (wrong)
        - NEW: Returns the ERC that is generated by the reaction's support (correct)
        
        Parameters
        ----------
        RN : ReactionNetwork
            The reaction network
        hierarchy : ERC_Hierarchy
            The hierarchy object (not used in this version, but kept for compatibility)
        reaction : Reaction
            The reaction whose generating ERC we want to find
            
        Returns
        -------
        ERC or None
            The unique ERC generated by this reaction, or None if not found
        """
        # Step 1: Get the support (reactants) of the reaction
        support = RN.get_supp_from_reactions(reaction)
        
        if not support:
            # Inflow reaction (no support), handle specially
            # These typically generate the empty/inflow ERC
            return None
        
        # Step 2: Compute the closure of the support
        # This is the definition: ERC(r) = closure(supp(r))
        reaction_closure = closure(RN, support)
        reaction_closure_names = set(species_list_to_names(reaction_closure))
        
        # Step 3: Find the ERC with exactly matching closure
        for erc in self.ercs:
            erc_closure_names = erc.get_closure_names(RN)
            
            # Must be exact match, not just containment
            if erc_closure_names == reaction_closure_names:
                return erc
        
        # Should not happen if ERCs were built correctly, but handle gracefully
        return None
    def plot_hierarchy(self, figsize=(10,10), title="ERC Hierarchy"):
        """
        Plot the ERC hierarchy.
        
        Parameters
        ----------
        figsize : tuple, optional
            Figure size (width, height)
        title : str, optional
            Plot title
            
        Returns
        -------
        networkx.DiGraph
            The hierarchy graph
        """
        if self.graph is None:
            raise ValueError("Hierarchy graph not built")
            
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)

        levels = ERC.get_node_levels(self.graph)  # Keep this as utility
        pos = {}
        level_nodes = defaultdict(list)

        for node, level in levels.items():
            level_nodes[level].append(node)

        # Position nodes by level
        for level, nodes in level_nodes.items():
            n_nodes = len(nodes)
            nodes.sort(key=lambda n: len(nx.ancestors(self.graph, n)), reverse=True)
            for i, node in enumerate(nodes):
                x = (i - (n_nodes-1)/2) * 2.0
                y = level * 2.0
                pos[node] = (x, y)

        # Draw graph components
        nodes = nx.draw_networkx_nodes(self.graph, pos,
                                     node_color='lightblue',
                                     node_size=1000)
        nx.draw_networkx_edges(self.graph, pos,
                             edge_color='gray',
                             arrows=True,
                             arrowsize=20)
        nx.draw_networkx_labels(self.graph, pos,
                              font_size=12,
                              font_weight='bold')

        # Add hover functionality using self.ercs and self.RN
        annot = ax.annotate("",
                           xy=(0,0), xytext=(20,20),
                           textcoords="offset points",
                           bbox=dict(boxstyle="round", fc="w", ec="0.5", alpha=0.9),
                           ha='center',
                           visible=False)

        def update_annot(ind):
            node = list(self.graph.nodes())[ind["ind"][0]]
            erc = next((e for e in self.ercs if e.label == node), None)
            if erc:
                closure_text = f"{species_list_to_names(erc.get_closure(self.RN))}"
            else:
                closure_text = "Hover data not available"
            pos_node = pos.get(node, (0,0))
            annot.xy = pos_node
            annot.set_text(closure_text)

        def hover(event):
            if event.inaxes == ax:
                cont, ind = nodes.contains(event)
                annot.set_visible(cont)
                if cont:
                    try:
                        update_annot(ind)
                        fig.canvas.draw_idle()
                    except Exception:
                        pass

        fig.canvas.mpl_connect("motion_notify_event", hover)

        plt.title(title)
        plt.axis('off')
        plt.show()
        return self.graph

__all__ = [
    'ERC',
    'ERC_Hierarchy',
    'generators',
    'species_list_to_names',
    'closure',
    'closures'
]