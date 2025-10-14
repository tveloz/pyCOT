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
    Generate all unique generators based on their closures, not just supports.
    This is critical for networks with reactions that have the same support 
    but produce different closures (e.g., a+b->c and a+b->d).
    """
    gen = [closure(RN, [])]
    
    # Track seen CLOSURES to avoid duplicates (not just supports!)
    seen_closures = set()
    
    # Add the empty set closure signature
    empty_closure_sig = tuple(sorted([sp.name for sp in gen[0]]))
    seen_closures.add(empty_closure_sig)
    
    for reaction in RN.reactions():
        support = RN.get_supp_from_reactions(reaction)
        
        if support:
            # Compute the closure of this support
            support_closure = closure(RN, support)
            
            # Create a hashable signature for the CLOSURE (not the support)
            closure_signature = tuple(sorted([sp.name for sp in support_closure]))
            
            # Only add if we haven't seen this closure before
            if closure_signature not in seen_closures:
                gen.append(support)
                seen_closures.add(closure_signature)
    
    return gen

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
        """Create ERCs with closure pre-computation for efficiency"""
        print("Computing ERCs with closure caching...")
        
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
            
            # Pre-compute and cache closure for efficiency
            erc._closure = closure_cache[id(min_gens[0])]
            erc._closure_names = set(species_list_to_names(erc._closure))
            erc._RN_reference = id(RN)
            
            ERC_list.append(erc)
        
        print(f"Created {len(ERC_list)} ERCs with pre-computed closures")
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

    @staticmethod
    def plot_hierarchy(RN, ercs=None , graph=None, figsize=(10,10), title="ERC Hierarchy"):
        """Plot the ERC hierarchy"""
        
        if ercs is None:
            ercs=ERC.ERCs(RN)
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
                                     node_size=1000)
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


class ERC_Hierarchy:
    """Helper class to manage ERC hierarchies and containment relationships"""
    
    def __init__(self, RN_or_ercs, RN=None):
        """
        Initialize ERC hierarchy.
        
        Parameters
        ----------
        RN_or_ercs : ReactionNetwork or list
            Either a ReactionNetwork object or a list of pre-computed ERCs
        RN : ReactionNetwork, optional
            Required when passing a list of ERCs without RN references
        """
        if isinstance(RN_or_ercs, list):
            self.ercs = RN_or_ercs
            # Try to get RN from ERCs first, fallback to provided RN
            self.RN = getattr(self.ercs[0], 'RN', None) if self.ercs else None
            if self.RN is None:
                if RN is None:
                    raise ValueError("RN must be provided when ERCs don't have RN reference")
                self.RN = RN
        else:
            if not isinstance(RN_or_ercs, ReactionNetwork):
                raise TypeError(f"Expected ReactionNetwork or list of ERCs, got {type(RN_or_ercs)}")
            self.RN = RN_or_ercs
            self.ercs = ERC.ERCs(self.RN)
        
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


__all__ = [
    'ERC',
    'ERC_Hierarchy',
    'generators',
    'species_list_to_names',
    'closure',
    'closures'
]