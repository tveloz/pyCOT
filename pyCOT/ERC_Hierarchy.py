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
    def find_minimal_synergies(base_erc, target_erc, erc_hierarchy, RN, visited_branch_pairs=None):
        """
        Find minimal synergetic ERCs by exploring hierarchy from top to bottom.
        Tracks visited branch pairs to avoid redundant exploration.
        
        Parameters
        ----------
        base_erc : ERC
            The base ERC to find synergies with
        target_erc : ERC
            The target ERC whose generators we want to cover
        erc_hierarchy : nx.DiGraph
            The hierarchy graph of ERCs
        RN : ReactionNetwork
            The reaction network
        visited_branch_pairs : set, optional
            Set of already explored branch pairs (base_branch, syn_branch, target)
            
        Returns
        -------
        dict
            Dictionary mapping minimal synergetic ERCs to tuples of:
            - list of uniquely covered generators
            - coverage ratio
        """
        if visited_branch_pairs is None:
            visited_branch_pairs = set()

        def get_branch_nodes(node):
            """Get all nodes in a branch (ancestors and descendants)"""
            ancestors = nx.ancestors(erc_hierarchy, node)
            descendants = nx.descendants(erc_hierarchy, node)
            return {node}.union(ancestors).union(descendants)

        def has_synergy(erc1, erc2, target_gens):
            """Check if two ERCs have synergy for covering target generators"""
            erc1_closure = set(species_list_to_names(erc1.get_closure(RN)))
            erc2_closure = set(species_list_to_names(erc2.get_closure(RN)))
            combined_closure = erc1_closure.union(erc2_closure)
            covered = []
            
            for gen in target_gens:
                gen_species = set(species_list_to_names(gen))
                if (gen_species.issubset(combined_closure) and
                    not gen_species.issubset(erc1_closure) and
                    not gen_species.issubset(erc2_closure) and
                    len(gen_species & erc1_closure - erc2_closure) > 0 and
                    len(gen_species & erc2_closure - erc1_closure) > 0):
                    covered.append(gen)
            
            return covered

        def explore_branch(current_erc, remaining_gens, minimal_synergies):
            """Recursively explore hierarchy branch for minimal synergies"""
            # Get current branch nodes
            current_branch = get_branch_nodes(current_erc.label)
            base_branch = get_branch_nodes(base_erc.label)
            
            # Check if this branch pair has been explored for this target
            branch_pair = frozenset([frozenset(base_branch), frozenset(current_branch)])
            branch_target_key = (branch_pair, target_erc.label)
            
            if branch_target_key in visited_branch_pairs:
                return
            
            visited_branch_pairs.add(branch_target_key)
            
            covered_gens = has_synergy(base_erc, current_erc, remaining_gens)
            
            if not covered_gens:
                return  # No synergy in this branch
                
            # Get successor nodes (contained ERCs)
            successors = list(erc_hierarchy.successors(current_erc.label))
            if not successors:
                # If leaf node with synergy, add to minimal synergies
                coverage_ratio = len(covered_gens) / len(target_erc.min_generators)
                minimal_synergies[current_erc] = (covered_gens, coverage_ratio)
                return
                
            # Try successors for better (more specific) synergies
            found_better = False
            for succ_label in successors:
                succ_erc = erc_hierarchy.nodes[succ_label]['erc']
                succ_covered = has_synergy(base_erc, succ_erc, covered_gens)
                if succ_covered:
                    found_better = True
                    explore_branch(succ_erc, succ_covered, minimal_synergies)
            
            # If no successor had synergy, current is minimal
            if not found_better:
                coverage_ratio = len(covered_gens) / len(target_erc.min_generators)
                minimal_synergies[current_erc] = (covered_gens, coverage_ratio)

        # Get top-level ERCs that aren't in containment with base or target
        top_ercs = [erc_hierarchy.nodes[node]['erc'] for node in erc_hierarchy.nodes()
                   if erc_hierarchy.in_degree(node) == 0 and
                   node not in [base_erc.label, target_erc.label]]
        
        minimal_synergies = {}
        # Explore each top-level ERC
        for top_erc in top_ercs:
            explore_branch(top_erc, target_erc.min_generators, minimal_synergies)
            
        return minimal_synergies
    
    @staticmethod
    def find_hierarchical_minimal_synergies(erc_hierarchy, RN):
        """Find minimal synergies by exploring hierarchy chains from bottom to top."""
        
        def get_ancestor_chains(node):
            """Get all possible chains from node to roots"""
            chains = []
            if erc_hierarchy.in_degree(node) == 0:  # Node is a root
                return [[node]]
            
            for parent in erc_hierarchy.predecessors(node):
                parent_chains = get_ancestor_chains(parent)
                for chain in parent_chains:
                    chains.append(chain + [node])
            
            return chains
        
        def has_synergy(erc1, erc2, target_gens):
            """Check if two ERCs have synergy for target generators"""
            erc1_closure = set(species_list_to_names(erc1.get_closure(RN)))
            erc2_closure = set(species_list_to_names(erc2.get_closure(RN)))
            combined_closure = erc1_closure.union(erc2_closure)
            covered = []
            
            for gen in target_gens:
                gen_species = set(species_list_to_names(gen))
                if (gen_species.issubset(combined_closure) and
                    not gen_species.issubset(erc1_closure) and
                    not gen_species.issubset(erc2_closure)):
                    covered.append(gen)
            
            return covered

        def explore_synergy_chain(base_erc, syn_erc, target_erc, visited_combinations):
            """Check if synergy exists and is minimal"""
            combo = frozenset([base_erc.label, syn_erc.label, target_erc.label])
            if combo in visited_combinations:
                return None
            
            visited_combinations.add(combo)
            covered_gens = has_synergy(base_erc, syn_erc, target_erc.min_generators)
            
            if covered_gens:
                return (syn_erc, covered_gens, len(covered_gens)/len(target_erc.min_generators))
            return None

        minimal_synergies = {}
        visited_combinations = set()
        
        # Get bottom ERCs (no outgoing edges)
        bottom_ercs = [erc_hierarchy.nodes[n]['erc'] for n in erc_hierarchy.nodes()
                      if erc_hierarchy.out_degree(n) == 0]
        
        # Get top ERCs (no incoming edges)
        top_ercs = [erc_hierarchy.nodes[n]['erc'] for n in erc_hierarchy.nodes()
                    if erc_hierarchy.in_degree(n) == 0]
        
        # Process each bottom ERC as potential base
        for base_erc in bottom_ercs:
            # Get all possible chains from this bottom ERC to roots
            base_chains = get_ancestor_chains(base_erc.label)
            base_chain_nodes = {node for chain in base_chains for node in chain}
            
            # For each top ERC not in any of base's chains
            for top_erc in top_ercs:
                if top_erc.label not in base_chain_nodes:
                    # Get all possible chains from this top ERC
                    top_chains = get_ancestor_chains(top_erc.label)
                    
                    # Try synergies with ERCs from each chain
                    for chain in top_chains:
                        chain_ercs = [erc_hierarchy.nodes[n]['erc'] for n in chain]
                        
                        for syn_erc in chain_ercs:
                            if syn_erc.label in {base_erc.label, top_erc.label}:
                                continue
                            
                            # Try synergy with chain ERCs from top to bottom
                            for target_erc in reversed(chain_ercs):
                                synergy = explore_synergy_chain(base_erc, syn_erc, target_erc, visited_combinations)
                                if synergy:
                                    key = (base_erc, target_erc)
                                    if key not in minimal_synergies:
                                        minimal_synergies[key] = []
                                    minimal_synergies[key].append(synergy)
                                    break  # Found minimal synergy in this chain
        
        return minimal_synergies


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

    @staticmethod
    def plot_hierarchy_synergies(ercs, RN, hierarchy_graph=None, synergies=None, figsize=(12, 8), title=None):
        """Plot ERC hierarchy with synergies as colored edges with hover information."""
        if hierarchy_graph is None:
            hierarchy_graph = ERC.build_hierarchy_graph(ercs, RN)
            
        fig, ax = plt.subplots(figsize=figsize)
        G = nx.DiGraph()
        
        # Add nodes with species information
        for erc in ercs:
            species_info = f"Species: {species_list_to_names(erc.get_closure(RN))}"
            G.add_node(erc.label, node_type="erc", species=species_info)
        
        # Add hierarchy edges (black)
        for u, v in hierarchy_graph.edges():
            G.add_edge(u, v, color='black', edge_type='hierarchy')
        
        # Use standardized get_node_levels
        levels = ERC.get_node_levels(G)
        
        # Position nodes by level
        pos = {}
        level_nodes = defaultdict(list)
        
        for node, level in levels.items():
            level_nodes[level].append(node)
        
        for level, nodes in level_nodes.items():
            n_nodes = len(nodes)
            nodes.sort(key=lambda n: len(nx.ancestors(G, n)), reverse=True)
            for i, node in enumerate(nodes):
                x = (i - (n_nodes-1)/2) * 2.0
                y = level * 2.0
                pos[node] = (x, y)

        # Add synergy edges without labels
        if synergies:
            for (base_erc, target_erc), syn_list in synergies.items():
                for syn_erc, _, _ in syn_list:
                    synergy_info = f"{base_erc.label}+{syn_erc.label}→{target_erc.label}"
                    G.add_edge(base_erc.label, target_erc.label,
                             color='red', edge_type='synergy',
                             synergy_info=synergy_info)

        # Draw the graph
        nodes = nx.draw_networkx_nodes(G, pos, node_color='lightblue', 
                                     node_size=500, ax=ax)
        
        # Draw edges with different styles
        hierarchy_edges = [(u, v) for u, v in G.edges() if G[u][v]['edge_type'] == 'hierarchy']
        synergy_edges = [(u, v) for u, v in G.edges() if G[u][v]['edge_type'] == 'synergy']
        
        nx.draw_networkx_edges(G, pos, edgelist=hierarchy_edges,
                             edge_color='black',
                             arrows=True,
                             arrowsize=20)
        
        nx.draw_networkx_edges(G, pos, edgelist=synergy_edges,
                             edge_color='red',
                             style='dashed',
                             arrows=True,
                             arrowsize=20,
                             connectionstyle="arc3,rad=0.3")
        
        nx.draw_networkx_labels(G, pos)

        # Create hover annotation
        annot = ax.annotate("", xy=(0,0), xytext=(10,10),
                          textcoords="offset points",
                          bbox=dict(boxstyle="round", fc="w", ec="0.5", alpha=0.9),
                          arrowprops=dict(arrowstyle="->"))
        annot.set_visible(False)

        def update_annot(obj_type, idx=None, edge=None):
            if obj_type == 'node':
                node = list(G.nodes())[idx]
                text = G.nodes[node]['species']
                pos_node = pos[node]
                annot.xy = pos_node
            else:  # edge
                text = G.edges[edge]['synergy_info']
                pos_src = pos[edge[0]]
                pos_dst = pos[edge[1]]
                annot.xy = ((pos_src[0] + pos_dst[0])/2, ((pos_src[1] + pos_dst[1])/2))
            annot.set_text(text)

        def hover(event):
            if event.inaxes != ax:
                return
            
            # Check for nodes
            cont, ind = nodes.contains(event)
            if cont:
                annot.set_visible(True)
                update_annot('node', ind['ind'][0])
                fig.canvas.draw_idle()
                return
            
            # Check for synergy edges
            for edge in synergy_edges:
                edge_pos = np.array([[pos[edge[0]][0], pos[edge[0]][1]],
                                   [pos[edge[1]][0], pos[edge[1]][1]]])
                dist = np.linalg.norm(np.cross(edge_pos[1]-edge_pos[0],
                                   edge_pos[0]-np.array([event.xdata, event.ydata])))/np.linalg.norm(edge_pos[1]-edge_pos[0])
                if dist < 0.1:
                    annot.set_visible(True)
                    update_annot('edge', edge=edge)
                    fig.canvas.draw_idle()
                    return
            
            annot.set_visible(False)
            fig.canvas.draw_idle()

        fig.canvas.mpl_connect("motion_notify_event", hover)
        if title is None:
            ax.set_title("ERC Hierarchy with Synergies\n(Hover over red edges to see synergy details)")
        else:
            ax.set_title(title)
        ax.set_axis_on()
        ax.grid(False)
        
        return fig

    @staticmethod
    def plot_hierarchy_html(ercs, RN, hierarchy_graph=None, figsize=(800, 600)):
        """Generate interactive HTML plot of ERC hierarchy."""
        if hierarchy_graph is None:
            hierarchy_graph = ERC.build_hierarchy_graph(ercs, RN)
        
        # Use get_node_levels for consistent level calculation
        levels = ERC.get_node_levels(hierarchy_graph)
        
        # Position nodes
        pos = {}
        level_nodes = defaultdict(list)
        
        for node, level in levels.items():
            level_nodes[level].append(node)
        
        for level, nodes in level_nodes.items():
            n_nodes = len(nodes)
            nodes.sort()
            for i, node in enumerate(nodes):
                x = (i - (n_nodes-1)/2) * 2.0
                y = level * 2.0
                pos[node] = (x, y)
        
        # Create edges
        edge_x = []
        edge_y = []
        for edge in hierarchy_graph.edges():
            x0, y0 = pos[edge[0]]
            x1, y1 = pos[edge[1]]
            edge_x.extend([x0, x1, None])
            edge_y.extend([y0, y1, None])
        
        # Create figure
        fig = go.Figure()
        
        # Add edges
        fig.add_trace(go.Scatter(
            x=edge_x, y=edge_y,
            line=dict(width=0.5, color='#888'),
            hoverinfo='none',
            mode='lines'))
        
        # Add nodes
        node_x = []
        node_y = []
        node_text = []
        
        for node in hierarchy_graph.nodes():
            x, y = pos[node]
            node_x.append(x)
            node_y.append(y)
            erc = hierarchy_graph.nodes[node]['erc']
            hover_text = f"ERC {erc.label}<br>Species: {species_list_to_names(erc.get_closure(RN))}"
            node_text.append(hover_text)
        
        fig.add_trace(go.Scatter(
            x=node_x, y=node_y,
            mode='markers+text',
            hoverinfo='text',
            hovertext=node_text,
            text=[f"ERC{node}" for node in hierarchy_graph.nodes()],
            textposition="top center",
            marker=dict(
                showscale=False,
                size=20,
                line_width=2)))
        
        fig.update_layout(
            title="ERC Hierarchy",
            showlegend=False,
            hovermode='closest',
            margin=dict(b=20,l=5,r=5,t=40),
            width=figsize[0],
            height=figsize[1],
            plot_bgcolor='white')
        
        return fig

    @staticmethod
    def plot_hierarchy_synergies_html(ercs, RN, hierarchy_graph=None, synergies=None, figsize=(800, 600)):
        """Generate interactive HTML plot of ERC hierarchy with synergies."""
        if hierarchy_graph is None:
            hierarchy_graph = ERC.build_hierarchy_graph(ercs, RN)
        
        # Get base hierarchy plot
        fig = ERC.plot_hierarchy_html(ercs, RN, hierarchy_graph, figsize)
        
        if synergies:
            # Get node positions (same as in plot_hierarchy_html)
            levels = ERC.get_node_levels(hierarchy_graph)
            pos = {}
            level_nodes = defaultdict(list)
            
            for node, level in levels.items():
                level_nodes[level].append(node)
            
            for level, nodes in level_nodes.items():
                n_nodes = len(nodes)
                nodes.sort()
                for i, node in enumerate(nodes):
                    x = (i - (n_nodes-1)/2) * 2.0
                    y = level * 2.0
                    pos[node] = (x, y)
            
            # Add synergy edges
            for (base_erc, target_erc), syn_list in synergies.items():
                for syn_erc, covered_gens, ratio in syn_list:
                    x0, y0 = pos[base_erc.label]
                    x1, y1 = pos[target_erc.label]
                    
                    hover_text = f"{base_erc.label}+{syn_erc.label}→{target_erc.label}"
                    
                    fig.add_trace(go.Scatter(
                        x=[x0, x1],
                        y=[y0, y1],
                        mode='lines',
                        line=dict(color='red', width=1, dash='dash'),
                        hoverinfo='text',
                        hovertext=hover_text))
        
        return fig
__all__ = [
    'ERC',
    'generators',
    'species_list_to_names',
    'closure',
    'closures'
]



