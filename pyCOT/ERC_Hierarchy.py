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
from collections import defaultdict
from collections import Counter # para el gráfico de cantidad de básicos vs repeticiones
from pyCOT.file_manipulation import *
from pyCOT.display import *
from pyCOT.reaction_network import *
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout

class ERC:
    def __init__(self, min_generators, label, all_generators=None):
        """
        Initialize an ERC object.
        
        :param min_generators: List of minimal generators (sets of species).
        :param species: List of species in the closure.
        """
        self.min_generators = min_generators
        self.all_generators =all_generators 
        self.label= label

    def get_closure(self,RN):
        return closure(RN,self.min_generators[0])
    
    def get_reacs(self,RN):
        specs=self.get_closure(RN)
        reacs=RN.get_reactions_from_species(specs)
        return reacs

class Hierarchy_ERC:
    def __init__(self, RN):
        """
        Initialize the Hierarchy_ERC object.

        :param ercs: List of ERC objects.
        """
        self.ercs =[]
        #Calculating the elements to build self.ercs
        ERCx=ERCs(RN)   
        for erc in ERCx:
            self.ercs.append(ERC(erc[0],erc[3],erc[2]))
        #Create containment matrix labelled
        #self.containments=
        
        #Get direct containments to make containment forest at the level of labels
        label_direct_containments= get_direct_containments(RN,ERCx)
        self.graph = nx.DiGraph()
        
        # Add nodes and edges
        for erc in self.ercs:
            self.graph.add_node(erc.label)
        for start, end in label_direct_containments:
            self.graph.add_edge(start, end)
    ########### Graph functions #######################
    def __str__(self):
        """
        Return a string representation of the Hierarchy_ERC object.

        :return: String showing the direct containment graph.
        """
        graph_str = "Hierarchy_ERC Containment Graph:\n"
        for edge in self.graph.edges:
            graph_str += f"{edge}\n"
        return graph_str

    ########## Getters from inner node conditions#######################
    def get_erc_from_label(self, label):
        """
        Retrieve the ERC object with the given label.

        :param label: The label of the ERC to retrieve.
        :return: The ERC object with the specified label, or None if not found.
        """
        for erc in self.ercs:
            if erc.label == label:
                return erc
        print("get_ERC was called with a wrong label= "+label)
        return None
    def get_erc_from_reaction(self,RN,r):   
        erc_species=closure(RN,RN.get_supp_from_reactions(r))
        for lst in self.ercs:
            if set(lst.get_closure(RN)) == set(erc_species):
                return lst
        print("get_erc_from_reaction for reaction "+str(r)+" did not find associated ERC")
        return None
         
    def get_erc_from_generator(self,species):
        erc_species=closure(RN,species)
        for lst in self.ercs:      
            if set(lst.get_closure) == erc_species:
                return lst
        print("get_erc_from_reaction for species "+str(species)+" did not find associated ERC")
        return None
    ############# Get ERCs from conditions comparing with other nodes###############

    def get_contained(self, erc, itself=False):
        """
        Retrieve the list of containments for a given ERC.

        :param label: The label of the ERC to query.
        :return: List of labels of directly contained ERCs.
        """
        label=erc.label
        if not self.graph.has_node(label):
            raise ValueError(f"The node {label} is not in the graph.")
    
        # Perform DFS or BFS to find all reachable nodes
        reachable_nodes = nx.descendants(self.graph, label)
        # Include the node itself
        if itself:
            reachable_nodes.add(label)
        
        return reachable_nodes
    
    def get_contain(self, erc, itself=False):
        """
        Get all nodes that can reach the given node transitively in a directed graph.
        
        :param graph: A NetworkX DiGraph object.
        :param node: The target node.
        :return: A set of all nodes that have paths to the given node.
        """
        label=erc.label
        if not self.graph.has_node(label):
            raise ValueError(f"The node {label} is not in the graph.")
        
        # Find all ancestors (nodes that have a path to the target node)
        reachable_nodes = nx.ancestors(self.graph, label)
        
        # Include the node itself
        if itself:
            reachable_nodes.add(label)
        
        return reachable_nodes
    
    def get_potential_synergies(self, RN, erc):
        """
        Get ERCs that have potential synergies with the given ERC.
        
        Parameters
        ----------
        erc : ERC
            The ERC to find potential synergies for
            
        Returns
        -------
        list[ERC]
            List of ERCs that have potential synergies with the given ERC
        """
        # Get the species closure for this ERC
        closure_species = erc.get_closure(RN)
        
        # Get reactions that partially intersect with this closure
        partial_reactions = RN.get_reactions_partially_intersecting_support(closure_species)
        
        # Get the ERCs associated with these reactions
        synergy_ercs = []
        for reaction in partial_reactions:
            erc_from_reaction = self.get_erc_from_reaction(RN, reaction)
            if erc_from_reaction is not None:
                synergy_ercs.append(erc_from_reaction)
                
        # Remove duplicates while preserving order
        return list(dict.fromkeys(synergy_ercs))
    
# Example usage
# Assume `get_direct_containments` is defined elsewhere and works as described.
# ercs = [
#     ERC(min_generators=[{'A', 'B'}], species=['A', 'B', 'C'], label="ERC1"),
#     ERC(min_generators=[{'B'}], species=['B', 'C'], label="ERC2"),
#     ERC(min_generators=[{'C'}], species=['C'], label="ERC3")
# ]

# hierarchy = Hierarchy_ERC(ercs)
# print(hierarchy)


########ERC constructor functions###################
def generators(RN):
    gen = []
    for r in [reaction.name() for reaction in RN.reactions()]:
        gen.append(RN.get_supp_from_reactions(r))
    return gen

def closure(RN,X):
    temp=X
    CL=list(set(temp).union(set(RN.get_prod_from_species(temp))))
    while RN.get_reactions_from_species(CL)!=RN.get_reactions_from_species(temp):
        temp=CL
        CL=list(set(temp).union(set(RN.get_prod_from_species(temp))))
    return CL

def closures(RN,ListX):
    L=[]
    for X in ListX:
        L.append(closure(RN,X))
    return L

def identify_minimal_generators(existing_erc, new_gen, new_closure):
    """
    Identifies if and how a new generator should be added to existing minimal generators
    Returns: (is_new_closure, updated_min_gens)
    """
    for erc in existing_erc:
        if set(new_closure) == set(erc[1]):
            # Same closure - check if we need to update minimal generators
            updated_min_gens = erc[0].copy()
            new_is_minimal = True
            
            # Check against existing minimal generators
            for i, min_gen in enumerate(erc[0]):
                if all(item in min_gen for item in new_gen):
                    new_is_minimal = False
                    if len(new_gen) < len(min_gen):
                        updated_min_gens[i] = new_gen
                
            if new_is_minimal:
                updated_min_gens.append(new_gen)
            return False, updated_min_gens
            
    return True, [new_gen]

def ERCs(RN):
    List_gen = generators(RN)
    List_closures = closures(RN, List_gen)
    ERC = []
    
    for i in range(len(List_gen)):
        is_new_closure, min_gens = identify_minimal_generators(ERC, List_gen[i], List_closures[i])
        
        if is_new_closure:
            ERC.append([min_gens, List_closures[i], List_gen[i]])
        else:
            # Update the minimal generators for existing closure
            for erc in ERC:
                if set(List_closures[i]) == set(erc[1]):
                    erc[0] = min_gens
                    break
    
    # Add labels
    for i in range(len(ERC)):
        ERC[i].append(f"E{i}")
    
    return ERC
############################# GEN Constructor Functions##################################
##############################Getting Direct Containments#################################
#define X contains Y
def set_containment(X,Y):
    return all(item in X for item in Y)
#define X strictly contains Y
def set_strict_containment(X,Y):
    return all(item in X for item in Y) and len(Y)<len(X)

#Get the containment relations among ERCs
def get_containments(RN,ERC):
    erc=ERC.copy()
    containments=[]
    erc_closures=[inner_list[1] for inner_list in erc]
    for i in range(len(erc)):
        erc_i=erc_closures[i]
        for j in range(i,len(erc)):
            erc_j=erc_closures[j]
            if set_strict_containment(erc_i, erc_j):  
                containments.append([i,j])
            if  set_strict_containment(erc_j, erc_i):
                containments.append([j,i])
    return containments

def get_direct_containments(RN,ERC):
    containments=get_containments(RN,ERC)
    direct_containments = []
    #indexes that must be eliminated must be as start and end in two different ERCs
    flattened_list = list(chain.from_iterable(containments))
    erc_dubious_indexes=list(set(flattened_list))
    
    #print("erc_dubious_indexes="+str(erc_dubious_indexes))    

    for pair in containments:
        start, end = pair
        to_add=True
        # Check if there is no intermediate element between start and end
        for mid in erc_dubious_indexes:
            if ([start, mid] in containments) and ([mid, end] in containments):
                to_add=False
                
        if to_add:
            direct_containments.append(pair)
    direct_label_containments=[]
    
    #Obtain the labels of the ERC containments to return it in the format to be used by Hierarchy_ERC class 
    for pair in direct_containments:
        start, end=pair
        direct_label_containments.append([f"E{start}", f"E{end}"])
    return direct_label_containments



