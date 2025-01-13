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
from pyCOT.ERC_Hierarchy import *
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout


class GEN:
    def __init__(self, RN):
        """
        Initialize a GEN object.
        
        :param RN: Reaction Network object used to compute ERCs and relations.
        """
        # Compute ERCs and their attributes
        self.ERCs = [
            ERC(min_generators=erc[0], species=erc[1], label="E"+str(erc[3]))
            for erc in ERCs(RN)
        ]

        # Create the containment graph
        #containments = get_containments(RN, ERCs(RN))
        direct_containments = get_direct_containments(RN, ERCs(RN))
        self.graph = nx.DiGraph()
        
        # Add nodes and edges
        for i, erc in enumerate(self.ERCs):
            self.graph.add_node(f"E{i}")
        for start, end in direct_containments:
            self.graph.add_edge(f"E{start}", f"E{end}")
    
    def get_erc_level(self, erc_index):
        node = f"E{erc_index}"
        if not self.graph.has_node(node):
            raise ValueError(f"Node {node} not found in the graph.")

        # Rename the graph for traversal
        rev_graph = self.graph

    # Perform DFS from the given node to calculate the longest path
        def dfs(node, visited):
            visited.add(node)
            max_depth = 0
            for neighbor in rev_graph.neighbors(node):
                if neighbor not in visited:
                    max_depth = max(max_depth, dfs(neighbor, visited) + 1)
            visited.remove(node)
            return max_depth

        # Call DFS starting from the given node
        return dfs(node, set())


########ERC constructor functions###################
def generators(RN):
    gen=[]
    for r in RN.RnStr:
        gen.append(RN.get_supp_from_reactions(r))
    return gen
#Compute closure of a set of species
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
###ERCs produces a list of four-tuple:
    #0th coordinate stores the list of minimal generators for a given ERC (sets of species)
    #1st coordinate stores the ERC as a set of species
    #2nd coordinate stores the reactions that generate same closure (equivalence class)
    #3rd coordinate stores the label

def ERCs(RN):
    List_gen=generators(RN).copy()
    List_closures=closures(RN,List_gen).copy()
    ERC=[]
    #print("To check "+str(len(List_gen))+" generators")
    for i in range(len(List_gen)):
        #print("Start Checking "+str(List_gen[i]))
        novel_min_gen=True
        ###Check that List_gen[i] needs to be added to min_gen###
        for j in range(len(ERC)):
        #Step1: Check List_gen[i] has the same closure than min_gen[j]    
            #print("Passing by min_gen "+str(j)+" whose closure is "+str(min_gen[j][1]))
            if set(List_closures[i])==set(ERC[j][1]):
                 #print("Closure same as generator "+str(j))
                 #Then ERC associated to List_gen[i] is already in min_gen, no need to append it as a new ERC
                 novel_min_gen=False
                 #Now we know we must check either replace current min_gen by a smaller one, add as new to the list or not add List_gen[i]
                 #new min_gen accounts for a local addition to current list of min_gen
                 new_min_gen=True
                 for k in range(len(ERC[j][0])):
                     #print("Starting comparison with element "+str(k)+ " of minimal generator list "+str(min_gen[j][0]))
                     #Step 2: Check if List_gen[i] is equal or smaller to some in for min_gen[j][0][k]
                     if all(item in ERC[j][0][k] for item in List_gen[i]):
                         #print("Generator is smaller than min_gen stored "+str(k)+" which is "+str(min_gen[j][0][k]))
                         #If so, we know there is no new_min_gen to add, but it can replace current ones  
                         new_min_gen=False
                         #verify containment is strict for replacing min_gen[j][0][k] by List_gen[i]
                         if len(List_gen[i])<len(ERC[j][0][k]):
                             ERC[j][0][k]=List_gen[i]
                        #If the latter condition does not meet, we skip List_gen[i] because is not minimal    
                 #In case the new_min_gen never changed to false, we know its minimal generator to keep tracK
                 if new_min_gen:        
                     #print("The generator "+str(List_gen[i])+" is new so far")
                     ERC[j][0].append(List_gen[i])
        #Once all has been checked if novel_min_gen remains True we must add because its a new closure
        if novel_min_gen:
            ERC.append([[List_gen[i]],List_closures[i],RN.get_reactions_from_species(List_closures[i])])
        #add the indexation to ERCs as fourth coordinate
    for i in range(len(ERC)):
        ERC[i].append("E"+str(i))
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
    return direct_containments

#########################ERC relations to compute Synergies efficiently######################33
def get_higher_ercs(gen, erc):
    """
    Get all nodes that transitively point to the erc in the gen.

    :param graph: A NetworkX DiGraph object.
    :param target_node: The target node in the graph.
    :return: A set of all nodes that transitively point to the target node.
    """
    # Reverse the graph to invert the direction of edges
    reversed_graph = gen.graph.reverse(copy=True)
    
    # Perform a depth-first search (or BFS) from the target node in the reversed graph
    predecessors = nx.dfs_predecessors(reversed_graph, erc)
    
    # Include the target node itself and all its reachable nodes
    result_set = set(predecessors.keys()) | {erc}
    
    return result_set
def get_lower_ercs(gen, erc):
    """
    Get all nodes that are transitively pointed from a given erc.

    :param graph: A NetworkX DiGraph object.
    :param target_node: The target node in the graph.
    :return: A set of all nodes that transitively point to the target node.
    """
    
    # Perform a depth-first search (or BFS) from the target node in the reversed graph
    successors = nx.dfs_predecessors(gen.graph, erc)
    
    # Include the target node itself and all its reachable nodes
    result_set = set(successors.keys()) | {erc}
    
    return result_set