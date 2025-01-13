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
    def __init__(self, min_generators, label):
        """
        Initialize an ERC object.
        
        :param min_generators: List of minimal generators (sets of species).
        :param species: List of species in the closure.
        """
        self.min_generators = min_generators
        self.label= label

    def get_closure(self,RN):
        return closure(RN,self.min_generators[0])
    
    def get_reacs(self,RN):
        specs=self.get_closure(RN)
        reacs=RN.get_reactions_from_species(specs)

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
            self.ercs.append(ERC(erc[0],erc[3]))
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
    def get_erc_from_reaction(self,RN,hierarchy,r):   
        erc_species=closure(RN,RN.get_supp_from_reactions(r))
        for lst in hierarchy.ercs:
            if set(lst.get_closure(RN)) == set(erc_species):
                return lst
        print("get_erc_from_reaction for reaction "+str(r)+" did not find associated ERC")
        return None
         
    def get_erc_from_generator(self,hierarchy,species):
        erc_species=closure(RN,species)
        for lst in hierarchy.ercs:      
            if set(lst.get_closure) == erc_species:
                return lst
        print("get_erc_from_reaction for species "+str(species)+" did not find associated ERC")
        return None

   

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
    def leveled_ERCs(self):
        ERCs=self.ercs
        fundamental=[]
        for erc in self.ercs:
            if len(get_contained(erc))==0:
                fundamental.append(erc)

        

    def __str__(self):
        """
        Return a string representation of the Hierarchy_ERC object.

        :return: String showing the direct containment graph.
        """
        graph_str = "Hierarchy_ERC Containment Graph:\n"
        for edge in self.graph.edges:
            graph_str += f"{edge}\n"
        return graph_str
    
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
    direct_label_containments=[]
    
    #Obtain the labels of the ERC containments to return it in the format to be used by Hierarchy_ERC class 
    for pair in direct_containments:
        start, end=pair
        direct_label_containments.append([f"E{start}", f"E{end}"])
    return direct_label_containments



