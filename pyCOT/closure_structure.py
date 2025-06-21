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
# from pyCOT.file_manipulation import *
# from pyCOT.display import *
#from pyCOT.reaction_network import *
import networkx as nx
# from networkx.drawing.nx_agraph import graphviz_layout
# from pyCOT.simulations import *  
from pyCOT.rn_types import *
# from pyCOT.rn_rustworkx import *
from pyCOT.rn_rustworkx import ReactionNetwork

###########FIRST APPROACH: COMPUTING ALL REACTIVE SEMIORGS BY BRUTE FORCE ON ERCS###################
def remove_duplicate_sublists(list_of_lists):
    seen = set()
    unique_list = []
    
    for sublist in list_of_lists:
        # Convert sublist to tuple to make it hashable
        tuple_sublist = tuple(sublist)
        if tuple_sublist not in seen:
            seen.add(tuple_sublist)
            unique_list.append(sublist)
    
    return unique_list

def list_of_lists_to_set(list_of_lists):
    return set(element for sublist in list_of_lists for element in sublist)

def power_list(list_of_sublists):
  # Convert each sublist to a set to avoid duplicate elements within sublists
    list_of_sets = [set(sublist) for sublist in list_of_sublists]
    
    # Generate all combinations of these sets
    combination_cases = list(chain.from_iterable(combinations(list_of_sets, r) for r in range(len(list_of_sets) + 1)))
    
    result = []
    for comb in combination_cases:
        # Merge sets within each combination to avoid duplicates across sublists
        merged_set = set().union(*comb)
        # Add as a sorted list to maintain consistency
        result.append(sorted(merged_set))
    
    # Remove duplicate lists by converting to tuples for hashable unique entries
    unique_result = [list(x) for x in set(tuple(x) for x in result)]
    return unique_result

def reactive_semi_orgs(RN):
     print("reactive_semi_orgs starting")
     ERC_list=ERCs(RN)
     ERCs_closures=[sublist[1] for sublist in ERC_list if len(sublist) > 1]
     for erc in ERC_list:
         print(erc)
     #print("###computing power sets###")
     power_list_ERCs=power_list(ERCs_closures)
     #print(power_list_ERCs)
     #print("###computing power set closures###")
     power_list_closures=remove_duplicate_sublists(closures(RN,power_list_ERCs))
     #print(power_list_closures)
    
     reactive_sorgs= [sorg for sorg in power_list_closures if RN.is_semi_self_maintaining(sorg)]
     print("reactive_semi_orgs ending")
     print("number of ERCs : "+str(len(ERCs_closures)))
     print("size of power list of ERCs with duplicates: "+str(len(power_list_ERCs)))
     print("closures found : "+str(len(power_list_closures)))
     print("semiorgs found: "+str(len(reactive_sorgs)))
     return(reactive_sorgs)     
###################END OF FIRST APPROACH: COMPUTING ALL REACTIVE SEMIORGS BY BRUTE FORCE ON ERCS###################


###################Second approach: Computing relevant semi-organisations using dynamical hierarchy algorithm



##################Calcular Self-maintainance####################################################
# Function for solve the Linear Programming Problem: S.v>=0, v>0
def minimize_sv(S, epsilon=1,method='highs'):
    """
    Minimizes c @ v subject to S @ v >= 0 and v > 0.
    
    Parameters:
    ----------
    S : numpy.ndarray
        Stoichiometric matrix of the reaction network (RN).
    epsilon : float, optional
        Small positive value to ensure all coordinates of v remain strictly positive (default is 1).

    Returns:
    -------
    numpy.ndarray
        Optimal process vector v with all strictly positive coordinates.

    Raises:
    ------
    ValueError
        If no feasible solution is found.
    """
    n_species, n_reactions = S.shape  # Dimensions of S
    c = np.ones(n_reactions)          # Objective function: minimize the sum of v
    A_ub = -S                         # Reformulate S @ v >= 0 as -S @ v <= 0
    b_ub = np.zeros(n_species)        # Inequality constraints
    bounds = [(epsilon, None)] * n_reactions  # v > 0 (avoiding exact zero values) 
    
    # Solve the linear programming problem: minimize c @ v subject to (A_ub @ v <= b_ub)
    result = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=bounds, 
                     method=method) # 'highs' uses the Dual Simplex method (good for many constraints)
                                     # 'highs-ipm' uses the Interior Point method (faster for large sparse problems)
        
    if result.success:
        return [True, result.x]
    else:
        print("minimize_sv: No feasible solution was found.")
        return [False, np.zeros(n_reactions)]

def is_self_maintaining(RN: ReactionNetwork,X):
    """
    Verifies X is a self-maintaining set w.r.t RN
    
    Parameters:
    ----------
    X : list of species
    RN: Reaction network (RN).
    
    Returns:
    -------
    [Boolean,np.array1,np.array2] 
    Boolean expresses if it is SelfMaintaining, 
    np.array1 delivers self-maintaining vector,
    np.array2 delivers production vector      
    """
    #matrix_data = universal_stoichiometric_matrix(RN) # Calculates the universal stoichiometric matrix associated with the reaction network
    # print("reactions to build S")
    # print(RN.get_reactions_from_species(X))
    sub_RN=RN.sub_reaction_network(X) # Creates a sub-reaction network with the required set of species and reactions
    S=sub_RN.stoichiometry_matrix()
    # print("Stoichiometric Matrix")
    # print(S)
    res=minimize_sv(S, epsilon=1,method='highs')
    xv=S @ res[1]
    res.append(xv)
    return res


#############################################################
def remove_duplicates_list(semi_org):
    """
    Remove duplicate sublists from a list of lists.

    Args:
        semi_org (list of list): A list of subsets (lists) that may contain duplicates.

    Returns:
        list of list: The input list with duplicate subsets removed, preserving original order.
    """
    unique_subsets = []
    seen_sets = []  # List to store already seen sets

    for sublist in semi_org:
        current_set = frozenset(sublist)  # Use frozenset to compare contents
        if current_set not in seen_sets:
            seen_sets.append(current_set)
            unique_subsets.append(sublist)  # Keep the original list
 
    semi_organisations = []
    for i, semi_org_new in enumerate(unique_subsets):
        semi_organisations.append(semi_org_new)

    return semi_organisations

#############################################################
# Function to find organisations in a reaction network 
def find_organisations(rn, semi_orgs):
    """
    Find organisations of a reaction network from the list of semi-organisations and check if they are self-maintaining.

    Args:
        rn: A reaction network object or structure required by `is_self_maintaining`.
        semi_orgs (list of list): A list of candidate semi-organisations.

    Returns:
        tuple:
            organisations (list of list): The self-maintaining organisations.
            vector_process (list of list): The corresponding process vectors.
            vector_production (list of list): The corresponding state vectors.
    """
    semi_orgs = remove_duplicates_list(semi_orgs)  # Remove duplicates
    semi_orgs.sort(key=lambda x: len(x))  # Sort by size
    print("Semi-organisations without duplicates and sorted by size:")
    for i, semi_org_new in enumerate(semi_orgs):
        print(f"S{i+1}:", semi_org_new)
    print("-" * 70)

    print("find_organisations starting")
    organisations = []
    vector_process = []
    vector_production = []
    for i, semi_organisation in enumerate(semi_orgs):
        if len(semi_organisation) == 0:
            print("Semi-organisation = []\n")
            print("     Empty semi-organisation.")
            continue

        print(f"\nSemi-organisation_{i+1} = {semi_organisation}")
        print(f"Semi-organisation_{i+1} size =", len(semi_organisation))
        res = is_self_maintaining(rn, X=semi_organisation)
        if res[0] == True:
            print("     Is self-maintaining:", res[0])
            organisations.append(semi_organisation)
            vector_process.append([x.tolist() for x in res[1]])
            vector_production.append([x.tolist() for x in res[2]])
        else:
            print("     Is self-maintaining: False")

    return organisations, vector_process, vector_production


################## END OF DISGRESSION ############################################


