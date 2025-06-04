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


#def closure_structure(RN,):
    #get generators
    #compute closures     
    #get minimal generators
    #get containments
    #get synergies
    #create second order RN    
    #obtain all pathways from the second order RN
    #create hierarchy
#    return 0

#Obtain generators of a reaction network
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
        ERC[i].append(i)
    return ERC


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

def get_direct_containments(RN,ERC,containment_list):
    containments=containment_list.copy()
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

def get_ERC_from_reaction(ERC,r):
    for erc in ERC:
        if r in erc[2]:
            return erc

########## Chain characterization in terms of levels #############
def build_chains(set_labels, containment_pairs):
    # Step 1: Build the adjacency list (graph)
    containment_graph = defaultdict(list)
    reverse_graph = defaultdict(list)
    
    for a, b in containment_pairs:
        containment_graph[a].append(b)
        reverse_graph[b].append(a)

    # Step 2: Identify root nodes (nodes that are not children of any other node)
    root_nodes = set(set_labels) - set(reverse_graph.keys())

    # Step 3: Function to recursively find chains via DFS
    def find_chains(node, current_chain, all_chains):
        current_chain.append(node)
        # If node has no children, it is the end of a chain
        if node not in containment_graph or not containment_graph[node]:
            all_chains.append(current_chain[:])  # Add a copy of the current chain
        else:
            # Recursively explore each child node
            for child in containment_graph[node]:
                find_chains(child, current_chain, all_chains)
        # Backtrack
        current_chain.pop()

    # Step 4: Find all chains starting from each root node
    all_chains = []
    for root in root_nodes:
        find_chains(root, [], all_chains)
    
    return all_chains, containment_graph

def plot_containment_forest(set_labels, containment_pairs):
    all_chains, containment_graph = build_chains(set_labels, containment_pairs)
    
    # Create a directed graph
    G = nx.DiGraph()
    
    # Add edges to the graph from the containment pairs
    G.add_edges_from(containment_pairs)
    
    # Generate positions using graphviz layout (hierarchical)
    pos = nx.nx_agraph.graphviz_layout(G, prog="dot")  # 'dot' layout is good for hierarchical graphs
    
    # Draw the graph
    plt.figure(figsize=(10, 10))
    nx.draw(G, pos, with_labels=True, node_size=2000, node_color="skyblue", font_size=10, font_weight="bold", arrows=True)
    
    plt.title("Containment Forest")
    plt.show()

# Calculate node levels
def calculate_levels(graph):
    levels = defaultdict(int)  # To store the level of each node
    visited = set()            # Keep track of visited nodes
    queue = [(node, 0) for node in graph if node not in reverse_graph]  # Start with root nodes (level 0)

    while queue:
        current, level = queue.pop(0)
        if current in visited:
            continue
        visited.add(current)
        levels[level] += 1  # Increment the count of nodes at this level

        # Add children to the queue with the next level
        for child in graph.get(current, []):
            queue.append((child, level + 1))

    return levels








def get_synergies(RN, ERCs):
    synergies=[]
    ERC=ERCs.copy()
    #Can this exploration be optimized? 
    for i in range(len(ERC)):
        erc_i=ERC[i]
        for j in range(i,len(ERC)):
            erc_j=ERC[j]
            #print("get_synergies, iteration "+str(i)+","+str(j))
            if not (set_containment(erc_i,erc_j) or set_containment(erc_j,erc_i)):
                compound=list(set(erc_i[1]).union(erc_j[1]))
                compound_reacs=RN.get_reactions_from_species(compound)
                union_reacs=list(set(erc_i[2]).union(erc_j[2]))                
                if len(compound_reacs)>len(union_reacs):
                    novel_reacs= [r for r in compound_reacs if r not in union_reacs]
                    for r in novel_reacs:
                        syn_index=get_ERC_from_reaction(ERC,r)[3]
                        if [i,j,syn_index] not in synergies:
                            synergies.append([i,j,syn_index])                 
    return synergies
def get_direct_synergies(RN,ERCs,direct_containment_list,synergies_list):
    direct_synergies = []
    containments = direct_containment_list.copy()
    #indexes to be eliminated must be start and end in two different transitions
    flattened_list = list(chain.from_iterable(containments))
    erc_dubious_indexes = list(set(flattened_list))
    
    for triplet in synergies_list:
        i, j, k = triplet
        to_add=True
        # Check if there is no intermediate element between start and end
        for mid in erc_dubious_indexes:
            if ([i, mid] in containments) and ([mid,j,k] in synergies_list or [j,mid,k] in synergies_list):
                #print("won't add "+str([i,j,k])+" beacuse "+str([i,mid])+" and "+str([mid,j,k])+" or "+str([j,mid,k]))
                to_add=False
                break
            elif ([j, mid] in containments) and ([mid,i,k] in synergies_list or [i,mid,k] in synergies_list):
                #print("won't add "+str([i,j,k])+" beacuse "+str([j,mid])+" and "+str([mid,i,k])+" or "+str([i,mid,k]))
                to_add=False
                break
            elif ([mid, k] in containments) and ([i,j,mid] in synergies_list or [j,i,mid] in synergies_list):
                #print("won't add "+str([i,j,k])+" beacuse "+str([k,mid])+" and "+str([i,j,mid])+" or "+str([j,i,mid]))
                to_add=False
                break
            elif [i,j,mid] in synergies_list and ([i,mid,k] in synergies_list or [mid,j,k] in synergies_list):
                #print("synergy "+str([i,j,k])+ "is not direct because "+str([i,j,mid])+" and "+str([i,mid,k])+" or "+str([mid,j,k]))
                to_add=False
                break
        if to_add:
            direct_synergies.append(triplet)    
    return direct_synergies


def second_order_network(RN,ERC,ERC_containments,ERC_synergies):
    ERC=ERCs(RN)
    erc=[]
    sorn_species=[]
    ERC_generative_str=ERC_containments.copy()
    for r in ERC_synergies:
        ERC_generative_str.append(r)
    for i in range(len(ERC)):
        sorn_species.append("E"+str(i))
        erc.append(ERC[i][1])
        
    sorn_reactions=[]
    for l in range(len(ERC_containments)):
        i=ERC_containments[l][0]
        j=ERC_containments[l][1]
        sorn_reactions.append("E"+str(i)+"=>"+"E"+str(j))
    
    for l in range(len(ERC_synergies)):
            
        i=ERC_synergies[l][0]
        j=ERC_synergies[l][1]
        k=ERC_synergies[l][2]
        sorn_reactions.append("E"+str(i)+"+"+"E"+str(j)+"=>"+"E"+str(k))
    # for l in range(len(ERC)):
    #     sorn_reactions.append("E"+str(l)+"=>"+"E"+str(l))                  
    return [[sorn_species,erc],sorn_reactions, ERC_generative_str]        

def MSORN(sorn):
    print("Creating minimal SORN")
    if len(sorn[1])>0:
        to_file=str(sorn[1])    
        to_file=to_file.replace("'","")
        to_file=to_file.replace("[","")
        to_file=to_file.replace("]","")
        to_file=to_file.replace(", ", "\n")
        #print(to_file)
        file_path = "MSORN.txt"
        with open(file_path, 'w') as file:
        # Write the string to the file
            for i, line in enumerate(to_file.split("\n"), start=0):
            # Add a prefix that includes the line number
                modified_line = f"r{i}: {line};\n"
                file.write(modified_line)
        msorn=load_pyCOT_from_file(file_path)        
        return msorn 
    else:
        ERC_species=sorn[0][0]
        ERC_Bt = bitarray(len(ERC_species))
        ERC_Bt.setall(True)
        return ReactionNetwork(SpBt=ERC_Bt, RnStr=[], RnBt= bitarray(), SpStr=ERC_species, RnMsupp=np.array([]),RnMprod=np.array([]))



#****************Compute Closed sets from terminal nodes of MSORN - Backward approach*********************

#Compute the ERCs that do not in the support of any reaction
def terminal_species(msorn):
    terminal=msorn.SpStr.copy()
    for r in msorn.RnStr:
        s=msorn.get_supp_from_reactions(r)
        if len(s)==1:
            sp=s[0]
            #print("support "+str(sp)+" is in "+str(r))
            if sp in terminal:
                terminal.remove(sp)
            #print("current terminal:")
            #print(str(terminal))
    return(terminal)

            
    
def Closed_paths_backward(msorn):
    print("Computing closed paths of MSORN")
    '''open_paths is the list of paths to be closed:
        open_paths[i,0] is a chain of ERCs that where open_paths[i,0,j+1] is the support of open_path[i,0,j]. 
        paths[i,1] keeps the last element to be generated in open_paths[i,0], as open_paths[i,0,-1] can contain a pair of ERCs. 
    Idea: The algorithm checks the MSORN backwards, storing each step as a closed_path. It keeps information about the last node to be generated. 
          If such last node is empty path is terminated, if last node is already in the open_path it means there is a loop, so the open_path is not built further.
    '''
    open_paths=[]
    ERC=msorn.SpStr    
    terminal_ERCs=terminal_species(msorn)
    if terminal_ERCs==[]:
        terminal_ERCs.append(ERC[0])
    for i in range(len(terminal_ERCs)):
        erc=terminal_ERCs[i]
        #print(str(erc))
        open_paths.append([[[erc]],[erc]])
    total=len(open_paths)
    itera=0    
    closed_paths=[]
    #print("terminated")
    while itera<10:
        paths_to_delete=[]
        new_open_paths=[]
        print("********** NEW ITERATION **********")        
        print("itera = "+str(itera)+", |total open paths| = "+str(len(open_paths))+" |finished closed_paths| = "+str(len(closed_paths)))
        for i in range(len(open_paths)):
            #This conditions checks whether the path still requires to add more nodes
            #print("working on open_path "+str(i)+" out of "+str(len(open_paths)))
            #print(str(open_paths[i]))
            if len(open_paths[i][1])==0:
                None
                #print("Path "+str(open_paths[i])+" completed")
            else:
                #print("Working on path "+str(i)+ "="+str(open_paths[i]))
                #tail is the last ERC (can be a single element or a pair of ERCs) added to the path
                tail=open_paths[i][1]          
                reactions_tail=[]
                reactions_tail=msorn.get_reactions_producing_species(tail)
                #print("for tail "+str(tail[k])+" the reactions_tail is "+str(reactions_tail))                       
                if len(reactions_tail)==0:
                    paths_to_delete.append(open_paths[i])
                else:
                    for r in reactions_tail:
                        reactants=msorn.get_supp_from_reactions([r])
                        if not isinstance(reactants,list):
                            reactants=[reactants]
                        #print("for reaction "+str(r)+" the reactants are "+str(reactants))                                
                        for s in reactants:
                            if s in flatten_list(open_paths[i][0]):
                                to_add=False
                                #print("species "+s+" already in "+str(open_paths[i][0]))
                                paths_to_delete.append(open_paths[i])
                            else:
                                to_add=True
                                #print("species "+s+" is new tail of "+str(open_paths[i][0]))
                            if to_add:
                                paths_to_delete.append(open_paths[i])
                                new_open_path=open_paths[i][0].copy()
                                new_open_path.append(reactants)
                                new_tail=s
                                new_open_paths.append([new_open_path,new_tail])
        paths_to_delete=remove_duplicates_pairs2(paths_to_delete)
        for path in paths_to_delete:
            closed_paths.append(path)
            open_paths=remove_sublist(open_paths,path)
        for path in new_open_paths:
            open_paths.append(path)
        #print("BEFORE DELETING DUPLICATES")
        #print("Current number of closed_paths = "+str(len(closed_paths)))
        #print("Current number of open_paths = "+str(len(open_paths)))
        #print("********* END OF ITERATION******************")
        #print("AFTER DELETING DUPLICATES")
        open_paths=remove_duplicates_pairs2(open_paths)
        closed_paths=remove_duplicates_pairs2(closed_paths)    
        #print("Current number of closed_paths = "+str(len(closed_paths)))
        #print("Current number of open_paths = "+str(len(open_paths)))
        #print("Total number of pending nodes = "+str(total))
        itera=itera+1
        if len(open_paths)==0:
            break
    return closed_paths

def remove_sublist(main_list, sublist_to_remove):
   return [sublist for sublist in main_list if sublist != sublist_to_remove]

def flatten_list(lst):
    flattened = []
    for i in lst:
        if isinstance(i, list):
            flattened.extend(flatten_list(i))
        else:
            flattened.append(i)
    return flattened


def Closed_paths_to_Sp(closed_paths,sorn):
    f_closed_paths=[]
    for i in range(len(closed_paths)):
        f_closed_paths.append(list(set(flatten_list(closed_paths[i]))))
    # Use a set to keep track of seen sublists
    seen = set()
    
    # Use a list comprehension to filter out duplicates
    erc_closed_sets= [x for x in f_closed_paths if tuple(x) not in seen and not seen.add(tuple(x))]
    erc = sorn[0]
    # print("the ERC list")
    # for i in range(len(erc[0])):
    #     print(str([erc[0][i],erc[1][i]]))
    sp_closed_sets=[[]]
    for i in range(len(erc_closed_sets)):
        #print("erc_closed "+str(i)+" = "+str(erc_closed_sets[i]))
        for E in erc_closed_sets[i]:
            for j in range(len(erc[0])):
                if E==erc[0][j]:
                    #print("Ei "+str(E)+" = "+str(erc[0][j]))                    
                    for s in erc[1][j]:   
                        sp_closed_sets[i].append(s)
            sp_closed_sets[i]=list(set(sp_closed_sets[i]))
        if i<len(erc_closed_sets)-1:
            sp_closed_sets.append([])      
    seen = set()    
    # Use a list comprehension to filter out duplicates
    unique_sp_closed_sets= [x for x in sp_closed_sets if tuple(x) not in seen and not seen.add(tuple(x))]

    return unique_sp_closed_sets

def remove_duplicates(input_list):
    unique_sublists = []
    seen_sets = set()

    for sublist in input_list:
        sublist_set = frozenset(sublist)
        if sublist_set not in seen_sets:
            unique_sublists.append(sublist)
            seen_sets.add(sublist_set)

    return unique_sublists



def remove_duplicates_pairs2(input_list):
    unique_sublists = []
    seen_list1 = []
    for i in range(len(input_list)):
        if len(seen_list1)==0:
            unique=True
        else: 
            unique=True            
            input1=frozenset(flatten_list(input_list[i][0]))
            for j in range(len(seen_list1)):                
                seen1=frozenset(flatten_list(seen_list1[j]))
                if  input1==seen1:
#                   seen_list1.append(input_list[i][0])
                    unique=False
        if unique==True:
            unique_sublists.append([input_list[i][0],input_list[i][1]])
            seen_list1.append(input_list[i][0])
    #print("unique removed "+str(len(input_list)-len(unique_sublists))+" elements")
    return unique_sublists