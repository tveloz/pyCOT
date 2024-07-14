#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 29 18:55:14 2023

@author: tveloz
"""
# import re
# import random as rm
# from itertools import combinations

import numpy as np
# import pandas as pd
# from bitarray import bitarray as bt
# from bitarray import frozenbitarray as fbt
# from scipy.optimize import linprog
# import matplotlib.pyplot as plt
# from pyvis.network import Network
# import networkx as nx

from Display import *
from pyCOT_constructor import *
from itertools import chain
from File_manipulation import *


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
###ERCs produces a list of triplets:
    #0th coordinate stores the list of minimal generators for a given ERC (sets of species)
    #1st coordinate stores the ERC
    #3rd coordinate stores the reactions that generate same closure

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
            if List_closures[i]==ERC[j][1]:
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
                 #In case the new_min_gen never changed to false, we know its minimal generator to keep track:
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

def get_synergies(RN, ERCs):
    synergies=[]
    ERC=ERCs.copy()
    for i in range(len(ERC)):
        erc_i=ERC[i]
        for j in range(i,len(ERC)):
            erc_j=ERC[j]
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
    print("Creating SORN with the following ERCs")
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

#****************Compute Closed sets from terminal nodes of MSORN

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

            
    



#***************Algorithms using the pyCOT object turned into a graph
# DEPRECATED

def find_closed_sets(graph):
    closed_sets = []
    visited = set()

    def dfs_species(species, current_set):
        current_set.add(species)
        visited.add(species)

        for reaction in graph.get_reactions_for_species(species):
            if all(reactant in current_set for reactant in graph.get_reactants_for_reaction(reaction)):
                for product in graph.get_products_for_reaction(reaction):
                    dfs_species(product, current_set)

        # Check if no new reactions can be activated
        if all(reactant in current_set for reactant in graph.get_reactants_for_species(species)):
            closed_sets.append(current_set.copy())
        current_set.remove(species)

    for species in graph.species_vertices():
        if species not in visited:
            dfs_species(species, set())

    return closed_sets








def lists_equivalent(list1, list2):
    Cond1=(len(list1) == len(list2)) and (all(sublist in list2 for sublist in list1))
    Cond2=list1[-1]==list2[-1]
    return Cond1 and Cond2 

def backward_chain_closed(msorn,erc):
    "computes the backward chain of ERCs from erc"
    condition=True
    pending_paths=[[erc]]
    finished_paths=[]
    back=[]
    new=[]
    pending_reacs=msorn.RnStr
    while condition: # se agregan paths           
        for path in pending_paths: # buscar si hay algo que agregar a paths
            #print("working on "+str(path))               
            tail=path[-1] #we identify the last element in the path
            candidate_reacs=msorn.get_reactions_producing_species(tail)#get all reactions that produce the tail
            #print("candidate reac "+str(candidate_reacs))
            for r in candidate_reacs:
                support=msorn.get_supp_from_reactions(r)#get the support of the reaction
                #print("candidate reac "+r+" has supp "+str(support))                    
                new=[]
                if len(support)==1:#check that it keeps the chain
                    condition=True #there is at least one new block to be added
                    new.append(support[0])#block will be created
                    #print("current new "+str(new))
                    for s in new:#time to create and add new blocks
                        #print("adding "+s+" to new block")
                        new_block=path.copy()#block creation
                        new_block.append(s)#adding s to new_block
                        #print("new block "+str(new_block))
                        pending_paths.append(new_block)#block addition                        
            finished_paths.append(path)
            #print("finished paths "+str(finished_paths))
            pending_paths.remove(path)
            #print("pending paths "+str(pending_paths))
            
            if len(pending_paths)==0:
                condition=False 
    return finished_paths


def Closed_paths_backward(msorn):
    print("Computing closed paths of MSORN")
    '''paths is the main computing variable:
        paths[0] is a list storing each closed set as a list of ERCs. 
        paths[0] is a list storing the list of ERCs that still need to be visited
    '''
    paths=[]
    ERC=msorn.SpStr    
    terminal_ERCs=terminal_species(msorn)
    if terminal_ERCs==[]:
        terminal_ERCs.append(ERC[0])
    for i in range(len(terminal_ERCs)):
        erc=terminal_ERCs[i]
        print(str(erc))
        backward_i=msorn.get_backward_connected_species_to_species([erc])
        backward_i.remove(erc)
        pair=([[erc]],backward_i)
        paths.append(pair)
    print("before computing closed paths, we have the paths and pending")
    for i in range(len(paths)):
        print(str(paths[i][0])+" has to for "+ str(len(paths[i][1]))+" nodes")
    total=0
    for i in range(len(paths)):
        total=total+len(paths[i][1])
    itera=0    
    closed_paths=[]
    print("terminated")
    while itera<30:
        #print("********** NEW ITERATION **********")        
        print("itera = "+str(itera)+", |total paths| = "+str(len(paths))+" |finished_paths| = "+str(len(closed_paths)))
        for i in range(len(paths)):
            #This conditions checks whether the path still requires to add more nodes
            if len(paths[i][1])==0:
                None
                #print("Path "+str(paths[i])+" completed")
            else:
                #print("Working on path "+str(i)+ "="+str(paths[i][0]))
                #print("w.r.t pending nodes to add "+str(paths[i][1]))                
                #tail is the last ERC (can be a single element or a pair of ERCs) added to the path
                tail=paths[i][0][-1].copy()
                #print("tail is "+str(tail))
                if len(tail)==0:
                    #print("Empty tail!")
                    None
                else:
                    if not isinstance(tail,list):
                        tail=[tail]  
                        #print("tail turned into a list!")
                    for j in range(len(paths[i][0])-1):    
                        for s in paths[i][0][j]:
                            if s in tail:
                                tail.remove(s)
                                #print("reactant "+str(s)+" removed from tail because already in "+str(paths[i][0][j])) 
                    reactions_tail=[]
                    for k in range(len(tail)):
                        reactions_tail=msorn.get_reactions_producing_species(tail[k])
                        #print("for tail "+str(tail[k])+" the reactions_tail is "+str(reactions_tail))                       
                        to_add=False
                        for r in reactions_tail:
                            reactants=msorn.get_supp_from_reactions([r])
                            if not isinstance(reactants,list):
                                reactants=[reactants]
                            #print("for reaction "+str(r)+" the reactants are "+str(reactants))                                
                            for s in reactants:
                                if s in paths[i][1]:
                                    to_add=True
                            if to_add:
                                new_path=paths[i][0].copy()
                                new_path.append(reactants)
                                to_explore=paths[i][1].copy()
                                for s in reactants:
                                    if s in to_explore:
                                        to_explore.remove(s)
                                paths.append([new_path,to_explore])
                                #print("new path = ["+str(new_path)+","+str(to_explore)+"]")
                #print("End of path "+str(i))
                paths[i]=[paths[i][0],[]]
                closed_paths.append(paths[i][0])
                
                #print("will remain now as")
                #print(paths[i])
        print("********* END OF ITERATION******************")
        #print("Path removal")
        #print("Total number before removal = "+str(len(paths)))
        paths=remove_duplicates_pairs2(paths)    
        #print("Total number after removal = "+str(len(paths)))
        total=0
        for i in range(len(paths)):
            total=total+len(paths[i][1])
            #print("path "+str(i)+" = " +str(paths[i]))
        print("Current number of paths = "+str(len(paths)))
        print("Total number of pending nodes = "+str(total))
        print("Iteration Number = "+str(itera))
        itera=itera+1
        if total==0:
            break
    return closed_paths


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
    print("the ERC list")
    for i in range(len(erc[0])):
        print(str([erc[0][i],erc[1][i]]))
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

def remove_duplicates_pairs(input_list):
    unique_sublists = []
    seen_sets1 = set()
    seen_sets2 = set()
    for sublist1,sublist2 in input_list:
        sublist_set1 = frozenset(flatten_list(sublist1))
        sublist_set2 = frozenset(flatten_list(sublist1))
        if (sublist_set1 not in seen_sets1) or (sublist_set2 not in seen_sets2):

            unique_sublists.append([sublist1,sublist2])
            seen_sets1.add(sublist_set1)
            seen_sets2.add(sublist_set2)

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
    print("unique removed "+str(len(input_list)-len(unique_sublists))+" elements")
    return unique_sublists