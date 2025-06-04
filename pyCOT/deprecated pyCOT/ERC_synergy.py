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

class ERC_Synergy:
    def __init__(self, reactants, product):
        """
        Initialize a Synergy object.
        :param reactants: List of ERC objects representing the reactants.
        :param product: An ERC object representing the product.
        """
        if not isinstance(reactants, list) or not all(isinstance(r, ERC) for r in reactants):
            raise ValueError("Reactants must be a list of ERC objects.")
        if not isinstance(product, ERC):
            raise ValueError("Product must be an ERC object.")
        
        self.reactants = reactants
        self.product = product
        self.rlabel=[r.label for r in reactants]
        self.plabel=product.label

    def __repr__(self):
        """
        Return a string representation of the Synergy object.
        """
        reactant_ids = [f"E{idx}" for idx, _ in enumerate(self.reactants)]
        product_id = f"E{self.product}"
        return f"Synergy(Reactants={reactant_ids}, Product={product_id})"
def Synergy(erc1,erc2,hierarchy,RN):
    """
        Returns the list of direct synergy (highest-level possible ERCs generated in the hierarchy) out of erc1 and erc2, if it exists.

        :param erc1,erc2: ercs to combine
        :param hierarchy: An Hierarchy_ERC object
        :param RN: The underlying RN
        """
    synergies=[]
    if erc1 in hierarchy.get_contain(erc2) or erc2 in hierarchy.get_contain(erc1):
        return synergies
    else:
        erc_1=erc1.get_closure(RN)
        erc_2=erc2.get_closure(RN)
        rn1=RN.get_reactions_from_species(erc_1)
        rn2=RN.get_reactions_from_species(erc_2)
        union_reacs=list(set(rn1).union(set(rn2)))            
        joint_closure=closure(RN,list(set(erc_1).union(erc_2)))          
        
        joint_closure_reacs=RN.get_reactions_from_species(joint_closure)
        if len(joint_closure_reacs)>len(union_reacs):   
                #print("jcr= "+str(joint_closure_reacs)+", union="+str(union_reacs))
                novel_reacs= [r for r in joint_closure_reacs if r not in union_reacs]  
                for r in novel_reacs:
                    syn=hierarchy.get_erc_from_reaction(RN,hierarchy,r)
                    add = True
                    for s in synergies:
                        #print("testing that "+str(syn.label)+" in "+str(hierarchy.get_contained(s.product)))
                        if syn in hierarchy.get_contained(s.product):
                            add=False
                            break
                        #print("testing that "+str(s.rlabel)+" is not "+str([erc1.label, erc2.label])+" and "+s.plabel+" is "+syn.label)
                        if set(s.rlabel) == set([erc1.label, erc2.label]) and s.plabel == syn.label:  
                            add=False
                            break
                    if add:
                        synergies.append(ERC_Synergy([erc1,erc2],syn))
                #Eliminate the still possibly redundant synergies, i.e 
                for i in range(len(synergies)):
                    for j in range(i,len(synergies)):
                        s1=synergies[i]
                        s2=synergies[j]
                        if s2.product in hierarchy.get_contained(s1.product):
                            print(s2.product+" contained in "+s1.product)
                            synergies.remove(s2)
                        elif s1.product in hierarchy.get_contained(s2.product):
                            print(s1.product+" contained in "+s2.product)
                            synergies.remove(s1)
                return synergies
        
def non_redundant_synergies(synergies,hierarchy):
    print("starting non-redudant-synergies")
    nr_synergies=synergies.copy()
    for i in range(len(nr_synergies)):
        for j in range(i+1,len(nr_synergies)):
            #Define ERCs in place
            s1=nr_synergies[i]
            s2=nr_synergies[j]
            s1r0=s1.reactants[0].label
            s1r1=s1.reactants[1].label
            s2r0=s2.reactants[0].label
            s2r1=s2.reactants[1].label
            #Obtain Containment matrix (possibly input, must have diagonal equals to 1)
            #CM[ij]
            #Obtain relevant containment coordinates
            #s1 =Ei+Ej->Ek and s2=Ei'+Ej'->Ek'
            #if ((CM[i,i'] and CM[j,j']) or (CM[i,j'] and CM[j,i'])) and CM[k',k]):  
            # Then reactants of s1 are same or lower than reactants of s2 and product of s1 are same or higher than product of s2 
            #eliminate s2
            #if ((CM[i',i] and CM[j',j]) or (CM[i',j] and CM[j',i])) and CM[k,k']):  
            # Then reactants of s2 are same or lower than reactants of s1 and product of s2 are same or higher than product of s1 
            #eliminate s1

            #Eliminate s2 if reactants of s1 are lower
            print("is "+s2r0+"+"+s2r1+"->"+s2.product.label)
            print("reduntant vs "+s1r0+"+"+s1r1+"->"+s1.product.label) 
            if ((s1r0 in hs2r0 and s1r1 in hs2r1 or (s1r1 in hs2r0 and s1r0 in hs2r1)) and s2.product in hs1p): 
                print("redundant!")
                nr_synergies.remove(s2)
            if ((s2r0 in hs1r0 and s2r1 in hs1r1 or (s2r1 in hs1r0 and s2r0 in hs1r1)) and s1.product in hs2p): 
                print("redundant!")
                nr_synergies.remove(s1)  
    return nr_synergies


# def get_all_synergies():
#      all_syn=[]
# for erc1 in hierarchy.ercs:
#     for erc2 in hierarchy.ercs:
#         #print("checking if there is a synergy between "+erc1.label+" and "+erc2.label)
#         syn=Synergy(erc1,erc2,hierarchy,RN)
#         #print(str(syn))
#         if syn!=None:
#             for s in syn:
#                 add=True
#                 #print(s.reactants[0].label+"+"+s.reactants[1].label+"->"+s.product.label)
#                 #print("testing that "+str(s.rlabel)+" is not "+str([erc1.label, erc2.label])+" and "+s.plabel+" is "+syn.label)
#                 for x in all_syn:
#                     if set(s.rlabel) == set(x.rlabel) and s.plabel == x.plabel:  
#                         add=False
#                         break
#                 if add:
#                     all_syn.append(s)