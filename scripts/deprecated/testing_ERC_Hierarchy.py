#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 15:28:44 2023

@author: tveloz
"""

import sys
import os

# Añadir el directorio raíz al PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pyCOT.reaction_network import ReactionNetwork
from pyCOT.ERC_Hierarchy import *
from  pyCOT.reactive_features import *
from  pyCOT.ERC_synergy import *
import networkx as nx
from pyCOT.file_manipulation import *
import matplotlib.pyplot as plt
import time
from collections import defaultdict
from collections import Counter # para el gráfico de cantidad de básicos vs repeticiones

# Create an instance of the HelloWorld class
# SpStr= ['a', 'b', 'c', 'd']  # Default: ['a', 'b', 'c', 'd']
# SpBt=bt([True,True,True,True])  # Default: [a, b, c, d]
# RnStr= ['r1', 'r2']  # Updated: ['r1', 'r2']
# RnBt= bt([True, True])  # Default: [r1, r2]
# RnVecS = np.array([[1, 0, 0, 0], [0, 1, 0, 0]])  # Updated for reactions
# RnVecP= np.array([[0, 2, 0, 0], [0, 1, 1, 0]])  # Updated for reactions

# testRN=pyCOT(SpStr,SpBt,RnStr,RnBt,RnVecS,RnVecP)
# print("specs")
# print(testRN.SpStr)
# print("specsBt")
# print(testRN.SpBt)
# print("reacsBt")
# print(testRN.RnBt)
# print("Reac")
# print(testRN.RnStr)
# print("SuppVec")
# print(testRN.RnVecS)
# print("ProdVec")
# print(testRN.RnVecP)
#Example usage:
#file_path = '../../networks/testing/in_out.txt'
# file_path = '../networks/testing/RedPDoSR01.txt'
# file_path = '../networks/testing/RN_AEP_03.txt'
# file_path = '../networks/testing/RN_AEP_04.txt'
# file_path = '../networks/testing/navarino.txt'
#file_path = 'networks/testing/ERCs_test.txt'
# file_path = '../networks/testing/RedPDoSR01.txt'
#file_path = '../networks/biomodels_all/BIOMD0000000011/BIOMD0000000011.xml'
#file_path = 'networks/biomodels_interesting/BIOMD0000000237_manyOrgs.txt'  #ERROR: ZeroDivisionError: division by zero


#file_path = 'networks/biomodels_interesting/BIOMD0000000652_manyOrgs.txt'

#file_path = '../networks/testing/ERCs_test.txt'
#file_path = '../../networks/testing/Synergy_test.txt'
# file_path = '../networks/testing/Farm.txt'
#file_path = '../networks/testing/RedPDoSR01.txt'
#file_path = 'networks/testing/MSORN_test1.txt'
file_path = 'networks/testing/autopoietic.txt'
#file_path = 'networks/Navarino/RN_IN_02_Py.COT.txt'
#file_path = 'networks/testing/Farm.txt'
#file_path = 'networks/testing/autopoietic_ext.txt'
file_path = 'networks/testing/Synergy_test.txt'

#file_path = 'networks/RandomAlife/RN_Ns_40_Norg_20_id_396.txt'
#file_path = 'networks/biomodels_interesting/central_ecoli.txt'


# Assuming load_pyCOT_from_file and reac_analysis are defined elsewhere in your code

RN = load_pyCOT_from_file(file_path)
hierarchy = Hierarchy_ERC(RN)
# Create a GEN object
# gen = GEN(RN)


for erc in hierarchy.ercs:
    print(erc.label)
    print(erc.min_generators)
    print(erc.get_closure(RN)) 
   
    print("##############")


print(hierarchy)

# for erc in hierarchy.ercs:
#     print(erc.label+ " is above "+str(hierarchy.get_contained(erc)))
#     print(erc.label+ " is below "+str(hierarchy.get_contain(erc)))
#erc1=hierarchy.ercs[1]
#erc2=hierarchy.ercs[2]

all_syn=[]
for erc1 in hierarchy.ercs:
    for erc2 in hierarchy.ercs:
        #print("checking if there is a synergy between "+erc1.label+" and "+erc2.label)
        syn=Synergy(erc1,erc2,hierarchy,RN)
        #print(str(syn))
        if syn!=None:
            for s in syn:
                add=True
                #print(s.reactants[0].label+"+"+s.reactants[1].label+"->"+s.product.label)
                #print("testing that "+str(s.rlabel)+" is not "+str([erc1.label, erc2.label])+" and "+s.plabel+" is "+syn.label)
                for x in all_syn:
                    if set(s.rlabel) == set(x.rlabel) and s.plabel == x.plabel:  
                        add=False
                        break
                if add:
                    all_syn.append(s)

all_syn=non_redundant_synergies(all_syn,hierarchy)

for s in all_syn:
    print(s.reactants[0].label+"+"+s.reactants[1].label+"->"+s.product.label)

# Get the level of a specific ERC
# for i, erc in enumerate(gen.ERCs):
#     level = gen.get_erc_level(i)  # Replace 0 with the desired ERC index
#     print(f"Level of ERC= E"+str(erc.species)+" : " +str(level))

# target = 3
# level = gen.get_erc_level(target)  # Replace 0 with the desired ERC index
# print(f"Level of target ERC= E"+str(erc.species)+" : " +str(level))

# result = get_higher_ercs(gen, target)
# print(f"ERCs that transitively point to node {target}: {result}")

