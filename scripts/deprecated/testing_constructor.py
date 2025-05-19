#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 21:08:29 2024

@author: tveloz
"""

import os
import sys
# Añadir el directorio raíz al PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pyCOT.reaction_network import *
from pyCOT.closure_structure import *
import networkx as nx
from pyCOT.file_manipulation import *
import matplotlib.pyplot as plt
import time

# testRN=pyCOT(SpStr,SpBt,RnStr,RnBt,RnVecS,RnVecP)
# print("specs")
# print(testRN.SpStr)
# print("specsBt")
# print(testRN.SpBt)
# print("reacsBt")
# print(testRN.RnBt)
# print("Reac")
# print(testRN.RnStr)
#Example usage:
#file_path = '../../networks/testing/in_out.txt'
#file_path = '../../networks/testing/RedPDoSR00.txt'
#file_path = '../../networks/biomodels_interesting/BIOMD0000000652_manyOrgs.xml'
#file_path = '../../networks/biomodels_interesting/Biomodel_724.txt'
file_path = '../../networks/testing/ERCs_test.txt'
file_path = '../../networks/testing/Synergy_test.txt'
file_path = '../../networks/testing/Farm.txt'
#file_path = '../../networks/testing/RedPDoSR01.txt'
#file_path = '../../networks/testing/MSORN_test1.txt'
#file_path = '../../networks/testing/MSORN_test2.txt'
#file_path = '../../networks/RandomAlife/RN_Ns_40_Norg_24_id_316.txt'
#file_path = '../../networks/biomodels_interesting/central_ecoli.txt'
#file_path = '../../networks/biomodels_all/BIOMD0000000086_url.xml'
#file_path = '../../networks/testing/MSORN_test2.txt'
#msorn= load_pyCOT_from_file(file_path)#load_pyCOT_from_Sbml(file_path)
#testRN2 = load_pyCOT_from_Sbml(file_path)
testRN = load_pyCOT_from_file(file_path)

print(testRN.RnStr)
print(testRN.get_reactions_from_species('farmer'))
rn_list=['R1','R4','R7']
print("Reacs = "+str(rn_list))

Ms=testRN.get_supp_matrix_from_reactions(rn_list)
Mp=testRN.get_prod_matrix_from_reactions(rn_list)
MS=Mp-Ms

for v in MS:
    print(v)