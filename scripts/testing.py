#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 15:28:44 2023

@author: tveloz
"""

import os
import sys
# Añadir el directorio raíz al PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# main.py (or another script)

from pyCOT.reaction_network import *
from pyCOT.closure_structure import *
import networkx as nx
from pyCOT.file_manipulation import *
import matplotlib.pyplot as plt
import time
import os

# Get the current working directory
#current_directory = os.getcwd()
#print("Current Directory:", current_directory)
#Create an instance of the HelloWorld class
# l1=["s1","s2"]
# l2=["s3","s2"]
# l=[l1,l2,["s3","s4"],["s1","s2"]]
# lp=power_list(l)
# for x in lp:
#     print(x)


# print(['s1','s2']==['s2','s1'])
# print(set(['s1','s2'])==set(['s2','s1']))

# SpStr= ['a', 'b', 'c', 'd']  # Default: ['a', 'b', 'c', 'd']
# SpBt=bt([True,True,True,True])  # Default: [a, b, c, d]
# RnStr= ['r1', 'r2']  # Updated: ['r1', 'r2']
# RnBt= bt([True, True])  # Default: [r1, r2]
# RnVecS = np.array([[1, 0, 0, 0], [0, 1, 0, 0]])  # Updated for reactions
# RnVecP= np.array([[0, 2, 0, 0], [0, 1, 1, 0]])  # Updated for reactions

file_path = r'C:\Users\tvelo\Dropbox\Public\AcademicWork\Europe\CLEA\Postdocs\TempletonPostdoc\sftw\pyCOT\networks\testing\Farm.txt'

testRN=load_pyCOT_from_file(file_path)
print("specs")
print(testRN.SpStr)
reactive_sorgs=reactive_semi_orgs(testRN)
print("printing results!")
# for sorgs in reactive_sorgs:
#    print(sorgs)

#print(testRN.get_reactions_from_bt(bt('0100')))
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

#file_path = '../../networks/testing/RedPDoSR00.txt'
#file_path = '../networks/biomodels_interesting/BIOMD0000000652_manyOrgs.xml'
#file_path = '../networks/biomodels_interesting/Biomodel_724.txt'
# file_path = '../networks/testing/ERCs_test.txt'
# file_path = '../networks/testing/Synergy_test.txt'
#file_path = '../networks/testing/Farm.txt'
# #file_path = '../networks/testing/RedPDoSR01.txt'
# #file_path = '../../networks/testing/MSORN_test1.txt'
# #file_path = '../../networks/testing/MSORN_test2.txt'
#file_path = '../networks/RandomAlife/RN_Ns_40_Norg_24_id_316.txt'
#file_path = '../networks/biomodels_interesting/central_ecoli.txt'
# file_path = '../networks/biomodels_all/BIOMD0000000086_url.xml'
#file_path = '../../networks/testing/MSORN_test2.txt'
#msorn= load_pyCOT_from_file(file_path)#load_pyCOT_from_Sbml(file_path)
#testRN2 = load_pyCOT_from_Sbml(file_path)
#testRN2 = load_pyCOT_from_file(file_path)
#print(str(testRN2.SpStr))
#print(str(testRN2.get_reactions_consuming_species("d")))
#print(str(msorn.get_connected_species_to_species("E4")))



