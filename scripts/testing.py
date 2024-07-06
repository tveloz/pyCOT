#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 15:28:44 2023

@author: tveloz
"""

# main.py (or another script)

from pyCOT.pyCOT_constructor import *
from pyCOT.pyCOT_closure_structure import *
import networkx as nx
from pyCOT.File_manipulation import *
import matplotlib.pyplot as plt
import time

#Create an instance of the HelloWorld class
SpStr= ['a', 'b', 'c', 'd']  # Default: ['a', 'b', 'c', 'd']
SpBt=bt([True,True,True,True])  # Default: [a, b, c, d]
RnStr= ['r1', 'r2']  # Updated: ['r1', 'r2']
RnBt= bt([True, True])  # Default: [r1, r2]
RnVecS = np.array([[1, 0, 0, 0], [0, 1, 0, 0]])  # Updated for reactions
RnVecP= np.array([[0, 2, 0, 0], [0, 1, 1, 0]])  # Updated for reactions

file_path = '/home/tveloz/Dropbox/Public/AcademicWork/Europe/CLEA/Postdocs/TempletonPostdoc/sftw/pyCOT/networks/testing/autocat_0_no_R.txt'

testRN=load_pyCOT_from_file(file_path)
print(testRN.get_reactions_from_bt(bt('0100')))
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



