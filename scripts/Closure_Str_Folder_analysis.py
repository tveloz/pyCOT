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
from pyCOT.closure_structure import *
from  pyCOT.reactive_features import *
import networkx as nx
from pyCOT.file_manipulation import *
import matplotlib.pyplot as plt
import time
from collections import defaultdict
from collections import Counter # para el gráfico de cantidad de básicos vs repeticiones
from pathlib import Path
import pickle
import glob
from pathlib import Path
from Closure_Str_RN_analysis import Closure_Str_RN_analysis 

# This script aims at performing the analysis made in Closure_Str_RN_analysis.py in a systematic way. It reads a folder of
# RN in .txt format. The script can be executed in multiple runs as it stores everything in .pkl files, one for each
# .txt RN and performs a check to identify where should start to continue computing the missing parts.

#folder_path = 'networks/Navarino'
#txt_files = glob.glob(os.path.join(folder_path, '*.txt'))
input_path = Path("networks/Navarino/RN_IN_02_PyCOT.txt")
input_path = Path('networks/biomodels_interesting/BIOMD0000000652_manyOrgs.txt')
output_file = input_path.with_name(input_path.stem + '.pkl')    
# Open the pickle file in read-binary mode
objects=[]
print("checking if "+str(output_file)+" exists")

if output_file.exists():
    print(str(output_file)+" already exists")
    with open(output_file, 'rb') as file:
        while True:
            try:
                # Attempt to load the next object
                obj = pickle.load(file)
                objects.append(obj)
            except EOFError:
                # End of file reached, break out of the loop
                break
        if len(objects)==21:
            print("file "+str(output_file)+ " already analyzed")
        else: 
            print("file "+str(output_file)+ " Analysis incomplete!")
            print(str(len(objects))+" of "+str(21)+"..."+"restarting calculation")
            Closure_Str_RN_analysis(str(input_path))
else:
    print("pkl for file "+str(input_path)+ " not created yet...starting calculation")
    Closure_Str_RN_analysis(str(input_path))

    
  