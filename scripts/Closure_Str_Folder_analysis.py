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
from Closure_Str_RN_Computing import Closure_Str_RN_Computing 


folder_path = 'networks/biomodels_all_txt'
txt_files = glob.glob(os.path.join(folder_path, '*.txt'))

txt_files_ok=[ff.replace("\\", "/") for ff in txt_files]
RN_size=[]
for file_path in txt_files_ok:
    input_path = Path(file_path)
    #input_path = Path('networks/biomodels_interesting/BIOMD0000000652_manyOrgs.txt')
    RN = load_pyCOT_from_file(input_path)
    number = re.findall(r'\d+', file_path)
    RN_size.append([number[0],len(RN.SpStr),len(RN.RnStr)])
cutted_list = [el for el in RN_size if el[2]>=70 and el[2]<1000] 
#sorted_list = sorted(cutted_list, key=lambda x: x[2])



sorted_list = sorted(cutted_list, key=lambda x: x[2])
print(sorted_list)
print([len(cutted_list),(len(cutted_list)/len(RN_size))])
plt.plot([el[2] for el in sorted_list])


# Add titles and labels
plt.title('RN Sorted by Size')
plt.xlabel('Size')
plt.ylabel('Network')

# Show the plot
plt.show()