import os       # Imports os to interact with the operating system
import sys      # Imports sys to manipulate interpreter parameters and functions
import warnings # Imports warnings to manage code warnings
warnings.filterwarnings("ignore", category=UserWarning, message="Pick support for FancyArrowPatch is missing.")

# Add the root directory to PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
 
from pyCOT.hierarchy_visualize import *  
from pyCOT.reaction_network import *
from pyCOT.closure_structure import * 
from pyCOT.file_manipulation import * 

#####################################################################################
'''Example manual'''
# input_data1 = [[], ['s1', 's2'], ['s1', 's3'], ['s1', 's2', 's3']]
# input_data1 = [['a'], ['a', 'b'], ['a', 'c', 'd'], ['a', 'b', 'x'], ['a', 'c', 'd', 'y', 'z'], ['a', 'b', 'c', 'd', 'x', 'y', 'z']]
# # Ejecutar hierarchy_visualize 
# hierarchy_visualize(input_data1)
#####################################################################################

##################################################################################### 
'''Example with semi-organisations'''
# Cargar ruta del archivo .txt
# file_path = 'Txt/autopoietic.txt'
# file_path = 'Txt/2019fig1.txt'
# file_path = 'Txt/2019fig2.txt'
# file_path = 'Txt/non_connected_example.txt'
file_path = 'Txt/Farm.txt'   

# Load the file using a function load_pyCOT_from_file()
testRN = load_pyCOT_from_file(file_path)
print(testRN.RnStr)

# Vector of Semi-Organisations calculated at brute force
input_data = reactive_semi_orgs(testRN)   
print("Semi-organisations:", input)
 
# Convert each subset to a set and then to a list to remove duplicates.
unique_subsets = []
for sublist in input_data:
    if set(sublist) not in [set(x) for x in unique_subsets]:
        unique_subsets.append(sublist)
print("Number of Semi-organisations without duplicates=", len(unique_subsets))
   
# Run hierarchy_visualize 
hierarchy_visualize(unique_subsets)