from pyvis.network import Network 

import networkx as nx
import matplotlib.pyplot as plt
import mplcursors
import webbrowser  # Allows opening URLs or local files in the system's default browser
import os  # For handling paths and checking file existence
from collections import defaultdict
import sys
sys.stdout.reconfigure(encoding='utf-8')
import tempfile

###############################################################
### Examples of use of organisations and reaction sub-networks
###############################################################
# Import necessary libraries and modules
import os # Import the os module to interact with the operating system
# os.system('clear') # limpiar la pantalla del terminal en macOS 

import sys # Import the sys module to manipulate the Python runtime environment
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))) # Add the root directory to the PYTHONPATH

# Import the necessary modules from pyCOT  
from pyCOT.io.functions import read_txt    
from pyCOT.rn_visualize import rn_get_visualization, rn_visualize_html
from pyCOT.rn_visualize import hierarchy_get_visualization_html, hierarchy_visualize_html

###################################################################################
# # Load the reaction network from a file
file_path = 'Txt/autopoietic.txt' 
# file_path = 'Txt/2019fig1.txt'
# file_path = 'Txt/2019fig2.txt'
# file_path = 'Txt/non_connected_example.txt' 
# file_path = 'Txt/PassiveUncomforableIndignated_problemsolution.txt' 
file_path = 'Txt/Farm.txt' 
file_path = 'networks/testing/Lake Chad Basin System_0.txt' 

# # Constructs a pyCOT ReactionNetwork object
rn = read_txt(file_path) # DIEGO: Solucionar errores en la lectura del archivo (JOE)

# # String representation of the reaction network
# rn_get_string(rn) 

##############################################################################
# # Example of use of the function rn_visualize_html
##############################################################################
# # Visualize the reaction network
# rn_get_visualization(rn, filename="reaction_network.html")

# # Open the HTML file in the default browser
rn_visualize_html(rn, filename="reaction_network.html")

##############################################################################
# # # Visualize the hierarchy of the reaction network
##############################################################################
from pyCOT.closure_structure import reactive_semi_orgs   # Utiliza file_manipulation.py 
from pyCOT.file_manipulation import load_pyCOT_from_file # Utiliza reaction_network.py

# Creates an object called testRN, which represents the RN with file_manipulation.
# testRN = load_pyCOT_from_file(file_path)
# semi_org = reactive_semi_orgs(testRN) 
# print("-"*70)  
# # for i, semi_org_new in enumerate(semi_org):
# #     print(f"S{i+1} =", semi_org_new) # Print the semi-organisations
# # print("-"*70)

# # hierarchy_get_visualization_html(semi_org,filename="hierarchy_org_autopoietic.html")
# hierarchy_visualize_html(semi_org,filename="hierarchy_org_autopoietic.html")