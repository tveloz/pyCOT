# Script 4: Metabolic Network Visualization with pyCOT

#############################################
# ======================================== source .venv/bin/activate
# 1. LIBRARY LOADING AND CONFIGURATION: 
# ======================================== Depurador de Python: Archivo actual (pyCOT)
# Import necessary libraries of python
import os
import sys 
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))) # Add the root directory to the PYTHONPATH
import pandas as pd

# Import the necessary modules from pyCOT  
from pyCOT.io.functions import read_txt, print_reaction_network
from pyCOT.rn_visualize import *

# ========================================
# 2. CREATING THE REACTION_NETWORK OBJECT
# ======================================== 
file_path = 'Txt/autopoietic.txt'  # Change this path according to the desired file
# file_path = 'Txt/Farm.txt' 
# file_path = 'Txt/GSM1.txt' 
# file_path = 'Txt/RN_IN_04.txt' 
# file_path = 'Txt/2025fig1.txt' 
# file_path = 'Txt/2025fig2.txt' 

# Create ReactionNetwork object from text file
rn = read_txt(file_path)   

# =========================================
# 3. CREATING BIPARTITE GRAPH FROM REACTION NETWORK 
# =========================================
graph, metab_nodes, rxn_nodes = create_bipartite_graph_from_rn(rn) 
print(metab_nodes)
print(rxn_nodes)

# =========================================
# 4. VISUALIZATION OF THE METABOLIC NETWORK
# =========================================
rn_visualize_html(rn,curvature='curvedCCW')

# Create bipartite graph visualization using Graphviz (png)
plot_graph_with_graphviz(graph,lst_color_spcs=[("red",["f1", "f2", "f3"]),("green",["y1", "y2", "y3"])],shape_reactions_node='rectangle')

# Create the metabolic network visualization using Pyvis (html)
visualize_metabolic_network(graph, shape_species_node='circle', curvature='curvedCCW') 