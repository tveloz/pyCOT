# Script 4: Metabolic Network Visualization with pyCOT

#############################################
# ======================================== source .venv/bin/activate
# 1. LIBRARY LOADING AND CONFIGURATION: 
# ======================================== Depurador de Python: Archivo actual (pyCOT)
# Import necessary libraries of python
import os
import sys 
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))) # Add the root directory to the PYTHONPATH

# Import the necessary modules from pyCOT  
from pyCOT.io.functions import read_txt 
from pyCOT.visualization.rn_visualize import rn_visualize_html, create_bipartite_graph_from_rn, rn_visualize_html_in_out, rn_visualize_png_in_out

# ========================================
# 2. CREATING THE REACTION_NETWORK OBJECT
# ======================================== 
file_path = 'Txt/autopoietic.txt'  # Change this path according to the desired file
# file_path = 'Txt/Farm.txt' 
# file_path = 'Txt/GSM1.txt' 
# file_path = 'Txt/RN_IN_04.txt' 
# file_path = 'Txt/2025fig1.txt' 
# file_path = 'Txt/2025fig2.txt' 
# file_path = 'Txt/2007Dittrich-Speroni_E.coli.txt'

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
# Create the visualization of reaction network using Pyvis (html)
rn_visualize_html(rn,curvature='curvedCCW')

# Create bipartite graph visualization using Graphviz (png)
rn_visualize_png_in_out(graph) 

# Create the metabolic network visualization using Pyvis (html) 
rn_visualize_html_in_out(graph, shape_species_node='circle', curvature=None) #'curvedCCW') #Curva las aristas