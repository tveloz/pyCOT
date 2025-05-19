# Import libraries
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt 
import re
from scipy.integrate import solve_ivp
from scipy.integrate import odeint
from matplotlib.widgets import Slider
from scipy.optimize import linprog  
import matplotlib.animation as animation 

def generate_subnetwork_txt(species, reactions, stoichiometric_matrix, folder_name="Txt_sub_network", file_name="sub_network.txt"):
    """
    Generates a .txt file with the species, reactions, and stoichiometric matrix corresponding to a subnetwork,
    saved inside a specified folder.
    
    Args:
    species (list): list of species names (e.g., ['s1', 's2', 's6'])
    reactions (list): list of reaction labels (e.g., ['R1', 'R3', 'R7'])
    stoichiometric_matrix (list of list or np.array): corresponding stoichiometric array (species x reactions)
    folder_name (str): name of the folder where the file will be saved
    file_name (str): name of the output file (e.g., 'sub_network.txt')
    
    Returns:
    str: absolute path to the saved file
    """
    # Create folder if it doesn't exist
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
    
    # Full file path
    file_path = os.path.join(folder_name, file_name)
    
    def format_reaction(reactants, products):
        left = '+'.join(reactants) if reactants else ''
        right = '+'.join(products) if products else ''
        return f"{left}=>{right}"
        
    with open(file_path, 'w') as f:
        for idx, reaction in enumerate(reactions):
            column = [row[idx] for row in stoichiometric_matrix]
            reactants = []
            products = []
            for coef, specie in zip(column, species):
                if coef < 0:
                    reactants.append(f"{-coef if abs(coef)!=1 else ''}{specie}")
                elif coef > 0:
                    products.append(f"{coef if coef!=1 else ''}{specie}")
            reaction_str = format_reaction(reactants, products)
            f.write(f"{reaction}:\t{reaction_str};\n")

    # Print absolute file path
    abs_path = os.path.abspath(file_path)
    print(f"\nFile saved at:\n{abs_path}")
    return abs_path

###############################################################
### Examples of use
###############################################################
# Import necessary libraries and modules
import os
import sys

# Add the root directory to the PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import the necessary modules from pyCOT  
from pyCOT.io.functions import read_txt  
from pyCOT.file_manipulation import load_pyCOT_from_file 
from pyCOT.rn_visualize import *

##################################################################
# Generate text file with subnetwork
##################################################################
species_subnetwork = ['s1', 's2']
reactions_subnetwork = ['R2', 'R3', 'R4', 'R5'] # Reactions of the sub-reaction network
S_subnetwork = [[1.00, 1.00, -1.00, 1.00],
                [-1.00, -1.00, 1.00, -1.00]]

# Generate the subnetwork text file ('Txt_sub_network/sub_network.txt')
generate_subnetwork_txt(species_subnetwork, reactions_subnetwork, S_subnetwork, 
                        file_name="sub_network.txt")

# Load the subnetwork from the text file
file_path2 = 'Txt_sub_network/sub_network.txt' 

# Load the subnetwork using rn_rustworkx and file_manipulation
rn2 = read_txt(file_path2) 
testRN2 = load_pyCOT_from_file(file_path2)

print("\n")
print("-"*70,"\nExample of a subnetwork:")
print("-"*70)
species2 = rn2.stoichiometry_matrix().species
print("Sub-network Species =",species2)   # Lista de especies

reactions2 = rn2.stoichiometry_matrix().reactions
print("Sub-network Reactions =",reactions2) # Lista de reacciones

stoichiometry_matrix_subnetwork = rn2.stoichiometry_matrix()
print("Sub-network Stoichiometry Matrix:\n",stoichiometry_matrix_subnetwork) # Imprime la representación de la matriz de estequiometría

rn_visualize_html(testRN2,filename="sub_network.html")