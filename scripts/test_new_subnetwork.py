###############################################################
### Examples of use of organisations and reaction sub-networks
###############################################################
# Import necessary libraries 
import os # Import the os module to interact with the operating system
os.system('clear') # limpiar la pantalla del terminal en macOS 

import sys # Import the sys module to manipulate the Python runtime environment
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))) # Add the root directory to the PYTHONPATH

# Import the necessary modules from pyCOT  
from pyCOT.io.functions import read_txt  
from pyCOT.io.functions import generate_subnetwork_txt
from pyCOT.file_manipulation import load_pyCOT_from_file
from pyCOT.closure_structure import reactive_semi_orgs
from pyCOT.closure_structure import find_organisations
from pyCOT.rn_visualize import *                         
from pyCOT.simulations import simulation
from pyCOT.plot_dynamics import plot_series_ode

###################################################################################
## 1. Load the reaction network from a file
# file_path = 'Txt/autopoietic.txt' 
# file_path = 'Txt/2019fig1.txt'
# file_path = 'Txt/2019fig2.txt'
# file_path = 'Txt/non_connected_example.txt' 
# file_path = 'Txt/PassiveUncomforableIndignated_problemsolution.txt' 
file_path = 'Txt/Farm.txt' 

###################################################################################
## 2. Constructs a pyCOT ReactionNetwork object with rn_rustworkx
rn = read_txt(file_path) # DIEGO: Solucionar errores en la lectura del archivo (JOE)

# Creates an object called testRN, which represents the RN with file_manipulation.
testRN = load_pyCOT_from_file(file_path)

###################################################################################
## 3. Calculate the semi-organisations 
print("-"*70) 
semi_org = reactive_semi_orgs(testRN)   

for i, semi_org_new in enumerate(semi_org):
    print(f"S{i+1} =", semi_org_new) # Print the semi-organisations
print("-"*70)

###################################################################################
## 4. Find organisations based on the semi-organisations and check if they are self-maintaining 
organisations, vector_process, vector_x = find_organisations(rn, semi_org)

print("\n")
print("-"*70,"\nOrganisations and their process and state vectors:")
print("-"*70)
print("Number of organisations =", len(organisations)) 
print("Organisations:")
for i, org in enumerate(organisations):
    print(f"Organisation_{i} =", org) # Print the organisations 
# CL={grass, 
# cows, 
# milk, 
# dung, 
# fertilizer, 
# worms, 
# grain, 
# straw, 
# water, 
# infrastructure, 
# farmer, 
# money}.

print("\nProcess vector:")
for i, vec in enumerate(vector_process):
    print(f"v_{i} =", vec) # Print the vector process of each organisation

print("\nProduction vector:")
for i, vec in enumerate(vector_x):
    print(f"x_{i} =", vec) # Print the vector production of each organisation

###################################################################################
# 5. Test: x_v=S.v
print("\n")
print("-"*70,"\nTEST: x_v = S.v ≥ 0")
print("-"*70) 

n = 1 # Index of the organisation to be analyzed 
organisation = organisations[n] # Organisation to be analyzed

print(f"Organisation_{n} = {organisation}\n") # Organisation to be analyzed

reactions_subnetwork = rn.sub_reaction_network(organisation).stoichiometry_matrix().reactions # Reactions of the sub-reaction network
print("Triggered Reactions =",reactions_subnetwork) # Reactions of the sub-reaction network

S_org_n = rn.sub_reaction_network(organisation).stoichiometry_matrix() # Sub-reaction network of the organisation 
print(f"\nSub Stoichiometric Matrix of Organisation_{n}\n", S_org_n) # Matriz de estequiometría de la organización  

print(f"\nVector process of organisation_{n}:\n", vector_process[n]) # Vector process of the organisation
print(f"Vector x of organisation_{n}:\n", vector_x[n]) # Vector x of the organisation
print("x_v = S.v =", S_org_n@vector_process[n]) # Multiplication of the stoichiometric matrix by the vector process
print("-"*70)

# ##################################################################
# # Example 1: Generate a text file of a subnetwork with 
# # a list of species and reactions predefined for file_path = 'Txt/autopoietic.txt'
# ##################################################################
# species_subnetwork = ['s1', 's2'] 
# reactions_subnetwork = ['R2', 'R3', 'R4', 'R5']  

# file_name = "sub_network_autopoietic.txt"
# generate_subnetwork_txt(species_subnetwork, reactions_subnetwork, rn, 
#                         file_name=file_name) 

# # Load the new reaction sub_network 
# file_path2 = f"Txt_sub_network/{file_name}" 

# rn2 = read_txt(file_path2)  

# print("\n")
# print("-"*70,"\nExample of a subnetwork:")
# print("-"*70)
# species2 = rn2.stoichiometry_matrix().species
# print("Sub-network Species =",species2)   # Lista de especies

# reactions2 = rn2.stoichiometry_matrix().reactions
# print("Sub-network Reactions =",reactions2) # Lista de reacciones

# S_subnetwork = rn2.stoichiometry_matrix()
# print("Sub-network Stoichiometry Matrix:\n",S_subnetwork) # Imprime la representación de la matriz de estequiometría
# print("\n")

# print("-"*70)
# print("Visualizing the reaction network and the subnetwork")
# print("-"*70)
# # Visualize the hierarchy of the reaction network
# hierarchy_visualize_html(organisations,filename="hierarchy_org_autopoietic.html")

# # Visualize the reaction network
# rn_visualize_html(rn,filename="rn1_autopoietic.html")

# # Visualize the subnetwork
# rn_visualize_html(rn2,filename="rn2_autopoietic.html")

##################################################################
# Example 2: Generate text file with Sub-network of the Organisation for file_path = 'Txt/Farm.txt' 
##################################################################
species_subnetwork = organisation
reactions_subnetwork = reactions_subnetwork

file_name = "sub_network_Farm.txt"
generate_subnetwork_txt(species_subnetwork, reactions_subnetwork, rn, 
                        file_name=file_name) 

# Load the new reaction sub_network 
file_path2 = f"Txt_sub_network/{file_name}" 

rn2 = read_txt(file_path2)  

print("\n")
print("-"*70,"\nExample of a subnetwork:")
print("-"*70)
species2 = rn2.stoichiometry_matrix().species
print("Sub-network Species =",species2)   # Lista de especies

reactions2 = rn2.stoichiometry_matrix().reactions
print("Sub-network Reactions =",reactions2) # Lista de reacciones

S_subnetwork = rn2.stoichiometry_matrix()
print("Sub-network Stoichiometry Matrix:\n",S_subnetwork) # Imprime la representación de la matriz de estequiometría
print("\n")

print("-"*70)
print("Visualizing the reaction network and the subnetwork")
print("-"*70)
# Visualize the hierarchy of the reaction network
hierarchy_visualize_html(organisations,filename="hierarchy_org_farm.html")

# Visualize the reaction network
rn_visualize_html(rn,filename="rn1_farm.html")

# Visualize the subnetwork
rn_visualize_html(rn2,filename="rn2_farm.html")

##################################################################
# Example 3: Simulation of the subnetwork
##################################################################
# # Plot the time series of the ODE of the subnetwork
# time_series, flux_vector = simulation(rn2, rate="mak") 
# print("Time series of ODE of the subnetwork:")
# print(time_series)
# plot_series_ode(time_series)
# print("Flux series of ODE of the subnetwork:")
# print(flux_vector)
# plot_series_ode(flux_vector)