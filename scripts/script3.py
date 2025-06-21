# Script 3: Exploring Organizations and Sub-networks with pyCOT

# ========================================
# 1. LIBRARY LOADING AND CONFIGURATION
# ========================================
import os
import sys

# Add the project root directory to the PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import pyCOT modules
from pyCOT.io.functions import read_txt, generate_subnetwork_txt
from pyCOT.rn_visualize import rn_visualize_html, hierarchy_visualize_html
from pyCOT.closure_structure import reactive_semi_orgs, find_organisations, remove_duplicates_list

# Older version of pyCOT: 0.1.0
from pyCOT.file_manipulation import load_pyCOT_from_file # Requires: from pyCOT.reaction_network import * # Requires: import pyCOT.data_conversion as dc

# ========================================
# 2. CREATING THE REACTION_NETWORK OBJECT
# ========================================
file_path = 'networks/testing/Farm.txt'  # Input file
file_path = 'networks/testing/Farm_milk_and_dung.txt'  # Input file
#file_path = 'Txt/Farm.txt'
rn = read_txt(file_path)    # Load using ReactionNetwork class
test_rn = load_pyCOT_from_file(file_path)  # Load using pyCOT object utilities

S= rn.stoichiometry_matrix()  # Stoichiometric matrix
# v1 = [5.3, 1.4, 1.4, 0.5, 1.0, 2.4, 0.5, 1.0, 0.5, 1.0, 1.9, 7.7, 0.5, 0.5, 0.4, 0.5, 0.5] # S.v1 < 0 (Falló)
# # # v1= [14.0, 4.0, 3.0, 1.0, 5.0, 7.0, 1.0, 2.0, 1.0, 5.0, 5.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0] # S.v2>=0

# # # Comprobation: S.v>=0
# xv10 = S @ v1
# print("x_v1 = S.v1 =",xv10)
# ========================================
# 3. COMPUTE SEMI-ORGANISATIONS
# ========================================
print("-" * 100)
semi_orgs = remove_duplicates_list(reactive_semi_orgs(test_rn))  # Remove duplicates

print("-" * 100)
print("Number of semi-organisations =", len(semi_orgs))

# ========================================
# 4. IDENTIFY ORGANISATIONS
# ========================================
organisations, vec_v, vec_x = find_organisations(rn, semi_orgs)

print("\n" + "-" * 70)
print("Organisations and their process and state vectors:")
print("-" * 70)
print("Number of organisations =", len(organisations))

for i, org in enumerate(organisations):
    print(f"Organisation_{i} =", org)

print("\nProcess vectors:")
for i, v in enumerate(vec_v):
    print(f"v_{i} =", v)

print("\nProduction vectors:")
for i, x in enumerate(vec_x):
    print(f"x_{i} =", x)

# ========================================
# 5. CREATION OF A SUB-NETWORK FOR A SELECTED ORGANISATION
# ========================================
n = 1  # Index of the organisation to analyze
org_n = organisations[n]
print(f"\nOrganisation_{n} =", org_n)

reactions_subnet = rn.sub_reaction_network(org_n).stoichiometry_matrix().reactions
print("Triggered Reactions =", reactions_subnet)

# ========================================
# 6. EXPORT SUB-NETWORK TO TEXT FILE
# ========================================
output_file = "sub_network_Farm.txt"
generate_subnetwork_txt(org_n, reactions_subnet, rn, file_name=output_file)

# Load the subnetwork from file
sub_file_path = f"Txt_sub_network/{output_file}"
rn_sub = read_txt(sub_file_path)

print("-" * 100)
print("Example of a subnetwork:")
print("-" * 100)

print("Sub-network Species =", rn_sub.stoichiometry_matrix().species)
print("Sub-network Reactions =", rn_sub.stoichiometry_matrix().reactions)
print("Stoichiometry Matrix:\n", rn_sub.stoichiometry_matrix())

# ========================================
# 7. VISUALIZATION
# ========================================
# 7.1 Visualize full reaction network with the selected organisation highlighted
print("-" * 100)
print("Visualizations:")
print("-" * 100)
rn_visualize_html(rn, lst_color_spcs=[("green", org_n)], filename="reaction_network.html")

# 7.2 Visualize hierarchy of semi-organisations
filtered_orgs = [org for org in organisations if set(org) != set(org_n)]
hierarchy_visualize_html(
    semi_orgs,
    lst_color_subsets=[
        ("green", organisations),
        ("yellow", filtered_orgs)
    ],
    filename="hierarchy_semiorg_farm.html"
)

# 7.3 Visualize the subnetwork
rn_visualize_html(rn_sub, global_species_color="green", filename="rn2_farm.html")