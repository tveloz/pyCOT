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
from pyCOT.Persistent_Modules import compute_all_organizations               

# ========================================
# 2. CREATING THE REACTION_NETWORK OBJECT
# ========================================
file_path = 'Txt/Farm.txt'
# file_path = 'Txt/2019influenza1.txt'
# file_path = 'Txt/2019influenza2.txt'
# file_path = 'networks/testing/Farm.txt'  # Input file
# file_path = 'networks/testing/Farm_milk_and_dung.txt'  # Input file
# file_path = 'networks/testing/Farm.txt'
# file_path = 'networks/RandomAlife/RN_Ns_20_Norg_4_id_12.txt'
# file_path = 'networks/RandomAlife/RN_Ns_40_Norg_12_id_358.txt'
# file_path = 'networks/Navarino/RN_IN_05.txt'

rn = read_txt(file_path) 

# ========================================
# 3. COMPUTE SEMI-ORGANIZATIONS AND ORGANIZATIONS
# ======================================== 
elementary_sos, elementary_organizations, all_organizations, all_semi_organizations, statistics, computation_data = compute_all_organizations(
    rn, 
    max_module_size=100,    # Small size for testing
    max_generator_size=100, # Small size for testing
    verbose=True            # Set to True to see detailed computation steps
)

# Set of all semi-organizations 
all_semi_organizations_sets = [set(sp.name for sp in so.combined_closure) for so in all_semi_organizations]
print("all_semi_organizations_sets =", all_semi_organizations_sets)

# Set of all organizations
orgs_sets = [set(sp.name for sp in org.closure_species) for org in elementary_organizations]
orgs_sets = sorted(orgs_sets, key=len) # Sort by size
print("\norgs_sets =", orgs_sets)

# ========================================
# 4. SELECTED ORGANIZATION
# ========================================  
n = 1                        # Index of the organization to analyze
org_n = [orgs_sets[n]]       # Organization with index n
print(f"\norg_{n} =", org_n) 

other_orgs = [s for s in orgs_sets if s != org_n[0]] # Set of other organizations, excluding the selected organization
print("\nother_orgs =", other_orgs)

# ========================================
# 5. CREATE SUBNETWORK FOR THE SELECTED ORGANIZATION
# ======================================== 
output_file = "sub_network_Farm.txt"
reactions_subnet = rn.sub_reaction_network(org_n[0]).stoichiometry_matrix().reactions
print("Triggered Reactions =", reactions_subnet)
generate_subnetwork_txt(org_n[0], reactions_subnet, rn, file_name=output_file)
sub_file_path = f"Txt_sub_network/{output_file}"
rn_sub = read_txt(sub_file_path)

# ========================================
# 6. VISUALIZATION
# ========================================
# Visualize the hierarchy of all semi-organizations, highlighting the selected organization and other organizations
hierarchy_visualize_html(all_semi_organizations_sets,
    lst_color_subsets=[ 
        ("cyan", all_semi_organizations_sets),  # All semi-organizations in yellow
        ("green", org_n),       # Organization selected in green 
        ("yellow", other_orgs)  # Other organizations in yellow
    ])

# Visualize the reaction network in HTML format with the selected organization highlighted in green
rn_visualize_html(rn, lst_color_spcs=[("green", org_n[0] )], filename="rn1_farm.html")

# Visualize the subnetwork in HTML format of the selected organization
rn_visualize_html(rn_sub, global_species_color="green", filename="rn2_farm.html")