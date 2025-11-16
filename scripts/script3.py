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
from pyCOT.Persistent_Modules_Generator import compute_all_organizations

# ========================================
# 2. CREATING THE REACTION_NETWORK OBJECT
# ========================================
file_path = 'Txt/Farm.txt'
# file_path = 'Txt/autopoietic.txt'
# file_path = 'Txt/2007Dittrich-Speroni_fixed_point.txt'
# file_path = 'Txt/2007Dittrich-Speroni_Ex_HIV.txt'
# file_path = 'Txt/2007Dittrich-Speroni_Ex_30.txt'
# file_path = 'Txt/2007Dittrich-Speroni_Ex_302.txt'
# file_path = 'Txt/2019influenza1.txt'
# file_path = 'Txt/2019influenza2.txt'
# file_path = 'Txt/RN_IN_04.txt' # No me corre
# file_path = 'Txt/2007Dittrich-Speroni_E.coli.txt'
# file_path = 'Txt/2010Veloz_Ex_4.txt'

# file_path = 'networks/testing/Farm.txt'  # Input file
# file_path = 'networks/testing/Farm_milk_and_dung.txt'  # Input file
# file_path = 'networks/testing/Farm.txt'
#file_path = 'networks/RandomAlife/RN_Ns_20_Norg_4_id_12.txt'
#file_path = 'networks/RandomAlife/RN_Ns_40_Norg_17_id_564.txt'
# file_path = 'networks/Navarino/RN_IN_05.txt'
# file_path = 'networks/Marine_Ecosystem/Las_Cruces_251021.txt'
# file_path = 'networks/Riverland_model/cause_driven_conflict.txt'



rn = read_txt(file_path) 

# ========================================
# 3. COMPUTE SEMI-ORGANIZATIONS AND ORGANIZATIONS
# ======================================== 
# elementary_sos, elementary_organizations, all_organizations, all_semi_organizations, statistics, computation_data 
results = compute_all_organizations(
    rn, 
    max_generator_size=8, max_organization_size=5,
    # max_module_size=100,    # Small size for testing
    # max_generator_size=100, # Small size for testing
    verbose=True            # Set to True to see detailed computation steps
) 

# === Semi-organizations ===
all_semi_organizations = results['elementary_sos']
all_semi_organizations_sets = [set(sp.name for sp in so.closure_species) for so in all_semi_organizations]
all_semi_organizations_sets = sorted(all_semi_organizations_sets, key=len) # Ordenar por tamaño 
print("\nSemi-organizations:\n", all_semi_organizations_sets) 

# === Organizations ===
all_organizations = results['elementary_organizations']
all_organizations_sets = []
for org in all_organizations:
    if hasattr(org, "combined_closure"):  # Es un Organization
        sset = set(sp.name for sp in org.combined_closure)
    elif hasattr(org, "closure_species"):  # Es un ElementarySO
        sset = set(sp.name for sp in org.closure_species)
    else:
        sset = set()  # fallback por seguridad
    all_organizations_sets.append(sset)
orgs_sets = sorted(all_organizations_sets, key=len)
print("\nOrganizations:\n")
for org in all_organizations_sets:
    print(org)
# ========================================
# 4. SELECTED ORGANIZATION
# ========================================  
# Verificar si hay organizaciones antes de procesar
if not orgs_sets:
    print("\nNo organizations were found in the network.")
    # print("The only organization is the empty set ∅.")
    
    # Visualización básica de la red principal aunque no haya organizaciones
    rn_visualize_html(rn, filename="rn1.html")
    
else: # Si hay organizaciones, proceder con el análisis de la organización seleccionada con el índice n
    n = 1                       # Index of the organization to analyze
    org_n = [orgs_sets[n]]       # Organization with index n
    print(f"\norg_{n} =", org_n) 

    other_orgs = [s for s in orgs_sets if s != org_n[0]] # Set of other organizations, excluding the selected organization
    print("\nother_orgs =", other_orgs)

    # ========================================
    # 5. CREATE SUBNETWORK FOR THE SELECTED ORGANIZATION
    # ======================================== 
    output_file = "sub_network.txt"
    reactions_subnet = rn.sub_reaction_network(org_n[0]).stoichiometry_matrix().reactions
    print("Triggered Reactions =", reactions_subnet)
    generate_subnetwork_txt(org_n[0], reactions_subnet, rn, file_name=output_file)
    sub_file_path = f"Txt_sub_network/{output_file}"
    rn_sub = read_txt(sub_file_path)

    # ========================================
    # 6. VISUALIZATION
    # ========================================
    # Visualize the hierarchy of all semi-organizations, highlighting the selected organization and other organizations
    # hierarchy_visualize_html(all_semi_organizations_sets,
    #     lst_color_subsets=[ 
    #         ("cyan", all_semi_organizations_sets),  # All semi-organizations in yellow
    #         ("green", org_n),       # Organization selected in green 
    #         ("yellow", other_orgs)  # Other organizations in yellow
    #     ])
    hierarchy_visualize_html(all_organizations_sets)
    # Visualize the reaction network in HTML format with the selected organization highlighted in green
    # rn_visualize_html(rn, lst_color_spcs=[("green", org_n[0] )], filename="rn1.html")

    # Visualize the subnetwork in HTML format of the selected organization
    #rn_visualize_html(rn_sub, global_species_color="green", filename="rn2.html")
