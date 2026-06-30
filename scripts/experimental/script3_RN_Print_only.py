# Script 3: Brute Force Organizations with pyCOT
# ========================================
# 1. LIBRARY LOADING AND CONFIGURATION
# ========================================
import os
import sys

project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
src_path = os.path.join(project_root, 'src')
sys.path.insert(0, src_path)

from pyCOT.io.functions import read_txt
from pyCOT.visualization.rn_visualize import rn_visualize_html, hierarchy_visualize_html
from pyCOT.analysis.Persistent_Modules_Generator import brute_force_organizations

# ========================================
# 2. CREATING THE REACTION_NETWORK OBJECT
# ========================================
file_path = 'projects/XCEPT_Article_1/scripts/data/Complex_example.txt'
# file_path = 'projects/XCEPT_Article_1/scripts/data/CE_2.txt'
file_path = 'networks/testing/AMF_RN_3.txt'
file_path = 'networks/testing/autopoietic_ext.txt'
file_path = 'data/biomodels/BioMD_other/BIOMD0000000237.txt'

file_path = os.path.join(project_root, file_path)
rn = read_txt(file_path)

# ========================================
# 3. SPECIES COLOUR GROUPS
# ========================================
# Add any number of (color, [species, ...]) pairs.
# When a species appears in multiple groups the LAST entry wins.
# Leave the list empty ([]) for the default uncoloured view.
print(f"Species: {[s.name for s in rn.species()]}")
species_color_groups = [
     ("red",    ["p"]),
     ("blue",   ["l","a","b"]),
     ("green",  ["c"]),
]

rn_visualize_html(rn, lst_color_spcs=species_color_groups, filename="rn1.html")

# ========================================
# 3. BRUTE FORCE COMPUTATION
# ========================================
# max_combo_size limits ERC combination depth.  None = exhaustive (slow for
# large networks).  Start with 4-5 to get all small organizations quickly.
bf = brute_force_organizations(rn, max_combo_size=None, verbose=True)

all_closures_frozen   = bf['closures']           # list of frozensets
semi_org_frozen       = bf['semi_organizations']  # list of frozensets (SSM)
org_frozen            = bf['organizations']        # list of frozensets (SM)

# Convert to sets for display / visualization
all_closures_sets = [set(s) for s in all_closures_frozen]
semi_org_sets     = [set(s) for s in semi_org_frozen]
org_sets          = [set(s) for s in org_frozen]

# ========================================
# 4. PRINTED REPORT
# ========================================
print("\n" + "=" * 70)
print("  BRUTE FORCE RESULTS")
print("=" * 70)
print(f"  Unique closures    : {len(all_closures_frozen)}")
print(f"  Semi-organizations : {len(semi_org_frozen)}")
print(f"  Organizations      : {len(org_frozen)}")

print("\n--- Closures ---")
for i, cl in enumerate(all_closures_frozen, 1):
    tag = ""
    if cl in org_frozen:
        tag = "  *** ORG ***"
    elif cl in semi_org_frozen:
        tag = "  (SSM)"
    print(f"  C{i:03d} ({len(cl):2d} sp)  {sorted(cl)}{tag}")

print("\n--- Semi-organizations ---")
for i, so in enumerate(semi_org_frozen, 1):
    tag = "  *** ORG ***" if so in org_frozen else ""
    print(f"  S{i:03d} ({len(so):2d} sp)  {sorted(so)}{tag}")

print("\n--- Organizations ---")
for i, org in enumerate(org_frozen, 1):
    print(f"  O{i:03d} ({len(org):2d} sp)  {sorted(org)}")

print("=" * 70)

# ========================================
# 5. HIERARCHY VISUALIZATION WITH FILTER
# ========================================
# Set relevant_species to a non-empty set to show ONLY nodes that contain
# at least one of those species.  Empty set = show everything.
# relevant_species = set()          # e.g. {"P_H", "H_Res"} or {"I", "P_I", "L"}
relevant_species = {'l'}  # Filter to show only sets containing 'mycellium'
relevant_species_basic = {'a','b'}
relevant_species_parasitic = {'p'}
relevant_species_coop = {'c'}

if True:
    #show_closures  = [s for s in all_closures_sets if relevant_species & s]
    show_semi_orgs = [s for s in semi_org_sets     if relevant_species & s]
    #show_orgs      = [s for s in org_sets          if relevant_species & s]
    show_basic = [s for s in org_sets if relevant_species_basic & s]
    show_parasitic = [s for s in org_sets if relevant_species_parasitic & s]
    show_coop = [s for s in org_sets if relevant_species_coop & s]
    show_all = [s for s in org_sets if len(s)==5]
    # print(f"\nFilter active — showing only sets containing {sorted(relevant_species)}")
    # print(f"  Nodes shown: {len(show_closures)} closures, "
    #       f"{len(show_semi_orgs)} SSMs, {len(show_orgs)} orgs")
else:
    show_closures  = all_closures_sets
    show_semi_orgs = semi_org_sets
    show_orgs      = org_sets
# Color scheme:
#   green  = full organization (SM)
#   yellow = semi-organization only (SSM but not SM)
#   cyan   = closure only (not SSM)
hierarchy_visualize_html(
    show_semi_orgs,
    lst_color_subsets=[
        ("blue", show_basic),
        ("red",  show_parasitic),
        ("green", show_coop),
        ("yellow", show_all)
    ],
    filename="hierarchy_bf.html"
)



