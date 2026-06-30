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

file_path = os.path.join(project_root, file_path)
rn = read_txt(file_path)
rn_visualize_html(rn, filename="rn1.html")

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
relevant_species = {'mycellium'}  # Filter to show only sets containing 'mycellium'

if relevant_species:
    show_closures  = [s for s in all_closures_sets if relevant_species & s]
    show_semi_orgs = [s for s in semi_org_sets     if relevant_species & s]
    show_orgs      = [s for s in org_sets          if relevant_species & s]
    print(f"\nFilter active — showing only sets containing {sorted(relevant_species)}")
    print(f"  Nodes shown: {len(show_closures)} closures, "
          f"{len(show_semi_orgs)} SSMs, {len(show_orgs)} orgs")
else:
    show_closures  = all_closures_sets
    show_semi_orgs = semi_org_sets
    show_orgs      = org_sets

# Color scheme:
#   green  = full organization (SM)
#   yellow = semi-organization only (SSM but not SM)
#   cyan   = closure only (not SSM)
hierarchy_visualize_html(
    show_closures,
    lst_color_subsets=[
        ("yellow", show_semi_orgs),
        ("green",  show_orgs),
    ],
    filename="hierarchy_bf.html"
)

# ========================================
# 6. RN VISUALIZATION WITH SPECIES COLORING
# ========================================
# rn_visualize_html accepts lst_color_spcs as a list of (color, [species, ...])
# tuples.  When a species belongs to multiple groups the LAST matching entry
# wins, so list higher-priority groups last.
#
# Three options are provided; uncomment the one you want.

# ---- helper: union of a list of sets ----
def _species_union(list_of_sets):
    result = set()
    for s in list_of_sets:
        result |= s
    return result

# ---- Option A: auto-derive groups from brute-force results ---------------
# orange = in at least one closure
# yellow = in at least one semi-organisation (overrides orange)
# green  = in at least one full organisation  (overrides yellow)
all_closure_species = sorted(_species_union(all_closures_sets))
all_ssm_species     = sorted(_species_union(semi_org_sets))
all_org_species     = sorted(_species_union(org_sets))

rn_visualize_html(
    rn,
    lst_color_spcs=[
        ("orange", all_closure_species),
        ("yellow", all_ssm_species),
        ("green",  all_org_species),
    ],
    filename="rn_colored_auto.html",
)

# ---- Option B: manual species groups -------------------------------------
# Define any number of named colour groups.  Edit species names and colours
# to match your network.
#
# custom_groups = [
#     ("red",    ["species_A", "species_B"]),
#     ("blue",   ["species_C", "species_D"]),
#     ("purple", ["species_E"]),
# ]
# rn_visualize_html(rn, lst_color_spcs=custom_groups,
#                   filename="rn_colored_custom.html")

# ---- Option C: highlight species of one specific organisation ------------
# Pick by index into org_sets (0-based).
#
# if org_sets:
#     idx = 0
#     rn_visualize_html(
#         rn,
#         lst_color_spcs=[("green", sorted(org_sets[idx]))],
#         filename=f"rn_org{idx}.html",
#     )
