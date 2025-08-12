# Script 3 Enhanced: Exploring Organizations and Sub-networks with Modern pyCOT

# ========================================
# 1. LIBRARY LOADING AND CONFIGURATION
# ========================================
import os
import sys

# Add the project root directory to the PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import modern pyCOT modules
from pyCOT.io.functions import read_txt, generate_subnetwork_txt
from pyCOT.rn_visualize import rn_visualize_html, hierarchy_visualize_html
from pyCOT.Persistent_Modules import compute_elementary_organizations

# ========================================
# 2. CREATING THE REACTION_NETWORK OBJECT
# ========================================
file_path = 'networks/testing/Farm.txt' 
file_path = 'networks/testing/autopoietic_ext2.txt'
file_path = 'networks/RandomAlife/RN_Ns_20_Norg_6_id_154.txt'  # Input file
#file_path = 'networks/Navarino/RN_IN_05.txt'  # Input file

print(f"Loading reaction network from {file_path}...")
rn = read_txt(file_path)    # Load using ReactionNetwork class

print(f"Reaction Network: {len(rn.species())} species, {len(rn.reactions())} reactions")

# ========================================
# 3. COMPUTE ORGANIZATIONS USING MODERN APPROACH
# ========================================
print("-" * 100)
print("Computing elementary semi-organizations and organizations...")

# Use the modern optimized approach
organizations, elementary_sos, org_statistics, flux_data = compute_elementary_organizations(
    rn, max_generator_size=8, verbose=True
)

print("-" * 100)
print(f"Found {len(elementary_sos)} elementary semi-organizations")
print(f"Found {len(organizations)} organizations")

if not organizations:
    print("No organizations found! Cannot proceed with visualization.")
    sys.exit(1)

# ========================================
# 4. SELECT LARGEST ORGANIZATION
# ========================================
print("\n" + "-" * 70)
print("Selecting largest organization...")
print("-" * 70)

# Find the largest organization(s) by closure size
organization_sizes = [(i, len(org.closure_species)) for i, org in enumerate(organizations)]
max_size = max(size for _, size in organization_sizes)
largest_orgs = [(i, size) for i, size in organization_sizes if size == max_size]

# Select the first of the largest organizations
selected_idx, selected_size = largest_orgs[0]
selected_org = organizations[selected_idx]

print(f"Selected Organization_{selected_idx} with {selected_size} species:")
print(f"Organization species: {selected_org.closure_species}")
print(f"Is P-ERC: {selected_org.is_p_erc}")
print(f"Number of constituent ERCs: {len(selected_org.constituent_ercs)}")

# ========================================
# 5. CREATION OF SUB-NETWORK FOR SELECTED ORGANIZATION
# ========================================
print(f"\nCreating sub-network for Organization_{selected_idx}...")

# Extract organization species as a list of names
org_species_names = [sp.name for sp in selected_org.closure_species]
print(f"Species in organization: {org_species_names}")

# Create sub-reaction network using species names
sub_rn_matrix = rn.sub_reaction_network(org_species_names)  # Changed to use names
triggered_reactions = sub_rn_matrix.stoichiometry_matrix().reactions

print(f"Triggered Reactions: {triggered_reactions}")

# ========================================
# 6. EXPORT AND LOAD SUB-NETWORK
# ========================================
output_file = "sub_network_Farm_enhanced.txt"
# Pass species names instead of Species objects
generate_subnetwork_txt(org_species_names, triggered_reactions, rn, file_name=output_file)

# Load and verify the subnetwork
sub_file_path = f"Txt_sub_network/{output_file}"
rn_sub = read_txt(sub_file_path)

# Verify loaded network
print("-" * 100)
print("Sub-network details:")
print("-" * 100)
print(f"Original species: {org_species_names}")
print(f"Loaded species: {[sp.name for sp in rn_sub.species()]}")
print(f"Reactions: {triggered_reactions}")
print(f"Matrix shape: {rn_sub.stoichiometry_matrix().shape}")

# ========================================
# 7. PREPARE DATA FOR VISUALIZATIONS
# ========================================
print("-" * 100)
print("Preparing data for visualizations...")

# Convert elementary SOs to old format (list of species lists) for hierarchy visualization
elementary_sos_old_format = []
for eso in elementary_sos:
    elementary_sos_old_format.append(list(eso.closure_species))

# Convert organizations to old format
organizations_old_format = []
for org in organizations:
    organizations_old_format.append(list(org.closure_species))

# Prepare filtered organizations (all except the selected one)
filtered_orgs = [list(org.closure_species) for i, org in enumerate(organizations) if i != selected_idx]

# ========================================
# 8. GENERATE HTML VISUALIZATIONS
# ========================================
print("Generating HTML visualizations...")
print("-" * 100)

# 8.1 Visualize full reaction network with the selected organization highlighted
print("Creating full network visualization with selected organization highlighted...")
# Use the already defined org_species_names instead of org_species_list
rn_visualize_html(
    rn, 
    lst_color_spcs=[("green", org_species_names)], 
    filename="reaction_network_enhanced.html"
)
print("✓ Full network visualization saved as: reaction_network_enhanced.html")

# 8.2 Visualize hierarchy of semi-organizations
print("Creating hierarchy visualization of semi-organizations...")
# Convert Species objects to string names before visualization
elementary_sos_old_format = []
for eso in elementary_sos:
    species_names = [sp.name for sp in eso.closure_species]
    elementary_sos_old_format.append(species_names)

# Convert organizations to string names
organizations_old_format = []
for org in organizations:
    species_names = [sp.name for sp in org.closure_species]
    organizations_old_format.append(species_names)

# Convert filtered organizations to string names
filtered_orgs = [
    [sp.name for sp in org.closure_species] 
    for i, org in enumerate(organizations) 
    if i != selected_idx
]

# Now call the visualization with string-based species lists
hierarchy_visualize_html(
    elementary_sos_old_format,
    lst_color_subsets=[
        ("green", organizations_old_format),
        ("yellow", filtered_orgs)
    ],
    filename="hierarchy_semiorg_farm_enhanced.html"
)
print("✓ Hierarchy visualization saved as: hierarchy_semiorg_farm_enhanced.html")

# 8.3 Visualize the subnetwork
print("Creating sub-network visualization...")
# Convert Species objects to names for visualization
sub_species_names = [sp.name for sp in rn_sub.species()]
rn_visualize_html(
    rn_sub, 
    lst_color_spcs=[("green", sub_species_names)],  # Changed parameter name and format
    filename="rn_subnetwork_enhanced.html"
)
print("✓ Sub-network visualization saved as: rn_subnetwork_enhanced.html")

# ========================================
# 9. SUMMARY REPORT
# ========================================
print("\n" + "=" * 80)
print("ANALYSIS COMPLETE - SUMMARY")
print("=" * 80)
print(f"Network: {len(rn.species())} species, {len(rn.reactions())} reactions")
print(f"Elementary Semi-Organizations found: {len(elementary_sos)}")
print(f"Organizations found: {len(organizations)}")
print(f"Selected organization: Organization_{selected_idx}")
print(f"  - Size: {selected_size} species")
print(f"  - Type: {'P-ERC' if selected_org.is_p_erc else 'Multi-ERC'}")
#display all species in the selected organization
org_species_list = selected_org.closure_species   
print(f"  - Species: {org_species_list}")  # Display all
print(f"Sub-network: {len(rn_sub.species())} species, {len(rn_sub.reactions())} reactions")

print(f"\nVisualizations generated:")
print("1. reaction_network_enhanced.html - Full network with selected organization")
print("2. hierarchy_semiorg_farm_enhanced.html - Semi-organization hierarchy")
print("3. rn_subnetwork_enhanced.html - Sub-network of selected organization")

print(f"\nConversion statistics:")
conversion_rate = len(organizations) / len(elementary_sos) * 100 if elementary_sos else 0
print(f"SO → Organization conversion rate: {conversion_rate:.1f}%")

# Display size distribution of organizations
if len(organizations) > 1:
    print(f"\nOrganization size distribution:")
    size_counts = {}
    for org in organizations:
        size = len(org.closure_species)
        size_counts[size] = size_counts.get(size, 0) + 1
    
    for size in sorted(size_counts.keys()):
        count = size_counts[size]
        print(f"  Size {size}: {count} organization{'s' if count > 1 else ''}")

print("\nAnalysis completed successfully!")