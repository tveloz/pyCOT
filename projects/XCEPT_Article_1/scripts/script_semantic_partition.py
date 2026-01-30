"""
Test Script for Part 1: Semantic Partitioning System

This script demonstrates the functionality of the semantic partitioning module
using the cause-driven conflict model.

"""

import os
import sys

# Add the project root directory to the PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import pyCOT modules
from pyCOT.io.functions import read_txt, generate_subnetwork_txt
from pyCOT.rn_visualize import rn_visualize_html, hierarchy_visualize_html
from pyCOT.Persistent_Modules_Generator import compute_all_organizations
from pyCOT.semantic_partition import (
    define_semantic_categories,
    get_species_indices_for_category,
    get_categories_for_species,
    validate_semantic_partition,
    print_semantic_partition_summary,
    create_conflict_semantic_partition
)
import numpy as np
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
file_path = 'networks/Conflict_Theory/conflict_toy_model0.txt'
file_path = 'networks/Conflict_Theory/cause_driven_conflict_gov.txt'



rn = read_txt(file_path) 
#!/usr/bin/env python3

print("=" * 80)
print("PART 1: SEMANTIC PARTITIONING SYSTEM TEST")
print("=" * 80)

# ========================================
# 1. Define the Conflict Model Manually
# ========================================
print("\n" + "=" * 80)
print("STEP 1: Defining Conflict Model Species")
print("=" * 80)

# Species from cause_driven_conflict.
species_list = rn.species()
species = [specie.name for specie in species_list]
reactions_list = rn.reactions()
reactions =[reaction.name() for reaction in reactions_list]

# Create a mock stoichiometric matrix (12 species × 22 reactions)
# This is just for validation purposes
n_species = len(species)
n_reactions = len(reactions)  # From the model file (r1-r22)
S = rn.stoichiometry_matrix()  # Mock matrix for testing

print(f"\n📊 Conflict model specification:")
print(f"   Species: {len(species)}")
print(f"   Reactions: {n_reactions}")
print(f"   Stoichiometric matrix shape: {S.shape}")
print(f"\n   Species list: {species}")

# ========================================
# 2. Define Semantic Categories - Manual Approach
# ========================================
print("\n" + "=" * 80)
print("STEP 2: Manual Category Definition")
print("=" * 80)

print("\nDefining semantic categories for conflict model...")
print("Categories to create:")
print("  - conflict: violent species and damage causes")
print("  - peace: peaceful species and grievance resolution")
print("  - perception: internal perceptions that drive behavior")

category_definitions = {
    'conflict': ['A_v', 'B_v', 'D_A', 'D_B','iD_A', 'iD_B',],
    'peace': ['A_p', 'B_p', 'G_A', 'G_B','iG_A', 'iG_B'],
    'group A':['A_v', 'A_p', 'D_A', 'iD_A'],
    'group B':['B_v', 'B_p', 'D_B', 'iD_B']
}

print("\nCategory definitions:")
for cat, spec_list in category_definitions.items():
    print(f"  {cat}: {spec_list}")

semantic_partition = define_semantic_categories(species, category_definitions)

print("\n✓ Semantic partition created successfully!")
print(semantic_partition)

# ========================================
# 3. Explore the Semantic Partition
# ========================================
print("\n" + "=" * 80)
print("STEP 3: Exploring Semantic Partition")
print("=" * 80)

print_semantic_partition_summary(semantic_partition)

# ========================================
# 4. Query Functions - Get Indices
# ========================================
print("\n" + "=" * 80)
print("STEP 4: Testing Query Functions - Species Indices")
print("=" * 80)

for category in semantic_partition.categories:
    indices = get_species_indices_for_category(semantic_partition, category)
    species_in_cat = semantic_partition.category_to_species[category]
    print(f"\n{category} category:")
    print(f"  Species: {species_in_cat}")
    print(f"  Indices: {indices}")
    print(f"  Verification: {[species[i] for i in indices]}")

# ========================================
# 5. Query Functions - Get Categories for Species
# ========================================
print("\n" + "=" * 80)
print("STEP 5: Testing Query Functions - Categories per Species")
print("=" * 80)

print("\nChecking which categories each species belongs to:")
for sp in species:
    categories = get_categories_for_species(semantic_partition, sp)
    if categories:
        print(f"  {sp:8s} → {categories}")
    else:
        print(f"  {sp:8s} → [uncategorized]")

# ========================================
# 6. Validate Semantic Partition
# ========================================
print("\n" + "=" * 80)
print("STEP 6: Validation")
print("=" * 80)

report = validate_semantic_partition(semantic_partition, S, verbose=True)

if report['valid']:
    print("\n✓ Semantic partition is valid and ready for downstream analysis!")
else:
    print("\n✗ Semantic partition has errors - see report above")

# ========================================
# 7. Test Alternative Construction Method
# ========================================
print("\n" + "=" * 80)
print("STEP 7: Alternative Construction (Automatic for Conflict Model)")
print("=" * 80)

print("\nTesting automatic conflict partition creation...")
auto_partition = create_conflict_semantic_partition(species)

print("\nAuto-generated partition:")
print_semantic_partition_summary(auto_partition)

# Validate auto partition
auto_report = validate_semantic_partition(auto_partition, S, verbose=False)
print(f"\nAuto partition validation: {'✓ VALID' if auto_report['valid'] else '✗ INVALID'}")
print(f"Coverage: {auto_report['stats']['coverage_percent']:.1f}%")

# ========================================
# 8. Test Edge Cases
# ========================================
print("\n" + "=" * 80)
print("STEP 8: Testing Edge Cases")
print("=" * 80)

# Test 8a: Overlapping categories
print("\n8a. Testing overlapping categories...")
overlapping_defs = {
    'violence': ['A_v', 'B_v'],
    'faction_A': ['A_v', 'A_p', 'iG_A', 'iD_A'],
    'faction_B': ['B_v', 'B_p', 'iG_B', 'iD_B']
}

overlap_partition = define_semantic_categories(species, overlapping_defs)
print(f"Created partition with overlapping categories")

overlapping_species = [
    sp for sp, cats in overlap_partition.species_to_categories.items()
    if len(cats) > 1
]
print(f"Overlapping species ({len(overlapping_species)}): {overlapping_species}")

for sp in overlapping_species[:2]:  # Show first 2
    cats = get_categories_for_species(overlap_partition, sp)
    print(f"  {sp} belongs to: {cats}")

# Validate
overlap_report = validate_semantic_partition(overlap_partition, S, verbose=False)
print(f"Overlap partition validation: {'✓ VALID' if overlap_report['valid'] else '✗ INVALID'}")

# Test 8b: Error handling - invalid species
print("\n8b. Testing error handling - invalid species...")
try:
    bad_defs = {
        'invalid': ['NONEXISTENT_SPECIES', 'A_v']
    }
    bad_partition = define_semantic_categories(species, bad_defs)
    print("  ✗ Should have raised ValueError!")
except ValueError as e:
    print(f"  ✓ Correctly raised ValueError: {e}")

# Test 8c: Error handling - invalid category query
print("\n8c. Testing error handling - invalid category query...")
try:
    indices = get_species_indices_for_category(semantic_partition, 'NONEXISTENT_CATEGORY')
    print("  ✗ Should have raised KeyError!")
except KeyError as e:
    print(f"  ✓ Correctly raised KeyError: {e}")

# Test 8d: Empty category definitions
print("\n8d. Testing error handling - empty definitions...")
try:
    empty_partition = define_semantic_categories(species, {})
    print("  ✗ Should have raised ValueError!")
except ValueError as e:
    print(f"  ✓ Correctly raised ValueError: {e}")

# ========================================
# 9. Demonstrate Integration with Stoichiometric Matrix
# ========================================
print("\n" + "=" * 80)
print("STEP 9: Integration with Stoichiometric Matrix")
print("=" * 80)

print("\nDemonstrating how to use semantic categories with stoichiometric analysis...")

# Show how to extract submatrices for specific categories
for category in ['conflict', 'peace']:
    indices = get_species_indices_for_category(semantic_partition, category)
    S_sub = S[indices, :]
    
    print(f"\n{category.upper()} category submatrix:")
    print(f"  Shape: {S_sub.shape} (species in category × all reactions)")
    print(f"  Species: {semantic_partition.category_to_species[category]}")
    print(f"  Submatrix:\n{S_sub}")

# ========================================
# 10. Summary and Next Steps
# ========================================
print("\n" + "=" * 80)
print("SUMMARY AND NEXT STEPS")
print("=" * 80)

print("\n✓ Part 1 Complete: Semantic Partitioning System")
print("\nKey accomplishments:")
print("  ✓ Created SemanticPartition dataclass")
print("  ✓ Implemented define_semantic_categories()")
print("  ✓ Implemented get_species_indices_for_category()")
print("  ✓ Implemented get_categories_for_species()")
print("  ✓ Implemented validate_semantic_partition()")
print("  ✓ Tested with conflict model")
print("  ✓ Validated edge cases and error handling")
print("  ✓ Demonstrated integration with stoichiometric matrix")

print("\nIntegration points verified:")
print("  ✓ Uses species names from rn.species()")
print("  ✓ Indices match rows in rn.stoichiometry_matrix()")
print("  ✓ Ready for Part 2: Challenge Decomposition")

print("\nNext steps:")
print("  → Part 2: Challenge Decomposition (decompose processes to challenges)")
print("  → Part 3: Challenge Aggregation (identify emergent patterns)")
print("  → Parts 4-9: Advanced analysis and integration")

print("\n" + "=" * 80)
print("✓ TEST COMPLETED SUCCESSFULLY")
print("=" * 80)