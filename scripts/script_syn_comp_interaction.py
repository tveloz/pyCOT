from pyCOT.rn_rustworkx import ReactionNetwork
from pyCOT.ERC_Hierarchy import ERC, ERC_Hierarchy, species_list_to_names
from pyCOT.ERC_Syn_Comp_Interaction import (
    analyze_interaction_patterns, 
    analyze_all_pairwise_interactions,
    analyze_complementarity_distribution,
    get_interaction_statistics,
    print_distribution_analysis
)
from pyCOT.io.functions import read_txt
import time
import pandas as pd

# Load network and compute ERCs
file_path = 'networks/testing/ERCs_test2.txt'
#file_path = 'networks/Navarino/RN_IN_05.txt'
#file_path = 'networks/biomodels_interesting/bigg_iAF692.txt'  # Adjust path as needed

print("Loading network and computing ERCs...")
RN = read_txt(file_path)
hierarchy = ERC_Hierarchy(RN)
hierarchy.build_hierarchy_graph()

print(f"\nNetwork statistics:")
print(f"- Reactions: {len(RN.reactions())}")
print(f"- Species: {len(RN.species())}")
print(f"- ERCs found: {len(hierarchy.ercs)}")

# Test interaction pattern analysis
print("\n" + "="*70)
print("TESTING INTERACTION PATTERN ANALYSIS")
print("="*70)

# First, test individual pair analysis
print("\nTesting individual pair interaction analysis...")

if len(hierarchy.ercs) >= 2:
    erc1, erc2 = hierarchy.ercs[0], hierarchy.ercs[1]
    print(f"\nAnalyzing interaction pattern between {erc1.label} and {erc2.label}:")
    
    start_time = time.time()
    pattern = analyze_interaction_patterns(erc1, erc2, hierarchy, RN)
    end_time = time.time()
    
    print(f"  Analysis time: {end_time - start_time:.4f} seconds")
    print(f"  Pattern: {pattern}")
    
    # Get detailed summary
    summary = pattern.get_interaction_summary()
    print(f"\n📊 Interaction Summary:")
    for key, value in summary.items():
        print(f"  {key}: {value}")
    
    # Show detailed info if there are interactions
    if pattern.has_synergy or pattern.has_complementarity:
        detailed_info = pattern.get_detailed_info()
        
        if detailed_info['synergies']:
            print(f"\n🔄 Synergies found:")
            for syn_info in detailed_info['synergies']:
                print(f"  {syn_info['description']}")
        
        if detailed_info['complementarities']:
            print(f"\n🔗 Complementarities found:")
            for comp_info in detailed_info['complementarities']:
                print(f"  {comp_info['description']}")

# Now test comprehensive analysis of all pairs
print(f"\n" + "="*70)
print("COMPREHENSIVE PAIRWISE ANALYSIS")
print("="*70)

print(f"\nAnalyzing all ERC pairs for interaction patterns...")
start_time = time.time()

# Limit analysis for large networks to avoid excessive computation
max_pairs = 50  # Limit for demonstration
total_possible = len(hierarchy.ercs) * (len(hierarchy.ercs) - 1) // 2

if total_possible <= max_pairs:
    interaction_df = analyze_all_pairwise_interactions(hierarchy, RN, verbose=True)
else:
    print(f"Network has {total_possible} possible pairs. Analyzing first {max_pairs} for demonstration...")
    # Create a smaller hierarchy for demonstration
    small_hierarchy = ERC_Hierarchy(RN)
    small_hierarchy.ercs = hierarchy.ercs[:10]  # Take first 10 ERCs
    small_hierarchy.graph = hierarchy.graph
    small_hierarchy._containment_cache = {}
    small_hierarchy._contained_cache = {}
    small_hierarchy._build_containment_caches()
    
    interaction_df = analyze_all_pairwise_interactions(small_hierarchy, RN, verbose=True)

end_time = time.time()

print(f"\nTotal analysis time: {end_time - start_time:.2f} seconds")
print(f"DataFrame shape: {interaction_df.shape}")

# Display basic DataFrame information
print(f"\n" + "="*70)
print("INTERACTION DATAFRAME OVERVIEW")
print("="*70)

print(f"\nDataFrame columns:")
for col in interaction_df.columns:
    print(f"  {col}")

print(f"\nSample of interaction data:")
print(interaction_df[['erc_pair', 'has_synergy', 'has_complementarity', 'interaction_category']].head(10))

# Get and display basic statistics
print(f"\n" + "="*70)
print("BASIC INTERACTION STATISTICS")
print("="*70)

stats = get_interaction_statistics(interaction_df)

print(f"\n📈 Overall Statistics:")
print(f"  Total pairs analyzed: {stats['total_pairs']}")
print(f"  Active interactions: {stats['active_interactions']}")
print(f"  Synergy only: {stats['synergy_only']}")
print(f"  Complementarity only: {stats['complementarity_only']}")
print(f"  Both synergy & complementarity: {stats['both_synergy_complementarity']}")

if stats['category_counts']:
    print(f"\n🏷️  Interaction Categories:")
    for category, count in sorted(stats['category_counts'].items()):
        print(f"  {category}: {count}")

if any(stats['synergy_type_counts'].values()):
    print(f"\n🔄 Synergy Type Totals:")
    for syn_type, count in stats['synergy_type_counts'].items():
        print(f"  {syn_type}: {count}")

if any(stats['complementarity_type_counts'].values()):
    print(f"\n🔗 Complementarity Type Totals:")
    for comp_type, count in stats['complementarity_type_counts'].items():
        print(f"  {comp_type}: {count}")

# MAIN ANALYSIS: Distribution analysis
print(f"\n" + "="*70)
print("DISTRIBUTION ANALYSIS")
print("="*70)

print(f"\nPerforming complementarity distribution analysis...")
distribution_results = analyze_complementarity_distribution(interaction_df)

# Print formatted distribution analysis
print_distribution_analysis(distribution_results)

# Additional detailed analysis with examples
print(f"\n" + "="*70)
print("DETAILED EXAMPLES AND INSIGHTS")
print("="*70)

# Show examples of different interaction patterns
active_interactions = interaction_df[(interaction_df['has_synergy']) | (interaction_df['has_complementarity'])]

if not active_interactions.empty:
    print(f"\n🔍 EXAMPLE INTERACTION PATTERNS:")
    
    # Examples of synergy + complementarity combinations
    syn_and_comp = active_interactions[(active_interactions['has_synergy']) & (active_interactions['has_complementarity'])]
    if not syn_and_comp.empty:
        print(f"\n  SYNERGY + COMPLEMENTARITY EXAMPLES:")
        for idx, row in syn_and_comp.head(3).iterrows():
            print(f"    {row['erc_pair']}: {row['interaction_category']}")
            print(f"      Synergies: Basic={row['basic_synergy_count']}, Maximal={row['maximal_synergy_count']}, Fundamental={row['fundamental_synergy_count']}")
            print(f"      Complementarities: Type1={row['comp_type1_count']}, Type2={row['comp_type2_count']}, Type3={row['comp_type3_count']}")
    
    # Examples of synergy only
    syn_only = active_interactions[(active_interactions['has_synergy']) & (~active_interactions['has_complementarity'])]
    if not syn_only.empty:
        print(f"\n  SYNERGY ONLY EXAMPLES:")
        for idx, row in syn_only.head(3).iterrows():
            print(f"    {row['erc_pair']}: {row['synergy_signature']}")
    
    # Examples of complementarity only
    comp_only = active_interactions[(~active_interactions['has_synergy']) & (active_interactions['has_complementarity'])]
    if not comp_only.empty:
        print(f"\n  COMPLEMENTARITY ONLY EXAMPLES:")
        for idx, row in comp_only.head(3).iterrows():
            print(f"    {row['erc_pair']}: {row['complementarity_signature']}")

# Key insights based on distribution analysis
print(f"\n💡 KEY INSIGHTS FROM DISTRIBUTION ANALYSIS:")

if 'synergy_complementarity_distribution' in distribution_results:
    syn_dist = distribution_results['synergy_complementarity_distribution']
    
    # Find synergy type with highest complementarity percentage
    max_comp_percentage = 0
    max_comp_type = None
    for syn_type, data in syn_dist.items():
        if data['percentage_with_any_comp'] > max_comp_percentage:
            max_comp_percentage = data['percentage_with_any_comp']
            max_comp_type = syn_type
    
    if max_comp_type:
        print(f"  - {max_comp_type.upper()} synergies are most likely to be complementary ({max_comp_percentage:.1f}%)")
    
    # Check for interesting patterns
    for syn_type, data in syn_dist.items():
        if data['percentage_with_type1'] > 50:
            print(f"  - {syn_type.upper()} synergies frequently involve Type 1 complementarity ({data['percentage_with_type1']:.1f}%)")
        if data['percentage_with_type2'] > 30:
            print(f"  - {syn_type.upper()} synergies often have Type 2 complementarity ({data['percentage_with_type2']:.1f}%)")

if 'non_synergy_complementarity_distribution' in distribution_results:
    non_syn_dist = distribution_results['non_synergy_complementarity_distribution']
    
    if non_syn_dist['percentage_of_non_syn_that_are_comp'] > 20:
        print(f"  - Significant complementarity exists without synergy ({non_syn_dist['percentage_of_non_syn_that_are_comp']:.1f}% of non-synergetic pairs)")
    
    # Identify dominant complementarity type in non-synergetic pairs
    max_type_percentage = max(non_syn_dist['type1_percentage'], non_syn_dist['type2_percentage'], non_syn_dist['type3_percentage'])
    if non_syn_dist['type1_percentage'] == max_type_percentage:
        print(f"  - Type 1 complementarity dominates in non-synergetic pairs ({non_syn_dist['type1_percentage']:.1f}%)")
    elif non_syn_dist['type2_percentage'] == max_type_percentage:
        print(f"  - Type 2 complementarity dominates in non-synergetic pairs ({non_syn_dist['type2_percentage']:.1f}%)")
    elif non_syn_dist['type3_percentage'] == max_type_percentage:
        print(f"  - Type 3 complementarity dominates in non-synergetic pairs ({non_syn_dist['type3_percentage']:.1f}%)")

overall = distribution_results['overall_statistics']
if overall['both_syn_and_comp'] > 0:
    percentage_overlap = (overall['both_syn_and_comp'] / overall['synergetic_pairs']) * 100
    print(f"  - {percentage_overlap:.1f}% of synergetic pairs also show complementarity")

if overall['complementarity_only'] > overall['synergy_only']:
    print(f"  - Network shows more pure complementarity than pure synergy")
elif overall['synergy_only'] > overall['complementarity_only']:
    print(f"  - Network shows more pure synergy than pure complementarity")

# Save results to CSV for further analysis
output_filename = f"interaction_analysis_results.csv"
interaction_df.to_csv(output_filename, index=False)
print(f"\n💾 Results saved to: {output_filename}")

print(f"\n✅ Comprehensive interaction pattern analysis complete!")
print(f"   Use the saved CSV file for detailed statistical analysis in external tools.")