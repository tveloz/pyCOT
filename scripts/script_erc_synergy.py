from pyCOT.rn_rustworkx import ReactionNetwork
from pyCOT.ERC_Hierarchy import ERC, ERC_Hierarchy, species_list_to_names
from pyCOT.ERC_Synergy_Complementarity import *
from pyCOT.io.functions import read_txt
import time
from itertools import combinations

# Load network and compute ERCs
file_path = 'networks/testing/ERCs_test2.txt'
file_path = 'networks/Navarino/RN_IN_05.txt'
#file_path = 'networks/testing/Farm.txt'
#file_path = 'networks/biomodels_interesting/bigg_iAF692.txt'
file_path= 'networks/RandomAlife/RN_Ns_40_Norg_20_id_396.txt'
#file_path = 'networks/biomodels_interesting/central_ecoli.txt'
print("Loading network and computing ERCs...")
RN = read_txt(file_path)
hierarchy = ERC_Hierarchy(RN)
hierarchy.build_hierarchy_graph()

print(f"\nNetwork statistics:")
print(f"- Reactions: {len(RN.reactions())}")
print(f"- Species: {len(RN.species())}")
print(f"- ERCs found: {len(hierarchy.ercs)}")

print("\n" + "="*60)
print("BRUTE FORCE vs EFFICIENT but incomplete ALGORITHM COMPARISON")
print("="*60)

# ============================================================================
# BRUTE FORCE APPROACH (pair by pair)
# ============================================================================
print("\n🐌 Running BRUTE FORCE algorithm...")
start_time = time.time()

brute_force_synergies = get_all_fundamental_synergies_brute_force(hierarchy.ercs, hierarchy, RN, verbose=False)

brute_force_time = time.time() - start_time

# ============================================================================
# EFFICIENT APPROACH (all at once)
# ============================================================================
print("🚀 Running EFFICIENT but incomplete  algorithm...")
start_time = time.time()

efficient_synergies = get_fundamental_synergies(hierarchy.ercs, hierarchy, RN, verbose=False)

efficient_time = time.time() - start_time

# ============================================================================
# RESULTS COMPARISON
# ============================================================================
print("\n" + "="*60)
print("COMPARISON RESULTS")
print("="*60)

print(f"\n⏱️  TIMING:")
print(f"  Brute Force: {brute_force_time:.3f} seconds")
print(f"  Efficient:   {efficient_time:.3f} seconds")
if efficient_time > 0:
    speedup = brute_force_time / efficient_time
    print(f"  Speedup:     {speedup:.1f}x faster")
else:
    print(f"  Speedup:     >1000x faster")

print(f"\n📊 SYNERGY COUNTS:")
print(f"  Brute Force: {len(brute_force_synergies)} synergies")
print(f"  Efficient:   {len(efficient_synergies)} synergies")

# ============================================================================
# RESULT VERIFICATION
# ============================================================================
def synergy_to_tuple(syn):
    """Convert synergy to comparable tuple"""
    return (tuple(sorted(syn.rlabel)), syn.plabel)

brute_force_set = {synergy_to_tuple(syn) for syn in brute_force_synergies}
efficient_set = {synergy_to_tuple(syn) for syn in efficient_synergies}

print(f"\n🔍 RESULT VERIFICATION:")
if brute_force_set == efficient_set:
    print("  ✅ IDENTICAL RESULTS - algorithms produce same synergies!")
else:
    print("  ❌ DIFFERENT RESULTS - algorithms disagree!")
    
    only_brute = brute_force_set - efficient_set
    only_efficient = efficient_set - brute_force_set
    
    if only_brute:
        print(f"     Only in brute force ({len(only_brute)}):")
        for syn_tuple in list(only_brute)[:3]:  # Show first 3
            print(f"       {'+'.join(syn_tuple[0])} → {syn_tuple[1]}")
        if len(only_brute) > 3:
            print(f"       ... and {len(only_brute)-3} more")
    
    if only_efficient:
        print(f"     Only in efficient ({len(only_efficient)}):")
        for syn_tuple in list(only_efficient)[:3]:  # Show first 3
            print(f"       {'+'.join(syn_tuple[0])} → {syn_tuple[1]}")
        if len(only_efficient) > 3:
            print(f"       ... and {len(only_efficient)-3} more")
    print("Synergies shared by both algorithms:")
    for syn_tuple in brute_force_set & efficient_set:
        print(f"       {'+'.join(syn_tuple[0])} → {syn_tuple[1]}")
# ============================================================================
# SUMMARY
# ============================================================================
print(f"\n📋 SUMMARY:")
if brute_force_set == efficient_set:
    print(f"  🎯 Efficient algorithm works correctly!")
    print(f"  🚀 Performance improvement: {speedup:.1f}x faster")
    print(f"  💾 Found {len(efficient_synergies)} fundamental synergies")
else:
    print(f"  ⚠️  Algorithms produce different results - need investigation!")
    print(f"  📈 Brute force: {len(brute_force_synergies)} synergies") 
    print(f"  📈 Efficient: {len(efficient_synergies)} synergies")

print(f"\n✅ Comparison complete!")