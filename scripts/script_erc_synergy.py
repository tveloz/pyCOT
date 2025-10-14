#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script for optimized fundamental synergy detection with detailed diagnostics
Run this after adding the optimization functions to ERC_Synergy_Complementarity.py

This version includes comprehensive analysis of differences between algorithms,
showing species sets and minimal generators for complete diagnosis.
"""

import time
from collections import defaultdict
import networkx as nx

# Import the necessary modules
from pyCOT.ERC_Hierarchy import (
    ERC,
    ERC_Hierarchy,
    closure,
    generators
)
from pyCOT.ERC_Synergy_Complementarity import (
    get_all_fundamental_synergies_brute_force,
    get_all_fundamental_synergies_optimized
)
from pyCOT.io.functions import read_txt


# ============================================================================
# HELPER FUNCTIONS FOR VERIFICATION
# ============================================================================

def target_is_contained_by(small_target, large_target, hierarchy):
    """
    Check if small_target is strictly contained by large_target.
    Returns True if small_target ⊂ large_target (strict containment).
    """
    if small_target == large_target:
        return False  # Not STRICT containment
    
    if not hierarchy.graph:
        return False
    
    # small_target is contained by large_target if there's a path from small to large
    return nx.has_path(hierarchy.graph, small_target, large_target)


def is_more_fundamental_than(pair1_labels, pair2_labels, hierarchy):
    """
    Check if pair1 is more fundamental than pair2.
    
    Pair1 is more fundamental if:
    - Both components of pair1 are contained by (or equal to) components of pair2
    - At least one is strictly contained
    
    Parameters
    ----------
    pair1_labels : tuple
        (label1, label2) for first pair
    pair2_labels : tuple
        (label1, label2) for second pair
    hierarchy : ERC_Hierarchy
        The hierarchy object
        
    Returns
    -------
    bool
        True if pair1 is more fundamental than pair2
    """
    if pair1_labels == pair2_labels:
        return False
    
    if not hierarchy.graph:
        return False
    
    # Unpack pairs
    p1_a, p1_b = pair1_labels
    p2_a, p2_b = pair2_labels
    
    # Check both possible orderings (since pairs are unordered)
    # Option 1: p1_a corresponds to p2_a, p1_b corresponds to p2_b
    option1_a_contained = (p1_a == p2_a or nx.has_path(hierarchy.graph, p1_a, p2_a))
    option1_b_contained = (p1_b == p2_b or nx.has_path(hierarchy.graph, p1_b, p2_b))
    option1_strict = (nx.has_path(hierarchy.graph, p1_a, p2_a) or 
                     nx.has_path(hierarchy.graph, p1_b, p2_b))
    
    if option1_a_contained and option1_b_contained and option1_strict:
        return True
    
    # Option 2: p1_a corresponds to p2_b, p1_b corresponds to p2_a
    option2_a_contained = (p1_a == p2_b or nx.has_path(hierarchy.graph, p1_a, p2_b))
    option2_b_contained = (p1_b == p2_a or nx.has_path(hierarchy.graph, p1_b, p2_a))
    option2_strict = (nx.has_path(hierarchy.graph, p1_a, p2_b) or 
                     nx.has_path(hierarchy.graph, p1_b, p2_a))
    
    if option2_a_contained and option2_b_contained and option2_strict:
        return True
    
    return False


def print_erc_details(label, erc, RN, indent="       "):
    """Print detailed information about an ERC."""
    species = sorted(erc.get_closure_names(RN))
    print(f"{indent}ERC {label}: {{{', '.join(species)}}}")
    print(f"{indent}Minimal generators:")
    for gen_idx, gen in enumerate(erc.min_generators, 1):
        gen_species_names = sorted([sp.name for sp in gen])
        print(f"{indent}  Gen {gen_idx}: {{{', '.join(gen_species_names)}}}")


def analyze_synergy_mismatch(pair, target, hierarchy, RN, context=""):
    """Analyze why a synergy might be found by one algorithm but not the other."""
    print(f"\n   {context}")
    print(f"   Synergy: {pair[0]} + {pair[1]} → {target}")
    print("   " + "-"*65)
    
    # Get the ERC objects
    erc1 = hierarchy.get_erc_by_label(pair[0])
    erc2 = hierarchy.get_erc_by_label(pair[1])
    target_erc = hierarchy.get_erc_by_label(target)
    
    if not (erc1 and erc2 and target_erc):
        print("   ERROR: Could not retrieve ERC objects")
        return
    
    # Show base ERCs
    print(f"\n   BASE ERCs:")
    erc1_species = sorted(erc1.get_closure_names(RN))
    erc2_species = sorted(erc2.get_closure_names(RN))
    print(f"     ERC {pair[0]}: {{{', '.join(erc1_species)}}}")
    print(f"       Minimal generators:")
    for gen_idx, gen in enumerate(erc1.min_generators, 1):
        gen_names = sorted([sp.name for sp in gen])
        print(f"         Gen {gen_idx}: {{{', '.join(gen_names)}}}")
    
    print(f"\n     ERC {pair[1]}: {{{', '.join(erc2_species)}}}")
    print(f"       Minimal generators:")
    for gen_idx, gen in enumerate(erc2.min_generators, 1):
        gen_names = sorted([sp.name for sp in gen])
        print(f"         Gen {gen_idx}: {{{', '.join(gen_names)}}}")
    
    # Show union
    union_species = sorted(set(erc1_species) | set(erc2_species))
    print(f"\n   UNION of base ERCs: {{{', '.join(union_species)}}}")
    
    # Show target ERC
    print(f"\n   TARGET ERC:")
    target_species = sorted(target_erc.get_closure_names(RN))
    print(f"     ERC {target}: {{{', '.join(target_species)}}}")
    print(f"       Minimal generators:")
    union_set = set(union_species)
    for gen_idx, gen in enumerate(target_erc.min_generators, 1):
        gen_names = set([sp.name for sp in gen])
        gen_names_sorted = sorted(gen_names)
        covered = gen_names.issubset(union_set)
        coverage_icon = "✓" if covered else "✗"
        print(f"         Gen {gen_idx}: {{{', '.join(gen_names_sorted)}}} [{coverage_icon}]")
    
    # Check reactions from union that are novel
    print(f"\n   REACTIONS:")
    erc1_reacs = erc1.get_reacs(RN)
    erc2_reacs = erc2.get_reacs(RN)
    joint_closure = closure(RN, list(set(erc1.get_closure(RN)) | set(erc2.get_closure(RN))))
    joint_reacs = RN.get_reactions_from_species(joint_closure)
    
    erc1_reac_names = set(r.name() for r in erc1_reacs)
    erc2_reac_names = set(r.name() for r in erc2_reacs)
    joint_reac_names = set(r.name() for r in joint_reacs)
    
    novel_reacs = joint_reac_names - (erc1_reac_names | erc2_reac_names)
    print(f"     Base 1 reactions: {len(erc1_reac_names)}")
    print(f"     Base 2 reactions: {len(erc2_reac_names)}")
    print(f"     Joint closure reactions: {len(joint_reac_names)}")
    print(f"     Novel reactions: {len(novel_reacs)}")
    if novel_reacs and len(novel_reacs) <= 5:
        print(f"       Novel: {', '.join(sorted(list(novel_reacs)))}")


def analyze_synergy_structure(ercs, hierarchy, RN, verbose=True):
    """
    Analyze the structure of maximal synergies to understand pruning potential.
    This helps us understand how much the optimization can help.
    """
    from pyCOT.ERC_Synergy_Complementarity import get_maximal_synergies
    from itertools import combinations
    
    if verbose:
        print("\n" + "="*60)
        print("SYNERGY STRUCTURE ANALYSIS")
        print("="*60)
    
    # Compute hierarchy levels
    if hierarchy.graph:
        levels = ERC.get_node_levels(hierarchy.graph)
    else:
        levels = {erc.label: 0 for erc in ercs}
    
    # Group ERCs by level
    synergies_by_level = defaultdict(list)
    target_to_producers = defaultdict(list)  # target -> [(pair, level)]
    pair_to_targets = defaultdict(list)      # pair -> [targets]
    
    total_maximal = 0
    pairs_with_synergies = 0
    
    for erc1, erc2 in combinations(ercs, 2):
        maximal_syns = get_maximal_synergies(erc1, erc2, hierarchy, RN)
        
        if maximal_syns:
            pairs_with_synergies += 1
            total_maximal += len(maximal_syns)
            
            pair_key = tuple(sorted([erc1.label, erc2.label]))
            level_sum = levels[erc1.label] + levels[erc2.label]
            
            for syn in maximal_syns:
                target = syn.plabel
                synergies_by_level[level_sum].append(syn)
                target_to_producers[target].append((pair_key, level_sum))
                pair_to_targets[pair_key].append(target)
    
    if verbose:
        print(f"\nMaximal synergy statistics:")
        print(f"  Total ERC pairs: {len(list(combinations(ercs, 2)))}")
        print(f"  Pairs with synergies: {pairs_with_synergies}")
        print(f"  Total maximal synergies: {total_maximal}")
        
        if pairs_with_synergies > 0:
            print(f"  Avg synergies per synergistic pair: {total_maximal/pairs_with_synergies:.2f}")
    
    # Analyze pair domination potential
    expected_fundamental_by_pair = 0
    for target, producers in target_to_producers.items():
        if len(producers) == 1:
            expected_fundamental_by_pair += 1
        else:
            expected_fundamental_by_pair += 1
    
    # Analyze target containment potential
    pairs_with_contained_targets = 0
    for pair_key, targets in pair_to_targets.items():
        if len(targets) > 1:
            has_containment = False
            for i, t1 in enumerate(targets):
                for t2 in targets[i+1:]:
                    if target_is_contained_by(t1, t2, hierarchy) or target_is_contained_by(t2, t1, hierarchy):
                        has_containment = True
                        break
                if has_containment:
                    break
            if has_containment:
                pairs_with_contained_targets += len(targets) - 1
    
    expected_fundamental_by_target = expected_fundamental_by_pair - pairs_with_contained_targets
    expected_fundamental = max(1, expected_fundamental_by_target)
    
    pruning_potential = total_maximal - expected_fundamental
    
    pruning_stats = {
        'total_maximal': total_maximal,
        'expected_fundamental': expected_fundamental,
        'pruning_potential': pruning_potential,
        'pruning_percentage': 100 * pruning_potential / max(1, total_maximal),
        'pair_domination_potential': total_maximal - expected_fundamental_by_pair,
        'target_containment_potential': pairs_with_contained_targets
    }
    
    if verbose:
        print(f"\nPruning potential:")
        print(f"  Total maximal synergies: {total_maximal}")
        print(f"  Expected fundamental (by pair domination): ~{expected_fundamental_by_pair}")
        print(f"  Expected fundamental (with target containment): ~{expected_fundamental}")
        print(f"  Potential pruning by pair domination: ~{pruning_stats['pair_domination_potential']}")
        print(f"  Potential pruning by target containment: ~{pruning_stats['target_containment_potential']}")
        print(f"  Total potential pruning: ~{pruning_potential} "
              f"({pruning_stats['pruning_percentage']:.1f}%)")
    
    return synergies_by_level, target_to_producers, pair_to_targets, pruning_stats


# ============================================================================
# VERIFICATION FUNCTIONS
# ============================================================================

def verify_target_containment_correctness(fundamental_synergies, hierarchy):
    """
    Verify that no synergy has a target that contains another target from the same pair.
    """
    print("\n" + "="*60)
    print("TARGET CONTAINMENT VERIFICATION")
    print("="*60)
    
    synergies_by_pair = defaultdict(list)
    for syn in fundamental_synergies:
        pair_key = tuple(sorted(syn.rlabel))
        synergies_by_pair[pair_key].append(syn)
    
    violations = []
    
    for pair_key, syns in synergies_by_pair.items():
        if len(syns) <= 1:
            continue
        
        for i, syn1 in enumerate(syns):
            for syn2 in syns[i+1:]:
                t1 = syn1.plabel
                t2 = syn2.plabel
                
                if target_is_contained_by(t1, t2, hierarchy):
                    violations.append({
                        'pair': pair_key,
                        'smaller_target': t1,
                        'larger_target': t2,
                        'issue': f"{t1} ⊂ {t2}, so {t2} should not be fundamental"
                    })
                
                if target_is_contained_by(t2, t1, hierarchy):
                    violations.append({
                        'pair': pair_key,
                        'smaller_target': t2,
                        'larger_target': t1,
                        'issue': f"{t2} ⊂ {t1}, so {t1} should not be fundamental"
                    })
    
    if violations:
        print(f"\n⚠️  Found {len(violations)} target containment violations!")
        for i, v in enumerate(violations[:3]):
            print(f"\n  Violation {i+1}:")
            print(f"    Pair: {v['pair'][0]} + {v['pair'][1]}")
            print(f"    Issue: {v['issue']}")
        if len(violations) > 3:
            print(f"\n  ... and {len(violations) - 3} more violations")
        return False
    else:
        print("\n✓ No target containment violations found!")
        return True


def verify_pair_domination_correctness(fundamental_synergies, hierarchy):
    """
    Verify that no synergy has a pair that's dominated by another pair for the same target.
    """
    print("\n" + "="*60)
    print("PAIR DOMINATION VERIFICATION")
    print("="*60)
    
    synergies_by_target = defaultdict(list)
    for syn in fundamental_synergies:
        synergies_by_target[syn.plabel].append(syn)
    
    violations = []
    
    for target, syns in synergies_by_target.items():
        if len(syns) <= 1:
            continue
        
        for i, syn1 in enumerate(syns):
            for syn2 in syns[i+1:]:
                pair1 = tuple(sorted(syn1.rlabel))
                pair2 = tuple(sorted(syn2.rlabel))
                
                if pair1 == pair2:
                    continue
                
                if is_more_fundamental_than(pair1, pair2, hierarchy):
                    violations.append({
                        'target': target,
                        'more_fundamental_pair': pair1,
                        'less_fundamental_pair': pair2
                    })
                
                elif is_more_fundamental_than(pair2, pair1, hierarchy):
                    violations.append({
                        'target': target,
                        'more_fundamental_pair': pair2,
                        'less_fundamental_pair': pair1
                    })
    
    if violations:
        print(f"\n⚠️  Found {len(violations)} pair domination violations!")
        for i, v in enumerate(violations[:3]):
            print(f"\n  Violation {i+1}:")
            print(f"    Target: {v['target']}")
            print(f"    More fundamental: {v['more_fundamental_pair']}")
            print(f"    Less fundamental: {v['less_fundamental_pair']}")
        if len(violations) > 3:
            print(f"\n  ... and {len(violations) - 3} more violations")
        return False
    else:
        print("\n✓ No pair domination violations found!")
        return True


# ============================================================================
# COMPARISON FUNCTION
# ============================================================================

def compare_methods(ercs, hierarchy, RN, run_brute_force=True):
    """Compare the brute force and optimized methods with detailed analysis."""
    print("\n" + "="*70)
    print(" "*20 + "SYNERGY DETECTION COMPARISON")
    print("="*70)
    
    # First, analyze structure
    print("\n[1] Analyzing synergy structure...")
    synergies_by_level, target_to_producers, pair_to_targets, pruning_stats = analyze_synergy_structure(
        ercs, hierarchy, RN, verbose=True
    )
    
    results = {}
    
    # Brute force method
    if run_brute_force:
        print("\n[2] Running brute force method...")
        start = time.time()
        brute_force_result = get_all_fundamental_synergies_brute_force(
            ercs, hierarchy, RN, verbose=False
        )
        brute_force_time = time.time() - start
        
        print(f"   ⏱️  Time: {brute_force_time:.3f}s")
        print(f"   ✓ Found: {len(brute_force_result)} fundamental synergies")
        
        results['brute_force'] = {
            'synergies': brute_force_result,
            'time': brute_force_time,
            'count': len(brute_force_result)
        }
        
        print("\n   Verifying brute force correctness...")
        bf_target_ok = verify_target_containment_correctness(brute_force_result, hierarchy)
        bf_pair_ok = verify_pair_domination_correctness(brute_force_result, hierarchy)
        
        if not (bf_target_ok and bf_pair_ok):
            print("\n   ⚠️  WARNING: Brute force has violations!")
    else:
        print("\n[2] Skipping brute force method (run_brute_force=False)")
        results['brute_force'] = None
    
    # Optimized method
    print("\n[3] Running optimized method...")
    start = time.time()
    optimized_result = get_all_fundamental_synergies_optimized(
        ercs, hierarchy, RN, verbose=True
    )
    optimized_time = time.time() - start
    
    print(f"   ⏱️  Time: {optimized_time:.3f}s")
    print(f"   ✓ Found: {len(optimized_result)} fundamental synergies")
    
    results['optimized'] = {
        'synergies': optimized_result,
        'time': optimized_time,
        'count': len(optimized_result)
    }
    
    print("\n   Verifying optimized correctness...")
    opt_target_ok = verify_target_containment_correctness(optimized_result, hierarchy)
    opt_pair_ok = verify_pair_domination_correctness(optimized_result, hierarchy)
    
    if not (opt_target_ok and opt_pair_ok):
        print("\n   ⚠️  WARNING: Optimized method has violations!")
    else:
        print("\n   ✓✓ Optimized method passed all correctness checks!")
    
    # Compare results
    print("\n" + "="*70)
    print(" "*25 + "COMPARISON RESULTS")
    print("="*70)
    
    if run_brute_force:
        print(f"\n⏱️  Performance:")
        print(f"   Brute force time:  {brute_force_time:.3f}s")
        print(f"   Optimized time:    {optimized_time:.3f}s")
        
        if optimized_time > 0:
            speedup = brute_force_time / optimized_time
            print(f"   Speedup:           {speedup:.2f}x")
    
    print(f"\n📊 Counts:")
    if run_brute_force:
        print(f"   Brute force found: {len(brute_force_result)}")
    print(f"   Optimized found:   {len(optimized_result)}")
    
    # Detailed comparison
    if run_brute_force:
        print(f"\n" + "="*70)
        print(" "*20 + "DETAILED DIFFERENCE ANALYSIS")
        print("="*70)
        
        brute_set = {(tuple(sorted(s.rlabel)), s.plabel) for s in brute_force_result}
        opt_set = {(tuple(sorted(s.rlabel)), s.plabel) for s in optimized_result}
        
        if brute_set == opt_set:
            print(f"\n✓✓✓ Results match perfectly!")
        else:
            print(f"\n⚠️⚠️⚠️ WARNING: Results differ!")
            
            only_brute = brute_set - opt_set
            only_opt = opt_set - brute_set
            
            if only_brute:
                print(f"\n{'='*70}")
                print(f"FOUND IN BRUTE FORCE BUT NOT IN OPTIMIZED: {len(only_brute)}")
                print(f"{'='*70}")
                
                for idx, (pair, target) in enumerate(only_brute, 1):
                    analyze_synergy_mismatch(
                        pair, target, hierarchy, RN,
                        context=f"[{idx}/{len(only_brute)}]"
                    )
            
            if only_opt:
                print(f"\n{'='*70}")
                print(f"FOUND IN OPTIMIZED BUT NOT IN BRUTE FORCE: {len(only_opt)}")
                print(f"{'='*70}")
                
                for idx, (pair, target) in enumerate(only_opt, 1):
                    analyze_synergy_mismatch(
                        pair, target, hierarchy, RN,
                        context=f"[{idx}/{len(only_opt)}]"
                    )
    
    # Show pruning effectiveness
    print(f"\n" + "="*70)
    print(f"🎯 PRUNING EFFECTIVENESS")
    print(f"="*70)
    print(f"   Total pruning potential: {pruning_stats['pruning_percentage']:.1f}%")
    print(f"   - By pair domination: {pruning_stats['pair_domination_potential']} synergies")
    print(f"   - By target containment: {pruning_stats['target_containment_potential']} synergies")
    
    print("\n" + "="*70 + "\n")
    
    return results


# ============================================================================
# MAIN TEST FUNCTION
# ============================================================================

def test_on_example_system(RN_or_ercs, run_brute_force=True):
    """Test the optimization on a reaction network or list of ERCs."""
    if isinstance(RN_or_ercs, tuple):
        ercs, hierarchy, RN = RN_or_ercs
    else:
        RN = RN_or_ercs
        print("Building ERCs from reaction network...")
        ercs = ERC.ERCs(RN)
        print(f"Found {len(ercs)} ERCs")
        print("Building hierarchy...")
        hierarchy = ERC_Hierarchy(ercs, RN)
        print(f"Hierarchy has {hierarchy.graph.number_of_nodes()} nodes and "
              f"{hierarchy.graph.number_of_edges()} edges")

    results = compare_methods(ercs, hierarchy, RN, run_brute_force=run_brute_force)
    
    return results


# ============================================================================
# RUN TESTS
# ============================================================================

if __name__ == "__main__":
    file_path = 'networks/RandomAlife/RN_Ns_40_Norg_10_id_568.txt'
    
    print("Loading reaction network from:", file_path)
    RN = read_txt(file_path)
    
    print(f"\nReaction network stats:")
    print(f"  Species: {len(RN.species())}")
    print(f"  Reactions: {len(RN.reactions())}")
    
    results = test_on_example_system(RN, run_brute_force=True)