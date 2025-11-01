#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script for ERC synergy and complementarity computation
Reports statistics and results from ERC_SORN computation
Can optionally generate a connected subnetwork for testing
"""

import time

# Import the necessary modules
from pyCOT.ERC_Hierarchy import (
    ERC,
    ERC_Hierarchy
)
from pyCOT.ERC_Synergy_Complementarity import (
    build_erc_sorn
)
from pyCOT.io.functions import read_txt


def analyze_synergies_and_complementarities(RN):
    """Compute and report synergies and complementarities for a reaction network."""
    
    print("\n" + "="*70)
    print(" "*15 + "SYNERGY & COMPLEMENTARITY ANALYSIS")
    print("="*70)
    
    # Build ERCs and hierarchy
    print("\n[1] Building ERCs...")
    start = time.time()
    ercs = ERC.ERCs(RN)
    erc_time = time.time() - start
    print(f"    Found {len(ercs)} ERCs in {erc_time:.3f}s")
    
    print("\n[2] Building ERC Hierarchy...")
    start = time.time()
    hierarchy = ERC_Hierarchy(RN, ercs)
    hierarchy_time = time.time() - start
    print(f"    Hierarchy has {hierarchy.graph.number_of_nodes()} nodes and "
          f"{hierarchy.graph.number_of_edges()} edges")
    print(f"    Built in {hierarchy_time:.3f}s")
    
    # Build SORN (computes synergies and complementarities)
    print("\n[3] Computing Synergies and Complementarities (ERC_SORN)...")
    sorn = build_erc_sorn(hierarchy, RN)
    
    # Get statistics
    stats = sorn.get_statistics()
    
    # Display results
    print("\n" + "="*70)
    print(" "*28 + "RESULTS")
    print("="*70)
    
    print(f"\n⏱️  Timing:")
    print(f"   ERCs:           {erc_time:.3f}s")
    print(f"   Hierarchy:      {hierarchy_time:.3f}s")
    print(f"   SORN:           {stats['build_time']:.3f}s")
    total_time = erc_time + hierarchy_time + stats['build_time']
    print(f"   Total:          {total_time:.3f}s")
    
    print(f"\n📊 ERC Pair Statistics:")
    print(f"   Total pairs checked:     {stats['total_pairs_checked']}")
    print(f"   Productive pairs:        {stats['productive_pairs']} "
          f"({100*stats['productive_pairs']/max(1,stats['total_pairs_checked']):.1f}%)")
    print(f"   - Synergistic pairs:     {stats['synergistic_pairs']}")
    print(f"   - Complementary pairs:   {stats['complementary_pairs']}")
    
    print(f"\n🔄 Synergy Statistics:")
    print(f"   Total synergies:         {stats['total_synergies']}")
    if stats['synergistic_pairs'] > 0:
        avg_syn = stats['total_synergies'] / stats['synergistic_pairs']
        print(f"   Avg per synergistic pair: {avg_syn:.2f}")
    
    print(f"\n🤝 Complementarity Statistics:")
    print(f"   Total complementarities: {stats['total_complementarities']}")
    if stats['complementary_pairs'] > 0:
        avg_comp = stats['total_complementarities'] / stats['complementary_pairs']
        print(f"   Avg per complementary pair: {avg_comp:.2f}")
    
    # Display sample synergies
    print(f"\n💡 Sample Synergies (first 10):")
    synergy_count = 0
    for erc1_label in list(sorn._erc_to_synergies.keys())[:5]:
        partners = sorn._erc_to_synergies[erc1_label]
        for erc2_label, synergies in partners:
            for syn in synergies[:2]:  # Show up to 2 per pair
                print(f"   {syn}")
                synergy_count += 1
                if synergy_count >= 10:
                    break
            if synergy_count >= 10:
                break
        if synergy_count >= 10:
            break
    
    if stats['total_synergies'] > 10:
        print(f"   ... and {stats['total_synergies'] - 10} more")
    
    # Display sample complementarities
    print(f"\n🤝 Sample Complementarities (first 10):")
    comp_count = 0
    for erc1_label in list(sorn._erc_to_complementarities.keys())[:5]:
        partners = sorn._erc_to_complementarities[erc1_label]
        for erc2_label, comps in partners:
            for comp in comps[:2]:  # Show up to 2 per pair
                print(f"   {comp}")
                comp_count += 1
                if comp_count >= 10:
                    break
            if comp_count >= 10:
                break
        if comp_count >= 10:
            break
    
    if stats['total_complementarities'] > 10:
        print(f"   ... and {stats['total_complementarities'] - 10} more")
    
    # Theoretical maximum synergies
    n_ercs = len(ercs)
    max_pairs = n_ercs * (n_ercs - 1) // 2
    max_synergies_upper_bound = max_pairs * n_ercs  # Each pair could produce synergies to all ERCs
    
    print(f"\n📈 Scaling Analysis:")
    print(f"   Number of ERCs:          {n_ercs}")
    print(f"   Possible ERC pairs:      {max_pairs}")
    print(f"   Theoretical upper bound: {max_synergies_upper_bound} (all pairs × all targets)")
    if stats['total_synergies'] > 0:
        efficiency = 100 * stats['total_synergies'] / max_synergies_upper_bound
        print(f"   Actual synergies:        {stats['total_synergies']} ({efficiency:.3f}% of upper bound)")
    
    print("\n" + "="*70 + "\n")
    
    return {
        'n_ercs': n_ercs,
        'erc_time': erc_time,
        'hierarchy_time': hierarchy_time,
        'sorn_time': stats['build_time'],
        'total_time': total_time,
        'stats': stats,
        'sorn': sorn
    }


if __name__ == "__main__":
    # Configuration
    USE_SUBNETWORK = True  # Set to True to generate a subnetwork
    SUBNETWORK_PERCENTAGE = 1.0  # Target percentage of reactions (1.0 = 100%)
    SUBNETWORK_SEED = 47  # For reproducibility
    
    # Test files - uncomment the one you want to use
    #file_path = 'networks/testing/Farm.txt'
    # file_path = 'networks/Navarino/RN_IN_05.txt'
    # file_path = 'networks/RandomAlife/RN_Ns_40_Norg_10_id_568.txt'
    # file_path = 'networks/RandomAlife/RN_Ns_40_Norg_20_id_396.txt'
    # file_path = 'networks/biomodels_interesting/central_ecoli.txt'
    # file_path = 'networks/biomodels_interesting/BIOMD0000000652_manyOrgs.txt'
    # file_path = 'networks/biomodels_interesting/bigg_iAF692.txt'
    file_path = 'networks/Marine_Ecosystem/Las_Cruces_251021.txt'

    print("="*70)
    print(" "*10 + "ERC SYNERGY & COMPLEMENTARITY TEST")
    print("="*70)
    print(f"\nLoading reaction network from: {file_path}")
    
    RN_full = read_txt(file_path)
    
    print(f"\nOriginal network stats:")
    print(f"  Species:   {len(RN_full.species())}")
    print(f"  Reactions: {len(RN_full.reactions())}")
    
    # Generate subnetwork if requested
    if USE_SUBNETWORK and SUBNETWORK_PERCENTAGE < 1.0:
        print("\n" + "="*70)
        print("GENERATING CONNECTED SUBNETWORK")
        print("="*70)
        RN = RN_full.generate_connected_subnetwork(
            target_percentage=SUBNETWORK_PERCENTAGE,
            seed=SUBNETWORK_SEED,
            verbose=True
        )
    else:
        RN = RN_full
    
    # Run analysis
    results = analyze_synergies_and_complementarities(RN)
    
    # Print final summary
    print("="*70)
    print(" "*25 + "SUMMARY")
    print("="*70)
    print(f"Network: {file_path}")
    print(f"Species: {len(RN.species())}, Reactions: {len(RN.reactions())}")
    print(f"ERCs: {results['n_ercs']}")
    print(f"Synergies: {results['stats']['total_synergies']}")
    print(f"Complementarities: {results['stats']['total_complementarities']}")
    print(f"Total time: {results['total_time']:.3f}s")
    print("="*70)