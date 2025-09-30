#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script for SORN_Generators.py

Tests the construction of irreducible generators through productive novelty
and the ERC_SORN functionality for efficient relationship lookup.
"""
import os
import sys
import time

# Add the root directory to the PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import required modules
from pyCOT.rn_rustworkx import ReactionNetwork
from pyCOT.io.functions import read_txt
from pyCOT.ERC_Hierarchy import ERC, ERC_Hierarchy, species_list_to_names
from pyCOT.ERC_Synergy_Complementarity import *
# Import our new SORN_Generators module
from pyCOT.SORN_Generators import (
    IrreducibleGenerator,
    ProductiveExtension,
    identify_p_ercs,
    build_erc_sorn,
    find_productive_extensions,
    build_irreducible_generators,
    analyze_generator_statistics
)

def test_sorn_generators():
    """Main test function for SORN_Generators module."""
    print("="*70)
    print("TESTING SORN_GENERATORS.PY")
    print("="*70)
    
    # Load reaction network
    print("Loading reaction network...")
    file_path = 'networks/testing/Farm.txt'
    #file_path = 'networks/Navarino/RN_IN_05.txt'
    # Alternative networks for testing:
    # file_path = 'networks/RandomAlife/RN_Ns_20_Norg_4_id_12.txt'
    #file_path = 'networks/RandomAlife/RN_Ns_40_Norg_12_id_358.txt'
    #file_path= 'networks/biomodels_interesting/bigg_iAF692.txt'  # Adjust path as needed


    RN = read_txt(file_path)
    
    # Verify we got a proper ReactionNetwork object
    if not isinstance(RN, ReactionNetwork):
        raise TypeError(f"read_txt returned {type(RN)}, expected ReactionNetwork")
    
    print(f"✅ Loaded network: {len(RN.species())} species, {len(RN.reactions())} reactions")
    
    # Create ERC hierarchy
    print("\nCreating ERC hierarchy...")
    start_time = time.time()
    ercs = ERC.ERCs(RN)
    hierarchy = ERC_Hierarchy(ercs, RN)
    hierarchy_time = time.time() - start_time
    print(f"✅ Created hierarchy: {len(hierarchy.ercs)} ERCs in {hierarchy_time:.2f}s")
    
    # Test 1: Identify P-ERCs
    print("\n" + "-"*50)
    print("TEST 1: Identifying P-ERCs")
    print("-"*50)
    
    start_time = time.time()
    p_ercs = identify_p_ercs(hierarchy, RN)
    p_erc_time = time.time() - start_time
    
    print(f"✅ Found {len(p_ercs)} P-ERCs in {p_erc_time:.3f}s")
    
    if p_ercs:
        print("P-ERCs found:")
        for i, p_erc in enumerate(p_ercs[:5]):  # Show first 5
            closure_species=species_list_to_names(p_erc.get_closure(RN)) 
            print(p_erc.label + ":" + str(closure_species)) 
            print(f"  Minimal Generators: {species_list_to_names(p_erc.min_generators)}")
            print(str(len(closure_species))+ " species in" f"  {i+1}")
        if len(p_ercs) > 10:
            print(f"  ... and {len(p_ercs) - 10} more")
    else:
        print("No P-ERCs found")
    
    # Test 2: Build ERC_SORN
    print("\n" + "-"*50)
    print("TEST 2: Building ERC_SORN")
    print("-"*50)
    
    start_time = time.time()
    erc_sorn = build_erc_sorn(hierarchy, RN)
    sorn_time = time.time() - start_time
    
    sorn_stats = erc_sorn.get_statistics()
    print(f"✅ ERC_SORN built in {sorn_time:.2f}s")
    print(f"   - Total pairs checked: {sorn_stats['total_pairs_checked']}")
    print(f"   - Productive pairs: {sorn_stats['productive_pairs']}")
    print(f"   - Synergistic pairs: {sorn_stats['synergistic_pairs']}")
    print(f"   - Complementary pairs: {sorn_stats['complementary_pairs']}")
    print(f"   - Total synergies: {sorn_stats['total_synergies']}")
    print(f"   - Total complementarities: {sorn_stats['total_complementarities']}")
    
    # Test some SORN queries
    if len(hierarchy.ercs) >= 2:
        erc1, erc2 = hierarchy.ercs[0], hierarchy.ercs[1]
        print(f"\nTesting SORN queries between {erc1.label} and {erc2.label}:")
        
        synergies = erc_sorn.get_synergies(erc1.label, erc2.label)
        complementarities = erc_sorn.get_complementarities(erc1.label, erc2.label)
        
        print(f"   - Synergies: {len(synergies)}")
        print(f"   - Complementarities: {len(complementarities)}")
        print(f"   - Has productive relationship: {erc_sorn.has_productive_relationship(erc1.label, erc2.label)}")
        
        # Show productive partners for first ERC
        partners = erc_sorn.get_productive_partners(erc1.label)
        print(f"   - {erc1.label} has {len(partners)} productive partners")
    
    # Test 3: Build Irreducible Generators
    print("\n" + "-"*50)
    print("TEST 3: Building Irreducible Generators")
    print("-"*50)
    
    start_time = time.time()
    generators, erc_sorn = build_irreducible_generators(hierarchy, RN, erc_sorn=erc_sorn, max_size=50, verbose=True)
    generator_time = time.time() - start_time
    
    print(f"\n✅ Built {len(generators)} irreducible generators in {generator_time:.2f}s")
    
    if generators:
        print("\nFirst few generators:")
        for i, gen in enumerate(generators[:5]):
            closure_size = len(gen.get_closure(RN))
            ssm_status = "SSM" if gen.check_ssm(RN) else "Non-SSM"
            erc_labels = gen.get_erc_labels()
            print(f"  {i+1}. {erc_labels} → {closure_size} species ({ssm_status})")
        if len(generators) > 5:
            print(f"  ... and {len(generators) - 5} more")
    
    # Test 4: Test productive extensions
    print("\n" + "-"*50)
    print("TEST 4: Testing Productive Extensions")
    print("-"*50)
    
    if generators:
        # Test with the first generator
        test_gen = generators[0]
        print(f"Testing extensions for generator: {test_gen.get_erc_labels()}")
        
        extensions = find_productive_extensions(test_gen, erc_sorn)
        print(f"✅ Found {len(extensions)} productive extensions")
        
        if extensions:
            print("Extension possibilities:")
            for i, ext in enumerate(extensions[:5]):
                print(f"  {i+1}. Add {ext.target_erc.label} via {ext.step_type}")
            if len(extensions) > 5:
                print(f"  ... and {len(extensions) - 5} more")
    
    # Test 5: Analyze Generator Statistics
    print("\n" + "-"*50)
    print("TEST 5: Analyzing Generator Statistics")
    print("-"*50)
    
    if generators:
        stats = analyze_generator_statistics(generators)
        
        print(f"✅ Generator Statistics:")
        print(f"   - Total generators: {stats['total_generators']}")
        print(f"   - Size distribution: {dict(stats['size_distribution'])}")
        print(f"   - Average size: {stats['summary']['avg_size']:.2f}")
        print(f"   - Maximum size: {stats['summary']['max_size']}")
        print(f"   - Average synergies per generator: {stats['summary']['avg_synergies']:.2f}")
        
        if stats['summary']['most_used_ercs']:
            print("   - Most used ERCs:")
            for erc_label, count in stats['summary']['most_used_ercs']:
                print(f"     • {erc_label}: used in {count} generators")
    
    # Test 6: Test IrreducibleGenerator class directly
    print("\n" + "-"*50)
    print("TEST 6: Testing IrreducibleGenerator Class")
    print("-"*50)
    
    if p_ercs:
        # Create a simple generator from a P-ERC
        test_gen = IrreducibleGenerator([p_ercs[2]])
        print(f"✅ Created generator: {test_gen}")
        
        closure = test_gen.get_closure(RN)
        print(f"   - Closure size: {len(closure)}")
        print(f"   - Is SSM: {test_gen.check_ssm(RN)}")
        print(f"   - Size: {test_gen.size()}")
        print(f"   - Synergy count: {test_gen.get_synergy_count()}")
        print(f"   - Complementarity counts: {test_gen.get_complementarity_counts()}")
        
        # Test adding a productive step (if we have extensions)
        extensions = find_productive_extensions(test_gen, erc_sorn)
        if extensions:
            ext = extensions[0]
            test_gen.add_productive_step(ext.target_erc, ext.step_type, ext.step_details)
            print(f"   - After extension: {test_gen}")
            print(f"   - New size: {test_gen.size()}")
    
    # Summary
    print("\n" + "="*70)
    print("TEST SUMMARY")
    print("="*70)
    print(f"✅ All tests completed successfully!")
    print(f"✅ Network: {len(RN.species())} species, {len(RN.reactions())} reactions")
    print(f"✅ ERCs: {len(hierarchy.ercs)}")
    print(f"✅ P-ERCs: {len(p_ercs)}")
    print(f"✅ Productive pairs: {sorn_stats['productive_pairs']}")
    print(f"✅ Irreducible generators: {len(generators)}")
    
    total_time = hierarchy_time + p_erc_time + sorn_time + generator_time
    print(f"✅ Total computation time: {total_time:.2f}s")
    
    return {
        'RN': RN,
        'hierarchy': hierarchy,
        'p_ercs': p_ercs,
        'erc_sorn': erc_sorn,
        'generators': generators,
        'statistics': stats if generators else {},
        'timing': {
            'hierarchy': hierarchy_time,
            'p_ercs': p_erc_time,
            'sorn': sorn_time,
            'generators': generator_time,
            'total': total_time
        }
    }

if __name__ == "__main__":
    try:
        results = test_sorn_generators()
        print("\n🎉 All tests passed!")
        
        # Optional: Save results for further analysis
        generators = results['generators']
        if generators:
            ssm_count = sum(1 for gen in generators if gen.check_ssm(results['RN']))
            print(f"\n📊 Additional insights:")
            print(f"   - SSM generators: {ssm_count}/{len(generators)}")
            print(f"   - SSM ratio: {ssm_count/len(generators)*100:.1f}%")
        
    except Exception as e:
        print(f"\n❌ Test failed with error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)