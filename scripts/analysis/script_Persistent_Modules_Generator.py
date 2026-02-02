#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script for Persistent_Modules_Generator.py

Tests the construction of elementary and non-elementary organizations
from irreducible generators, including self-maintenance verification
and hierarchical organization construction.
"""
import os
import sys
import time
import numpy as np

from pyCOT.tests.reaction_network.test_add_from_reaction_string import rn

# Add the root directory to the PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import required modules
from pyCOT.core.rn_rustworkx import ReactionNetwork
from pyCOT.io.functions import read_txt
from pyCOT.analysis.ERC_Hierarchy import ERC, ERC_Hierarchy
from pyCOT.analysis.ERC_Synergy_Complementarity import build_erc_sorn
# Import SORN_Generators for building generators
from pyCOT.analysis.SORN_Generators import (
    build_irreducible_generators
)

# Import our new Persistent_Modules_Generator module
from pyCOT.analysis.Persistent_Modules_Generator import (
    ElementarySO,
    Organization, 
    OrganizationHierarchy,
    compute_elementary_sos,
    check_self_maintenance,
    find_productive_so_extensions,
    build_non_elementary_organizations,
    compute_all_organizations
)

def test_persistent_modules_generator():
    """Main test function for Persistent_Modules_Generator module."""
    print("="*80)
    print("TESTING PERSISTENT_MODULES_GENERATOR.PY")
    print("="*80)
    
    # Load reaction network
    print("Loading reaction network...")
    file_path = 'networks/testing/Farm.txt'
    file_path = 'networks/Marine_Ecosystem/Las_Cruces_251021.txt'  # Specify the network file

    # Alternative networks for testing:
    #file_path = 'networks/RandomAlife/RN_Ns_40_Norg_12_id_358.txt'
    #file_path = 'networks/Navarino/RN_IN_05.txt'
    #file_path = 'networks/biomodels_interesting/bigg_iAF692.txt'
    RN = read_txt(file_path)
    
    # Verify we got a proper ReactionNetwork object
    if not isinstance(RN, ReactionNetwork): 
        raise TypeError(f"read_txt returned {type(RN)}, expected ReactionNetwork")
    
    print(f"✅ Loaded network: {len(RN.species())} species, {len(RN.reactions())} reactions")
    
    # Create ERC hierarchy and build generators (prerequisite)
    print("\nSetting up prerequisites...")
    start_time = time.time()
    ercs = ERC.ERCs(RN)
    hierarchy = ERC_Hierarchy(RN, ercs)
    setup_time = time.time() - start_time
    print(f"✅ Created hierarchy: {len(hierarchy.ercs)} ERCs in {setup_time:.2f}s")
        # Test 2: Build ERC_SORN
    print("\n" + "-"*50)
    print("TEST 2: Building ERC_SORN")
    print("-"*50)
    
    start_time = time.time()
    erc_sorn = build_erc_sorn(hierarchy, RN)
    sorn_time = time.time() - start_time

    # Build generators using SORN_Generators
    print("Building irreducible generators...")
    start_time = time.time()
    generators, erc_sorn = build_irreducible_generators(hierarchy, RN, erc_sorn=erc_sorn, max_size=50, verbose=False)
    generators_time = time.time() - start_time
    print(f"✅ Built {len(generators)} generators in {generators_time:.2f}s")
    
    if not generators:
        print("❌ No generators found. Cannot proceed with organization tests.")
        return None
    
    # Count SSM generators for context
    ssm_count = sum(1 for gen in generators if gen.check_ssm(RN))
    print(f"   - SSM generators: {ssm_count}/{len(generators)}")
    
    # Test 1: Compute Elementary Semi-Organizations
    print("\n" + "-"*60)
    print("TEST 1: Computing Elementary Semi-Organizations")
    print("-"*60)
    
    start_time = time.time()
    elementary_sos = compute_elementary_sos(generators, RN, hierarchy)
    esos_time = time.time() - start_time
    
    print(f"✅ Found {len(elementary_sos)} elementary SOs in {esos_time:.3f}s")
    
    if elementary_sos:
        print("Elementary SOs found:")
        for i, eso in enumerate(elementary_sos[:5]):  # Show first 5
            stats = eso.get_generation_statistics()
            print(f"  {i+1}. {eso.id}:")
            print(f"     - Species: {len(eso.closure_species)}")
            print(f"     - ERCs: {len(eso.constituent_ercs)}")
            print(f"     - P-ERC: {eso.is_p_erc}")
            print(f"     - Pathways: {stats['pathway_count']}")
        if len(elementary_sos) > 5:
            print(f"  ... and {len(elementary_sos) - 5} more")
    else:
        print("No elementary SOs found")
        return None
    
    # Test 2: Check Self-Maintenance
    print("\n" + "-"*60)
    print("TEST 2: Checking Self-Maintenance")
    print("-"*60)
    
    elementary_organizations = []
    sm_times = []
    
    print("Testing self-maintenance for elementary SOs...")
    for i, eso in enumerate(elementary_sos):
        start_time = time.time()
        is_org, flux = eso.check_self_maintenance()
        sm_time = time.time() - start_time
        sm_times.append(sm_time)
        
        if is_org:
            elementary_organizations.append(eso)
            status = "✅ Organization"
        else:
            status = "❌ Semi-organization only"
        
        if i < 5:  # Show details for first 5
            print(f"  {i+1}. {eso.id}: {status} ({sm_time:.3f}s)")
    
    if len(elementary_sos) > 5:
        print(f"  ... and {len(elementary_sos) - 5} more tested")
    
    total_sm_time = sum(sm_times)
    print(f"✅ Self-maintenance check completed in {total_sm_time:.3f}s")
    print(f"   - Elementary organizations: {len(elementary_organizations)}")
    print(f"   - Semi-organizations only: {len(elementary_sos) - len(elementary_organizations)}")
    
    if elementary_sos:
        conversion_rate = len(elementary_organizations) / len(elementary_sos) * 100
        print(f"   - Conversion rate: {conversion_rate:.1f}%")
    
    # Test 3: Test Self-Maintenance Function Directly
    print("\n" + "-"*60)
    print("TEST 3: Testing Self-Maintenance Function Directly")
    print("-"*60)
    
    if elementary_sos:
        test_eso = elementary_sos[0]
        print(f"Testing self-maintenance for {test_eso.id}...")
        
        start_time = time.time()
        is_sm, flux_vector, production_vector = check_self_maintenance(test_eso.closure_species, RN)
        direct_sm_time = time.time() - start_time
        
        print(f"✅ Direct self-maintenance test completed in {direct_sm_time:.4f}s")
        print(f"   - Is self-maintaining: {is_sm}")
        
        if is_sm and flux_vector is not None:
            print(f"   - Flux vector shape: {flux_vector.shape}")
            print(f"   - Flux vector sum: {np.sum(flux_vector):.4f}")
            print(f"   - Production vector shape: {production_vector.shape}")
            print(f"   - Min production: {np.min(production_vector):.6f}")
            print(f"   - Max production: {np.max(production_vector):.6f}")
    
    # Test 4: Find Productive SO Extensions
    print("\n" + "-"*60)
    print("TEST 4: Finding Productive SO Extensions")
    print("-"*60)
    
    if len(elementary_organizations) >= 2:
        print(f"Finding productive extensions for {len(elementary_organizations)} elementary organizations...")
        
        # Build ERC_SORN for extension finding
        start_time = time.time()
        erc_sorn = build_erc_sorn(hierarchy, RN)
        sorn_time = time.time() - start_time
        print(f"✅ Built ERC_SORN in {sorn_time:.2f}s")
        
        start_time = time.time()
        productive_combinations = find_productive_so_extensions(elementary_organizations, erc_sorn)
        extensions_time = time.time() - start_time
        
        print(f"✅ Found {len(productive_combinations)} productive combinations in {extensions_time:.3f}s")
        
        if productive_combinations:
            print("Productive combinations found:")
            for i, ((so1_id, so2_id), relationships) in enumerate(list(productive_combinations.items())[:3]):
                print(f"  {i+1}. {so1_id} + {so2_id}:")
                print(f"     - {len(relationships)} productive relationships")
                # Show relationship types
                rel_types = {}
                for rel in relationships:
                    rel_type = rel['type']
                    rel_types[rel_type] = rel_types.get(rel_type, 0) + 1
                for rel_type, count in rel_types.items():
                    print(f"       • {count} {rel_type}(s)")
            if len(productive_combinations) > 3:
                print(f"  ... and {len(productive_combinations) - 3} more")
    else:
        print("❌ Need at least 2 elementary organizations to find productive extensions")
        productive_combinations = {}
        erc_sorn = build_erc_sorn(hierarchy, RN)
    
    # Test 5: Build Non-Elementary Organizations
    print("\n" + "-"*60)
    print("TEST 5: Building Non-Elementary Organizations")
    print("-"*60)
    
    if len(elementary_organizations) >= 2:
        start_time = time.time()
        non_elementary_orgs = build_non_elementary_organizations(
            elementary_organizations, erc_sorn, max_size=50, verbose=True)
        non_elem_time = time.time() - start_time
        
        print(f"✅ Built {len(non_elementary_orgs)} non-elementary organizations in {non_elem_time:.3f}s")
        
        if non_elementary_orgs:
            print("Non-elementary organizations:")
            for i, org in enumerate(non_elementary_orgs[:3]):
                print(f"  {i+1}. {org.id}:")
                print(f"     - Elementary SOs: {len(org.elementary_sos)}")
                print(f"     - Total species: {len(org.combined_closure)}")
                print(f"     - Total ERCs: {len(org.constituent_ercs)}")
                print(f"     - Construction method: {org.construction_method}")
            if len(non_elementary_orgs) > 3:
                print(f"  ... and {len(non_elementary_orgs) - 3} more")
    else:
        print("❌ Need at least 2 elementary organizations to build non-elementary ones")
        non_elementary_orgs = []
        non_elem_time = 0
    
    # Test 6: Test Organization Hierarchy
    print("\n" + "-"*60)
    print("TEST 6: Testing Organization Hierarchy")
    print("-"*60)
    
    print("Building organization hierarchy...")
    start_time = time.time()
    
    # Create hierarchy object
    org_hierarchy = OrganizationHierarchy()
    
    # Add all organizations
    for eso in elementary_organizations:
        org_hierarchy.add_organization(eso)
    for org in non_elementary_orgs:
        org_hierarchy.add_organization(org)
    
    # Build containment relationships
    org_hierarchy.build_containment_relationships()
    hierarchy_build_time = time.time() - start_time
    
    print(f"✅ Built organization hierarchy in {hierarchy_build_time:.3f}s")
    
    # Get hierarchy statistics
    hierarchy_stats = org_hierarchy.get_statistics()
    print("Hierarchy statistics:")
    for key, value in hierarchy_stats.items():
        print(f"  - {key}: {value}")
    
    # Find maximal organizations
    maximal_orgs = org_hierarchy.find_maximal_organizations()
    print(f"✅ Found {len(maximal_orgs)} maximal organizations")
    
    if maximal_orgs:
        print("Maximal organizations:")
        for i, org in enumerate(maximal_orgs[:3]):
            # Get species count depending on organization type
            if hasattr(org, 'combined_closure'):
                species_count = len(org.combined_closure)
            else:  # ElementarySO
                species_count = len(org.closure_species)
            print(f"  {i+1}. {org.id}: {species_count} species")
    
    # Test 7: Complete Pipeline Test
    print("\n" + "-"*60)
    print("TEST 7: Complete Pipeline (compute_all_organizations)")
    print("-"*60)
    
    print("Running complete organization computation pipeline...")
    start_time = time.time()
    
    complete_results = compute_all_organizations(
        RN, 
        max_generator_size=50,  # Smaller for testing
        max_organization_size=10,
        verbose=True
    )
    
    complete_time = time.time() - start_time
    
    print(f"✅ Complete pipeline finished in {complete_time:.2f}s")
    
    # Show pipeline results
    print("\nPipeline Results Summary:")
    stats = complete_results['statistics']
    timing = complete_results['timing']
    
    print(f"  - Total generators: {stats['total_generators']}")
    print(f"  - Elementary SOs: {stats['total_elementary_sos']}")
    print(f"  - Elementary organizations: {stats['total_elementary_organizations']}")
    print(f"  - Non-elementary organizations: {stats['total_non_elementary_organizations']}")
    print(f"  - Total organizations: {stats['total_organizations']}")
    print(f"  - Conversion rate: {stats['conversion_rate_sos_to_orgs']*100:.1f}%")
    
    print("\nTiming breakdown:")
    for phase, time_val in timing.items():
        print(f"  - {phase}: {time_val:.3f}s")
    
    # Test 8: Test Individual Classes
    print("\n" + "-"*60)
    print("TEST 8: Testing Individual Classes")
    print("-"*60)
    
    if elementary_sos:
        # Test ElementarySO class
        test_eso = elementary_sos[0]
        print(f"Testing ElementarySO class with {test_eso.id}:")
        print(f"  - Size: {test_eso.size()}")
        print(f"  - Is P-ERC: {test_eso.is_p_erc}")
        
        # Test statistics
        stats = test_eso.get_generation_statistics()
        print(f"  - Generation statistics keys: {list(stats.keys())}")
        print(f"  - Pathway count: {stats['pathway_count']}")
        
        # Test self-maintenance
        is_org, flux = test_eso.check_self_maintenance()
        print(f"  - Self-maintaining: {is_org}")
    
    if non_elementary_orgs:
        # Test Organization class
        test_org = non_elementary_orgs[0]
        print(f"\nTesting Organization class with {test_org.id}:")
        print(f"  - Size: {test_org.size()}")
        print(f"  - Elementary SOs: {[eso.id for eso in test_org.elementary_sos]}")
        print(f"  - Combined closure size: {len(test_org.combined_closure)}")
        print(f"  - Construction method: {test_org.construction_method}")
        
        # Test self-maintenance for non-elementary
        is_org, flux = test_org.check_self_maintenance(RN)
        print(f"  - Self-maintaining: {is_org}")
    
    # Summary
    print("\n" + "="*80)
    print("TEST SUMMARY")
    print("="*80)
    print(f"✅ All tests completed successfully!")
    print(f"✅ Network: {len(RN.species())} species, {len(RN.reactions())} reactions")
    print(f"✅ Generators: {len(generators)}")
    print(f"✅ Elementary SOs: {len(elementary_sos)}")
    print(f"✅ Elementary organizations: {len(elementary_organizations)}")
    print(f"✅ Non-elementary organizations: {len(non_elementary_orgs)}")
    print(f"✅ Total organizations: {len(elementary_organizations) + len(non_elementary_orgs)}")
    
    total_test_time = (setup_time + generators_time + esos_time + 
                      total_sm_time + extensions_time + non_elem_time + 
                      hierarchy_build_time + complete_time)
    print(f"✅ Total test time: {total_test_time:.2f}s")
    
    print("Elementary Organizations")
    for org in elementary_organizations:
        species = [specie.name for specie in org.closure_species]
        print(org.id+" =" +str(species))
            

    print("Non-Elementary Organizations")
    for org in non_elementary_orgs:
        species = [specie.name for specie in org.combined_closure]
        print(org.id+" =" + str(species))
    return {
        'RN': RN,
        'hierarchy': hierarchy,
        'generators': generators,
        'elementary_sos': elementary_sos,
        'elementary_organizations': elementary_organizations,
        'non_elementary_organizations': non_elementary_orgs,
        'organization_hierarchy': org_hierarchy,
        'complete_results': complete_results,
        'timing': {
            'setup': setup_time,
            'generators': generators_time,
            'elementary_sos': esos_time,
            'self_maintenance': total_sm_time,
            'extensions': extensions_time,
            'non_elementary': non_elem_time,
            'hierarchy_build': hierarchy_build_time,
            'complete_pipeline': complete_time,
            'total': total_test_time
        }
    }

def analyze_results(results):
    """Analyze and display additional insights from the test results."""
    if not results:
        return
    
    print("\n" + "="*80)
    print("ADDITIONAL ANALYSIS")
    print("="*80)
    
    elementary_sos = results['elementary_sos']
    elementary_orgs = results['elementary_organizations']
    non_elementary_orgs = results['non_elementary_organizations']
    
    # Analyze size distributions
    if elementary_sos:
        eso_sizes = [len(eso.closure_species) for eso in elementary_sos]
        print(f"📊 Elementary SO size distribution:")
        print(f"   - Min species: {min(eso_sizes)}")
        print(f"   - Max species: {max(eso_sizes)}")
        print(f"   - Average species: {sum(eso_sizes)/len(eso_sizes):.1f}")
    
    # Analyze P-ERC vs multi-ERC
    if elementary_sos:
        p_erc_count = sum(1 for eso in elementary_sos if eso.is_p_erc)
        multi_erc_count = len(elementary_sos) - p_erc_count
        print(f"📊 Elementary SO composition:")
        print(f"   - P-ERCs: {p_erc_count}")
        print(f"   - Multi-ERC SOs: {multi_erc_count}")
    
    # Analyze organization vs semi-organization ratio
    if elementary_sos:
        org_ratio = len(elementary_orgs) / len(elementary_sos)
        print(f"📊 Organization analysis:")
        print(f"   - Self-maintaining ratio: {org_ratio:.1%}")
        print(f"   - Organizations: {len(elementary_orgs)}")
        print(f"   - Semi-organizations only: {len(elementary_sos) - len(elementary_orgs)}")
    
    # Timing analysis
    timing = results['timing']
    total_time = timing['total']
    print(f"📊 Performance analysis:")
    for phase, time_val in timing.items():
        if phase != 'total':
            percentage = (time_val / total_time) * 100
            print(f"   - {phase}: {time_val:.3f}s ({percentage:.1f}%)")

if __name__ == "__main__":
    try:
        results = test_persistent_modules_generator()
        
        if results:
            print("\n🎉 All tests passed!")
            analyze_results(results)
        else:
            print("\n⚠️  Tests completed but no meaningful results obtained.")
        
    except Exception as e:
        print(f"\n❌ Test failed with error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)