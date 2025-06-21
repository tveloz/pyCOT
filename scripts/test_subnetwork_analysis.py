#!/usr/bin/env python3
"""
Comprehensive comparison script for connectivity analysis and synergy algorithms.
Validates correctness and measures performance improvements.

ALGORITHM CONTROL:
- Set RUN_CONNECTIVITY = True/False to enable/disable connectivity analysis
- Set RUN_SYNERGY_ORIGINAL = True/False to enable/disable original synergy algorithm  
- Set RUN_SYNERGY_OPTIMIZED = True/False to enable/disable optimized synergy algorithm
- Set COMPARE_SYNERGY_ALGORITHMS = True/False to enable/disable synergy algorithm comparison
"""

import time
import os
import sys
import glob
import random
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict, deque

# Add the project root directory to the PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import pyCOT modules
from pyCOT.io.functions import read_txt, generate_subnetwork_txt
from pyCOT.rn_visualize import rn_visualize_html
from pyCOT.ERC_Hierarchy import ERC

# Try to import synergy algorithms
try:
    from pyCOT.Dynamical_Hierarchy_Synergy_Algorithm import EfficientFundamentalSynergyFinder
    DYNAMICAL_AVAILABLE = True
except ImportError as e:
    print(f"Warning: Dynamical synergy algorithm not available: {e}")
    DYNAMICAL_AVAILABLE = False

try:
    from pyCOT.Opt_Dynamical_Hierarchy_Synergy import OptimizedFundamentalSynergyFinder  
    OPTIMIZED_AVAILABLE = True
except ImportError as e:
    print(f"Warning: Optimized synergy algorithm not available: {e}")
    OPTIMIZED_AVAILABLE = False

# ALGORITHM CONTROL FLAGS
# ================================================================
RUN_CONNECTIVITY = True                    # Run connectivity analysis
RUN_SYNERGY_ORIGINAL = True and DYNAMICAL_AVAILABLE    # Run original synergy algorithm
RUN_SYNERGY_OPTIMIZED = True and OPTIMIZED_AVAILABLE   # Run optimized synergy algorithm
COMPARE_SYNERGY_ALGORITHMS = RUN_SYNERGY_ORIGINAL and RUN_SYNERGY_OPTIMIZED  # Compare algorithms

# Configuration
SAMPLING_CONFIG = {
    'num_reactions_to_sample': 10,      # Number of reactions to sample
    'num_samplings': 5,               # Number of different samplings to try
    'random_seed': 42,
    'min_component_size': 3,           # Minimum species in a component to analyze
    'save_subnetworks': False,          # Save generated subnetworks to files
}


def get_closure(rn, initial_species):
    """
    Get the closure of a set of species (all species reachable from them).
    
    Parameters:
    - rn: ReactionNetwork
    - initial_species: List of Species objects
    
    Returns:
    - Set of species names in the closure
    """
    if not initial_species:
        return set()
    
    # Start with initial species
    current_species = set(initial_species)
    closure = set(initial_species)
    
    # Keep expanding until no new species are found
    iteration = 0
    while True:
        # Get reactions that can be triggered by current species
        reactions = rn.get_reactions_from_species(list(current_species))
        
        # Filter to reactions where ALL reactants are available
        available_species_names = {sp.name for sp in current_species}
        valid_reactions = []
        
        for reaction in reactions:
            reactants = rn.get_supp_from_reactions([reaction])
            reactant_names = {sp.name for sp in reactants}
            
            # Check if all reactants are available
            if reactant_names.issubset(available_species_names):
                valid_reactions.append(reaction)
        
        # Get products from valid reactions
        if valid_reactions:
            products = rn.get_prod_from_reactions(valid_reactions)
            new_species = set(products) - closure
            
            if not new_species:
                break  # No new species found
            
            closure.update(new_species)
            current_species = closure
        else:
            break  # No valid reactions
        
        iteration += 1
    
    return {sp.name for sp in closure}


def find_connected_components_in_closure(rn, closure_species_names):
    """
    Find connected components within a closure.
    Two species are connected if they participate in the same reaction.
    
    Parameters:
    - rn: ReactionNetwork
    - closure_species_names: Set of species names in the closure
    
    Returns:
    - List of components, where each component is a set of species names
    """
    if not closure_species_names:
        return []
    
    # Build adjacency graph
    adjacency = defaultdict(set)
    
    # Get all species objects for the closure
    closure_species = [rn.get_species(name) for name in closure_species_names]
    
    # Get all reactions in the closure
    all_reactions = rn.get_reactions_from_species(closure_species)
    
    # Filter to reactions where ALL species are in the closure
    valid_reactions = []
    for reaction in all_reactions:
        reaction_species = rn.get_species_from_reactions([reaction])
        reaction_species_names = {sp.name for sp in reaction_species}
        
        if reaction_species_names.issubset(closure_species_names):
            valid_reactions.append(reaction)
    
    # Build adjacency from reactions
    for reaction in valid_reactions:
        # Get all species in this reaction
        reaction_species = rn.get_species_from_reactions([reaction])
        species_names = [sp.name for sp in reaction_species]
        
        # Connect all pairs of species in this reaction
        for i, sp1 in enumerate(species_names):
            for sp2 in species_names[i+1:]:
                adjacency[sp1].add(sp2)
                adjacency[sp2].add(sp1)
    
    # Find connected components using BFS
    visited = set()
    components = []
    
    for species_name in closure_species_names:
        if species_name in visited:
            continue
        
        # BFS from this species
        component = set()
        queue = deque([species_name])
        
        while queue:
            current = queue.popleft()
            if current in visited:
                continue
            
            visited.add(current)
            component.add(current)
            
            # Add neighbors
            for neighbor in adjacency[current]:
                if neighbor not in visited and neighbor in closure_species_names:
                    queue.append(neighbor)
        
        if component:
            components.append(component)
    
    # Handle isolated species (not in any reaction)
    for species_name in closure_species_names:
        if species_name not in visited:
            components.append({species_name})
    
    return components


def sample_and_analyze(rn, config, sample_id=1, verbose=True):
    """
    Sample reactions, compute closure, find components, and analyze.
    
    Returns:
    - Dictionary with analysis results
    """
    if verbose:
        print(f"\n{'='*60}")
        print(f"SAMPLING {sample_id}: {config['num_reactions_to_sample']} reactions")
        print(f"{'='*60}")
    
    # Sample random reactions
    all_reactions = rn.reactions()
    num_to_sample = min(config['num_reactions_to_sample'], len(all_reactions))
    sampled_reactions = random.sample(all_reactions, num_to_sample)
    sampled_reaction_names = [r.name() for r in sampled_reactions]
    
    if verbose:
        print(f"Sampled reactions: {sampled_reaction_names}")
    
    # Get initial species from sampled reactions
    initial_species = rn.get_species_from_reactions(sampled_reactions)
    initial_species_names = [sp.name for sp in initial_species]
    
    if verbose:
        print(f"Initial species ({len(initial_species_names)}): {initial_species_names}")
    
    # Compute closure
    start_time = time.time()
    closure_species_names = get_closure(rn, initial_species)
    closure_time = time.time() - start_time
    
    if verbose:
        print(f"Closure computed in {closure_time:.3f}s")
        print(f"Closure size: {len(closure_species_names)} species")
    
    # Find connected components within the closure
    components = find_connected_components_in_closure(rn, closure_species_names)
    
    # Sort components by size (largest first)
    components = sorted(components, key=len, reverse=True)
    
    if verbose:
        print(f"Found {len(components)} connected components:")
        for i, comp in enumerate(components):
            print(f"  Component {i+1}: {len(comp)} species")
    
    # Filter components by minimum size
    valid_components = [comp for comp in components if len(comp) >= config.get('min_component_size', 1)]
    
    if verbose and len(valid_components) < len(components):
        print(f"Filtered to {len(valid_components)} components with >= {config.get('min_component_size', 1)} species")
    
    # Get reactions for each component
    component_info = []
    for i, component in enumerate(valid_components):
        # Get species objects
        component_species = [rn.get_species(name) for name in component]
        
        # Get all reactions for these species
        all_component_reactions = rn.get_reactions_from_species(component_species)
        
        # Filter to reactions where ALL species are in the component
        valid_reactions = []
        for reaction in all_component_reactions:
            reaction_species = rn.get_species_from_reactions([reaction])
            reaction_species_names = {sp.name for sp in reaction_species}
            
            if reaction_species_names.issubset(component):
                valid_reactions.append(reaction)
        
        component_info.append({
            'id': f'sample{sample_id}_component{i+1}',
            'species': sorted(list(component)),
            'reactions': [r.name() for r in valid_reactions],
            'n_species': len(component),
            'n_reactions': len(valid_reactions),
            'from_reactions': sampled_reaction_names,
            'initial_species': initial_species_names
        })
    
    return {
        'sample_id': sample_id,
        'sampled_reactions': sampled_reaction_names,
        'initial_species': initial_species_names,
        'closure_size': len(closure_species_names),
        'num_components': len(components),
        'components': component_info,
        'closure_time': closure_time
    }


def run_synergy_analysis_on_components(rn, components, verbose=True):
    """Run synergy analysis on each component."""
    if not components:
        print("No components available for synergy analysis")
        return {}
    
    results = {}
    
    for i, comp_info in enumerate(components):
        comp_id = comp_info['id']
        
        if verbose:
            print(f"\n{'='*60}")
            print(f"SYNERGY ANALYSIS: {comp_id} ({i+1}/{len(components)})")
            print(f"{'='*60}")
            print(f"Component: {comp_info['n_species']} species, {comp_info['n_reactions']} reactions")
            print(f"From reactions: {comp_info['from_reactions']}")
        
        try:
            # Create sub-reaction network
            species_objects = [rn.get_species(name) for name in comp_info['species']]
            sub_rn = rn.sub_reaction_network(species_objects)
            
            # Verify the subnetwork
            sub_species_count = len(sub_rn.species())
            sub_reaction_count = len(sub_rn.reactions())
            
            if verbose:
                print(f"Created sub-reaction network: {sub_species_count} species, {sub_reaction_count} reactions")
            
            # Skip if too small
            if sub_species_count < 2 or sub_reaction_count < 1:
                if verbose:
                    print(f"Skipping - network too small for meaningful analysis")
                continue
            
            # Compute ERCs and hierarchy
            if verbose:
                print("Computing ERCs and hierarchy...")
            
            sub_ercs = ERC.ERCs(sub_rn)
            sub_hierarchy = ERC.build_hierarchy_graph(sub_ercs, sub_rn)
            
            if verbose:
                print(f"Component has {len(sub_ercs)} ERCs")
            
            if len(sub_ercs) == 0:
                if verbose:
                    print(f"Skipping - no ERCs found")
                continue
            
            # Run synergy algorithms
            synergy_results = {}
            
            # Original algorithm
            if RUN_SYNERGY_ORIGINAL:
                if verbose:
                    print(f"\nRunning ORIGINAL algorithm...")
                
                try:
                    original_finder = EfficientFundamentalSynergyFinder(sub_rn, sub_ercs, sub_hierarchy)
                    start_time = time.time()
                    original_results = original_finder.find_fundamental_synergies(verbose=False)
                    original_time = time.time() - start_time
                    
                    synergy_results['original'] = {
                        'synergies': original_results,
                        'time': original_time,
                        'count': len(original_results)
                    }
                    
                    if verbose:
                        print(f"  Original: {original_time:.3f}s, {len(original_results)} synergies")
                        
                except Exception as e:
                    print(f"  Error in original algorithm: {e}")
                    synergy_results['original'] = {'error': str(e)}
            
            # Optimized algorithm
            if RUN_SYNERGY_OPTIMIZED:
                if verbose:
                    print(f"Running OPTIMIZED algorithm...")
                
                try:
                    optimized_finder = OptimizedFundamentalSynergyFinder(sub_rn, sub_ercs, sub_hierarchy)
                    start_time = time.time()
                    optimized_results = optimized_finder.find_fundamental_synergies(verbose=False)
                    optimized_time = time.time() - start_time
                    
                    synergy_results['optimized'] = {
                        'synergies': optimized_results,
                        'time': optimized_time,
                        'count': len(optimized_results)
                    }
                    
                    if verbose:
                        print(f"  Optimized: {optimized_time:.3f}s, {len(optimized_results)} synergies")
                        
                except Exception as e:
                    print(f"  Error in optimized algorithm: {e}")
                    synergy_results['optimized'] = {'error': str(e)}
            
            # Compare results
            if (COMPARE_SYNERGY_ALGORITHMS and 
                'original' in synergy_results and 
                'optimized' in synergy_results and
                'error' not in synergy_results['original'] and 
                'error' not in synergy_results['optimized']):
                
                def synergy_to_tuple(syn):
                    return tuple(sorted([syn['base1'].label, syn['base2'].label]) + [syn['target'].label])
                
                original_set = set(synergy_to_tuple(s) for s in synergy_results['original']['synergies'])
                optimized_set = set(synergy_to_tuple(s) for s in synergy_results['optimized']['synergies'])
                
                results_match = original_set == optimized_set
                speedup = (synergy_results['original']['time'] / 
                          synergy_results['optimized']['time'] 
                          if synergy_results['optimized']['time'] > 0 else float('inf'))
                
                synergy_results['comparison'] = {
                    'results_match': results_match,
                    'speedup': speedup,
                    'original_only': original_set - optimized_set,
                    'optimized_only': optimized_set - original_set
                }
                
                if verbose:
                    print(f"\nComparison:")
                    print(f"  Results match: {'✓' if results_match else '✗'}")
                    if not results_match:
                        if synergy_results['comparison']['original_only']:
                            print(f"  Missing in optimized: {synergy_results['comparison']['original_only']}")
                        if synergy_results['comparison']['optimized_only']:
                            print(f"  Extra in optimized: {synergy_results['comparison']['optimized_only']}")
                    print(f"  Speedup: {speedup:.2f}x")
            
            results[comp_id] = {
                'info': comp_info,
                'n_ercs': len(sub_ercs),
                'synergy_results': synergy_results
            }
            
        except Exception as e:
            print(f"Error analyzing {comp_id}: {e}")
            results[comp_id] = {
                'info': comp_info,
                'error': str(e)
            }
    
    return results


def save_components_to_files(rn, sample_results, output_dir="generated_components"):
    """Save components to text files."""
    if not sample_results or not sample_results.get('components'):
        return
    
    os.makedirs(output_dir, exist_ok=True)
    
    saved_count = 0
    for comp_info in sample_results['components']:
        comp_id = comp_info['id']
        output_file = f"{comp_id}.txt"
        
        try:
            generate_subnetwork_txt(
                comp_info['species'], 
                comp_info['reactions'], 
                rn, 
                folder_name=output_dir,
                file_name=output_file
            )
            saved_count += 1
        except Exception as e:
            print(f"Error saving {comp_id}: {e}")
    
    if saved_count > 0:
        print(f"Saved {saved_count} components to {output_dir}/")


def analyze_file_with_sampling(file_path, config, verbose=True):
    """Analyze a file using multiple samplings."""
    if verbose:
        print(f"\n{'='*70}")
        print(f"ANALYZING: {os.path.basename(file_path)}")
        print(f"{'='*70}")
    
    # Set random seed
    if config.get('random_seed') is not None:
        random.seed(config['random_seed'])
    
    # Load reaction network
    try:
        rn = read_txt(file_path)
        print(f"Loaded network: {len(rn.species())} species, {len(rn.reactions())} reactions")
    except Exception as e:
        print(f"Error loading file: {e}")
        return None
    
    # Perform multiple samplings
    all_results = []
    all_components = []
    
    for sample_id in range(1, config['num_samplings'] + 1):
        sample_result = sample_and_analyze(rn, config, sample_id, verbose)
        all_results.append(sample_result)
        all_components.extend(sample_result['components'])
        
        # Save components if requested
        if config.get('save_subnetworks', False):
            save_components_to_files(rn, sample_result)
    
    # Remove duplicate components (by species set)
    unique_components = []
    seen_species_sets = set()
    
    for comp in all_components:
        species_tuple = tuple(sorted(comp['species']))
        if species_tuple not in seen_species_sets:
            seen_species_sets.add(species_tuple)
            unique_components.append(comp)
    
    if verbose:
        print(f"\n{'='*70}")
        print("SAMPLING SUMMARY")
        print(f"{'='*70}")
        print(f"Total samplings: {config['num_samplings']}")
        print(f"Total components found: {len(all_components)}")
        print(f"Unique components: {len(unique_components)}")
        
        # Component size distribution
        sizes = [comp['n_species'] for comp in unique_components]
        if sizes:
            print(f"Component sizes: min={min(sizes)}, max={max(sizes)}, avg={np.mean(sizes):.1f}")
    
    # Run synergy analysis on unique components
    if RUN_SYNERGY_ORIGINAL or RUN_SYNERGY_OPTIMIZED:
        if verbose:
            print(f"\n{'='*70}")
            print(f"RUNNING SYNERGY ANALYSIS ON {len(unique_components)} UNIQUE COMPONENTS")
            print(f"{'='*70}")
        
        synergy_results = run_synergy_analysis_on_components(rn, unique_components, verbose)
        
        # Print summary
        if synergy_results and verbose:
            print(f"\n{'='*70}")
            print("SYNERGY COMPARISON SUMMARY")
            print(f"{'='*70}")
            
            all_speedups = []
            all_matches = []
            
            for comp_id, result in synergy_results.items():
                if 'error' in result:
                    continue
                
                info = result['info']
                synergy = result.get('synergy_results', {})
                
                print(f"\n{comp_id}:")
                print(f"  Size: {info['n_species']} species, {info['n_reactions']} reactions, {result.get('n_ercs', 0)} ERCs")
                
                if 'comparison' in synergy:
                    comp = synergy['comparison']
                    print(f"  Match: {'✓' if comp['results_match'] else '✗'}")
                    print(f"  Speedup: {comp['speedup']:.2f}x")
                    
                    all_matches.append(comp['results_match'])
                    if comp['speedup'] != float('inf'):
                        all_speedups.append(comp['speedup'])
            
            # Overall statistics
            if all_matches:
                print(f"\n{'='*70}")
                print("OVERALL STATISTICS")
                print(f"{'='*70}")
                print(f"Components analyzed: {len(all_matches)}")
                print(f"Results matching: {sum(all_matches)}/{len(all_matches)} ({100*sum(all_matches)/len(all_matches):.1f}%)")
                
                if all_speedups:
                    print(f"Speedup statistics:")
                    print(f"  Average: {np.mean(all_speedups):.2f}x")
                    print(f"  Median: {np.median(all_speedups):.2f}x")
                    print(f"  Range: {min(all_speedups):.2f}x - {max(all_speedups):.2f}x")
    
    return {
        'file': os.path.basename(file_path),
        'network_size': {'species': len(rn.species()), 'reactions': len(rn.reactions())},
        'sampling_results': all_results,
        'unique_components': unique_components,
        'synergy_results': synergy_results if (RUN_SYNERGY_ORIGINAL or RUN_SYNERGY_OPTIMIZED) else None
    }


def main():
    """Main function."""
    print("="*70)
    print("REACTION NETWORK SAMPLING AND SYNERGY COMPARISON")
    print("="*70)
    print(f"Configuration:")
    print(f"  Reactions to sample: {SAMPLING_CONFIG['num_reactions_to_sample']}")
    print(f"  Number of samplings: {SAMPLING_CONFIG['num_samplings']}")
    print(f"  Minimum component size: {SAMPLING_CONFIG['min_component_size']}")
    print(f"  Save components: {'✓' if SAMPLING_CONFIG['save_subnetworks'] else '✗'}")
    print(f"\nAlgorithms:")
    print(f"  Original synergy: {'✓' if RUN_SYNERGY_ORIGINAL else '✗'}")
    print(f"  Optimized synergy: {'✓' if RUN_SYNERGY_OPTIMIZED else '✗'}")
    print(f"  Comparison: {'✓' if COMPARE_SYNERGY_ALGORITHMS else '✗'}")
    
    # Analyze file
    single_file = 'networks/biomodels_interesting/bigg_iAF692.txt'
    if os.path.exists(single_file):
        result = analyze_file_with_sampling(single_file, SAMPLING_CONFIG, verbose=True)
    else:
        print(f"\nFile {single_file} not found")
    
    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)


if __name__ == "__main__":
    main()