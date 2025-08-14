#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Complete ERC Synergy and Complementarity Analysis Script

Standalone comprehensive analysis script with specialized visualizations.
Designed for IDE execution with full pipeline from network loading to visualization.

Features:
- Complete synergy and complementarity analysis
- 4 specialized research-focused plots
- Additional comprehensive analysis plots
- Batch processing for multiple networks
- All-in-one pipeline for IDE execution

Usage (in IDE):
    - Set input_path to your network file or directory
    - Run the script directly
    
Author: Based on theoretical work by Tomas Veloz et al.
"""

import os
import sys
import time
import glob
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict, Counter
from itertools import combinations
from pathlib import Path

# Import required modules
from pyCOT.io.functions import read_txt
from pyCOT.ERC_Hierarchy import ERC, ERC_Hierarchy, closure

# Import core synergy/complementarity functions
from pyCOT.ERC_synergy import (
    ERC_Synergy, ERC_Complementarity,
    get_basic_synergies, get_maximal_synergies, get_fundamental_synergies,
    get_complementarity, can_interact
)

# ============================================================================
# CORE ANALYSIS FUNCTIONS
# ============================================================================

def analyze_all_synergies_and_complementarities(hierarchy, RN):
    """
    Analyze all synergies and complementarities in the ERC hierarchy.
    """
    print("Starting comprehensive synergy and complementarity analysis...")
    
    results = {
        'basic_synergies': [],
        'maximal_synergies': [],
        'fundamental_synergies': [],
        'complementarities_type1': [],
        'complementarities_type2': [],
        'complementarities_type3': [],
        'statistics': {}
    }
    
    total_pairs = 0
    interacting_pairs = 0
    
    print(f"Analyzing {len(hierarchy.ercs)} ERCs...")
    
    # Analyze all ERC pairs
    for i, erc1 in enumerate(hierarchy.ercs):
        if i % 10 == 0:  # Progress indicator
            print(f"  Progress: {i}/{len(hierarchy.ercs)} ERCs processed")
            
        for j, erc2 in enumerate(hierarchy.ercs):
            if i >= j:  # Avoid duplicates and self-comparison
                continue
                
            total_pairs += 1
            
            try:
                if can_interact(erc1, erc2, hierarchy):
                    interacting_pairs += 1
                    
                    # Get synergies
                    basic_syn = get_basic_synergies(erc1, erc2, hierarchy, RN)
                    maximal_syn = get_maximal_synergies(erc1, erc2, hierarchy, RN)
                    fundamental_syn = get_fundamental_synergies(erc1, erc2, hierarchy, RN)
                    
                    results['basic_synergies'].extend(basic_syn)
                    results['maximal_synergies'].extend(maximal_syn)
                    results['fundamental_synergies'].extend(fundamental_syn)
                    
                    # Get complementarities
                    complementarities = get_complementarity(erc1, erc2, hierarchy, RN)
                    
                    for comp in complementarities:
                        if comp.comp_type == 1:
                            results['complementarities_type1'].append(comp)
                        elif comp.comp_type == 2:
                            results['complementarities_type2'].append(comp)
                        elif comp.comp_type == 3:
                            results['complementarities_type3'].append(comp)
                            
            except Exception as e:
                print(f"⚠️  Warning: Error analyzing pair ({erc1.label}, {erc2.label}): {e}")
                continue
    
    # Calculate statistics
    results['statistics'] = {
        'total_erc_pairs': total_pairs,
        'interacting_pairs': interacting_pairs,
        'non_interacting_pairs': total_pairs - interacting_pairs,
        'basic_synergies_count': len(results['basic_synergies']),
        'maximal_synergies_count': len(results['maximal_synergies']),
        'fundamental_synergies_count': len(results['fundamental_synergies']),
        'complementarity_type1_count': len(results['complementarities_type1']),
        'complementarity_type2_count': len(results['complementarities_type2']),
        'complementarity_type3_count': len(results['complementarities_type3']),
        'synergy_efficiency': len(results['fundamental_synergies']) / len(results['basic_synergies']) if results['basic_synergies'] else 0,
        'interaction_efficiency': interacting_pairs / total_pairs if total_pairs > 0 else 0
    }
    
    print("Analysis completed!")
    return results

def prepare_participation_data(hierarchy, RN, results):
    """Calculate which ERCs participate in each type of synergy/complementarity."""
    # Initialize participation tracking for all ERCs
    all_erc_labels = [erc.label for erc in hierarchy.ercs]
    erc_participation = {
        'basic_synergy': {label: 0 for label in all_erc_labels},
        'maximal_synergy': {label: 0 for label in all_erc_labels}, 
        'fundamental_synergy': {label: 0 for label in all_erc_labels},
        'type1_complementarity': {label: 0 for label in all_erc_labels},
        'type2_complementarity': {label: 0 for label in all_erc_labels},
        'type3_complementarity': {label: 0 for label in all_erc_labels}
    }
    
    # Count synergy participations
    for syn in results['basic_synergies']:
        for erc in syn.reactants:
            if erc.label in erc_participation['basic_synergy']:
                erc_participation['basic_synergy'][erc.label] += 1
        if syn.product.label in erc_participation['basic_synergy']:
            erc_participation['basic_synergy'][syn.product.label] += 1
    
    for syn in results['maximal_synergies']:
        for erc in syn.reactants:
            if erc.label in erc_participation['maximal_synergy']:
                erc_participation['maximal_synergy'][erc.label] += 1
        if syn.product.label in erc_participation['maximal_synergy']:
            erc_participation['maximal_synergy'][syn.product.label] += 1
            
    for syn in results['fundamental_synergies']:
        for erc in syn.reactants:
            if erc.label in erc_participation['fundamental_synergy']:
                erc_participation['fundamental_synergy'][erc.label] += 1
        if syn.product.label in erc_participation['fundamental_synergy']:
            erc_participation['fundamental_synergy'][syn.product.label] += 1
    
    # Count complementarity participations
    for comp in results['complementarities_type1']:
        if comp.erc1.label in erc_participation['type1_complementarity']:
            erc_participation['type1_complementarity'][comp.erc1.label] += 1
        if comp.erc2.label in erc_participation['type1_complementarity']:
            erc_participation['type1_complementarity'][comp.erc2.label] += 1
        
    for comp in results['complementarities_type2']:
        if comp.erc1.label in erc_participation['type2_complementarity']:
            erc_participation['type2_complementarity'][comp.erc1.label] += 1
        if comp.erc2.label in erc_participation['type2_complementarity']:
            erc_participation['type2_complementarity'][comp.erc2.label] += 1
        
    for comp in results['complementarities_type3']:
        if comp.erc1.label in erc_participation['type3_complementarity']:
            erc_participation['type3_complementarity'][comp.erc1.label] += 1
        if comp.erc2.label in erc_participation['type3_complementarity']:
            erc_participation['type3_complementarity'][comp.erc2.label] += 1
    
    # Convert to lists of participation counts (preserving order)
    participation_data = {}
    for interaction_type, erc_counts in erc_participation.items():
        participation_counts = [erc_counts[label] for label in all_erc_labels]
        participation_data[interaction_type] = participation_counts
    
    return participation_data

def prepare_hierarchy_level_data(hierarchy, RN, results):
    """Calculate synergy distribution across hierarchy levels."""
    if not hierarchy.graph or len(hierarchy.graph.nodes()) == 0:
        return {'level_pair_counts': {}, 'level_pair_proportions': {}, 'total_synergies': 0, 'max_level': 0}
    
    # Get ERC hierarchy levels using topological sorting
    try:
        import networkx as nx
        if hierarchy.graph.is_directed():
            dag = hierarchy.graph
        else:
            dag = nx.DiGraph()
            dag.add_edges_from(hierarchy.graph.edges())
        
        # Calculate levels for each ERC
        erc_levels = {}
        if nx.is_directed_acyclic_graph(dag) and len(dag.nodes()) > 0:
            try:
                levels = list(nx.topological_generations(dag))
                for level_idx, level_nodes in enumerate(levels):
                    for node in level_nodes:
                        erc_levels[node] = level_idx
            except:
                # Fallback: assign levels based on in-degree
                for node in dag.nodes():
                    erc_levels[node] = dag.in_degree(node)
        else:
            # Fallback: assign all to level 0
            for erc in hierarchy.ercs:
                erc_levels[erc.label] = 0
                
    except Exception as e:
        print(f"Warning: Could not calculate hierarchy levels: {e}")
        # Fallback: assign all ERCs to level 0
        erc_levels = {erc.label: 0 for erc in hierarchy.ercs}
    
    # Count synergies by level pairs
    level_pair_counts = defaultdict(int)
    total_synergies = 0
    
    # Process fundamental synergies (most important)
    for syn in results['fundamental_synergies']:
        if len(syn.reactants) >= 2:
            level1 = erc_levels.get(syn.reactants[0].label, 0)
            level2 = erc_levels.get(syn.reactants[1].label, 0)
            # Create sorted pair to avoid (0,1) vs (1,0) duplication
            level_pair = tuple(sorted([level1, level2]))
            level_pair_counts[level_pair] += 1
            total_synergies += 1
    
    # Convert to proportions
    level_pair_proportions = {}
    for level_pair, count in level_pair_counts.items():
        level_pair_proportions[level_pair] = count / total_synergies if total_synergies > 0 else 0
    
    return {
        'level_pair_counts': {str(k): v for k, v in level_pair_counts.items()},  # Convert tuples to strings for JSON compatibility
        'level_pair_proportions': {str(k): v for k, v in level_pair_proportions.items()},
        'total_synergies': total_synergies,
        'max_level': max(erc_levels.values()) if erc_levels else 0
    }

def analyze_basic_hierarchy_structure(graph, n_ercs):
    """Basic hierarchy structure analysis (fallback when enhanced benchmark unavailable)."""
    n_nodes = len(graph.nodes()) if graph else 0
    n_edges = len(graph.edges()) if graph else 0
    
    return {
        'n_edges': n_edges,
        'hierarchy_density': n_edges / (n_nodes * (n_nodes - 1) / 2) if n_nodes > 1 else 0,
        'avg_degree': 2 * n_edges / n_nodes if n_nodes > 0 else 0,
        'n_levels': 1,  # Default
        'avg_level_size': n_nodes,
        'max_level_size': n_nodes,
        'level_variance': 0,
        'leaf_nodes': 0,
        'root_nodes': 0,
        'leaf_ratio': 0,
        'root_ratio': 0,
        'max_chain_length': 0,
        'hierarchy_depth_ratio': 0,
        'avg_branching': 0,
        'max_branching': 0,
        'n_cliques': 1,
        'avg_closure_size': 0,
        'ercs_per_species': 0
    }

def analyze_single_network_complete(network_file):
    """Complete analysis of a single network including all data preparation."""
    print(f"\n{'='*80}")
    print(f"COMPLETE ANALYSIS: {os.path.basename(network_file)}")
    print(f"{'='*80}")
    
    try:
        # Load network
        print("Loading reaction network...")
        RN = read_txt(network_file)
        print(f"Network loaded: {len(RN.species())} species, {len(RN.reactions())} reactions")
        
        # Create hierarchy using optimized ERC_Hierarchy class
        print("Creating ERC hierarchy...")
        start_time = time.time()
        hierarchy = ERC_Hierarchy(RN)
        hierarchy.build_hierarchy_graph3()
        hierarchy_time = time.time() - start_time
        
        # Analyze hierarchy structure
        try:
            from enhanced_hierarchy_benchmark import analyze_hierarchy_structure
            hierarchy_props = analyze_hierarchy_structure(hierarchy.graph, hierarchy.ercs, RN)
        except ImportError:
            print("⚠️  Warning: Using basic structure analysis")
            hierarchy_props = analyze_basic_hierarchy_structure(hierarchy.graph, len(hierarchy.ercs))
        
        print(f"Hierarchy created: {len(hierarchy.ercs)} ERCs in {hierarchy_time:.4f}s")
        
        # Run synergy and complementarity analysis
        print("Running synergy and complementarity analysis...")
        start_time = time.time()
        results = analyze_all_synergies_and_complementarities(hierarchy, RN)
        analysis_time = time.time() - start_time
        print(f"Analysis completed in {analysis_time:.4f}s")
        
        # Prepare detailed data for plotting
        participation_data = prepare_participation_data(hierarchy, RN, results)
        hierarchy_level_data = prepare_hierarchy_level_data(hierarchy, RN, results)
        
        # Compile comprehensive result
        stats = results['statistics']
        network_data = {
            'network_file': os.path.basename(network_file),
            'n_species': len(RN.species()),
            'n_reactions': len(RN.reactions()),
            'n_ercs': len(hierarchy.ercs),
            'hierarchy_time': hierarchy_time,
            'analysis_time': analysis_time,
            'total_time': hierarchy_time + analysis_time,
            'basic_synergies': stats['basic_synergies_count'],
            'maximal_synergies': stats['maximal_synergies_count'],
            'fundamental_synergies': stats['fundamental_synergies_count'],
            'type1_complementarity': stats['complementarity_type1_count'],
            'type2_complementarity': stats['complementarity_type2_count'],
            'type3_complementarity': stats['complementarity_type3_count'],
            'total_complementarity': (stats['complementarity_type1_count'] + 
                                    stats['complementarity_type2_count'] + 
                                    stats['complementarity_type3_count']),
            'synergy_efficiency': stats['synergy_efficiency'],
            'interaction_efficiency': stats['interaction_efficiency'],
            'total_productive': (stats['fundamental_synergies_count'] + 
                               stats['complementarity_type1_count'] + 
                               stats['complementarity_type2_count'] + 
                               stats['complementarity_type3_count']),
            'productivity_rate': ((stats['fundamental_synergies_count'] + 
                                 stats['complementarity_type1_count'] + 
                                 stats['complementarity_type2_count'] + 
                                 stats['complementarity_type3_count']) / 
                                 max(stats['interacting_pairs'], 1)),
            'synergy_density': stats['fundamental_synergies_count'] / len(hierarchy.ercs) if len(hierarchy.ercs) > 0 else 0,
            'complementarity_density': (stats['complementarity_type1_count'] + 
                                      stats['complementarity_type2_count'] + 
                                      stats['complementarity_type3_count']) / len(hierarchy.ercs) if len(hierarchy.ercs) > 0 else 0,
            'participation_data': participation_data,
            'hierarchy_level_data': hierarchy_level_data
        }
        
        # Add hierarchy properties
        network_data.update(hierarchy_props)
        
        return network_data, results
        
    except Exception as e:
        print(f"❌ Error analyzing {network_file}: {e}")
        import traceback
        traceback.print_exc()
        return None, None

# ============================================================================
# SPECIALIZED PLOTTING FUNCTIONS
# ============================================================================

def create_synergy_scaling_plot(df, ax):
    """Plot 1: x-axis number ERCs, y-axis 3 types of synergy, 3 curves"""
    df_sorted = df.sort_values('n_ercs')
    
    # Plot three synergy curves
    ax.plot(df_sorted['n_ercs'], df_sorted['basic_synergies'], 'o-', 
           label='Basic Synergies', alpha=0.7, color='lightcoral', linewidth=2)
    ax.plot(df_sorted['n_ercs'], df_sorted['maximal_synergies'], 's-', 
           label='Maximal Synergies', alpha=0.7, color='gold', linewidth=2)
    ax.plot(df_sorted['n_ercs'], df_sorted['fundamental_synergies'], '^-', 
           label='Fundamental Synergies', alpha=0.7, color='lightgreen', linewidth=2)
    
    ax.set_xlabel('Number of ERCs')
    ax.set_ylabel('Number of Synergies')
    ax.set_title('Synergy Types vs Network Size')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_xscale('log')
    ax.set_yscale('log')

def create_participation_distribution_plot(df, detailed_data, ax):
    """Plot 2: ERC participation distribution across interaction types (6 subplots)"""
    ax.clear()
    ax.axis('off')
    
    participation_types = [
        ('basic_synergy', 'Basic Synergy'),
        ('maximal_synergy', 'Maximal Synergy'), 
        ('fundamental_synergy', 'Fundamental Synergy'),
        ('type1_complementarity', 'Type 1 Comp.'),
        ('type2_complementarity', 'Type 2 Comp.'),
        ('type3_complementarity', 'Type 3 Comp.')
    ]
    
    # Create 2x3 subplots within the axis
    fig = plt.gcf()
    
    for i, (interaction_type, title) in enumerate(participation_types):
        # Position calculation for 2x3 grid
        row = i // 3
        col = i % 3
        
        # Create subplot position [left, bottom, width, height]
        left = 0.1 + col * 0.28
        bottom = 0.55 - row * 0.4
        width = 0.25
        height = 0.35
        
        ax_sub = fig.add_axes([left, bottom, width, height])
        
        # Collect participation counts across all networks
        all_participation_counts = []
        
        for data in detailed_data:
            participation_data = data.get('participation_data', {})
            counts = participation_data.get(interaction_type, [])
            
            if counts:  # Only if we have data
                all_participation_counts.extend(counts)
        
        if all_participation_counts:
            # Create histogram of participation counts
            max_participation = max(all_participation_counts)
            bins = list(range(0, min(max_participation + 2, 11)))  # Cap at 10+ for readability
            
            # Create histogram
            counts, bins_edges, patches = ax_sub.hist(all_participation_counts, bins=bins, 
                                                    alpha=0.7, edgecolor='black', color='skyblue')
            
            # Calculate and display statistics
            total_ercs = len(all_participation_counts)
            participating_ercs = sum(1 for c in all_participation_counts if c > 0)
            participation_rate = participating_ercs / total_ercs if total_ercs > 0 else 0
            avg_participation = np.mean(all_participation_counts)
            
            # Add statistics text box
            stats_text = f'Total: {total_ercs}\nActive: {participating_ercs}\nRate: {participation_rate:.2f}\nAvg: {avg_participation:.1f}'
            ax_sub.text(0.98, 0.98, stats_text,
                       transform=ax_sub.transAxes, verticalalignment='top', horizontalalignment='right',
                       bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
                       fontsize=8)
        else:
            # No data case
            ax_sub.text(0.5, 0.5, 'No Data', ha='center', va='center', 
                       transform=ax_sub.transAxes, fontsize=12, color='gray')
        
        # Format subplot
        ax_sub.set_title(title, fontsize=10, fontweight='bold')
        ax_sub.set_xlabel('Participation Count', fontsize=8)
        ax_sub.set_ylabel('# ERCs', fontsize=8)
        ax_sub.grid(True, alpha=0.3)
        ax_sub.tick_params(labelsize=8)
    
    # Add main title to the entire plot
    ax.text(0.5, 0.95, 'ERC Participation Distribution by Interaction Type', 
           ha='center', va='top', transform=ax.transAxes, fontsize=14, fontweight='bold')
    ax.text(0.5, 0.05, 'Hub ERCs (high participation) vs Peripheral ERCs (low participation)', 
           ha='center', va='bottom', transform=ax.transAxes, fontsize=10, style='italic')

def create_complementarity_scaling_plot(df, ax):
    """Plot 3: x-axis number ERCs, y-axis 3 types of complementarity, 3 curves"""
    df_sorted = df.sort_values('n_ercs')
    
    # Plot three complementarity curves
    ax.plot(df_sorted['n_ercs'], df_sorted['type1_complementarity'], 'o-', 
           label='Type 1 (Requirement Reduction)', alpha=0.7, color='skyblue', linewidth=2)
    ax.plot(df_sorted['n_ercs'], df_sorted['type2_complementarity'], 's-', 
           label='Type 2 (Requirement Change)', alpha=0.7, color='orange', linewidth=2)
    ax.plot(df_sorted['n_ercs'], df_sorted['type3_complementarity'], '^-', 
           label='Type 3 (Product Expansion)', alpha=0.7, color='lightpink', linewidth=2)
    
    ax.set_xlabel('Number of ERCs')
    ax.set_ylabel('Number of Complementarities')
    ax.set_title('Complementarity Types vs Network Size')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_xscale('log')
    ax.set_yscale('log')

def create_hierarchy_level_synergy_plot(df, detailed_data, ax):
    """Plot 4: x-axis ERC size quartiles, y-axis all level pairs, color: synergy proportion"""
    
    # Group networks by ERC size quartiles (handle edge cases)
    try:
        df['size_quartile'] = pd.qcut(df['n_ercs'], q=4, labels=['Q1', 'Q2', 'Q3', 'Q4'], duplicates='drop')
    except ValueError:
        # If not enough unique values for quartiles, use simple binning
        df['size_quartile'] = pd.cut(df['n_ercs'], bins=4, labels=['Q1', 'Q2', 'Q3', 'Q4'])
    
    # Handle case where quartile creation fails
    if df['size_quartile'].isna().all():
        # Fallback: assign all to Q1
        df['size_quartile'] = 'Q1'
    
    # Determine all possible level pairs from the data
    all_level_pairs = set()
    max_level = 0
    
    for data in detailed_data:
        hierarchy_data = data.get('hierarchy_level_data', {})
        level_pairs = hierarchy_data.get('level_pair_counts', {})
        for level_pair_str in level_pairs.keys():
            # Parse level pair string like "(0, 1)" back to tuple
            try:
                level_pair = eval(level_pair_str) if isinstance(level_pair_str, str) else level_pair_str
                if isinstance(level_pair, tuple) and len(level_pair) == 2:
                    all_level_pairs.add(level_pair)
                    max_level = max(max_level, max(level_pair))
            except:
                continue
    
    # If no level pairs found, create default ones
    if not all_level_pairs:
        max_level = 3  # Default assumption
        for i in range(max_level + 1):
            for j in range(i, max_level + 1):
                all_level_pairs.add((i, j))
    
    # Convert to sorted list for consistent ordering
    level_pairs = sorted(list(all_level_pairs))
    
    # Get actual quartiles present in the data
    actual_quartiles = [q for q in ['Q1', 'Q2', 'Q3', 'Q4'] if q in df['size_quartile'].values]
    if not actual_quartiles:
        actual_quartiles = ['Q1']  # Fallback
    
    # Initialize synergy proportion matrix
    synergy_matrix = np.zeros((len(actual_quartiles), len(level_pairs)))
    
    # Calculate actual synergy proportions by quartile and level pair
    for q_idx, quartile in enumerate(actual_quartiles):
        quartile_networks = df[df['size_quartile'] == quartile]
        quartile_detailed = [data for data in detailed_data 
                           if data['network_file'] in quartile_networks['network_file'].values]
        
        # Aggregate synergy counts by level pair for this quartile
        total_synergies_in_quartile = 0
        level_pair_synergies = defaultdict(int)
        
        for data in quartile_detailed:
            hierarchy_data = data.get('hierarchy_level_data', {})
            level_pair_counts = hierarchy_data.get('level_pair_counts', {})
            total_synergies = hierarchy_data.get('total_synergies', 0)
            
            total_synergies_in_quartile += total_synergies
            
            for level_pair_str, count in level_pair_counts.items():
                try:
                    level_pair = eval(level_pair_str) if isinstance(level_pair_str, str) else level_pair_str
                    if isinstance(level_pair, tuple) and len(level_pair) == 2:
                        level_pair_synergies[level_pair] += count
                except:
                    continue
        
        # Calculate proportions for this quartile
        for lp_idx, level_pair in enumerate(level_pairs):
            if total_synergies_in_quartile > 0:
                proportion = level_pair_synergies[level_pair] / total_synergies_in_quartile
                synergy_matrix[q_idx, lp_idx] = proportion
    
    # Create heatmap
    im = ax.imshow(synergy_matrix, cmap='Reds', aspect='auto', vmin=0, vmax=1)
    
    # Set labels
    ax.set_xticks(range(len(level_pairs)))
    ax.set_xticklabels([f'({l1},{l2})' for l1, l2 in level_pairs], rotation=45, fontsize=10)
    ax.set_yticks(range(len(actual_quartiles)))
    ax.set_yticklabels(actual_quartiles)
    
    ax.set_xlabel('Hierarchy Level Pairs')
    ax.set_ylabel('ERC Size Quartiles')
    ax.set_title('Synergy Proportion by Hierarchy Levels\n(White = 0, Red = 1)')
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Proportion of Total Synergies')
    
    # Add text annotations for values > 0.01 (avoid cluttering)
    for i in range(len(actual_quartiles)):
        for j in range(len(level_pairs)):
            if synergy_matrix[i, j] > 0.01:
                text = ax.text(j, i, f'{synergy_matrix[i, j]:.2f}',
                             ha="center", va="center", 
                             color="white" if synergy_matrix[i, j] > 0.5 else "black", 
                             fontsize=8)

def create_additional_analysis_plots(df, detailed_data, output_dir):
    """Additional comprehensive analysis plots."""
    plt.style.use('seaborn-v0_8' if 'seaborn-v0_8' in plt.style.available else 'default')
    
    fig = plt.figure(figsize=(20, 15))
    gs = fig.add_gridspec(3, 4, hspace=0.3, wspace=0.3)
    
    # 1. 3D plot: ERCs vs Synergies vs Complementarities
    ax1 = fig.add_subplot(gs[0, 0:2], projection='3d')
    ax1.scatter(df['n_ercs'], df['fundamental_synergies'], df['total_complementarity'],
               c=df['productivity_rate'], cmap='plasma', alpha=0.7)
    ax1.set_xlabel('Number of ERCs')
    ax1.set_ylabel('Fundamental Synergies')
    ax1.set_zlabel('Total Complementarities')
    ax1.set_title('3D Generative Capacity Analysis')
    
    # 2. Correlation heatmap
    ax2 = fig.add_subplot(gs[0, 2:4])
    correlation_cols = ['n_ercs', 'fundamental_synergies', 'total_complementarity', 
                       'productivity_rate', 'synergy_efficiency', 'interaction_efficiency']
    correlation_cols = [col for col in correlation_cols if col in df.columns]
    
    if len(correlation_cols) > 1:
        corr_matrix = df[correlation_cols].corr()
        sns.heatmap(corr_matrix, annot=True, cmap='coolwarm', center=0, ax=ax2)
        ax2.set_title('Correlation Matrix: Key Metrics')
    
    # 3. Productivity distribution pie chart
    ax3 = fig.add_subplot(gs[1, :2])
    df['efficiency_category'] = pd.cut(df['productivity_rate'], 
                                     bins=[0, 0.1, 0.3, 0.5, 1.0], 
                                     labels=['Low', 'Medium', 'High', 'Very High'])
    
    efficiency_counts = df['efficiency_category'].value_counts()
    colors = ['red', 'orange', 'lightgreen', 'darkgreen']
    wedges, texts, autotexts = ax3.pie(efficiency_counts.values, labels=efficiency_counts.index, 
                                       autopct='%1.1f%%', colors=colors)
    ax3.set_title('Network Productivity Distribution')
    
    # 4. Synergy efficiency vs complementarity density
    ax4 = fig.add_subplot(gs[1, 2:])
    ax4.scatter(df['synergy_efficiency'], df['complementarity_density'], 
               alpha=0.6, c=df['n_ercs'], cmap='viridis', s=60)
    ax4.set_xlabel('Synergy Efficiency (Fundamental/Basic)')
    ax4.set_ylabel('Complementarity Density')
    ax4.set_title('Synergy Efficiency vs Complementarity Density')
    
    scatter = ax4.collections[0]
    cbar = plt.colorbar(scatter, ax=ax4)
    cbar.set_label('Number of ERCs')
    
    # 5. Time complexity analysis
    ax5 = fig.add_subplot(gs[2, :2])
    valid_data = df[df['total_time'].notna()]
    if len(valid_data) > 0:
        ax5.scatter(valid_data['n_ercs'], valid_data['total_time'], 
                   c=valid_data['productivity_rate'], cmap='viridis', alpha=0.7)
        ax5.set_xlabel('Number of ERCs')
        ax5.set_ylabel('Computation Time (s)')
        ax5.set_title('Computational Complexity vs Network Size')
        ax5.set_xscale('log')
        ax5.set_yscale('log')
        
        scatter = ax5.collections[0]
        cbar = plt.colorbar(scatter, ax=ax5)
        cbar.set_label('Productivity Rate')
    
    # 6. Summary statistics table
    ax6 = fig.add_subplot(gs[2, 2:])
    ax6.axis('tight')
    ax6.axis('off')
    
    summary_stats = {
        'Metric': ['Networks Analyzed', 'Avg ERCs', 'Avg Synergies', 'Avg Complementarities', 
                   'Avg Productivity', 'Avg Efficiency', 'Avg Computation Time'],
        'Value': [len(df), 
                 f"{df['n_ercs'].mean():.1f}",
                 f"{df['fundamental_synergies'].mean():.1f}",
                 f"{df['total_complementarity'].mean():.1f}",
                 f"{df['productivity_rate'].mean():.3f}",
                 f"{df['synergy_efficiency'].mean():.3f}",
                 f"{df['total_time'].mean():.3f}s" if 'total_time' in df.columns else "N/A"]
    }
    
    table = ax6.table(cellText=list(zip(summary_stats['Metric'], summary_stats['Value'])),
                     colLabels=['Metric', 'Value'],
                     cellLoc='center',
                     loc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 1.5)
    ax6.set_title('Summary Statistics', pad=20)
    
    plt.suptitle('Additional Comprehensive Analysis', fontsize=16, fontweight='bold')
    
    plot_file = os.path.join(output_dir, "additional_analysis.png")
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"📊 Additional analysis plots saved to: {plot_file}")
    plt.show()

def create_all_plots(df, detailed_data, output_dir):
    """Create all visualization plots."""
    os.makedirs(output_dir, exist_ok=True)
    
    # Main specialized plots (4 requested plots)
    print("Creating specialized research plots...")
    plt.style.use('seaborn-v0_8' if 'seaborn-v0_8' in plt.style.available else 'default')
    
    fig = plt.figure(figsize=(20, 15))
    gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.25)
    
    # Plot 1: Synergy Types vs Network Size
    ax1 = fig.add_subplot(gs[0, 0])
    create_synergy_scaling_plot(df, ax1)
    
    # Plot 2: ERC Participation Distribution  
    ax2 = fig.add_subplot(gs[0, 1])
    create_participation_distribution_plot(df, detailed_data, ax2)
    
    # Plot 3: Complementarity Types vs Network Size
    ax3 = fig.add_subplot(gs[1, 0])
    create_complementarity_scaling_plot(df, ax3)
    
    # Plot 4: Hierarchy Level Synergy Distribution
    ax4 = fig.add_subplot(gs[1, 1])
    create_hierarchy_level_synergy_plot(df, detailed_data, ax4)
    
    plt.suptitle('Specialized Synergy and Complementarity Analysis', fontsize=16, fontweight='bold')
    
    plot_file = os.path.join(output_dir, "specialized_synergy_analysis.png")
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"📊 Specialized synergy plots saved to: {plot_file}")
    plt.show()
    
    # Additional comprehensive plots
    create_additional_analysis_plots(df, detailed_data, output_dir)

# ============================================================================
# BATCH ANALYSIS FUNCTIONS
# ============================================================================

def analyze_directory_complete(directory_path, output_dir="complete_synergy_results"):
    """Complete analysis of all networks in a directory."""
    print(f"\n{'='*80}")
    print(f"COMPLETE BATCH ANALYSIS: {directory_path}")
    print(f"{'='*80}")
    
    # Find network files
    network_files = []
    for ext in ['*.txt', '*.net', '*.rn']:
        network_files.extend(glob.glob(os.path.join(directory_path, ext)))
    
    if not network_files:
        print(f"❌ No network files found in {directory_path}")
        return None
    
    print(f"Found {len(network_files)} network files")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Analyze each network
    summary_data = []
    detailed_data = []
    successful_analyses = 0
    
    for i, network_file in enumerate(network_files, 1):
        print(f"\n{'='*20} NETWORK {i}/{len(network_files)} {'='*20}")
        
        network_data, results = analyze_single_network_complete(network_file)
        
        if network_data:
            successful_analyses += 1
            summary_data.append({k: v for k, v in network_data.items() 
                               if k not in ['participation_data', 'hierarchy_level_data']})
            detailed_data.append(network_data)
    
    if not summary_data:
        print("❌ No successful analyses to process.")
        return None
    
    # Create DataFrame
    df = pd.DataFrame(summary_data)
    
    # Save detailed results
    summary_file = os.path.join(output_dir, "complete_synergy_analysis.csv")
    df.to_csv(summary_file, index=False)
    
    # Create all visualizations
    create_all_plots(df, detailed_data, output_dir)
    
    # Print summary
    print_complete_summary(df)
    
    print(f"\n💾 Complete results saved to: {summary_file}")
    print(f"📊 All visualizations saved in: {output_dir}")
    
    return df, detailed_data

def print_complete_summary(df):
    """Print complete analysis summary."""
    print("\n" + "="*80)
    print("COMPLETE SYNERGY AND COMPLEMENTARITY ANALYSIS SUMMARY")
    print("="*80)
    
    print(f"\nDataset Overview:")
    print(f"  Networks analyzed: {len(df)}")
    print(f"  ERC size range: {df['n_ercs'].min()}-{df['n_ercs'].max()}")
    print(f"  Total ERCs: {df['n_ercs'].sum()}")
    print(f"  Total fundamental synergies: {df['fundamental_synergies'].sum()}")
    print(f"  Total complementarities: {df['total_complementarity'].sum()}")
    
    print(f"\nGenerative Capacity Statistics:")
    print(f"  Average productivity rate: {df['productivity_rate'].mean():.3f}")
    print(f"  Average synergy efficiency: {df['synergy_efficiency'].mean():.3f}")
    print(f"  Average interaction efficiency: {df['interaction_efficiency'].mean():.3f}")
    print(f"  Networks with high productivity (>0.3): {len(df[df['productivity_rate'] > 0.3])}")
    
    print(f"\nScaling Relationships:")
    synergy_corr = df[['n_ercs', 'fundamental_synergies']].corr().iloc[0, 1]
    comp_corr = df[['n_ercs', 'total_complementarity']].corr().iloc[0, 1]
    print(f"  ERCs vs Fundamental Synergies: r = {synergy_corr:.3f}")
    print(f"  ERCs vs Total Complementarities: r = {comp_corr:.3f}")
    
    print(f"\nComputational Performance:")
    if 'total_time' in df.columns:
        print(f"  Average analysis time: {df['total_time'].mean():.3f}s")
        print(f"  Time per ERC: {(df['total_time'] / df['n_ercs']).mean():.4f}s")

# ============================================================================
# MAIN EXECUTION FOR IDE
# ============================================================================

def run_complete_analysis():
    """
    Main function for complete synergy analysis - designed for IDE execution.
    Modify the input_path variable below to point to your network(s).
    """
    
    # ===== CONFIGURATION SECTION - MODIFY THESE PATHS =====
    
    # Option 1: Single network file
    # input_path = "path/to/your/network.txt"
    
    # Option 2: Directory with multiple networks
    input_path = "C:/Users/tvelo/Dropbox/Public/AcademicWork/Europe/CLEA/Postdocs/TempletonPostdoc/sftw/pyCOT/networks/RandomAlife"
    
    # Output directory for results and plots
    output_dir = "complete_synergy_results_random_alife"
    
    # ===== END CONFIGURATION SECTION =====
    
    # Validate input path
    input_path = os.path.abspath(input_path)
    if input_path==None:
        print(f"❌ Path not found: {input_path}")
        return None
    
    print(f"🚀 Starting Complete Synergy Analysis")
    print(f"📁 Input: {input_path}")
    print(f"📊 Output: {output_dir}")
    
    try:
        if os.path.isfile(input_path):
            # Single file analysis
            print(f"\n🔬 Analyzing single network...")
            network_data, results = analyze_single_network_complete(input_path)
            
            if network_data:
                print(f"\n✅ Single network analysis completed successfully!")
                
                # Create simple summary for single network
                print(f"\nSingle Network Summary:")
                print(f"  ERCs: {network_data['n_ercs']}")
                print(f"  Fundamental Synergies: {network_data['fundamental_synergies']}")
                print(f"  Total Complementarities: {network_data['total_complementarity']}")
                print(f"  Productivity Rate: {network_data['productivity_rate']:.3f}")
                print(f"  Analysis Time: {network_data['total_time']:.3f}s")
                
                return network_data, results
            else:
                print(f"\n❌ Single network analysis failed!")
                return None
                
        elif os.path.isdir(input_path):
            # Directory analysis with complete visualization
            print(f"\n📊 Analyzing directory of networks...")
            result = analyze_directory_complete(input_path, output_dir)
            
            if result is not None:
                df, detailed_data = result
                print(f"\n✅ Complete batch analysis successful!")
                print(f"📈 Analyzed {len(df)} networks with specialized visualizations")
                return df, detailed_data
            else:
                print(f"\n❌ Batch analysis failed!")
                return None
        else:
            print(f"❌ Invalid path type: {input_path}")
            return None
            
    except Exception as e:
        print(f"❌ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        return None

# ============================================================================
# IDE EXECUTION
# ============================================================================

if __name__ == "__main__":
    # Run the complete analysis
    results = run_complete_analysis()
    
    if results:
        print(f"\n🎉 Analysis completed successfully!")
        print(f"📊 Check the output directory for detailed results and visualizations")
        print(f"🔬 The following plots were generated:")
        print(f"   1. Synergy scaling with network size (3 curves)")
        print(f"   2. ERC participation distributions (6 subplots)")
        print(f"   3. Complementarity scaling with network size (3 curves)")
        print(f"   4. Hierarchy level synergy patterns (quartiles heatmap)")
        print(f"   5. Additional comprehensive analysis plots")
    else:
        print(f"\n💥 Analysis failed - check the configuration and error messages above")