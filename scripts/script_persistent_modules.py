#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Optimized Analysis Script for Persistent Modules

This script tests the optimized Persistent_Modules library and generates comprehensive
statistics and visualizations for elementary semi-organizations and their
irreducible generators, with special focus on the computational gains from
closure-based deduplication.

Usage:
    python analyze_persistent_modules_optimized.py

Author: Based on theoretical work by Tomas Veloz et al.
Enhanced with closure-based optimization
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter, defaultdict
import time
import json
from pathlib import Path

# Import the reaction network library and persistent modules
from pyCOT.io.functions import *
from pyCOT.ERC_Hierarchy import ERC, ERC_Hierarchy
from pyCOT.ERC_Synergy_Complementarity import build_erc_sorn
from pyCOT.Persistent_Modules import (
    compute_persistent_modules,
    build_irreducible_generators,
    IrreducibleGenerator, 
    ElementarySO,
    identify_p_ercs,
    analyze_generator_statistics,
    compute_elementary_sos
)

# ============================================================================
# ENHANCED ANALYSIS FUNCTIONS FOR OPTIMIZATION METRICS
# ============================================================================

def analyze_reduction_gains(build_stats):
    """
    Analyze the reduction gains from the optimized algorithm.
    
    Parameters
    ----------
    build_stats : dict
        Build statistics from the optimized algorithm
        
    Returns
    -------
    dict
        Detailed reduction analysis
    """
    reduction_analysis = {
        'total_explored': build_stats['total_generators_explored'],
        'unique_closures': build_stats['unique_closures_found'],
        'redundant_pruned': build_stats['redundant_generators_pruned'],
        'overall_reduction_percent': 0,
        'reduction_by_size': {},
        'efficiency_factor': 0,
        'p_ercs': build_stats.get('p_ercs', 0),
        'fundamental_synergies': build_stats.get('fundamental_synergies', 0)
    }
    
    # Calculate overall reduction
    if reduction_analysis['total_explored'] > 0:
        reduction_analysis['overall_reduction_percent'] = (
            100 * reduction_analysis['redundant_pruned'] / reduction_analysis['total_explored']
        )
        reduction_analysis['efficiency_factor'] = (
            reduction_analysis['total_explored'] / reduction_analysis['unique_closures']
        )
    
    # Analyze reduction by size
    for size, data in build_stats.get('reduction_by_size', {}).items():
        reduction_analysis['reduction_by_size'][size] = {
            'explored': data['explored'],
            'unique': data['unique'],
            'redundant': data['redundant'],
            'reduction_percent': data['reduction_ratio'] * 100
        }
    
    return reduction_analysis

def create_optimization_visualizations(build_stats, reduction_analysis, output_dir="analysis_output"):
    """
    Create visualizations specifically for optimization metrics.
    
    Parameters
    ----------
    build_stats : dict
        Build statistics from the optimized algorithm
    reduction_analysis : dict
        Reduction analysis results
    output_dir : str
        Directory to save visualizations
    """
    Path(output_dir).mkdir(exist_ok=True)
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # 1. Overall reduction pie chart
    ax1 = axes[0, 0]
    unique = reduction_analysis['unique_closures']
    redundant = reduction_analysis['redundant_pruned']
    
    if redundant > 0:
        ax1.pie([unique, redundant], 
                labels=['Unique Closures', 'Redundant Eliminated'],
                autopct='%1.1f%%',
                colors=['#2ecc71', '#e74c3c'],
                startangle=90)
        ax1.set_title(f'Generator Efficiency\n{reduction_analysis["efficiency_factor"]:.1f}x reduction')
    else:
        ax1.text(0.5, 0.5, 'No redundancy detected', 
                horizontalalignment='center', verticalalignment='center')
        ax1.set_title('Generator Efficiency')
    
    # 2. Reduction by size
    ax2 = axes[0, 1]
    if reduction_analysis['reduction_by_size']:
        sizes = sorted(reduction_analysis['reduction_by_size'].keys())
        explored = [reduction_analysis['reduction_by_size'][s]['explored'] for s in sizes]
        unique = [reduction_analysis['reduction_by_size'][s]['unique'] for s in sizes]
        
        x = np.arange(len(sizes))
        width = 0.35
        
        bars1 = ax2.bar(x - width/2, explored, width, label='Explored', alpha=0.7, color='#e74c3c')
        bars2 = ax2.bar(x + width/2, unique, width, label='Unique', alpha=0.7, color='#2ecc71')
        
        ax2.set_xlabel('Generator Size')
        ax2.set_ylabel('Count')
        ax2.set_title('Generators Explored vs Unique by Size')
        ax2.set_xticks(x)
        ax2.set_xticklabels([f'Size {s}' for s in sizes])
        ax2.legend()
        
        # Add value labels on bars
        for bars in [bars1, bars2]:
            for bar in bars:
                height = bar.get_height()
                if height > 0:
                    ax2.text(bar.get_x() + bar.get_width()/2., height,
                            f'{int(height)}', ha='center', va='bottom')
    
    # 3. Reduction percentage by size
    ax3 = axes[1, 0]
    if reduction_analysis['reduction_by_size']:
        sizes = sorted(reduction_analysis['reduction_by_size'].keys())
        reduction_percents = [reduction_analysis['reduction_by_size'][s]['reduction_percent'] for s in sizes]
        
        bars = ax3.bar([f'Size {s}' for s in sizes], reduction_percents, color='#3498db')
        ax3.set_ylabel('Reduction (%)')
        ax3.set_title('Redundancy Reduction by Generator Size')
        ax3.set_ylim(0, 100)
        
        # Add percentage labels
        for bar, pct in zip(bars, reduction_percents):
            if pct > 0:
                ax3.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 1,
                        f'{pct:.1f}%', ha='center', va='bottom')
    
    # 4. Construction phases
    ax4 = axes[1, 1]
    phases = ['P-ERCs', 'Fund. Synergies', 'Extensions']
    values = [
        build_stats.get('p_ercs', 0),
        build_stats.get('fundamental_synergies', 0),
        reduction_analysis['unique_closures'] - build_stats.get('p_ercs', 0) - build_stats.get('fundamental_synergies', 0)
    ]
    
    ax4.bar(phases, values, color=['#9b59b6', '#f39c12', '#1abc9c'])
    ax4.set_title('Unique Closures by Construction Phase')
    ax4.set_ylabel('Count')
    
    for i, (phase, value) in enumerate(zip(phases, values)):
        ax4.text(i, value + 0.5, str(value), ha='center', va='bottom')
    
    plt.suptitle('Closure-Based Optimization Analysis', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'{output_dir}/optimization_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Optimization visualizations saved to {output_dir}/")

def profile_irreducible_generators(generators):
    """
    Profile irreducible generators in terms of synergy and complementarity usage.
    
    Parameters
    ----------
    generators : list of IrreducibleGenerator
        List of irreducible generators to profile
        
    Returns
    -------
    dict
        Detailed profiling statistics
    """
    profile_data = {
        'total_generators': len(generators),
        'size_distribution': Counter(),
        'synergy_distribution': Counter(),
        'complementarity_distribution': {'type1': Counter(), 'type2': Counter(), 'type3': Counter()},
        'construction_patterns': [],
        'ssm_generators': 0,
        'pathway_details': []
    }
    
    for i, gen in enumerate(generators):
        size = gen.size()
        synergy_count = gen.get_synergy_count()
        comp_counts = gen.get_complementarity_counts()
        
        # Size distribution
        profile_data['size_distribution'][size] += 1
        
        # Synergy distribution
        profile_data['synergy_distribution'][synergy_count] += 1
        
        # Complementarity distribution
        for comp_type in ['type1', 'type2', 'type3']:
            profile_data['complementarity_distribution'][comp_type][comp_counts[comp_type]] += 1
        
        # Construction pattern
        pattern = []
        for step in gen.construction_path:
            if step['step_type'] == 'synergy':
                pattern.append('S')
            else:
                pattern.append(f"C{step['details'].get('comp_type', 1)}")
        
        profile_data['construction_patterns'].append(''.join(pattern))
        
        # Pathway details
        pathway_detail = {
            'generator_id': i,
            'size': size,
            'synergies': synergy_count,
            'complementarities': comp_counts,
            'pattern': ''.join(pattern),
            'erc_labels': [erc.label for erc in gen.erc_sequence]
        }
        profile_data['pathway_details'].append(pathway_detail)
    
    return profile_data

def analyze_elementary_so_structure(elementary_sos):
    """
    Analyze the structure of elementary semi-organizations.
    
    Parameters
    ----------
    elementary_sos : list of ElementarySO
        Elementary semi-organizations to analyze
        
    Returns
    -------
    dict
        Structural analysis results
    """
    structure_data = {
        'total_elementary_sos': len(elementary_sos),
        'p_erc_analysis': {'count': 0, 'labels': []},
        'multi_erc_analysis': {'count': 0, 'size_distribution': Counter()},
        'erc_frequency': Counter(),
        'generator_diversity': [],
        'size_vs_generators': [],
        'synergy_vs_complementarity': {'synergy_only': 0, 'complementarity_only': 0, 'mixed': 0, 'neither': 0}
    }
    
    for so in elementary_sos:
        # P-ERC analysis
        if so.is_p_erc:
            structure_data['p_erc_analysis']['count'] += 1
            structure_data['p_erc_analysis']['labels'].append(so.constituent_ercs[0].label)
        else:
            structure_data['multi_erc_analysis']['count'] += 1
            structure_data['multi_erc_analysis']['size_distribution'][len(so.constituent_ercs)] += 1
        
        # ERC frequency
        for erc in so.constituent_ercs:
            structure_data['erc_frequency'][erc.label] += 1
        
        # Generator diversity
        pathway_count = len(so.generating_pathways)
        structure_data['generator_diversity'].append(pathway_count)
        
        # Size vs generators relationship
        structure_data['size_vs_generators'].append({
            'closure_size': len(so.closure_species),
            'erc_count': len(so.constituent_ercs),
            'pathway_count': pathway_count,
            'is_p_erc': so.is_p_erc
        })
        
        # Synergy vs complementarity usage
        so_stats = so.get_generation_statistics()
        has_synergy = so_stats.get('synergy_stats', {}).get('max', 0) > 0
        has_comp1 = any(count > 0 for count in so_stats.get('complementarity_stats', {}).get('type1', {}).get('counts', []))
        has_comp2 = any(count > 0 for count in so_stats.get('complementarity_stats', {}).get('type2', {}).get('counts', []))
        has_comp3 = any(count > 0 for count in so_stats.get('complementarity_stats', {}).get('type3', {}).get('counts', []))
        has_complementarity = has_comp1 or has_comp2 or has_comp3
        
        if has_synergy and has_complementarity:
            structure_data['synergy_vs_complementarity']['mixed'] += 1
        elif has_synergy:
            structure_data['synergy_vs_complementarity']['synergy_only'] += 1
        elif has_complementarity:
            structure_data['synergy_vs_complementarity']['complementarity_only'] += 1
        else:
            structure_data['synergy_vs_complementarity']['neither'] += 1
    
    return structure_data

def create_standard_visualizations(elementary_sos, generators, statistics, output_dir="analysis_output"):
    """Create comprehensive visualizations of the analysis results."""
    Path(output_dir).mkdir(exist_ok=True)
    
    # Set style
    plt.style.use('seaborn-v0_8')
    sns.set_palette("husl")
    
    # Figure 1: Generator Size Distribution
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # Generator size distribution
    size_dist = Counter(gen.size() for gen in generators)
    axes[0, 0].bar(size_dist.keys(), size_dist.values())
    axes[0, 0].set_title('Generator Size Distribution')
    axes[0, 0].set_xlabel('Generator Size (# ERCs)')
    axes[0, 0].set_ylabel('Count')
    
    # Synergy usage distribution
    synergy_usage = statistics.get('synergy_usage', [])
    if synergy_usage:
        max_synergies = max(synergy_usage)
        bins = max(1, max_synergies)
        axes[0, 1].hist(synergy_usage, bins=bins, alpha=0.7)
    else:
        axes[0, 1].text(0.5, 0.5, 'No synergy data available',
                       ha='center', va='center',
                       transform=axes[0, 1].transAxes)
    axes[0, 1].set_title('Synergy Usage Distribution')
    axes[0, 1].set_xlabel('Number of Synergies per Generator')
    axes[0, 1].set_ylabel('Count')
    
    # Complementarity type usage
    comp_data = []
    comp_labels = []
    for comp_type in ['type1', 'type2', 'type3']:
        comp_usage = statistics.get('complementarity_usage', {}).get(comp_type, [])
        comp_data.append(len([x for x in comp_usage if x > 0]))
        comp_labels.append(f'Type {comp_type[-1]}')
    
    if any(comp_data):
        axes[1, 0].bar(comp_labels, comp_data)
    else:
        axes[1, 0].text(0.5, 0.5, 'No complementarity data available',
                       ha='center', va='center',
                       transform=axes[1, 0].transAxes)
    axes[1, 0].set_title('Complementarity Type Usage')
    axes[1, 0].set_ylabel('Count of Generators Using Type')
    
    # P-ERC vs Multi-ERC distribution
    p_erc_count = statistics.get('p_erc_count', 0)
    multi_erc_count = statistics.get('multi_erc_count', 0)
    
    if p_erc_count > 0 or multi_erc_count > 0:
        axes[1, 1].pie([max(p_erc_count, 0.1), max(multi_erc_count, 0.1)],
                      labels=['P-ERCs', 'Multi-ERC SOs'],
                      autopct='%1.1f%%',
                      startangle=90)
    else:
        axes[1, 1].text(0.5, 0.5, 'No ERC data available',
                       ha='center', va='center',
                       transform=axes[1, 1].transAxes)
    axes[1, 1].set_title('Elementary SO Composition')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/generator_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Figure 2: ERC Distribution and Pathways
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # Most common ERCs
    erc_dist = statistics.get('erc_distribution', Counter())
    most_common_ercs = erc_dist.most_common(10)
    if most_common_ercs:
        labels, counts = zip(*most_common_ercs)
        axes[0, 0].barh(range(len(labels)), counts)
        axes[0, 0].set_yticks(range(len(labels)))
        axes[0, 0].set_yticklabels(labels)
    else:
        axes[0, 0].text(0.5, 0.5, 'No ERC frequency data available',
                       ha='center', va='center',
                       transform=axes[0, 0].transAxes)
    axes[0, 0].set_title('Most Frequent ERCs in Elementary SOs')
    axes[0, 0].set_xlabel('Frequency')
    
    # Generator count per SO
    gen_counts = statistics.get('generator_counts', [])
    if gen_counts:
        max_count = max(gen_counts)
        bins = max(1, max_count)
        axes[0, 1].hist(gen_counts, bins=bins, alpha=0.7)
    else:
        axes[0, 1].text(0.5, 0.5, 'No generator count data available',
                       ha='center', va='center',
                       transform=axes[0, 1].transAxes)
    axes[0, 1].set_title('Generators per Elementary SO')
    axes[0, 1].set_xlabel('Number of Generators')
    axes[0, 1].set_ylabel('Number of SOs')
    
    # Size ratio distribution
    size_ratios = statistics.get('size_ratios', [])
    if size_ratios:
        axes[1, 0].hist(size_ratios, bins=20, alpha=0.7)
    else:
        axes[1, 0].text(0.5, 0.5, 'No size ratio data available',
                       ha='center', va='center',
                       transform=axes[1, 0].transAxes)
    axes[1, 0].set_title('Generator Size / SO Size Ratio')
    axes[1, 0].set_xlabel('Size Ratio')
    axes[1, 0].set_ylabel('Count')
    
    # Pathway diversity
    pathway_div = statistics.get('pathway_diversity', [])
    if pathway_div:
        max_div = max(pathway_div)
        bins = max(1, max_div)
        axes[1, 1].hist(pathway_div, bins=bins, alpha=0.7)
    else:
        axes[1, 1].text(0.5, 0.5, 'No pathway diversity data available',
                       ha='center', va='center',
                       transform=axes[1, 1].transAxes)
    axes[1, 1].set_title('Pathway Diversity per SO')
    axes[1, 1].set_xlabel('Number of Different Pathways')
    axes[1, 1].set_ylabel('Number of SOs')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/pathway_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Standard visualizations saved to {output_dir}/")

def generate_enhanced_report(elementary_sos, generators, statistics, profile_data, structure_data, 
                           build_stats, reduction_analysis, 
                           erc_sorn=None, output_dir="analysis_output"):
    """Generate an enhanced analysis report with optimization metrics."""
    Path(output_dir).mkdir(exist_ok=True)
    
    # Initialize default values for all potentially missing statistics
    statistics = statistics or {}
    profile_data = profile_data or {}
    structure_data = structure_data or {
        'total_elementary_sos': 0,
        'p_erc_analysis': {'count': 0, 'labels': []},
        'multi_erc_analysis': {'count': 0, 'size_distribution': {}},
        'synergy_vs_complementarity': {'synergy_only': 0, 'complementarity_only': 0, 'mixed': 0, 'neither': 0}
    }
    build_stats = build_stats or {}
    reduction_analysis = reduction_analysis or {
        'total_explored': 0,
        'unique_closures': 0,
        'redundant_pruned': 0,
        'overall_reduction_percent': 0.0,
        'efficiency_factor': 1.0
    }

    # Get summary with defaults
    summary = statistics.get('summary', {})
    
    report = []
    report.append("=" * 80)
    report.append("PERSISTENT MODULES ANALYSIS REPORT (OPTIMIZED)")
    report.append("=" * 80)
    report.append("")
    
    # Optimization Summary
    report.append("CLOSURE-BASED OPTIMIZATION SUMMARY")
    report.append("-" * 40)
    report.append(f"Total generators explored: {reduction_analysis.get('total_explored', 0):,}")
    report.append(f"Unique closures found: {reduction_analysis.get('unique_closures', 0):,}")
    report.append(f"Redundant generators eliminated: {reduction_analysis.get('redundant_pruned', 0):,}")
    report.append(f"Overall reduction: {reduction_analysis.get('overall_reduction_percent', 0.0):.1f}%")
    report.append(f"Efficiency factor: {reduction_analysis.get('efficiency_factor', 1.0):.1f}x")
    report.append("")
    
    # Summary Statistics
    report.append("SUMMARY STATISTICS")
    report.append("-" * 40)
    report.append(f"Total Elementary Semi-Organizations: {structure_data.get('total_elementary_sos', 0)}")
    report.append(f"  - P-ERCs: {structure_data.get('p_erc_analysis', {}).get('count', 0)}")
    report.append(f"  - Multi-ERC SOs: {structure_data.get('multi_erc_analysis', {}).get('count', 0)}")
    report.append(f"Total Unique Generators (after deduplication): {len(generators)}")
    report.append(f"Total Generation Pathways: {summary.get('total_pathways', 0)}")
    report.append("")
    
    # Generator Analysis
    report.append("GENERATOR ANALYSIS")
    report.append("-" * 40)
    
    # Get summary statistics with defaults
    avg_generators = summary.get('avg_generators_per_so', float('nan'))
    avg_synergies = summary.get('avg_synergies_per_path', float('nan'))
    avg_size_ratio = summary.get('avg_size_ratio', float('nan'))
    
    report.append(f"Average Generators per SO: {avg_generators if not np.isnan(avg_generators) else 'N/A'}")
    report.append(f"Average Synergies per Pathway: {avg_synergies if not np.isnan(avg_synergies) else 'N/A'}")
    report.append(f"Average Size Ratio (Generator/SO): {avg_size_ratio if not np.isnan(avg_size_ratio) else 'N/A'}")
    report.append("")
    
    # Size Distribution
    report.append("SIZE DISTRIBUTION")
    report.append("-" * 40)
    size_dist = profile_data.get('size_distribution', Counter())
    for size, count in sorted(size_dist.items()):
        report.append(f"Size {size}: {count} generators")
    report.append("")
    
    # Construction Patterns
    report.append("CONSTRUCTION PATTERNS (Top 10)")
    report.append("-" * 40)
    patterns = profile_data.get('construction_patterns', [])
    pattern_counts = Counter(patterns)
    for pattern, count in pattern_counts.most_common(10):
        report.append(f"{pattern}: {count} generators")
    report.append("(S=Synergy, C1=Complementarity Type 1, C2=Type 2, C3=Type 3)")
    report.append("")
    
    # ERC Frequency
    report.append("MOST FREQUENT ERCs")
    report.append("-" * 40)
    most_common = summary.get('most_common_ercs', [])
    if most_common:
        for erc_label, frequency in most_common:
            report.append(f"{erc_label}: appears in {frequency} elementary SOs")
    else:
        report.append("No ERC frequency data available")
    report.append("")
    
    # Synergy vs Complementarity
    report.append("SYNERGY VS COMPLEMENTARITY USAGE")
    report.append("-" * 40)
    svc = structure_data.get('synergy_vs_complementarity', {
        'synergy_only': 0,
        'complementarity_only': 0,
        'mixed': 0,
        'neither': 0
    })
    total = sum(svc.values())
    if total > 0:
        for category, count in svc.items():
            percentage = (count/total*100) if total > 0 else 0
            report.append(f"{category.replace('_', ' ').title()}: {count} ({percentage:.1f}%)")
    else:
        report.append("No synergy/complementarity usage data available")
    report.append("")
    
    # P-ERC Analysis
    report.append("P-ERC ANALYSIS")
    report.append("-" * 40)
    p_erc_analysis = structure_data.get('p_erc_analysis', {'count': 0, 'labels': []})
    report.append(f"Number of P-ERCs: {p_erc_analysis.get('count', 0)}")
    if p_erc_analysis.get('labels'):
        report.append("P-ERC Labels: " + ", ".join(p_erc_analysis['labels']))
    report.append("")
    
    # Multi-ERC SO Analysis
    report.append("MULTI-ERC SO ANALYSIS")
    report.append("-" * 40)
    multi_erc = structure_data.get('multi_erc_analysis', {'count': 0, 'size_distribution': {}})
    report.append(f"Number of Multi-ERC SOs: {multi_erc.get('count', 0)}")
    report.append("ERC Count Distribution:")
    for erc_count, so_count in sorted(multi_erc.get('size_distribution', {}).items()):
        report.append(f"  {erc_count} ERCs: {so_count} SOs")
    report.append("")
    
    # Save report and data
    with open(f'{output_dir}/analysis_report_optimized.txt', 'w') as f:
        f.write('\n'.join(report))
    
    # Prepare detailed data with safe defaults
    detailed_data = {
        'statistics': statistics,
        'profile_data': profile_data,
        'structure_data': structure_data,
        'build_stats': build_stats,
        'reduction_analysis': reduction_analysis,
        'optimization_gains': {
            'total_reduction_percent': reduction_analysis.get('overall_reduction_percent', 0.0),
            'efficiency_factor': reduction_analysis.get('efficiency_factor', 1.0),
            'unique_closures': reduction_analysis.get('unique_closures', 0),
            'redundant_eliminated': reduction_analysis.get('redundant_pruned', 0)
        },
        'pathway_details': profile_data.get('pathway_details', [])
    }
    
    # Add ERC_SORN data if available
    if erc_sorn is not None:
        try:
            detailed_data['erc_sorn_stats'] = erc_sorn.get_statistics()
            detailed_data['erc_sorn_analysis'] = erc_sorn.analyze_relationship_patterns()
        except AttributeError:
            detailed_data['erc_sorn_stats'] = {}
            detailed_data['erc_sorn_analysis'] = {}
    
    with open(f'{output_dir}/detailed_data_optimized.json', 'w') as f:
        json.dump(detailed_data, f, indent=2, default=str)
    
    print(f"Enhanced report saved to {output_dir}/analysis_report_optimized.txt")
    print(f"Detailed data saved to {output_dir}/detailed_data_optimized.json")

# ============================================================================
# MAIN ANALYSIS SCRIPT (ENHANCED)
# ============================================================================

def analyze_reaction_network_optimized(rn_file=None, RN=None, max_generator_size=8, 
                                     compare_algorithms_flag=False, output_dir="analysis_output"):
    """Main function to analyze a reaction network with optimized algorithm."""
    print("=" * 80)
    print("PERSISTENT MODULES ANALYSIS (OPTIMIZED)")
    print("=" * 80)
    
    # Initialize default statistics structure
    default_statistics = {
        'summary': {
            'avg_generators_per_so': float('nan'),
            'avg_synergies_per_path': float('nan'),
            'avg_size_ratio': float('nan'),
            'total_pathways': 0,
            'most_common_ercs': [],
            'p_erc_count': 0,
            'multi_erc_count': 0
        },
        'build_stats': {
            'total_generators_explored': 0,
            'unique_closures_found': 0,
            'redundant_generators_pruned': 0,
            'p_ercs': 0,
            'fundamental_synergies': 0,
            'reduction_by_size': {}
        },
        'synergy_usage': [],
        'complementarity_usage': {'type1': [], 'type2': [], 'type3': []},
        'generator_counts': [],
        'size_ratios': [],
        'pathway_diversity': [],
        'erc_distribution': Counter()
    }
    
    # Load or use provided reaction network
    if RN is None:
        if rn_file is None:
            print("Error: Must provide either rn_file or RN parameter")
            return
        print(f"Loading reaction network from {rn_file}...")
        RN = read_txt(rn_file)
    
    print(f"Reaction Network: {len(RN.species())} species, {len(RN.reactions())} reactions")
    
    # Perform optimized computation
    print("\nRunning optimized computation...")
    start_time = time.time()
    
    elementary_sos, generators, statistics, erc_sorn = compute_persistent_modules(
        RN, max_generator_size=max_generator_size, use_optimized=True
    )
    
    # Merge with defaults
    for key in default_statistics:
        if key not in statistics:
            statistics[key] = default_statistics[key]
        elif isinstance(default_statistics[key], dict):
            for subkey in default_statistics[key]:
                if subkey not in statistics[key]:
                    statistics[key][subkey] = default_statistics[key][subkey]
    
    computation_time = time.time() - start_time
    
    # Extract build statistics with defaults
    build_stats = statistics.get('build_stats', default_statistics['build_stats'])
    
    print(f"\nComputation completed in {computation_time:.2f} seconds")
    print(f"Found {len(elementary_sos)} elementary semi-organizations")
    print(f"Generated {len(generators)} unique generators (after deduplication)")
    
    # Analyze reduction gains
    print("\nAnalyzing reduction gains...")
    reduction_analysis = analyze_reduction_gains(build_stats)
    
    print(f"Reduction achieved: {reduction_analysis['overall_reduction_percent']:.1f}%")
    print(f"Efficiency factor: {reduction_analysis['efficiency_factor']:.1f}x")
    
    # Profile generators
    print("\nProfiling generators...")
    profile_data = profile_irreducible_generators(generators)
    
    # Analyze SO structure
    print("Analyzing elementary SO structure...")
    structure_data = analyze_elementary_so_structure(elementary_sos)
    
    # Create visualizations
    print("\nCreating visualizations...")
    create_standard_visualizations(elementary_sos, generators, statistics, output_dir)
    create_optimization_visualizations(build_stats, reduction_analysis, output_dir)
    
    # Generate enhanced report
    print("Generating enhanced report...")
    generate_enhanced_report(
        elementary_sos, generators, statistics, 
        profile_data, structure_data, build_stats, 
        reduction_analysis,  
        erc_sorn, output_dir
    )
    
    # Print summary using safe dictionary access
    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)
    print(f"Elementary SOs: {len(elementary_sos)}")
    print(f"Unique Generators: {len(generators)}")
    print(f"Redundant Eliminated: {reduction_analysis.get('redundant_pruned', 0)}")
    print(f"Overall Reduction: {reduction_analysis.get('overall_reduction_percent', 0.0):.1f}%")
    print(f"Average Generators per SO: {statistics.get('summary', {}).get('avg_generators_per_so', float('nan')):.2f}")
    
    print(f"\nMost common construction patterns:")
    pattern_counts = Counter(profile_data['construction_patterns'])
    for pattern, count in pattern_counts.most_common(3):
        print(f"  {pattern}: {count} generators")
    
    print(f"\nComputation time: {computation_time:.2f} seconds")
    print(f"Results saved to: {output_dir}/")
    
    return elementary_sos, generators, statistics, profile_data, structure_data, reduction_analysis, erc_sorn

# ============================================================================
# EXAMPLE USAGE AND TESTING
# ============================================================================

def create_test_reaction_network():
    """
    Create a simple test reaction network for demonstration.
    
    Returns
    -------
    ReactionNetwork
        A test reaction network
    """
    from pyCOT.rn_rustworkx import ReactionNetwork
    
    RN = ReactionNetwork()
    
    # Add some test reactions that will show redundancy
    test_reactions = [
        # Basic metabolism with multiple paths
        "r1: f1 + s1 => 2*s1 + y1",
        "r2: f2 + s2 => 2*s2 + y2", 
        "r3: f3 + s3 => 2*s3 + y3",
        "r4: y1 + s4 => 2*s4 + f1",
        "r5: y1 + y2 => f1",
        "r6: y2 + y3 => f1",
        "r7: => f1",  # Inflow (P-ERC)
        "r8: => f2",  # Inflow (P-ERC)
        "r9: => f3",  # Inflow (P-ERC)
        # Additional reactions creating redundant paths
        "r10: s1 + s2 => p1",
        "r11: s2 + s3 => p2",
        "r12: p1 + p2 => f1",
        "r13: s1 + s3 => p3",
        "r14: p3 + y1 => f1",
        # Example of a non-P-ERC that might be part of a larger SO
        "r15: a + b => c", # ERC for r15 might not be P-ERC
        "r16: c + d => a",
        "r17: => d" # P-ERC
    ]
    
    for reaction in test_reactions:
        RN.add_from_reaction_string(reaction)
    
    return RN

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Analyze persistent modules with optimized algorithm')
    parser.add_argument('--network', type=str, help='Path to reaction network file')
    parser.add_argument('--max-size', type=int, default=8, help='Maximum generator size')
    parser.add_argument('--compare', action='store_true', help='Compare original vs optimized algorithms')
    parser.add_argument('--output', type=str, default='analysis_output', help='Output directory')
    parser.add_argument('--test', action='store_true', help='Run with test network')
    
    args = parser.parse_args()
    
    if args.test:
        # Run with test network
        print("Creating test reaction network...")
        test_RN = create_test_reaction_network()
        results = analyze_reaction_network_optimized(
            RN=test_RN, 
            max_generator_size=args.max_size,
            compare_algorithms_flag=args.compare,
            output_dir=args.output
        )
    elif args.network:
        # Run with provided network file
        results = analyze_reaction_network_optimized(
            rn_file=args.network,
            max_generator_size=args.max_size,
            compare_algorithms_flag=args.compare,
            output_dir=args.output
        )
    else:
        # Default: run with example networks
        print("Running with example network...")
        file_path = 'networks/Navarino/RN_IN_05.txt'
        file_path = 'networks/testing/Farm.txt'
        file_path = 'networks/biomodels_interesting/central_ecoli.txt'
        file_path = 'networks/biomodels_interesting/BIOMD0000000652_manyOrgs.txt'

        results = analyze_reaction_network_optimized(
            rn_file=file_path,
            max_generator_size=args.max_size,
            compare_algorithms_flag=True,  # Always compare for demo
            output_dir=args.output
        )
    
    print("\nAnalysis complete!")
    print(f"Check '{args.output}/' directory for detailed results.")

