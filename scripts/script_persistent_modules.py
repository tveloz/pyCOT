#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Enhanced Analysis Script for Organizations and Persistent Modules

This script tests both the optimized Persistent_Modules library and the new
organization detection functionality, generating comprehensive statistics and 
visualizations for elementary semi-organizations, organizations, and their
irreducible generators.

Usage:
    python analyze_organizations_enhanced.py

Author: Based on theoretical work by Tomas Veloz et al.
Enhanced with closure-based optimization and organization detection
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

# Import the new organization detection functions
try:
    from pyCOT.Persistent_Modules import (
        compute_elementary_organizations,
        check_elementary_sos_are_organizations,
        analyze_organization_patterns,
        is_self_maintaining
    )
    ORGANIZATION_FUNCTIONS_AVAILABLE = True
except ImportError:
    print("Warning: Organization detection functions not available. Please ensure organization_extensions.py is in the path.")
    ORGANIZATION_FUNCTIONS_AVAILABLE = False

# ============================================================================
# ENHANCED ANALYSIS FUNCTIONS FOR OPTIMIZATION AND ORGANIZATIONS
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

def analyze_organization_conversion(elementary_sos, organizations, flux_data=None):
    """
    Analyze the conversion from elementary SOs to organizations.
    
    Parameters
    ----------
    elementary_sos : list of ElementarySO
        All elementary semi-organizations
    organizations : list of ElementarySO
        Organizations (self-maintaining elementary SOs)
    flux_data : tuple, optional
        (flux_vectors, production_vectors) from organization computation
        
    Returns
    -------
    dict
        Conversion analysis results
    """
    conversion_analysis = {
        'total_elementary_sos': len(elementary_sos),
        'total_organizations': len(organizations),
        'conversion_rate': len(organizations) / len(elementary_sos) if elementary_sos else 0,
        'p_erc_analysis': {
            'so_count': 0,
            'org_count': 0,
            'conversion_rate': 0
        },
        'multi_erc_analysis': {
            'so_count': 0,
            'org_count': 0,
            'conversion_rate': 0
        },
        'size_distribution_comparison': {
            'sos': Counter(),
            'orgs': Counter()
        },
        'failed_sos': [],  # SOs that are not organizations
        'organization_characteristics': {
            'avg_size': 0,
            'size_range': (0, 0),
            'most_common_ercs': [],
            'flux_statistics': {}
        }
    }
    
    # Analyze P-ERC vs Multi-ERC conversion
    p_erc_sos = [so for so in elementary_sos if so.is_p_erc]
    p_erc_orgs = [org for org in organizations if org.is_p_erc]
    multi_erc_sos = [so for so in elementary_sos if not so.is_p_erc]
    multi_erc_orgs = [org for org in organizations if not org.is_p_erc]
    
    conversion_analysis['p_erc_analysis'] = {
        'so_count': len(p_erc_sos),
        'org_count': len(p_erc_orgs),
        'conversion_rate': len(p_erc_orgs) / len(p_erc_sos) if p_erc_sos else 0
    }
    
    conversion_analysis['multi_erc_analysis'] = {
        'so_count': len(multi_erc_sos),
        'org_count': len(multi_erc_orgs),
        'conversion_rate': len(multi_erc_orgs) / len(multi_erc_sos) if multi_erc_sos else 0
    }
    
    # Size distribution comparison
    for so in elementary_sos:
        size = len(so.closure_species)
        conversion_analysis['size_distribution_comparison']['sos'][size] += 1
    
    for org in organizations:
        size = len(org.closure_species)
        conversion_analysis['size_distribution_comparison']['orgs'][size] += 1
    
    # Identify failed SOs (those that are not organizations)
    org_ids = {id(org) for org in organizations}
    conversion_analysis['failed_sos'] = [so for so in elementary_sos if id(so) not in org_ids]
    
    # Organization characteristics
    if organizations:
        org_sizes = [len(org.closure_species) for org in organizations]
        conversion_analysis['organization_characteristics'] = {
            'avg_size': np.mean(org_sizes),
            'size_range': (min(org_sizes), max(org_sizes)),
            'most_common_ercs': Counter([erc.label for org in organizations for erc in org.constituent_ercs]).most_common(5),
            'flux_statistics': {}
        }
        
        # Analyze flux data if available
        if flux_data and len(flux_data) >= 2:
            flux_vectors, production_vectors = flux_data
            if flux_vectors:
                flux_norms = [np.linalg.norm(fv) if fv is not None else 0 for fv in flux_vectors]
                conversion_analysis['organization_characteristics']['flux_statistics'] = {
                    'avg_flux_norm': np.mean(flux_norms),
                    'flux_range': (min(flux_norms), max(flux_norms)),
                    'active_reactions_avg': np.mean([np.sum(fv > 1e-6) if fv is not None else 0 for fv in flux_vectors])
                }
    
    return conversion_analysis

def create_organization_visualizations(elementary_sos, organizations, conversion_analysis, output_dir="analysis_output"):
    """
    Create visualizations specifically for organization analysis.
    
    Parameters
    ----------
    elementary_sos : list of ElementarySO
        All elementary semi-organizations
    organizations : list of ElementarySO
        Organizations (self-maintaining elementary SOs)
    conversion_analysis : dict
        Results from analyze_organization_conversion
    output_dir : str
        Directory to save visualizations
    """
    Path(output_dir).mkdir(exist_ok=True)
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    
    # 1. Overall conversion pie chart
    ax1 = axes[0, 0]
    orgs = conversion_analysis['total_organizations']
    failed = conversion_analysis['total_elementary_sos'] - orgs
    
    if orgs > 0 or failed > 0:
        ax1.pie([orgs, failed], 
                labels=['Organizations', 'Failed SOs'],
                autopct='%1.1f%%',
                colors=['#2ecc71', '#e74c3c'],
                startangle=90)
        conversion_rate = conversion_analysis['conversion_rate'] * 100
        ax1.set_title(f'SO → Organization Conversion\n{conversion_rate:.1f}% success rate')
    else:
        ax1.text(0.5, 0.5, 'No data available', ha='center', va='center')
        ax1.set_title('SO → Organization Conversion')
    
    # 2. P-ERC vs Multi-ERC conversion rates
    ax2 = axes[0, 1]
    p_erc_rate = conversion_analysis['p_erc_analysis']['conversion_rate'] * 100
    multi_erc_rate = conversion_analysis['multi_erc_analysis']['conversion_rate'] * 100
    
    categories = ['P-ERCs', 'Multi-ERC SOs']
    rates = [p_erc_rate, multi_erc_rate]
    colors = ['#9b59b6', '#3498db']
    
    bars = ax2.bar(categories, rates, color=colors)
    ax2.set_ylabel('Conversion Rate (%)')
    ax2.set_title('Conversion Rate by SO Type')
    ax2.set_ylim(0, 100)
    
    # Add value labels
    for bar, rate in zip(bars, rates):
        ax2.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 1,
                f'{rate:.1f}%', ha='center', va='bottom')
    
    # 3. Size distribution comparison
    ax3 = axes[0, 2]
    size_dist_sos = conversion_analysis['size_distribution_comparison']['sos']
    size_dist_orgs = conversion_analysis['size_distribution_comparison']['orgs']
    
    if size_dist_sos or size_dist_orgs:
        all_sizes = sorted(set(size_dist_sos.keys()) | set(size_dist_orgs.keys()))
        so_counts = [size_dist_sos.get(size, 0) for size in all_sizes]
        org_counts = [size_dist_orgs.get(size, 0) for size in all_sizes]
        
        x = np.arange(len(all_sizes))
        width = 0.35
        
        bars1 = ax3.bar(x - width/2, so_counts, width, label='SOs', alpha=0.7, color='#e74c3c')
        bars2 = ax3.bar(x + width/2, org_counts, width, label='Organizations', alpha=0.7, color='#2ecc71')
        
        ax3.set_xlabel('Closure Size (# Species)')
        ax3.set_ylabel('Count')
        ax3.set_title('Size Distribution: SOs vs Organizations')
        ax3.set_xticks(x)
        ax3.set_xticklabels([f'{s}' for s in all_sizes])
        ax3.legend()
    else:
        ax3.text(0.5, 0.5, 'No size data available', ha='center', va='center')
        ax3.set_title('Size Distribution: SOs vs Organizations')
    
    # 4. Conversion rate by size
    ax4 = axes[1, 0]
    if size_dist_sos and size_dist_orgs:
        sizes = sorted(set(size_dist_sos.keys()) | set(size_dist_orgs.keys()))
        conversion_rates_by_size = []
        
        for size in sizes:
            so_count = size_dist_sos.get(size, 0)
            org_count = size_dist_orgs.get(size, 0)
            rate = (org_count / so_count * 100) if so_count > 0 else 0
            conversion_rates_by_size.append(rate)
        
        bars = ax4.bar([f'{s}' for s in sizes], conversion_rates_by_size, color='#f39c12')
        ax4.set_ylabel('Conversion Rate (%)')
        ax4.set_title('Conversion Rate by Closure Size')
        ax4.set_ylim(0, 100)
        
        # Add value labels
        for bar, rate in zip(bars, conversion_rates_by_size):
            if rate > 0:
                ax4.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 1,
                        f'{rate:.0f}%', ha='center', va='bottom')
    else:
        ax4.text(0.5, 0.5, 'No size data available', ha='center', va='center')
        ax4.set_title('Conversion Rate by Closure Size')
    
    # 5. Organization characteristics
    ax5 = axes[1, 1]
    org_chars = conversion_analysis['organization_characteristics']
    
    if organizations:
        # Show flux statistics if available
        flux_stats = org_chars.get('flux_statistics', {})
        if flux_stats:
            metrics = ['Avg Flux Norm', 'Avg Active Reactions']
            values = [
                flux_stats.get('avg_flux_norm', 0),
                flux_stats.get('active_reactions_avg', 0)
            ]
            
            ax5.bar(metrics, values, color=['#1abc9c', '#e67e22'])
            ax5.set_title('Organization Flux Characteristics')
            ax5.set_ylabel('Average Value')
            
            for i, (metric, value) in enumerate(zip(metrics, values)):
                ax5.text(i, value + max(values) * 0.01, f'{value:.2f}', 
                        ha='center', va='bottom')
        else:
            ax5.text(0.5, 0.5, f'Organizations found: {len(organizations)}\nAvg size: {org_chars["avg_size"]:.1f} species',
                    ha='center', va='center', transform=ax5.transAxes)
            ax5.set_title('Organization Summary')
    else:
        ax5.text(0.5, 0.5, 'No organizations found', ha='center', va='center')
        ax5.set_title('Organization Summary')
    
    # 6. Most common ERCs in organizations
    ax6 = axes[1, 2]
    most_common_ercs = org_chars.get('most_common_ercs', [])
    
    if most_common_ercs:
        labels, counts = zip(*most_common_ercs)
        ax6.barh(range(len(labels)), counts, color='#34495e')
        ax6.set_yticks(range(len(labels)))
        ax6.set_yticklabels(labels)
        ax6.set_xlabel('Frequency')
        ax6.set_title('Most Common ERCs in Organizations')
    else:
        ax6.text(0.5, 0.5, 'No ERC data available', ha='center', va='center')
        ax6.set_title('Most Common ERCs in Organizations')
    
    plt.suptitle('Organization Analysis', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'{output_dir}/organization_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Organization visualizations saved to {output_dir}/")

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
                           build_stats, reduction_analysis, organizations=None, conversion_analysis=None,
                           erc_sorn=None, output_dir="analysis_output"):
    """Generate an enhanced analysis report with optimization and organization metrics."""
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
    organizations = organizations or []
    conversion_analysis = conversion_analysis or {
        'total_organizations': 0,
        'conversion_rate': 0.0,
        'p_erc_analysis': {'conversion_rate': 0.0},
        'multi_erc_analysis': {'conversion_rate': 0.0}
    }

    # Get summary with defaults
    summary = statistics.get('summary', {})
    
    report = []
    report.append("=" * 80)
    report.append("ENHANCED PERSISTENT MODULES & ORGANIZATIONS ANALYSIS REPORT")
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
    report.append(f"Total Organizations: {len(organizations)}")
    report.append(f"  - P-ERC Organizations: {sum(1 for org in organizations if org.is_p_erc)}")
    report.append(f"  - Multi-ERC Organizations: {sum(1 for org in organizations if not org.is_p_erc)}")
    report.append(f"SO → Organization Conversion Rate: {conversion_analysis.get('conversion_rate', 0.0)*100:.1f}%")
    report.append(f"Total Unique Generators (after deduplication): {len(generators)}")
    report.append(f"Total Generation Pathways: {summary.get('total_pathways', 0)}")
    report.append("")
    
    # Organization Analysis
    if ORGANIZATION_FUNCTIONS_AVAILABLE:
        report.append("ORGANIZATION ANALYSIS")
        report.append("-" * 40)
        report.append(f"Organizations found: {len(organizations)}")
        report.append(f"Overall conversion rate: {conversion_analysis.get('conversion_rate', 0.0)*100:.1f}%")
        report.append(f"P-ERC conversion rate: {conversion_analysis.get('p_erc_analysis', {}).get('conversion_rate', 0.0)*100:.1f}%")
        report.append(f"Multi-ERC conversion rate: {conversion_analysis.get('multi_erc_analysis', {}).get('conversion_rate', 0.0)*100:.1f}%")
        
        if organizations:
            org_sizes = [len(org.closure_species) for org in organizations]
            report.append(f"Average organization size: {np.mean(org_sizes):.1f} species")
            report.append(f"Organization size range: {min(org_sizes)} - {max(org_sizes)} species")
        
        failed_count = structure_data.get('total_elementary_sos', 0) - len(organizations)
        report.append(f"Failed SOs (not self-maintaining): {failed_count}")
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
    
    # Save report and data - Changed to use UTF-8 encoding
    with open(f'{output_dir}/analysis_report_enhanced.txt', 'w', encoding='utf-8') as f:
        f.write('\n'.join(report))
    
    # Prepare detailed data with safe defaults
    detailed_data = {
        'statistics': statistics,
        'profile_data': profile_data,
        'structure_data': structure_data,
        'build_stats': build_stats,
        'reduction_analysis': reduction_analysis,
        'organization_data': {
            'total_organizations': len(organizations),
            'organizations': [{'closure_size': len(org.closure_species), 'is_p_erc': org.is_p_erc} for org in organizations],
            'conversion_analysis': conversion_analysis
        },
        'optimization_gains': {
            'total_reduction_percent': reduction_analysis.get('overall_reduction_percent', 0.0),
            'efficiency_factor': reduction_analysis.get('efficiency_factor', 1.0),
            'unique_closures': reduction_analysis.get('unique_closures', 0),
            'redundant_eliminated': reduction_analysis.get('redundant_pruned', 0)
        },
        'pathway_details': profile_data.get('pathway_details', [])
    }
    
    # Save JSON with UTF-8 encoding
    with open(f'{output_dir}/detailed_data_enhanced.json', 'w', encoding='utf-8') as f:
        json.dump(detailed_data, f, indent=2, default=str, ensure_ascii=False)
    
    print(f"Enhanced report saved to {output_dir}/analysis_report_enhanced.txt")
    print(f"Detailed data saved to {output_dir}/detailed_data_enhanced.json")

# ============================================================================
# MAIN ANALYSIS SCRIPT (ENHANCED WITH ORGANIZATIONS)
# ============================================================================

def analyze_reaction_network_enhanced(rn_file=None, RN=None, max_generator_size=8, 
                                    compare_algorithms_flag=False, output_dir="analysis_output"):
    """Main function to analyze a reaction network with organizations and optimized algorithm."""
    print("=" * 80)
    print("ENHANCED PERSISTENT MODULES & ORGANIZATIONS ANALYSIS")
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
    
    # ENHANCED: Compute organizations
    organizations = []
    conversion_analysis = {}
    flux_data = None
    
    if ORGANIZATION_FUNCTIONS_AVAILABLE:
        print("\nComputing organizations...")
        org_start_time = time.time()
        
        organizations, elementary_sos_unused, org_statistics, flux_data = compute_elementary_organizations(
            RN, max_generator_size=max_generator_size, verbose=True
        )
        
        org_computation_time = time.time() - org_start_time
        
        print(f"Organization computation completed in {org_computation_time:.2f} seconds")
        print(f"Found {len(organizations)} organizations from {len(elementary_sos)} elementary SOs")
        
        # Analyze conversion from SOs to organizations
        print("\nAnalyzing SO → Organization conversion...")
        conversion_analysis = analyze_organization_conversion(elementary_sos, organizations, flux_data)
        
        print(f"Conversion rate: {conversion_analysis['conversion_rate']*100:.1f}%")
        print(f"P-ERC conversion rate: {conversion_analysis['p_erc_analysis']['conversion_rate']*100:.1f}%")
        print(f"Multi-ERC conversion rate: {conversion_analysis['multi_erc_analysis']['conversion_rate']*100:.1f}%")
    else:
        print("\nWarning: Organization detection functions not available")
        print("Skipping organization computation...")
    
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
    
    if ORGANIZATION_FUNCTIONS_AVAILABLE and organizations:
        create_organization_visualizations(elementary_sos, organizations, conversion_analysis, output_dir)
    
    # Generate enhanced report
    print("Generating enhanced report...")
    generate_enhanced_report(
        elementary_sos, generators, statistics, 
        profile_data, structure_data, build_stats, 
        reduction_analysis, organizations, conversion_analysis,
        erc_sorn, output_dir
    )
    
    # Print summary using safe dictionary access
    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)
    print(f"Elementary SOs: {len(elementary_sos)}")
    if ORGANIZATION_FUNCTIONS_AVAILABLE:
        print(f"Organizations: {len(organizations)}")
        print(f"Conversion Rate: {conversion_analysis.get('conversion_rate', 0)*100:.1f}%")
    print(f"Unique Generators: {len(generators)}")
    print(f"Redundant Eliminated: {reduction_analysis.get('redundant_pruned', 0)}")
    print(f"Overall Reduction: {reduction_analysis.get('overall_reduction_percent', 0.0):.1f}%")
    print(f"Average Generators per SO: {statistics.get('summary', {}).get('avg_generators_per_so', float('nan')):.2f}")
    
    print(f"\nMost common construction patterns:")
    pattern_counts = Counter(profile_data['construction_patterns'])
    for pattern, count in pattern_counts.most_common(3):
        print(f"  {pattern}: {count} generators")
    
    total_computation_time = computation_time + (org_computation_time if ORGANIZATION_FUNCTIONS_AVAILABLE else 0)
    print(f"\nTotal computation time: {total_computation_time:.2f} seconds")
    print(f"Results saved to: {output_dir}/")
    
    return elementary_sos, generators, statistics, profile_data, structure_data, reduction_analysis, organizations, conversion_analysis, erc_sorn

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
    
    # Add some test reactions that will show redundancy and organizations
    test_reactions = [
        # Basic metabolism with multiple paths
        "r1: f1 + s1 => 2*s1 + y1",
        "r2: f2 + s2 => 2*s2 + y2", 
        "r3: f3 + s3 => 2*s3 + y3",
        "r4: y1 + s4 => 2*s4 + f1",
        "r5: y1 + y2 => f1",
        "r6: y2 + y3 => f1",
        "r7: => f1",  # Inflow (P-ERC, should be organization)
        "r8: => f2",  # Inflow (P-ERC, should be organization)
        "r9: => f3",  # Inflow (P-ERC, should be organization)
        # Additional reactions creating redundant paths
        "r10: s1 + s2 => p1",
        "r11: s2 + s3 => p2",
        "r12: p1 + p2 => f1",
        "r13: s1 + s3 => p3",
        "r14: p3 + y1 => f1",
        # Example of a potentially self-maintaining cycle
        "r15: a + b => c",
        "r16: c + d => a",
        "r17: => d",  # P-ERC (should be organization)
        "r18: => b"   # P-ERC (should be organization)
    ]
    
    for reaction in test_reactions:
        RN.add_from_reaction_string(reaction)
    
    return RN

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Analyze persistent modules and organizations with enhanced features')
    parser.add_argument('--network', type=str, help='Path to reaction network file')
    parser.add_argument('--max-size', type=int, default=8, help='Maximum generator size')
    parser.add_argument('--compare', action='store_true', help='Compare original vs optimized algorithms')
    parser.add_argument('--output', type=str, default='analysis_output_enhanced', help='Output directory')
    parser.add_argument('--test', action='store_true', help='Run with test network')
    
    args = parser.parse_args()
    
    if args.test:
        # Run with test network
        print("Creating test reaction network...")
        test_RN = create_test_reaction_network()
        results = analyze_reaction_network_enhanced(
            RN=test_RN, 
            max_generator_size=args.max_size,
            compare_algorithms_flag=args.compare,
            output_dir=args.output
        )
    elif args.network:
        # Run with provided network file
        results = analyze_reaction_network_enhanced(
            rn_file=args.network,
            max_generator_size=args.max_size,
            compare_algorithms_flag=args.compare,
            output_dir=args.output
        )
    else:
        # Default: run with example networks
        print("Running with example network...")
        #file_path = 'networks/Navarino/RN_IN_05.txt'
        file_path = 'networks/testing/Farm.txt'
        #file_path = 'networks/biomodels_interesting/central_ecoli.txt'
        #file_path = 'networks/biomodels_interesting/BIOMD0000000652_manyOrgs.txt'

        results = analyze_reaction_network_enhanced(
            rn_file=file_path,
            max_generator_size=args.max_size,
            compare_algorithms_flag=True,  # Always compare for demo
            output_dir=args.output
        )
    
    print("\nEnhanced analysis complete!")
    print(f"Check '{args.output}/' directory for detailed results.")
    
    if ORGANIZATION_FUNCTIONS_AVAILABLE:
        print("\nOrganization analysis included:")
        print("- SO → Organization conversion analysis")
        print("- Self-maintenance verification using linear programming")
        print("- Flux vector analysis for organizations")
    else:
        print("\nNote: Organization analysis not available.")
        print("To enable, ensure organization_extensions.py is in the Python path.")