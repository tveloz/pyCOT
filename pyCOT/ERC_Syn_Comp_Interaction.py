#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ERC Synergy-Complementarity Interaction Patterns

Enhanced interaction pattern tracking and analysis for combined synergy and 
complementarity relationships between ERCs.

This module focuses on:
- Combined interaction pattern analysis
- Distribution analysis of complementarities over synergies
- Distribution analysis of complementarities over non-synergetic pairs
- DataFrame-based result storage for analysis

Classes:
    ERC_Interaction_Pattern: Represents combined synergy and complementarity interaction

Functions:
    analyze_interaction_patterns: Analyze complete interaction patterns between ERCs
    analyze_all_pairwise_interactions: Analyze all ERC pairs and return DataFrame
    analyze_complementarity_distribution: Analyze how complementarities distribute

Author: Based on theoretical work by Tomas Veloz et al.
"""

import pandas as pd
from pyCOT.ERC_Synergy_Complementarity import (
    get_basic_synergies, 
    get_maximal_synergies, 
    get_fundamental_synergies,
    get_complementarity
)

# ============================================================================
# ENHANCED INTERACTION PATTERN TRACKING
# ============================================================================

class ERC_Interaction_Pattern:
    """
    Represents combined synergy and complementarity interaction between two ERCs.
    Tracks synergy and complementarity relationships for distribution analysis.
    """
    
    def __init__(self, erc1, erc2, synergies, complementarities):
        """
        Initialize an Interaction Pattern object.
        
        Parameters
        ----------
        erc1, erc2 : ERC
            ERCs involved in the interaction
        synergies : list of ERC_Synergy
            List of synergy objects found between the ERCs
        complementarities : list of ERC_Complementarity
            List of complementarity objects found between the ERCs
        """
        self.erc1 = erc1
        self.erc2 = erc2
        self.synergies = synergies  # List of ERC_Synergy objects
        self.complementarities = complementarities  # List of ERC_Complementarity objects
        
        # Analyze interaction patterns
        self.has_synergy = len(synergies) > 0
        self.has_complementarity = len(complementarities) > 0
        self.synergy_types = self._get_synergy_types()
        self.complementarity_types = self._get_complementarity_types()
        
        # Create signatures for analysis
        self.synergy_signature = self._get_synergy_signature()
        self.complementarity_signature = self._get_complementarity_signature()
        self.combined_signature = f"{self.synergy_signature}×{self.complementarity_signature}"
        self.interaction_category = self._classify_interaction()
    
    def _get_synergy_types(self):
        """Get set of synergy types present."""
        return {s.synergy_type for s in self.synergies}
    
    def _get_complementarity_types(self):
        """Get set of complementarity types present."""
        return {c.comp_type for c in self.complementarities}
    
    def _get_synergy_signature(self):
        """Create synergy signature for pattern analysis."""
        basic = sum(1 for s in self.synergies if s.synergy_type == "regular")
        maximal = sum(1 for s in self.synergies if s.synergy_type == "maximal") 
        fundamental = sum(1 for s in self.synergies if s.synergy_type == "fundamental")
        return f"S[{basic},{maximal},{fundamental}]"
    
    def _get_complementarity_signature(self):
        """Create complementarity signature for pattern analysis."""
        type1 = sum(1 for c in self.complementarities if c.comp_type == 1)
        type2 = sum(1 for c in self.complementarities if c.comp_type == 2)
        type3 = sum(1 for c in self.complementarities if c.comp_type == 3)
        return f"C[{type1},{type2},{type3}]"
    
    def _classify_interaction(self):
        """Classify the interaction type."""
        if self.has_synergy and self.has_complementarity:
            # Further classify by types
            if 1 in self.complementarity_types:
                return "super_productive"  # Type 1 + Synergy
            elif 3 in self.complementarity_types:
                return "expansion_driven"  # Type 3 + Synergy
            else:
                return "synergistic_complementary"
        elif self.has_synergy:
            if "fundamental" in self.synergy_types:
                return "fundamental_synergistic"
            else:
                return "synergistic_only"
        elif self.has_complementarity:
            if 2 in self.complementarity_types:
                return "pure_metabolic"  # Type 2 alone (no synergy)
            else:
                return "complementary_only"
        else:
            return "neutral"
    
    def get_interaction_summary(self):
        """Get a detailed summary of the interaction."""
        summary = {
            'erc_pair': f"{self.erc1.label} ⇔ {self.erc2.label}",
            'has_synergy': self.has_synergy,
            'has_complementarity': self.has_complementarity,
            'synergy_count': len(self.synergies),
            'complementarity_count': len(self.complementarities),
            'synergy_types': list(self.synergy_types),
            'complementarity_types': list(self.complementarity_types),
            'interaction_category': self.interaction_category,
            'combined_signature': self.combined_signature
        }
        return summary
    
    def get_detailed_info(self):
        """Get detailed information about all synergies and complementarities."""
        info = {
            'synergies': [],
            'complementarities': []
        }
        
        for syn in self.synergies:
            info['synergies'].append({
                'type': syn.synergy_type,
                'reactants': syn.rlabel,
                'product': syn.plabel,
                'description': str(syn)
            })
        
        for comp in self.complementarities:
            info['complementarities'].append({
                'type': comp.comp_type,
                'ercs': [comp.erc1.label, comp.erc2.label],
                'reduction_info': comp.reduction_info,
                'description': str(comp)
            })
        
        return info
    
    def __repr__(self):
        """String representation of the Interaction Pattern object."""
        return f"{self.interaction_category}: {self.combined_signature}"

def analyze_interaction_patterns(erc1, erc2, hierarchy, RN):
    """
    Analyze complete interaction patterns between two ERCs.
    
    This function computes all types of synergies and complementarities
    between two ERCs and returns a comprehensive interaction pattern analysis.
    
    Parameters
    ----------
    erc1, erc2 : ERC
        ERCs to analyze for interaction patterns
    hierarchy : ERC_Hierarchy
        Hierarchy containing the ERCs
    RN : ReactionNetwork
        The reaction network
    
    Returns
    -------
    ERC_Interaction_Pattern
        Complete interaction analysis including synergies and complementarities
    """
    # Get all synergies
    basic_syn = get_basic_synergies(erc1, erc2, hierarchy, RN)
    maximal_syn = get_maximal_synergies(erc1, erc2, hierarchy, RN)
    fundamental_syn = get_fundamental_synergies(erc1, erc2, hierarchy, RN)
    all_synergies = basic_syn + maximal_syn + fundamental_syn
    
    # Get all complementarities  
    complementarities = get_complementarity(erc1, erc2, hierarchy, RN)
    
    # Create interaction pattern
    return ERC_Interaction_Pattern(erc1, erc2, all_synergies, complementarities)

def analyze_all_pairwise_interactions(hierarchy, RN, verbose=False):
    """
    Analyze interaction patterns for all pairs of ERCs in the hierarchy.
    
    Parameters
    ----------
    hierarchy : ERC_Hierarchy
        Hierarchy containing the ERCs
    RN : ReactionNetwork
        The reaction network
    verbose : bool, optional
        If True, print progress information
    
    Returns
    -------
    pandas.DataFrame
        DataFrame with detailed interaction data for all ERC pairs
    """
    interaction_data = []
    total_pairs = len(hierarchy.ercs) * (len(hierarchy.ercs) - 1) // 2
    processed = 0
    
    if verbose:
        print(f"Analyzing interaction patterns for {total_pairs} ERC pairs...")
    
    for i, erc1 in enumerate(hierarchy.ercs):
        for j, erc2 in enumerate(hierarchy.ercs[i+1:], i+1):
            pattern = analyze_interaction_patterns(erc1, erc2, hierarchy, RN)
            
            # Extract detailed information for DataFrame
            row_data = {
                'erc1_label': erc1.label,
                'erc2_label': erc2.label,
                'erc_pair': f"{erc1.label}_{erc2.label}",
                'has_synergy': pattern.has_synergy,
                'has_complementarity': pattern.has_complementarity,
                'synergy_count': len(pattern.synergies),
                'complementarity_count': len(pattern.complementarities),
                'has_basic_synergy': 'regular' in pattern.synergy_types,
                'has_maximal_synergy': 'maximal' in pattern.synergy_types,
                'has_fundamental_synergy': 'fundamental' in pattern.synergy_types,
                'has_comp_type1': 1 in pattern.complementarity_types,
                'has_comp_type2': 2 in pattern.complementarity_types,
                'has_comp_type3': 3 in pattern.complementarity_types,
                'basic_synergy_count': sum(1 for s in pattern.synergies if s.synergy_type == "regular"),
                'maximal_synergy_count': sum(1 for s in pattern.synergies if s.synergy_type == "maximal"),
                'fundamental_synergy_count': sum(1 for s in pattern.synergies if s.synergy_type == "fundamental"),
                'comp_type1_count': sum(1 for c in pattern.complementarities if c.comp_type == 1),
                'comp_type2_count': sum(1 for c in pattern.complementarities if c.comp_type == 2),
                'comp_type3_count': sum(1 for c in pattern.complementarities if c.comp_type == 3),
                'interaction_category': pattern.interaction_category,
                'synergy_signature': pattern.synergy_signature,
                'complementarity_signature': pattern.complementarity_signature,
                'combined_signature': pattern.combined_signature
            }
            
            interaction_data.append(row_data)
            
            processed += 1
            if verbose and processed % 10 == 0:
                print(f"  Processed {processed}/{total_pairs} pairs...")
    
    df = pd.DataFrame(interaction_data)
    
    if verbose:
        active_interactions = len(df[(df['has_synergy']) | (df['has_complementarity'])])
        print(f"Found {active_interactions} active interactions out of {total_pairs} possible pairs")
    
    return df

def analyze_complementarity_distribution(df):
    """
    Analyze how complementarities distribute over synergies and non-synergetic pairs.
    
    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame from analyze_all_pairwise_interactions
    
    Returns
    -------
    dict
        Distribution analysis results
    """
    results = {}
    
    # 1. Distribution of complementarities over synergies
    synergetic_pairs = df[df['has_synergy']]
    
    if len(synergetic_pairs) > 0:
        results['synergy_complementarity_distribution'] = {}
        
        # For each synergy type, what percentage are complementary of each type
        for syn_type in ['basic', 'maximal', 'fundamental']:
            col_name = f'has_{syn_type}_synergy'
            syn_subset = synergetic_pairs[synergetic_pairs[col_name]]
            
            if len(syn_subset) > 0:
                comp_dist = {
                    'total_pairs': len(syn_subset),
                    'with_any_complementarity': len(syn_subset[syn_subset['has_complementarity']]),
                    'with_type1': len(syn_subset[syn_subset['has_comp_type1']]),
                    'with_type2': len(syn_subset[syn_subset['has_comp_type2']]),
                    'with_type3': len(syn_subset[syn_subset['has_comp_type3']]),
                    'percentage_with_any_comp': len(syn_subset[syn_subset['has_complementarity']]) / len(syn_subset) * 100,
                    'percentage_with_type1': len(syn_subset[syn_subset['has_comp_type1']]) / len(syn_subset) * 100,
                    'percentage_with_type2': len(syn_subset[syn_subset['has_comp_type2']]) / len(syn_subset) * 100,
                    'percentage_with_type3': len(syn_subset[syn_subset['has_comp_type3']]) / len(syn_subset) * 100
                }
                results['synergy_complementarity_distribution'][syn_type] = comp_dist
    
    # 2. Distribution of complementarities over non-synergetic pairs
    non_synergetic_pairs = df[~df['has_synergy']]
    
    if len(non_synergetic_pairs) > 0:
        complementary_non_syn = non_synergetic_pairs[non_synergetic_pairs['has_complementarity']]
        
        if len(complementary_non_syn) > 0:
            results['non_synergy_complementarity_distribution'] = {
                'total_non_synergetic_pairs': len(non_synergetic_pairs),
                'complementary_non_synergetic_pairs': len(complementary_non_syn),
                'percentage_of_non_syn_that_are_comp': len(complementary_non_syn) / len(non_synergetic_pairs) * 100,
                'type1_count': len(complementary_non_syn[complementary_non_syn['has_comp_type1']]),
                'type2_count': len(complementary_non_syn[complementary_non_syn['has_comp_type2']]),
                'type3_count': len(complementary_non_syn[complementary_non_syn['has_comp_type3']]),
                'type1_percentage': len(complementary_non_syn[complementary_non_syn['has_comp_type1']]) / len(complementary_non_syn) * 100,
                'type2_percentage': len(complementary_non_syn[complementary_non_syn['has_comp_type2']]) / len(complementary_non_syn) * 100,
                'type3_percentage': len(complementary_non_syn[complementary_non_syn['has_comp_type3']]) / len(complementary_non_syn) * 100
            }
    
    # 3. Overall statistics
    results['overall_statistics'] = {
        'total_pairs': len(df),
        'synergetic_pairs': len(df[df['has_synergy']]),
        'complementary_pairs': len(df[df['has_complementarity']]),
        'both_syn_and_comp': len(df[(df['has_synergy']) & (df['has_complementarity'])]),
        'synergy_only': len(df[(df['has_synergy']) & (~df['has_complementarity'])]),
        'complementarity_only': len(df[(~df['has_synergy']) & (df['has_complementarity'])]),
        'neither': len(df[(~df['has_synergy']) & (~df['has_complementarity'])]),
        'percentage_synergetic': len(df[df['has_synergy']]) / len(df) * 100,
        'percentage_complementary': len(df[df['has_complementarity']]) / len(df) * 100
    }
    
    return results

def get_interaction_statistics(df):
    """
    Get basic statistics about interaction patterns from DataFrame.
    
    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame from analyze_all_pairwise_interactions
    
    Returns
    -------
    dict
        Basic statistics about the interaction patterns
    """
    # Filter for active interactions only
    active_df = df[(df['has_synergy']) | (df['has_complementarity'])]
    
    stats = {
        'total_pairs': len(df),
        'active_interactions': len(active_df),
        'synergy_only': len(df[(df['has_synergy']) & (~df['has_complementarity'])]),
        'complementarity_only': len(df[(~df['has_synergy']) & (df['has_complementarity'])]),
        'both_synergy_complementarity': len(df[(df['has_synergy']) & (df['has_complementarity'])]),
        'category_counts': df['interaction_category'].value_counts().to_dict(),
        'synergy_type_counts': {
            'basic': df['basic_synergy_count'].sum(),
            'maximal': df['maximal_synergy_count'].sum(),
            'fundamental': df['fundamental_synergy_count'].sum()
        },
        'complementarity_type_counts': {
            'type1': df['comp_type1_count'].sum(),
            'type2': df['comp_type2_count'].sum(),
            'type3': df['comp_type3_count'].sum()
        }
    }
    
    return stats

def print_distribution_analysis(distribution_results):
    """
    Print formatted distribution analysis results.
    
    Parameters
    ----------
    distribution_results : dict
        Results from analyze_complementarity_distribution
    """
    print("="*80)
    print("COMPLEMENTARITY DISTRIBUTION ANALYSIS")
    print("="*80)
    
    # Overall statistics
    overall = distribution_results['overall_statistics']
    print(f"\n📊 OVERALL STATISTICS:")
    print(f"  Total ERC pairs analyzed: {overall['total_pairs']}")
    print(f"  Synergetic pairs: {overall['synergetic_pairs']} ({overall['percentage_synergetic']:.1f}%)")
    print(f"  Complementary pairs: {overall['complementary_pairs']} ({overall['percentage_complementary']:.1f}%)")
    print(f"  Both synergy & complementarity: {overall['both_syn_and_comp']}")
    print(f"  Synergy only: {overall['synergy_only']}")
    print(f"  Complementarity only: {overall['complementarity_only']}")
    print(f"  Neither: {overall['neither']}")
    
    # Distribution over synergies
    if 'synergy_complementarity_distribution' in distribution_results:
        print(f"\n🔄 COMPLEMENTARITY DISTRIBUTION OVER SYNERGIES:")
        syn_dist = distribution_results['synergy_complementarity_distribution']
        
        for syn_type, data in syn_dist.items():
            print(f"\n  {syn_type.upper()} SYNERGIES ({data['total_pairs']} pairs):")
            print(f"    With any complementarity: {data['with_any_complementarity']} ({data['percentage_with_any_comp']:.1f}%)")
            print(f"    With Type 1 complementarity: {data['with_type1']} ({data['percentage_with_type1']:.1f}%)")
            print(f"    With Type 2 complementarity: {data['with_type2']} ({data['percentage_with_type2']:.1f}%)")
            print(f"    With Type 3 complementarity: {data['with_type3']} ({data['percentage_with_type3']:.1f}%)")
    
    # Distribution over non-synergetic pairs
    if 'non_synergy_complementarity_distribution' in distribution_results:
        print(f"\n🔗 COMPLEMENTARITY DISTRIBUTION OVER NON-SYNERGETIC PAIRS:")
        non_syn_dist = distribution_results['non_synergy_complementarity_distribution']
        
        print(f"  Total non-synergetic pairs: {non_syn_dist['total_non_synergetic_pairs']}")
        print(f"  Complementary non-synergetic pairs: {non_syn_dist['complementary_non_synergetic_pairs']}")
        print(f"  Percentage of non-synergetic that are complementary: {non_syn_dist['percentage_of_non_syn_that_are_comp']:.1f}%")
        
        print(f"\n  DISTRIBUTION OF COMPLEMENTARITY TYPES IN NON-SYNERGETIC PAIRS:")
        print(f"    Type 1: {non_syn_dist['type1_count']} ({non_syn_dist['type1_percentage']:.1f}%)")
        print(f"    Type 2: {non_syn_dist['type2_count']} ({non_syn_dist['type2_percentage']:.1f}%)")
        print(f"    Type 3: {non_syn_dist['type3_count']} ({non_syn_dist['type3_percentage']:.1f}%)")