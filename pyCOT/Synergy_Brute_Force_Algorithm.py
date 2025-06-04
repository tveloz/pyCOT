#!/usr/bin/env python3
"""
Brute Force Synergy Algorithm Module

This module introduces various algorithms for computing synergies.
It finds pairs of base1 ERCs that together can fully cover target generators, 
with eventually proper minimality enforcement.

@author: ERC Synergy Framework
"""

##########################Brute Force Synergy Calculator##########################
"""
Unified Brute Force Synergy Calculator

This module provides a unified interface for finding synergies in ERC hierarchies
using brute force approach with optional minimal filtering.

@author: Integrated Synergy Algorithm
"""

import numpy as np
import networkx as nx
from collections import defaultdict, deque
from itertools import combinations
from typing import List, Dict, Set, Tuple, Optional, Union, FrozenSet
from dataclasses import dataclass

# Import required modules from your existing codebase
from pyCOT.ERC_Hierarchy import ERC, species_list_to_names, closure


class BruteForceSynergyCalculator:
    """Unified class for finding synergies using hierarchy-aware brute force approach."""
    
    def __init__(self, reaction_network, ercs=None, hierarchy_graph=None):
        """Initialize synergy calculator with reaction network and ERCs."""
        self.rn = reaction_network
        self.ercs = ercs or ERC.ERCs(reaction_network)
        self.hierarchy_graph = hierarchy_graph or ERC.build_hierarchy_graph(self.ercs, reaction_network)
        
        # Results storage
        self.all_synergies = []
        self.minimal_synergies = []
        self.checked_pairs = set()
        
        # Statistics
        self.stats = {
            'total_checks': 0,
            'synergies_found': 0,
            'minimal_synergies_found': 0,
            'symmetric_avoided': 0,
            'minimality_analysis': {}
        }
        
    def get_erc_by_label(self, label: str) -> ERC:
        """Get ERC object by its label."""
        for erc in self.ercs:
            if erc.label == label:
                return erc
        raise ValueError(f"ERC with label {label} not found")
    
    def compute_combined_closure(self, base1: ERC, base2: ERC) -> Set[str]:
        """
        Properly compute the closure of combined generators from two ERCs.
        This is the key fix for correct synergy detection.
        """
        # Get all minimal generators from both ERCs
        combined_generators = []
        
        # Add generators from base1
        for gen in base1.min_generators:
            combined_generators.extend(gen)
        
        # Add generators from base2
        for gen in base2.min_generators:
            combined_generators.extend(gen)
        
        # Remove duplicates while preserving order
        unique_generators = []
        seen = set()
        for species in combined_generators:
            if species not in seen:
                unique_generators.append(species)
                seen.add(species)
        
        # Compute proper closure from combined generators
        combined_closure_species = closure(self.rn, unique_generators)
        combined_closure_names = set(species_list_to_names(combined_closure_species))
        
        return combined_closure_names
    
    def validates_synergy(self, base1: ERC, base2: ERC, target: ERC, verbose: bool = False) -> Tuple[bool, Dict]:
        """
        Check if base1 and base2 form a valid synergy for target.
        Uses proper closure computation for combined generators.
        
        Returns:
            Tuple[bool, Dict]: (is_synergy, details_dict)
        """
        base1_closure = set(species_list_to_names(base1.get_closure(self.rn)))
        base2_closure = set(species_list_to_names(base2.get_closure(self.rn)))
        target_closure = set(species_list_to_names(target.get_closure(self.rn)))
        
        # Compute proper combined closure
        combined = self.compute_combined_closure(base1, base2)
        
        if verbose:
            print(f"    Checking synergy {base1.label} + {base2.label} -> {target.label}")
            print(f"      Base1 closure: {sorted(base1_closure)}")
            print(f"      Base2 closure: {sorted(base2_closure)}")
            print(f"      Combined closure: {sorted(combined)}")
            print(f"      Target closure: {sorted(target_closure)}")
        
        # Basic synergy conditions
        basic_conditions = (
            combined.issuperset(target_closure) and
            not base1_closure.issuperset(target_closure) and
            not base2_closure.issuperset(target_closure)
        )
        
        if not basic_conditions:
            if verbose:
                print(f"      ✗ Basic conditions failed")
            return False, {}
        
        # Check for meaningful contribution from both bases at generator level
        covered_generators = []
        synergistic_generators = []
        
        for gen in target.min_generators:
            gen_species = set(species_list_to_names(gen))
            
            if gen_species.issubset(combined):
                covered_generators.append(gen)
                
                # Check if this generator requires true synergy (both bases contribute)
                if (not gen_species.issubset(base1_closure) and
                    not gen_species.issubset(base2_closure) and
                    len(gen_species & base1_closure - base2_closure) > 0 and
                    len(gen_species & base2_closure - base1_closure) > 0):
                    synergistic_generators.append(gen)
                    if verbose:
                        print(f"      ✓ Synergistic generator found: {species_list_to_names(gen)}")
        
        if not synergistic_generators:
            if verbose:
                print(f"      ✗ No synergistic generators found")
            return False, {}
        
        # Calculate synergy details
        coverage_ratio = len(covered_generators) / len(target.min_generators)
        synergy_ratio = len(synergistic_generators) / len(target.min_generators)
        
        details = {
            'covered_generators': covered_generators,
            'synergistic_generators': synergistic_generators,
            'coverage_ratio': coverage_ratio,
            'synergy_ratio': synergy_ratio,
            'base1_closure_size': len(base1_closure),
            'base2_closure_size': len(base2_closure),
            'target_closure_size': len(target_closure),
            'combined_closure_size': len(combined),
            'base1_closure': base1_closure,
            'base2_closure': base2_closure,
            'target_closure': target_closure,
            'combined_closure': combined
        }
        
        if verbose:
            print(f"      ✓ SYNERGY VALIDATED!")
        
        return True, details
    
    def find_all_synergies(self, verbose: bool = False) -> List[Dict]:
        """
        Find all synergies using brute force approach.
        
        Args:
            verbose: If True, print detailed information during search
            
        Returns:
            List of synergy dictionaries
        """
        if verbose:
            print("="*70)
            print("STARTING BRUTE FORCE SYNERGY SEARCH")
            print("="*70)
        
        self.all_synergies = []
        self.checked_pairs = set()
        
        # Get all possible combinations of ERCs for synergy testing
        for target in self.ercs:
            if verbose:
                print(f"\nChecking target ERC: {target.label}")
            
            # Get all other ERCs as potential bases
            other_ercs = [erc for erc in self.ercs if erc != target]
            
            # Check all pairs of other ERCs as potential synergy bases
            for i, base1 in enumerate(other_ercs):
                for j, base2 in enumerate(other_ercs[i+1:], i+1):
                    # Create ordered pair to avoid symmetric checks
                    if base1.label < base2.label:
                        ordered_base1, ordered_base2 = base1, base2
                    else:
                        ordered_base1, ordered_base2 = base2, base1
                    
                    # Create unique identifier for this combination
                    combo_id = (ordered_base1.label, ordered_base2.label, target.label)
                    
                    if combo_id in self.checked_pairs:
                        self.stats['symmetric_avoided'] += 1
                        continue
                    
                    self.checked_pairs.add(combo_id)
                    self.stats['total_checks'] += 1
                    
                    # Check synergy
                    is_synergy, details = self.validates_synergy(ordered_base1, ordered_base2, target, verbose)
                    
                    if is_synergy:
                        synergy_info = {
                            'base1': ordered_base1,
                            'base2': ordered_base2,
                            'target': target,
                            'details': details
                        }
                        self.all_synergies.append(synergy_info)
                        self.stats['synergies_found'] += 1
                        
                        if verbose:
                            print(f"  ✓ Synergy found: {ordered_base1.label} + {ordered_base2.label} → {target.label}")
        
        if verbose:
            print(f"\n" + "="*70)
            print("BRUTE FORCE SYNERGY SEARCH COMPLETE")
            print("="*70)
            print(f"Total synergies found: {len(self.all_synergies)}")
            print(f"Total combinations checked: {self.stats['total_checks']}")
            print(f"Symmetric combinations avoided: {self.stats['symmetric_avoided']}")
        
        return self.all_synergies
    
    def get_successors(self, erc_label: str) -> Set[str]:
        """Get all successors (more specific ERCs) of given ERC in hierarchy."""
        return set(self.hierarchy_graph.successors(erc_label))
    
    def get_predecessors(self, erc_label: str) -> Set[str]:
        """Get all predecessors (more general ERCs) of given ERC in hierarchy."""
        return set(self.hierarchy_graph.predecessors(erc_label))
    
    def check_minimal_base1_criterion(self, synergy: Dict, verbose: bool = False) -> Tuple[bool, List[str]]:
        """
        Check if base1 is minimal (no successor base1' such that base1'+base2->target).
        
        Returns:
            Tuple[bool, List[str]]: (is_minimal, list_of_violating_successors)
        """
        base1_label = synergy['base1'].label
        base2_label = synergy['base2'].label
        target_label = synergy['target'].label
        
        successors = self.get_successors(base1_label)
        violating_successors = []
        
        if verbose:
            print(f"    Checking base1 minimality for {base1_label}")
            print(f"      Successors to check: {list(successors)}")
        
        for successor in successors:
            successor_erc = self.get_erc_by_label(successor)
            base2_erc = self.get_erc_by_label(base2_label)
            target_erc = self.get_erc_by_label(target_label)
            
            is_synergy, _ = self.validates_synergy(successor_erc, base2_erc, target_erc)
            if is_synergy:
                violating_successors.append(successor)
                if verbose:
                    print(f"      ✗ Successor {successor} also forms synergy with {base2_label} -> {target_label}")
        
        is_minimal = len(violating_successors) == 0
        if verbose:
            print(f"      Result: {'✓ Minimal' if is_minimal else '✗ Not minimal'}")
        
        return is_minimal, violating_successors
    
    def check_minimal_base2_criterion(self, synergy: Dict, verbose: bool = False) -> Tuple[bool, List[str]]:
        """
        Check if base2 is minimal (no successor base2' such that base1+base2'->target).
        
        Returns:
            Tuple[bool, List[str]]: (is_minimal, list_of_violating_successors)
        """
        base1_label = synergy['base1'].label
        base2_label = synergy['base2'].label
        target_label = synergy['target'].label
        
        successors = self.get_successors(base2_label)
        violating_successors = []
        
        if verbose:
            print(f"    Checking base2 minimality for {base2_label}")
            print(f"      Successors to check: {list(successors)}")
        
        for successor in successors:
            base1_erc = self.get_erc_by_label(base1_label)
            successor_erc = self.get_erc_by_label(successor)
            target_erc = self.get_erc_by_label(target_label)
            
            is_synergy, _ = self.validates_synergy(base1_erc, successor_erc, target_erc)
            if is_synergy:
                violating_successors.append(successor)
                if verbose:
                    print(f"      ✗ Successor {successor} also forms synergy with {base1_label} -> {target_label}")
        
        is_minimal = len(violating_successors) == 0
        if verbose:
            print(f"      Result: {'✓ Minimal' if is_minimal else '✗ Not minimal'}")
        
        return is_minimal, violating_successors
    
    def check_maximal_target_criterion(self, synergy: Dict, verbose: bool = False) -> Tuple[bool, List[str]]:
        """
        Check if target is maximal (no predecessor target' such that base1+base2->target').
        
        Returns:
            Tuple[bool, List[str]]: (is_maximal, list_of_violating_predecessors)
        """
        base1_label = synergy['base1'].label
        base2_label = synergy['base2'].label
        target_label = synergy['target'].label
        
        predecessors = self.get_predecessors(target_label)
        violating_predecessors = []
        
        if verbose:
            print(f"    Checking target maximality for {target_label}")
            print(f"      Predecessors to check: {list(predecessors)}")
        
        for predecessor in predecessors:
            base1_erc = self.get_erc_by_label(base1_label)
            base2_erc = self.get_erc_by_label(base2_label)
            predecessor_erc = self.get_erc_by_label(predecessor)
            
            is_synergy, _ = self.validates_synergy(base1_erc, base2_erc, predecessor_erc)
            if is_synergy:
                violating_predecessors.append(predecessor)
                if verbose:
                    print(f"      ✗ Predecessor {predecessor} also achievable with {base1_label} + {base2_label}")
        
        is_maximal = len(violating_predecessors) == 0
        if verbose:
            print(f"      Result: {'✓ Maximal' if is_maximal else '✗ Not maximal'}")
        
        return is_maximal, violating_predecessors
    
    def filter_minimal_synergies(self, verbose: bool = False) -> List[Dict]:
        """
        Filter synergies to keep only minimal ones based on three criteria.
        
        Args:
            verbose: If True, print detailed information during filtering
            
        Returns:
            List of minimal synergy dictionaries
        """
        if verbose:
            print("="*70)
            print("FILTERING MINIMAL SYNERGIES")
            print("="*70)
        
        if not self.all_synergies:
            if verbose:
                print("No synergies to filter.")
            return []
        
        self.minimal_synergies = []
        minimality_analysis = {}
        
        if verbose:
            print(f"Checking {len(self.all_synergies)} synergies for minimality...")
        
        for i, synergy in enumerate(self.all_synergies):
            base1_label = synergy['base1'].label
            base2_label = synergy['base2'].label
            target_label = synergy['target'].label
            
            if verbose:
                print(f"\nSynergy {i+1}: {base1_label} + {base2_label} → {target_label}")
            
            # Check all three criteria
            criterion1, base1_violations = self.check_minimal_base1_criterion(synergy, verbose)
            criterion2, base2_violations = self.check_minimal_base2_criterion(synergy, verbose)
            criterion3, target_violations = self.check_maximal_target_criterion(synergy, verbose)
            
            # Add minimality info to synergy
            synergy_with_criteria = synergy.copy()
            synergy_with_criteria['minimality'] = {
                'base1_minimal': criterion1,
                'base2_minimal': criterion2,
                'target_maximal': criterion3,
                'is_minimal': criterion1 and criterion2 and criterion3,
                'base1_violations': base1_violations,
                'base2_violations': base2_violations,
                'target_violations': target_violations
            }
            
            # Store analysis for statistics
            analysis_key = f"{base1_label}+{base2_label}->{target_label}"
            minimality_analysis[analysis_key] = synergy_with_criteria['minimality']
            
            if criterion1 and criterion2 and criterion3:
                self.minimal_synergies.append(synergy_with_criteria)
                self.stats['minimal_synergies_found'] += 1
                if verbose:
                    print("  ✓ MINIMAL SYNERGY")
            else:
                if verbose:
                    reasons = []
                    if not criterion1:
                        reasons.append(f"base1 not minimal (violations: {base1_violations})")
                    if not criterion2:
                        reasons.append(f"base2 not minimal (violations: {base2_violations})")
                    if not criterion3:
                        reasons.append(f"target not maximal (violations: {target_violations})")
                    print(f"  ✗ Not minimal: {', '.join(reasons)}")
        
        self.stats['minimality_analysis'] = minimality_analysis
        
        if verbose:
            print(f"\n" + "="*70)
            print("MINIMAL SYNERGY FILTERING COMPLETE")
            print("="*70)
            print(f"Minimal synergies: {len(self.minimal_synergies)}")
            print(f"Filtered out: {len(self.all_synergies) - len(self.minimal_synergies)}")
            print(f"Total analyzed: {len(self.all_synergies)}")
        
        return self.minimal_synergies
    
    def brute_force(self, minimal: bool = False, verbose: bool = False) -> Union[List[Dict], Tuple[List[Dict], List[Dict]]]:
        """
        Main interface function for brute force synergy finding.
        
        Args:
            minimal: If True, also compute minimal synergies
            verbose: If True, print detailed information during computation
            
        Returns:
            If minimal=False: List of all synergies
            If minimal=True: Tuple of (all_synergies, minimal_synergies)
        """
        # Find all synergies
        all_synergies = self.find_all_synergies(verbose=verbose)
        
        if not minimal:
            return all_synergies
        
        # Also find minimal synergies
        minimal_synergies = self.filter_minimal_synergies(verbose=verbose)
        
        return all_synergies, minimal_synergies
    
    def print_synergy_summary(self, include_minimal: bool = True):
        """Print a comprehensive summary of found synergies."""
        print(f"\n" + "="*70)
        print("SYNERGY ANALYSIS SUMMARY")
        print("="*70)
        
        print(f"Total synergies found: {len(self.all_synergies)}")
        if include_minimal and self.minimal_synergies:
            print(f"Minimal synergies: {len(self.minimal_synergies)}")
            print(f"Minimality rate: {len(self.minimal_synergies)/len(self.all_synergies)*100:.1f}%" if self.all_synergies else "N/A")
        
        print(f"\nComputational Statistics:")
        print(f"  Total combinations checked: {self.stats['total_checks']}")
        print(f"  Symmetric combinations avoided: {self.stats['symmetric_avoided']}")
        print(f"  Synergy hit rate: {self.stats['synergies_found']/self.stats['total_checks']*100:.2f}%" if self.stats['total_checks'] > 0 else "N/A")
        
        if self.all_synergies:
            print(f"\nAll Synergies:")
            for i, synergy in enumerate(self.all_synergies, 1):
                base1 = synergy['base1']
                base2 = synergy['base2']
                target = synergy['target']
                details = synergy['details']
                
                print(f"  {i}. {base1.label} + {base2.label} → {target.label}")
                print(f"     Coverage: {details['coverage_ratio']:.1%}, Synergy: {details['synergy_ratio']:.1%}")
                
                if 'minimality' in synergy:
                    minimality = synergy['minimality']
                    status = "MINIMAL" if minimality['is_minimal'] else "NOT MINIMAL"
                    print(f"     Status: {status}")
        
        if include_minimal and self.minimal_synergies:
            print(f"\nMinimal Synergies Only:")
            for i, synergy in enumerate(self.minimal_synergies, 1):
                base1 = synergy['base1']
                base2 = synergy['base2']
                target = synergy['target']
                print(f"  {i}. {base1.label} + {base2.label} → {target.label}")


# Utility function for easy import
def brute_force_synergies(reaction_network, ercs=None, hierarchy_graph=None, 
                         minimal=False, verbose=False):
    """
    Convenience function for finding synergies using brute force approach.
    
    Args:
        reaction_network: The reaction network object
        ercs: List of ERCs (computed automatically if None)
        hierarchy_graph: ERC hierarchy graph (computed automatically if None)
        minimal: If True, also compute minimal synergies
        verbose: If True, print detailed information
        
    Returns:
        If minimal=False: List of all synergies
        If minimal=True: Tuple of (all_synergies, minimal_synergies)
    """
    calculator = BruteForceSynergyCalculator(reaction_network, ercs, hierarchy_graph)
    return calculator.brute_force(minimal=minimal, verbose=verbose)
