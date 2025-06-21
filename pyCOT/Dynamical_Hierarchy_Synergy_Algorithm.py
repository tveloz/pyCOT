#!/usr/bin/env python3
"""
Efficient Fundamental Synergy Finder - OPTIMIZED VERSION

This module implements the optimized three-step algorithm for finding fundamental
synergies by leveraging generator-based partial overlap and target-centric exploration.

The key insight: Base pairs are NOT generated globally. They are generated PER TARGET
using only bases that have partial overlap with that target's minimal generators.

Algorithm Steps:
1. Build base-target relationships via generator partial overlap
2. For each target, generate pairs from its relevant bases only  
3. Explore targets top-down with dynamic constraint propagation

@author: Optimized Efficient Synergy Algorithm
"""

import numpy as np
import networkx as nx
from collections import defaultdict, deque
from itertools import combinations
from typing import List, Dict, Set, Tuple, Optional, FrozenSet
import time

# Import required modules
from pyCOT.ERC_Hierarchy import ERC, species_list_to_names, closure
from pyCOT.Synergy_Brute_Force_Algorithm import BruteForceSynergyCalculator


class EfficientFundamentalSynergyFinder:
    """
    Optimized algorithm using the three-step approach:
    
    1. Build base-target relationships: For each (base, target) pair, check if base
       has partial overlap with any minimal generator of target.
       
    2. Organize target-centric structure: For each target, generate base pairs ONLY
       from bases that passed the generator overlap test for that specific target.
       
    3. Explore with dynamic pruning: Process targets top-down, checking only relevant
       base pairs and propagating constraints as synergies are discovered.
    """
    
    def __init__(self, reaction_network, ercs=None, hierarchy_graph=None):
        """Initialize the efficient synergy finder."""
        self.rn = reaction_network
        self.ercs = ercs or ERC.ERCs(reaction_network)
        self.hierarchy_graph = hierarchy_graph or ERC.build_hierarchy_graph(self.ercs, reaction_network)
        
        # Build lookup structures
        self.erc_by_label = {erc.label: erc for erc in self.ercs}
        self.levels = ERC.get_node_levels(self.hierarchy_graph)
        
        # Results and statistics
        self.fundamental_synergies = []
        self.stats = {
            'partial_overlap_checks': 0,
            'valid_base_target_pairs': 0,
            'base_pairs_generated': 0,
            'synergies_found': 0,
            'within_target_pruning': 0,
            'cross_target_pruning': 0,
            'verification_checks': 0,
            'final_fundamental': 0
        }
        
        # Cache for closures
        self.closure_cache = {}
        
    def get_closure_set(self, erc: ERC) -> Set[str]:
        """Get closure as a set of species names (cached)."""
        if erc.label not in self.closure_cache:
            self.closure_cache[erc.label] = set(species_list_to_names(erc.get_closure(self.rn)))
        return self.closure_cache[erc.label]
    
    def has_partial_overlap_with_generators(self, base: ERC, target: ERC) -> bool:
        """
        Check if base has partial overlap with at least one minimal generator of target.
        This is the KEY filtering step that makes the algorithm maximally efficient.
        
        A base can only contribute to a synergy if it partially covers at least one
        of the target's minimal generators (not empty intersection, not full coverage).
        """
        self.stats['partial_overlap_checks'] += 1
        
        base_closure = self.get_closure_set(base)
        
        # Check each minimal generator of the target
        for generator in target.min_generators:
            gen_species = set(species_list_to_names(generator))
            
            # Check for partial overlap with this generator
            intersection = base_closure & gen_species
            
            # Partial overlap means: has intersection but doesn't fully cover
            if intersection and not gen_species.issubset(base_closure):
                return True
        
        return False
    
    def build_base_target_relationships(self, verbose: bool = False) -> Dict[str, List[str]]:
        """
        STEP 1: Build base-target relationships through generator partial overlap analysis.
        
        For each (base, target) pair, check if base has partial overlap with at least
        one minimal generator of target. This identifies which bases can potentially
        contribute to synergies for each target.
        
        Returns a dict: target_label -> list of base_labels that can contribute.
        """
        if verbose:
            print("\nBuilding base-target relationships through generator partial overlap analysis...")
        
        base_target_map = defaultdict(list)
        
        # For each potential target
        for target in self.ercs:
            # Skip impossible cases based on hierarchy
            # Get ancestors (ERCs that contain target)
            ancestors = nx.ancestors(self.hierarchy_graph, target.label)
            
            # Check each potential base
            for base in self.ercs:
                # Skip if base is the target itself
                if base.label == target.label:
                    continue
                
                # Skip if base contains target (impossible synergy)
                if base.label in ancestors:
                    continue
                
                # Check partial overlap with target's minimal generators
                if self.has_partial_overlap_with_generators(base, target):
                    base_target_map[target.label].append(base.label)
                    self.stats['valid_base_target_pairs'] += 1
        
        if verbose:
            total_possible = len(self.ercs) * (len(self.ercs) - 1)
            print(f"  Partial overlap checks: {self.stats['partial_overlap_checks']}")
            print(f"  Valid base-target pairs: {self.stats['valid_base_target_pairs']}")
            print(f"  Reduction: {(1 - self.stats['valid_base_target_pairs']/total_possible)*100:.1f}%")
            
            # Show details for each target
            for target, bases in base_target_map.items():
                target_erc = self.erc_by_label[target]
                gens = [species_list_to_names(g) for g in target_erc.min_generators]
                print(f"  Target {target} (generators: {gens}) has {len(bases)} valid bases: {bases}")
        
        return base_target_map
    
    def generate_target_centric_structure(self, base_target_map: Dict[str, List[str]], 
                                        verbose: bool = False) -> Dict[str, List[Tuple[str, str]]]:
        """
        Organize the target-centric structure and generate base pairs.
        Key insight: Only generate pairs from bases that are relevant to each specific target.
        This is NOT generating all base pairs - it's generating them per target!
        """
        if verbose:
            print("\nOrganizing target-centric structure and generating relevant base pairs...")
        
        target_pairs = defaultdict(list)
        
        for target, potential_bases in base_target_map.items():
            if len(potential_bases) < 2:
                continue  # Need at least 2 bases for a synergy
            
            # Generate pairs ONLY from bases that can contribute to THIS target
            pairs_for_this_target = 0
            for i, base1 in enumerate(potential_bases):
                for base2 in potential_bases[i+1:]:
                    # Check that bases don't contain each other
                    if nx.has_path(self.hierarchy_graph, base1, base2):
                        continue
                    if nx.has_path(self.hierarchy_graph, base2, base1):
                        continue
                    
                    # Store in canonical order
                    if base1 < base2:
                        pair = (base1, base2)
                    else:
                        pair = (base2, base1)
                    
                    target_pairs[target].append(pair)
                    pairs_for_this_target += 1
                    self.stats['base_pairs_generated'] += 1
            
            if verbose and pairs_for_this_target > 0:
                print(f"  Target {target}: {pairs_for_this_target} pairs from {len(potential_bases)} relevant bases")
        
        # Sort pairs for each target by combined level (bottom-up)
        for target, pairs in target_pairs.items():
            pairs.sort(key=lambda p: self.levels[p[0]] + self.levels[p[1]])
        
        if verbose:
            print(f"  Total base pairs generated: {self.stats['base_pairs_generated']}")
            print(f"  Average pairs per target: {self.stats['base_pairs_generated']/len(target_pairs):.1f}")
        
        return target_pairs
    
    def check_synergy(self, base1: ERC, base2: ERC, target: ERC) -> bool:
        """Check if base1 and base2 form a synergy for target."""
        base1_closure = self.get_closure_set(base1)
        base2_closure = self.get_closure_set(base2)
        target_closure = self.get_closure_set(target)
        
        # Compute combined closure properly
        combined_species = list(set(base1.get_closure(self.rn)) | set(base2.get_closure(self.rn)))
        combined_closure = set(species_list_to_names(closure(self.rn, combined_species)))
        
        # Basic synergy conditions
        if not combined_closure.issuperset(target_closure):
            return False
        if base1_closure.issuperset(target_closure):
            return False
        if base2_closure.issuperset(target_closure):
            return False
        
        # Check for synergistic generators
        for gen in target.min_generators:
            gen_species = set(species_list_to_names(gen))
            
            if (gen_species.issubset(combined_closure) and
                not gen_species.issubset(base1_closure) and
                not gen_species.issubset(base2_closure) and
                len(gen_species & base1_closure - base2_closure) > 0 and
                len(gen_species & base2_closure - base1_closure) > 0):
                return True
        
        return False
    
    def apply_constraints_from_synergy(self, base1: str, base2: str, target: str,
                                     target_pairs: Dict[str, List[Tuple[str, str]]],
                                     constraints: Dict[str, Set[Tuple[str, str]]],
                                     base_target_map: Dict[str, List[str]]) -> int:
        """Apply constraints based on discovered synergy with hierarchy awareness."""
        total_pruned = 0
        
        # Get ancestors of bases (for non-minimality pruning)
        ancestors_base1 = nx.ancestors(self.hierarchy_graph, base1)
        ancestors_base2 = nx.ancestors(self.hierarchy_graph, base2)
        
        # 1. Within-target pruning: remove pairs with non-minimal bases
        target_bases = set(base_target_map[target])
        
        # Prune pairs where one or both bases are ancestors of the found bases
        for pair in target_pairs[target]:
            b1, b2 = pair
            if pair in constraints[target]:
                continue
                
            # Check if this pair should be pruned
            should_prune = False
            
            # Case 1: b1 is ancestor of base1, b2 is base2
            if b1 in ancestors_base1 and b2 == base2:
                should_prune = True
            # Case 2: b1 is base1, b2 is ancestor of base2
            elif b1 == base1 and b2 in ancestors_base2:
                should_prune = True
            # Case 3: b1 is ancestor of base1, b2 is ancestor of base2
            elif b1 in ancestors_base1 and b2 in ancestors_base2:
                should_prune = True
            # Case 4: Reverse cases due to commutativity
            elif b2 in ancestors_base1 and b1 == base2:
                should_prune = True
            elif b2 == base1 and b1 in ancestors_base2:
                should_prune = True
            elif b2 in ancestors_base1 and b1 in ancestors_base2:
                should_prune = True
            elif b1 in ancestors_base2 and b2 in ancestors_base1:
                should_prune = True
                
            if should_prune:
                constraints[target].add(pair)
                self.stats['within_target_pruning'] += 1
                total_pruned += 1
        
        # 2. Cross-target pruning: for targets contained in current target
        descendants = nx.descendants(self.hierarchy_graph, target)
        
        for desc_target in descendants:
            if desc_target not in target_pairs:
                continue
                
            # Verify the synergy works for contained target
            self.stats['verification_checks'] += 1
            if self.check_synergy(self.erc_by_label[base1],
                                self.erc_by_label[base2],
                                self.erc_by_label[desc_target]):
                
                # Add constraints for this target
                # The found pair itself
                pair = (base1, base2) if base1 < base2 else (base2, base1)
                if pair not in constraints[desc_target]:
                    constraints[desc_target].add(pair)
                    self.stats['cross_target_pruning'] += 1
                    total_pruned += 1
                
                # Also prune non-minimal versions
                desc_bases = set(base_target_map[desc_target])
                
                for test_pair in target_pairs[desc_target]:
                    if test_pair in constraints[desc_target]:
                        continue
                        
                    b1, b2 = test_pair
                    should_prune = False
                    
                    # Check all non-minimal cases
                    if ((b1 == base1 or b1 in ancestors_base1) and
                        (b2 == base2 or b2 in ancestors_base2)):
                        should_prune = True
                    elif ((b1 == base2 or b1 in ancestors_base2) and
                          (b2 == base1 or b2 in ancestors_base1)):
                        should_prune = True
                        
                    if should_prune:
                        constraints[desc_target].add(test_pair)
                        self.stats['cross_target_pruning'] += 1
                        total_pruned += 1
        
        return total_pruned
    
    def explore_targets(self, target_pairs: Dict[str, List[Tuple[str, str]]],
                       base_target_map: Dict[str, List[str]],
                       verbose: bool = True) -> List[Dict]:
        """Explore targets from top to bottom with dynamic pruning."""
        if verbose:
            print("\nExploring targets with dynamic pruning...")
        
        # Sort targets by level (top to bottom)
        sorted_targets = sorted(target_pairs.keys(),
                              key=lambda t: self.levels[t],
                              reverse=True)
        
        constraints = defaultdict(set)
        discovered_synergies = []
        
        for target_label in sorted_targets:
            pairs = target_pairs[target_label]
            if not pairs:
                continue
            
            if verbose:
                print(f"\nTarget {target_label} (level {self.levels[target_label]})")
                print(f"  Initial pairs: {len(pairs)}")
            
            # Filter out constrained pairs
            valid_pairs = [p for p in pairs if p not in constraints[target_label]]
            
            if verbose:
                print(f"  Valid pairs after constraints: {len(valid_pairs)}")
            
            # Check each valid pair
            for base1_label, base2_label in valid_pairs:
                base1 = self.erc_by_label[base1_label]
                base2 = self.erc_by_label[base2_label]
                target = self.erc_by_label[target_label]
                
                if self.check_synergy(base1, base2, target):
                    self.stats['synergies_found'] += 1
                    
                    synergy = {
                        'base1': base1,
                        'base2': base2,
                        'target': target
                    }
                    discovered_synergies.append(synergy)
                    
                    if verbose:
                        print(f"  ✓ Synergy: {base1_label} + {base2_label} → {target_label}")
                    
                    # Apply constraints
                    pruned = self.apply_constraints_from_synergy(
                        base1_label, base2_label, target_label,
                        target_pairs, constraints, base_target_map)
                    
                    if verbose and pruned > 0:
                        print(f"    Pruned {pruned} pairs")
        
        return discovered_synergies
    
    def verify_fundamental_synergies(self, synergies: List[Dict], verbose: bool = True) -> List[Dict]:
        """Verify that discovered synergies are truly fundamental."""
        if verbose:
            print(f"\nVerifying {len(synergies)} synergies for fundamentality...")
        
        fundamental = []
        
        for i, synergy in enumerate(synergies):
            is_fundamental = True
            base1 = synergy['base1']
            base2 = synergy['base2']
            target = synergy['target']
            
            # Check against all other synergies
            for j, other in enumerate(synergies):
                if i == j:
                    continue
                
                other_base1 = other['base1']
                other_base2 = other['base2']
                other_target = other['target']
                
                # Check if other has more minimal bases for same target
                if target.label == other_target.label:
                    # Check all cases accounting for commutativity
                    if ((nx.has_path(self.hierarchy_graph, base1.label, other_base1.label) and
                         base2.label == other_base2.label) or
                        (nx.has_path(self.hierarchy_graph, base1.label, other_base2.label) and
                         base2.label == other_base1.label) or
                        (nx.has_path(self.hierarchy_graph, base2.label, other_base1.label) and
                         base1.label == other_base2.label) or
                        (nx.has_path(self.hierarchy_graph, base2.label, other_base2.label) and
                         base1.label == other_base1.label)):
                        is_fundamental = False
                        if verbose:
                            print(f"  {base1.label}+{base2.label}→{target.label} not fundamental")
                            print(f"    Violated by {other_base1.label}+{other_base2.label}→{other_target.label}")
                        break
                
                # Check if other has same/smaller bases but higher target
                bases_match = ((base1.label == other_base1.label and base2.label == other_base2.label) or
                              (base1.label == other_base2.label and base2.label == other_base1.label))
                
                if bases_match and nx.has_path(self.hierarchy_graph, other_target.label, target.label):
                    is_fundamental = False
                    if verbose:
                        print(f"  {base1.label}+{base2.label}→{target.label} not fundamental")
                        print(f"    Target not maximal, violated by {other_target.label}")
                    break
            
            if is_fundamental:
                fundamental.append(synergy)
                self.stats['final_fundamental'] += 1
        
        return fundamental
    
    def find_fundamental_synergies(self, verbose: bool = True) -> List[Dict]:
        """
        Main method to find all fundamental synergies efficiently.
        
        The algorithm works in three steps:
        1. Identify which bases can contribute to which targets (base-target relationships)
        2. For each target, generate pairs ONLY from its relevant bases
        3. Explore targets with dynamic pruning
        
        This is fundamentally different from generating all base pairs upfront!
        """
        if verbose:
            print("="*70)
            print("OPTIMIZED EFFICIENT FUNDAMENTAL SYNERGY FINDER")
            print("="*70)
            print(f"Network: {len(self.ercs)} ERCs, {len(self.rn.species())} species")
            print(f"Hierarchy: {self.hierarchy_graph.number_of_edges()} containment relations")
        
        start_time = time.time()
        
        # Step 1: Build base-target relationships through generator partial overlap
        base_target_map = self.build_base_target_relationships(verbose)
        
        # Step 2: Organize target-centric structure (pairs will be generated lazily)
        target_pairs = self.generate_target_centric_structure(base_target_map, verbose)
        
        # Step 3: Explore targets with dynamic pruning
        discovered_synergies = self.explore_targets(target_pairs, base_target_map, verbose)
        
        # Step 4: Verify fundamentality
        fundamental_synergies = self.verify_fundamental_synergies(discovered_synergies, verbose)
        
        end_time = time.time()
        
        self.fundamental_synergies = fundamental_synergies
        
        if verbose:
            print(f"\n" + "="*70)
            print("RESULTS")
            print("="*70)
            print(f"Synergies discovered: {len(discovered_synergies)}")
            print(f"Fundamental synergies: {len(fundamental_synergies)}")
            print(f"Time elapsed: {end_time - start_time:.3f} seconds")
            print(f"\nStatistics:")
            for key, value in self.stats.items():
                print(f"  {key}: {value}")
            
            # Calculate efficiency metrics
            n_ercs = len(self.ercs)
            max_possible = n_ercs * (n_ercs - 1) * (n_ercs - 2) // 2
            actual_checks = self.stats['partial_overlap_checks'] + self.stats['synergies_found']
            
            print(f"\nEfficiency Metrics:")
            print(f"  Maximum possible synergy checks: {max_possible}")
            print(f"  Actual operations performed: {actual_checks}")
            print(f"  Reduction: {(1 - actual_checks/max_possible)*100:.1f}%")
            print(f"\nKey optimization: Checking partial overlap with MINIMAL GENERATORS")
            print(f"rather than full closures provides maximum search space reduction.")
            
            if fundamental_synergies:
                print(f"\nFundamental Synergies:")
                for i, synergy in enumerate(fundamental_synergies, 1):
                    print(f"  {i}. {synergy['base1'].label} + {synergy['base2'].label} → {synergy['target'].label}")
        
        return fundamental_synergies


def efficient_fundamental_synergies(reaction_network, ercs=None, hierarchy_graph=None, verbose=True):
    """Convenience function for finding fundamental synergies efficiently."""
    finder = EfficientFundamentalSynergyFinder(reaction_network, ercs, hierarchy_graph)
    return finder.find_fundamental_synergies(verbose=verbose)