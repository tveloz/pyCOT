#!/usr/bin/env python3
"""
Optimized Efficient Fundamental Synergy Finder - CORRECTED VERSION

Fixed issues:
1. Constraint propagation only applies to pairs that exist for each target
2. Proper handling of cross-target synergy verification
3. Complete exploration without premature termination
4. Correct fundamentality verification

@author: Corrected Optimized Algorithm
"""

import numpy as np
import networkx as nx
from collections import defaultdict, deque
from itertools import combinations
from typing import List, Dict, Set, Tuple, Optional, FrozenSet
import time
from dataclasses import dataclass

# Import required modules
from pyCOT.ERC_Hierarchy import ERC, species_list_to_names, closure


@dataclass(frozen=True)
class BasePair:
    """Immutable base pair representation for efficient hashing."""
    base1: str
    base2: str
    
    def __post_init__(self):
        # Ensure canonical ordering
        if self.base1 > self.base2:
            object.__setattr__(self, 'base1', self.base2)
            object.__setattr__(self, 'base2', self.base1)
    
    def __hash__(self):
        return hash((self.base1, self.base2))


class OptimizedFundamentalSynergyFinder:
    """
    Optimized algorithm with correct constraint management and complete exploration.
    """
    
    def __init__(self, reaction_network, ercs=None, hierarchy_graph=None):
        """Initialize the optimized synergy finder."""
        self.rn = reaction_network
        self.ercs = ercs or ERC.ERCs(reaction_network)
        self.hierarchy_graph = hierarchy_graph or ERC.build_hierarchy_graph(self.ercs, reaction_network)
        
        # Build lookup structures
        self.erc_by_label = {erc.label: erc for erc in self.ercs}
        self.levels = ERC.get_node_levels(self.hierarchy_graph)
        
        # Pre-compute hierarchy relationships
        self._precompute_hierarchy()
        
        # Results and statistics
        self.fundamental_synergies = []
        self.stats = defaultdict(int)
        
        # Caches
        self.closure_cache = {}
        self.synergy_check_cache = {}
        
    def _precompute_hierarchy(self):
        """Pre-compute hierarchy relationships for O(1) lookup."""
        self.ancestors = {}
        self.descendants = {}
        
        for node in self.hierarchy_graph.nodes():
            self.ancestors[node] = nx.ancestors(self.hierarchy_graph, node)
            self.descendants[node] = nx.descendants(self.hierarchy_graph, node)
    
    def get_closure_set(self, erc: ERC) -> Set[str]:
        """Get closure as a set of species names (cached)."""
        if erc.label not in self.closure_cache:
            self.closure_cache[erc.label] = set(species_list_to_names(erc.get_closure(self.rn)))
        return self.closure_cache[erc.label]
    
    def has_partial_overlap_with_generators(self, base: ERC, target: ERC) -> bool:
        """Check if base has partial overlap with at least one minimal generator of target."""
        self.stats['partial_overlap_checks'] += 1
        
        base_closure = self.get_closure_set(base)
        
        for generator in target.min_generators:
            gen_species = set(species_list_to_names(generator))
            intersection = base_closure & gen_species
            
            if intersection and not gen_species.issubset(base_closure):
                return True
        
        return False
    
    def build_base_target_relationships(self, verbose: bool = False) -> Dict[str, List[str]]:
        """Build base-target relationships through generator partial overlap analysis."""
        if verbose:
            print("\nBuilding base-target relationships...")
        
        base_target_map = defaultdict(list)
        
        for target in self.ercs:
            for base in self.ercs:
                # Skip if base is the target itself
                if base.label == target.label:
                    continue
                
                # Skip if base contains target
                if base.label in self.ancestors.get(target.label, set()):
                    continue
                
                # Check partial overlap with target's minimal generators
                if self.has_partial_overlap_with_generators(base, target):
                    base_target_map[target.label].append(base.label)
                    self.stats['valid_base_target_pairs'] += 1
        
        if verbose:
            print(f"  Valid base-target pairs: {self.stats['valid_base_target_pairs']}")
        
        return base_target_map
    
    def generate_target_pairs(self, base_target_map: Dict[str, List[str]], 
                            verbose: bool = False) -> Dict[str, List[Tuple[str, str]]]:
        """Generate base pairs for each target."""
        if verbose:
            print("\nGenerating target-centric base pairs...")
        
        target_pairs = defaultdict(list)
        
        for target, potential_bases in base_target_map.items():
            if len(potential_bases) < 2:
                continue
            
            for i, base1 in enumerate(potential_bases):
                for j in range(i + 1, len(potential_bases)):
                    base2 = potential_bases[j]
                    
                    # Check that bases don't contain each other
                    if base1 in self.ancestors.get(base2, set()):
                        continue
                    if base2 in self.ancestors.get(base1, set()):
                        continue
                    
                    # Store in canonical order
                    pair = (base1, base2) if base1 < base2 else (base2, base1)
                    target_pairs[target].append(pair)
                    self.stats['base_pairs_generated'] += 1
        
        # Sort pairs for each target by combined level (bottom-up)
        for target, pairs in target_pairs.items():
            pairs.sort(key=lambda p: self.levels[p[0]] + self.levels[p[1]])
        
        if verbose:
            print(f"  Total base pairs generated: {self.stats['base_pairs_generated']}")
        
        return target_pairs
    
    def check_synergy(self, base1: ERC, base2: ERC, target: ERC) -> bool:
        """Check if base1 and base2 form a synergy for target."""
        self.stats['synergy_checks'] += 1
        
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
        """Apply constraints based on discovered synergy - CORRECTED VERSION."""
        total_pruned = 0
        
        # Get ancestors of bases
        ancestors_base1 = self.ancestors.get(base1, set())
        ancestors_base2 = self.ancestors.get(base2, set())
        
        # 1. Within-target pruning
        for pair in target_pairs[target]:
            b1, b2 = pair
            if pair in constraints[target]:
                continue
            
            should_prune = False
            
            # Check all non-minimal cases
            if ((b1 in ancestors_base1 and b2 == base2) or
                (b1 == base1 and b2 in ancestors_base2) or
                (b1 in ancestors_base1 and b2 in ancestors_base2) or
                (b2 in ancestors_base1 and b1 == base2) or
                (b2 == base1 and b1 in ancestors_base2) or
                (b2 in ancestors_base1 and b1 in ancestors_base2) or
                (b1 in ancestors_base2 and b2 in ancestors_base1)):
                should_prune = True
            
            if should_prune:
                constraints[target].add(pair)
                self.stats['within_target_pruning'] += 1
                total_pruned += 1
        
        # 2. Cross-target pruning for descendants
        descendants = self.descendants.get(target, set())
        
        for desc_target in descendants:
            if desc_target not in target_pairs:
                continue
            
            # Check if the base pair exists for this descendant target
            base_pair = (base1, base2) if base1 < base2 else (base2, base1)
            
            # Only process if this pair was actually generated for the descendant
            if base_pair in target_pairs[desc_target]:
                # Verify the synergy works for descendant
                self.stats['verification_checks'] += 1
                if self.check_synergy(self.erc_by_label[base1],
                                    self.erc_by_label[base2],
                                    self.erc_by_label[desc_target]):
                    
                    # Constrain this pair
                    if base_pair not in constraints[desc_target]:
                        constraints[desc_target].add(base_pair)
                        self.stats['cross_target_pruning'] += 1
                        total_pruned += 1
                    
                    # Also prune non-minimal versions FOR THIS TARGET
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
            
            # Filter out constrained pairs
            valid_pairs = [p for p in pairs if p not in constraints[target_label]]
            
            if verbose and len(valid_pairs) > 0:
                print(f"\nTarget {target_label} (level {self.levels[target_label]})")
                print(f"  Initial pairs: {len(pairs)}")
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
                    if ((base1.label in self.ancestors.get(other_base1.label, set()) and
                         base2.label == other_base2.label) or
                        (base1.label in self.ancestors.get(other_base2.label, set()) and
                         base2.label == other_base1.label) or
                        (base2.label in self.ancestors.get(other_base1.label, set()) and
                         base1.label == other_base2.label) or
                        (base2.label in self.ancestors.get(other_base2.label, set()) and
                         base1.label == other_base1.label)):
                        is_fundamental = False
                        if verbose:
                            print(f"  {base1.label}+{base2.label}→{target.label} not fundamental")
                            print(f"    Violated by {other_base1.label}+{other_base2.label}→{other_target.label}")
                        break
                
                # Check if other has same/smaller bases but higher target
                bases_match = ((base1.label == other_base1.label and base2.label == other_base2.label) or
                              (base1.label == other_base2.label and base2.label == other_base1.label))
                
                if bases_match and other_target.label in self.ancestors.get(target.label, set()):
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
        """Main method to find all fundamental synergies with optimizations."""
        if verbose:
            print("="*70)
            print("OPTIMIZED FUNDAMENTAL SYNERGY FINDER (CORRECTED)")
            print("="*70)
            print(f"Network: {len(self.ercs)} ERCs")
        
        start_time = time.time()
        
        # Step 1: Build base-target relationships
        base_target_map = self.build_base_target_relationships(verbose)
        
        # Step 2: Generate target-centric pairs
        target_pairs = self.generate_target_pairs(base_target_map, verbose)
        
        # Step 3: Explore targets with pruning
        discovered_synergies = self.explore_targets(target_pairs, base_target_map, verbose)
        
        # Step 4: Verify fundamentality
        fundamental_synergies = self.verify_fundamental_synergies(discovered_synergies, verbose)
        
        end_time = time.time()
        
        self.fundamental_synergies = fundamental_synergies
        
        if verbose:
            print(f"\n" + "="*70)
            print("RESULTS")
            print("="*70)
            print(f"Fundamental synergies: {len(fundamental_synergies)}")
            print(f"Time elapsed: {end_time - start_time:.3f} seconds")
            print(f"\nStatistics:")
            for key, value in sorted(self.stats.items()):
                print(f"  {key}: {value}")
            
            if fundamental_synergies:
                print(f"\nFundamental Synergies:")
                for i, synergy in enumerate(fundamental_synergies, 1):
                    print(f"  {i}. {synergy['base1'].label} + {synergy['base2'].label} → {synergy['target'].label}")
        
        return fundamental_synergies


def optimized_fundamental_synergies(reaction_network, ercs=None, hierarchy_graph=None, verbose=True):
    """Convenience function for finding fundamental synergies with optimizations."""
    finder = OptimizedFundamentalSynergyFinder(reaction_network, ercs, hierarchy_graph)
    return finder.find_fundamental_synergies(verbose=verbose)