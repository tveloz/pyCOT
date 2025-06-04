#!/usr/bin/env python3
"""
Dynamic Synergy Explorer with Propagation-Based Pruning

This module implements the corrected strategy for fundamental synergy discovery:
1. Build complete target-base structure (no premature filtering)
2. Systematic exploration: highest targets first, lowest base combinations first
3. Dynamic propagation: discoveries prune exploration space optimally
4. Guarantees fundamentality by construction

@author: Corrected Synergy Framework
"""

import networkx as nx
import numpy as np
from collections import defaultdict, deque
from dataclasses import dataclass, field
from typing import List, Set, Dict, Tuple, Optional, FrozenSet
from itertools import combinations

from pyCOT.ERC_Hierarchy import ERC, species_list_to_names, closure


@dataclass
class BaseInfo:
    """Information about a base ERC associated with a target."""
    base_erc: 'ERC'
    hierarchy_level: int
    signature: frozenset
    status: str = 'active'  # 'active' | 'removed'


@dataclass
class TargetExplorationNode:
    """Complete exploration information for a target ERC."""
    target_erc: 'ERC'
    hierarchy_level: int
    closure: Set[str]
    generators: List[List[str]]
    associated_bases: Dict[str, BaseInfo] = field(default_factory=dict)
    valid_combinations: List[Tuple[str, str]] = field(default_factory=list)
    removed_combinations: Set[Tuple[str, str]] = field(default_factory=set)
    discovered_synergies: List[Tuple[str, str]] = field(default_factory=list)


@dataclass
class FundamentalSynergy:
    """A discovered fundamental synergy."""
    base1_erc: 'ERC'
    base2_erc: 'ERC' 
    target_erc: 'ERC'
    discovery_order: int
    base1_signature: frozenset
    base2_signature: frozenset


class DynamicSynergyExplorer:
    """
    Dynamic synergy explorer implementing the corrected strategy with propagation-based pruning.
    
    Key principles:
    1. Fundamentality is a property of (base1, base2, target) triples, not individual components
    2. No premature target filtering - build complete structure first
    3. Systematic exploration ensures fundamentality by construction
    4. Dynamic propagation optimally reduces search space after each discovery
    """
    
    def __init__(self, reaction_network, ercs=None, hierarchy_graph=None):
        """Initialize dynamic synergy explorer."""
        self.RN = reaction_network
        self.ercs = ercs if ercs is not None else ERC.ERCs(reaction_network)
        self.hierarchy = hierarchy_graph if hierarchy_graph is not None else ERC.build_hierarchy_graph(self.ercs, reaction_network)
        
        # Exploration state
        self.exploration_structure = {}
        self.fundamental_synergies = []
        self.discovery_count = 0
        self.levels = ERC.get_node_levels(self.hierarchy)
        
        # Statistics
        self.stats = {
            'total_targets': 0,
            'total_base_associations': 0,
            'total_combinations_initial': 0,
            'combinations_tested': 0,
            'synergies_discovered': 0,
            'combinations_pruned': 0
        }
    
    # =================== PHASE 1: COMPLETE STRUCTURE BUILDING ===================
    
    def compute_overlap_signature(self, base_erc, target_erc):
        """Compute overlap signature between base and target."""
        base_closure = set(species_list_to_names(base_erc.get_closure(self.RN)))
        signature_elements = []
        
        for i, generator in enumerate(target_erc.min_generators):
            generator_species = set(species_list_to_names(generator))
            intersection = base_closure & generator_species
            
            # Check for strict partial overlap
            if (intersection and 
                not generator_species.issubset(base_closure) and 
                not base_closure.issubset(generator_species)):
                
                element = (i, frozenset(intersection))
                signature_elements.append(element)
        
        return frozenset(signature_elements)
    
    def build_complete_structure(self, verbose=False):
        """
        Build complete exploration structure without any target filtering.
        Creates target-centric representation with associated base hierarchies.
        """
        if verbose:
            print("🏗️  BUILDING COMPLETE EXPLORATION STRUCTURE")
            print("=" * 60)
        
        # First, identify all potential base-target associations
        base_target_associations = []
        
        for base_erc in self.ercs:
            base_closure = set(species_list_to_names(base_erc.get_closure(self.RN)))
            
            # Get ERCs contained by base_erc (should be excluded as targets)
            base_contained_labels = set(nx.descendants(self.hierarchy, base_erc.label))
            
            for target_erc in self.ercs:
                if base_erc == target_erc:
                    continue
                
                # CRITICAL FIX: Exclude targets contained by base_erc
                if target_erc.label in base_contained_labels:
                    continue
                
                # Check if base has strict partial overlap with any target generator
                has_overlap = False
                for generator in target_erc.min_generators:
                    generator_species = set(species_list_to_names(generator))
                    intersection = base_closure & generator_species
                    
                    if (intersection and 
                        not generator_species.issubset(base_closure) and 
                        not base_closure.issubset(generator_species)):
                        has_overlap = True
                        break
                
                if has_overlap:
                    signature = self.compute_overlap_signature(base_erc, target_erc)
                    base_target_associations.append((base_erc, target_erc, signature))
        
        if verbose:
            print(f"Found {len(base_target_associations)} base-target associations")
        
        # Build target-centric exploration nodes
        target_bases = defaultdict(list)
        for base_erc, target_erc, signature in base_target_associations:
            target_bases[target_erc].append((base_erc, signature))
        
        for target_erc, base_list in target_bases.items():
            # Create target exploration node
            target_closure = set(species_list_to_names(target_erc.get_closure(self.RN)))
            target_generators = [species_list_to_names(gen) for gen in target_erc.min_generators]
            
            node = TargetExplorationNode(
                target_erc=target_erc,
                hierarchy_level=self.levels[target_erc.label],
                closure=target_closure,
                generators=target_generators
            )
            
            # Add associated bases
            for base_erc, signature in base_list:
                base_info = BaseInfo(
                    base_erc=base_erc,
                    hierarchy_level=self.levels[base_erc.label],
                    signature=signature,
                    status='active'
                )
                node.associated_bases[base_erc.label] = base_info
            
            # Generate valid combinations (non-containing pairs)
            active_bases = [info.base_erc for info in node.associated_bases.values()]
            node.valid_combinations = self.get_non_containing_pairs(active_bases)
            
            self.exploration_structure[target_erc.label] = node
        
        # Update statistics
        self.stats['total_targets'] = len(self.exploration_structure)
        self.stats['total_base_associations'] = len(base_target_associations)
        self.stats['total_combinations_initial'] = sum(
            len(node.valid_combinations) for node in self.exploration_structure.values()
        )
        
        if verbose:
            print(f"Built structure:")
            print(f"  Targets: {self.stats['total_targets']}")
            print(f"  Base associations: {self.stats['total_base_associations']}")
            print(f"  Initial combinations: {self.stats['total_combinations_initial']}")
        
        return self.exploration_structure
    
    def get_non_containing_pairs(self, bases):
        """Get all pairs of bases where neither contains the other."""
        valid_pairs = []
        
        for base1, base2 in combinations(bases, 2):
            # Check containment relationships
            base1_ancestors = set(nx.ancestors(self.hierarchy, base1.label))
            base2_ancestors = set(nx.ancestors(self.hierarchy, base2.label))
            
            # Neither should contain the other
            base1_contains_base2 = base2.label in nx.descendants(self.hierarchy, base1.label)
            base2_contains_base1 = base1.label in nx.descendants(self.hierarchy, base2.label)
            
            if not base1_contains_base2 and not base2_contains_base1:
                # Order consistently (lower hierarchy level first)
                if self.levels[base1.label] <= self.levels[base2.label]:
                    valid_pairs.append((base1.label, base2.label))
                else:
                    valid_pairs.append((base2.label, base1.label))
        
        return valid_pairs
    
    # =================== PHASE 2: SYSTEMATIC FUNDAMENTAL DISCOVERY ===================
    
    def validates_synergy(self, base1_erc, base2_erc, target_erc):
        """Check if base1 + base2 forms a valid synergy for target."""
        # Compute combined closure properly
        base1_generators = []
        for gen in base1_erc.min_generators:
            base1_generators.extend(gen)
        
        base2_generators = []
        for gen in base2_erc.min_generators:
            base2_generators.extend(gen)
        
        # Combine and remove duplicates
        combined_generators = []
        seen = set()
        for species in base1_generators + base2_generators:
            if species not in seen:
                combined_generators.append(species)
                seen.add(species)
        
        # Compute proper closure
        combined_closure_species = closure(self.RN, combined_generators)
        combined_closure = set(species_list_to_names(combined_closure_species))
        
        base1_closure = set(species_list_to_names(base1_erc.get_closure(self.RN)))
        base2_closure = set(species_list_to_names(base2_erc.get_closure(self.RN)))
        target_closure = set(species_list_to_names(target_erc.get_closure(self.RN)))
        
        # Basic synergy conditions
        if not (combined_closure.issuperset(target_closure) and
                not base1_closure.issuperset(target_closure) and
                not base2_closure.issuperset(target_closure)):
            return False
        
        # Check for meaningful synergy at generator level
        for generator in target_erc.min_generators:
            gen_species = set(species_list_to_names(generator))
            
            if (gen_species.issubset(combined_closure) and
                not gen_species.issubset(base1_closure) and
                not gen_species.issubset(base2_closure) and
                len(gen_species & base1_closure - base2_closure) > 0 and
                len(gen_species & base2_closure - base1_closure) > 0):
                return True
        
        return False
    
    def discover_fundamental_synergies(self, verbose=False):
        """
        Systematic discovery of fundamental synergies with dynamic propagation.
        
        Explores targets from highest to lowest level, and within each target
        tries base combinations starting from lowest levels. Fundamentality 
        is guaranteed by construction.
        """
        if not self.exploration_structure:
            self.build_complete_structure(verbose)
        
        if verbose:
            print("\n🔍 SYSTEMATIC FUNDAMENTAL SYNERGY DISCOVERY")
            print("=" * 60)
        
        # Sort targets by hierarchy level (highest first)
        sorted_targets = sorted(
            self.exploration_structure.values(),
            key=lambda node: node.hierarchy_level,
            reverse=True
        )
        
        for target_node in sorted_targets:
            if verbose:
                print(f"\n🎯 Exploring target: {target_node.target_erc.label} (Level {target_node.hierarchy_level})")
            
            # Get current valid combinations for this target
            current_combinations = [
                combo for combo in target_node.valid_combinations
                if combo not in target_node.removed_combinations
            ]
            
            if verbose:
                print(f"   Valid combinations: {len(current_combinations)}")
            
            # Sort combinations by base hierarchy levels (lowest first)
            def combination_priority(combo):
                base1_level = self.levels[combo[0]]
                base2_level = self.levels[combo[1]]
                return (base1_level + base2_level, min(base1_level, base2_level))
            
            current_combinations.sort(key=combination_priority)
            
            # Try combinations in order
            for base1_label, base2_label in current_combinations:
                base1_erc = self.get_erc_by_label(base1_label)
                base2_erc = self.get_erc_by_label(base2_label)
                
                self.stats['combinations_tested'] += 1
                
                if self.validates_synergy(base1_erc, base2_erc, target_node.target_erc):
                    # FUNDAMENTAL synergy discovered!
                    fundamental_synergy = FundamentalSynergy(
                        base1_erc=base1_erc,
                        base2_erc=base2_erc,
                        target_erc=target_node.target_erc,
                        discovery_order=self.discovery_count,
                        base1_signature=target_node.associated_bases[base1_label].signature,
                        base2_signature=target_node.associated_bases[base2_label].signature
                    )
                    
                    self.fundamental_synergies.append(fundamental_synergy)
                    target_node.discovered_synergies.append((base1_label, base2_label))
                    self.discovery_count += 1
                    self.stats['synergies_discovered'] += 1
                    
                    if verbose:
                        print(f"   ✅ FUNDAMENTAL SYNERGY: {base1_label} + {base2_label} → {target_node.target_erc.label}")
                    
                    # PROPAGATE discovery to prune exploration space
                    self.propagate_discovery(fundamental_synergy, verbose)
        
        if verbose:
            print(f"\n📊 DISCOVERY COMPLETE")
            print(f"   Fundamental synergies found: {len(self.fundamental_synergies)}")
            print(f"   Combinations tested: {self.stats['combinations_tested']}")
            print(f"   Combinations pruned: {self.stats['combinations_pruned']}")
        
        return self.fundamental_synergies
    
    # =================== PHASE 3: DYNAMIC PROPAGATION ===================
    
    def propagate_discovery(self, synergy, verbose=False):
        """
        Propagate fundamental synergy discovery to optimally prune exploration space.
        
        Two types of propagation:
        1. Base hierarchy: Remove higher bases with same signatures for current target
        2. Target hierarchy: Remove combinations involving discovered bases for contained targets
        """
        if verbose:
            print(f"   🔄 Propagating discovery...")
        
        pruned_count = 0
        
        # 1. BASE HIERARCHY PROPAGATION
        pruned_count += self.propagate_base_hierarchy(synergy, verbose)
        
        # 2. TARGET HIERARCHY PROPAGATION  
        pruned_count += self.propagate_target_hierarchy(synergy, verbose)
        
        self.stats['combinations_pruned'] += pruned_count
        
        if verbose:
            print(f"      Pruned {pruned_count} combinations")
    
    def propagate_base_hierarchy(self, synergy, verbose=False):
        """Remove higher bases with same signatures for the discovered target."""
        target_label = synergy.target_erc.label
        target_node = self.exploration_structure[target_label]
        
        pruned_count = 0
        
        # Get signatures of fundamental bases
        base1_sig = synergy.base1_signature
        base2_sig = synergy.base2_signature
        base1_level = self.levels[synergy.base1_erc.label]
        base2_level = self.levels[synergy.base2_erc.label]
        
        # Remove higher bases with same signatures
        bases_to_remove = []
        
        for base_label, base_info in target_node.associated_bases.items():
            if base_info.status == 'removed':
                continue
                
            # Check if this base should be removed
            should_remove = False
            
            if (base_info.signature == base1_sig and 
                base_info.hierarchy_level > base1_level):
                should_remove = True
                
            elif (base_info.signature == base2_sig and 
                  base_info.hierarchy_level > base2_level):
                should_remove = True
            
            if should_remove:
                bases_to_remove.append(base_label)
        
        # Remove identified bases and their combinations
        for base_label in bases_to_remove:
            target_node.associated_bases[base_label].status = 'removed'
            
            # Remove combinations involving this base
            combinations_to_remove = [
                combo for combo in target_node.valid_combinations
                if base_label in combo and combo not in target_node.removed_combinations
            ]
            
            for combo in combinations_to_remove:
                target_node.removed_combinations.add(combo)
                pruned_count += 1
        
        return pruned_count
    
    def propagate_target_hierarchy(self, synergy, verbose=False):
        """
        Remove combinations for contained targets where both bases contain discovered bases.
        
        CRITICAL FIX: Remove any combination (base1', base2') for contained targets where:
        - base1' contains (or equals) base1 AND base2' contains (or equals) base2
        
        Rationale: If base1 + base2 → target_high is fundamental, then for any target_low ⊆ target_high,
        combinations using base1' ⊇ base1 and base2' ⊇ base2 would be non-fundamental
        (more resources for lesser achievement).
        """
        target_label = synergy.target_erc.label
        base1_label = synergy.base1_erc.label
        base2_label = synergy.base2_erc.label
        
        pruned_count = 0
        
        # Get all targets contained in the discovered target
        contained_targets = nx.descendants(self.hierarchy, target_label)
        
        # Get bases that contain the discovered bases (including themselves)
        base1_containing_bases = {base1_label} | nx.ancestors(self.hierarchy, base1_label)
        base2_containing_bases = {base2_label} | nx.ancestors(self.hierarchy, base2_label)
        
        for contained_target_label in contained_targets:
            if contained_target_label not in self.exploration_structure:
                continue
            
            contained_node = self.exploration_structure[contained_target_label]
            
            # Remove combinations where both bases contain the discovered bases
            combinations_to_remove = []
            
            for combo in contained_node.valid_combinations:
                if combo in contained_node.removed_combinations:
                    continue
                    
                combo_base1, combo_base2 = combo
                
                # Check if both bases in combination contain the discovered bases
                # (in either order since combinations are symmetric)
                condition1 = (combo_base1 in base1_containing_bases and 
                             combo_base2 in base2_containing_bases)
                condition2 = (combo_base1 in base2_containing_bases and 
                             combo_base2 in base1_containing_bases)
                
                if condition1 or condition2:
                    combinations_to_remove.append(combo)
            
            # Remove identified combinations
            for combo in combinations_to_remove:
                contained_node.removed_combinations.add(combo)
                pruned_count += 1
        
        return pruned_count
    
    # =================== UTILITIES ===================
    
    def get_erc_by_label(self, label):
        """Get ERC object by label."""
        for erc in self.ercs:
            if erc.label == label:
                return erc
        raise ValueError(f"ERC with label {label} not found")
    
    def print_exploration_statistics(self):
        """Print comprehensive statistics about the exploration process."""
        print("\n" + "=" * 80)
        print("📊 DYNAMIC EXPLORATION STATISTICS")
        print("=" * 80)
        
        print(f"Structure Building:")
        print(f"  Total targets: {self.stats['total_targets']}")
        print(f"  Base associations: {self.stats['total_base_associations']}")
        print(f"  Initial combinations: {self.stats['total_combinations_initial']}")
        
        print(f"\nDiscovery Process:")
        print(f"  Combinations tested: {self.stats['combinations_tested']}")
        print(f"  Fundamental synergies found: {self.stats['synergies_discovered']}")
        print(f"  Success rate: {self.stats['synergies_discovered']/self.stats['combinations_tested']*100:.1f}%" if self.stats['combinations_tested'] > 0 else "N/A")
        
        print(f"\nPruning Efficiency:")
        print(f"  Combinations pruned: {self.stats['combinations_pruned']}")
        print(f"  Final combinations remaining: {self.stats['total_combinations_initial'] - self.stats['combinations_pruned']}")
        if self.stats['total_combinations_initial'] > 0:
            print(f"  Pruning efficiency: {self.stats['combinations_pruned']/self.stats['total_combinations_initial']*100:.1f}%")
    
    def print_fundamental_synergies(self):
        """Print discovered fundamental synergies."""
        print("\n" + "=" * 80)
        print("🎯 DISCOVERED FUNDAMENTAL SYNERGIES")
        print("=" * 80)
        
        if not self.fundamental_synergies:
            print("No fundamental synergies discovered.")
            return
        
        for i, synergy in enumerate(self.fundamental_synergies, 1):
            print(f"\n{i}. {synergy.base1_erc.label} + {synergy.base2_erc.label} → {synergy.target_erc.label}")
            print(f"   Base1 Level: {self.levels[synergy.base1_erc.label]}")
            print(f"   Base2 Level: {self.levels[synergy.base2_erc.label]}")
            print(f"   Target Level: {self.levels[synergy.target_erc.label]}")
            print(f"   Discovery Order: {synergy.discovery_order}")
            print(f"   Base1 Closure: {species_list_to_names(synergy.base1_erc.get_closure(self.RN))}")
            print(f"   Base2 Closure: {species_list_to_names(synergy.base2_erc.get_closure(self.RN))}")
            print(f"   Target Closure: {species_list_to_names(synergy.target_erc.get_closure(self.RN))}")
    
    def run_complete_analysis(self, verbose=True):
        """Run the complete dynamic synergy analysis."""
        if verbose:
            print("🚀 DYNAMIC FUNDAMENTAL SYNERGY DISCOVERY")
            print("=" * 80)
        
        # Phase 1: Build complete structure
        self.build_complete_structure(verbose)
        
        # Phase 2 & 3: Discover with propagation
        self.discover_fundamental_synergies(verbose)
        
        # Results
        if verbose:
            self.print_exploration_statistics()
            self.print_fundamental_synergies()
        
        return self.fundamental_synergies