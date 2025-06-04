#!/usr/bin/env python3
"""
ERC Exploration Module

This module provides advanced analysis capabilities for ERC hierarchies including:
- Target identification for potential synergies
- Target-centric synergy organization
- Hierarchy structural analysis
- Synergy filtering and optimization

@author: ERC Analysis Framework
"""

import networkx as nx
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict, Counter
from dataclasses import dataclass
from typing import List, Set, Dict, Tuple, Optional
from pathlib import Path

from pyCOT.ERC_Hierarchy import ERC, species_list_to_names


@dataclass
class OverlapInfo:
    """Store information about partial overlap between base1 and target generators."""
    generator_index: int
    generator_species: Set[str]
    overlapping_species: Set[str]
    base1_unique: Set[str]  # Species in base1 but not in generator
    generator_unique: Set[str]  # Species in generator but not in base1


@dataclass
class Base1SynergyInfo:
    """Store synergy information for a base1 ERC with a target."""
    base1_erc: 'ERC'
    hierarchy_level: int
    overlaps: List[OverlapInfo]
    overlap_signature: frozenset  # Unique signature for deduplication


class ERCExplorer:
    """
    Advanced ERC exploration and synergy analysis class.
    
    Provides methods for:
    - Target identification for potential synergies
    - Target-centric synergy organization with filtering
    - Hierarchy structural analysis
    - Comprehensive reporting and visualization
    """
    
    def __init__(self, reaction_network, ercs=None, hierarchy_graph=None):
        """
        Initialize ERC Explorer.
        
        Parameters
        ----------
        reaction_network : ReactionNetwork
            The reaction network
        ercs : list[ERC], optional
            List of ERCs. If None, will be computed automatically
        hierarchy_graph : nx.DiGraph, optional
            The ERC hierarchy graph. If None, will be built automatically
        """
        self.RN = reaction_network
        self.ercs = ercs if ercs is not None else ERC.ERCs(reaction_network)
        self.hierarchy = hierarchy_graph if hierarchy_graph is not None else ERC.build_hierarchy_graph(self.ercs, reaction_network)
        
        # Cache for analysis results
        self._potential_synergies = None
        self._target_centric_synergies = None
        self._hierarchy_analysis = None
        
    # =================== TARGET IDENTIFICATION ===================
    
    def identify_potential_synergy_targets(self, force_recompute=False):
        """
        Identify potential synergy targets using bottom-up base1 traversal 
        and top-down target search with recursive successor exploration.
        
        Parameters
        ----------
        force_recompute : bool
            Force recomputation even if cached results exist
        
        Returns
        -------
        dict
            Dictionary mapping base1 ERCs to lists of potential target ERCs
            {base1_erc: [target1, target2, ...], ...}
        """
        if self._potential_synergies is not None and not force_recompute:
            return self._potential_synergies
        
        # Get hierarchy levels - higher level means contains more ERCs
        levels = ERC.get_node_levels(self.hierarchy)
        
        # Sort ERCs by hierarchy level for base1 (bottom-up: ascending levels)
        base1_candidates = sorted(self.ercs, key=lambda e: levels[e.label])
        
        # Dictionary to store results: base1 -> list of potential targets
        potential_synergies = {}
        
        def has_strict_partial_overlap(base1_closure, generator_species):
            """Check if generator has strict partial overlap with base1 closure."""
            intersection = base1_closure & generator_species
            return (bool(intersection) and 
                    not generator_species.issubset(base1_closure) and 
                    not base1_closure.issubset(generator_species))
        
        def is_potential_target(target_erc, base1, base1_closure):
            """Check if target ERC is a potential synergy target for base1."""
            if target_erc == base1:
                return False
                
            for generator in target_erc.min_generators:
                generator_species = set(species_list_to_names(generator))
                if has_strict_partial_overlap(base1_closure, generator_species):
                    return True
            return False
        
        def explore_successors_recursively(current_erc, base1, base1_closure, found_targets, visited):
            """Recursively explore successors of a potential target."""
            # Get all successors (ERCs contained by current_erc)
            for successor_label in self.hierarchy.successors(current_erc.label):
                if successor_label in visited:
                    continue
                    
                successor_erc = next(erc for erc in self.ercs if erc.label == successor_label)
                visited.add(successor_label)
                
                # Check if successor is also a potential target
                if is_potential_target(successor_erc, base1, base1_closure):
                    found_targets.append(successor_erc)
                    # Recursively check successors of this successor
                    explore_successors_recursively(successor_erc, base1, base1_closure, found_targets, visited)
        
        # Main algorithm: for each base1 (bottom-up traversal)
        for base1 in base1_candidates:
            base1_closure = set(species_list_to_names(base1.get_closure(self.RN)))
            potential_targets = []
            
            # Get all ERCs as potential targets (excluding base1 and all ERCs contained by base1)
            base1_contained_labels = set(nx.descendants(self.hierarchy, base1.label))
            target_candidates = [erc for erc in self.ercs 
                               if erc != base1 and erc.label not in base1_contained_labels]
            
            # Sort target candidates by level (descending for top-down traversal)
            target_candidates.sort(key=lambda e: levels[e.label], reverse=True)
            
            # Track which ERCs we've already processed to avoid redundant checks
            visited = set()
            
            # Check each target candidate from top to bottom
            for target_candidate in target_candidates:
                if target_candidate.label in visited:
                    continue
                    
                visited.add(target_candidate.label)
                
                # Check if this candidate is a potential target
                if is_potential_target(target_candidate, base1, base1_closure):
                    potential_targets.append(target_candidate)
                    
                    # Recursively explore all successors of this potential target
                    explore_successors_recursively(target_candidate, base1, base1_closure, 
                                                 potential_targets, visited)
            
            # Store results if any potential targets were found
            if potential_targets:
                potential_synergies[base1] = potential_targets
        
        self._potential_synergies = potential_synergies
        return potential_synergies
    
    # =================== TARGET-CENTRIC ORGANIZATION ===================
    
    def compute_overlap_info(self, base1_erc, target_erc):
        """Compute detailed overlap information between base1 and target."""
        base1_closure = set(species_list_to_names(base1_erc.get_closure(self.RN)))
        overlaps = []
        
        for i, generator in enumerate(target_erc.min_generators):
            generator_species = set(species_list_to_names(generator))
            intersection = base1_closure & generator_species
            
            # Check for strict partial overlap
            if (intersection and 
                not generator_species.issubset(base1_closure) and 
                not base1_closure.issubset(generator_species)):
                
                overlap_info = OverlapInfo(
                    generator_index=i,
                    generator_species=generator_species,
                    overlapping_species=intersection,
                    base1_unique=base1_closure - generator_species,
                    generator_unique=generator_species - base1_closure
                )
                overlaps.append(overlap_info)
        
        return overlaps
    
    def create_overlap_signature(self, overlaps):
        """Create a unique signature for the overlap pattern."""
        signature_elements = []
        
        for overlap in overlaps:
            element = (
                overlap.generator_index,
                frozenset(overlap.overlapping_species)
            )
            signature_elements.append(element)
        
        return frozenset(signature_elements)
    
    def filter_redundant_base1s(self, base1_infos):
        """Filter out redundant base1 ERCs that provide identical overlap patterns."""
        # Group base1s by their overlap signature
        signature_groups = defaultdict(list)
        for base1_info in base1_infos:
            signature_groups[base1_info.overlap_signature].append(base1_info)
        
        filtered_base1s = []
        
        # Process each signature group
        for signature, group in signature_groups.items():
            if len(group) == 1:
                # Only one base1 with this signature, keep it
                filtered_base1s.extend(group)
            else:
                # Multiple base1s with same signature, keep only lowest hierarchy level
                min_level = min(info.hierarchy_level for info in group)
                min_level_base1s = [info for info in group if info.hierarchy_level == min_level]
                filtered_base1s.extend(min_level_base1s)
        
        return filtered_base1s
    
    def reorganize_to_target_centric(self, potential_synergies=None, force_recompute=False):
        """
        Reorganize potential synergies to target-centric structure with filtering.
        
        Parameters
        ----------
        potential_synergies : dict, optional
            Potential synergies. If None, will be computed automatically
        force_recompute : bool
            Force recomputation even if cached results exist
            
        Returns
        -------
        dict
            Target-centric structure: {target: [Base1SynergyInfo]}
        """
        if self._target_centric_synergies is not None and not force_recompute:
            return self._target_centric_synergies
        
        if potential_synergies is None:
            potential_synergies = self.identify_potential_synergy_targets(force_recompute)
        
        # Get hierarchy levels
        levels = ERC.get_node_levels(self.hierarchy)
        
        # First pass: collect all base1-target relationships with overlap info
        target_to_base1s = defaultdict(list)
        
        for base1_erc, target_list in potential_synergies.items():
            for target_erc in target_list:
                # Compute detailed overlap information
                overlaps = self.compute_overlap_info(base1_erc, target_erc)
                
                if overlaps:  # Only if there are actual overlaps
                    overlap_signature = self.create_overlap_signature(overlaps)
                    
                    base1_info = Base1SynergyInfo(
                        base1_erc=base1_erc,
                        hierarchy_level=levels[base1_erc.label],
                        overlaps=overlaps,
                        overlap_signature=overlap_signature
                    )
                    
                    target_to_base1s[target_erc].append(base1_info)
        
        # Second pass: filter redundant base1s for each target
        filtered_target_synergies = {}
        
        for target_erc, base1_infos in target_to_base1s.items():
            filtered_base1s = self.filter_redundant_base1s(base1_infos)
            if filtered_base1s:
                filtered_target_synergies[target_erc] = filtered_base1s
        
        self._target_centric_synergies = filtered_target_synergies
        return filtered_target_synergies
    
    # =================== HIERARCHY ANALYSIS ===================
    
    def analyze_hierarchy_structure(self, force_recompute=False):
        """Perform comprehensive structural analysis of the hierarchy."""
        if self._hierarchy_analysis is not None and not force_recompute:
            return self._hierarchy_analysis
        
        G = self.hierarchy
        levels = ERC.get_node_levels(G)
        
        analysis = {}
        
        # Basic metrics
        analysis['basic_metrics'] = {
            'total_nodes': G.number_of_nodes(),
            'total_edges': G.number_of_edges(),
            'max_level': max(levels.values()) if levels else 0,
            'min_level': min(levels.values()) if levels else 0,
            'hierarchy_depth': max(levels.values()) - min(levels.values()) + 1 if levels else 0
        }
        
        # Level distribution
        level_distribution = Counter(levels.values())
        analysis['level_distribution'] = dict(level_distribution)
        
        # Nodes per level
        nodes_per_level = defaultdict(list)
        for node, level in levels.items():
            nodes_per_level[level].append(node)
        analysis['nodes_per_level'] = dict(nodes_per_level)
        
        # Root and leaf analysis
        root_nodes = [node for node in G.nodes() if G.in_degree(node) == 0]
        leaf_nodes = [node for node in G.nodes() if G.out_degree(node) == 0]
        
        analysis['root_leaf_analysis'] = {
            'root_nodes': root_nodes,
            'leaf_nodes': leaf_nodes,
            'num_roots': len(root_nodes),
            'num_leaves': len(leaf_nodes)
        }
        
        self._hierarchy_analysis = analysis
        return analysis
    
    # =================== REPORTING AND VISUALIZATION ===================
    
    def print_potential_synergies(self, potential_synergies=None):
        """Print the identified potential synergies in a readable format."""
        if potential_synergies is None:
            potential_synergies = self.identify_potential_synergy_targets()
        
        if not potential_synergies:
            print("No potential synergies identified.")
            return
        
        print("=" * 80)
        print("POTENTIAL SYNERGY TARGETS IDENTIFIED")
        print("=" * 80)
        
        for base1, targets in potential_synergies.items():
            print(f"\nBASE1: {base1.label}")
            print(f"  Closure: {species_list_to_names(base1.get_closure(self.RN))}")
            print(f"  Potential Targets ({len(targets)}):")
            
            for target in targets:
                print(f"    • {target.label}")
                print(f"      Closure: {species_list_to_names(target.get_closure(self.RN))}")
                
                # Show which generators have strict partial overlap
                overlapping_gens = []
                base1_closure = set(species_list_to_names(base1.get_closure(self.RN)))
                
                for i, generator in enumerate(target.min_generators):
                    generator_species = set(species_list_to_names(generator))
                    intersection = base1_closure & generator_species
                    if (intersection and 
                        not generator_species.issubset(base1_closure) and 
                        not base1_closure.issubset(generator_species)):
                        overlapping_gens.append(f"Gen{i+1}:{species_list_to_names(generator)}")
                
                print(f"      Overlapping generators: {', '.join(overlapping_gens)}")
            
            print("-" * 60)
    
    def print_target_centric_synergies(self, target_synergies=None):
        """Print the target-centric synergy structure in a readable format."""
        if target_synergies is None:
            target_synergies = self.reorganize_to_target_centric()
        
        if not target_synergies:
            print("No target-centric synergies found.")
            return
        
        print("=" * 80)
        print("TARGET-CENTRIC SYNERGY ANALYSIS")
        print("=" * 80)
        
        for target_erc, base1_infos in target_synergies.items():
            print(f"\nTARGET: {target_erc.label}")
            print(f"  Closure: {species_list_to_names(target_erc.get_closure(self.RN))}")
            print(f"  Generators: {[species_list_to_names(gen) for gen in target_erc.min_generators]}")
            print(f"  Potential Base1 ERCs ({len(base1_infos)}):")
            
            for i, base1_info in enumerate(base1_infos, 1):
                print(f"\n    [{i}] BASE1: {base1_info.base1_erc.label} (Level {base1_info.hierarchy_level})")
                print(f"        Closure: {species_list_to_names(base1_info.base1_erc.get_closure(self.RN))}")
                print(f"        Overlapping Generators ({len(base1_info.overlaps)}):")
                
                for overlap in base1_info.overlaps:
                    print(f"          Generator {overlap.generator_index + 1}: {list(overlap.generator_species)}")
                    print(f"            Overlap: {list(overlap.overlapping_species)}")
                    print(f"            Base1 unique: {list(overlap.base1_unique) if overlap.base1_unique else 'None'}")
                    print(f"            Generator unique: {list(overlap.generator_unique) if overlap.generator_unique else 'None'}")
            
            print("-" * 60)
    
    def analyze_filtering_statistics(self, original_synergies=None, filtered_synergies=None):
        """Analyze and print statistics about the filtering process."""
        if original_synergies is None:
            original_synergies = self.identify_potential_synergy_targets()
        if filtered_synergies is None:
            filtered_synergies = self.reorganize_to_target_centric()
        
        # Count original base1-target pairs
        original_pairs = sum(len(targets) for targets in original_synergies.values())
        
        # Count filtered target-base1 pairs
        filtered_pairs = sum(len(base1_infos) for base1_infos in filtered_synergies.values())
        
        print("=" * 80)
        print("FILTERING STATISTICS")
        print("=" * 80)
        print(f"Original base1-target pairs: {original_pairs}")
        print(f"Filtered target-base1 pairs: {filtered_pairs}")
        print(f"Reduction: {original_pairs - filtered_pairs} pairs removed")
        if original_pairs > 0:
            print(f"Efficiency gain: {(1 - filtered_pairs/original_pairs)*100:.1f}% reduction")
        
        print(f"\nTargets with potential synergies: {len(filtered_synergies)}")
        
        # Show distribution of base1s per target
        base1_counts = [len(base1_infos) for base1_infos in filtered_synergies.values()]
        if base1_counts:
            print(f"Average base1s per target: {sum(base1_counts)/len(base1_counts):.1f}")
            print(f"Max base1s for a target: {max(base1_counts)}")
            print(f"Min base1s for a target: {min(base1_counts)}")
    
    def generate_comprehensive_report(self):
        """Generate a comprehensive exploration report."""
        print("\n" + "="*80)
        print("🎯 ERC EXPLORATION COMPREHENSIVE REPORT")
        print("="*80)
        
        # Basic network info
        print(f"\n📊 NETWORK OVERVIEW:")
        print(f"  Species: {len(self.RN.species())}")
        print(f"  Reactions: {len(self.RN.reactions())}")
        print(f"  ERCs: {len(self.ercs)}")
        
        # Hierarchy analysis
        print(f"\n🏗️  HIERARCHY ANALYSIS:")
        hierarchy_analysis = self.analyze_hierarchy_structure()
        basic = hierarchy_analysis['basic_metrics']
        print(f"  Hierarchy depth: {basic['hierarchy_depth']} levels")
        print(f"  Containment relations: {basic['total_edges']}")
        print(f"  Root ERCs: {hierarchy_analysis['root_leaf_analysis']['num_roots']}")
        print(f"  Leaf ERCs: {hierarchy_analysis['root_leaf_analysis']['num_leaves']}")
        
        # Target identification
        print(f"\n🎯 TARGET IDENTIFICATION:")
        potential_synergies = self.identify_potential_synergy_targets()
        print(f"  Base1 ERCs with targets: {len(potential_synergies)}")
        total_potential_pairs = sum(len(targets) for targets in potential_synergies.values())
        print(f"  Total potential synergy pairs: {total_potential_pairs}")
        
        # Target-centric organization
        print(f"\n🎛️  TARGET-CENTRIC ORGANIZATION:")
        target_synergies = self.reorganize_to_target_centric()
        print(f"  Targets with potential synergies: {len(target_synergies)}")
        
        # Filtering statistics
        self.analyze_filtering_statistics(potential_synergies, target_synergies)
        
        # Detailed results
        print(f"\n" + "-" * 60)
        self.print_target_centric_synergies(target_synergies)
        
        print("\n" + "="*80)
        print("✅ ERC EXPLORATION REPORT COMPLETE")
        print("="*80)
        
        return {
            'potential_synergies': potential_synergies,
            'target_synergies': target_synergies,
            'hierarchy_analysis': hierarchy_analysis
        }