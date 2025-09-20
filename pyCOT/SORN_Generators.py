#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SORN_Generators.py

Second Order Reaction Network and Irreducible Generator Construction

This module handles the construction of irreducible generators through productive 
novelty steps (synergy and complementarity). It provides efficient lookup of 
productive relationships between ERCs via the ERC_SORN class.

Classes:
    IrreducibleGenerator: Represents ERC sequences built through productive novelty
    ERC_SORN: Pre-computed network of productive relationships between ERCs
    ProductiveExtension: Helper class for extension possibilities
    
Functions:
    identify_p_ercs: Find persistent single ERCs
    build_erc_sorn: Factory function for ERC_SORN construction
    find_productive_extensions: Find productive extensions for generators
    build_irreducible_generators: Main generator construction function
    analyze_generator_statistics: Analysis utilities

Author: Based on theoretical work by Tomas Veloz et al.
"""

import time
from itertools import combinations
import networkx as nx
from collections import defaultdict, Counter
from pyCOT.ERC_Hierarchy import ERC, ERC_Hierarchy, closure, species_list_to_names
from pyCOT.ERC_Synergy_Complementarity import *

# ============================================================================
# HELPER CLASSES
# ============================================================================

class ProductiveExtension:
    """
    Represents a single productive extension possibility for a generator.
    
    This helper class encapsulates the information needed to extend a generator
    through productive novelty (synergy or complementarity).
    """
    
    def __init__(self, target_erc, step_type, step_details):
        """
        Initialize a productive extension.
        
        Parameters
        ----------
        target_erc : ERC
            The ERC that can be added to extend the generator
        step_type : str
            Type of productive step ('synergy' or 'complementarity')
        step_details : dict
            Details about the productive relationship
        """
        self.target_erc = target_erc
        self.step_type = step_type
        self.step_details = step_details
        
    def __repr__(self):
        return f"ProductiveExtension({self.target_erc.label}, {self.step_type})"

# ============================================================================
# CORE CLASSES
# ============================================================================

class IrreducibleGenerator:
    """
    Represents a sequence of ERCs that generates a closed set through 
    exclusively productive novelty steps (synergy and complementarity).
    
    This class tracks the construction path and provides methods to compute
    closures and check semi-self-maintenance properties.
    """
    
    def __init__(self, initial_erc_sequence=None):
        """
        Initialize an IrreducibleGenerator.
        
        Parameters
        ----------
        initial_erc_sequence : list of ERC, optional
            Initial sequence of ERCs
        """
        self.erc_sequence = initial_erc_sequence or []
        self.construction_path = []  # List of dicts with step details
        self.is_irreducible = True   # Built exclusively through productive novelty
        self._cached_closure = None
        self._cached_ssm = None
        
    def add_productive_step(self, new_erc, step_type, step_details):
        """
        Add an ERC through a productive novelty step.
        
        Parameters
        ----------
        new_erc : ERC
            ERC to add to the sequence
        step_type : str
            'synergy' or 'complementarity'
        step_details : dict
            Details about the productive step
        """
        self.erc_sequence.append(new_erc)
        self.construction_path.append({
            'erc': new_erc,
            'step_type': step_type,
            'details': step_details
        })
        # Clear caches when generator changes
        self._cached_closure = None
        self._cached_ssm = None
        
    def get_closure(self, RN):
        """
        Get closure using RN.generated_closure method with caching.
        
        Parameters
        ----------
        RN : ReactionNetwork
            The reaction network
            
        Returns
        -------
        list of Species
            The closure of all ERCs in the sequence
        """
        if self._cached_closure is None:
            if not self.erc_sequence:
                self._cached_closure = []
            else:
                # Collect all species from ERC closures
                all_species = []
                for erc in self.erc_sequence:
                    all_species.extend(erc.get_closure(RN))
                
                # Remove duplicates while preserving Species objects
                unique_species_dict = {sp.name: sp for sp in all_species}
                unique_species = list(unique_species_dict.values())
                
                # Use RN's closure method
                self._cached_closure = list(RN.generated_closure(unique_species))
            
        return self._cached_closure
    
    def check_ssm(self, RN):
        """
        Check if generator produces semi-self-maintaining set.
        
        Parameters
        ----------
        RN : ReactionNetwork
            The reaction network
            
        Returns
        -------
        bool
            True if the generated closure is semi-self-maintaining
        """
        if self._cached_ssm is None:
            closure_species = self.get_closure(RN)
            self._cached_ssm = is_semi_self_maintaining(RN, closure_species)
            
        return self._cached_ssm
    
    def size(self):
        """Get the number of ERCs in the generator."""
        return len(self.erc_sequence)
    
    def get_synergy_count(self):
        """Count synergistic steps in construction."""
        return sum(1 for step in self.construction_path if step['step_type'] == 'synergy')
    
    def get_complementarity_counts(self):
        """Count complementarity steps by type."""
        comp_counts = {'type1': 0, 'type2': 0, 'type3': 0}
        for step in self.construction_path:
            if step['step_type'] == 'complementarity':
                comp_type = step['details'].get('comp_type', 1)
                comp_counts[f'type{comp_type}'] += 1
        return comp_counts
    
    def get_erc_labels(self):
        """Get list of ERC labels in the generator."""
        return [erc.label for erc in self.erc_sequence]
    
    def __repr__(self):
        labels = self.get_erc_labels()
        return f"IrreducibleGenerator({labels}, size={self.size()})"

class ERC_SORN:
    """
    Second Order Reaction Network: Pre-computed network of relationships between ERCs.
    
    Stores all synergies and complementarities between ERC pairs for efficient lookup
    during generator construction, eliminating the need for repeated calculations.
    """
    
    def __init__(self, hierarchy, RN):
        """
        Initialize ERC_SORN by pre-computing all productive relationships.
        
        Parameters
        ----------
        hierarchy : ERC_Hierarchy
            The ERC hierarchy containing all ERCs
        RN : ReactionNetwork
            The reaction network
        """
        self.hierarchy = hierarchy
        self.RN = RN
        self.ercs = hierarchy.ercs

        # Core storage structures for O(1) lookup
        self._synergies = {}  # (a,b) -> list[ERC_Synergy]
        self._complementarities = {}  # (a,b) -> list[ERC_Complementarity]

        # Indexes for efficient queries
        self._erc_to_synergies = defaultdict(list)
        self._erc_to_complementarities = defaultdict(list)
        self._productive_partners = defaultdict(set)

        # Statistics
        self.computation_stats = {
            'total_pairs_checked': 0,
            'productive_pairs': 0,
            'total_synergies': 0,
            'total_complementarities': 0,
            'synergistic_pairs': 0,
            'complementary_pairs': 0,
            'build_time': 0.0
        }
        
        # Build the SORN
        self._build_sorn()
    
    def _build_sorn(self):
        """Build the SORN by computing all productive relationships."""
        start = time.time()
        print(f"Building ERC_SORN for {len(self.ercs)} ERCs...")

        # Compute all fundamental synergies
        all_synergies = []
        # for erc1 in self.ercs:
        #     for erc2 in self.ercs:
        #         if erc1.label >= erc2.label:
        #             continue
        #         synergies = get_fundamental_synergies_brute_force(erc1, erc2, self.hierarchy, self.RN)
        #         all_synergies.extend(synergies)
        synergies = get_all_fundamental_synergies_brute_force(self.ercs,self.hierarchy, self.RN, verbose=False)
        all_synergies.extend(synergies)

        # Map synergies by pair key for O(1) per-pair lookup
        synergies_by_pair = defaultdict(list)
        for syn in all_synergies:
            a_label = syn.rlabel[0]
            b_label = syn.rlabel[1]
            key = tuple(sorted((a_label, b_label)))
            synergies_by_pair[key].append(syn)

        # Store synergies
        for pair_key, syn_list in synergies_by_pair.items():
            a_label, b_label = pair_key
            self._store_synergies(a_label, b_label, syn_list)

        # Compute complementarities for all pairs
        for erc1, erc2 in combinations(self.ercs, 2):
            self.computation_stats['total_pairs_checked'] += 1

            if not can_interact(erc1, erc2, self.hierarchy):
                continue

            key = tuple(sorted((erc1.label, erc2.label)))
            pair_has_productive = False

            # Compute complementarities
            comps = get_complementarity(erc1, erc2, self.hierarchy, self.RN)
            if comps:
                self._store_complementarities(erc1.label, erc2.label, comps)
                pair_has_productive = True
                self.computation_stats['total_complementarities'] += len(comps)
                self.computation_stats['complementary_pairs'] += 1

            # Check if this pair had synergies
            syns = synergies_by_pair.get(key, [])
            if syns:
                pair_has_productive = True
                self.computation_stats['total_synergies'] += len(syns)
                self.computation_stats['synergistic_pairs'] += 1

            if pair_has_productive:
                self._productive_partners[erc1.label].add(erc2.label)
                self._productive_partners[erc2.label].add(erc1.label)
                self.computation_stats['productive_pairs'] += 1

        # Finalize stats
        self.computation_stats['build_time'] = time.time() - start

        print(f"ERC_SORN built: {self.computation_stats['productive_pairs']} productive pairs, "
              f"{self.computation_stats['total_synergies']} synergies, "
              f"{self.computation_stats['total_complementarities']} complementarities, "
              f"time={self.computation_stats['build_time']:.2f}s")
    
    def _store_synergies(self, erc1_label, erc2_label, synergies):
        """Store synergies with bidirectional lookup."""
        key1 = (erc1_label, erc2_label)
        key2 = (erc2_label, erc1_label)
        
        self._synergies[key1] = synergies
        self._synergies[key2] = synergies  # Bidirectional
        
        # Update indexes
        self._erc_to_synergies[erc1_label].append((erc2_label, synergies))
        self._erc_to_synergies[erc2_label].append((erc1_label, synergies))
    
    def _store_complementarities(self, erc1_label, erc2_label, complementarities):
        """Store complementarities with bidirectional lookup."""
        key1 = (erc1_label, erc2_label)
        key2 = (erc2_label, erc1_label)
        
        self._complementarities[key1] = complementarities
        self._complementarities[key2] = complementarities  # Bidirectional
        
        # Update indexes
        self._erc_to_complementarities[erc1_label].append((erc2_label, complementarities))
        self._erc_to_complementarities[erc2_label].append((erc1_label, complementarities))
    
    # ============================================================================
    # PUBLIC QUERY METHODS
    # ============================================================================
    
    def get_synergies(self, erc1_label, erc2_label):
        """Get all fundamental synergies between two ERCs."""
        return self._synergies.get((erc1_label, erc2_label), [])
    
    def get_complementarities(self, erc1_label, erc2_label):
        """Get all complementarities between two ERCs."""
        return self._complementarities.get((erc1_label, erc2_label), [])
    
    def has_synergy(self, erc1_label, erc2_label):
        """Check if two ERCs have any fundamental synergy."""
        return len(self.get_synergies(erc1_label, erc2_label)) > 0
    
    def has_complementarity(self, erc1_label, erc2_label):
        """Check if two ERCs have any complementarity."""
        return len(self.get_complementarities(erc1_label, erc2_label)) > 0
    
    def has_productive_relationship(self, erc1_label, erc2_label):
        """Check if two ERCs have any productive relationship."""
        return self.has_synergy(erc1_label, erc2_label) or self.has_complementarity(erc1_label, erc2_label)
    
    def get_productive_partners(self, erc_label):
        """Get all ERCs that have productive relationships with the given ERC."""
        return self._productive_partners.get(erc_label, set())
    
    def get_all_synergistic_partners(self, erc_label):
        """Get all ERCs that have synergistic relationships with the given ERC."""
        return self._erc_to_synergies.get(erc_label, [])
    
    def get_all_complementary_partners(self, erc_label):
        """Get all ERCs that have complementary relationships with the given ERC."""
        return self._erc_to_complementarities.get(erc_label, [])
    
    def get_productive_extensions_for_generator(self, erc_labels_in_generator):
        """
        Get all possible productive extensions for a generator (set of ERCs).
        
        Parameters
        ----------
        erc_labels_in_generator : list of str
            Labels of ERCs currently in the generator
            
        Returns
        -------
        list of ProductiveExtension
            List of productive extension possibilities
        """
        current_erc_set = set(erc_labels_in_generator)
        extensions = []
        candidate_extensions = {}  # candidate_label -> list of extension details
        
        # For each ERC in the current generator, find its productive partners
        for erc_label in erc_labels_in_generator:
            # Get synergistic partners
            for partner_label, synergies in self._erc_to_synergies.get(erc_label, []):
                if partner_label not in current_erc_set:
                    if partner_label not in candidate_extensions:
                        candidate_extensions[partner_label] = []
                    for synergy in synergies:
                        candidate_extensions[partner_label].append(('synergy', {
                            'synergy_type': 'fundamental',
                            'synergy_object': synergy,
                            'with_erc': erc_label
                        }))
            
            # Get complementary partners
            for partner_label, complementarities in self._erc_to_complementarities.get(erc_label, []):
                if partner_label not in current_erc_set:
                    if partner_label not in candidate_extensions:
                        candidate_extensions[partner_label] = []
                    for comp in complementarities:
                        candidate_extensions[partner_label].append(('complementarity', {
                            'comp_type': comp.comp_type,
                            'comp_object': comp,
                            'with_erc': erc_label
                        }))
        
        # Convert to ProductiveExtension objects
        label_to_erc = {erc.label: erc for erc in self.ercs}
        for candidate_label, extension_list in candidate_extensions.items():
            candidate_erc = label_to_erc[candidate_label]
            # For now, take the first extension found (could be enhanced)
            step_type, step_details = extension_list[0]
            extensions.append(ProductiveExtension(candidate_erc, step_type, step_details))
        
        return extensions
    
    def get_statistics(self):
        """Get statistics about the SORN."""
        return self.computation_stats.copy()

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def is_semi_self_maintaining(RN, species_set):
    """
    Check if a set of species is semi-self-maintaining.
    
    Parameters
    ----------
    RN : ReactionNetwork
        The reaction network
    species_set : list of Species
        Set of species to check
        
    Returns
    -------
    bool
        True if semi-self-maintaining (requirements = ∅)
    """
    # Get reactions from the species set
    reactions = RN.get_reactions_from_species(species_set)
    
    # Calculate consumed and produced species
    consumed = set()
    produced = set()
    
    for reaction in reactions:
        for edge in reaction.edges:
            if edge.type == "reactant":
                consumed.add(edge.species_name)
            elif edge.type == "product":
                produced.add(edge.species_name)
    
    # Semi-self-maintaining if no species are required (consumed but not produced)
    required = consumed - produced
    return len(required) == 0

# ============================================================================
# MAIN FUNCTIONS
# ============================================================================

def identify_p_ercs(hierarchy, RN):
    """
    Identify persistent ERCs (single ERCs that are semi-self-maintaining).
    
    Parameters
    ----------
    hierarchy : ERC_Hierarchy
        The ERC hierarchy
    RN : ReactionNetwork
        The reaction network
        
    Returns
    -------
    list of ERC
        List of persistent ERCs
    """
    p_ercs = []
    
    for erc in hierarchy.ercs:
        erc_closure = erc.get_closure(RN)
        if is_semi_self_maintaining(RN, erc_closure):
            p_ercs.append(erc)
    
    return p_ercs

def build_erc_sorn(hierarchy, RN):
    """
    Factory function to build an ERC_SORN.
    
    Parameters
    ----------
    hierarchy : ERC_Hierarchy
        The ERC hierarchy
    RN : ReactionNetwork
        The reaction network
        
    Returns
    -------
    ERC_SORN
        Pre-computed second-order reaction network
    """
    return ERC_SORN(hierarchy, RN)

def find_productive_extensions(generator, erc_sorn):
    """
    Find all ERCs that can extend a generator through productive novelty.
    
    Uses pre-computed ERC_SORN for O(1) lookup instead of expensive recalculation.
    
    Parameters
    ----------
    generator : IrreducibleGenerator
        Current generator to extend
    erc_sorn : ERC_SORN
        Pre-computed second-order reaction network
        
    Returns
    -------
    list of ProductiveExtension
        List of productive extension possibilities
    """
    # Get ERC labels in current generator
    current_erc_labels = generator.get_erc_labels()
    
    # Use ERC_SORN's method to find extensions
    return erc_sorn.get_productive_extensions_for_generator(current_erc_labels)

def build_irreducible_generators(hierarchy, RN, max_size=10, verbose=True):
    """
    Build irreducible generators maintaining only one per unique closure.
    
    This version dramatically reduces redundant computation by recognizing that multiple
    generators producing the same closure will explore identical productive extensions.
    
    Parameters
    ----------
    hierarchy : ERC_Hierarchy
        The ERC hierarchy
    RN : ReactionNetwork
        The reaction network
    max_size : int
        Maximum size of generators to build
    verbose : bool
        Whether to print detailed statistics
        
    Returns
    -------
    list of IrreducibleGenerator
        List of irreducible generators
    """
    if verbose:
        print("=" * 80)
        print("BUILDING IRREDUCIBLE GENERATORS")
        print("=" * 80)
    
    # Build ERC_SORN once at the beginning for efficient lookup
    erc_sorn = build_erc_sorn(hierarchy, RN)
    
    # Statistics tracking
    stats = {
        'total_generators_explored': 0,
        'unique_closures_found': 0,
        'redundant_generators_pruned': 0,
        'generators_by_size': {},
        'p_ercs': 0,
        'fundamental_synergies': 0
    }
    
    # Map from closure signature -> representative generator
    closure_to_generator = {}
    
    # Phase 1: Initialize with P-ERCs
    p_ercs = identify_p_ercs(hierarchy, RN)
    stats['p_ercs'] = len(p_ercs)
    
    if verbose:
        print(f"\nPhase 1: Found {len(p_ercs)} P-ERCs")
    
    for p_erc in p_ercs:
        gen = IrreducibleGenerator([p_erc])
        closure = gen.get_closure(RN)
        closure_sig = tuple(sorted(sp.name for sp in closure))
        closure_to_generator[closure_sig] = gen
        stats['total_generators_explored'] += 1
    
    # Phase 2: Initialize with fundamental synergies
    if verbose:
        print("\nPhase 2: Building fundamental synergy generators using ERC_SORN...")
    
    synergy_count = 0
    redundant_synergies = 0
    
    for erc1 in hierarchy.ercs:
        synergistic_partners = erc_sorn.get_all_synergistic_partners(erc1.label)
        
        for partner_label, synergies in synergistic_partners:
            # Avoid duplicate pairs
            if erc1.label >= partner_label:
                continue
                
            erc2 = next(erc for erc in hierarchy.ercs if erc.label == partner_label)
            
            for synergy in synergies:
                gen = IrreducibleGenerator([erc1, erc2])
                gen.add_productive_step(erc2, 'synergy', {
                    'synergy_type': 'fundamental',
                    'synergy_object': synergy,
                    'with_erc': erc1.label
                })
                
                # Check if this closure already exists
                closure = gen.get_closure(RN)
                closure_sig = tuple(sorted(sp.name for sp in closure))
                
                stats['total_generators_explored'] += 1
                
                if closure_sig not in closure_to_generator:
                    closure_to_generator[closure_sig] = gen
                    synergy_count += 1
                else:
                    redundant_synergies += 1
                    stats['redundant_generators_pruned'] += 1
    
    stats['fundamental_synergies'] = synergy_count
    
    if verbose:
        print(f"  - Found {synergy_count} unique fundamental synergy closures")
        print(f"  - Pruned {redundant_synergies} redundant synergy generators")
    
    # Phase 3: Iterative extension with closure-based deduplication
    if verbose:
        print("\nPhase 3: Extending generators through productive novelty...")
        print("         (maintaining one generator per unique closure)")
        print("-" * 80)
    
    current_size = 1 if p_ercs else 2
    
    while current_size < max_size:
        # Get current generators of this size
        current_gens = [(sig, gen) for sig, gen in closure_to_generator.items() 
                       if gen.size() == current_size]
        
        if not current_gens:
            break
        
        new_closure_to_generator = {}
        generators_explored_this_size = 0
        redundant_this_size = 0
        
        if verbose:
            print(f"\nSize {current_size}: Processing {len(current_gens)} unique closures")
        
        for closure_sig, gen in current_gens:
            # Find productive extensions using ERC_SORN
            extensions = find_productive_extensions(gen, erc_sorn)
            
            for extension in extensions:
                # Create new generator
                new_gen = IrreducibleGenerator(gen.erc_sequence.copy())
                new_gen.construction_path = gen.construction_path.copy()
                new_gen.add_productive_step(extension.target_erc, extension.step_type, extension.step_details)
                
                # Check if this creates a new closure
                new_closure = new_gen.get_closure(RN)
                new_closure_sig = tuple(sorted(sp.name for sp in new_closure))
                
                generators_explored_this_size += 1
                stats['total_generators_explored'] += 1
                
                # Only keep if this produces a NEW closure
                if (new_closure_sig not in closure_to_generator and 
                    new_closure_sig not in new_closure_to_generator):
                    new_closure_to_generator[new_closure_sig] = new_gen
                else:
                    redundant_this_size += 1
                    stats['redundant_generators_pruned'] += 1
        
        # Update main dictionary with new closures
        closure_to_generator.update(new_closure_to_generator)
        
        # Update statistics
        stats['generators_by_size'][current_size] = generators_explored_this_size
        
        if verbose:
            print(f"  └─ Size {current_size} → {current_size + 1}:")
            print(f"     - Explored: {generators_explored_this_size} potential generators")
            print(f"     - Found: {len(new_closure_to_generator)} NEW unique closures")
            print(f"     - Pruned: {redundant_this_size} redundant generators")
            if generators_explored_this_size > 0:
                reduction = 100 * redundant_this_size / generators_explored_this_size
                print(f"     - Reduction: {reduction:.1f}% redundancy eliminated")
        
        current_size += 1
    
    # Convert to list of generators
    all_generators = list(closure_to_generator.values())
    stats['unique_closures_found'] = len(all_generators)
    
    if verbose:
        print("\n" + "=" * 80)
        print("SUMMARY:")
        print(f"  Total generators explored: {stats['total_generators_explored']}")
        print(f"  Unique closures found: {stats['unique_closures_found']}")
        print(f"  Redundant generators pruned: {stats['redundant_generators_pruned']}")
        total_reduction = 100 * stats['redundant_generators_pruned'] / stats['total_generators_explored']
        print(f"  Overall reduction: {total_reduction:.1f}%")
        print("=" * 80)
    
    return all_generators

def analyze_generator_statistics(generators):
    """
    Analyze statistics of irreducible generators.
    
    Parameters
    ----------
    generators : list of IrreducibleGenerator
        Generators to analyze
        
    Returns
    -------
    dict
        Comprehensive statistics about the generators
    """
    if not generators:
        return {
            'total_generators': 0,
            'size_distribution': {},
            'construction_statistics': {},
            'error': 'No generators provided'
        }
    
    stats = {
        'total_generators': len(generators),
        'size_distribution': Counter(),
        'synergy_usage': [],
        'complementarity_usage': {'type1': [], 'type2': [], 'type3': []},
        'erc_usage': Counter(),
        'construction_paths': []
    }
    
    for gen in generators:
        # Size distribution
        stats['size_distribution'][gen.size()] += 1
        
        # Synergy usage
        stats['synergy_usage'].append(gen.get_synergy_count())
        
        # Complementarity usage
        comp_counts = gen.get_complementarity_counts()
        for comp_type in ['type1', 'type2', 'type3']:
            stats['complementarity_usage'][comp_type].append(comp_counts[comp_type])
        
        # ERC usage
        for erc_label in gen.get_erc_labels():
            stats['erc_usage'][erc_label] += 1
        
        # Construction paths
        stats['construction_paths'].append(len(gen.construction_path))
    
    # Compute summary statistics
    stats['summary'] = {
        'avg_size': sum(stats['synergy_usage']) / len(generators) if generators else 0,
        'max_size': max(stats['size_distribution'].keys()) if stats['size_distribution'] else 0,
        'most_used_ercs': stats['erc_usage'].most_common(5),
        'avg_synergies': sum(stats['synergy_usage']) / len(generators) if generators else 0
    }
    
    return stats