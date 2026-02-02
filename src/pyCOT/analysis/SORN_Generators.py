#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SORN_Generators.py

Irreducible Generator Construction

This module handles the construction of irreducible generators through productive 
novelty steps (synergy and complementarity). It provides efficient lookup of 
productive relationships between ERCs via the ERC_SORN class.

Classes:
    IrreducibleGenerator: Represents ERC sequences built through productive novelty
    ProductiveExtension: Helper class for extension possibilities
    
Functions:
    identify_p_ercs: Find persistent single ERCs
    find_productive_extensions: Find productive extensions for generators
    build_irreducible_generators: Main generator construction function
    analyze_generator_statistics: Analysis utilities

Author: Based on theoretical work by Tomas Veloz et al.
"""

import time
from itertools import combinations
import networkx as nx
from collections import defaultdict, Counter
from pyCOT.analysis.ERC_Hierarchy import ERC, ERC_Hierarchy, closure, species_list_to_names
from pyCOT.analysis.ERC_Synergy_Complementarity import (
    ERC_SORN,
    build_erc_sorn,
    get_fundamental_synergies_brute_force
)

# ============================================================================
# HELPER CLASSES
# ============================================================================

class ProductiveExtension:
    """
    Represents a single productive extension possibility for a generator.
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
        self.construction_path = []
        self.is_irreducible = True
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
                all_species = []
                for erc in self.erc_sequence:
                    all_species.extend(erc.get_closure(RN))
                
                unique_species_dict = {sp.name: sp for sp in all_species}
                unique_species = list(unique_species_dict.values())
                
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
    reactions = RN.get_reactions_from_species(species_set)
    
    consumed = set()
    produced = set()
    
    for reaction in reactions:
        for edge in reaction.edges:
            if edge.type == "reactant":
                consumed.add(edge.species_name)
            elif edge.type == "product":
                produced.add(edge.species_name)
    
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
    current_erc_labels = generator.get_erc_labels()
    
    extension_tuples = erc_sorn.get_productive_extensions_for_generator(current_erc_labels)
    
    extensions = []
    for candidate_label, step_type, step_details in extension_tuples:
        target_erc = next(erc for erc in erc_sorn.ercs if erc.label == candidate_label)
        extensions.append(ProductiveExtension(target_erc, step_type, step_details))
    
    return extensions

def build_irreducible_generators(hierarchy, RN, erc_sorn=None, max_size=10, verbose=True):
    """
    Build irreducible generators maintaining only one per unique closure.
    
    Parameters
    ----------
    hierarchy : ERC_Hierarchy
        The ERC hierarchy
    RN : ReactionNetwork
        The reaction network
    erc_sorn : ERC_SORN, optional
        Pre-computed SORN (will be built if not provided)
    max_size : int
        Maximum size of generators to build
    verbose : bool
        Whether to print detailed statistics
        
    Returns
    -------
    tuple
        (list of IrreducibleGenerator, ERC_SORN)
    """
    if verbose:
        print("=" * 80)
        print("BUILDING IRREDUCIBLE GENERATORS")
        print("=" * 80)
    
    if erc_sorn is None:
        erc_sorn = build_erc_sorn(hierarchy, RN)
    
    stats = {
        'total_generators_explored': 0,
        'unique_closures_found': 0,
        'redundant_generators_pruned': 0,
        'generators_by_size': {},
        'p_ercs': 0,
        'fundamental_synergies': 0
    }
    
    closure_to_generator = {}
    
    # Phase 1: P-ERCs
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
    
    # Phase 2: Fundamental synergies
    if verbose:
        print("\nPhase 2: Building fundamental synergy generators using ERC_SORN...")
    
    synergy_count = 0
    redundant_synergies = 0
    
    for erc1 in hierarchy.ercs:
        synergistic_partners = erc_sorn.get_all_synergistic_partners(erc1.label)
        
        for partner_label, synergies in synergistic_partners:
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
    
    # Phase 3: Iterative extension
    if verbose:
        print("\nPhase 3: Extending generators through productive novelty...")
        print("         (maintaining one generator per unique closure)")
        print("-" * 80)
    
    current_size = 1 if p_ercs else 2
    
    while current_size < max_size:
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
            extensions = find_productive_extensions(gen, erc_sorn)
            
            for extension in extensions:
                new_gen = IrreducibleGenerator(gen.erc_sequence.copy())
                new_gen.construction_path = gen.construction_path.copy()
                new_gen.add_productive_step(extension.target_erc, extension.step_type, extension.step_details)
                
                new_closure = new_gen.get_closure(RN)
                new_closure_sig = tuple(sorted(sp.name for sp in new_closure))
                
                generators_explored_this_size += 1
                stats['total_generators_explored'] += 1
                
                if (new_closure_sig not in closure_to_generator and 
                    new_closure_sig not in new_closure_to_generator):
                    new_closure_to_generator[new_closure_sig] = new_gen
                else:
                    redundant_this_size += 1
                    stats['redundant_generators_pruned'] += 1
        
        closure_to_generator.update(new_closure_to_generator)
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
    
    return all_generators, erc_sorn

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
        stats['size_distribution'][gen.size()] += 1
        stats['synergy_usage'].append(gen.get_synergy_count())
        
        comp_counts = gen.get_complementarity_counts()
        for comp_type in ['type1', 'type2', 'type3']:
            stats['complementarity_usage'][comp_type].append(comp_counts[comp_type])
        
        for erc_label in gen.get_erc_labels():
            stats['erc_usage'][erc_label] += 1
        
        stats['construction_paths'].append(len(gen.construction_path))
    
    stats['summary'] = {
        'avg_size': sum(stats['synergy_usage']) / len(generators) if generators else 0,
        'max_size': max(stats['size_distribution'].keys()) if stats['size_distribution'] else 0,
        'most_used_ercs': stats['erc_usage'].most_common(5),
        'avg_synergies': sum(stats['synergy_usage']) / len(generators) if generators else 0
    }
    
    return stats