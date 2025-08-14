#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Persistent Modules Library

Implementation of relevant generators theorem for computing elementary 
semi-organizations through irreducible ERC sequences exhibiting synergy 
and complementarity.

Based on theoretical framework by Tomas Veloz et al:
"Synergy and Complementarity: The Generative Basis of Chemical Organizations"

Classes:
    IrreducibleGenerator: Represents ERC sequences built through productive novelty
    ElementarySO: Represents elementary semi-organizations
    
Functions:
    identify_p_ercs: Find persistent single ERCs
    build_irreducible_generators: Construct all irreducible generators
    compute_elementary_sos: Generate elementary semi-organizations
    analyze_generator_statistics: Compute analysis metrics

Author: Based on theoretical work by Tomas Veloz et al.
"""

from itertools import combinations
from collections import defaultdict, Counter
from pyCOT.ERC_Hierarchy import ERC, ERC_Hierarchy, closure, species_list_to_names
from pyCOT.ERC_Synergy_Complementarity import (
    get_fundamental_synergies, get_complementarity, 
    ERC_Synergy, ERC_Complementarity, build_erc_sorn, ERC_SORN
)

class IrreducibleGenerator:
    """
    Represents a sequence of ERCs that generates a closed set through 
    exclusively productive novelty steps (synergy and complementarity).
    """
    
    def __init__(self, initial_erc_sequence=None, RN=None):
        """
        Initialize an IrreducibleGenerator.
        
        Parameters
        ----------
        initial_erc_sequence : list of ERC, optional
            Initial sequence of ERCs
        RN : ReactionNetwork, optional
            The reaction network
        """
        self.erc_sequence = initial_erc_sequence or []
        self.construction_path = []  # List of dicts with step details
        self.is_irreducible = True   # Built exclusively through productive novelty
        self._cached_closure = None
        self._cached_ssm = None
        self._RN_reference = None
        
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
        # Clear caches
        self._cached_closure = None
        self._cached_ssm = None
        
    def get_closure(self, RN):
        """
        Get closure using existing RN.generated_closure method with caching.
        
        Parameters
        ----------
        RN : ReactionNetwork
            The reaction network
            
        Returns
        -------
        list of Species
            The closure of all ERCs in the sequence
        """
        if self._cached_closure is None or self._RN_reference != id(RN):
            if not self.erc_sequence:
                self._cached_closure = []
            else:
                # Collect all species from ERC closures (using cached ERC closures)
                all_species = []
                for erc in self.erc_sequence:
                    all_species.extend(erc.get_closure(RN))
                
                # Remove duplicates while preserving Species objects
                unique_species_dict = {sp.name: sp for sp in all_species}
                unique_species = list(unique_species_dict.values())
                
                # Use existing RN.generated_closure method
                self._cached_closure = list(RN.generated_closure(unique_species))
            
            self._RN_reference = id(RN)
            self._cached_ssm = None  # SSM depends on closure
            
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
        if self._cached_ssm is None or self._RN_reference != id(RN):
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
    
    def __repr__(self):
        labels = [erc.label for erc in self.erc_sequence]
        return f"IrreducibleGenerator({labels}, irreducible={self.is_irreducible})"

class ElementarySO:
    """
    Represents an elementary semi-organization (persistent module 
    containing no smaller persistent modules).
    """
    
    def __init__(self, closure_species, constituent_ercs, RN):
        """
        Initialize an ElementarySO.
        
        Parameters
        ----------
        closure_species : list of Species
            Species in the elementary SO
        constituent_ercs : list of ERC
            ERCs contained in this elementary SO
        RN : ReactionNetwork
            The reaction network
        """
        self.closure_species = closure_species
        self.constituent_ercs = constituent_ercs
        self.generating_pathways = []  # List of IrreducibleGenerators
        # is_p_erc is true if it's a single ERC and that ERC is SSM
        # Corrected: Use the global is_semi_self_maintaining function
        self.is_p_erc = len(constituent_ercs) == 1 and \
                        is_semi_self_maintaining(RN, constituent_ercs[0].get_closure(RN))
        
        # Cache species names for easy comparison
        self.species_names = set(sp.name for sp in closure_species)
        
        # Statistics
        self._synergy_stats = None
        self._complementarity_stats = None
        self._size_stats = None
        
    def add_generator(self, generator):
        """Add an irreducible generator that produces this SO."""
        self.generating_pathways.append(generator)
        # Clear cached statistics
        self._synergy_stats = None
        self._complementarity_stats = None
        self._size_stats = None
    
    def get_generation_statistics(self):
        """
        Analyze the productive novelty patterns used to generate this SO.
        
        Returns
        -------
        dict
            Statistics about synergy and complementarity usage
        """
        if not self.generating_pathways:
            return {
                'synergy_stats': {
                    'mean': float('nan'),
                    'max': float('nan'),
                    'min': float('nan'),
                    'counts': []
                },
                'complementarity_stats': {
                    'type1': {'mean': float('nan'), 'counts': []},
                    'type2': {'mean': float('nan'), 'counts': []},
                    'type3': {'mean': float('nan'), 'counts': []}
                },
                'size_stats': {
                    'generator_sizes': [],
                    'closure_size': len(self.closure_species),
                    'size_ratios': []
                },
                'pathway_count': 0,
                'is_p_erc': self.is_p_erc
            }
            
        if self._synergy_stats is None:
            self._compute_statistics()
            
        return {
            'synergy_stats': self._synergy_stats,
            'complementarity_stats': self._complementarity_stats,
            'size_stats': self._size_stats,
            'pathway_count': len(self.generating_pathways),
            'is_p_erc': self.is_p_erc
        }
    
    def _compute_statistics(self):
        """Compute cached statistics."""
        synergy_counts = []
        comp_type1_counts = []
        comp_type2_counts = []
        comp_type3_counts = []
        generator_sizes = []
        
        for gen in self.generating_pathways:
            synergy_counts.append(gen.get_synergy_count())
            comp_counts = gen.get_complementarity_counts()
            comp_type1_counts.append(comp_counts['type1'])
            comp_type2_counts.append(comp_counts['type2'])
            comp_type3_counts.append(comp_counts['type3'])
            generator_sizes.append(gen.size())
        
        # Handle empty lists with NaN values
        self._synergy_stats = {
            'mean': float('nan') if not synergy_counts else sum(synergy_counts) / len(synergy_counts),
            'max': float('nan') if not synergy_counts else max(synergy_counts),
            'min': float('nan') if not synergy_counts else min(synergy_counts),
            'counts': synergy_counts
        }
        
        self._complementarity_stats = {
            'type1': {
                'mean': float('nan') if not comp_type1_counts else sum(comp_type1_counts) / len(comp_type1_counts),
                'counts': comp_type1_counts
            },
            'type2': {
                'mean': float('nan') if not comp_type2_counts else sum(comp_type2_counts) / len(comp_type2_counts),
                'counts': comp_type2_counts
            },
            'type3': {
                'mean': float('nan') if not comp_type3_counts else sum(comp_type3_counts) / len(comp_type3_counts),
                'counts': comp_type3_counts
            }
        }
        
        closure_size = len(self.closure_species)
        self._size_stats = {
            'generator_sizes': generator_sizes,
            'closure_size': closure_size,
            'size_ratios': (
                [gen_size / closure_size for gen_size in generator_sizes]
                if closure_size > 0 and generator_sizes 
                else []
            )
        }
    
    def __repr__(self):
        return f"ElementarySO(species={len(self.closure_species)}, ERCs={len(self.constituent_ercs)}, P-ERC={self.is_p_erc})"

# ============================================================================
# CORE FUNCTIONS
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
        True if semi-self-maintaining (reqs = ∅)
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
    OPTIMIZED VERSION: Find all ERCs that can extend a generator through productive novelty.
    
    Uses pre-computed ERC_SORN for O(1) lookup instead of expensive recalculation.
    
    Parameters
    ----------
    generator : IrreducibleGenerator
        Current generator to extend
    erc_sorn : ERC_SORN
        Pre-computed second-order reaction network
        
    Returns
    -------
    list of tuple
        List of (erc_object, step_type, step_details) for productive extensions
    """
    # Get ERC labels in current generator
    current_erc_labels = [erc.label for erc in generator.erc_sequence]
    
    # Use ERC_SORN's optimized method to find extensions
    extensions_data = erc_sorn.get_productive_extensions_for_generator(current_erc_labels)
    
    # Convert labels back to ERC objects and format for return
    extensions = []
    label_to_erc = {erc.label: erc for erc in erc_sorn.ercs}
    
    for candidate_label, step_type, step_details in extensions_data:
        candidate_erc = label_to_erc[candidate_label]
        extensions.append((candidate_erc, step_type, step_details))
    
    return extensions

def build_irreducible_generators(hierarchy, RN, max_size=10, verbose=True):
    """
    OPTIMIZED VERSION: Build irreducible generators maintaining only one per unique closure.
    
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
    tuple
        (generators, stats) where generators is list of IrreducibleGenerator
        and stats is a dict with reduction statistics
    """
    if verbose:
        print("=" * 80)
        print("Building ERC_SORN (pre-computing all productive relationships)...")
    
    # Build ERC_SORN once at the beginning
    erc_sorn = build_erc_sorn(hierarchy, RN)
    
    # Statistics tracking
    stats = {
        'total_generators_explored': 0,
        'unique_closures_found': 0,
        'redundant_generators_pruned': 0,
        'generators_by_size': {},
        'closures_by_size': {},
        'reduction_by_size': {},
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
        gen = IrreducibleGenerator([p_erc], RN)
        closure = gen.get_closure(RN)
        closure_sig = tuple(sorted(sp.name for sp in closure))
        closure_to_generator[closure_sig] = gen
        stats['total_generators_explored'] += 1
    
    # Phase 2: Initialize with fundamental synergies (OPTIMIZED)
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
                gen = IrreducibleGenerator([erc1, erc2], RN)
                gen.add_productive_step(erc2, 'synergy', {
                    'synergy_type': 'fundamental',
                    'synergy_object': synergy,
                    'with_erc': erc1.label
                })
                
                # Check if this closure already exists
                closure = gen.get_closure(RN)
                closure_sig = tuple(sorted(sp.name for sp in closure))
                
                #SHOULD CACHE THESE CLOSURES??



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
            
            for new_erc, step_type, step_details in extensions:
                # Create new generator
                new_gen = IrreducibleGenerator(gen.erc_sequence.copy(), RN)
                new_gen.construction_path = gen.construction_path.copy()
                new_gen.add_productive_step(new_erc, step_type, step_details)
                
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
        stats['closures_by_size'][current_size] = len(new_closure_to_generator)
        stats['reduction_by_size'][current_size] = {
            'explored': generators_explored_this_size,
            'unique': len(new_closure_to_generator),
            'redundant': redundant_this_size,
            'reduction_ratio': redundant_this_size / generators_explored_this_size if generators_explored_this_size > 0 else 0
        }
        
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
    
    return all_generators, stats

def compute_elementary_sos(generators, RN, hierarchy=None):
    """
    Compute elementary semi-organizations from irreducible generators.
    
    Based on the theoretical framework, the irreducible generators, when they
    produce a semi-organization, inherently produce an elementary one.
    Therefore, no explicit post-filtering for elementarity is needed.
    
    Parameters
    ----------
    generators : list of IrreducibleGenerator
        All irreducible generators
    RN : ReactionNetwork
        The reaction network
        
    Returns
    -------
    list of ElementarySO
        Elementary semi-organizations found
    """
    print("Computing elementary semi-organizations (relying on theoretical construction for elementarity)...")
    if hierarchy is None:
        print("WARNING: No ERC hierarchy provided. Creating new ERC_Hierarchy (closures will be (re?)computed)")
        hierarchy = ERC_Hierarchy(RN)  # This recreates ERCs and loses caches!
    else:
        print("✓ Using existing hierarchy (preserving cached closures)")
    # Filter for SSM generators: These are all semi-organizations (closed and self-maintaining)
    ssm_generators = []
    for gen in generators:
        if gen.check_ssm(RN):
            ssm_generators.append(gen)
    
    print(f"Found {len(ssm_generators)} SSM generators.")
    
    # Group generators by their closure (same closure = same SO)
    closure_groups = defaultdict(list)
    for gen in ssm_generators:
        closure_species = gen.get_closure(RN)
        closure_signature = tuple(sorted(sp.name for sp in closure_species))
        closure_groups[closure_signature].append(gen)
    
    print(f"Found {len(closure_groups)} unique closures (potential semi-organizations).")
    
    # Create ElementarySO objects for ALL identified semi-organizations.
    # Based on the theoretical guarantees, these are already the elementary ones.
    elementary_sos = []
    hierarchy = ERC_Hierarchy(RN) # Re-initialize hierarchy for consistency if needed
    
    for closure_signature, gens in closure_groups.items():
        representative_gen = gens[0]
        closure_species = representative_gen.get_closure(RN)
        
        # Collect all unique constituent ERCs for this semi-organization
        constituent_ercs = []
        unique_ercs_in_closure = set() 
        for gen in gens:
            for erc in gen.erc_sequence:
                if erc not in unique_ercs_in_closure:
                    unique_ercs_in_closure.add(erc)
        
        closure_names_set = set(sp.name for sp in closure_species)
        for erc in unique_ercs_in_closure:
            erc_closure_names = set(sp.name for sp in erc.get_closure(RN))
            if erc_closure_names.issubset(closure_names_set):
                constituent_ercs.append(erc)
        
        elementary_so = ElementarySO(closure_species, constituent_ercs, RN)
        for gen in gens:
            elementary_so.add_generator(gen)
        elementary_sos.append(elementary_so)
    
    print(f"Identified {len(elementary_sos)} elementary semi-organizations (by construction).")
    return elementary_sos

def analyze_generator_statistics(elementary_sos):
    """
    Analyze statistics of generators and elementary SOs.
    
    Parameters
    ----------
    elementary_sos : list of ElementarySO
        Elementary semi-organizations to analyze
        
    Returns
    -------
    dict
        Comprehensive statistics about the generation process
    """
    stats = {
        'total_elementary_sos': len(elementary_sos),
        'p_erc_count': sum(1 for so in elementary_sos if so.is_p_erc),
        'multi_erc_count': sum(1 for so in elementary_sos if not so.is_p_erc),
        'generator_counts': [],
        'synergy_usage': [],
        'complementarity_usage': {'type1': [], 'type2': [], 'type3': []},
        'size_ratios': [],
        'erc_distribution': Counter(),
        'pathway_diversity': []
    }
    
    for so in elementary_sos:
        so_stats = so.get_generation_statistics()
        
        # Generator count per SO
        stats['generator_counts'].append(so_stats.get('pathway_count', 0))
        
        # Synergy usage
        if 'synergy_stats' in so_stats:
            stats['synergy_usage'].extend(so_stats['synergy_stats'].get('counts', []))
        
        # Complementarity usage
        if 'complementarity_stats' in so_stats:
            for comp_type in ['type1', 'type2', 'type3']:
                if comp_type in so_stats['complementarity_stats']:
                    stats['complementarity_usage'][comp_type].extend(
                        so_stats['complementarity_stats'][comp_type].get('counts', [])
                    )
        
        # Size ratios
        if 'size_stats' in so_stats:
            stats['size_ratios'].extend(so_stats['size_stats'].get('size_ratios', []))
        
        # ERC distribution
        for erc in so.constituent_ercs:
            stats['erc_distribution'][erc.label] += 1
        
        # Pathway diversity
        stats['pathway_diversity'].append(len(so.generating_pathways))
    
    # Compute summary statistics
    stats['summary'] = {
        'avg_generators_per_so': sum(stats['generator_counts']) / len(stats['generator_counts']) if stats['generator_counts'] else 0,
        'avg_synergies_per_path': sum(stats['synergy_usage']) / len(stats['synergy_usage']) if stats['synergy_usage'] else 0,
        'avg_size_ratio': sum(stats['size_ratios']) / len(stats['size_ratios']) if stats['size_ratios'] else 0,
        'most_common_ercs': stats['erc_distribution'].most_common(5),
        'total_pathways': sum(stats['generator_counts'])
    }
    
    return stats

# ============================================================================
# MAIN COMPUTATION FUNCTION
# ============================================================================

def compute_persistent_modules(RN, max_generator_size=8, use_optimized=True):
    """
    OPTIMIZED VERSION: Compute all elementary persistent modules with closure-based deduplication.
    
    Parameters
    ----------
    RN : ReactionNetwork
        The reaction network to analyze
    max_generator_size : int
        Maximum size of generators to build
    use_optimized : bool
        Whether to use the optimized algorithm (True) or original (False) for comparison
        
    Returns
    -------
    tuple
        (elementary_sos, generators, statistics, erc_sorn)
    """
    print("=" * 80)
    if use_optimized:
        print("Computing Persistent Modules (OPTIMIZED with closure deduplication)")
    else:
        print("Computing Persistent Modules (ORIGINAL algorithm)")
    print("=" * 80)
    
    # Build ERC hierarchy
    print("\nBuilding ERC hierarchy...")
    hierarchy = ERC_Hierarchy(RN)
    hierarchy.build_hierarchy_graph()
    
    print(f"Network statistics:")
    print(f"  - Species: {len(RN.species())}")
    print(f"  - Reactions: {len(RN.reactions())}")  
    print(f"  - ERCs: {len(hierarchy.ercs)}")
    
    # Build irreducible generators
    print(f"\nBuilding irreducible generators (max size: {max_generator_size})...")
    
    if use_optimized:
        generators, build_stats = build_irreducible_generators(
            hierarchy, RN, max_generator_size
        )
    else:
        # Placeholder for original algorithm if needed for comparison
        # For this fix, we are only concerned with the optimized path
        raise NotImplementedError("Original algorithm not implemented for this comparison. Please use use_optimized=True.")
    
    # Compute elementary SOs
    # This function now includes the P-ERC filtering logic
    print("\nComputing elementary semi-organizations...")
    elementary_sos = compute_elementary_sos(generators, RN, hierarchy)
    
    # Analyze statistics
    print("\nAnalyzing statistics...")
    analysis_stats = analyze_generator_statistics(elementary_sos)
    
    # Combine all statistics
    statistics = {
        'build_stats': build_stats,
        'analysis_stats': analysis_stats,
        'algorithm': 'optimized' if use_optimized else 'original',
        'summary': {
            'total_ercs': len(hierarchy.ercs),
            'generators_explored': build_stats.get('total_generators_explored', 
                                                  build_stats.get('total_generators_built', 0)),
            'unique_closures': build_stats['unique_closures_found'],
            'redundancy_eliminated': build_stats.get('redundant_generators_pruned', 0),
            'elementary_sos': len(elementary_sos),
            'reduction_percentage': 0  # Will calculate below
        }
    }
    
    # Calculate reduction percentage
    if use_optimized and 'redundant_generators_pruned' in build_stats:
        total = build_stats['total_generators_explored']
        pruned = build_stats['redundant_generators_pruned']
        statistics['summary']['reduction_percentage'] = 100 * pruned / total if total > 0 else 0
    
    # Build ERC_SORN for additional analysis
    erc_sorn = build_erc_sorn(hierarchy, RN)
    
    # Print summary
    print("\n" + "=" * 80)
    print("COMPUTATION COMPLETE")
    print("=" * 80)
    print(f"Algorithm: {statistics['algorithm'].upper()}")
    print(f"ERCs: {statistics['summary']['total_ercs']}")
    print(f"Generators explored: {statistics['summary']['generators_explored']}")
    print(f"Unique closures found: {statistics['summary']['unique_closures']}")
    if use_optimized:
        print(f"Redundant generators eliminated: {statistics['summary']['redundancy_eliminated']}")
        print(f"Reduction: {statistics['summary']['reduction_percentage']:.1f}%")
    print(f"Elementary semi-organizations: {statistics['summary']['elementary_sos']}")
    print("=" * 80)
    
    return elementary_sos, generators, statistics, erc_sorn

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Organization Detection Extensions for Persistent_Modules.py

Functions to verify which elementary semi-organizations are organizations
by checking the self-maintaining condition using linear programming.
"""

import numpy as np
from scipy.optimize import linprog
from collections import defaultdict


def minimize_sv(S, epsilon=0.01, method='highs'):
    """
    Solve the linear programming problem to find self-maintaining vector.
    
    We want to find v such that:
    - S @ v >= 0 (net production >= 0 for each species)
    - v >= epsilon for all reactions (all reactions must be active)
    
    Parameters:
    ----------
    S : numpy.ndarray or StoichiometryMatrix
        Stoichiometric matrix (species x reactions)
    epsilon : float
        Minimum value for active reactions
    method : str
        LP solver method
        
    Returns:
    -------
    list: [is_feasible (bool), solution_vector (np.array or None)]
    """
    # Convert to numpy array if needed
    if hasattr(S, '__array__'):
        S_array = np.asarray(S)
    else:
        S_array = S
    
    n_species, n_reactions = S_array.shape
    
    if n_reactions == 0:
        return [False, None]
    
    # Objective: minimize sum of reaction rates (find feasible solution)
    c = np.ones(n_reactions)
    
    # Inequality constraints: -S @ v <= 0  (i.e., S @ v >= 0)
    A_ub = -S_array
    b_ub = np.zeros(n_species)
    
    # Bounds: v >= epsilon for all reactions (all must be active)
    bounds = [(epsilon, None) for _ in range(n_reactions)]
    
    try:
        # Solve linear program
        result = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method=method)
        
        if result.success:
            return [True, result.x]
        else:
            return [False, None]
    except Exception as e:
        print(f"Linear programming error: {e}")
        return [False, None]


def is_self_maintaining(RN, X,epsilon=1e-6):
    """
    Verifies if X is a self-maintaining set w.r.t RN using linear programming.
    
    A set is self-maintaining if there exists a flux vector v such that:
    - All reactions in the closure of X are active (v[i] > 0)
    - The net production of each species is non-negative (S @ v >= 0)
    
    Parameters:
    ----------
    RN : ReactionNetwork
        The reaction network
    X : list of Species objects or species names
        Set of species to check
        
    Returns:
    -------
    list: [Boolean, np.array or None, np.array or None] 
        Boolean: True if self-maintaining
        np.array1: self-maintaining flux vector (if exists)
        np.array2: production vector S @ v (if exists)
    """
    try:
        # Handle different input types
        if not X:
            return [False, None, None]
            
        # Convert to species objects if needed
        if isinstance(X[0], str):
            species_objects = [RN.get_species(name) for name in X]
        else:
            species_objects = X
        
        # Create sub-reaction network with species and their activated reactions
        sub_RN = RN.sub_reaction_network(species_objects)
        
        # Get stoichiometric matrix
        S = sub_RN.stoichiometry_matrix()
        
        # Check if there are any reactions
        if S.shape[1] == 0:
            # No reactions means not self-maintaining (unless it's empty set)
            return [False, None, None]
        
        # Solve for self-maintaining vector
        res = minimize_sv(S, epsilon=epsilon, method='highs')
        print(f"Self-maintaining check for {[sp.name for sp in species_objects]}: {res[0]}")
        if res[0]:  # If successful
            flux_vector = res[1]
            production_vector = np.asarray(S) @ flux_vector
            return [True, flux_vector, production_vector]
        else:
            return [False, None, None]
            
    except Exception as e:
        print(f"Error in is_self_maintaining: {e}")
        return [False, None, None]


def check_elementary_sos_are_organizations(elementary_sos, RN, verbose=True):
    """
    Check which elementary semi-organizations are actually organizations.
    
    Parameters:
    ----------
    elementary_sos : list of ElementarySO
        List of elementary semi-organizations to check
    RN : ReactionNetwork
        The reaction network
    verbose : bool
        Whether to print detailed information
        
    Returns:
    -------
    tuple: (organizations, flux_vectors, production_vectors)
        organizations: list of ElementarySO that are self-maintaining
        flux_vectors: corresponding flux vectors
        production_vectors: corresponding production vectors
    """
    if verbose:
        print("="*70)
        print("CHECKING ELEMENTARY SEMI-ORGANIZATIONS FOR SELF-MAINTENANCE")
        print("="*70)
        print(f"Total elementary SOs to check: {len(elementary_sos)}")
    
    organizations = []
    flux_vectors = []
    production_vectors = []
    
    for i, elem_so in enumerate(elementary_sos):
        if verbose:
            print(f"\nElementary SO {i+1}:")
            print(f"  Species count: {len(elem_so.closure_species)}")
            print(f"  ERCs count: {len(elem_so.constituent_ercs)}")
            print(f"  Is P-ERC: {elem_so.is_p_erc}")
            if len(elem_so.closure_species) <= 20:  # Only show species if not too many
                species_names = [sp.name for sp in elem_so.closure_species]
                print(f"  Species: {species_names}")
        
        # Check self-maintenance
        result = is_self_maintaining(RN, elem_so.closure_species)
        
        if result[0]:  # Is self-maintaining
            if verbose:
                print(f"  ✓ ORGANIZATION: Self-maintaining")
            organizations.append(elem_so)
            flux_vectors.append(result[1])
            production_vectors.append(result[2])
        else:
            if verbose:
                print(f"  ✗ Not self-maintaining")
    
    if verbose:
        print("\n" + "="*70)
        print("RESULTS SUMMARY")
        print("="*70)
        print(f"Elementary SOs checked: {len(elementary_sos)}")
        print(f"Organizations found: {len(organizations)}")
        print(f"Conversion rate: {len(organizations)/len(elementary_sos)*100:.1f}%")
        
        # Analyze by type
        p_erc_count = sum(1 for so in elementary_sos if so.is_p_erc)
        p_erc_orgs = sum(1 for org in organizations if org.is_p_erc)
        print(f"P-ERCs: {p_erc_orgs}/{p_erc_count} are organizations")
        
        multi_erc_count = len(elementary_sos) - p_erc_count
        multi_erc_orgs = len(organizations) - p_erc_orgs
        print(f"Multi-ERC SOs: {multi_erc_orgs}/{multi_erc_count} are organizations")
    
    return organizations, flux_vectors, production_vectors


def compute_elementary_organizations(RN, max_generator_size=8, verbose=True):
    """
    Compute all elementary organizations by finding elementary SOs and checking self-maintenance.
    
    This is the main function that combines the irreducible generator approach
    with linear programming verification.
    
    Parameters:
    ----------
    RN : ReactionNetwork
        The reaction network to analyze
    max_generator_size : int
        Maximum size of generators to build
    verbose : bool
        Whether to print detailed information
        
    Returns:
    -------
    tuple: (organizations, elementary_sos, statistics, flux_data)
        organizations: list of ElementarySO that are self-maintaining
        elementary_sos: all elementary semi-organizations found
        statistics: computation statistics
        flux_data: tuple of (flux_vectors, production_vectors)
    """
    # Import here to avoid circular imports
    from pyCOT.Persistent_Modules import compute_persistent_modules
    
    if verbose:
        print("="*80)
        print("COMPUTING ALL ORGANIZATIONS VIA ELEMENTARY SEMI-ORGANIZATIONS")
        print("="*80)
    
    # Step 1: Compute elementary semi-organizations using irreducible generators
    elementary_sos, generators, pm_statistics, erc_sorn = compute_persistent_modules(
        RN, max_generator_size, use_optimized=True
    )
    
    if verbose:
        print(f"\nStep 1 Complete: Found {len(elementary_sos)} elementary semi-organizations")
    
    # Step 2: Check which elementary SOs are self-maintaining (organizations)
    organizations, flux_vectors, production_vectors = check_elementary_sos_are_organizations(
        elementary_sos, RN, verbose
    )
    
    # Step 3: Compile comprehensive statistics
    final_statistics = {
        'persistent_module_stats': pm_statistics,
        'elementary_sos_found': len(elementary_sos),
        'organizations_found': len(organizations),
        'conversion_rate': len(organizations)/len(elementary_sos) if elementary_sos else 0,
        'p_erc_organizations': sum(1 for org in organizations if org.is_p_erc),
        'multi_erc_organizations': sum(1 for org in organizations if not org.is_p_erc),
        'average_org_size': np.mean([len(org.closure_species) for org in organizations]) if organizations else 0,
        'size_distribution': {
            'min': min(len(org.closure_species) for org in organizations) if organizations else 0,
            'max': max(len(org.closure_species) for org in organizations) if organizations else 0,
            'median': np.median([len(org.closure_species) for org in organizations]) if organizations else 0
        }
    }
    
    flux_data = (flux_vectors, production_vectors)
    
    if verbose:
        print("\n" + "="*80)
        print("FINAL RESULTS")
        print("="*80)
        print(f"Total organizations found: {len(organizations)}")
        print(f"P-ERC organizations: {final_statistics['p_erc_organizations']}")
        print(f"Multi-ERC organizations: {final_statistics['multi_erc_organizations']}")
        print(f"Average organization size: {final_statistics['average_org_size']:.1f} species")
        print(f"Size range: {final_statistics['size_distribution']['min']} - {final_statistics['size_distribution']['max']} species")
        print("="*80)
    
    return organizations, elementary_sos, final_statistics, flux_data


def analyze_organization_patterns(organizations, RN, verbose=True):
    """
    Analyze patterns in the found organizations.
    
    Parameters:
    ----------
    organizations : list of ElementarySO
        Organizations to analyze
    RN : ReactionNetwork
        The reaction network
    verbose : bool
        Whether to print analysis
        
    Returns:
    -------
    dict: Analysis results
    """
    if not organizations:
        return {'total': 0}
    
    analysis = {
        'total': len(organizations),
        'by_size': defaultdict(int),
        'by_type': {'p_erc': 0, 'multi_erc': 0},
        'species_usage': defaultdict(int),
        'erc_usage': defaultdict(int),
        'generation_patterns': {
            'synergy_usage': [],
            'complementarity_usage': [],
            'pathway_diversity': []
        }
    }
    
    for org in organizations:
        # Size distribution
        size = len(org.closure_species)
        analysis['by_size'][size] += 1
        
        # Type classification
        if org.is_p_erc:
            analysis['by_type']['p_erc'] += 1
        else:
            analysis['by_type']['multi_erc'] += 1
        
        # Species usage
        for species in org.closure_species:
            analysis['species_usage'][species.name] += 1
        
        # ERC usage
        for erc in org.constituent_ercs:
            analysis['erc_usage'][erc.label] += 1
        
        # Generation patterns
        gen_stats = org.get_generation_statistics()
        if 'synergy_stats' in gen_stats:
            analysis['generation_patterns']['synergy_usage'].extend(
                gen_stats['synergy_stats'].get('counts', [])
            )
        if 'complementarity_stats' in gen_stats:
            for comp_type in ['type1', 'type2', 'type3']:
                if comp_type in gen_stats['complementarity_stats']:
                    analysis['generation_patterns']['complementarity_usage'].extend(
                        gen_stats['complementarity_stats'][comp_type].get('counts', [])
                    )
        analysis['generation_patterns']['pathway_diversity'].append(
            gen_stats.get('pathway_count', 0)
        )
    
    if verbose:
        print("\nORGANIZATION ANALYSIS")
        print("-" * 50)
        print(f"Total organizations: {analysis['total']}")
        print(f"P-ERC organizations: {analysis['by_type']['p_erc']}")
        print(f"Multi-ERC organizations: {analysis['by_type']['multi_erc']}")
        
        print("\nSize distribution:")
        for size in sorted(analysis['by_size'].keys()):
            print(f"  {size} species: {analysis['by_size'][size]} organizations")
        
        print(f"\nMost common species (top 5):")
        sorted_species = sorted(analysis['species_usage'].items(), 
                              key=lambda x: x[1], reverse=True)
        for species, count in sorted_species[:5]:
            print(f"  {species}: {count} organizations")
        
        print(f"\nMost common ERCs (top 5):")
        sorted_ercs = sorted(analysis['erc_usage'].items(), 
                           key=lambda x: x[1], reverse=True)
        for erc, count in sorted_ercs[:5]:
            print(f"  {erc}: {count} organizations")
    
    return analysis


# Add to Persistent_Modules.py __all__ list
__all__ = [
    'minimize_sv',
    'is_self_maintaining', 
    'check_elementary_sos_are_organizations',
    'compute_all_organizations',
    'analyze_organization_patterns'
]   

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extension to Persistent_Modules.py: Hierarchical Organization Computation

Implementation of hierarchical organization detection using elementary semi-organizations
as building blocks for constructing larger organizations through productive novelty.
"""

from collections import defaultdict
from pyCOT.ERC_Hierarchy import ERC, ERC_Hierarchy, closure, species_list_to_names
from pyCOT.ERC_Synergy_Complementarity import build_erc_sorn
from pyCOT.Persistent_Modules import (
    compute_elementary_organizations, is_semi_self_maintaining, 
    is_self_maintaining, ElementarySO
)

"""
Extension to Persistent_Modules.py: Hierarchical Organization Computation

Implementation of hierarchical organization detection using elementary semi-organizations
as building blocks for constructing larger organizations through productive novelty.
"""

from collections import defaultdict
from pyCOT.ERC_Hierarchy import ERC, ERC_Hierarchy, closure, species_list_to_names
from pyCOT.ERC_Synergy_Complementarity import build_erc_sorn
from pyCOT.Persistent_Modules import (
    compute_elementary_organizations, is_semi_self_maintaining, 
    is_self_maintaining, ElementarySO
)
class CompositeModule:
    """
    Represents a module composed of multiple elementary semi-organizations.
    
    This class tracks both the constituent elementary SOs and the combined
    closure that results from their union, enabling efficient extension
    and deduplication during hierarchical construction.
    """
    
    def __init__(self, elementary_sos, combined_closure, constituent_ercs, RN):
        """
        Initialize a CompositeModule.
        
        Parameters
        ----------
        elementary_sos : list of ElementarySO
            Elementary SOs that compose this module
        combined_closure : list of Species
            Combined closure of all elementary SOs
        constituent_ercs : list of ERC
            All unique ERCs from all elementary SOs
        RN : ReactionNetwork
            The reaction network
        """
        self.elementary_sos = elementary_sos
        self.combined_closure = combined_closure
        self.constituent_ercs = constituent_ercs
        self.closure_signature = tuple(sorted(sp.name for sp in combined_closure))
        
        # Caching for SSM check
        self._cached_ssm = None
        self._RN_reference = None
        
        # Track species names for efficient lookup
        self.species_names = set(sp.name for sp in combined_closure)
        
    def check_ssm(self, RN):
        """Check if this module is semi-self-maintaining."""
        if self._cached_ssm is None or self._RN_reference != id(RN):
            self._cached_ssm = is_semi_self_maintaining(RN, self.combined_closure)
            self._RN_reference = id(RN)
        return self._cached_ssm
    
    def size(self):
        """Get the number of elementary SOs in this module."""
        return len(self.elementary_sos)
    
    def contains_elementary_so(self, elementary_so):
        """Check if this module already contains the given elementary SO."""
        return elementary_so in self.elementary_sos
    
    def get_elementary_so_labels(self):
        """Get string representation of constituent elementary SOs."""
        return [f"ESO_{i}" for i in range(len(self.elementary_sos))]
    
    def __repr__(self):
        eso_labels = self.get_elementary_so_labels()
        return f"CompositeModule({eso_labels}, species={len(self.combined_closure)})"


def has_productive_novelty_with_elementary_so(current_module, candidate_elementary_so, erc_sorn):
    """
    Check if adding candidate_elementary_so to current_module would introduce productive novelty.
    
    Productive novelty exists if any ERC in the current module has a synergistic or
    complementary relationship with any ERC in the candidate elementary SO.
    
    Parameters
    ----------
    current_module : CompositeModule
        Current module to potentially extend
    candidate_elementary_so : ElementarySO
        Elementary SO candidate for extension
    erc_sorn : ERC_SORN
        Pre-computed second-order reaction network
        
    Returns
    -------
    bool
        True if productive novelty exists
    """
    # Check if any ERC in current module has productive relationship 
    # with any ERC in candidate elementary SO
    for current_erc in current_module.constituent_ercs:
        for candidate_erc in candidate_elementary_so.constituent_ercs:
            if erc_sorn.has_productive_relationship(current_erc.label, candidate_erc.label):
                return True
    return False


def extend_module_with_elementary_so(current_module, elementary_so, RN):
    """
    Create a new CompositeModule by extending current_module with elementary_so.
    
    This function combines the elementary SOs, merges their ERCs, and computes
    the new combined closure using the existing closure computation logic.
    
    Parameters
    ----------
    current_module : CompositeModule
        Module to extend
    elementary_so : ElementarySO
        Elementary SO to add
    RN : ReactionNetwork
        The reaction network
        
    Returns
    -------
    CompositeModule
        New extended module
    """
    # Combine elementary SOs
    new_elementary_sos = current_module.elementary_sos + [elementary_so]
    
    # Combine ERCs (ensure uniqueness)
    new_constituent_ercs = list(current_module.constituent_ercs)
    for erc in elementary_so.constituent_ercs:
        if erc not in new_constituent_ercs:
            new_constituent_ercs.append(erc)
    
    # Compute new combined closure
    # Start with current closure species, add elementary SO species, then compute closure
    current_species = current_module.combined_closure
    additional_species = elementary_so.closure_species
    
    # Combine and remove duplicates while preserving Species objects
    all_species_dict = {sp.name: sp for sp in current_species}
    for sp in additional_species:
        all_species_dict[sp.name] = sp
    
    union_species = list(all_species_dict.values())
    
    # Use existing closure function to get final closure
    new_combined_closure = closure(RN, union_species)
    
    return CompositeModule(new_elementary_sos, new_combined_closure, new_constituent_ercs, RN)


def compute_all_organizations(RN, max_module_size=10, max_generator_size=8, 
                            elementary_sos=None, elementary_organizations=None, 
                            verbose=True):
    """
    Compute all organizations by extending elementary semi-organizations through productive novelty.
    
    This function implements the hierarchical construction approach where elementary SOs 
    are used as building blocks to create larger semi-organizations, which are then 
    filtered for self-maintenance to identify organizations.
    
    The algorithm works as follows:
    1. Compute elementary SOs and organizations (if not provided)
    2. Build ERC_SORN for efficient productive relationship lookup
    3. Initialize each elementary SO as a size-1 CompositeModule
    4. Iteratively extend modules by adding elementary SOs with productive novelty
    5. Use closure-based deduplication to avoid redundant computation
    6. Filter all generated modules for self-maintenance (organizations)
    
    Parameters
    ----------
    RN : ReactionNetwork
        The reaction network to analyze
    max_module_size : int
        Maximum number of elementary SOs per composite module
    max_generator_size : int
        Maximum size of generators when computing elementary SOs (if not provided)
    elementary_sos : list of ElementarySO, optional
        Pre-computed elementary semi-organizations
    elementary_organizations : list of ElementarySO, optional
        Pre-computed elementary organizations (subset of elementary_sos)
    verbose : bool
        Whether to print detailed progress information
        
    Returns
    -------
    tuple
        (all_organizations, all_semi_organizations, statistics, computation_data)
        all_organizations: list of CompositeModule that are self-maintaining
        all_semi_organizations: list of all CompositeModule generated
        statistics: comprehensive computation statistics
        computation_data: tuple of (elementary_sos, elementary_organizations, erc_sorn, flux_data)
    """
    if verbose:
        print("=" * 80)
        print("COMPUTING ALL ORGANIZATIONS VIA HIERARCHICAL EXTENSION")
        print("=" * 80)
    
    # Step 1: Ensure we have elementary SOs and organizations
    if elementary_sos is None or elementary_organizations is None:
        if verbose:
            print("Computing elementary semi-organizations and organizations...")
        
        elementary_organizations, elementary_sos, elem_stats, elem_flux_data = compute_elementary_organizations(
            RN, max_generator_size, verbose=verbose
        )
        
        if verbose:
            print(f"Found {len(elementary_sos)} elementary SOs, {len(elementary_organizations)} elementary organizations")
    else:
        if verbose:
            print(f"Using provided elementary SOs: {len(elementary_sos)} SOs, {len(elementary_organizations)} organizations")
        elem_stats = None
        elem_flux_data = None
    
    # Step 2: Build ERC_SORN for efficient productive relationship lookup
    if verbose:
        print("\nBuilding ERC_SORN for efficient productive relationship lookup...")
    
    # We need the hierarchy for SORN construction
    hierarchy = ERC_Hierarchy(RN)
    hierarchy.build_hierarchy_graph()
    erc_sorn = build_erc_sorn(hierarchy, RN)
    
    if verbose:
        sorn_stats = erc_sorn.get_statistics()
        print(f"SORN built: {sorn_stats['synergistic_pairs']} synergistic pairs, "
              f"{sorn_stats['complementary_pairs']} complementary pairs")
    
    # Step 3: Initialize with elementary SOs as size-1 modules
    if verbose:
        print(f"\nInitializing with {len(elementary_sos)} elementary SOs as size-1 modules...")
    
    # Map from closure signature to CompositeModule for deduplication
    modules_by_signature = {}
    
    # Convert each elementary SO to a CompositeModule
    for elem_so in elementary_sos:
        module = CompositeModule([elem_so], elem_so.closure_species, elem_so.constituent_ercs, RN)
        modules_by_signature[module.closure_signature] = module
    
    if verbose:
        print(f"Initialized {len(modules_by_signature)} unique size-1 modules")
    
    # Statistics tracking
    stats = {
        'elementary_sos_input': len(elementary_sos),
        'elementary_organizations_input': len(elementary_organizations),
        'modules_by_size': {1: len(modules_by_signature)},
        'extensions_explored_by_size': {},
        'duplicates_pruned_by_size': {},
        'total_modules_generated': len(modules_by_signature),
        'total_extensions_explored': 0,
        'total_duplicates_pruned': 0
    }
    
    # Step 4: Iterative extension phase
    if verbose:
        print(f"\nIteratively extending modules up to size {max_module_size}...")
        print("-" * 60)
    
    current_size = 1
    
    while current_size < max_module_size:
        # Get all modules of current size
        current_size_modules = [
            module for module in modules_by_signature.values()
            if module.size() == current_size
        ]
        
        if not current_size_modules:
            if verbose:
                print(f"No modules of size {current_size} to extend. Stopping.")
            break
        
        if verbose:
            print(f"\nSize {current_size} → {current_size + 1}:")
            print(f"  Processing {len(current_size_modules)} modules of size {current_size}")
        
        new_modules = {}
        extensions_explored = 0
        duplicates_pruned = 0
        
        # Try to extend each current module with each elementary SO
        for current_module in current_size_modules:
            for candidate_elem_so in elementary_sos:
                # Skip if this elementary SO is already in the module
                if current_module.contains_elementary_so(candidate_elem_so):
                    continue
                
                extensions_explored += 1
                
                # Check for productive novelty
                if has_productive_novelty_with_elementary_so(current_module, candidate_elem_so, erc_sorn):
                    # Create extended module
                    try:
                        extended_module = extend_module_with_elementary_so(current_module, candidate_elem_so, RN)
                        
                        # Check for duplicates using closure signature
                        if (extended_module.closure_signature not in modules_by_signature and
                            extended_module.closure_signature not in new_modules):
                            new_modules[extended_module.closure_signature] = extended_module
                        else:
                            duplicates_pruned += 1
                            
                    except Exception as e:
                        if verbose:
                            print(f"    Warning: Error extending module: {e}")
                        continue
        
        # Add new modules to main collection
        modules_by_signature.update(new_modules)
        
        # Update statistics
        stats['modules_by_size'][current_size + 1] = len(new_modules)
        stats['extensions_explored_by_size'][current_size] = extensions_explored
        stats['duplicates_pruned_by_size'][current_size] = duplicates_pruned
        stats['total_modules_generated'] += len(new_modules)
        stats['total_extensions_explored'] += extensions_explored
        stats['total_duplicates_pruned'] += duplicates_pruned
        
        if verbose:
            print(f"  Extensions explored: {extensions_explored}")
            print(f"  New unique modules: {len(new_modules)}")
            print(f"  Duplicates pruned: {duplicates_pruned}")
            if extensions_explored > 0:
                efficiency = 100 * duplicates_pruned / extensions_explored
                print(f"  Pruning efficiency: {efficiency:.1f}%")
        
        current_size += 1
    
    # Step 5: Collect all generated modules
    all_semi_organizations = list(modules_by_signature.values())
    
    if verbose:
        print(f"\n" + "=" * 60)
        print("EXTENSION PHASE COMPLETE")
        print("=" * 60)
        print(f"Total modules generated: {len(all_semi_organizations)}")
        print("Distribution by size:")
        for size in sorted(stats['modules_by_size'].keys()):
            count = stats['modules_by_size'][size]
            print(f"  Size {size}: {count} modules")
    
    # Step 6: Filter for organizations (self-maintaining modules)
    if verbose:
        print(f"\nFiltering {len(all_semi_organizations)} modules for self-maintenance (organizations)...")
    
    all_organizations = []
    organization_flux_vectors = []
    organization_production_vectors = []
    
    for i, module in enumerate(all_semi_organizations):
        if verbose and (i + 1) % 50 == 0:
            print(f"  Checking module {i+1}/{len(all_semi_organizations)}...")
        
        # Check self-maintenance using linear programming
        result = is_self_maintaining(RN, module.combined_closure)
        
        if result[0]:  # Is self-maintaining
            all_organizations.append(module)
            organization_flux_vectors.append(result[1])
            organization_production_vectors.append(result[2])
    
    # Step 7: Compile final statistics
    final_statistics = {
        'elementary_computation': elem_stats,
        'extension_phase': stats,
        'final_results': {
            'total_semi_organizations': len(all_semi_organizations),
            'total_organizations': len(all_organizations),
            'conversion_rate': len(all_organizations) / len(all_semi_organizations) if all_semi_organizations else 0,
            'size_distribution_sos': {
                size: sum(1 for m in all_semi_organizations if m.size() == size)
                for size in range(1, max_module_size + 1)
            },
            'size_distribution_orgs': {
                size: sum(1 for m in all_organizations if m.size() == size)
                for size in range(1, max_module_size + 1)
            }
        }
    }
    
    # Calculate organization rates by size
    org_rates_by_size = {}
    for size in range(1, max_module_size + 1):
        sos_count = final_statistics['final_results']['size_distribution_sos'].get(size, 0)
        orgs_count = final_statistics['final_results']['size_distribution_orgs'].get(size, 0)
        org_rates_by_size[size] = orgs_count / sos_count if sos_count > 0 else 0
    
    final_statistics['final_results']['organization_rates_by_size'] = org_rates_by_size
    
    flux_data = (organization_flux_vectors, organization_production_vectors)
    computation_data = (elementary_sos, elementary_organizations, erc_sorn, flux_data)
    
    if verbose:
        print("\n" + "=" * 80)
        print("COMPUTATION COMPLETE")
        print("=" * 80)
        print(f"Total semi-organizations found: {len(all_semi_organizations)}")
        print(f"Total organizations found: {len(all_organizations)}")
        print(f"Overall conversion rate: {final_statistics['final_results']['conversion_rate']:.1%}")
        
        print(f"\nOrganization rates by module size:")
        for size in sorted(org_rates_by_size.keys()):
            if stats['modules_by_size'].get(size, 0) > 0:
                rate = org_rates_by_size[size]
                count = final_statistics['final_results']['size_distribution_orgs'].get(size, 0)
                total = final_statistics['final_results']['size_distribution_sos'].get(size, 0)
                print(f"  Size {size}: {count}/{total} ({rate:.1%})")
        
        print("=" * 80)
    
    return elementary_sos, elementary_organizations, all_organizations, all_semi_organizations, final_statistics, computation_data


def analyze_hierarchical_organization_patterns(all_organizations, all_semi_organizations, 
                                             elementary_sos, verbose=True):
    """
    Analyze patterns in the hierarchically constructed organizations.
    
    Parameters
    ----------
    all_organizations : list of CompositeModule
        All organizations found
    all_semi_organizations : list of CompositeModule
        All semi-organizations generated
    elementary_sos : list of ElementarySO
        Elementary semi-organizations used as building blocks
    verbose : bool
        Whether to print analysis results
        
    Returns
    -------
    dict
        Analysis results
    """
    analysis = {
        'total_organizations': len(all_organizations),
        'total_semi_organizations': len(all_semi_organizations),
        'elementary_sos_used': len(elementary_sos),
        'organization_size_stats': {},
        'semi_organization_size_stats': {},
        'elementary_so_usage': defaultdict(int),
        'most_productive_elementary_sos': [],
        'size_efficiency': {}
    }
    
    # Analyze by size
    org_sizes = [org.size() for org in all_organizations]
    so_sizes = [so.size() for so in all_semi_organizations]
    
    if org_sizes:
        analysis['organization_size_stats'] = {
            'min': min(org_sizes),
            'max': max(org_sizes),
            'mean': sum(org_sizes) / len(org_sizes),
            'distribution': defaultdict(int)
        }
        for size in org_sizes:
            analysis['organization_size_stats']['distribution'][size] += 1
    
    if so_sizes:
        analysis['semi_organization_size_stats'] = {
            'min': min(so_sizes),
            'max': max(so_sizes),
            'mean': sum(so_sizes) / len(so_sizes),
            'distribution': defaultdict(int)
        }
        for size in so_sizes:
            analysis['semi_organization_size_stats']['distribution'][size] += 1
    
    # Track elementary SO usage
    for org in all_organizations:
        for elem_so in org.elementary_sos:
            analysis['elementary_so_usage'][id(elem_so)] += 1
    
    # Find most productive elementary SOs
    usage_counts = list(analysis['elementary_so_usage'].values())
    if usage_counts:
        most_used_count = max(usage_counts)
        analysis['most_productive_elementary_sos'] = [
            (elem_so_id, count) for elem_so_id, count in analysis['elementary_so_usage'].items()
            if count == most_used_count
        ]
    
    # Calculate organization efficiency by size
    for size in range(1, max(so_sizes) + 1 if so_sizes else 1):
        org_count = analysis['organization_size_stats']['distribution'].get(size, 0) if org_sizes else 0
        so_count = analysis['semi_organization_size_stats']['distribution'].get(size, 0) if so_sizes else 0
        analysis['size_efficiency'][size] = org_count / so_count if so_count > 0 else 0
    
    if verbose:
        print("\nHIERARCHICAL ORGANIZATION ANALYSIS")
        print("-" * 50)
        print(f"Total organizations: {analysis['total_organizations']}")
        print(f"Total semi-organizations: {analysis['total_semi_organizations']}")
        print(f"Elementary SOs used as building blocks: {analysis['elementary_sos_used']}")
        
        if org_sizes:
            print(f"\nOrganization size statistics:")
            print(f"  Range: {analysis['organization_size_stats']['min']} - {analysis['organization_size_stats']['max']} elementary SOs")
            print(f"  Average: {analysis['organization_size_stats']['mean']:.1f} elementary SOs")
            print(f"  Distribution:")
            for size in sorted(analysis['organization_size_stats']['distribution'].keys()):
                count = analysis['organization_size_stats']['distribution'][size]
                print(f"    Size {size}: {count} organizations")
        
        print(f"\nOrganization efficiency by size:")
        for size in sorted(analysis['size_efficiency'].keys()):
            if analysis['semi_organization_size_stats']['distribution'].get(size, 0) > 0:
                efficiency = analysis['size_efficiency'][size]
                org_count = analysis['organization_size_stats']['distribution'].get(size, 0)
                so_count = analysis['semi_organization_size_stats']['distribution'].get(size, 0)
                print(f"  Size {size}: {org_count}/{so_count} ({efficiency:.1%})")
    
    return analysis


# Add to __all__ for proper module export
__all__ = [
    'CompositeModule',
    'compute_all_organizations',
    'analyze_hierarchical_organization_patterns',
    'has_productive_novelty_with_elementary_so',
    'extend_module_with_elementary_so'
]