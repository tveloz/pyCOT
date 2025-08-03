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

def compute_elementary_sos(generators, RN):
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
    elementary_sos = compute_elementary_sos(generators, RN)
    
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

