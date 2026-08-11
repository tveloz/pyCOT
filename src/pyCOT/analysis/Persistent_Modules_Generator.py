#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Persistent_Modules_Generator.py

*** For new work, use pyCOT.analysis.organizations.compute_organizations. ***

That module is now the recommended way this library computes organizations:
genome-scale ERC -> fundamental relations -> fundamental hierarchy -> EPM/ESPM
exploration, followed by the same LP self-maintenance verification this file
originated (minimize_sv, re-exported below unchanged). It replaces the
ElementarySO/Organization/OrganizationHierarchy construction and the
brute_force_organizations/compute_all_organizations enumerators kept in this
file below, which do not scale past small networks (brute_force_organizations
took ~52 minutes on a 92-species, 23-ERC network; compute_organizations
solves the same network — and reproduces the same published organizations,
species-for-species — in under 2 seconds; validated this way against all
five scenarios of Centler et al. 2006's E. coli sugar-metabolism model,
see projects/COT_fundamental_Generators/scripts/reproduce_centler2006.py).

Everything below this point is kept for backward compatibility (existing
scripts that import ElementarySO / Organization / OrganizationHierarchy /
compute_all_organizations / brute_force_organizations directly) and as a
reference implementation, not because it is still the preferred path.

    minimize_sv / check_self_maintenance / diagnose_self_maintenance
        Re-exported unchanged from pyCOT.analysis.organizations.self_maintenance
        (the validated LP kernel originally developed in this file — moved
        there so it has one canonical home, importable from either path).

LEGACY (this file, below) — Elementary and Non-Elementary Organization
Construction via irreducible generators and productive novelty:

    ElementarySO: Elementary semi-organization from irreducible generators
    Organization: Non-elementary organization containing multiple ElementarySOs
    OrganizationHierarchy: Manages containment relationships between organizations
    compute_elementary_sos: Convert SSM generators to ElementarySOs
    find_productive_so_extensions: Extend SOs through productive novelty
    build_non_elementary_organizations: Hierarchical organization construction
    compute_all_organizations: Old orchestrator (pairwise-limited, not exhaustive)
    brute_force_organizations: Old exhaustive orchestrator (exponential, exact)

Author: Based on theoretical work by Tomas Veloz et al.
"""

import time
import numpy as np
from scipy.optimize import linprog
from itertools import combinations
from collections import defaultdict, Counter
from pyCOT.analysis.ERC_Hierarchy import ERC, ERC_Hierarchy, closure, species_list_to_names

# Import from our SORN_Generators module
from pyCOT.analysis.SORN_Generators import (
    IrreducibleGenerator, ERC_SORN, ProductiveExtension,
    identify_p_ercs, build_erc_sorn, find_productive_extensions,
    build_irreducible_generators, is_semi_self_maintaining
)

# The validated LP kernel now lives in pyCOT.analysis.organizations.self_maintenance;
# re-exported here so `from pyCOT.analysis.Persistent_Modules_Generator import
# minimize_sv` (used directly by projects/Decomposition_Theorem/decomp/circuits.py
# and .witness.py) keeps working without any change on their part.
from pyCOT.analysis.organizations.self_maintenance import (  # noqa: F401
    minimize_sv, check_self_maintenance, diagnose_self_maintenance,
)

# ============================================================================
# LEGACY — CORE CLASSES (see module docstring: prefer pyCOT.analysis.organizations)
# ============================================================================

class ElementarySO:
    """
    Elementary Semi-Organization built from irreducible generators.
    
    Represents a semi-organization that cannot be decomposed into smaller
    semi-organizations, constructed from SSM irreducible generators.
    """
    
    def __init__(self, closure_species, constituent_ercs, RN):
        """
        Initialize an elementary semi-organization.
        
        Parameters
        ----------
        closure_species : list of Species
            Species in the closure of this SO
        constituent_ercs : list of ERC
            ERCs that contribute to this SO
        RN : ReactionNetwork
            The reaction network
        """
        self.closure_species = closure_species
        self.constituent_ercs = constituent_ercs
        self.RN = RN
        self.generating_pathways = []  # List of IrreducibleGenerator objects
        
        # Check if this is a P-ERC (primary ERC)
        self.is_p_erc = (len(self.constituent_ercs) == 1 and 
                        len(self.closure_species) == len(self.constituent_ercs[0].get_closure(RN)))
        
        # Generate unique ID from sorted species names (truncated for readability)
        species_names = sorted([sp.name for sp in self.closure_species])
        if len(species_names) <= 3:
            name_part = '-'.join(species_names)
        else:
            name_part = '-'.join(species_names[:3]) + f'-etc{len(species_names)-3}'
        self.id = f"ESO_{name_part}"
        
        # Cached statistics and self-maintenance result
        self._synergy_stats = None
        self._complementarity_stats = None
        self._size_stats = None
        self._self_maintenance_result = None
    
    def add_generator(self, generator):
        """Add an irreducible generator that produces this SO."""
        if not isinstance(generator, IrreducibleGenerator):
            raise TypeError(f"Expected IrreducibleGenerator, got {type(generator)}")
        self.generating_pathways.append(generator)
        # Clear cached statistics
        self._synergy_stats = None
        self._complementarity_stats = None
        self._size_stats = None
    
    def size(self):
        """Return the number of constituent ERCs."""
        return len(self.constituent_ercs)
    
    def check_self_maintenance(self):
        """
        Check if this elementary SO is self-maintaining (i.e., an organization).
        
        Returns
        -------
        tuple
            (is_organization: bool, flux_vector: np.array or None)
        """
        if self._self_maintenance_result is None:
            result = check_self_maintenance(self.closure_species, self.RN)
            is_org = result[0]
            flux_vector = result[1] if is_org else None
            self._self_maintenance_result = (is_org, flux_vector)
        
        return self._self_maintenance_result
    
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
                'synergy_stats': {'mean': 0, 'counts': []},
                'complementarity_stats': {
                    'type1': {'mean': 0, 'counts': []},
                    'type2': {'mean': 0, 'counts': []},
                    'type3': {'mean': 0, 'counts': []}
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
        
        self._synergy_stats = {
            'mean': sum(synergy_counts) / len(synergy_counts) if synergy_counts else 0,
            'max': max(synergy_counts) if synergy_counts else 0,
            'min': min(synergy_counts) if synergy_counts else 0,
            'counts': synergy_counts
        }
        
        self._complementarity_stats = {
            'type1': {
                'mean': sum(comp_type1_counts) / len(comp_type1_counts) if comp_type1_counts else 0,
                'counts': comp_type1_counts
            },
            'type2': {
                'mean': sum(comp_type2_counts) / len(comp_type2_counts) if comp_type2_counts else 0,
                'counts': comp_type2_counts
            },
            'type3': {
                'mean': sum(comp_type3_counts) / len(comp_type3_counts) if comp_type3_counts else 0,
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
        return f"ElementarySO(id={self.id}, species={len(self.closure_species)}, ERCs={len(self.constituent_ercs)}, P-ERC={self.is_p_erc})"

class Organization:
    """
   Organization containing one or multiple ElementarySOs.
    
    Built through productive combination of elementary semi-organizations,
    representing higher-order organizational structures.
    """
    
    def __init__(self, elementary_sos, construction_method="productive_combination"):
        """
        Initialize a non-elementary organization.
        
        Parameters
        ----------
        elementary_sos : list of ElementarySO
            Elementary SOs that compose this organization
        construction_method : str
            Method used to construct this organization
        """
        self.elementary_sos = elementary_sos
        self.construction_method = construction_method
        
        # Compute combined closure
        all_species = []
        for eso in elementary_sos:
            all_species.extend(eso.closure_species)
        
        # Remove duplicates while preserving Species objects
        unique_species_dict = {sp.name: sp for sp in all_species}
        self.combined_closure = list(unique_species_dict.values())
        
        # Collect all constituent ERCs
        all_ercs = []
        for eso in elementary_sos:
            all_ercs.extend(eso.constituent_ercs)
        unique_ercs_dict = {erc.label: erc for erc in all_ercs}
        self.constituent_ercs = list(unique_ercs_dict.values())
        
        # Generate ID
        eso_ids = sorted([eso.id for eso in elementary_sos])
        if len(eso_ids) <= 2:
            id_part = '+'.join(eso_ids)
        else:
            id_part = '+'.join(eso_ids[:2]) + f'+{len(eso_ids)-2}more'
        self.id = f"ORG_{id_part}"
        
        # Cached self-maintenance result
        self._self_maintenance_result = None
    
    def size(self):
        """Return the number of elementary SOs."""
        return len(self.elementary_sos)
    
    def check_self_maintenance(self, RN):
        """
        Check if this organization is self-maintaining.
        
        Parameters
        ----------
        RN : ReactionNetwork
            The reaction network
            
        Returns
        -------
        tuple
            (is_organization: bool, flux_vector: np.array or None)
        """
        if self._self_maintenance_result is None:
            result = check_self_maintenance(self.combined_closure, RN)
            is_org = result[0]
            flux_vector = result[1] if is_org else None
            self._self_maintenance_result = (is_org, flux_vector)
        
        return self._self_maintenance_result
    
    def __repr__(self):
        return f"Organization(id={self.id}, ESOs={len(self.elementary_sos)}, species={len(self.combined_closure)})"

class OrganizationHierarchy:
    """
    Manages containment relationships between organizations.
    
    Tracks which organizations contain others and provides methods for
    hierarchical analysis and maximal organization identification.
    """
    
    def __init__(self):
        """Initialize an empty organization hierarchy."""
        self.organizations = []  # All organizations (elementary and non-elementary)
        self.elementary_organizations = []  # Only elementary organizations
        self.containment_graph = None  # NetworkX graph of containment relationships
        self._containment_matrix = None
    
    def add_organization(self, organization):
        """Add an organization to the hierarchy."""
        self.organizations.append(organization)
        if isinstance(organization, ElementarySO):
            self.elementary_organizations.append(organization)
        # Clear cached structures
        self.containment_graph = None
        self._containment_matrix = None
    
    def build_containment_relationships(self):
        """Build the containment graph between organizations."""
        import networkx as nx
        
        self.containment_graph = nx.DiGraph()
        
        # Add all organizations as nodes
        for org in self.organizations:
            self.containment_graph.add_node(org.id, organization=org)
        
        # Add containment edges
        for i, org1 in enumerate(self.organizations):
            for j, org2 in enumerate(self.organizations):
                if i != j:
                    # Check if org1 contains org2 (org2's closure is subset of org1's closure)
                    # Get closure for org1
                    if hasattr(org1, 'combined_closure'):
                        closure1_names = set(sp.name for sp in org1.combined_closure)
                    else:  # ElementarySO
                        closure1_names = set(sp.name for sp in org1.closure_species)
                    
                    # Get closure for org2  
                    if hasattr(org2, 'combined_closure'):
                        closure2_names = set(sp.name for sp in org2.combined_closure)
                    else:  # ElementarySO
                        closure2_names = set(sp.name for sp in org2.closure_species)
                    
                    if closure2_names.issubset(closure1_names) and closure1_names != closure2_names:
                        self.containment_graph.add_edge(org1.id, org2.id)
    
    def find_maximal_organizations(self):
        """Find organizations that are not contained by any other organization."""
        if self.containment_graph is None:
            self.build_containment_relationships()
        
        maximal_orgs = []
        for org in self.organizations:
            # An organization is maximal if it has no incoming edges (not contained by others)
            if self.containment_graph.in_degree(org.id) == 0:
                maximal_orgs.append(org)
        
        return maximal_orgs
    
    def get_containment_tree(self):
        """Get the containment relationships as a tree structure."""
        if self.containment_graph is None:
            self.build_containment_relationships()
        
        return self.containment_graph
    
    def get_statistics(self):
        """Get statistics about the organization hierarchy."""
        elementary_count = len(self.elementary_organizations)
        non_elementary_count = len(self.organizations) - elementary_count
        
        # Count organizations vs semi-organizations
        orgs_count = 0
        sos_count = 0
        
        for org in self.organizations:
            if hasattr(org, 'check_self_maintenance'):
                # For ElementarySO, we already computed this
                if isinstance(org, ElementarySO):
                    is_org, _ = org.check_self_maintenance()
                else:
                    # For non-elementary, need to pass RN (simplified for now)
                    is_org = True  # Assume for now, should be checked properly
                
                if is_org:
                    orgs_count += 1
                else:
                    sos_count += 1
        
        return {
            'total_organizations': len(self.organizations),
            'elementary_organizations': elementary_count,
            'non_elementary_organizations': non_elementary_count,
            'self_maintaining_count': orgs_count,
            'semi_organizations_count': sos_count,
            'maximal_organizations': len(self.find_maximal_organizations())
        }

# (minimize_sv / check_self_maintenance / diagnose_self_maintenance now live in
#  pyCOT.analysis.organizations.self_maintenance -- re-exported at top of this file)

# ============================================================================
# MAIN FUNCTIONS
# ============================================================================

def compute_elementary_sos(generators, RN, hierarchy=None):
    """
    Compute elementary semi-organizations from irreducible generators.
    
    Parameters
    ----------
    generators : list of IrreducibleGenerator
        All irreducible generators
    RN : ReactionNetwork
        The reaction network
    hierarchy : ERC_Hierarchy, optional
        ERC hierarchy for finding constituent ERCs
        
    Returns
    -------
    list of ElementarySO
        Elementary semi-organizations found
    """
    print("Computing elementary semi-organizations...")
    
    if not generators:
        print("No generators provided.")
        return []
    
    # Filter for SSM generators: These are all semi-organizations
    ssm_generators = []
    for gen in generators:
        if gen.check_ssm(RN):
            ssm_generators.append(gen)
    
    print(f"Found {len(ssm_generators)} SSM generators from {len(generators)} total generators.")
    
    # Group generators by their closure (same closure = same SO)
    closure_groups = defaultdict(list)
    for gen in ssm_generators:
        closure_species = gen.get_closure(RN)
        closure_signature = tuple(sorted(sp.name for sp in closure_species))
        closure_groups[closure_signature].append(gen)
    
    print(f"Found {len(closure_groups)} unique closures (potential semi-organizations).")
    
    # Create ElementarySO objects
    elementary_sos = []
    
    for closure_signature, gens in closure_groups.items():
        representative_gen = gens[0]
        closure_species = representative_gen.get_closure(RN)
        
        # Find constituent ERCs for this closure
        constituent_ercs = []
        if hierarchy:
            closure_names_set = set(sp.name for sp in closure_species)
            for erc in hierarchy.ercs:
                erc_closure_names = set(sp.name for sp in erc.get_closure(RN))
                if erc_closure_names.issubset(closure_names_set):
                    constituent_ercs.append(erc)
        else:
            # If no hierarchy provided, use ERCs from generators
            unique_ercs = {}
            for gen in gens:
                for erc in gen.erc_sequence:
                    unique_ercs[erc.label] = erc
            constituent_ercs = list(unique_ercs.values())
        
        elementary_so = ElementarySO(closure_species, constituent_ercs, RN)
        for gen in gens:
            elementary_so.add_generator(gen)
        elementary_sos.append(elementary_so)
    
    print(f"Created {len(elementary_sos)} elementary semi-organizations.")
    return elementary_sos

def find_productive_so_extensions(elementary_sos, erc_sorn):
    """
    Find productive extensions for elementary semi-organizations.
    
    This adapts the generator extension logic to work with complete SOs,
    looking for ways to productively combine elementary SOs.
    
    Parameters
    ----------
    elementary_sos : list of ElementarySO
        Elementary SOs to find extensions for
    erc_sorn : ERC_SORN
        Pre-computed relationship network
        
    Returns
    -------
    dict
        Mapping from SO pairs to their productive relationship details
    """
    productive_combinations = {}
    
    print(f"Finding productive extensions for {len(elementary_sos)} elementary SOs...")
    
    # Check all pairs of elementary SOs for productive relationships
    combinations_checked = 0
    productive_found = 0
    
    for so1, so2 in combinations(elementary_sos, 2):
        combinations_checked += 1
        
        # Check if any ERCs from so1 have productive relationships with ERCs from so2
        productive_relationships = []
        
        for erc1 in so1.constituent_ercs:
            for erc2 in so2.constituent_ercs:
                # Check for synergies
                synergies = erc_sorn.get_synergies(erc1.label, erc2.label)
                if synergies:
                    for synergy in synergies:
                        productive_relationships.append({
                            'type': 'synergy',
                            'erc1': erc1.label,
                            'erc2': erc2.label,
                            'details': synergy
                        })
                
                # Check for complementarities
                complementarities = erc_sorn.get_complementarities(erc1.label, erc2.label)
                if complementarities:
                    for comp in complementarities:
                        productive_relationships.append({
                            'type': 'complementarity',
                            'erc1': erc1.label,
                            'erc2': erc2.label,
                            'details': comp
                        })
        
        if productive_relationships:
            productive_combinations[(so1.id, so2.id)] = productive_relationships
            productive_found += 1
    
    print(f"Found {productive_found} productive SO combinations from {combinations_checked} pairs checked.")
    return productive_combinations

def build_non_elementary_organizations(elementary_sos, erc_sorn, max_size=5, verbose=True):
    """
    Build non-elementary organizations through productive combination of elementary SOs.
    
    Parameters
    ----------
    elementary_sos : list of ElementarySO
        Elementary semi-organizations to combine
    erc_sorn : ERC_SORN
        Pre-computed relationship network
    max_size : int
        Maximum number of elementary SOs to combine
    verbose : bool
        Whether to print progress information
        
    Returns
    -------
    list of Organization
        Non-elementary organizations found
    """
    if verbose:
        print(f"Building non-elementary organizations from {len(elementary_sos)} elementary SOs...")
    
    # Find productive combinations
    productive_combinations = find_productive_so_extensions(elementary_sos, erc_sorn)
    
    if not productive_combinations:
        print("No productive combinations found.")
        return []
    
    # Build organizations of increasing size, tracking unique closures
    organizations = []
    seen_closures = set()  # Track unique closures to avoid duplicates

    # Pre-populate with elementary organization closures to avoid duplicating them
    for eso in elementary_sos:
        closure_signature = frozenset(sp.name for sp in eso.closure_species)
        seen_closures.add(closure_signature)

    current_combinations = list(productive_combinations.keys())

    # Start with size 2 (pairs)
    for (so1_id, so2_id) in current_combinations:
        so1 = next(so for so in elementary_sos if so.id == so1_id)
        so2 = next(so for so in elementary_sos if so.id == so2_id)

        org = Organization([so1, so2], "productive_combination")

        # Create a hashable closure signature (frozenset of species names)
        closure_signature = frozenset(sp.name for sp in org.combined_closure)

        # Only add if this closure is new (not in elementary or already seen)
        if closure_signature not in seen_closures:
            seen_closures.add(closure_signature)
            organizations.append(org)

    if verbose:
        print(f"Built {len(organizations)} unique size-2 organizations (from {len(current_combinations)} pairs).")
    
    # For larger sizes, we could implement hierarchical extension
    # For now, we focus on pairwise combinations
    
    return organizations

def compute_all_organizations(RN, max_generator_size=8, max_organization_size=5, verbose=True):
    """
    Compute all organizations (elementary and non-elementary) for a reaction network.
    
    This is the main orchestrator function that coordinates the entire computation.
    
    Parameters
    ----------
    RN : ReactionNetwork
        The reaction network to analyze
    max_generator_size : int
        Maximum size for irreducible generators
    max_organization_size : int
        Maximum size for non-elementary organizations
    verbose : bool
        Whether to print detailed progress
        
    Returns
    -------
    dict
        Complete results including all organizations and statistics
    """
    if verbose:
        print("="*80)
        print("COMPUTING ALL ORGANIZATIONS")
        print("="*80)
        print(f"Network: {len(RN.species())} species, {len(RN.reactions())} reactions")
    
    start_time = time.time()
    
    # Step 1: Create ERC hierarchy
    if verbose:
        print("\nStep 1: Creating ERC hierarchy...")
    hierarchy_start = time.time()
    ercs = ERC.ERCs(RN)
    hierarchy = ERC_Hierarchy(RN,ercs)
    hierarchy_time = time.time() - hierarchy_start
    if verbose:
        print(f"✅ Created hierarchy: {len(hierarchy.ercs)} ERCs in {hierarchy_time:.2f}s")
    
    # Step 2: Build ERC_SORN
    if verbose:
        print("\nStep 2: Building ERC_SORN...")
    sorn_start = time.time()
    erc_sorn = build_erc_sorn(hierarchy, RN)
    sorn_time = time.time() - sorn_start
    if verbose:
        print(f"✅ Built SORN in {sorn_time:.2f}s")
    
    # Step 3: Build irreducible generators
    if verbose:
        print("\nStep 3: Building irreducible generators...")
    generators_start = time.time()
    generators, _ = build_irreducible_generators(hierarchy, RN, erc_sorn, max_generator_size, verbose)
    generators_time = time.time() - generators_start
    if verbose:
        print(f"✅ Built {len(generators)} generators in {generators_time:.2f}s")

    # Step 4: Compute elementary semi-organizations
    if verbose:
        print("\nStep 4: Computing elementary semi-organizations...")
    esos_start = time.time()
    elementary_sos = compute_elementary_sos(generators, RN, hierarchy)
    esos_time = time.time() - esos_start
    if verbose:
        print(f"✅ Found {len(elementary_sos)} elementary SOs in {esos_time:.2f}s")
    
    # Step 5: Check self-maintenance for elementary SOs
    if verbose:
        print("\nStep 5: Checking self-maintenance...")
    sm_start = time.time()
    elementary_organizations = []
    for eso in elementary_sos:
        is_org, flux = eso.check_self_maintenance()
        if is_org:
            elementary_organizations.append(eso)
    sm_time = time.time() - sm_start
    if verbose:
        print(f"✅ Found {len(elementary_organizations)} elementary organizations in {sm_time:.2f}s")
    
    # Step 6: Build non-elementary organizations (if we have elementary ones)
    non_elementary_organizations = []
    if elementary_organizations and len(elementary_organizations) >= 2:
        if verbose:
            print("\nStep 6: Building non-elementary organizations...")
        non_elem_start = time.time()
        # Use the SORN from generators construction
        erc_sorn = build_erc_sorn(hierarchy, RN)  # Reuse if possible
        non_elementary_organizations = build_non_elementary_organizations(
            elementary_organizations, erc_sorn, max_organization_size, verbose)
        non_elem_time = time.time() - non_elem_start
        if verbose:
            print(f"✅ Built {len(non_elementary_organizations)} non-elementary organizations in {non_elem_time:.2f}s")
    else:
        non_elem_time = 0
        if verbose:
            print("\nStep 5: Skipping non-elementary organizations (need ≥2 elementary organizations)")
    
    # Step 7: Create organization hierarchy
    if verbose:
        print("\nStep 7: Creating organization hierarchy...")
    hierarchy_obj = OrganizationHierarchy()
    
    # Add all organizations to hierarchy
    for eso in elementary_organizations:
        hierarchy_obj.add_organization(eso)
    for org in non_elementary_organizations:
        hierarchy_obj.add_organization(org)
    
    # Build containment relationships
    hierarchy_obj.build_containment_relationships()
    
    total_time = time.time() - start_time
    
    # Compile results
    results = {
        'elementary_sos': elementary_sos,
        'elementary_organizations': elementary_organizations,
        'non_elementary_organizations': non_elementary_organizations,
        'all_organizations': elementary_organizations + non_elementary_organizations,
        'organization_hierarchy': hierarchy_obj,
        'generators': generators,
        'erc_hierarchy': hierarchy,
        'statistics': {
            'total_generators': len(generators),
            'total_elementary_sos': len(elementary_sos),
            'total_elementary_organizations': len(elementary_organizations),
            'total_non_elementary_organizations': len(non_elementary_organizations),
            'total_organizations': len(elementary_organizations) + len(non_elementary_organizations),
            'conversion_rate_sos_to_orgs': (
                len(elementary_organizations) / len(elementary_sos) 
                if elementary_sos else 0
            ),
            'hierarchy_statistics': hierarchy_obj.get_statistics()
        },
        'timing': {
            'hierarchy_creation': hierarchy_time,
            'generator_construction': generators_time,
            'elementary_sos': esos_time,
            'self_maintenance_check': sm_time,
            'non_elementary_construction': non_elem_time,
            'total_time': total_time
        }
    }
    
    if verbose:
        print(f"\n" + "="*80)
        print("FINAL RESULTS")
        print("="*80)
        print(f"✅ Elementary semi-organizations: {len(elementary_sos)}")
        print(f"✅ Elementary organizations: {len(elementary_organizations)}")
        print(f"✅ Non-elementary organizations: {len(non_elementary_organizations)}")
        print(f"✅ Total organizations: {len(elementary_organizations) + len(non_elementary_organizations)}")
        if elementary_sos:
            conversion_rate = len(elementary_organizations) / len(elementary_sos) * 100
            print(f"✅ Conversion rate (SOs → Organizations): {conversion_rate:.1f}%")
        print(f"✅ Total computation time: {total_time:.2f}s")
        print("="*80)
    
    return results


# ============================================================================
# BRUTE FORCE ORGANIZATION CALCULATOR
# ============================================================================

def brute_force_organizations(RN, max_combo_size=None, early_stop=True, verbose=True):
    """
    Compute all closures, semi-organizations, and organizations by exhaustively
    enumerating all combinations of ERCs in lexicographic order.

    For each combination, the closure is computed; if new, it is stored and
    checked for semi-self-maintenance (structural) and self-maintenance (LP).

    Parameters
    ----------
    RN : ReactionNetwork
        The reaction network.
    max_combo_size : int or None
        Maximum number of ERCs per combination.  None = all sizes.
    early_stop : bool
        Stop increasing combo size once a full size-k pass produces no new
        closures (default True).  Keeps runtime practical for most networks.
    verbose : bool
        Print progress.

    Returns
    -------
    dict with keys:
        'closures'          : list of frozenset (species names), sorted by size
        'semi_organizations': list of frozenset  (SSM closures)
        'organizations'     : list of frozenset  (SM closures)
    """
    # ── Step 1: get ERCs, sort by label for deterministic ordering ──────────
    ercs = ERC.ERCs(RN)
    ercs.sort(key=lambda e: e.label)
    n = len(ercs)
    cap = min(max_combo_size, n) if max_combo_size is not None else n

    if verbose:
        print(f"\n{'='*60}")
        print(f"  BRUTE FORCE ORGANIZATION CALCULATOR")
        print(f"{'='*60}")
        print(f"  ERCs found : {n}  |  max combo size : {cap}")

    all_closures_set = set()   # frozensets already processed
    closures     = []          # all unique closures (frozensets)
    semi_orgs    = []          # SSM closures
    orgs         = []          # SM closures

    checked  = 0
    new_found = 0

    # ── Step 2: enumerate combinations in lexicographic order ───────────────
    for size in range(1, cap + 1):
        new_this_size = 0
        for combo_idx in combinations(range(n), size):
            checked += 1

            # union of species names from the selected ERC closures
            sp_names_union = set()
            for i in combo_idx:
                sp_names_union.update(ercs[i]._closure_names)

            # compute closure of the union
            sp_objects = [sp for sp in RN.species() if sp.name in sp_names_union]
            cl_objects = closure(RN, sp_objects)
            cl_frozen  = frozenset(sp.name for sp in cl_objects)

            if cl_frozen in all_closures_set:
                continue

            all_closures_set.add(cl_frozen)
            closures.append(cl_frozen)
            new_found  += 1
            new_this_size += 1

            cl_sp = [sp for sp in RN.species() if sp.name in cl_frozen]

            # ── SSM check (structural) ───────────────────────────────────
            if not is_semi_self_maintaining(RN, cl_sp):
                continue
            semi_orgs.append(cl_frozen)

            # ── SM check: structural pre-check + LP ──────────────────────
            sub_RN = RN.sub_reaction_network(cl_sp)
            S_obj  = sub_RN.stoichiometry_matrix()
            S      = np.asarray(S_obj, dtype=float)

            # species with no producer are always blocking — skip LP
            structurally_blocked = any(
                np.all(S[i, :] <= 0) and np.any(S[i, :] < 0)
                for i in range(S.shape[0])
            )
            if structurally_blocked:
                continue

            sm_result = check_self_maintenance(cl_sp, RN)
            if sm_result[0]:
                orgs.append(cl_frozen)

        if verbose:
            print(f"  size {size:2d}: {new_this_size} new closures  "
                  f"(total {new_found}, checked {checked})")

        # early stop: if no new closures at this size, larger combos won't add any
        if early_stop and new_this_size == 0:
            if verbose:
                print(f"  → no new closures at size {size}; stopping early.")
            break

    _key = lambda s: (len(s), sorted(s))
    closures  = sorted(closures,  key=_key)
    semi_orgs = sorted(semi_orgs, key=_key)
    orgs      = sorted(orgs,      key=_key)

    if verbose:
        print(f"\n  Combinations checked : {checked}")
        print(f"  Unique closures      : {new_found}")
        print(f"  Semi-organizations   : {len(semi_orgs)}")
        print(f"  Organizations        : {len(orgs)}")
        print(f"{'='*60}\n")

    return {
        'closures'          : closures,
        'semi_organizations': semi_orgs,
        'organizations'     : orgs,
    }