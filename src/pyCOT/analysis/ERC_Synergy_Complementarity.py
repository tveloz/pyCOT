#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Core Synergy and Complementarity Classes and Functions

Organization:
1. Core Classes (ERC_Synergy, ERC_Complementarity)
2. Helper Functions
3. Synergy Detection Functions (basic → maximal → fundamental)
4. Complementarity Detection Functions
5. ERC_SORN Class
6. Factory Functions
"""
import time
from itertools import combinations
import networkx as nx
from collections import defaultdict
from pyCOT.analysis.ERC_Hierarchy import ERC, ERC_Hierarchy, closure, species_list_to_names

# ============================================================================
# 1. CORE CLASSES
# ============================================================================

class ERC_Synergy:
    """
    Represents a synergistic relationship between ERCs.
    
    A synergy occurs when two ERCs together can generate a target ERC
    that neither could generate alone.
    """
    def __init__(self, reactants, product, synergy_type="regular"):
        self.reactants = reactants
        self.product = product
        self.synergy_type = synergy_type
        self.rlabel = [r.label for r in reactants]
        self.plabel = product.label

    def __str__(self):
        return f"{'+'.join(self.rlabel)} → {self.plabel} ({self.synergy_type})"

    def __eq__(self, other):
        if not isinstance(other, ERC_Synergy):
            return False
        return (set(self.rlabel) == set(other.rlabel) and 
                self.plabel == other.plabel and 
                self.synergy_type == other.synergy_type)

    def __hash__(self):
        return hash((tuple(sorted(self.rlabel)), self.plabel, self.synergy_type))


class ERC_Complementarity:
    """
    Represents a complementary relationship between ERCs.
    
    Complementarity types (from theory):
    - Type 1: Requirement reduction - |req(X')| < |req(X)| + |req(E)|
    - Type 2: Requirement change - |req(X')| = |req(X)| + |req(E)| but req(X') ≠ req(X)
    - Type 3: Product expansion - req(X') = req(X) but prod(X') ≠ prod(X)
    """
    def __init__(self, erc1, erc2, comp_type, info=None):
        self.erc1 = erc1
        self.erc2 = erc2
        self.comp_type = comp_type
        self.info = info or {}
        
    def __str__(self):
        return f"{self.erc1.label} ⊕ {self.erc2.label} (Type {self.comp_type})"


# ============================================================================
# 2. HELPER FUNCTIONS
# ============================================================================

def can_interact(erc1, erc2, hierarchy):
    """
    Check if two ERCs can interact (not in containment relationship).
    
    Two ERCs cannot productively interact if one contains the other.
    """
    if erc1 in hierarchy.get_contain(erc2) or erc2 in hierarchy.get_contain(erc1):
        return False
    return True


def has_partial_overlap_with_generators(base_erc, target_erc, RN):
    """
    Check if base ERC can partially contribute to covering target's generators.
    
    Returns True if base_erc covers some (but not all) species in any generator.
    """
    base_closure = base_erc.get_closure_names(RN)
    
    for generator in target_erc.min_generators:
        gen_species = set(species_list_to_names(generator))
        intersection = base_closure & gen_species
        
        # Partial overlap: some intersection but not complete coverage
        if intersection and not gen_species.issubset(base_closure):
            return True
    
    return False


def has_generator_coverage(base1_erc, base2_erc, target_erc, RN):
    """
    Check if combined bases cover at least one complete generator.
    
    Returns True if the union of two ERCs covers a generator that
    neither could cover alone.
    """
    base1_closure = base1_erc.get_closure_names(RN)
    base2_closure = base2_erc.get_closure_names(RN)
    combined_closure = base1_closure | base2_closure
    
    for generator in target_erc.min_generators:
        gen_species = set(species_list_to_names(generator))
        
        # Check if joint closure covers this generator
        if gen_species.issubset(combined_closure):
            # Verify neither base alone can cover it (true synergy)
            if (not gen_species.issubset(base1_closure) and 
                not gen_species.issubset(base2_closure)):
                return True
    
    return False


# ============================================================================
# 3. SYNERGY DETECTION FUNCTIONS
# ============================================================================

def get_basic_synergies(erc1, erc2, hierarchy, RN):
    """
    Find all basic synergies between two ERCs.
    
    A basic synergy occurs when erc1 + erc2 together can cover at least
    one generator of a target ERC that neither could cover alone.
    
    Returns:
        List of ERC_Synergy objects
    """
    if not can_interact(erc1, erc2, hierarchy):
        return []
    
    synergies = []
    
    erc1_closure_names = erc1.get_closure_names(RN)
    erc2_closure_names = erc2.get_closure_names(RN)
    joint_closure_names = erc1_closure_names | erc2_closure_names
    
    # Check each ERC in the hierarchy as a potential target
    for target_erc in hierarchy.ercs:
        # Skip if target is contained by either base ERC
        if target_erc in hierarchy.get_contained(erc1) or target_erc in hierarchy.get_contained(erc2):
            continue
        
        # Check generator coverage conditions
        has_covered_generator = False
        erc1_alone_covers = False
        erc2_alone_covers = False
        
        for generator in target_erc.min_generators:
            gen_species = set(species_list_to_names(generator))
            
            # Check if this generator is fully covered by joint closure
            if gen_species.issubset(joint_closure_names):
                has_covered_generator = True
                
                # Check if either base alone can cover this generator
                if gen_species.issubset(erc1_closure_names):
                    erc1_alone_covers = True
                    break  # Not a synergy if erc1 alone covers a generator
                
                if gen_species.issubset(erc2_closure_names):
                    erc2_alone_covers = True
                    break  # Not a synergy if erc2 alone covers a generator
        
        # This is a synergy if:
        # - At least one generator is covered by the joint closure
        # - Neither base ERC alone can cover any generator
        if has_covered_generator and not erc1_alone_covers and not erc2_alone_covers:
            synergies.append(ERC_Synergy([erc1, erc2], target_erc, "basic"))
    
    return synergies


def get_maximal_synergies(erc1, erc2, hierarchy, RN):
    """
    Find maximal synergies between two ERCs.
    
    A maximal synergy is a basic synergy where the target is not contained
    in any other target of a basic synergy with the same reactants.
    
    Returns:
        List of ERC_Synergy objects with synergy_type="maximal"
    """
    basic_synergies = get_basic_synergies(erc1, erc2, hierarchy, RN)
    if not basic_synergies:
        return []
    
    maximal_synergies = []
    
    for syn in basic_synergies:
        is_maximal = True
        
        # Check if any other synergy with same reactants has a target that contains this one
        for other_syn in basic_synergies:
            if syn != other_syn and set(syn.rlabel) == set(other_syn.rlabel):
                # If syn.product is contained in other_syn.product, syn is not maximal
                if syn.product in hierarchy.get_contained(other_syn.product):
                    is_maximal = False
                    break
        
        if is_maximal:
            maximal_synergies.append(ERC_Synergy(syn.reactants, syn.product, "maximal"))
    
    return maximal_synergies


def get_fundamental_synergies_brute_force(erc1, erc2, hierarchy, RN, maximal_synergies=None, verbose=False):
    """
    Find fundamental synergies between two ERCs using brute force method.
    
    A fundamental synergy is a maximal synergy that cannot be obtained
    by combining contained ERCs (descendants in the hierarchy).
    
    Args:
        erc1, erc2: The two ERCs to check
        hierarchy: The ERC hierarchy
        RN: The reaction network
        maximal_synergies: Pre-computed maximal synergies (optional)
        verbose: Print debug info
        
    Returns:
        List of ERC_Synergy objects with synergy_type="fundamental"
    """
    if maximal_synergies is None:
        maximal_synergies = get_maximal_synergies(erc1, erc2, hierarchy, RN)
    
    if not maximal_synergies:
        return []
    
    fundamental_synergies = []
    
    # Get descendants (contained ERCs) for each ERC
    if hierarchy.graph:
        erc1_successors = set(nx.descendants(hierarchy.graph, erc1.label))
        erc1_successors.add(erc1.label)
        
        erc2_successors = set(nx.descendants(hierarchy.graph, erc2.label))  
        erc2_successors.add(erc2.label)
    else:
        erc1_successors = {erc1.label}
        erc2_successors = {erc2.label}
    
    # Check each maximal synergy
    for syn in maximal_synergies:
        is_fundamental = True
        target_label = syn.plabel
        
        # Check if any pair of descendants can produce the same target
        for other_erc1_label in erc1_successors:
            for other_erc2_label in erc2_successors:
                # Skip if same ERC
                if other_erc1_label == other_erc2_label:
                    continue
                    
                # Skip the original pair
                if (set([other_erc1_label, other_erc2_label]) == 
                    set([erc1.label, erc2.label])):
                    continue
                
                # Get the other ERCs
                other_erc1 = hierarchy.get_erc_by_label(other_erc1_label)
                other_erc2 = hierarchy.get_erc_by_label(other_erc2_label)
                
                # Check if this pair has a synergy to the same target
                other_maximal = get_maximal_synergies(other_erc1, other_erc2, hierarchy, RN)
                
                for other_syn in other_maximal:
                    if other_syn.plabel == target_label:
                        is_fundamental = False
                        break
                
                if not is_fundamental:
                    break
            if not is_fundamental:
                break
        
        if is_fundamental:
            fundamental_synergies.append(ERC_Synergy(syn.reactants, syn.product, "fundamental"))
    
    return fundamental_synergies


def get_all_fundamental_synergies_brute_force(ercs, hierarchy, RN, verbose=False):
    """
    Find all fundamental synergies by exhaustive pairwise checking.
    
    This is a utility function for comprehensive analysis.
    """
    all_fundamental_synergies = []
    
    total_pairs = len(list(combinations(ercs, 2)))
    
    for i, (erc1, erc2) in enumerate(combinations(ercs, 2)):
        if verbose:
            print(f"Checking pair {i+1}/{total_pairs}: {erc1.label} + {erc2.label}")
        
        pair_synergies = get_fundamental_synergies_brute_force(erc1, erc2, hierarchy, RN)
        all_fundamental_synergies.extend(pair_synergies)
        
        if verbose and pair_synergies:
            print(f"  Found {len(pair_synergies)} fundamental synergies")
    
    return all_fundamental_synergies


# ============================================================================
# 4. COMPLEMENTARITY DETECTION FUNCTIONS
# ============================================================================

def get_complementarity(erc1, erc2, hierarchy, RN):
    """
    Detect complementarity relationships between two ERCs.
    
    Checks for three types of complementarity:
    - Type 1: Requirement reduction (fewer required_species together)
    - Type 2: Requirement change (same size but different required_species)
    - Type 3: Product expansion (same required_species but more products)
    
    Args:
        erc1, erc2: The two ERCs to check
        hierarchy: The ERC hierarchy
        RN: The reaction network
        
    Returns:
        List of ERC_Complementarity objects
    """
    if not can_interact(erc1, erc2, hierarchy):
        return []
    
    complementarities = []
    
    # Get required_species and products
    req1 = erc1.get_required_species(RN)   # Correct
    req2 = erc2.get_required_species(RN) 
    prod1 = erc1.get_produced_species(RN)
    prod2 = erc2.get_produced_species(RN)
    
    # Combined required_species (union of what each needs)
    combined_req_union = req1 | req2
    
    # Try to compute the closure of erc1 + erc2
    # This represents what they can produce together
    closure1 = erc1.get_closure_names(RN)
    closure2 = erc2.get_closure_names(RN)
    joint_closure = closure1 | closure2
    
    # Compute actual required_species for the joint system
    # (species needed that aren't produced by either)
    joint_req = combined_req_union - joint_closure
    
    # Type 1: Requirement reduction
    # |req(X')| < |req(X)| + |req(E)|
    if len(joint_req) < len(req1) + len(req2):
        info = {
            'req1': req1,
            'req2': req2,
            'joint_req': joint_req,
            'reduction': len(req1) + len(req2) - len(joint_req)
        }
        complementarities.append(ERC_Complementarity(erc1, erc2, comp_type=1, info=info))
    
    # Type 2: Requirement change (same size but different)
    # |req(X')| = |req(X)| + |req(E)| and req(X') ≠ req(X)
    elif len(joint_req) == len(req1) + len(req2) and joint_req != (req1 | req2):
        info = {
            'req1': req1,
            'req2': req2,
            'joint_req': joint_req,
            'changed_species': (req1 | req2) ^ joint_req
        }
        complementarities.append(ERC_Complementarity(erc1, erc2, comp_type=2, info=info))
    
    # Type 3: Product expansion
    # req(X') = req(X) and prod(X') ≠ prod(X)
    # Check if required_species stay the same but products expand
    joint_prod = prod1 | prod2
    if joint_req == req1 and joint_prod != prod1:
        info = {
            'req1': req1,
            'prod1': prod1,
            'joint_prod': joint_prod,
            'new_products': joint_prod - prod1
        }
        complementarities.append(ERC_Complementarity(erc1, erc2, comp_type=3, info=info))
    
    return complementarities


# ============================================================================
# 5. ERC_SORN CLASS
# ============================================================================

class ERC_SORN:
    """
    Second Order Reaction Network (SORN).
    
    Pre-computes and stores all fundamental synergies and complementarities
    between ERCs for efficient querying.
    
    This is the main data structure for analyzing productive relationships
    between ERCs.
    """
    
    def __init__(self, hierarchy, RN):
        """
        Initialize ERC_SORN by computing all productive relationships.
        
        Args:
            hierarchy: ERC_Hierarchy object
            RN: Reaction network
        """
        self.hierarchy = hierarchy
        self.RN = RN
        self.ercs = hierarchy.ercs

        # Storage dictionaries
        self._synergies = {}
        self._complementarities = {}
        self._erc_to_synergies = defaultdict(list)
        self._erc_to_complementarities = defaultdict(list)
        self._productive_partners = defaultdict(set)

        # Statistics tracking
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
        """
        Build SORN by computing synergies and complementarities for all ERC pairs.
        
        This is called automatically during initialization.
        """
        start = time.time()
        print(f"Building ERC_SORN for {len(self.ercs)} ERCs...")

        # Check all pairs of ERCs
        for erc1, erc2 in combinations(self.ercs, 2):
            self.computation_stats['total_pairs_checked'] += 1

            # Skip pairs that can't interact
            if not can_interact(erc1, erc2, self.hierarchy):
                continue

            pair_has_productive = False

            # Compute fundamental synergies
            synergies = get_fundamental_synergies_brute_force(
                erc1, erc2, self.hierarchy, self.RN
            )
            
            if synergies:
                self._store_synergies(erc1.label, erc2.label, synergies)
                pair_has_productive = True
                self.computation_stats['total_synergies'] += len(synergies)
                self.computation_stats['synergistic_pairs'] += 1

            # Compute complementarities
            comps = get_complementarity(erc1, erc2, self.hierarchy, self.RN)
            
            if comps:
                self._store_complementarities(erc1.label, erc2.label, comps)
                pair_has_productive = True
                self.computation_stats['total_complementarities'] += len(comps)
                self.computation_stats['complementary_pairs'] += 1

            # Track productive partnerships
            if pair_has_productive:
                self._productive_partners[erc1.label].add(erc2.label)
                self._productive_partners[erc2.label].add(erc1.label)
                self.computation_stats['productive_pairs'] += 1

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
        self._synergies[key2] = synergies
        
        self._erc_to_synergies[erc1_label].append((erc2_label, synergies))
        self._erc_to_synergies[erc2_label].append((erc1_label, synergies))
    
    def _store_complementarities(self, erc1_label, erc2_label, complementarities):
        """Store complementarities with bidirectional lookup."""
        key1 = (erc1_label, erc2_label)
        key2 = (erc2_label, erc1_label)
        
        self._complementarities[key1] = complementarities
        self._complementarities[key2] = complementarities
        
        self._erc_to_complementarities[erc1_label].append((erc2_label, complementarities))
        self._erc_to_complementarities[erc2_label].append((erc1_label, complementarities))
    
    # ========================================================================
    # Query methods
    # ========================================================================
    
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
        Get all possible productive extensions for a generator.
        
        Given a set of ERCs already in a generator, find all ERCs that
        could be added via synergy or complementarity.
        """
        current_erc_set = set(erc_labels_in_generator)
        candidate_extensions = {}
        
        for erc_label in erc_labels_in_generator:
            # Check for synergistic partners
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
            
            # Check for complementary partners
            for partner_label, comps in self._erc_to_complementarities.get(erc_label, []):
                if partner_label not in current_erc_set:
                    if partner_label not in candidate_extensions:
                        candidate_extensions[partner_label] = []
                    for comp in comps:
                        candidate_extensions[partner_label].append(('complementarity', {
                            'comp_type': comp.comp_type,
                            'comp_object': comp,
                            'with_erc': erc_label
                        }))
        
        # Format extensions
        extensions = []
        for candidate_label, extension_list in candidate_extensions.items():
            step_type, step_details = extension_list[0]
            extensions.append((candidate_label, step_type, step_details))
        
        return extensions
    
    def get_statistics(self):
        """Get statistics about the SORN."""
        return self.computation_stats.copy()
    
    def __repr__(self):
        return (f"ERC_SORN({len(self.ercs)} ERCs, "
                f"{self.computation_stats['synergistic_pairs']} synergistic pairs, "
                f"{self.computation_stats['complementary_pairs']} complementary pairs)")


# ============================================================================
# 6. FACTORY FUNCTIONS
# ============================================================================

def build_erc_sorn(hierarchy, RN):
    """
    Factory function to build an ERC_SORN.
    
    Args:
        hierarchy: ERC_Hierarchy object
        RN: Reaction network
        
    Returns:
        ERC_SORN object with all relationships computed
    """
    return ERC_SORN(hierarchy, RN)