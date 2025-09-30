#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Core Synergy and Complementarity Classes and Functions
"""
import time
from itertools import combinations
import networkx as nx
from collections import defaultdict
from pyCOT.ERC_Hierarchy import ERC, ERC_Hierarchy, closure, species_list_to_names

# ============================================================================
# CORE CLASSES
# ============================================================================

class ERC_Synergy:
    """Base synergy class"""
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
    """Base complementarity class"""
    def __init__(self, erc1, erc2, comp_type, info=None):
        self.erc1 = erc1
        self.erc2 = erc2
        self.comp_type = comp_type
        self.info = info or {}
        
    def __str__(self):
        return f"{self.erc1.label} ⊕ {self.erc2.label} (Type {self.comp_type})"

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def can_interact(erc1, erc2, hierarchy):
    """Check if ERCs can interact"""
    if erc1 in hierarchy.get_contain(erc2) or erc2 in hierarchy.get_contain(erc1):
        return False
    return True

def has_partial_overlap_with_generators(base_erc, target_erc, RN):
    """Check if base can contribute to covering target's generators."""
    base_closure = base_erc.get_closure_names(RN)
    
    for generator in target_erc.min_generators:
        gen_species = set(species_list_to_names(generator))
        intersection = base_closure & gen_species
        
        if intersection and not gen_species.issubset(base_closure):
            return True
    
    return False

def has_generator_coverage(base1_erc, base2_erc, target_erc, RN):
    """Check if combined bases cover at least one complete generator."""
    base1_closure = base1_erc.get_closure_names(RN)
    base2_closure = base2_erc.get_closure_names(RN)
    combined_closure = base1_closure | base2_closure
    
    for generator in target_erc.min_generators:
        gen_species = set(species_list_to_names(generator))
        
        if gen_species.issubset(combined_closure):
            if (not gen_species.issubset(base1_closure) and 
                not gen_species.issubset(base2_closure)):
                return True
    
    return False

# ============================================================================
# SYNERGY DETECTION
# ============================================================================

def get_basic_synergies(erc1, erc2, hierarchy, RN):
    """Generate all synergetic pairs by brute force"""
    if not can_interact(erc1, erc2, hierarchy):
        return []
    
    synergies = []
    
    erc1_closure = erc1.get_closure(RN)
    erc2_closure = erc2.get_closure(RN)
    
    rn1 = erc1.get_reacs(RN)
    rn2 = erc2.get_reacs(RN)
    
    rn1_names = set()
    rn2_names = set()
    
    for r in rn1:
        if hasattr(r, 'name'):
            rn1_names.add(r.name())
        else:
            rn1_names.add(str(r))
    
    for r in rn2:
        if hasattr(r, 'name'):
            rn2_names.add(r.name())
        else:
            rn2_names.add(str(r))
    
    union_reac_names = rn1_names.union(rn2_names)
    
    union_species = list(set(erc1_closure).union(set(erc2_closure)))
    joint_closure = closure(RN, union_species)
    
    joint_closure_reacs = RN.get_reactions_from_species(joint_closure)
    
    joint_closure_names = set()
    for r in joint_closure_reacs:
        if hasattr(r, 'name'):
            joint_closure_names.add(r.name())
        else:
            joint_closure_names.add(str(r))
    
    if len(joint_closure_names) > len(union_reac_names):   
        novel_reac_names = joint_closure_names - union_reac_names
        
        target_ercs_found = set()
        
        for r in joint_closure_reacs:
            r_name = r.name() if hasattr(r, 'name') else str(r)
            if r_name in novel_reac_names:
                syn_erc = hierarchy.get_erc_from_reaction(RN, hierarchy, r)
                if syn_erc is None:
                    continue
                
                if syn_erc.label not in target_ercs_found:
                    target_ercs_found.add(syn_erc.label)
                    synergies.append(ERC_Synergy([erc1, erc2], syn_erc, "regular"))
    
    return synergies

def get_maximal_synergies(erc1, erc2, hierarchy, RN):
    """Find maximal synergies between two ERCs"""
    basic_synergies = get_basic_synergies(erc1, erc2, hierarchy, RN)
    if not basic_synergies:
        return []
    
    maximal_synergies = []
    
    for syn in basic_synergies:
        is_maximal = True
        
        for other_syn in basic_synergies:
            if syn != other_syn and set(syn.rlabel) == set(other_syn.rlabel):
                if syn.product in hierarchy.get_contained(other_syn.product):
                    is_maximal = False
                    break
        
        if is_maximal:
            maximal_synergies.append(ERC_Synergy(syn.reactants, syn.product, "maximal"))
    
    return maximal_synergies

def get_fundamental_synergies_brute_force(erc1, erc2, hierarchy, RN, maximal_synergies=None, verbose=False):
    """Find fundamental synergies between two ERCs"""
    if maximal_synergies == None:
        maximal_synergies = get_maximal_synergies(erc1, erc2, hierarchy, RN)
    if not maximal_synergies:
        return []
    
    fundamental_synergies = []
    
    if hierarchy.graph:
        erc1_successors = set(nx.descendants(hierarchy.graph, erc1.label))
        erc1_successors.add(erc1.label)
        
        erc2_successors = set(nx.descendants(hierarchy.graph, erc2.label))  
        erc2_successors.add(erc2.label)
    else:
        erc1_successors = {erc1.label}
        erc2_successors = {erc2.label}
    
    for syn in maximal_synergies:
        is_fundamental = True
        target_label = syn.plabel
        
        for other_erc1_label in erc1_successors:
            for other_erc2_label in erc2_successors:
                if other_erc1_label == other_erc2_label:
                    continue
                    
                if (set([other_erc1_label, other_erc2_label]) == 
                    set([erc1.label, erc2.label])):
                    continue
                
                other_erc1 = hierarchy.get_erc_by_label(other_erc1_label)
                other_erc2 = hierarchy.get_erc_by_label(other_erc2_label)
                
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
    """Find all fundamental synergies by exhaustive pairwise checking."""
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
# COMPLEMENTARITY DETECTION
# ============================================================================

def is_complementary_type1(erc1, erc2, req1, req2, joint_req, joint_consumed,
                           joint_produced, prod1, prod2, hierarchy, RN):
    """Check Type 1 complementarity: |reqs(X')|<|reqs(X)|+|reqs(E)|"""
    reduction = len(req1) + len(req2) - len(joint_req)
        
    if reduction > 0:
        return True, {
            'req1': req1,
            'req2': req2, 
            'joint_req': joint_req,
            'joint_consumed': joint_consumed,
            'joint_produced': joint_produced,
            'reduction': reduction,
            'satisfied_by_1': req2 & prod1,
            'satisfied_by_2': req1 & prod2,
            'satisfied_by_synergy': (req1 | req2) & joint_produced - (prod1 | prod2)
        }
    else:
        return False, {}

def is_complementary_type2(erc1, erc2, req1, req2, joint_req, joint_consumed,
                           joint_produced, prod1, prod2, hierarchy, RN):
    """Check Type 2 complementarity"""
    if not can_interact(erc1, erc2, hierarchy):
        return False, {}
        
    try:
        if len(joint_req) == len(req1) + len(req2) and joint_req != req1:
            return True, {
                'req1': req1,
                'req2': req2,
                'joint_req': joint_req,
                'joint_consumed': joint_consumed,
                'joint_produced': joint_produced,
                'requirement_change': joint_req - req1,
                'requirement_shift': (req1 | req2) - joint_req
            }
    except Exception as e:
        print(f"⚠️  Warning: Error in type2 complementarity for {erc1.label}, {erc2.label}: {e}")
    
    return False, {}

def is_complementary_type3(erc1, erc2, req1, req2, joint_req, joint_consumed,
                           joint_produced, prod1, prod2, hierarchy, RN):
    """Check Type 3 complementarity"""
    if not can_interact(erc1, erc2, hierarchy):
        return False, {}
        
    try:
        novel_products = set(joint_produced) - (prod1 | prod2)
        total_products = joint_produced
        
        if joint_req == req1 and total_products != prod1 | prod2:
            return True, {
                'req1': req1,
                'joint_req': joint_req,
                'prod1': prod1,
                'joint_produced': joint_produced,
                'total_products': total_products,
                'novel_products': joint_produced - (prod1 | prod2),
                'has_synergy': len(joint_produced - (prod1 | prod2)) > 0
            }
    except Exception as e:
        print(f"⚠️  Warning: Error in type3 complementarity for {erc1.label}, {erc2.label}: {e}")
    
    return False, {}

def get_complementarity(erc1, erc2, hierarchy, RN):
    """Get all complementarity relationships between two ERCs."""
    if not can_interact(erc1, erc2, hierarchy):
        return []
    
    complementarities = []
    
    req1 = erc1.get_required_species(RN)
    req2 = erc2.get_required_species(RN)
    prod1 = erc1.get_produced_species(RN)
    prod2 = erc2.get_produced_species(RN)
    
    erc1_closure = erc1.get_closure(RN)
    erc2_closure = erc2.get_closure(RN)
    union_species = list(set(erc1_closure).union(set(erc2_closure)))
    joint_closure = closure(RN, union_species)
    
    joint_closure_reacs = RN.get_reactions_from_species(joint_closure)
    
    joint_consumed = set()
    joint_produced = set()
    
    for reaction in joint_closure_reacs:
        for edge in reaction.edges:
            if edge.type == "reactant":
                joint_consumed.add(edge.species_name)
            elif edge.type == "product":
                joint_produced.add(edge.species_name)
    
    joint_req = joint_consumed - joint_produced

    is_type1, info1 = is_complementary_type1(erc1, erc2, req1, req2, joint_req, joint_consumed,
                           joint_produced, prod1, prod2, hierarchy, RN)
    if is_type1:
        complementarities.append(ERC_Complementarity(erc1, erc2, 1, info1))
    
    is_type2, info2 = is_complementary_type2(erc1, erc2, req1, req2, joint_req, joint_consumed,
                           joint_produced, prod1, prod2, hierarchy, RN)
    if is_type2:
        complementarities.append(ERC_Complementarity(erc1, erc2, 2, info2))
    
    is_type3, info3 = is_complementary_type3(erc1, erc2, req1, req2, joint_req, joint_consumed,
                           joint_produced, prod1, prod2, hierarchy, RN)
    if is_type3:
        complementarities.append(ERC_Complementarity(erc1, erc2, 3, info3))
    
    return complementarities

# ============================================================================
# ERC_SORN CLASS
# ============================================================================

class ERC_SORN:
    """
    Second Order Reaction Network: Pre-computed network of relationships between ERCs.
    """
    
    def __init__(self, hierarchy, RN):
        """Initialize ERC_SORN by pre-computing all productive relationships."""
        self.hierarchy = hierarchy
        self.RN = RN
        self.ercs = hierarchy.ercs

        self._synergies = {}
        self._complementarities = {}
        self._erc_to_synergies = defaultdict(list)
        self._erc_to_complementarities = defaultdict(list)
        self._productive_partners = defaultdict(set)

        self.computation_stats = {
            'total_pairs_checked': 0,
            'productive_pairs': 0,
            'total_synergies': 0,
            'total_complementarities': 0,
            'synergistic_pairs': 0,
            'complementary_pairs': 0,
            'build_time': 0.0
        }
        
        self._build_sorn()
    
    def _build_sorn(self):
        """Build SORN with single-pass computation of synergies and complementarities."""
        start = time.time()
        print(f"Building ERC_SORN for {len(self.ercs)} ERCs...")

        for erc1, erc2 in combinations(self.ercs, 2):
            self.computation_stats['total_pairs_checked'] += 1

            if not can_interact(erc1, erc2, self.hierarchy):
                continue

            key = tuple(sorted((erc1.label, erc2.label)))
            pair_has_productive = False

            synergies = get_fundamental_synergies_brute_force(
                erc1, erc2, self.hierarchy, self.RN
            )
            
            if synergies:
                self._store_synergies(erc1.label, erc2.label, synergies)
                pair_has_productive = True
                self.computation_stats['total_synergies'] += len(synergies)
                self.computation_stats['synergistic_pairs'] += 1

            comps = get_complementarity(erc1, erc2, self.hierarchy, self.RN)
            
            if comps:
                self._store_complementarities(erc1.label, erc2.label, comps)
                pair_has_productive = True
                self.computation_stats['total_complementarities'] += len(comps)
                self.computation_stats['complementary_pairs'] += 1

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
    
    # Query methods
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
        """Get all possible productive extensions for a generator."""
        current_erc_set = set(erc_labels_in_generator)
        candidate_extensions = {}
        
        for erc_label in erc_labels_in_generator:
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
# FACTORY FUNCTIONS
# ============================================================================

def build_erc_sorn(hierarchy, RN):
    """Factory function to build an ERC_SORN."""
    return ERC_SORN(hierarchy, RN)