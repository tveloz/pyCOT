#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Core Synergy and Complementarity Classes and Functions

Pure mathematical definitions and detection functions for ERC synergies and 
complementarities based on the theoretical framework by Tomas Veloz et al.

This module contains ONLY core synergy and complementarity functions and classes.
NO analysis, statistics, plotting, or interaction pattern functions.

Classes:
    ERC_Synergy: Represents a synergy between ERCs
    ERC_Complementarity: Represents a complementarity between ERCs
    ERC_SORN: Pre-computed second-order reaction network for efficient lookup



Functions:
    get_basic_synergies: Detect basic synergies between ERCs
    get_maximal_synergies: Detect maximal synergies between ERCs  
    get_fundamental_synergies: Detect fundamental synergies between ERCs
    get_complementarity: Detect complementarity relationships between ERCs
    is_complementary_type1/2/3: Check specific complementarity types
    can_interact: Check if two ERCs can have productive interactions
    build_erc_sorn: Factory function to create ERC_SORN
    get_all_productive_pairs: Get all ERC pairs with productive relationships

Author: Based on theoretical work by Tomas Veloz et al.
"""

from itertools import combinations
from pyCOT.ERC_Hierarchy import ERC, ERC_Hierarchy, closure
from collections import defaultdict  # (if not already imported)

# ============================================================================
# CORE CLASSES
# ============================================================================

class ERC_Synergy:
    """
    Represents a synergy between ERCs.
    
    A synergy occurs when combining ERCs activates novel reactions that 
    cannot be activated by the individual ERCs alone.
    """
    
    def __init__(self, reactants, product, synergy_type="regular"):
        """
        Initialize a Synergy object.
        
        Parameters
        ----------
        reactants : list of ERC
            List of ERC objects representing the reactants
        product : ERC
            ERC object representing the product
        synergy_type : str
            Type of synergy ("regular", "maximal", "fundamental")
        """
        if not isinstance(reactants, list) or not all(isinstance(r, ERC) for r in reactants):
            raise ValueError("Reactants must be a list of ERC objects.")
        if not isinstance(product, ERC):
            raise ValueError("Product must be an ERC object.")
        
        self.reactants = reactants
        self.product = product
        self.rlabel = [r.label for r in reactants]
        self.plabel = product.label
        self.synergy_type = synergy_type

    def __repr__(self):
        """String representation of the Synergy object."""
        reactant_labels = "+".join(self.rlabel)
        return f"{self.synergy_type.capitalize()}Synergy({reactant_labels}→{self.plabel})"

class ERC_Complementarity:
    """
    Represents a complementarity between ERCs.
    
    Complementarity occurs when combining ERCs improves productive capacity
    through requirement reduction, requirement change, or product expansion.
    """
    
    def __init__(self, erc1, erc2, comp_type, reduction_info=None):
        """
        Initialize a Complementarity object.
        
        Parameters
        ----------
        erc1, erc2 : ERC
            ERC objects that are complementary
        comp_type : int
            Type of complementarity (1, 2, or 3)
        reduction_info : dict, optional
            Dictionary with details about the reduction/change
        """
        self.erc1 = erc1
        self.erc2 = erc2
        self.comp_type = comp_type
        self.reduction_info = reduction_info or {}
        
    def __repr__(self):
        """String representation of the Complementarity object."""
        return f"Complementarity_Type{self.comp_type}({self.erc1.label}⇔{self.erc2.label})"

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def can_interact(erc1, erc2, hierarchy):
    """
    Check if two ERCs can have synergy or complementarity.
    
    Returns False if one contains the other (efficiency optimization).
    
    Parameters
    ----------
    erc1, erc2 : ERC
        ERCs to check for interaction potential
    hierarchy : ERC_Hierarchy
        Hierarchy containing the ERCs
        
    Returns
    -------
    bool
        True if ERCs can interact, False if one contains the other
    """
    # Check if one ERC contains the other using cached containment relationships
    if erc1 in hierarchy.get_contain(erc2) or erc2 in hierarchy.get_contain(erc1):
        return False
    return True

# ============================================================================
# SYNERGY DETECTION FUNCTIONS
# ============================================================================

def get_basic_synergies(erc1, erc2, hierarchy, RN):
    """
    Get all basic synergies between two ERCs using cached closures.
    
    Parameters
    ----------
    erc1, erc2 : ERC
        ERCs to check for synergies
    hierarchy : ERC_Hierarchy
        Hierarchy containing the ERCs
    RN : ReactionNetwork
        The reaction network
        
    Returns
    -------
    list of ERC_Synergy
        List of basic synergies found
    """
    if not can_interact(erc1, erc2, hierarchy):
        return []
    
    synergies = []
    
    # Use cached closures - NO recomputation needed
    erc1_closure = erc1.get_closure(RN)
    erc2_closure = erc2.get_closure(RN)
    
    # Get cached reactions
    rn1 = erc1.get_reacs(RN)
    rn2 = erc2.get_reacs(RN)
    
    # Convert reactions to hashable format (reaction names/labels)
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
    
    # Compute joint closure using the union of cached closures
    union_species = list(set(erc1_closure).union(set(erc2_closure)))
    joint_closure = closure(RN, union_species)
    
    joint_closure_reacs = RN.get_reactions_from_species(joint_closure)
    
    # Convert joint closure reactions to names
    joint_closure_names = set()
    for r in joint_closure_reacs:
        if hasattr(r, 'name'):
            joint_closure_names.add(r.name())
        else:
            joint_closure_names.add(str(r))
    
    if len(joint_closure_names) > len(union_reac_names):   
        novel_reac_names = joint_closure_names - union_reac_names
        for r in joint_closure_reacs:
            r_name = r.name() if hasattr(r, 'name') else str(r)
            if r_name in novel_reac_names:
                syn_erc = hierarchy.get_erc_from_reaction(RN, hierarchy, r)
                if syn_erc is None:
                    continue
                    
                # Check for redundancy
                add = True
                for s in synergies:
                    if syn_erc in hierarchy.get_contained(s.product):
                        add = False
                        break
                    if set(s.rlabel) == set([erc1.label, erc2.label]) and s.plabel == syn_erc.label:  
                        add = False
                        break
                if add:
                    synergies.append(ERC_Synergy([erc1, erc2], syn_erc, "regular"))
    
    return synergies

def get_maximal_synergies(erc1, erc2, hierarchy, RN):
    """
    Get maximal synergies between two ERCs.
    
    A synergy E1+E2→E3 is maximal if there's no other synergy E1+E2→E3' 
    such that E3'→E3 (E3' contains E3).
    
    Parameters
    ----------
    erc1, erc2 : ERC
        ERCs to check for maximal synergies
    hierarchy : ERC_Hierarchy
        Hierarchy containing the ERCs
    RN : ReactionNetwork
        The reaction network
        
    Returns
    -------
    list of ERC_Synergy
        List of maximal synergies found
    """
    basic_synergies = get_basic_synergies(erc1, erc2, hierarchy, RN)
    if not basic_synergies:
        return []
    
    maximal_synergies = []
    
    for syn in basic_synergies:
        is_maximal = True
        
        # Check if there's another synergy with the same reactants but larger product
        for other_syn in basic_synergies:
            if syn != other_syn and set(syn.rlabel) == set(other_syn.rlabel):
                # Check if other_syn.product contains syn.product using cached containment
                if syn.product in hierarchy.get_contained(other_syn.product):
                    is_maximal = False
                    break
        
        if is_maximal:
            maximal_synergies.append(ERC_Synergy(syn.reactants, syn.product, "maximal"))
    
    return maximal_synergies

def get_fundamental_synergies(erc1, erc2, hierarchy, RN):
    """
    Get fundamental synergies between two ERCs.
    
    A maximal synergy E1+E2→E3 is fundamental if there's no synergy E1'+E2→E3 
    or E1+E2'→E3 with Ei→Ei' for i=1,2.
    
    Parameters
    ----------
    erc1, erc2 : ERC
        ERCs to check for fundamental synergies
    hierarchy : ERC_Hierarchy
        Hierarchy containing the ERCs
    RN : ReactionNetwork
        The reaction network
        
    Returns
    -------
    list of ERC_Synergy
        List of fundamental synergies found
    """
    maximal_synergies = get_maximal_synergies(erc1, erc2, hierarchy, RN)
    if not maximal_synergies:
        return []
    
    fundamental_synergies = []
    
    for syn in maximal_synergies:
        is_fundamental = True
        
        # Check all other ERC pairs to see if there's a more fundamental synergy
        for other_erc1 in hierarchy.ercs:
            for other_erc2 in hierarchy.ercs:
                if other_erc1 == other_erc2:
                    continue
                    
                # Skip if it's the same pair
                if set([other_erc1.label, other_erc2.label]) == set([erc1.label, erc2.label]):
                    continue
                
                other_maximal = get_maximal_synergies(other_erc1, other_erc2, hierarchy, RN)
                
                for other_syn in other_maximal:
                    # Check if other synergy produces the same result with more fundamental reactants
                    if other_syn.plabel == syn.plabel:
                        # Check if other reactants are more fundamental using cached containment
                        if ((other_erc1 in hierarchy.get_contained(erc1) and other_erc2.label == erc2.label) or
                            (other_erc1.label == erc1.label and other_erc2 in hierarchy.get_contained(erc2)) or
                            (other_erc1 in hierarchy.get_contained(erc2) and other_erc2.label == erc1.label) or
                            (other_erc1.label == erc2.label and other_erc2 in hierarchy.get_contained(erc1))):
                            is_fundamental = False
                            break
                
                if not is_fundamental:
                    break
            if not is_fundamental:
                break
        
        if is_fundamental:
            fundamental_synergies.append(ERC_Synergy(syn.reactants, syn.product, "fundamental"))
    
    return fundamental_synergies

# ============================================================================
# COMPLEMENTARITY DETECTION FUNCTIONS
# ============================================================================

def is_complementary_type1(erc1, erc2, hierarchy, RN):
    """
    Check Type 1 complementarity: |reqs(X')|<|reqs(X)|+|reqs(E)|
    
    CORRECTED: Uses closure operator (∨) to calculate joint requirements.
    This means the combination reduces the total required species.
    
    Parameters
    ----------
    erc1, erc2 : ERC
        ERCs to check for Type 1 complementarity
    hierarchy : ERC_Hierarchy
        Hierarchy containing the ERCs
    RN : ReactionNetwork
        The reaction network
        
    Returns
    -------
    tuple
        (is_complementary: bool, info: dict)
    """
    if not can_interact(erc1, erc2, hierarchy):
        return False, {}
    
    try:
        # Get individual ERC requirements and products (cached)
        req1 = erc1.get_required_species(RN)
        req2 = erc2.get_required_species(RN)
        prod1 = erc1.get_produced_species(RN)
        prod2 = erc2.get_produced_species(RN)
        
        # CRITICAL FIX: Calculate joint closure using ∨ operator (like synergies)
        erc1_closure = erc1.get_closure(RN)
        erc2_closure = erc2.get_closure(RN)
        union_species = list(set(erc1_closure).union(set(erc2_closure)))
        joint_closure = closure(RN, union_species)  # ∨ operator!
        
        # Get reactions from the joint closure
        joint_closure_reacs = RN.get_reactions_from_species(joint_closure)
        
        # Calculate joint requirements from the FULL closure (including synergies)
        joint_consumed = set()
        joint_produced = set()
        
        for reaction in joint_closure_reacs:
            # Get reactants (support) from reaction edges
            for edge in reaction.edges:
                if edge.type == "reactant":
                    joint_consumed.add(edge.species_name)
                elif edge.type == "product":
                    joint_produced.add(edge.species_name)
        
        # Joint requirements = what's consumed but not produced in the joint closure
        joint_req = joint_consumed - joint_produced
        
        # Type 1: Reduction in total requirements
        reduction = len(req1) + len(req2) - len(joint_req)
        
        if reduction > 0:
            return True, {
                'req1': req1,
                'req2': req2, 
                'joint_req': joint_req,
                'joint_consumed': joint_consumed,
                'joint_produced': joint_produced,
                'reduction': reduction,
                'satisfied_by_1': req2 & prod1,  # What erc2 needs that erc1 produces
                'satisfied_by_2': req1 & prod2,  # What erc1 needs that erc2 produces
                'satisfied_by_synergy': (req1 | req2) & joint_produced - (prod1 | prod2)  # Novel satisfaction
            }
    except Exception as e:
        print(f"⚠️  Warning: Error in type1 complementarity for {erc1.label}, {erc2.label}: {e}")
    
    return False, {}

def is_complementary_type2(erc1, erc2, hierarchy, RN):
    """
    Check Type 2 complementarity: |reqs(X')|=|reqs(X)|+|reqs(E)| and reqs(X')≠reqs(X)
    
    CORRECTED: Uses closure operator (∨) to calculate joint requirements.
    No reduction in total requirements but different requirement set.
    
    Parameters
    ----------
    erc1, erc2 : ERC
        ERCs to check for Type 2 complementarity
    hierarchy : ERC_Hierarchy
        Hierarchy containing the ERCs
    RN : ReactionNetwork
        The reaction network
        
    Returns
    -------
    tuple
        (is_complementary: bool, info: dict)
    """
    if not can_interact(erc1, erc2, hierarchy):
        return False, {}
        
    try:
        # Get individual ERC requirements and products (cached)
        req1 = erc1.get_required_species(RN)
        req2 = erc2.get_required_species(RN)
        prod1 = erc1.get_produced_species(RN)
        prod2 = erc2.get_produced_species(RN)
        
        # CRITICAL FIX: Calculate joint closure using ∨ operator
        erc1_closure = erc1.get_closure(RN)
        erc2_closure = erc2.get_closure(RN)
        union_species = list(set(erc1_closure).union(set(erc2_closure)))
        joint_closure = closure(RN, union_species)  # ∨ operator!
        
        # Get reactions from the joint closure
        joint_closure_reacs = RN.get_reactions_from_species(joint_closure)
        
        # Calculate joint requirements from the FULL closure
        joint_consumed = set()
        joint_produced = set()
        
        for reaction in joint_closure_reacs:
            for edge in reaction.edges:
                if edge.type == "reactant":
                    joint_consumed.add(edge.species_name)
                elif edge.type == "product":
                    joint_produced.add(edge.species_name)
        
        joint_req = joint_consumed - joint_produced
        
        # Type 2: Same total requirements but different set
        if len(joint_req) == len(req1) + len(req2) and joint_req != req1:
            return True, {
                'req1': req1,
                'req2': req2,
                'joint_req': joint_req,
                'joint_consumed': joint_consumed,
                'joint_produced': joint_produced,
                'requirement_change': joint_req - req1,
                'requirement_shift': (req1 | req2) - joint_req  # What was eliminated
            }
    except Exception as e:
        print(f"⚠️  Warning: Error in type2 complementarity for {erc1.label}, {erc2.label}: {e}")
    
    return False, {}

def is_complementary_type3(erc1, erc2, hierarchy, RN):
    """
    Check Type 3 complementarity: reqs(X')=reqs(X) and prods(X')≠prods(X)
    
    CORRECTED: Uses closure operator (∨) to calculate joint requirements and products.
    Same requirements but different products (implies synergy).
    
    Parameters
    ----------
    erc1, erc2 : ERC
        ERCs to check for Type 3 complementarity
    hierarchy : ERC_Hierarchy
        Hierarchy containing the ERCs
    RN : ReactionNetwork
        The reaction network
        
    Returns
    -------
    tuple
        (is_complementary: bool, info: dict)
    """
    if not can_interact(erc1, erc2, hierarchy):
        return False, {}
        
    try:
        # Get individual ERC requirements and products (cached)
        req1 = erc1.get_required_species(RN)
        req2 = erc2.get_required_species(RN)
        prod1 = erc1.get_produced_species(RN)
        prod2 = erc2.get_produced_species(RN)
        
        # CRITICAL FIX: Calculate joint closure using ∨ operator
        erc1_closure = erc1.get_closure(RN)
        erc2_closure = erc2.get_closure(RN)
        union_species = list(set(erc1_closure).union(set(erc2_closure)))
        joint_closure = closure(RN, union_species)  # ∨ operator!
        
        # Get reactions from the joint closure
        joint_closure_reacs = RN.get_reactions_from_species(joint_closure)
        
        # Calculate joint requirements and products from the FULL closure
        joint_consumed = set()
        joint_produced = set()
        
        for reaction in joint_closure_reacs:
            for edge in reaction.edges:
                if edge.type == "reactant":
                    joint_consumed.add(edge.species_name)
                elif edge.type == "product":
                    joint_produced.add(edge.species_name)
        
        joint_req = joint_consumed - joint_produced
        
        # Check if there are novel products from synergies
        basic_synergies = get_basic_synergies(erc1, erc2, hierarchy, RN)
        novel_products_from_synergies = set()
        
        for syn in basic_synergies:
            syn_prod = syn.product.get_produced_species(RN)
            novel_products_from_synergies |= syn_prod
        
        # Total products including synergistic products
        total_products = joint_produced
        
        # Type 3: Same requirements as one ERC but different products
        if joint_req == req1 and total_products != prod1 | prod2:
            return True, {
                'req1': req1,
                'joint_req': joint_req,
                'prod1': prod1,
                'joint_produced': joint_produced,
                'total_products': total_products,
                'novel_products_from_synergies': novel_products_from_synergies,
                'novel_products_direct': joint_produced - (prod1 | prod2),  # Direct novel products
                'has_synergy': len(basic_synergies) > 0
            }
    except Exception as e:
        print(f"⚠️  Warning: Error in type3 complementarity for {erc1.label}, {erc2.label}: {e}")
    
    return False, {}

def get_complementarity(erc1, erc2, hierarchy, RN):
    """
    Get all complementarity relationships between two ERCs.
    
    Parameters
    ----------
    erc1, erc2 : ERC
        ERCs to check for complementarity
    hierarchy : ERC_Hierarchy
        Hierarchy containing the ERCs
    RN : ReactionNetwork
        The reaction network
        
    Returns
    -------
    list of ERC_Complementarity
        List of complementarity relationships found
    """
    if not can_interact(erc1, erc2, hierarchy):
        return []
    
    complementarities = []
    
    # Check Type 1
    is_type1, info1 = is_complementary_type1(erc1, erc2, hierarchy, RN)
    if is_type1:
        complementarities.append(ERC_Complementarity(erc1, erc2, 1, info1))
    
    # Check Type 2
    is_type2, info2 = is_complementary_type2(erc1, erc2, hierarchy, RN)
    if is_type2:
        complementarities.append(ERC_Complementarity(erc1, erc2, 2, info2))
    
    # Check Type 3
    is_type3, info3 = is_complementary_type3(erc1, erc2, hierarchy, RN)
    if is_type3:
        complementarities.append(ERC_Complementarity(erc1, erc2, 3, info3))
    
    return complementarities

# ============================================================================
"""

ERC_SORN (Second Order Reaction Network) Class
The SORN represents the network of productive relationships BETWEEN ERCs,
enabling O(1) lookup of synergies and complementarities instead of expensive
recalculation during generator construction.


"""

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
        self._synergies = {}  # (erc1_label, erc2_label) -> list of ERC_Synergy
        self._complementarities = {}  # (erc1_label, erc2_label) -> list of ERC_Complementarity
        
        # Additional indexes for efficient queries
        self._erc_to_synergies = defaultdict(list)  # erc_label -> list of (partner_label, synergies)
        self._erc_to_complementarities = defaultdict(list)  # erc_label -> list of (partner_label, complementarities)
        self._productive_partners = defaultdict(set)  # erc_label -> set of partner labels with ANY productive relationship
        
        # Performance tracking
        self.computation_stats = {
            'total_pairs_checked': 0,
            'synergistic_pairs': 0,
            'complementary_pairs': 0,
            'total_synergies': 0,
            'total_complementarities': 0
        }
        
        # Build the SORN
        self._build_sorn()
    
    def _build_sorn(self):
        """Pre-compute all synergies and complementarities between ERC pairs."""
        print(f"Building ERC_SORN for {len(self.ercs)} ERCs...")
        
        # Check all pairs of ERCs
        for erc1, erc2 in combinations(self.ercs, 2):
            self.computation_stats['total_pairs_checked'] += 1
            
            # Skip if ERCs cannot interact (one contains the other)
            if not can_interact(erc1, erc2, self.hierarchy):
                continue
            
            # Get fundamental synergies
            synergies = get_fundamental_synergies(erc1, erc2, self.hierarchy, self.RN)
            if synergies:
                self._store_synergies(erc1.label, erc2.label, synergies)
                self.computation_stats['synergistic_pairs'] += 1
                self.computation_stats['total_synergies'] += len(synergies)
            
            # Get complementarities
            complementarities = get_complementarity(erc1, erc2, self.hierarchy, self.RN)
            if complementarities:
                self._store_complementarities(erc1.label, erc2.label, complementarities)
                self.computation_stats['complementary_pairs'] += 1
                self.computation_stats['total_complementarities'] += len(complementarities)
            
            # Update productive partners index
            if synergies or complementarities:
                self._productive_partners[erc1.label].add(erc2.label)
                self._productive_partners[erc2.label].add(erc1.label)
        
        print(f"ERC_SORN built: {self.computation_stats['synergistic_pairs']} synergistic pairs, "
              f"{self.computation_stats['complementary_pairs']} complementary pairs")
        print(f"Total relationships: {self.computation_stats['total_synergies']} synergies, "
              f"{self.computation_stats['total_complementarities']} complementarities")
    
    def _store_synergies(self, erc1_label, erc2_label, synergies):
        """Store synergies with bidirectional lookup."""
        # Canonical ordering for consistent lookup
        key1 = (erc1_label, erc2_label)
        key2 = (erc2_label, erc1_label)
        
        self._synergies[key1] = synergies
        self._synergies[key2] = synergies  # Bidirectional
        
        # Update indexes
        self._erc_to_synergies[erc1_label].append((erc2_label, synergies))
        self._erc_to_synergies[erc2_label].append((erc1_label, synergies))
    
    def _store_complementarities(self, erc1_label, erc2_label, complementarities):
        """Store complementarities with bidirectional lookup."""
        # Canonical ordering for consistent lookup
        key1 = (erc1_label, erc2_label)
        key2 = (erc2_label, erc1_label)
        
        self._complementarities[key1] = complementarities
        self._complementarities[key2] = complementarities  # Bidirectional
        
        # Update indexes
        self._erc_to_complementarities[erc1_label].append((erc2_label, complementarities))
        self._erc_to_complementarities[erc2_label].append((erc1_label, complementarities))
    
    # ============================================================================
    # PUBLIC QUERY METHODS (for use by Persistent_Modules.py)
    # ============================================================================
    
    def get_synergies(self, erc1_label, erc2_label):
        """
        Get all fundamental synergies between two ERCs.
        
        Parameters
        ----------
        erc1_label, erc2_label : str
            Labels of the ERCs to check
            
        Returns
        -------
        list of ERC_Synergy
            List of fundamental synergies (empty if none exist)
        """
        return self._synergies.get((erc1_label, erc2_label), [])
    
    def get_complementarities(self, erc1_label, erc2_label):
        """
        Get all complementarities between two ERCs.
        
        Parameters
        ----------
        erc1_label, erc2_label : str
            Labels of the ERCs to check
            
        Returns
        -------
        list of ERC_Complementarity
            List of complementarities (empty if none exist)
        """
        return self._complementarities.get((erc1_label, erc2_label), [])
    
    def has_synergy(self, erc1_label, erc2_label):
        """Check if two ERCs have any fundamental synergy."""
        return len(self.get_synergies(erc1_label, erc2_label)) > 0
    
    def has_complementarity(self, erc1_label, erc2_label):
        """Check if two ERCs have any complementarity."""
        return len(self.get_complementarities(erc1_label, erc2_label)) > 0
    
    def has_productive_relationship(self, erc1_label, erc2_label):
        """Check if two ERCs have any productive relationship (synergy or complementarity)."""
        return self.has_synergy(erc1_label, erc2_label) or self.has_complementarity(erc1_label, erc2_label)
    
    def get_productive_partners(self, erc_label):
        """
        Get all ERCs that have productive relationships with the given ERC.
        
        Parameters
        ----------
        erc_label : str
            Label of the ERC to query
            
        Returns
        -------
        set of str
            Set of ERC labels that have productive relationships with the given ERC
        """
        return self._productive_partners.get(erc_label, set())
    
    def get_all_synergistic_partners(self, erc_label):
        """
        Get all ERCs that have synergistic relationships with the given ERC.
        
        Parameters
        ----------
        erc_label : str
            Label of the ERC to query
            
        Returns
        -------
        list of tuple
            List of (partner_label, synergies) tuples
        """
        return self._erc_to_synergies.get(erc_label, [])
    
    def get_all_complementary_partners(self, erc_label):
        """
        Get all ERCs that have complementary relationships with the given ERC.
        
        Parameters
        ----------
        erc_label : str
            Label of the ERC to query
            
        Returns
        -------
        list of tuple
            List of (partner_label, complementarities) tuples
        """
        return self._erc_to_complementarities.get(erc_label, [])
    
    def get_productive_extensions_for_erc(self, erc_label, exclude_labels=None):
        """
        Get all possible productive extensions for a given ERC.
        
        Optimized method specifically for use in find_productive_extensions.
        
        Parameters
        ----------
        erc_label : str
            Label of the ERC to find extensions for
        exclude_labels : set of str, optional
            Set of ERC labels to exclude from results
            
        Returns
        -------
        list of tuple
            List of (partner_erc_label, step_type, step_details) tuples
        """
        if exclude_labels is None:
            exclude_labels = set()
        
        extensions = []
        
        # Add synergistic extensions
        for partner_label, synergies in self._erc_to_synergies.get(erc_label, []):
            if partner_label not in exclude_labels:
                for synergy in synergies:
                    extensions.append((partner_label, 'synergy', {
                        'synergy_type': 'fundamental',
                        'synergy_object': synergy,
                        'with_erc': erc_label
                    }))
        
        # Add complementary extensions
        for partner_label, complementarities in self._erc_to_complementarities.get(erc_label, []):
            if partner_label not in exclude_labels:
                for comp in complementarities:
                    extensions.append((partner_label, 'complementarity', {
                        'comp_type': comp.comp_type,
                        'comp_object': comp,
                        'with_erc': erc_label
                    }))
        
        return extensions
    
    def get_productive_extensions_for_generator(self, erc_labels_in_generator):
        """
        Get all possible productive extensions for a generator (set of ERCs).
        
        This is the optimized version of find_productive_extensions for use
        in Persistent_Modules.py.
        
        Parameters
        ----------
        erc_labels_in_generator : list of str
            Labels of ERCs currently in the generator
            
        Returns
        -------
        list of tuple
            List of (candidate_erc_label, step_type, step_details) tuples
        """
        current_erc_set = set(erc_labels_in_generator)
        extensions = []
        candidate_extensions = {}  # candidate_label -> list of extension details
        
        # For each ERC in the current generator, find its productive partners
        for erc_label in erc_labels_in_generator:
            erc_extensions = self.get_productive_extensions_for_erc(erc_label, current_erc_set)
            
            for candidate_label, step_type, step_details in erc_extensions:
                if candidate_label not in candidate_extensions:
                    candidate_extensions[candidate_label] = []
                candidate_extensions[candidate_label].append((step_type, step_details))
        
        # Convert to the expected format
        for candidate_label, extension_list in candidate_extensions.items():
            # For now, just take the first extension found
            # (could be enhanced to handle multiple relationships)
            step_type, step_details = extension_list[0]
            extensions.append((candidate_label, step_type, step_details))
        
        return extensions
    
    # ============================================================================
    # ANALYSIS AND UTILITY METHODS
    # ============================================================================
    
    def get_statistics(self):
        """Get statistics about the SORN."""
        return self.computation_stats.copy()
    
    def get_erc_productivity_ranking(self):
        """
        Get ERCs ranked by their productivity (number of productive relationships).
        
        Returns
        -------
        list of tuple
            List of (erc_label, productivity_score) sorted by productivity
        """
        productivity_scores = []
        
        for erc in self.ercs:
            synergy_count = len(self._erc_to_synergies.get(erc.label, []))
            comp_count = len(self._erc_to_complementarities.get(erc.label, []))
            total_score = synergy_count + comp_count
            productivity_scores.append((erc.label, total_score))
        
        return sorted(productivity_scores, key=lambda x: x[1], reverse=True)
    
    def analyze_relationship_patterns(self):
        """
        Analyze patterns in synergy and complementarity relationships.
        
        Returns
        -------
        dict
            Analysis of relationship patterns
        """
        analysis = {
            'synergy_only_pairs': 0,
            'complementarity_only_pairs': 0,
            'both_relationships_pairs': 0,
            'complementarity_type_distribution': {'type1': 0, 'type2': 0, 'type3': 0},
            'avg_synergies_per_pair': 0,
            'avg_complementarities_per_pair': 0
        }
        
        checked_pairs = set()
        total_synergies = 0
        total_complementarities = 0
        synergistic_pairs = 0
        complementary_pairs = 0
        
        # Analyze each unique pair
        for (erc1_label, erc2_label), synergies in self._synergies.items():
            if (erc2_label, erc1_label) in checked_pairs:
                continue
            checked_pairs.add((erc1_label, erc2_label))
            
            has_synergy = len(synergies) > 0
            complementarities = self._complementarities.get((erc1_label, erc2_label), [])
            has_complementarity = len(complementarities) > 0
            
            if has_synergy and has_complementarity:
                analysis['both_relationships_pairs'] += 1
            elif has_synergy:
                analysis['synergy_only_pairs'] += 1
            elif has_complementarity:
                analysis['complementarity_only_pairs'] += 1
            
            if has_synergy:
                synergistic_pairs += 1
                total_synergies += len(synergies)
            
            if has_complementarity:
                complementary_pairs += 1
                total_complementarities += len(complementarities)
                
                # Count complementarity types
                for comp in complementarities:
                    comp_type_key = f'type{comp.comp_type}'
                    analysis['complementarity_type_distribution'][comp_type_key] += 1
        
        # Calculate averages
        if synergistic_pairs > 0:
            analysis['avg_synergies_per_pair'] = total_synergies / synergistic_pairs
        if complementary_pairs > 0:
            analysis['avg_complementarities_per_pair'] = total_complementarities / complementary_pairs
        
        return analysis
    
    def __repr__(self):
        """String representation of the ERC_SORN."""
        return (f"ERC_SORN({len(self.ercs)} ERCs, "
                f"{self.computation_stats['synergistic_pairs']} synergistic pairs, "
                f"{self.computation_stats['complementary_pairs']} complementary pairs)")

# ============================================================================
# FACTORY FUNCTIONS AND UTILITIES
# ============================================================================

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

def get_all_productive_pairs(erc_sorn):
    """
    Get all pairs of ERCs that have productive relationships.
    
    Parameters
    ----------
    erc_sorn : ERC_SORN
        The second-order reaction network
        
    Returns
    -------
    list of tuple
        List of (erc1_label, erc2_label, relationship_types) tuples
    """
    productive_pairs = []
    checked_pairs = set()
    
    for erc in erc_sorn.ercs:
        for partner_label in erc_sorn.get_productive_partners(erc.label):
            pair = tuple(sorted([erc.label, partner_label]))
            if pair in checked_pairs:
                continue
            checked_pairs.add(pair)
            
            relationship_types = []
            if erc_sorn.has_synergy(erc.label, partner_label):
                relationship_types.append('synergy')
            if erc_sorn.has_complementarity(erc.label, partner_label):
                relationship_types.append('complementarity')
            
            productive_pairs.append((pair[0], pair[1], relationship_types))
    
    return productive_pairs

def analyze_sorn_efficiency(erc_sorn, verbose=True):
    """
    Analyze the efficiency gains provided by the ERC_SORN.
    
    Parameters
    ----------
    erc_sorn : ERC_SORN
        The second-order reaction network
    verbose : bool
        Whether to print analysis results
        
    Returns
    -------
    dict
        Efficiency analysis results
    """
    n_ercs = len(erc_sorn.ercs)
    total_possible_pairs = n_ercs * (n_ercs - 1) // 2
    
    analysis = {
        'total_ercs': n_ercs,
        'total_possible_pairs': total_possible_pairs,
        'pairs_with_relationships': erc_sorn.computation_stats['synergistic_pairs'] + erc_sorn.computation_stats['complementary_pairs'],
        'relationship_density': 0,
        'synergy_density': 0,
        'complementarity_density': 0,
        'efficiency_factor': 0  # How many repeated calculations we avoid
    }
    
    if total_possible_pairs > 0:
        analysis['relationship_density'] = analysis['pairs_with_relationships'] / total_possible_pairs
        analysis['synergy_density'] = erc_sorn.computation_stats['synergistic_pairs'] / total_possible_pairs
        analysis['complementarity_density'] = erc_sorn.computation_stats['complementary_pairs'] / total_possible_pairs
    
    # Estimate efficiency factor (how many calculations we avoid during generation)
    # Rough estimate: each generator extension might check O(n) pairs, and we might build O(n^k) generators
    analysis['efficiency_factor'] = total_possible_pairs  # Conservative estimate
    
    if verbose:
        print("ERC_SORN Efficiency Analysis:")
        print(f"  Total ERCs: {analysis['total_ercs']}")
        print(f"  Possible pairs: {analysis['total_possible_pairs']}")
        print(f"  Productive pairs: {analysis['pairs_with_relationships']} ({analysis['relationship_density']:.2%})")
        print(f"  Synergistic pairs: {erc_sorn.computation_stats['synergistic_pairs']} ({analysis['synergy_density']:.2%})")
        print(f"  Complementary pairs: {erc_sorn.computation_stats['complementary_pairs']} ({analysis['complementarity_density']:.2%})")
        print(f"  Estimated efficiency factor: ~{analysis['efficiency_factor']}x fewer calculations")
    
    return analysis