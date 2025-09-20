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
    # Check if one ERC contains the other using cached containment relationships
    if erc1 in hierarchy.get_contain(erc2) or erc2 in hierarchy.get_contain(erc1):
        return False
    return True

def has_partial_overlap_with_generators(base_erc, target_erc, RN):
    """
    OPTIMIZED: Check if base can contribute to covering target's generators.
    
    This is the necessary condition for a base to participate in synergies
    with the target. A base can contribute if it has partial overlap with
    at least one generator (intersection but not full coverage).
    """
    base_closure = base_erc.get_closure_names(RN)
    
    # Check each minimal generator of the target
    for generator in target_erc.min_generators:
        gen_species = set(species_list_to_names(generator))
        
        # Check for partial overlap with this generator
        intersection = base_closure & gen_species
        
        # Partial overlap means: has intersection but doesn't fully cover
        if intersection and not gen_species.issubset(base_closure):
            return True
    
    return False
def has_generator_coverage(base1_erc, base2_erc, target_erc, RN):
    """
    OPTIMIZED: Check if combined bases cover at least one complete generator.
    
    This is the NECESSARY AND SUFFICIENT condition for fundamental synergy.
    The combined closures must cover at least one generator completely,
    but neither base alone can cover that generator.
    """
    base1_closure = base1_erc.get_closure_names(RN)
    base2_closure = base2_erc.get_closure_names(RN)
    combined_closure = base1_closure | base2_closure
    
    # Check if combined bases cover at least one generator completely
    for generator in target_erc.min_generators:
        gen_species = set(species_list_to_names(generator))
        
        # Check if combined closure covers this generator
        if gen_species.issubset(combined_closure):
            # Ensure it's actually a synergy (neither base alone covers it)
            if (not gen_species.issubset(base1_closure) and 
                not gen_species.issubset(base2_closure)):
                return True
    
    return False
# ============================================================================
# SYNERGY DETECTION: BRUTE FORCE FASHION ONLY FOR TESTING
# ============================================================================

def get_basic_synergies(erc1, erc2, hierarchy, RN):
    """
    Generate all synergetic pairs by brute force
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
        
        # CORRECTED: Track all target )s separately
        target_ercs_found = set()  # Track which targets we've found
        
        for r in joint_closure_reacs:
            r_name = r.name() if hasattr(r, 'name') else str(r)
            if r_name in novel_reac_names:
                syn_erc = hierarchy.get_erc_from_reaction(RN, hierarchy, r)
                if syn_erc is None:
                    continue
                
                # CORRECTED: Only check for exact duplicates, not containment
                if syn_erc.label not in target_ercs_found:
                    target_ercs_found.add(syn_erc.label)
                    synergies.append(ERC_Synergy([erc1, erc2], syn_erc, "regular"))
    
    return synergies

def get_maximal_synergies(erc1, erc2, hierarchy, RN):
    """
    A synergy ERC1 + ERC2 → Target is maximal if there's no other synergy
    ERC1 + ERC2 → LargerTarget where LargerTarget contains Target.
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

def get_fundamental_synergies_brute_force(erc1, erc2, hierarchy, RN, maximal_synergies=None, verbose=False):
    """    
    A maximal synergy ERC1 + ERC2 → Target is fundamental if there's no synergy 
    MoreFundamental1 + MoreFundamental2 → Target where the reactants are more fundamental.  
    CRITICAL: Fundamentality is checked SEPARATELY for each target!
    """
    if maximal_synergies == None:
        maximal_synergies = get_maximal_synergies(erc1, erc2, hierarchy, RN)
    if not maximal_synergies:
        return []
    
    fundamental_synergies = []
    
    # Get successors (more fundamental ERCs) for each reactant
    if hierarchy.graph:
        erc1_successors = set(nx.descendants(hierarchy.graph, erc1.label))
        erc1_successors.add(erc1.label)  # Include self
        
        erc2_successors = set(nx.descendants(hierarchy.graph, erc2.label))  
        erc2_successors.add(erc2.label)  # Include self
    else:
        # No hierarchy - only self is considered
        erc1_successors = {erc1.label}
        erc2_successors = {erc2.label}
    
    for syn in maximal_synergies:
        is_fundamental = True
        target_label = syn.plabel
        
        # CORRECTED: Only check more fundamental pairs
        for other_erc1_label in erc1_successors:
            for other_erc2_label in erc2_successors:
                # Skip if it's the same ERC (already handled)
                if other_erc1_label == other_erc2_label:
                    continue
                    
                # Skip if it's the exact same pair (no improvement in fundamentality)
                if (set([other_erc1_label, other_erc2_label]) == 
                    set([erc1.label, erc2.label])):
                    continue
                
                # Get the ERC objects
                other_erc1 = hierarchy.get_erc_by_label(other_erc1_label)
                other_erc2 = hierarchy.get_erc_by_label(other_erc2_label)
                
                # Get maximal synergies for this more fundamental pair
                other_maximal = get_maximal_synergies(other_erc1, other_erc2, hierarchy, RN)
                
                for other_syn in other_maximal:
                    # CRITICAL: Only compare synergies that produce the SAME target
                    if other_syn.plabel == target_label:
                        # Found a more fundamental pair that produces same target
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
    BRUTE FORCE WRAPPER: Find all fundamental synergies by exhaustive pairwise checking.
    """
    all_fundamental_synergies = []
    
    total_pairs = len(list(combinations(ercs, 2)))
    
    for i, (erc1, erc2) in enumerate(combinations(ercs, 2)):
        if verbose:
            print(f"Checking pair {i+1}/{total_pairs}: {erc1.label} + {erc2.label}")
        
        # Call the PAIRWISE brute force function (note: plural "synergies")
        pair_synergies = get_fundamental_synergies_brute_force(erc1, erc2, hierarchy, RN)
        all_fundamental_synergies.extend(pair_synergies)
        
        if verbose and pair_synergies:
            print(f"  Found {len(pair_synergies)} fundamental synergies")
    
    return all_fundamental_synergies


# ============================================================================
# SYNERGY COMPUTATION EFFICIENT ALGORITHM...still not working!
# ============================================================================
def update_pruning_structures(synergy, global_pruned_pairs, level_combinations_to_skip, hierarchy, verbose=False):
    """
    Update pruning structures when a fundamental synergy is discovered.
    
    Prunes non-fundamental combinations based on hierarchy relationships:
    1. Same pair for descendant targets (less maximal)
    2. Pairs using less fundamental bases (ancestors) for same/descendant targets
    3. Bulk level combinations when patterns are clear
    """
    base1_label = synergy.rlabel[0] 
    base2_label = synergy.rlabel[1]
    target_label = synergy.plabel
    
    if not hierarchy.graph:
        return 0  # No pruning possible without hierarchy
    
    initial_pruned_count = len(global_pruned_pairs)
    
    # Get hierarchy information
    levels = ERC.get_node_levels(hierarchy.graph)
    
    # Get ancestors (containers) of bases - these are LESS fundamental
    base1_ancestors = set(nx.ancestors(hierarchy.graph, base1_label)) if base1_label in hierarchy.graph else set()
    base2_ancestors = set(nx.ancestors(hierarchy.graph, base2_label)) if base2_label in hierarchy.graph else set()
    
    # Get descendants (contained by) of target - these are LESS maximal  
    target_descendants = set(nx.descendants(hierarchy.graph, target_label)) if target_label in hierarchy.graph else set()
    
    current_pair = tuple(sorted([base1_label, base2_label]))
    target_level = levels.get(target_label, 0)
    base1_level = levels.get(base1_label, 0)
    base2_level = levels.get(base2_label, 0)
    
    if verbose:
        print(f"        Pruning for synergy: {base1_label}(L{base1_level}) + {base2_label}(L{base2_level}) → {target_label}(L{target_level})")
    
    # 1. INDIVIDUAL PAIR PRUNING
    
    # Mark the discovered pair as processed
    global_pruned_pairs.add(current_pair)
    
    # Prune less fundamental base combinations
    # Pattern: If base1_ancestor → base1, then base1_ancestor + base2 is less fundamental than base1 + base2
    
    for ancestor1 in base1_ancestors:
        # ancestor1 + base2 (less fundamental than base1 + base2)
        less_fundamental_pair = tuple(sorted([ancestor1, base2_label]))
        global_pruned_pairs.add(less_fundamental_pair)
        
        # ancestor1 + ancestor2 (less fundamental than base1 + base2)
        for ancestor2 in base2_ancestors:
            less_fundamental_pair = tuple(sorted([ancestor1, ancestor2]))
            global_pruned_pairs.add(less_fundamental_pair)
    
    # Symmetric case: base1 + ancestor2
    for ancestor2 in base2_ancestors:
        less_fundamental_pair = tuple(sorted([base1_label, ancestor2]))
        global_pruned_pairs.add(less_fundamental_pair)
    
    # 2. BULK LEVEL COMBINATION PRUNING
    
    # Conservative bulk pruning: Skip level combinations that are guaranteed non-fundamental
    # Only skip if ALL pairs at those levels would be less fundamental
    
    # If we found a fundamental synergy at (base1_level, base2_level),
    # then combinations (ancestor_level, base2_level) where ancestor_level > base1_level are less fundamental
    # BUT only if there are no other fundamental ERCs at ancestor_level that aren't ancestors of base1
    
    # For now, implement conservative individual-pair pruning only
    # Bulk pruning can be added later as an optimization if needed
    
    # 3. CROSS-TARGET PRUNING 
    
    # For descendant targets (less maximal), the same base combinations are less fundamental
    for desc_target in target_descendants:
        desc_level = levels.get(desc_target, 0)
        if desc_level < target_level:  # Only prune truly less maximal targets
            # The exact same pair for descendant target is less fundamental
            global_pruned_pairs.add(current_pair)  # Already added, but explicit
            
            # Pairs with ancestors for descendant targets are also less fundamental
            for ancestor1 in base1_ancestors:
                less_fundamental_pair = tuple(sorted([ancestor1, base2_label]))
                global_pruned_pairs.add(less_fundamental_pair)
                
                for ancestor2 in base2_ancestors:
                    less_fundamental_pair = tuple(sorted([ancestor1, ancestor2]))
                    global_pruned_pairs.add(less_fundamental_pair)
            
            for ancestor2 in base2_ancestors:
                less_fundamental_pair = tuple(sorted([base1_label, ancestor2]))
                global_pruned_pairs.add(less_fundamental_pair)
    
    # Calculate how many pairs were pruned by this update
    final_pruned_count = len(global_pruned_pairs)
    pairs_added = final_pruned_count - initial_pruned_count
    
    if verbose and pairs_added > 0:
        print(f"          → Pruned {pairs_added} non-fundamental pair combinations")
    
    return pairs_added
def get_fundamental_synergies(ercs, hierarchy, RN, verbose=True):
    """
    HIERARCHICAL CLOSURE-BASED: Find fundamental synergies using systematic 
    hierarchical exploration with actual closure computation.
    
    Guarantees fundamentality by construction through bottom-up base exploration
    and top-down target discovery with dynamic pruning.
    """
    if verbose:
        print("="*70)
        print("HIERARCHICAL CLOSURE-BASED SYNERGY FINDER")
        print("="*70)   
        print(f"Network: {len(ercs)} ERCs, {len(RN.species())} species")
        if hierarchy.graph:
            print(f"Hierarchy: {hierarchy.graph.number_of_edges()} containment relations")
    
    import time
    start_time = time.time()
    
    # Group ERCs by hierarchy level
    ercs_by_level = group_ercs_by_level(ercs, hierarchy)
    max_level = max(ercs_by_level.keys()) if ercs_by_level else 0
    
    # Pruning structures
    global_pruned_pairs = set()
    level_combinations_to_skip = set()
    discovered_synergies = []
    
    # Statistics
    stats = {
        'pairs_checked': 0,
        'pairs_skipped_pruning': 0,
        'pairs_skipped_containment': 0,
        'synergies_found': 0,
        'level_combinations_processed': 0,
        'level_combinations_skipped': 0
    }
    
    if verbose:
        print(f"Hierarchy depth: {max_level} levels")
        for level in range(max_level + 1):
            level_ercs = ercs_by_level.get(level, [])
            print(f"  Level {level}: {len(level_ercs)} ERCs")
    
    # HIERARCHICAL EXPLORATION: Level combinations from most fundamental to least
    for level_i in range(max_level + 1):
        for level_j in range(level_i, max_level + 1):
            
            # Check if this entire level combination can be skipped
            if (level_i, level_j) in level_combinations_to_skip:
                stats['level_combinations_skipped'] += 1
                continue
            
            stats['level_combinations_processed'] += 1
            
            if verbose:
                print(f"\nProcessing level combination ({level_i}, {level_j})")
            
            # Generate pairs for these specific levels
            ercs_level_i = ercs_by_level.get(level_i, [])
            ercs_level_j = ercs_by_level.get(level_j, [])
            
            level_pairs = []
            for erc1 in ercs_level_i:
                for erc2 in ercs_level_j:
                    # Avoid duplicates within same level
                    if level_i == level_j and erc1.label >= erc2.label:
                        continue
                    
                    # Check containment (no synergy possible)
                    if hierarchy.graph and (
                        hierarchy.graph.has_edge(erc1.label, erc2.label) or 
                        hierarchy.graph.has_edge(erc2.label, erc1.label)):
                        stats['pairs_skipped_containment'] += 1
                        continue
                    
                    # Check global pruning
                    pair_key = tuple(sorted([erc1.label, erc2.label]))
                    if pair_key in global_pruned_pairs:
                        stats['pairs_skipped_pruning'] += 1
                        continue
                    
                    level_pairs.append((erc1, erc2))
            
            if verbose:
                print(f"  Valid pairs to check: {len(level_pairs)}")
            
            # Process pairs at this level combination
            for base1, base2 in level_pairs:
                stats['pairs_checked'] += 1
                
                # CLOSURE-BASED SYNERGY DETECTION
                synergies = detect_synergies_by_closure(base1, base2, hierarchy, RN)
                
                for synergy in synergies:
                    discovered_synergies.append(synergy)
                    stats['synergies_found'] += 1
                    
                    if verbose:
                        print(f"    FUNDAMENTAL SYNERGY: {synergy.rlabel[0]} + {synergy.rlabel[1]} → {synergy.plabel}")
                    
                    # DYNAMIC PRUNING
                    pruned_count = update_pruning_structures(
                        synergy, global_pruned_pairs, level_combinations_to_skip, hierarchy, verbose)
    
    end_time = time.time()
    
    if verbose:
        print(f"\n" + "="*70)
        print("HIERARCHICAL CLOSURE-BASED RESULTS")
        print("="*70)
        print(f"Level combinations processed: {stats['level_combinations_processed']}")
        print(f"Level combinations skipped: {stats['level_combinations_skipped']}")
        print(f"Pairs checked: {stats['pairs_checked']}")
        print(f"Pairs skipped by pruning: {stats['pairs_skipped_pruning']}")
        print(f"Pairs skipped by containment: {stats['pairs_skipped_containment']}")
        print(f"Fundamental synergies found: {stats['synergies_found']}")
        print(f"Time elapsed: {end_time - start_time:.3f} seconds")
        
        if discovered_synergies:
            print(f"\nFundamental Synergies:")
            for i, synergy in enumerate(discovered_synergies, 1):
                print(f"  {i}. {synergy.rlabel[0]} + {synergy.rlabel[1]} → {synergy.plabel}")
    
    return discovered_synergies


def group_ercs_by_level(ercs, hierarchy):
    """Group ERCs by their hierarchy level."""
    if hierarchy.graph:
        levels = ERC.get_node_levels(hierarchy.graph)
    else:
        levels = {erc.label: 0 for erc in ercs}
    
    ercs_by_level = defaultdict(list)
    for erc in ercs:
        level = levels.get(erc.label, 0)
        ercs_by_level[level].append(erc)
    
    return ercs_by_level


def detect_synergies_by_closure(base1, base2, hierarchy, RN):
    """
    Detect synergies using actual closure computation (like brute force approach).
    
    Returns list of fundamental synergies discovered from this base pair.
    """
    synergies = []
    
    # Compute individual closures
    base1_closure = base1.get_closure(RN)
    base2_closure = base2.get_closure(RN)
    
    # Compute joint closure
    union_species = list(set(base1_closure).union(set(base2_closure)))
    joint_closure = closure(RN, union_species)
    
    # Get reactions from closures
    base1_reactions = RN.get_reactions_from_species(base1_closure)
    base2_reactions = RN.get_reactions_from_species(base2_closure)
    joint_reactions = RN.get_reactions_from_species(joint_closure)
    
    # Find novel reactions
    base1_reaction_names = {r.name() for r in base1_reactions}
    base2_reaction_names = {r.name() for r in base2_reactions}
    union_reaction_names = base1_reaction_names.union(base2_reaction_names)
    
    joint_reaction_names = {r.name() for r in joint_reactions}
    
    if len(joint_reaction_names) > len(union_reaction_names):
        novel_reaction_names = joint_reaction_names - union_reaction_names
        
        # Map novel reactions to target ERCs
        target_ercs_found = set()
        
        for reaction in joint_reactions:
            reaction_name = reaction.name()
            if reaction_name in novel_reaction_names:
                # Find which ERC this novel reaction belongs to
                target_erc = hierarchy.get_erc_from_reaction(RN, hierarchy, reaction)
                if target_erc is None:
                    continue
                
                # Avoid duplicate targets
                if target_erc.label not in target_ercs_found:
                    target_ercs_found.add(target_erc.label)
                    synergies.append(ERC_Synergy([base1, base2], target_erc, "fundamental"))
    
    return synergies


# ============================================================================
# COMPLEMENTARITY DETECTION FUNCTIONS
# ============================================================================

def is_complementary_type1(erc1, erc2, req1, req2, joint_req, joint_consumed,
                           joint_produced, prod1, prod2, hierarchy, RN):
    """
    Check Type 1 complementarity: |reqs(X')|<|reqs(X)|+|reqs(E)|
    Parameters
    ----------
    erc1, erc2 : ERC
        ERCs to check for Type 1 complementarity
    hierarchy : ERC_Hierarchy
        Hierarchy containing the ERCs
    RN : ReactionNetwork
        The reaction network
    req1, req2, prod1, prod2 : set
        Sets of requirements and produced species for ERC1 and ERC2
    joint_req, joint_consumed, joint_produced : set
        Set of joint requirements, consumed, and produced species for both ERCs
    Returns
    -------
    tuple
        (is_complementary: bool, info: dict)
    """

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
    else:
        return False, {}

def is_complementary_type2(erc1, erc2, req1, req2, joint_req, joint_consumed,
                           joint_produced, prod1, prod2, hierarchy, RN):
    """
    Check Type 2 complementarity: |reqs(X')|=|reqs(X)|+|reqs(E)| and reqs(X')≠reqs(X)
    
    CORRECTED: Uses closure operator (∨) to calculate joint requirements.
    No reduction in total requirements but different requirement set.
    
    Parameters
    ----------
    erc1, erc2 : ERC
        ERCs to check for Type 1 complementarity
    hierarchy : ERC_Hierarchy
        Hierarchy containing the ERCs
    RN : ReactionNetwork
        The reaction network
    req1, req2, prod1, prod2 : set
        Sets of requirements and produced species for ERC1 and ERC2
    joint_req, joint_consumed, joint_produced : set
        Set of joint requirements, consumed, and produced species for both ERCs
    Returns
    -------
    tuple
        (is_complementary: bool, info: dict)
    """
    if not can_interact(erc1, erc2, hierarchy):
        return False, {}
        
    try:
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

def is_complementary_type3(erc1, erc2, req1, req2, joint_req, joint_consumed,
                           joint_produced, prod1, prod2, hierarchy, RN):
    """
    Check Type 3 complementarity: reqs(X')=reqs(X) and prods(X')≠prods(X)
    
    CORRECTED: Uses closure operator (∨) to calculate joint requirements and products.
    Same requirements but different products (implies synergy).
    
    Parameters
    ----------
    erc1, erc2 : ERC
        ERCs to check for Type 1 complementarity
    hierarchy : ERC_Hierarchy
        Hierarchy containing the ERCs
    RN : ReactionNetwork
        The reaction network
    req1, req2, prod1, prod2 : set
        Sets of requirements and produced species for ERC1 and ERC2
    joint_req, joint_consumed, joint_produced : set
        Set of joint requirements, consumed, and produced species for both ERCs
    Returns
    -------
    tuple
        (is_complementary: bool, info: dict)
    """
    if not can_interact(erc1, erc2, hierarchy):
        return False, {}
        
    try:
        # Get individual ERC requirements and products (cached)

        
        # Check if joint_prod contains elements not in the prod1 or prod2
        novel_products = set(joint_produced) - (prod1 | prod2)

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
                'novel_products': joint_produced - (prod1 | prod2),
                'has_synergy': len(joint_produced - (prod1 | prod2)) >0
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

    # Check Type 1
    is_type1, info1 = is_complementary_type1(erc1, erc2, req1, req2, joint_req, joint_consumed,
                           joint_produced, prod1, prod2, hierarchy, RN)
    if is_type1:
        complementarities.append(ERC_Complementarity(erc1, erc2, 1, info1))
    
    # Check Type 2
    is_type2, info2 = is_complementary_type2(erc1, erc2, req1, req2, joint_req, joint_consumed,
                           joint_produced, prod1, prod2, hierarchy, RN)
    if is_type2:
        complementarities.append(ERC_Complementarity(erc1, erc2, 2, info2))
    
    # Check Type 3
    is_type3, info3 = is_complementary_type3(erc1, erc2, req1, req2, joint_req, joint_consumed,
                           joint_produced, prod1, prod2, hierarchy, RN)
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
        self._synergies = {}  # (a,b) -> list[ERC_Synergy]
        self._complementarities = {}  # (a,b) -> list[ERC_Complementarity]

        # Indexes
        self._erc_to_synergies = defaultdict(list)
        self._erc_to_complementarities = defaultdict(list)
        self._productive_partners = defaultdict(set)

        # Minimal stats
        self.computation_stats = {
            'total_pairs_checked': 0,
            'productive_pairs': 0,
            'total_synergies': 0,
            'total_complementarities': 0,
            'synergistic_pairs': 0,
            'build_time': 0.0
        }
        
        # Build the SORN
        self._build_sorn()
    
    def _build_sorn(self):
        start = time.time()
        print(f"Building ERC_SORN for {len(self.ercs)} ERCs...")

        # Fix: Pass pairs of ERCs and RN correctly
        all_synergies = []
        for erc1 in self.ercs:
            for erc2 in self.ercs:
                if erc1.label >= erc2.label:
                    continue
                # Pass arguments in correct order: erc1, erc2, hierarchy, RN
                synergies = get_fundamental_synergies_brute_force(erc1, erc2, self.hierarchy, self.RN)
                all_synergies.extend(synergies)

        # Map synergies by pair key for O(1) per-pair lookup
        synergies_by_pair = defaultdict(list)
        for syn in all_synergies:
            a_label = syn.rlabel[0]
            b_label = syn.rlabel[1]
            key = tuple(sorted((a_label, b_label)))
            synergies_by_pair[key].append(syn)

        # Update internal stores with synergies_by_pair
        for pair_key, syn_list in synergies_by_pair.items():
            a_label, b_label = pair_key
            self._store_synergies(a_label, b_label, syn_list)
            # productivity indexes updated in _store_synergies

        # STEP B: iterate all pairs to compute complementarities (one pass)
        for erc1, erc2 in combinations(self.ercs, 2):
            self.computation_stats['total_pairs_checked'] += 1

            if not can_interact(erc1, erc2, self.hierarchy):
                continue

            key = tuple(sorted((erc1.label, erc2.label)))
            pair_has_productive = False

            # Complementarities (compute once per pair)
            comps = get_complementarity(erc1, erc2, self.hierarchy, self.RN)
            if comps:
                self._store_complementarities(erc1.label, erc2.label, comps)
                pair_has_productive = True
                self.computation_stats['total_complementarities'] += len(comps)

            # Check if this pair had synergies (from the precomputed map)
            syns = synergies_by_pair.get(key, [])
            if syns:
                # already stored by _store_synergies above; but mark pair as productive
                pair_has_productive = True
                self.computation_stats['total_synergies'] += len(syns)
                self.computation_stats['synergistic_pairs'] += 1
            if pair_has_productive:
                self._productive_partners[erc1.label].add(erc2.label)
                self._productive_partners[erc2.label].add(erc1.label)
                self.computation_stats['productive_pairs'] += 1

        # Finalize stats
        self.computation_stats['build_time'] = time.time() - start

        print(f"ERC_SORN built: productive_pairs={self.computation_stats['productive_pairs']}, "
            f"synergies={self.computation_stats['total_synergies']}, "
            f"complementarities={self.computation_stats['total_complementarities']}, "
            f"time={self.computation_stats['build_time']:.2f}s")
    
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
                #synergistic_pairs += 1
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