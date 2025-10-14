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
# OPTIMIZED FUNDAMENTAL SYNERGY DETECTION WITH TARGET-SPECIFIC PRUNING
# ============================================================================
# Replace the get_all_fundamental_synergies_optimized function in 
# ERC_Synergy_Complementarity.py with this improved version

def get_all_fundamental_synergies_optimized(ercs, hierarchy, RN, verbose=False):
    """
    Find all fundamental synergies using aggressive target-specific pruning.
    
    Key improvements:
    1. Fundamentality is per-triplet (E1, E2, E3), not per-pair
    2. Once E1 + E2 ↠ E3 is found fundamental:
       - Prune all (E1', E2, E3) where E1' ⊇ E1 (pair domination for same target)
       - Prune all (E1, E2', E3) where E2' ⊇ E2 (pair domination for same target)
       - Prune all (E1, E2, E4) where E3 → E4 (target containment - higher targets)
    3. Process ERCs bottom-up (smallest first) to maximize early pruning
    
    Parameters
    ----------
    ercs : list
        List of ERC objects
    hierarchy : ERC_Hierarchy
        The hierarchy object containing ERC containment relationships
    RN : ReactionNetwork
        The reaction network
    verbose : bool, optional
        If True, print detailed progress information
        
    Returns
    -------
    list
        List of ERC_Synergy objects with synergy_type="fundamental"
    """
    
    if verbose:
        print("\n" + "="*70)
        print("OPTIMIZED FUNDAMENTAL SYNERGY COMPUTATION")
        print("="*70)
        print("Computing ERC levels for bottom-up processing...")
    
    # Compute hierarchy levels for bottom-up processing
    from pyCOT.ERC_Hierarchy import ERC
    if hierarchy.graph:
        levels = ERC.get_node_levels(hierarchy.graph)
    else:
        levels = {erc.label: 0 for erc in ercs}
    
    if verbose:
        distinct_levels = len(set(levels.values()))
        print(f"  Found {distinct_levels} distinct levels in hierarchy")
        print(f"  Level 0 (leaves): {sum(1 for v in levels.values() if v == 0)} ERCs")
        print(f"  Max level: {max(levels.values())}")
    
    # Pre-compute containment relationships for fast lookup
    descendants_cache = {}  # ERC -> set of ERCs it directly contains
    ancestors_cache = {}    # ERC -> set of ERCs that contain it
    
    for erc in ercs:
        if hierarchy.graph:
            descendants_cache[erc.label] = set(nx.descendants(hierarchy.graph, erc.label))
            ancestors_cache[erc.label] = set(nx.ancestors(hierarchy.graph, erc.label))
        else:
            descendants_cache[erc.label] = set()
            ancestors_cache[erc.label] = set()
    
    # Sort pairs by combined level (LOWER first = smaller ERCs = BOTTOM-UP)
    def pair_level(pair):
        erc1, erc2 = pair
        level_sum = levels[erc1.label] + levels[erc2.label]
        level_max = max(levels[erc1.label], levels[erc2.label])
        return (level_sum, level_max)
    
    all_pairs = list(combinations(ercs, 2))
    sorted_pairs = sorted(all_pairs, key=pair_level)
    
    if verbose:
        print(f"\nTotal ERC pairs to consider: {len(sorted_pairs)}")
        print("Starting fundamental synergy search with dynamic pruning...\n")
    
    # Track pruned triplets: (erc1_label, erc2_label, target_label) -> reason
    pruned_triplets = {}
    
    # Track found fundamental synergies: target_label -> [(pair_labels, level)]
    fundamental_producers_per_target = defaultdict(list)
    
    # Track all targets a pair can produce
    pair_to_targets = defaultdict(list)
    
    # Results
    fundamental_synergies = []
    
    # Statistics
    pairs_checked = 0
    pairs_with_synergies = 0
    total_maximal = 0
    synergies_pruned_by_pair = 0
    synergies_pruned_by_target = 0
    targets_seen = set()
    
    # Main loop: process pairs bottom-up
    for erc1, erc2 in sorted_pairs:
        pairs_checked += 1
        
        current_pair_labels = tuple(sorted([erc1.label, erc2.label]))
        current_level_sum = levels[erc1.label] + levels[erc2.label]
        
        # Get maximal synergies for this pair
        maximal_synergies = get_maximal_synergies(erc1, erc2, hierarchy, RN)
        
        if not maximal_synergies:
            continue
        
        pairs_with_synergies += 1
        total_maximal += len(maximal_synergies)
        
        if verbose and pairs_checked <= 10:
            print(f"Pair {pairs_checked}: {erc1.label}+{erc2.label} "
                  f"(levels: {levels[erc1.label]},{levels[erc2.label]})")
            print(f"  Found {len(maximal_synergies)} maximal synergies")
        
        # Check each maximal synergy for fundamentality
        for syn in maximal_synergies:
            target_label = syn.plabel
            targets_seen.add(target_label)
            triplet_key = (erc1.label, erc2.label, target_label)
            
            # Check if this triplet was pruned
            if triplet_key in pruned_triplets:
                if verbose and pairs_checked <= 10:
                    print(f"    PRUNED: →{target_label} ({pruned_triplets[triplet_key]})")
                continue
            
            # Check fundamentality: is there a more fundamental pair for this target?
            is_fundamental = True
            pruned_reason = None
            
            # Check if a more fundamental pair (with smaller ERCs) already produces this target
            for (other_pair, other_level) in fundamental_producers_per_target[target_label]:
                other_erc1_label, other_erc2_label = other_pair
                
                # Check if current pair is dominated by the other pair
                # (other pair uses smaller or equal ERCs in containment hierarchy)
                erc1_dominated = (erc1.label in ancestors_cache[other_erc1_label] or 
                                 erc1.label == other_erc1_label)
                erc2_dominated = (erc2.label in ancestors_cache[other_erc2_label] or 
                                 erc2.label == other_erc2_label)
                
                # If both ERCs are dominated and at least one is strictly dominated
                if erc1_dominated and erc2_dominated:
                    if (erc1.label in ancestors_cache[other_erc1_label] or 
                        erc2.label in ancestors_cache[other_erc2_label]):
                        is_fundamental = False
                        pruned_reason = f"Dominated by ({other_erc1_label},{other_erc2_label})"
                        synergies_pruned_by_pair += 1
                        break
            
            # Check if this pair already produces a smaller target (target containment)
            if is_fundamental:
                for other_target in pair_to_targets[current_pair_labels]:
                    # If other_target → target_label (other is smaller, contains target)
                    if target_label in ancestors_cache[other_target]:
                        is_fundamental = False
                        pruned_reason = f"Pair produces smaller target {other_target}"
                        synergies_pruned_by_target += 1
                        break
            
            # If fundamental, add it and prune non-fundamental variants
            if is_fundamental:
                if verbose and pairs_checked <= 10:
                    print(f"    FUNDAMENTAL: →{target_label}")
                
                fundamental_synergies.append(
                    ERC_Synergy(syn.reactants, syn.product, "fundamental")
                )
                
                # Track this fundamental synergy
                fundamental_producers_per_target[target_label].append(
                    (current_pair_labels, current_level_sum)
                )
                pair_to_targets[current_pair_labels].append(target_label)
                
                # PRUNING STEP 1: Prune descendants for SAME target
                # Any pair with larger ERCs producing this same target is non-fundamental
                for ancestor1_label in ancestors_cache[erc1.label]:
                    pruned_key = (ancestor1_label, erc2.label, target_label)
                    if pruned_key not in pruned_triplets:
                        pruned_triplets[pruned_key] = f"Pair dominated: {erc1.label}→{ancestor1_label}"
                
                for ancestor2_label in ancestors_cache[erc2.label]:
                    pruned_key = (erc1.label, ancestor2_label, target_label)
                    if pruned_key not in pruned_triplets:
                        pruned_triplets[pruned_key] = f"Pair dominated: {erc2.label}→{ancestor2_label}"
                
                # Cross products: both ERCs larger
                for ancestor1_label in ancestors_cache[erc1.label]:
                    for ancestor2_label in ancestors_cache[erc2.label]:
                        pruned_key = (ancestor1_label, ancestor2_label, target_label)
                        if pruned_key not in pruned_triplets:
                            pruned_triplets[pruned_key] = "Both dominated"
                
                # PRUNING STEP 2: Prune SAME pair for higher targets
                # If E1 + E2 ↠ E3, then E1 + E2 → E4 cannot be fundamental when E3 → E4
                for higher_target_label in ancestors_cache[target_label]:
                    pruned_key = (erc1.label, erc2.label, higher_target_label)
                    if pruned_key not in pruned_triplets:
                        pruned_triplets[pruned_key] = f"Lower target found: {target_label}→{higher_target_label}"
                    
                    # Also prune descendants of this pair for higher targets
                    for ancestor1_label in ancestors_cache[erc1.label]:
                        pruned_key = (ancestor1_label, erc2.label, higher_target_label)
                        if pruned_key not in pruned_triplets:
                            pruned_triplets[pruned_key] = f"Dominated pair + higher target"
                    
                    for ancestor2_label in ancestors_cache[erc2.label]:
                        pruned_key = (erc1.label, ancestor2_label, higher_target_label)
                        if pruned_key not in pruned_triplets:
                            pruned_triplets[pruned_key] = f"Dominated pair + higher target"
            
            elif verbose and pairs_checked <= 10:
                print(f"    NON-FUNDAMENTAL: →{target_label} ({pruned_reason})")
    
    if verbose:
        print(f"\n{'='*70}")
        print("OPTIMIZATION RESULTS")
        print(f"{'='*70}")
        print(f"Total pairs checked: {pairs_checked}")
        print(f"Pairs with synergies: {pairs_with_synergies}")
        print(f"Total maximal synergies: {total_maximal}")
        print(f"Synergies pruned by pair domination: {synergies_pruned_by_pair}")
        print(f"Synergies pruned by target containment: {synergies_pruned_by_target}")
        total_pruned = synergies_pruned_by_pair + synergies_pruned_by_target
        pruning_pct = 100 * total_pruned / max(1, total_maximal)
        print(f"Total synergies pruned: {total_pruned} ({pruning_pct:.1f}%)")
        print(f"Distinct targets seen: {len(targets_seen)}")
        print(f"Candidate fundamental synergies found: {len(fundamental_synergies)}")
        
        # Additional statistics
        synergies_per_target = defaultdict(int)
        synergies_per_pair = defaultdict(int)
        for syn in fundamental_synergies:
            synergies_per_target[syn.plabel] += 1
            pair_key = tuple(sorted(syn.rlabel))
            synergies_per_pair[pair_key] += 1
        
        if synergies_per_target:
            max_syn_target = max(synergies_per_target.values())
            avg_syn_target = sum(synergies_per_target.values()) / len(synergies_per_target)
            print(f"Average fundamental synergies per target: {avg_syn_target:.2f}")
            print(f"Max fundamental synergies to one target: {max_syn_target}")
        
        if synergies_per_pair:
            max_syn_pair = max(synergies_per_pair.values())
            avg_syn_pair = sum(synergies_per_pair.values()) / len(synergies_per_pair)
            print(f"Average fundamental synergies per pair: {avg_syn_pair:.2f}")
            print(f"Max fundamental synergies from one pair: {max_syn_pair}")
        print(f"{'='*70}\n")
    
    # CRITICAL FIX: Final verification pass
    # The pruning heuristics are aggressive and may miss some cases
    # Apply strict fundamentality check to all candidates
    if verbose:
        print("Running final verification pass to ensure strict fundamentality...")
    
    verified_fundamental = []
    for syn in fundamental_synergies:
        if verify_synergy_is_fundamental(syn, hierarchy, verbose=False):
            verified_fundamental.append(syn)
        elif verbose:
            print(f"  Filtered out: {syn.rlabel[0]}+{syn.rlabel[1]}→{syn.plabel} (not truly fundamental)")
    
    if verbose:
        filtered_count = len(fundamental_synergies) - len(verified_fundamental)
        if filtered_count > 0:
            print(f"Filtered out {filtered_count} false positives")
        print(f"Final verified fundamental synergies: {len(verified_fundamental)}")
    
    return verified_fundamental


def verify_synergy_is_fundamental(synergy, hierarchy, verbose=False):
    """
    Verify that a synergy is truly fundamental according to Definition 19.
    
    A synergy E1 + E2 → E3 is fundamental iff there is NO synergy:
    - E1' + E2 → E3 where E1' ⊂ E1 (E1' is contained by E1)
    - E1 + E2' → E3 where E2' ⊂ E2 (E2' is contained by E2)
    
    Parameters
    ----------
    synergy : ERC_Synergy
        The synergy to verify
    hierarchy : ERC_Hierarchy
        The hierarchy object
    verbose : bool
        Print debug info
        
    Returns
    -------
    bool
        True if the synergy is truly fundamental
    """
    erc1_label, erc2_label = synergy.rlabel
    target_label = synergy.plabel
    
    # Get descendants (contained ERCs) for both reactants
    if hierarchy.graph:
        erc1_descendants = set(nx.descendants(hierarchy.graph, erc1_label))
        erc2_descendants = set(nx.descendants(hierarchy.graph, erc2_label))
    else:
        return True  # No hierarchy, so it's fundamental by default
    
    # Check if any descendant of erc1 (with erc2) produces the same target
    for desc1_label in erc1_descendants:
        desc1 = hierarchy.get_erc_by_label(desc1_label)
        erc2 = hierarchy.get_erc_by_label(erc2_label)
        
        # Get maximal synergies for this pair
        desc_maximal = get_maximal_synergies(desc1, erc2, hierarchy, hierarchy.RN)
        
        for desc_syn in desc_maximal:
            if desc_syn.plabel == target_label:
                if verbose:
                    print(f"Not fundamental: {desc1_label}+{erc2_label}→{target_label} exists")
                return False
    
    # Check if any descendant of erc2 (with erc1) produces the same target
    for desc2_label in erc2_descendants:
        erc1 = hierarchy.get_erc_by_label(erc1_label)
        desc2 = hierarchy.get_erc_by_label(desc2_label)
        
        # Get maximal synergies for this pair
        desc_maximal = get_maximal_synergies(erc1, desc2, hierarchy, hierarchy.RN)
        
        for desc_syn in desc_maximal:
            if desc_syn.plabel == target_label:
                if verbose:
                    print(f"Not fundamental: {erc1_label}+{desc2_label}→{target_label} exists")
                return False
    
    return True


def get_fundamental_synergies_optimized(erc1, erc2, hierarchy, RN, 
                                       maximal_synergies=None, verbose=False):
    """
    Find fundamental synergies between two specific ERCs with target-aware pruning.
    
    This is a helper function for the optimized algorithm, but can also be used
    standalone for checking a specific pair.
    
    Parameters
    ----------
    erc1, erc2 : ERC
        The two ERCs to check for fundamental synergies
    hierarchy : ERC_Hierarchy
        The hierarchy object
    RN : ReactionNetwork
        The reaction network
    maximal_synergies : list, optional
        Pre-computed maximal synergies. If None, will compute them.
    verbose : bool, optional
        Print debug information
        
    Returns
    -------
    list
        List of ERC_Synergy objects with synergy_type="fundamental"
    """
    if maximal_synergies is None:
        maximal_synergies = get_maximal_synergies(erc1, erc2, hierarchy, RN)
    
    if not maximal_synergies:
        return []
    
    fundamental_synergies = []
    
    # Build descendant sets for both ERCs
    if hierarchy.graph:
        erc1_descendants = set(nx.descendants(hierarchy.graph, erc1.label))
        erc1_descendants.add(erc1.label)
        
        erc2_descendants = set(nx.descendants(hierarchy.graph, erc2.label))
        erc2_descendants.add(erc2.label)
    else:
        erc1_descendants = {erc1.label}
        erc2_descendants = {erc2.label}
    
    # For each maximal synergy, check if it's fundamental
    for syn in maximal_synergies:
        target_label = syn.plabel
        is_fundamental = True
        
        # Check all combinations of descendants
        for desc1_label in erc1_descendants:
            for desc2_label in erc2_descendants:
                # Skip if it's the same pair
                if (set([desc1_label, desc2_label]) == set([erc1.label, erc2.label])):
                    continue
                
                # Skip if both labels are the same
                if desc1_label == desc2_label:
                    continue
                
                # Get the descendant ERCs
                desc1 = hierarchy.get_erc_by_label(desc1_label)
                desc2 = hierarchy.get_erc_by_label(desc2_label)
                
                if desc1 is None or desc2 is None:
                    continue
                
                # Check if descendants produce the same target
                desc_maximal = get_maximal_synergies(desc1, desc2, hierarchy, RN)
                
                for desc_syn in desc_maximal:
                    if desc_syn.plabel == target_label:
                        is_fundamental = False
                        if verbose:
                            print(f"  Non-fundamental: {erc1.label}+{erc2.label}→{target_label} "
                                  f"because {desc1_label}+{desc2_label}→{target_label}")
                        break
                
                if not is_fundamental:
                    break
            
            if not is_fundamental:
                break
        
        if is_fundamental:
            fundamental_synergies.append(
                ERC_Synergy(syn.reactants, syn.product, "fundamental")
            )
    
    return fundamental_synergies


# ============================================================================
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