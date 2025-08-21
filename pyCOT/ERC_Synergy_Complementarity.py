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

def get_fundamental_synergies_brute_force(erc1, erc2, hierarchy, RN, maximal_synergies=None):
    """    
    A maximal synergy ERC1 + ERC2 → Target is fundamental if there's no synergy 
    MoreFundamental1 + MoreFundamental2 → Target where the reactants are more fundamental.  
    CRITICAL: Fundamentality is checked SEPARATELY for each target!
    """
    if maximal_synergies==None:
        maximal_synergies = get_maximal_synergies(erc1, erc2, hierarchy, RN)
    if not maximal_synergies:
        return []
    
    fundamental_synergies = []
    
    for syn in maximal_synergies:
        is_fundamental = True
        target_label = syn.plabel
        
        # CORRECTED: Check fundamentality specifically for THIS target
        # Look for any other pair that can produce the SAME target with more fundamental bases
        for other_erc1 in hierarchy.ercs:
            for other_erc2 in hierarchy.ercs:
                if other_erc1 == other_erc2:
                    continue
                    
                # Skip if it's the same pair
                if set([other_erc1.label, other_erc2.label]) == set([erc1.label, erc2.label]):
                    continue
                
                # Get maximal synergies for this other pair
                other_maximal = get_maximal_synergies(other_erc1, other_erc2, hierarchy, RN)
                
                for other_syn in other_maximal:
                    # CRITICAL: Only compare synergies that produce the SAME target
                    if other_syn.plabel == target_label:
                        # Check if other reactants are more fundamental using cached containment
                        more_fundamental = False
                        
                        # Case 1: other_erc1 is more fundamental than erc1, other_erc2 same as erc2
                        if (other_erc1 in hierarchy.get_contained(erc1) and other_erc2.label == erc2.label):
                            more_fundamental = True
                        # Case 2: other_erc1 same as erc1, other_erc2 is more fundamental than erc2
                        elif (other_erc1.label == erc1.label and other_erc2 in hierarchy.get_contained(erc2)):
                            more_fundamental = True
                        # Case 3: other_erc1 is more fundamental than erc2, other_erc2 same as erc1 (commutativity)
                        elif (other_erc1 in hierarchy.get_contained(erc2) and other_erc2.label == erc1.label):
                            more_fundamental = True
                        # Case 4: other_erc1 same as erc2, other_erc2 is more fundamental than erc1 (commutativity)
                        elif (other_erc1.label == erc2.label and other_erc2 in hierarchy.get_contained(erc1)):
                            more_fundamental = True
                        
                        if more_fundamental:
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
# SYNERGY COMPUTATION EFFICIENT ALGORITHM
# ============================================================================

def build_base_target_relationships_efficient(ercs, hierarchy, RN, verbose=False):
    """
    CORRECTED HIERARCHICAL: Build base-target relationships with proper level-based ordering.
    
    Key improvements:
    1. Process targets by hierarchy level (top-down: highest level first)
    2. For each target level, explore bases bottom-up (lowest level first) 
    3. Apply fundamentality optimization: if base1->base2 and both have same overlap with target, keep only base2
    4. Proper transitive pruning respecting hierarchy levels
    """
    if verbose:
        print("Building base-target relationships with CORRECTED hierarchical optimization...")
    
    base_target_map = defaultdict(list)
    partial_overlap_checks = 0
    valid_pairs = 0
    fundamentality_pruned = 0
    
    # Get hierarchy levels - critical for proper ordering
    if hierarchy.graph:
        levels = ERC.get_node_levels(hierarchy.graph)
        max_level = max(levels.values()) if levels else 0
        
        # Group ERCs by level for systematic processing
        ercs_by_level = defaultdict(list)
        for erc in ercs:
            level = levels.get(erc.label, 0)
            ercs_by_level[level].append(erc)
    else:
        levels = {erc.label: 0 for erc in ercs}
        max_level = 0
        ercs_by_level = {0: ercs}
    
    if verbose:
        print(f"  Hierarchy levels: 0 to {max_level}")
        for level in range(max_level + 1):
            level_ercs = ercs_by_level.get(level, [])
            print(f"    Level {level}: {len(level_ercs)} ERCs = {[erc.label for erc in level_ercs]}")
    
    # Track fundamentality relationships for optimization
    fundamentality_cache = {}  # (base1, base2, target) -> keep_base2_only
    
    # CORRECTED ALGORITHM: Process targets TOP-DOWN (highest level first)
    for target_level in range(max_level, -1, -1):
        targets_at_level = ercs_by_level.get(target_level, [])
        
        if not targets_at_level:
            continue
            
        if verbose:
            print(f"\n  🎯 Processing TARGET LEVEL {target_level}: {[t.label for t in targets_at_level]}")
        
        for target in targets_at_level:
            if verbose:
                print(f"\n    Target: {target.label} (level {target_level})")
            
            # Get ancestors (ERCs that contain target) for efficiency
            ancestors = set()
            if hierarchy.graph and target.label in hierarchy.graph:
                ancestors = set(nx.ancestors(hierarchy.graph, target.label))
            
            # Track valid bases and their overlaps for fundamentality optimization
            target_bases_with_overlap = []
            
            # CORRECTED: Explore bases BOTTOM-UP (lowest level first, most fundamental)
            for base_level in range(max_level + 1):
                bases_at_level = ercs_by_level.get(base_level, [])
                
                if verbose and bases_at_level:
                    print(f"      Checking bases at level {base_level}: {[b.label for b in bases_at_level]}")
                
                for base in bases_at_level:
                    # Skip if base is the target itself
                    if base.label == target.label:
                        continue
                    
                    # Skip if base contains target (impossible synergy)
                    if base.label in ancestors:
                        continue
                    
                    # Check partial overlap with target's minimal generators
                    partial_overlap_checks += 1
                    if has_partial_overlap_with_generators(base, target, RN):
                        
                        # Store base with its overlap info for fundamentality check
                        base_overlap_info = {
                            'base': base,
                            'level': base_level,
                            'overlap_signature': get_overlap_signature(base, target, RN)
                        }
                        target_bases_with_overlap.append(base_overlap_info)
                        
                    else:
                        if verbose:
                            print(f"        ✗ {base.label}: no partial overlap")
            
            # FUNDAMENTALITY OPTIMIZATION: Remove non-fundamental bases
            if verbose and len(target_bases_with_overlap) > 1:
                print(f"      Applying fundamentality optimization to {len(target_bases_with_overlap)} potential bases...")
            
            fundamental_bases = apply_fundamentality_optimization(
                target_bases_with_overlap, hierarchy, verbose)
                
            fundamentality_pruned += len(target_bases_with_overlap) - len(fundamental_bases)
            
            # Add fundamental bases to result
            for base_info in fundamental_bases:
                base_target_map[target.label].append(base_info['base'].label)
                valid_pairs += 1
                
                if verbose:
                    print(f"        ✓ {base_info['base'].label} (level {base_info['level']}): fundamental base")
    
    if verbose:
        total_possible = len(ercs) * (len(ercs) - 1)
        print(f"\n  📊 CORRECTED ALGORITHM STATISTICS:")
        print(f"    Partial overlap checks: {partial_overlap_checks}")
        print(f"    Valid base-target pairs: {valid_pairs}")
        print(f"    Fundamentality pruned: {fundamentality_pruned}")
        print(f"    Reduction: {(1 - valid_pairs/total_possible)*100:.1f}%")
        
        # Show level distribution
        level_stats = defaultdict(int)
        for target_label, bases in base_target_map.items():
            target_level = levels.get(target_label, 0)
            level_stats[target_level] += len(bases)
        
        print(f"    Base-target pairs by target level:")
        for level in sorted(level_stats.keys(), reverse=True):
            print(f"      Level {level}: {level_stats[level]} pairs")
    
    return base_target_map


def get_overlap_signature(base, target, RN):
    """
    Get a signature representing the overlap pattern between base and target.
    
    This is used for the fundamentality optimization: if two bases have the 
    same overlap signature with a target, we can keep only the more fundamental one.
    """
    base_closure = base.get_closure_names(RN)
    
    overlap_patterns = []
    
    for i, generator in enumerate(target.min_generators):
        gen_species = set(species_list_to_names(generator))
        intersection = base_closure & gen_species
        
        # Create a signature for this generator overlap
        if intersection:
            # Sort for consistency
            overlap_pattern = tuple(sorted(intersection))
            coverage_ratio = len(intersection) / len(gen_species)
            overlap_patterns.append((i, overlap_pattern, coverage_ratio))
    
    # Return a tuple that represents the overlap signature
    return tuple(sorted(overlap_patterns))


def apply_fundamentality_optimization(bases_with_overlap, hierarchy, verbose=False):
    """
    Apply fundamentality optimization: if base1->base2 and both have the same 
    overlap with target, keep only base2 (more fundamental).
    
    Returns the list of fundamental bases (after pruning non-fundamental ones).
    """
    if len(bases_with_overlap) <= 1:
        return bases_with_overlap
    
    # Group bases by overlap signature
    overlap_groups = defaultdict(list)
    for base_info in bases_with_overlap:
        signature = base_info['overlap_signature']
        overlap_groups[signature].append(base_info)
    
    fundamental_bases = []
    
    for signature, bases_group in overlap_groups.items():
        if len(bases_group) == 1:
            # Only one base with this signature - keep it
            fundamental_bases.extend(bases_group)
        else:
            # Multiple bases with same overlap - keep only most fundamental ones
            if verbose:
                base_labels = [b['base'].label for b in bases_group]
                print(f"          Same overlap signature: {base_labels}")
            
            # Find most fundamental bases (those not contained by others in the group)
            most_fundamental = []
            
            for base_info in bases_group:
                base_erc = base_info['base']
                is_fundamental = True
                
                # Check if this base is contained by any other base in the group
                for other_info in bases_group:
                    if other_info == base_info:
                        continue
                    
                    other_erc = other_info['base']
                    
                    # Check containment using hierarchy
                    if hierarchy.graph and hierarchy.graph.has_edge(other_erc.label, base_erc.label):
                        # other_erc contains base_erc, so base_erc is more fundamental
                        continue
                    elif hierarchy.graph and hierarchy.graph.has_edge(base_erc.label, other_erc.label):
                        # base_erc contains other_erc, so base_erc is less fundamental
                        is_fundamental = False
                        break
                    else:
                        # No direct containment - compare by level (lower level = more fundamental)
                        if other_info['level'] < base_info['level']:
                            is_fundamental = False
                            break
                
                if is_fundamental:
                    most_fundamental.append(base_info)
            
            if verbose and len(most_fundamental) < len(bases_group):
                pruned_labels = [b['base'].label for b in bases_group if b not in most_fundamental]
                kept_labels = [b['base'].label for b in most_fundamental]
                print(f"          Pruned non-fundamental: {pruned_labels}")
                print(f"          Kept fundamental: {kept_labels}")
            
            fundamental_bases.extend(most_fundamental)
    
    return fundamental_bases

def generate_target_pairs_efficient(base_target_map, hierarchy, verbose=False):
    """
    OPTIMIZED: Generate target pairs with enhanced pruning.
    """
    if verbose:
        print("Generating target-centric base pairs with optimization...")
    
    target_pairs = defaultdict(list)
    total_pairs = 0
    hierarchy_pruned = 0
    
    for target, potential_bases in base_target_map.items():
        if len(potential_bases) < 2:
            continue  # Need at least 2 bases for a synergy
        
        # OPTIMIZATION: Order bases by hierarchy level for systematic exploration
        if hierarchy.graph:
            levels = ERC.get_node_levels(hierarchy.graph)
            potential_bases = sorted(potential_bases, key=lambda b: levels.get(b, 0))
        
        # Generate pairs ONLY from bases that can contribute to THIS target
        pairs_for_target = 0
        for i, base1 in enumerate(potential_bases):
            for base2 in potential_bases[i+1:]:
                # Check that bases don't contain each other using hierarchy
                if hierarchy.graph:
                    if (hierarchy.graph.has_edge(base1, base2) or 
                        hierarchy.graph.has_edge(base2, base1)):
                        hierarchy_pruned += 1
                        continue
                
                # Store in canonical order
                pair = (base1, base2) if base1 < base2 else (base2, base1)
                target_pairs[target].append(pair)
                pairs_for_target += 1
                total_pairs += 1
        
        if verbose and pairs_for_target > 0:
            print(f"  Target {target}: {pairs_for_target} pairs from {len(potential_bases)} bases")
    
    if verbose:
        print(f"  Total base pairs generated: {total_pairs}")
        print(f"  Pairs pruned by hierarchy: {hierarchy_pruned}")
        avg_pairs = total_pairs / len(target_pairs) if target_pairs else 0
        print(f"  Average pairs per target: {avg_pairs:.1f}")
    
    return target_pairs

def apply_synergy_constraints(base1_label, base2_label, target_label, 
                            target_pairs, constraints, base_target_map, hierarchy):
    """
    OPTIMIZED: Apply constraints with enhanced fundamentality pruning.
    """
    total_pruned = 0
    
    if not hierarchy.graph:
        return total_pruned
    
    # Get hierarchy levels for systematic pruning
    levels = ERC.get_node_levels(hierarchy.graph)
    
    # Get ancestors of bases (for non-minimality pruning)
    ancestors_base1 = set()
    ancestors_base2 = set()
    if base1_label in hierarchy.graph:
        ancestors_base1 = set(nx.ancestors(hierarchy.graph, base1_label))
    if base2_label in hierarchy.graph:
        ancestors_base2 = set(nx.ancestors(hierarchy.graph, base2_label))
    
    # 1. ENHANCED within-target pruning
    for pair in list(target_pairs[target_label]):
        if pair in constraints[target_label]:
            continue
            
        b1, b2 = pair
        should_prune = False
        
        # More systematic non-minimality checking
        b1_level = levels.get(b1, 0)
        b2_level = levels.get(b2, 0)
        base1_level = levels.get(base1_label, 0)
        base2_level = levels.get(base2_label, 0)
        
        # Prune if we found a more fundamental combination
        if ((b1 in ancestors_base1 and b2 == base2_label) or
            (b1 == base1_label and b2 in ancestors_base2) or
            (b1 in ancestors_base1 and b2 in ancestors_base2) or
            (b2 in ancestors_base1 and b1 == base2_label) or
            (b2 == base1_label and b1 in ancestors_base2) or
            (b2 in ancestors_base1 and b1 in ancestors_base2)):
            should_prune = True
            
        if should_prune:
            constraints[target_label].add(pair)
            total_pruned += 1
    
    # 2. ENHANCED cross-target pruning with level awareness
    descendants = set()
    if target_label in hierarchy.graph:
        descendants = set(nx.descendants(hierarchy.graph, target_label))
    
    # Sort descendants by level for systematic processing
    if descendants:
        desc_by_level = sorted(descendants, key=lambda d: levels.get(d, 0), reverse=True)
        
        for desc_target in desc_by_level:
            if desc_target not in target_pairs:
                continue
            
            # Add the found pair as constraint for descendant targets
            pair = (base1_label, base2_label) if base1_label < base2_label else (base2_label, base1_label)
            if pair not in constraints[desc_target]:
                constraints[desc_target].add(pair)
                total_pruned += 1
            
            # Aggressively prune non-minimal versions in descendant targets
            for test_pair in list(target_pairs[desc_target]):
                if test_pair in constraints[desc_target]:
                    continue
                    
                b1, b2 = test_pair
                should_prune = False
                
                # Enhanced fundamentality check using levels
                if (((b1 == base1_label or b1 in ancestors_base1) and
                     (b2 == base2_label or b2 in ancestors_base2)) or
                    ((b1 == base2_label or b1 in ancestors_base2) and
                     (b2 == base1_label or b2 in ancestors_base1))):
                    should_prune = True
                    
                if should_prune:
                    constraints[desc_target].add(test_pair)
                    total_pruned += 1
    
    return total_pruned

def explore_targets_with_proper_hierarchy(target_pairs, base_target_map, ercs, hierarchy, RN, verbose=True):
    """
    CORRECTED: Proper hierarchical exploration respecting level-based ordering.
    
    Key fixes:
    1. Base pairs explored systematically by hierarchy levels (bottom-up)
    2. Targets explored by hierarchy levels (top-down) 
    3. No arbitrary sorting that breaks hierarchy logic
    """
    if verbose:
        print("Exploring with PROPER hierarchical ordering...")
    
    # Create label to ERC mapping
    erc_by_label = {erc.label: erc for erc in ercs}
    
    # Get hierarchy levels
    if hierarchy.graph:
        levels = ERC.get_node_levels(hierarchy.graph)
        max_level = max(levels.values()) if levels else 0
    else:
        levels = {erc.label: 0 for erc in ercs}
        max_level = 0
    
    constraints = defaultdict(set)
    discovered_synergies = []
    
    # Statistics
    stats = {
        'pairs_checked': 0,
        'pairs_skipped_constraint': 0,
        'pairs_skipped_generator': 0,
        'synergies_found': 0,
        'level_combinations_checked': 0
    }
    
    if verbose:
        print(f"  Hierarchy depth: {max_level} levels")
        print(f"  Target pairs available for: {list(target_pairs.keys())}")
    
    # CORRECTED: Systematic exploration by hierarchy levels
    # Targets: TOP-DOWN (highest level first - most maximal)
    for target_level in range(max_level, -1, -1):
        
        # Get targets at this level
        targets_at_level = [label for label in target_pairs.keys() 
                           if levels.get(label, 0) == target_level]
        
        if not targets_at_level:
            continue
            
        if verbose:
            print(f"\n  🎯 TARGET LEVEL {target_level}: {targets_at_level}")
        
        for target_label in targets_at_level:
            pairs = target_pairs[target_label]
            if not pairs:
                continue
            
            if verbose:
                print(f"\n    Target {target_label} (level {target_level})")
                print(f"      Initial pairs: {len(pairs)}")
            
            # Filter out constrained pairs
            valid_pairs = [p for p in pairs if p not in constraints[target_label]]
            stats['pairs_skipped_constraint'] += len(pairs) - len(valid_pairs)
            
            if verbose:
                print(f"      Valid pairs after constraints: {len(valid_pairs)}")
            
            # CORRECTED: Systematic base pair exploration by levels
            # Group pairs by their combined hierarchy levels
            pair_levels = {}
            for pair in valid_pairs:
                base1_level = levels.get(pair[0], 0)
                base2_level = levels.get(pair[1], 0)
                combined_level = base1_level + base2_level  # Lower is more fundamental
                
                if combined_level not in pair_levels:
                    pair_levels[combined_level] = []
                pair_levels[combined_level].append(pair)
            
            # Process pairs from most fundamental to least fundamental
            for combined_level in sorted(pair_levels.keys()):
                level_pairs = pair_levels[combined_level]
                
                if verbose:
                    print(f"        Processing {len(level_pairs)} pairs at combined level {combined_level}")
                
                stats['level_combinations_checked'] += 1
                
                # Check each pair at this level
                for base1_label, base2_label in level_pairs:
                    stats['pairs_checked'] += 1
                    
                    base1 = erc_by_label[base1_label]
                    base2 = erc_by_label[base2_label]
                    target_erc = erc_by_label[target_label]
                    
                    # Generator coverage check (fast pre-filter)
                    if not has_generator_coverage(base1, base2, target_erc, RN):
                        stats['pairs_skipped_generator'] += 1
                        continue
                    
                    
                # Check if any synergy produces this specific target
                # Optional debug flag at the top of your function or as parameter
                verify_with_brute_force = False  # Set True only for testing correctness

                if verify_with_brute_force:
                    # Debug-only: confirm using brute force fundamental synergy method
                    brute_synergies = get_fundamental_synergies_brute_force(base1, base2, hierarchy, RN)
                    for synergy in brute_synergies:
                        if synergy.plabel == target_label:
                            target_synergy_found = True
                            discovered_synergies.append(synergy)
                            stats['synergies_found'] += 1
                            break
                else:
                    # Use current algorithm's result directly (construct ERC_Synergy here)
                    synergy_obj = ERC_Synergy([base1, base2], hierarchy.get_erc_by_label(target_label), "fundamental")
                    discovered_synergies.append(synergy_obj)
                    stats['synergies_found'] += 1
                    target_synergy_found = True

                if target_synergy_found:
                    if verbose:
                        base1_level = levels.get(base1_label, 0)
                        base2_level = levels.get(base2_label, 0)
                        print(f"        ✓ SYNERGY: {base1_label}(L{base1_level}) + {base2_label}(L{base2_level}) → {target_label}(L{target_level})")

                    # Apply constraints with proper hierarchy awareness
                    pruned = apply_hierarchical_constraints(
                        base1_label, base2_label, target_label,
                        target_pairs, constraints, base_target_map, hierarchy, levels)

                    if verbose and pruned > 0:
                        print(f"          Pruned {pruned} non-fundamental combinations")

                    # Skip higher-level combinations for this target
                    remaining_levels = [lv for lv in pair_levels.keys() if lv > combined_level]
                    if remaining_levels and verbose:
                        remaining_pairs = sum(len(pair_levels[lv]) for lv in remaining_levels)
                        print(f"          Skipping {remaining_pairs} less fundamental pairs")
                    break
                    
    
    if verbose:
        print(f"\n📊 Hierarchical Exploration Statistics:")
        print(f"  Level combinations checked: {stats['level_combinations_checked']}")
        print(f"  Pairs checked: {stats['pairs_checked']}")
        print(f"  Pairs skipped by constraints: {stats['pairs_skipped_constraint']}")
        print(f"  Pairs skipped by generator check: {stats['pairs_skipped_generator']}")
        print(f"  Synergies found: {stats['synergies_found']}")
        
        total_possible = sum(len(pairs) for pairs in target_pairs.values())
        if total_possible > 0:
            efficiency = (1 - stats['pairs_checked']/total_possible) * 100
            print(f"  Efficiency: {efficiency:.1f}% reduction in checks")
    
    return discovered_synergies

def apply_hierarchical_constraints(base1_label, base2_label, target_label, 
                                 target_pairs, constraints, base_target_map, hierarchy, levels):
    """
    CORRECTED: Apply constraints with proper hierarchy level awareness.
    """
    total_pruned = 0
    
    if not hierarchy.graph:
        return total_pruned
    
    # Get hierarchy information
    base1_level = levels.get(base1_label, 0)
    base2_level = levels.get(base2_label, 0)
    target_level = levels.get(target_label, 0)
    
    # Get ancestors for containment-based pruning
    ancestors_base1 = set()
    ancestors_base2 = set()
    if base1_label in hierarchy.graph:
        ancestors_base1 = set(nx.ancestors(hierarchy.graph, base1_label))
    if base2_label in hierarchy.graph:
        ancestors_base2 = set(nx.ancestors(hierarchy.graph, base2_label))
    
    # 1. Within-target pruning: remove non-fundamental base combinations
    for pair in list(target_pairs[target_label]):
        if pair in constraints[target_label]:
            continue
            
        b1, b2 = pair
        b1_level = levels.get(b1, 0)
        b2_level = levels.get(b2, 0)
        
        should_prune = False
        
        # Prune if this pair uses less fundamental bases than the found synergy
        if ((b1 in ancestors_base1 or b1_level > base1_level) and 
            (b2 == base2_label or b2_level >= base2_level)):
            should_prune = True
        elif ((b2 in ancestors_base1 or b2_level > base1_level) and 
              (b1 == base2_label or b1_level >= base2_level)):
            should_prune = True
        elif ((b1 == base1_label or b1_level >= base1_level) and 
              (b2 in ancestors_base2 or b2_level > base2_level)):
            should_prune = True
        elif ((b2 == base1_label or b2_level >= base1_level) and 
              (b1 in ancestors_base2 or b1_level > base2_level)):
            should_prune = True
            
        if should_prune:
            constraints[target_label].add(pair)
            total_pruned += 1
    
    # 2. Cross-target pruning: for targets contained in current target
    descendants = set()
    if target_label in hierarchy.graph:
        descendants = set(nx.descendants(hierarchy.graph, target_label))
    
    for desc_target in descendants:
        if desc_target not in target_pairs:
            continue
        
        desc_level = levels.get(desc_target, 0)
        
        # Only prune if the descendant is at a lower level (less maximal)
        if desc_level < target_level:
            pair = (base1_label, base2_label) if base1_label < base2_label else (base2_label, base1_label)
            if pair not in constraints[desc_target]:
                constraints[desc_target].add(pair)
                total_pruned += 1
            
            # Also prune non-fundamental versions in descendant targets
            for test_pair in list(target_pairs[desc_target]):
                if test_pair in constraints[desc_target]:
                    continue
                    
                tb1, tb2 = test_pair
                tb1_level = levels.get(tb1, 0)
                tb2_level = levels.get(tb2, 0)
                
                # Prune if this uses less fundamental bases
                if (((tb1 == base1_label or tb1_level >= base1_level) and
                     (tb2 == base2_label or tb2_level >= base2_level)) or
                    ((tb1 == base2_label or tb1_level >= base2_level) and
                     (tb2 == base1_label or tb2_level >= base1_level))):
                    constraints[desc_target].add(test_pair)
                    total_pruned += 1
    
    return total_pruned
def get_fundamental_synergies(ercs, hierarchy, RN, verbose=True):
    """
    CORRECTED: Efficient algorithm with proper hierarchical exploration.
    """
    if verbose:
        print("="*70)
        print("CORRECTED HIERARCHICAL SYNERGY FINDER")
        print("="*70)   
        print(f"Network: {len(ercs)} ERCs, {len(RN.species())} species")
        if hierarchy.graph:
            print(f"Hierarchy: {hierarchy.graph.number_of_edges()} containment relations")
    
    import time
    start_time = time.time()
    
    # Step 1: Build base-target relationships (keep existing logic)
    base_target_map = build_base_target_relationships_efficient(ercs, hierarchy, RN, verbose)
    
    # Step 2: Generate target-centric structure (keep existing logic)
    target_pairs = generate_target_pairs_efficient(base_target_map, hierarchy, verbose)
    
    # Step 3: CORRECTED hierarchical exploration
    fundamental_synergies = explore_targets_with_proper_hierarchy(
        target_pairs, base_target_map, ercs, hierarchy, RN, verbose)
    
    end_time = time.time()
    
    if verbose:
        print(f"\n" + "="*70)
        print("CORRECTED ALGORITHM RESULTS")
        print("="*70)
        print(f"Fundamental synergies found: {len(fundamental_synergies)}")
        print(f"Time elapsed: {end_time - start_time:.3f} seconds")
        
        if fundamental_synergies:
            print(f"\nFundamental Synergies:")
            for i, synergy in enumerate(fundamental_synergies, 1):
                print(f"  {i}. {synergy.rlabel[0]} + {synergy.rlabel[1]} → {synergy.plabel}")
    
    return fundamental_synergies


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