"""
Semantic Partitioning Module for Reaction Network Analysis

This module provides the N-category semantic partitioning infrastructure
that enables semantic-level analysis of reaction networks.

Part 1 of the Challenge-Centered Semantic Pipeline Implementation
"""

import numpy as np
from typing import Dict, List, Set, Optional, Tuple
from dataclasses import dataclass, field


@dataclass
class SemanticPartition:
    """
    A semantic partition of a reaction network's species into meaningful categories.
    
    Attributes:
        categories: List of category names
        species_to_categories: Mapping from species name to list of categories it belongs to
        category_to_species: Mapping from category name to list of species in it
        species_indices: Mapping from species name to index in stoichiometric matrix
        category_indices: Mapping from category name to list of indices in stoichiometric matrix
    """
    categories: List[str]
    species_to_categories: Dict[str, List[str]]
    category_to_species: Dict[str, List[str]]
    species_indices: Dict[str, int]
    category_indices: Dict[str, List[int]]
    
    def __repr__(self):
        return (f"SemanticPartition(\n"
                f"  categories={self.categories},\n"
                f"  n_species={len(self.species_indices)},\n"
                f"  overlapping_species={sum(1 for cats in self.species_to_categories.values() if len(cats) > 1)}\n"
                f")")


def define_semantic_categories(
    species_list: List[str],
    category_definitions: Dict[str, List[str]]
) -> SemanticPartition:
    """
    Create a semantic partition from species list and category definitions.
    
    This function maps species into semantic categories, allowing species to belong
    to multiple categories (overlapping semantics). It creates both forward and
    reverse mappings for efficient queries.
    
    Args:
        species_list: List of species names from the reaction network (e.g., from rn.species())
        category_definitions: Dictionary mapping category names to lists of species
            Example: {
                'conflict': ['A_v', 'B_v', 'D_A', 'D_B'],
                'peace': ['A_p', 'B_p', 'G_A', 'G_B'],
                'perception': ['iD_A', 'iD_B', 'iG_A', 'iG_B']
            }
    
    Returns:
        SemanticPartition object containing all mappings and indices
        
    Raises:
        ValueError: If a species in category_definitions is not in species_list
        ValueError: If category_definitions is empty
        
    Example:
        >>> species = ['A_v', 'A_p', 'B_v', 'B_p', 'iG_A', 'iG_B']
        >>> categories = {
        ...     'conflict': ['A_v', 'B_v'],
        ...     'peace': ['A_p', 'B_p'],
        ...     'perception': ['iG_A', 'iG_B']
        ... }
        >>> partition = define_semantic_categories(species, categories)
        >>> print(partition.category_to_species['conflict'])
        ['A_v', 'B_v']
    """
    if not category_definitions:
        raise ValueError("category_definitions cannot be empty")
    
    # Create species index mapping (matches order in stoichiometric matrix)
    species_indices = {species: idx for idx, species in enumerate(species_list)}
    
    # Validate that all species in categories exist
    for category, species_in_cat in category_definitions.items():
        for species in species_in_cat:
            if species not in species_indices:
                raise ValueError(
                    f"Species '{species}' in category '{category}' "
                    f"not found in species_list"
                )
    
    # Build forward mapping: species -> categories
    species_to_categories = {species: [] for species in species_list}
    for category, species_in_cat in category_definitions.items():
        for species in species_in_cat:
            species_to_categories[species].append(category)
    
    # Build reverse mapping: category -> species (copy to preserve immutability)
    category_to_species = {
        category: list(species_in_cat)
        for category, species_in_cat in category_definitions.items()
    }
    
    # Build category -> indices mapping
    category_indices = {
        category: [species_indices[species] for species in species_in_cat]
        for category, species_in_cat in category_definitions.items()
    }
    
    # Get list of categories
    categories = list(category_definitions.keys())
    
    return SemanticPartition(
        categories=categories,
        species_to_categories=species_to_categories,
        category_to_species=category_to_species,
        species_indices=species_indices,
        category_indices=category_indices
    )


def get_species_indices_for_category(
    semantic_partition: SemanticPartition,
    category: str
) -> List[int]:
    """
    Get the list of stoichiometric matrix indices for species in a category.
    
    Args:
        semantic_partition: A SemanticPartition object
        category: Name of the category
        
    Returns:
        List of integer indices corresponding to species in this category
        
    Raises:
        KeyError: If category is not in the partition
        
    Example:
        >>> indices = get_species_indices_for_category(partition, 'conflict')
        >>> print(indices)
        [0, 2]  # Indices of A_v and B_v
    """
    if category not in semantic_partition.category_indices:
        raise KeyError(f"Category '{category}' not found in semantic partition")
    
    return semantic_partition.category_indices[category]


def get_categories_for_species(
    semantic_partition: SemanticPartition,
    species: str
) -> List[str]:
    """
    Get the list of categories that a species belongs to.
    
    Args:
        semantic_partition: A SemanticPartition object
        species: Name of the species
        
    Returns:
        List of category names this species belongs to (may be empty if uncategorized)
        
    Raises:
        KeyError: If species is not in the partition
        
    Example:
        >>> categories = get_categories_for_species(partition, 'A_v')
        >>> print(categories)
        ['conflict']
    """
    if species not in semantic_partition.species_to_categories:
        raise KeyError(f"Species '{species}' not found in semantic partition")
    
    return semantic_partition.species_to_categories[species]


def validate_semantic_partition(
    semantic_partition: SemanticPartition,
    S: np.ndarray,
    verbose: bool = True
) -> Dict[str, any]:
    """
    Validate a semantic partition against a stoichiometric matrix.
    
    Performs several checks:
    1. All indices are within valid bounds for S
    2. No duplicate indices within a category
    3. Consistency between species mappings
    4. Coverage statistics
    
    Args:
        semantic_partition: A SemanticPartition object
        S: Stoichiometric matrix (n_species × n_reactions)
        verbose: If True, print validation report
        
    Returns:
        Dictionary containing validation results:
            - 'valid': bool, overall validation status
            - 'errors': list of error messages
            - 'warnings': list of warning messages
            - 'stats': dictionary of statistics
            
    Example:
        >>> report = validate_semantic_partition(partition, S)
        >>> if report['valid']:
        ...     print("Partition is valid!")
    """
    errors = []
    warnings = []
    stats = {}
    
    n_species = S.shape[0]
    
    # Check 1: Validate indices are in bounds
    for category, indices in semantic_partition.category_indices.items():
        for idx in indices:
            if idx < 0 or idx >= n_species:
                errors.append(
                    f"Category '{category}': index {idx} out of bounds "
                    f"(valid range: 0-{n_species-1})"
                )
    
    # Check 2: Check for duplicate indices within categories
    for category, indices in semantic_partition.category_indices.items():
        if len(indices) != len(set(indices)):
            duplicates = [idx for idx in set(indices) if indices.count(idx) > 1]
            errors.append(
                f"Category '{category}': duplicate indices found: {duplicates}"
            )
    
    # Check 3: Verify consistency between mappings
    for species, categories in semantic_partition.species_to_categories.items():
        for category in categories:
            if species not in semantic_partition.category_to_species[category]:
                errors.append(
                    f"Inconsistency: species '{species}' lists category '{category}' "
                    f"but category doesn't list species"
                )
    
    # Check 4: Coverage statistics
    all_categorized_species = set()
    for species_list in semantic_partition.category_to_species.values():
        all_categorized_species.update(species_list)
    
    n_categorized = len(all_categorized_species)
    n_total = len(semantic_partition.species_indices)
    coverage_pct = 100 * n_categorized / n_total if n_total > 0 else 0
    
    stats['n_species'] = n_total
    stats['n_categorized'] = n_categorized
    stats['n_uncategorized'] = n_total - n_categorized
    stats['coverage_percent'] = coverage_pct
    stats['n_categories'] = len(semantic_partition.categories)
    
    # Check for uncategorized species
    uncategorized = []
    for species in semantic_partition.species_indices:
        if not semantic_partition.species_to_categories[species]:
            uncategorized.append(species)
    
    if uncategorized:
        warnings.append(
            f"Found {len(uncategorized)} uncategorized species: {uncategorized}"
        )
    
    # Check for overlapping species
    overlapping = [
        species for species, cats in semantic_partition.species_to_categories.items()
        if len(cats) > 1
    ]
    stats['n_overlapping'] = len(overlapping)
    if overlapping:
        stats['overlapping_species'] = overlapping
    
    # Compute category sizes
    category_sizes = {
        cat: len(species_list)
        for cat, species_list in semantic_partition.category_to_species.items()
    }
    stats['category_sizes'] = category_sizes
    
    # Overall validation
    valid = len(errors) == 0
    
    # Print report if requested
    if verbose:
        print("=" * 70)
        print("SEMANTIC PARTITION VALIDATION REPORT")
        print("=" * 70)
        print(f"\nStatus: {'✓ VALID' if valid else '✗ INVALID'}")
        
        if errors:
            print(f"\n❌ Errors ({len(errors)}):")
            for error in errors:
                print(f"  - {error}")
        
        if warnings:
            print(f"\n⚠️  Warnings ({len(warnings)}):")
            for warning in warnings:
                print(f"  - {warning}")
        
        print(f"\n📊 Statistics:")
        print(f"  Total species: {stats['n_species']}")
        print(f"  Categorized species: {stats['n_categorized']} ({stats['coverage_percent']:.1f}%)")
        print(f"  Uncategorized species: {stats['n_uncategorized']}")
        print(f"  Number of categories: {stats['n_categories']}")
        print(f"  Overlapping species: {stats['n_overlapping']}")
        
        print(f"\n📦 Category Sizes:")
        for category, size in category_sizes.items():
            species_list = semantic_partition.category_to_species[category]
            print(f"  {category}: {size} species {species_list}")
        
        print("=" * 70)
    
    return {
        'valid': valid,
        'errors': errors,
        'warnings': warnings,
        'stats': stats
    }


def get_semantic_impact_template(semantic_partition: SemanticPartition) -> Dict[str, Dict]:
    """
    Create a template dictionary for recording semantic impacts.
    
    This is a utility function for downstream analysis (Part 2).
    Creates a nested dictionary structure pre-populated with categories.
    
    Args:
        semantic_partition: A SemanticPartition object
        
    Returns:
        Dictionary with structure:
        {
            category_name: {
                'consumed': [],
                'produced': [],
                'net': 0.0
            }
        }
        
    Example:
        >>> template = get_semantic_impact_template(partition)
        >>> template['conflict']['consumed'] = [0, 2]  # indices
        >>> template['conflict']['net'] = -1.5
    """
    template = {}
    for category in semantic_partition.categories:
        template[category] = {
            'consumed': [],
            'produced': [],
            'net': 0.0
        }
    return template


def print_semantic_partition_summary(semantic_partition: SemanticPartition):
    """
    Print a human-readable summary of a semantic partition.
    
    Args:
        semantic_partition: A SemanticPartition object
    """
    print("=" * 70)
    print("SEMANTIC PARTITION SUMMARY")
    print("=" * 70)
    print(f"\nCategories ({len(semantic_partition.categories)}):")
    
    for category in semantic_partition.categories:
        species = semantic_partition.category_to_species[category]
        indices = semantic_partition.category_indices[category]
        print(f"\n  {category}:")
        print(f"    Species: {species}")
        print(f"    Indices: {indices}")
    
    # Show overlapping species if any
    overlapping = [
        (species, cats)
        for species, cats in semantic_partition.species_to_categories.items()
        if len(cats) > 1
    ]
    
    if overlapping:
        print(f"\n⚠️  Overlapping Species ({len(overlapping)}):")
        for species, categories in overlapping:
            print(f"    {species}: {categories}")
    
    print("=" * 70)


# Utility function for testing
def create_conflict_semantic_partition(species_list: List[str]) -> SemanticPartition:
    """
    Create a standard semantic partition for the conflict model.
    
    This is a convenience function that defines the standard categories
    for the cause-driven conflict model:
    - conflict: violent species and damage
    - peace: peaceful species and grievances
    - perception: internal perception 
    
    Args:
        species_list: List of all species names from the conflict model
        
    Returns:
        SemanticPartition configured for conflict analysis
    """
    category_definitions = {
        'conflict': [s for s in species_list if '_v' in s or s.startswith('D_')],
        'peace': [s for s in species_list if '_p' in s or s.startswith('G_')],
        'perception': [s for s in species_list if s.startswith('i')]
    }
    
    return define_semantic_categories(species_list, category_definitions)