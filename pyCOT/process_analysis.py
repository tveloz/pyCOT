"""
Process Analysis Library for Reaction Networks (Extended)

This module extends process_analysis.py with semantic-aware analysis and 
process decomposition capabilities.

Main components:
1. Process Classification - classify process vectors based on stoichiometric effects
2. Time Series Processing - rescaling and batch classification
3. Cone Analysis - nullspace, feasible regions, and comprehensive cone analysis
4. Utilities - helper functions for interval extraction and data export
5. Semantic Process Analysis - semantic-aware classification and decomposition (NEW)
"""

import numpy as np
import pandas as pd
import sympy as sp
from typing import List, Dict, Tuple, Optional, Union
from collections import Counter
import os


######################################################################################
# 1. PROCESS CLASSIFICATION
######################################################################################

def classify_process_mode(v, S, species_subset=None, species_list=None, tol=1e-3):
    """
    Clasifica un vector de proceso v con respecto a la matriz estequiométrica S.
    
    La clasificación puede hacerse sobre todas las especies o sobre un subconjunto específico.
    Esto permite analizar cómo un proceso afecta diferentes grupos semánticos de especies.

    Categorías base:
    - Stationary Mode: El proceso no produce cambios netos en las concentraciones (Sv=0).
    - Problem: El proceso consume especies pero no produce ninguna (Sv <= 0 y al menos un Sv < 0).
    - Challenge: El proceso consume al menos una especie (al menos un Sv < 0) y produce al menos una (Sv > 0).
    - Overproduction Mode: El proceso produce especies sin consumir (Sv >= 0 y al menos un Sv > 0).
    
    Args:
        v (np.ndarray): El vector de proceso a clasificar.
        S (np.ndarray): La matriz estequiométrica.
        species_subset (list, optional): Subconjunto de especies a considerar.
            Puede ser:
            - Lista de nombres de especies (strings): ['A_v', 'B_v', 'D_A']
            - Lista de índices (integers): [0, 1, 4]
            Si es None, considera todas las especies.
        species_list (list, optional): Lista completa de nombres de especies.
            Requerido si species_subset contiene nombres (strings).
        tol (float): Tolerancia para comparaciones numéricas.
        
    Returns:
        list[str]: Una lista de las categorías a las que pertenece el proceso v
                   con respecto al subconjunto de especies especificado.
                   
    Example:
        >>> # Clasificar sobre todas las especies
        >>> classify_process_mode(v, S)
        ['Challenge', 'Incomplete Process']
        
        >>> # Clasificar usando nombres de especies
        >>> classify_process_mode(v, S, 
        ...                      species_subset=['A_v', 'B_v', 'D_A'],
        ...                      species_list=species)
        ['Overproduction Mode', 'Incomplete Process']
        
        >>> # Clasificar usando índices (backend/performance)
        >>> classify_process_mode(v, S, species_subset=[0, 1, 4])
        ['Overproduction Mode', 'Incomplete Process']
    """
    
    if not isinstance(v, np.ndarray) or not isinstance(S, np.ndarray):
        raise TypeError("Los inputs 'v' y 'S' deben ser arrays de NumPy.")

    # Compute full stoichiometric effect
    Sv = S @ v
    
    # Convert species_subset to indices if needed
    if species_subset is not None:
        if len(species_subset) == 0:
            raise ValueError("species_subset no puede estar vacío")
        
        # Check if species_subset contains strings (species names)
        if isinstance(species_subset[0], str):
            if species_list is None:
                raise ValueError("species_list es requerido cuando species_subset contiene nombres de especies")
            
            # Convert names to indices
            species_indices = {name: idx for idx, name in enumerate(species_list)}
            
            # Validate that all species exist
            for name in species_subset:
                if name not in species_indices:
                    raise ValueError(f"Especie '{name}' no encontrada en species_list")
            
            species_subset_indices = [species_indices[name] for name in species_subset]
        else:
            # Already indices
            species_subset_indices = species_subset
        
        Sv_relevant = Sv[species_subset_indices]
    else:
        Sv_relevant = Sv
    
    classifications = []

    # Classify based on relevant species only
    is_stationary = np.all((-tol <= Sv_relevant) & (Sv_relevant <= tol))
    is_overproduced = np.all(Sv_relevant >= -tol) and np.any(Sv_relevant > tol)
    is_challenge = np.any(Sv_relevant < -tol) and np.any(Sv_relevant > tol)
    is_problem = np.all(Sv_relevant <= tol) and np.any(Sv_relevant <= -tol)
    is_complete = np.all(v > tol) 
    
    if is_stationary:
        classifications = ["Stationary Mode"]
    elif is_overproduced:
        classifications = ["Overproduction Mode"]
    elif is_challenge:
        classifications = ["Challenge"]
    elif is_problem:
        classifications = ["Problem"]
    else:
        classifications = ["Other"]
    
    if is_complete:
        classifications.append("Complete Process")
    else:
        classifications.append("Incomplete Process")
    
    return classifications


def is_cognitive_domain(v, S, species_subset=None, species_list=None, tol=1e-3):
    """
    Verifica si un proceso v está en el dominio cognitivo.
    
    Un proceso está en el dominio cognitivo si:
    - Es Stationary Mode o Overproduction Mode
    - Es un proceso completo (v > 0 para todas las reacciones)
    
    Args:
        v: Proceso a verificar
        S: Matriz estequiométrica
        species_subset: Opcional, subconjunto de especies (nombres o índices)
        species_list: Opcional, lista de nombres de especies (requerido si species_subset usa nombres)
        tol: Tolerancia
    """
    v_class = classify_process_mode(v, S, species_subset=species_subset, species_list=species_list, tol=tol) 
    if (v_class[0] == "Stationary Mode" or v_class[0] == "Overproduction Mode") and v_class[1] == "Complete Process":
        return True 
    else:
        return False


######################################################################################
# 5. SEMANTIC PROCESS ANALYSIS
######################################################################################

def decompose_process(v, S=None, tol=1e-3):
    """
    Decompose a process vector into single reaction components.
    
    Args:
        v: Process vector
        S: Stoichiometric matrix (optional, only used if Sv_component needed)
        tol: Tolerance for considering reaction active
        
    Returns:
        List of dicts, each containing:
            'reaction_idx': int
            'rate': float
            'Sv_component': array (if S provided)
    """
    if not isinstance(v, np.ndarray):
        raise TypeError("v must be numpy array")
    
    components = []
    active_reactions = np.where(np.abs(v) > tol)[0]
    
    for reaction_idx in active_reactions:
        rate = v[reaction_idx]
        component = {
            'reaction_idx': int(reaction_idx),
            'rate': float(rate)
        }
        
        if S is not None:
            component['Sv_component'] = S[:, reaction_idx] * rate
        
        components.append(component)
    
    return components


def compute_semantic_impact(Sv, semantic_partition, tol=1e-3):
    """
    Compute net effect on each semantic category.
    
    Args:
        Sv: Net stoichiometric effect (S @ v)
        semantic_partition: SemanticPartition object
        tol: Tolerance
        
    Returns:
        Dict mapping category to net effect
    """
    net_effects = {}
    
    for category in semantic_partition.categories:
        indices = semantic_partition.category_indices[category]
        net_effect = np.sum(Sv[indices])
        net_effects[category] = net_effect
    
    return net_effects


def classify_semantic_relation(Sv, semantic_partition, tol=1e-3):
    """
    Classify the semantic relation of a process.
    
    Args:
        Sv: Net stoichiometric effect
        semantic_partition: SemanticPartition object
        tol: Tolerance
        
    Returns:
        Dict with keys:
            'type': relation type
            'source': source category
            'target': target category
            'net_effects': dict of net effects per category
            'interpretation': human-readable string
    """
    net_effects = compute_semantic_impact(Sv, semantic_partition, tol)
    
    depleting_categories = [cat for cat, net in net_effects.items() if net < -tol]
    amplifying_categories = [cat for cat, net in net_effects.items() if net > tol]
    
    n_depleting = len(depleting_categories)
    n_amplifying = len(amplifying_categories)
    
    # Determine relation type
    if n_depleting == 0 and n_amplifying == 0:
        relation_type = 'PROCESS_BALANCES'
        source = None
        target = None
        interpretation = "Balanced process (steady state)"
    
    elif n_depleting == 0 and n_amplifying > 0:
        if n_amplifying == 1:
            relation_type = 'PROCESS_AMPLIFIES'
            source = amplifying_categories[0]
            target = amplifying_categories[0]
            interpretation = f"Amplifies {amplifying_categories[0]}"
        else:
            relation_type = 'PROCESS_AMPLIFIES_MULTIPLE'
            source = amplifying_categories
            target = amplifying_categories
            interpretation = f"Amplifies {amplifying_categories}"
    
    elif n_depleting > 0 and n_amplifying == 0:
        if n_depleting == 1:
            relation_type = 'PROCESS_DEPLETES'
            source = depleting_categories[0]
            target = None
            interpretation = f"Depletes {depleting_categories[0]}"
        else:
            relation_type = 'PROCESS_DEPLETES_MULTIPLE'
            source = depleting_categories
            target = None
            interpretation = f"Depletes {depleting_categories}"
    
    elif n_depleting == 1 and n_amplifying == 1:
        relation_type = 'PROCESS_GENERATES'
        source = depleting_categories[0]
        target = amplifying_categories[0]
        interpretation = f"Generates {amplifying_categories[0]} from {depleting_categories[0]}"
    
    else:
        relation_type = 'PROCESS_TRANSFORMS'
        source = depleting_categories
        target = amplifying_categories
        interpretation = f"Transforms {depleting_categories} into {amplifying_categories}"
    
    return {
        'type': relation_type,
        'source': source,
        'target': target,
        'net_effects': net_effects,
        'interpretation': interpretation
    }


def analyze_process_semantic(v, S, semantic_partition, species_subset=None, species_list=None, tol=1e-3):
    """
    Analyze a process with both mode and semantic classifications.
    
    Args:
        v: Process vector
        S: Stoichiometric matrix
        semantic_partition: SemanticPartition object
        species_subset: Optional, subset of species (names or indices) for mode classification
        species_list: Optional, list of species names (required if species_subset uses names)
        tol: Tolerance
        
    Returns:
        Dict with:
            'Sv': net stoichiometric effect
            'mode': process mode classification
            'relation': semantic relation classification
    """
    Sv = S @ v
    mode = classify_process_mode(v, S, species_subset=species_subset, species_list=species_list, tol=tol)
    relation = classify_semantic_relation(Sv, semantic_partition, tol)
    
    return {
        'Sv': Sv,
        'mode': mode[0],
        'relation': relation
    }


def random_decompose_process(v, n, strategy='mixed', seed=None):
    """
    Randomly decompose a process into n sub-processes.
    
    Supports two decomposition strategies:
    - 'coordinate': Assign each active reaction entirely to one part
    - 'value': Split each active reaction's value across multiple parts
    - 'mixed': Randomly choose strategy per reaction
    
    Args:
        v: Process vector
        n: Number of parts (must be >= 1)
        strategy: Decomposition strategy
        seed: Random seed for reproducibility
        
    Returns:
        List of n process vectors [v1, v2, ..., vn] where sum(vi) = v
    """
    if seed is not None:
        np.random.seed(seed)
    
    if not isinstance(v, np.ndarray):
        raise TypeError("v must be numpy array")
    
    active_reactions = np.where(v > 1e-10)[0]
    k = len(active_reactions)
    
    if n < 1:
        raise ValueError("n must be at least 1")
    
    if n > k:
        raise ValueError(f"n={n} cannot exceed number of active reactions k={k}")
    
    # Initialize n empty process vectors
    parts = [np.zeros_like(v) for _ in range(n)]
    
    for reaction_idx in active_reactions:
        rate = v[reaction_idx]
        
        # Decide strategy for this reaction
        if strategy == 'coordinate':
            # Assign entire reaction to one random part
            part_idx = np.random.randint(0, n)
            parts[part_idx][reaction_idx] = rate
            
        elif strategy == 'value':
            # Split value across all n parts using Dirichlet
            alphas = np.ones(n)
            fractions = np.random.dirichlet(alphas)
            for i in range(n):
                parts[i][reaction_idx] = rate * fractions[i]
                
        elif strategy == 'mixed':
            # Random choice: coordinate or value split
            if np.random.random() < 0.5:
                # Coordinate assignment
                part_idx = np.random.randint(0, n)
                parts[part_idx][reaction_idx] = rate
            else:
                # Value split
                alphas = np.ones(n)
                fractions = np.random.dirichlet(alphas)
                for i in range(n):
                    parts[i][reaction_idx] = rate * fractions[i]
        else:
            raise ValueError(f"Unknown strategy: {strategy}")
    
    return parts


def analyze_decomposition_statistics(decompositions, S, semantic_partition=None, species_subset=None, species_list=None, tol=1e-3):
    """
    Analyze statistical properties of process decompositions.
    
    Args:
        decompositions: List of decomposition results, each is list of process vectors
        S: Stoichiometric matrix
        semantic_partition: Optional SemanticPartition object
        species_subset: Optional list of species (names or indices) to focus classification on
        species_list: Optional list of species names (required if species_subset uses names)
        tol: Tolerance
        
    Returns:
        Dict with:
            'mode_distribution': Counter of process modes
            'relation_distribution': Counter of semantic relations (if semantic_partition provided)
            'all_parts': list of all part analyses
    """
    mode_counts = Counter()
    relation_counts = Counter()
    all_parts = []
    
    for decomposition in decompositions:
        for part in decomposition:
            # Skip zero vectors
            if np.all(np.abs(part) < tol):
                continue
            
            # Analyze mode (with optional species subset)
            mode = classify_process_mode(part, S, species_subset=species_subset, species_list=species_list, tol=tol)
            mode_counts[mode[0]] += 1
            
            part_analysis = {
                'v': part,
                'mode': mode[0]
            }
            
            # Analyze semantic relation if partition provided
            if semantic_partition is not None:
                Sv = S @ part
                relation = classify_semantic_relation(Sv, semantic_partition, tol)
                relation_counts[relation['type']] += 1
                part_analysis['relation'] = relation['type']
                part_analysis['relation_full'] = relation
            
            all_parts.append(part_analysis)
    
    result = {
        'mode_distribution': dict(mode_counts),
        'all_parts': all_parts,
        'n_parts': len(all_parts)
    }
    
    if semantic_partition is not None:
        result['relation_distribution'] = dict(relation_counts)
    
    return result


######################################################################################
# EXPORTS
######################################################################################

__all__ = [
    # Process Classification
    'classify_process_mode',
    'is_cognitive_domain',
    
    # Semantic Analysis
    'decompose_process',
    'compute_semantic_impact',
    'classify_semantic_relation',
    'analyze_process_semantic',
    
    # Decomposition
    'random_decompose_process',
    'analyze_decomposition_statistics'
]