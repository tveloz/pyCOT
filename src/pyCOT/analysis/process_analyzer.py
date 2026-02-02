"""
Process Analysis Module - Semantic Interpretation of Process Structures

This module interprets process structures (from process_structure.py) through the
lens of semantic categories (from semantic_partition.py). It answers: "What do 
structural operations mean for the categories we care about?"

Provides three analytical modes:
1. One-Process-Many-Decompositions: Pathway variability analysis
2. Many-Processes-Rolling-Time: Temporal scale analysis  
3. Many-Processes-Many-Decompositions-And-Rollings: Full multi-scale robustness

Part of the Challenge-Centered Semantic Pipeline (refactored architecture)
"""


import numpy as np
import pandas as pd
from typing import List, Dict, Tuple, Optional, Union, Callable
from dataclasses import dataclass, field
from collections import Counter
import warnings
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# Local imports

from pyCOT.core.semantic_partition import SemanticPartition
from pyCOT.analysis.process_structure import (
    ProcessAggregation, ProcessDisaggregation,
    disaggregation_random_sampling, rolling_window_aggregation,
    decompose_to_single_reactions
)


######################################################################################
# DATA STRUCTURES
######################################################################################

@dataclass
class CategoryBehavior:
    """
    Describes how a semantic category behaves under a process.
    
    Attributes:
        category: Category name
        net_effect: Net production (+) or consumption (-)
        total_production: Total amount produced
        total_consumption: Total amount consumed
        balance_state: 'overproduced', 'depleted', 'balanced'
        balance_quality: How close to zero (0=perfect balance, higher=more imbalanced)
        species_affected: Dict of species indices to their effects
    """
    category: str
    net_effect: float
    total_production: float
    total_consumption: float
    balance_state: str
    balance_quality: float
    species_affected: Dict[int, float] = field(default_factory=dict)
    
    def __repr__(self):
        return (f"CategoryBehavior({self.category}: net={self.net_effect:.3f}, "
                f"state={self.balance_state})")


@dataclass
class SemanticPattern:
    """
    Identified semantic pattern in process behavior.
    
    Attributes:
        pattern_type: Type of pattern
        categories_involved: List of category names
        strength: Magnitude of pattern
        description: Human-readable description
        metadata: Additional pattern-specific information
    """
    pattern_type: str
    categories_involved: List[str]
    strength: float
    description: str
    metadata: Dict = field(default_factory=dict)
    
    def __repr__(self):
        return f"SemanticPattern({self.pattern_type}: {self.description})"


@dataclass
class ThresholdDetection:
    """
    Detected threshold in category behavior across scales.
    
    Attributes:
        scale_value: Scale at which threshold occurs (window size or n_components)
        category: Category affected (or None for system-level)
        threshold_type: Type of threshold
        feature_before: Feature value before threshold
        feature_after: Feature value after threshold
        confidence: Statistical confidence (0-1)
        interpretation: Human-readable explanation
    """
    scale_value: float
    category: Optional[str]
    threshold_type: str
    feature_before: float
    feature_after: float
    confidence: float
    interpretation: str
    metadata: Dict = field(default_factory=dict)
    
    def __repr__(self):
        cat_str = self.category if self.category else "system"
        return (f"ThresholdDetection(scale={self.scale_value}, {cat_str}: "
                f"{self.threshold_type})")


######################################################################################
# PROCESS MODE CLASSIFICATION
######################################################################################

def classify_process_mode(
    v: np.ndarray,
    S: np.ndarray,
    species_subset: Optional[Union[List[str], List[int]]] = None,
    species_list: Optional[List[str]] = None,
    tol: float = 1e-3
) -> Tuple[str, str]:
    """
    Classify operational mode of a process (context-dependent).
    
    Process modes based on stoichiometric effect Sv:
    - Stationary Mode: Sv â‰ˆ 0 (balanced production/consumption)
    - Overproduction Mode: Sv â‰¥ 0 with some > 0 (net production)
    - Challenge Mode: Sv has mixed signs (some up, some down)
    - Problem Mode: Sv â‰¤ 0 with some < 0 (net consumption)
    
    Context-dependency: Classification is relative to species_subset if provided,
    allowing category-specific analysis (e.g., "overproduction for peace category").
    
    Args:
        v: Process vector
        S: Stoichiometric matrix (n_species x n_reactions)
        species_subset: Optional subset of species (names or indices) for classification
        species_list: Full list of species names (required if species_subset uses names)
        tol: Tolerance for numerical zero
        
    Returns:
        Tuple of (mode, completeness):
            mode: Process mode classification
            completeness: Whether process is complete (all reactions fire)
    """
    if not isinstance(v, np.ndarray):
        raise TypeError("v must be numpy array")
    
    if not isinstance(S, np.ndarray):
        raise TypeError("S must be numpy array")
    
    # Compute stoichiometric effect
    Sv = S @ v
    
    # Handle species subset (context-dependent classification)
    if species_subset is not None:
        if isinstance(species_subset[0], str):
            # Names provided - need species_list to convert to indices
            if species_list is None:
                raise ValueError("species_list required when species_subset contains names")
            indices = [species_list.index(name) for name in species_subset]
        else:
            # Indices provided directly
            indices = species_subset
        
        # Restrict to subset
        Sv = Sv[indices]
    
    # Classify based on sign pattern of Sv
    positive = Sv > tol
    negative = Sv < -tol
    zero = np.abs(Sv) <= tol
    
    n_positive = np.sum(positive)
    n_negative = np.sum(negative)
    n_zero = np.sum(zero)
    
    # Determine mode
    if n_zero == len(Sv):
        mode = "Stationary Mode"
    elif n_negative == 0 and n_positive > 0:
        mode = "Overproduction Mode"
    elif n_positive == 0 and n_negative > 0:
        mode = "Problem Mode"
    else:
        mode = "Challenge Mode"
    
    # Determine completeness
    n_active = np.sum(np.abs(v) > tol)
    if n_active == len(v):
        completeness = "Complete Process"
    else:
        completeness = "Incomplete Process"
    
    return (mode, completeness)


def is_cognitive_domain(
    v: np.ndarray,
    S: np.ndarray,
    species_subset: Optional[Union[List[str], List[int]]] = None,
    species_list: Optional[List[str]] = None,
    tol: float = 1e-3
) -> bool:
    """
    Check if process is in cognitive domain (self-maintaining).
    
    Cognitive domain = Stationary or Overproduction mode with complete process.
    These are processes that can maintain or grow organizations.
    
    Args:
        v: Process vector
        S: Stoichiometric matrix
        species_subset: Optional species subset for context
        species_list: Full species names list
        tol: Tolerance
        
    Returns:
        True if process is in cognitive domain
    """
    mode, completeness = classify_process_mode(v, S, species_subset, species_list, tol)
    
    is_maintaining = mode in ["Stationary Mode", "Overproduction Mode"]
    is_complete = completeness == "Complete Process"
    
    return is_maintaining and is_complete
def classify_response_to_disturbance(v, S, v_pert=None, x_pert=None, tol=1e-3, verbose=False):
    """
    Clasifica la respuesta de un proceso a una perturbaciÃ³n.
    
    Analiza cÃ³mo un proceso v responde a una perturbaciÃ³n, que puede ser:
    - v_pert: Una perturbaciÃ³n en el espacio de procesos
    - x_pert: Una perturbaciÃ³n en el espacio de estados
    
    Args:
        v (np.ndarray): El vector de proceso de respuesta.
        S (np.ndarray): La matriz estequiomÃ©trica.
        v_pert (np.ndarray, optional): PerturbaciÃ³n en el espacio de procesos.
        x_pert (np.ndarray, optional): PerturbaciÃ³n en el espacio de estados.
        tol (float): Tolerancia para comparaciones numÃ©ricas.
        verbose (bool): Si True, imprime informaciÃ³n de debug.
        
    Returns:
        list[str]: ClasificaciÃ³n de la respuesta a la perturbaciÃ³n.
        
    Raises:
        TypeError: Si los tipos de entrada no son correctos.
        ValueError: Si las dimensiones no coinciden.
    """
    if v_pert is not None and not isinstance(v_pert, np.ndarray):
        raise TypeError("El input 'v_pert', si se proporciona, debe ser un array de NumPy.")
    
    if v_pert is not None and v_pert.shape[0] != S.shape[1]:
        raise ValueError("El input 'v_pert' debe tener la misma dimensiÃ³n que el nÃºmero de reacciones en 'S'.")
    
    if x_pert is not None and not isinstance(x_pert, np.ndarray):
        raise TypeError("El input 'x_pert', si se proporciona, debe ser un array de NumPy.")
    
    if x_pert is not None and x_pert.shape[0] != S.shape[0]:
        raise ValueError("El input 'x_pert' debe tener la misma dimensiÃ³n que el nÃºmero de especies en 'S'.")

    v_mode = classify_process_mode(v, S, tol=tol)
    
    if v_pert is None and x_pert is None:
        if verbose:
            print("Null disturbance")
        return v_mode

    elif v_pert is not None and x_pert is None:
        if verbose:
            print("Process disturbance")
        v_pert_class = classify_process_mode(v_pert, S, tol=tol)
        v_cognitive_domain = is_cognitive_domain(v, S, tol=tol)
        v_combined = (v + v_pert)
        v_combined_class = classify_process_mode(v_combined, S, tol=tol)
        v_combined_cognitive_domain = is_cognitive_domain(v_combined, S, tol=tol)
        
        # Cognitive Control Situations
        if v_cognitive_domain:
            if verbose:
                print("In Cognitive Domain")
            # Cognitive Domain controla Challenge
            if v_pert_class[0] == "Challenge":
                if v_combined_cognitive_domain:
                    return ["Cognitive Controls Challenge"]
                else:
                    return ["Cognitive Breakdown by Challenge"]
            # Cognitive Domain controla Problem
            if v_pert_class[0] == "Problem": 
                if v_combined_cognitive_domain:
                    return ["Cognitive Controls Problem"]
                else:
                    return ["Cognitive Breakdown by Problem"]
            # Cognitive Domain se sostiene por proceso automantenido o sobreproducido
            if (v_pert_class[0] == "Stationary Mode" or v_pert_class[0] == "Overproduction Mode"): 
                if v_combined_cognitive_domain:
                    return ["Cognitive Domain Sustained"]
            else:
                if verbose:
                    print("Check Cognitive Breakdown by Glitch!!" + str(v_mode) + " + " + str(v_pert_class) + " =>" + str(v_combined_class))
                    return ["Cognitive Breakdown by Glitch"]
        else:
            if verbose:
                print("Out of Cognitive Domain")
                print(f"v_mode: {v_mode}, v_pert_class: {v_pert_class}, v_combined_class: {v_combined_class}")
                print(f"v_combined_cognitive_domain: {v_combined_cognitive_domain}")
            if v_combined_cognitive_domain:
                return ["Cognitive Domain Recovered"]
            else:
                return ["Cognitive Breakdown Sustained"]
                
    elif v_pert is None and x_pert is not None:
        # State Disturbance and response result
        v_cognitive_domain = is_cognitive_domain(v, S, tol=tol)
        x_next = x_pert + S @ v
        
        if v_cognitive_domain:
            if np.any(x_pert < -tol) and np.any(x_pert > tol):
                if verbose:
                    print("State disturbance is a challenge")
                if np.all(x_next >= -tol):
                    return ["Cognitive Controls Challenge"]
                else:
                    if np.any(x_next < -tol):
                        return ["Cognitive Breakdown by Challenge"]
            elif np.all(x_pert <= -tol) and np.any(x_pert < tol):
                if verbose:
                    print("State disturbance is a problem")
                if np.all(x_next >= -tol):
                    return ["Cognitive Controls Problem"]
                else:
                    if np.any(x_next < -tol):
                        return ["Cognitive Breakdown by Problem"]
            else:
                if verbose:
                    print("State disturbance is resources incoming")
                if np.all(x_next >= -tol):
                    return ["Cognitive Domain Sustained"]
                else:
                    return ["Cognitive Breakdown by Glitch"]
        else:
            if verbose:
                print("Out of Cognitive Domain")
                print(f"v_mode: {v_cognitive_domain} x_next: {x_next}")
            if np.all(x_next >= -tol):
                return ["Cognitive Domain Recovered"]
            else:
                return ["Cognitive Breakdown Sustained"]

######################################################################################
# CATEGORY BEHAVIOR ANALYSIS
######################################################################################

def analyze_category_behavior(
    process: np.ndarray,
    S: np.ndarray,
    semantic_partition: SemanticPartition,
    tol: float = 1e-3
) -> Dict[str, CategoryBehavior]:
    """
    Analyze how each semantic category behaves under a process.
    
    For each category, computes:
    - Net effect (production - consumption)
    - Total production and consumption
    - Balance state classification
    - Balance quality metric
    
    Args:
        process: Process vector
        S: Stoichiometric matrix
        semantic_partition: SemanticPartition object
        tol: Tolerance for balance detection
        
    Returns:
        Dictionary mapping category name to CategoryBehavior
        
    Example:
        >>> behaviors = analyze_category_behavior(v, S, semantic_partition)
        >>> print(behaviors['conflict'].balance_state)
        'overproduced'
    """
    Sv = S @ process
    
    category_behaviors = {}
    
    for category in semantic_partition.categories:
        indices = semantic_partition.category_indices[category]
        
        # Get effects on this category's species
        category_effects = Sv[indices]
        
        # Compute aggregates
        production = np.sum(category_effects[category_effects > tol])
        consumption = -np.sum(category_effects[category_effects < -tol])
        net = production - consumption
        
        # Classify balance state
        if abs(net) <= tol:
            balance_state = 'balanced'
        elif net > tol:
            balance_state = 'overproduced'
        else:
            balance_state = 'depleted'
        
        # Balance quality: distance from perfect balance
        balance_quality = abs(net)
        
        # Species-level details
        species_affected = {idx: float(Sv[idx]) for idx in indices if abs(Sv[idx]) > tol}
        
        behavior = CategoryBehavior(
            category=category,
            net_effect=float(net),
            total_production=float(production),
            total_consumption=float(consumption),
            balance_state=balance_state,
            balance_quality=float(balance_quality),
            species_affected=species_affected
        )
        
        category_behaviors[category] = behavior
    
    return category_behaviors


def track_category_sequence(
    process_sequence: List[np.ndarray],
    S: np.ndarray,
    semantic_partition: SemanticPartition,
    tol: float = 1e-3
) -> Dict[str, List[CategoryBehavior]]:
    """
    Track category behavior across a sequence of processes.
    
    This analyzes how categories evolve through time or through decomposition steps.
    
    Args:
        process_sequence: List of process vectors (temporal or decomposition)
        S: Stoichiometric matrix
        semantic_partition: SemanticPartition object
        tol: Tolerance
        
    Returns:
        Dictionary mapping category -> list of behaviors over sequence
    """
    # Initialize storage
    category_sequences = {cat: [] for cat in semantic_partition.categories}
    
    for process in process_sequence:
        behaviors = analyze_category_behavior(process, S, semantic_partition, tol)
        
        for category, behavior in behaviors.items():
            category_sequences[category].append(behavior)
    
    return category_sequences


def compute_category_statistics(
    category_sequence: List[CategoryBehavior]
) -> Dict:
    """
    Compute statistics for a category across a sequence.
    
    Args:
        category_sequence: List of CategoryBehavior objects for same category
        
    Returns:
        Dictionary with statistical measures
    """
    if not category_sequence:
        return {}
    
    net_effects = [b.net_effect for b in category_sequence]
    balance_qualities = [b.balance_quality for b in category_sequence]
    balance_states = [b.balance_state for b in category_sequence]
    
    stats = {
        'n_steps': len(category_sequence),
        'net_effect': {
            'mean': float(np.mean(net_effects)),
            'std': float(np.std(net_effects)),
            'min': float(np.min(net_effects)),
            'max': float(np.max(net_effects)),
            'final': net_effects[-1] if net_effects else None
        },
        'balance_quality': {
            'mean': float(np.mean(balance_qualities)),
            'std': float(np.std(balance_qualities)),
            'min': float(np.min(balance_qualities)),
            'max': float(np.max(balance_qualities))
        },
        'balance_state_distribution': dict(Counter(balance_states)),
        'transitions': _count_state_transitions(balance_states)
    }
    
    return stats


def _count_state_transitions(states: List[str]) -> Dict:
    """Count transitions between balance states."""
    transitions = Counter()
    for i in range(len(states) - 1):
        transition = f"{states[i]}->{states[i+1]}"
        transitions[transition] += 1
    return dict(transitions)


######################################################################################
# PATTERN RECOGNITION
######################################################################################

def identify_semantic_patterns(
    category_behaviors: Dict[str, CategoryBehavior],
    tol: float = 1e-3
) -> List[SemanticPattern]:
    """
    Identify semantic patterns from category behaviors.
    
    Patterns detected:
    - Balanced Transformation: One category depletes, another grows proportionally
    - Reinforcing Cycle: Category overproduced with self-amplification
    - Cascade: Multiple categories changing together
    - Stabilizing Balance: Multiple categories near steady state
    
    Args:
        category_behaviors: Dictionary of category behaviors
        tol: Tolerance for pattern detection
        
    Returns:
        List of identified SemanticPattern objects
    """
    patterns = []
    
    depleting = [cat for cat, b in category_behaviors.items() 
                 if b.balance_state == 'depleted']
    amplifying = [cat for cat, b in category_behaviors.items() 
                  if b.balance_state == 'overproduced']
    balanced = [cat for cat, b in category_behaviors.items() 
                if b.balance_state == 'balanced']
    
    # Pattern: Balanced Transformation
    if len(depleting) == 1 and len(amplifying) == 1:
        depl_cat = depleting[0]
        ampl_cat = amplifying[0]
        depl_mag = abs(category_behaviors[depl_cat].net_effect)
        ampl_mag = abs(category_behaviors[ampl_cat].net_effect)
        
        # Check if roughly balanced
        if abs(depl_mag - ampl_mag) / max(depl_mag, ampl_mag) < 0.3:
            pattern = SemanticPattern(
                pattern_type='BALANCED_TRANSFORMATION',
                categories_involved=[depl_cat, ampl_cat],
                strength=float((depl_mag + ampl_mag) / 2),
                description=f"Balanced transformation from {depl_cat} to {ampl_cat}",
                metadata={
                    'source': depl_cat,
                    'target': ampl_cat,
                    'balance_ratio': float(min(depl_mag, ampl_mag) / max(depl_mag, ampl_mag))
                }
            )
            patterns.append(pattern)
    
    # Pattern: Reinforcing Cycle
    for cat in amplifying:
        behavior = category_behaviors[cat]
        # Simple heuristic: strong overproduction suggests self-reinforcement
        if behavior.net_effect > 2 * tol:
            pattern = SemanticPattern(
                pattern_type='REINFORCING_CYCLE',
                categories_involved=[cat],
                strength=float(behavior.net_effect),
                description=f"Reinforcing cycle amplifying {cat}",
                metadata={'category': cat, 'net_effect': behavior.net_effect}
            )
            patterns.append(pattern)
    
    # Pattern: Cascade
    if len(depleting) > 1:
        total_depletion = sum(abs(category_behaviors[cat].net_effect) 
                             for cat in depleting)
        pattern = SemanticPattern(
            pattern_type='CASCADE_COLLAPSE',
            categories_involved=depleting,
            strength=float(total_depletion),
            description=f"Cascade affecting {', '.join(depleting)}",
            metadata={'affected_categories': depleting}
        )
        patterns.append(pattern)
    
    # Pattern: Stabilizing Balance
    if len(balanced) >= 2:
        avg_balance_quality = np.mean([category_behaviors[cat].balance_quality 
                                      for cat in balanced])
        pattern = SemanticPattern(
            pattern_type='STABILIZING_BALANCE',
            categories_involved=balanced,
            strength=float(1.0 / (1.0 + avg_balance_quality)),  # Inverse quality
            description=f"Stabilizing balance across {', '.join(balanced)}",
            metadata={'balanced_categories': balanced}
        )
        patterns.append(pattern)
    
    return patterns


######################################################################################
# MODE 1: ONE-PROCESS-MANY-DECOMPOSITIONS
######################################################################################

def analyze_pathway_variability(
    process: np.ndarray,
    S: np.ndarray,
    semantic_partition: SemanticPartition,
    n_samples: int = 100,
    n_components: int = 10,
    strategy: str = 'mixed',
    seed: Optional[int] = None,
    tol: float = 1e-3
) -> Dict:
    """
    Mode 1: Analyze pathway variability for single process.
    
    Takes a single aggregate process and generates multiple decomposition
    hypotheses, analyzing how semantic categories behave across different
    construction pathways.
    
    Args:
        process: Process vector to analyze
        S: Stoichiometric matrix
        semantic_partition: SemanticPartition object
        n_samples: Number of decomposition trials (N1)
        n_components: Number of components per decomposition (N2)
        strategy: Disaggregation strategy
        seed: Random seed
        tol: Tolerance
        
    Returns:
        Dictionary with:
            'decompositions': List of ProcessDisaggregation objects
            'category_sequences': Category behavior across samples
            'category_statistics': Statistical summary per category
            'patterns_by_sample': Patterns identified in each sample
            'pattern_robustness': Pattern frequency across samples
    """
    # Generate decompositions
    disaggregations = disaggregation_random_sampling(
        process, n_samples, n_components, strategy, seed
    )
    
    # Analyze each decomposition
    all_category_sequences = []
    all_patterns = []
    
    for disagg in disaggregations:
        # Track categories through this decomposition
        cat_seq = track_category_sequence(
            disagg.components, S, semantic_partition, tol
        )
        all_category_sequences.append(cat_seq)
        
        # Identify patterns in final state
        final_behaviors = analyze_category_behavior(
            disagg.original_process, S, semantic_partition, tol
        )
        patterns = identify_semantic_patterns(final_behaviors, tol)
        all_patterns.append(patterns)
    
    # Compute statistics across samples
    category_statistics = {}
    for category in semantic_partition.categories:
        # Collect all behavior sequences for this category
        category_data = [seq[category] for seq in all_category_sequences]
        
        # Aggregate statistics
        all_nets = []
        all_balance_qualities = []
        for behaviors in category_data:
            nets = [b.net_effect for b in behaviors]
            quals = [b.balance_quality for b in behaviors]
            all_nets.extend(nets)
            all_balance_qualities.extend(quals)
        
        category_statistics[category] = {
            'net_effect': {
                'mean': float(np.mean(all_nets)),
                'std': float(np.std(all_nets)),
                'min': float(np.min(all_nets)),
                'max': float(np.max(all_nets))
            },
            'balance_quality': {
                'mean': float(np.mean(all_balance_qualities)),
                'std': float(np.std(all_balance_qualities))
            },
            'pathway_variance': float(np.var(all_nets))  # Key metric!
        }
    
    # Pattern robustness
    pattern_counts = Counter()
    for pattern_list in all_patterns:
        for pattern in pattern_list:
            pattern_counts[pattern.pattern_type] += 1
    
    pattern_robustness = {
        ptype: count / n_samples 
        for ptype, count in pattern_counts.items()
    }
    
    return {
        'decompositions': disaggregations,
        'category_sequences': all_category_sequences,
        'category_statistics': category_statistics,
        'patterns_by_sample': all_patterns,
        'pattern_robustness': pattern_robustness,
        'n_samples': n_samples,
        'n_components': n_components
    }


######################################################################################
# MODE 2: MANY-PROCESSES-ROLLING-TIME
######################################################################################

def analyze_temporal_scale(
    process_series: Union[List[np.ndarray], np.ndarray],
    S: np.ndarray,
    semantic_partition: SemanticPartition,
    window_sizes: List[int],
    stride: int = 1,
    tol: float = 1e-3
) -> Dict:
    """
    Mode 2: Analyze semantic behavior across temporal scales.
    
    Applies rolling window aggregation at multiple scales and analyzes
    how categories behave at each temporal resolution.
    
    Args:
        process_series: Sequence of process vectors
        S: Stoichiometric matrix
        semantic_partition: SemanticPartition object
        window_sizes: List of aggregation windows to analyze
        stride: Step size for rolling windows
        tol: Tolerance
        
    Returns:
        Dictionary with:
            'scales': Window sizes analyzed
            'aggregations_by_scale': ProcessAggregation objects per scale
            'category_behaviors_by_scale': Category behaviors at each scale
            'category_statistics_by_scale': Statistical summaries per scale
            'patterns_by_scale': Patterns identified at each scale
    """
    from pyCOT.analysis.process_structure import multi_scale_aggregation
    
    # Perform multi-scale aggregation
    aggregations_by_scale = multi_scale_aggregation(
        process_series, window_sizes, stride=stride
    )
    
    results_by_scale = {}
    
    for window_size in window_sizes:
        aggs = aggregations_by_scale.get(window_size, [])
        
        if not aggs:
            continue
        
        # Analyze each aggregated process at this scale
        all_behaviors = []
        all_patterns = []
        
        for agg in aggs:
            behaviors = analyze_category_behavior(
                agg.aggregated_process, S, semantic_partition, tol
            )
            all_behaviors.append(behaviors)
            
            patterns = identify_semantic_patterns(behaviors, tol)
            all_patterns.append(patterns)
        
        # Compute statistics across windows at this scale
        category_stats = {}
        for category in semantic_partition.categories:
            cat_behaviors = [b[category] for b in all_behaviors]
            category_stats[category] = compute_category_statistics(cat_behaviors)
        
        # Pattern frequency at this scale
        pattern_counts = Counter()
        for pattern_list in all_patterns:
            for pattern in pattern_list:
                pattern_counts[pattern.pattern_type] += 1
        
        pattern_freq = {
            ptype: count / len(aggs) if len(aggs) > 0 else 0
            for ptype, count in pattern_counts.items()
        }
        
        results_by_scale[window_size] = {
            'aggregations': aggs,
            'category_behaviors': all_behaviors,
            'category_statistics': category_stats,
            'patterns': all_patterns,
            'pattern_frequency': pattern_freq,
            'n_windows': len(aggs)
        }
    
    return {
        'scales': window_sizes,
        'results_by_scale': results_by_scale
    }


######################################################################################
# MODE 3: MANY-PROCESSES-MANY-DECOMPOSITIONS-AND-ROLLINGS
######################################################################################

def analyze_full_multiscale_robustness(
    process_series: Union[List[np.ndarray], np.ndarray],
    S: np.ndarray,
    semantic_partition: SemanticPartition,
    window_sizes: List[int],
    n_decomp_samples: int = 50,
    n_decomp_components: int = 5,
    stride: int = 1,
    strategy: str = 'mixed',
    seed: Optional[int] = None,
    tol: float = 1e-3
) -> Dict:
    """
    Mode 3: Full multi-scale robustness analysis.
    
    Combines temporal scale (Mode 2) with pathway variability (Mode 1).
    For each temporal scale, analyzes whether observed patterns are robust
    to different construction pathways.
    
    Args:
        process_series: Sequence of process vectors
        S: Stoichiometric matrix
        semantic_partition: SemanticPartition object
        window_sizes: Temporal scales to analyze
        n_decomp_samples: Number of decomposition samples per aggregated process
        n_decomp_components: Components per decomposition
        stride: Rolling window stride
        strategy: Disaggregation strategy
        seed: Random seed
        tol: Tolerance
        
    Returns:
        Dictionary with comprehensive multi-scale robustness analysis
    """
    from pyCOT.analysis.process_structure import multi_scale_aggregation
    
    # First, perform temporal aggregation at all scales
    aggregations_by_scale = multi_scale_aggregation(
        process_series, window_sizes, stride=stride
    )
    
    results_by_scale = {}
    
    for window_size in window_sizes:
        aggs = aggregations_by_scale.get(window_size, [])
        
        if not aggs:
            continue
        
        scale_results = []
        
        # For each aggregated process at this scale
        for i, agg in enumerate(aggs):
            # Analyze temporal pathway (actual sequence)
            temporal_behaviors = track_category_sequence(
                agg.component_processes, S, semantic_partition, tol
            )
            
            # Analyze decomposition pathways (sampled alternatives)
            decomp_analysis = analyze_pathway_variability(
                agg.aggregated_process, S, semantic_partition,
                n_samples=n_decomp_samples,
                n_components=n_decomp_components,
                strategy=strategy,
                seed=None if seed is None else seed + i,
                tol=tol
            )
            
            # Compare temporal vs. decomposition statistics
            comparison = _compare_temporal_vs_decomposition(
                temporal_behaviors,
                decomp_analysis['category_statistics'],
                semantic_partition
            )
            
            scale_results.append({
                'window_id': i,
                'aggregation': agg,
                'temporal_behaviors': temporal_behaviors,
                'decomposition_analysis': decomp_analysis,
                'temporal_vs_decomposition': comparison
            })
        
        # Aggregate statistics across all windows at this scale
        scale_summary = _summarize_scale_results(scale_results, semantic_partition)
        
        results_by_scale[window_size] = {
            'window_results': scale_results,
            'scale_summary': scale_summary,
            'n_windows': len(aggs)
        }
    
    return {
        'scales': window_sizes,
        'results_by_scale': results_by_scale,
        'n_decomp_samples': n_decomp_samples,
        'n_decomp_components': n_decomp_components
    }


def _compare_temporal_vs_decomposition(
    temporal_seq: Dict,
    decomp_stats: Dict,
    semantic_partition: SemanticPartition
) -> Dict:
    """Compare temporal sequence statistics with decomposition statistics."""
    comparison = {}
    
    for category in semantic_partition.categories:
        temp_behaviors = temporal_seq[category]
        temp_nets = [b.net_effect for b in temp_behaviors]
        
        decomp_net_mean = decomp_stats[category]['net_effect']['mean']
        decomp_net_std = decomp_stats[category]['net_effect']['std']
        
        temp_mean = np.mean(temp_nets)
        temp_std = np.std(temp_nets)
        
        # Z-score: how many standard deviations is temporal from decomposition mean?
        if decomp_net_std > 0:
            z_score = (temp_mean - decomp_net_mean) / decomp_net_std
        else:
            z_score = 0.0
        
        comparison[category] = {
            'temporal_mean': float(temp_mean),
            'decomposition_mean': float(decomp_net_mean),
            'difference': float(temp_mean - decomp_net_mean),
            'z_score': float(z_score),
            'temporal_within_decomp_range': bool(abs(z_score) < 2.0)
        }
    
    return comparison


def _summarize_scale_results(scale_results: List[Dict], semantic_partition: SemanticPartition) -> Dict:
    """Summarize results across all windows at a scale."""
    if not scale_results:
        return {}
    
    summary = {
        'n_windows': len(scale_results)
    }
    
    # Aggregate comparisons
    for category in semantic_partition.categories:
        z_scores = [
            result['temporal_vs_decomposition'][category]['z_score']
            for result in scale_results
        ]
        
        summary[category] = {
            'mean_z_score': float(np.mean(z_scores)),
            'std_z_score': float(np.std(z_scores)),
            'robustness': float(np.mean([abs(z) < 2.0 for z in z_scores]))
        }
    
    return summary


######################################################################################
# THRESHOLD DETECTION
######################################################################################

def detect_thresholds(
    scale_results: Dict,
    feature_name: str = 'net_effect',
    method: str = 'derivative',
    threshold_sensitivity: float = 0.5,
    min_confidence: float = 0.6
) -> List[ThresholdDetection]:
    """
    Detect thresholds in category behavior across scales.
    
    Identifies scales where category properties change qualitatively.
    
    Args:
        scale_results: Results from Mode 2 or Mode 3 analysis
        feature_name: Category feature to analyze ('net_effect', 'balance_quality', etc.)
        method: Detection method ('derivative', 'sign_change', 'variance_spike')
        threshold_sensitivity: Sensitivity parameter (method-dependent)
        min_confidence: Minimum confidence for reporting threshold
        
    Returns:
        List of ThresholdDetection objects
    """
    thresholds = []
    
    if 'results_by_scale' not in scale_results:
        return thresholds
    
    scales = sorted(scale_results['results_by_scale'].keys())
    
    if len(scales) < 3:
        warnings.warn("Need at least 3 scales for threshold detection")
        return thresholds
    
    # Extract feature values across scales for each category
    # This is simplified - full implementation would be more sophisticated
    for category in scale_results['results_by_scale'][scales[0]]['category_statistics'].keys():
        feature_values = []
        
        for scale in scales:
            scale_data = scale_results['results_by_scale'][scale]
            if 'category_statistics' in scale_data:
                cat_stat = scale_data['category_statistics'].get(category, {})
                feature_dict = cat_stat.get(feature_name, {})
                value = feature_dict.get('mean', 0.0)
                feature_values.append(value)
            else:
                feature_values.append(0.0)
        
        # Apply detection method
        if method == 'derivative':
            detected = _detect_threshold_derivative(
                scales, feature_values, threshold_sensitivity
            )
        elif method == 'sign_change':
            detected = _detect_threshold_sign_change(
                scales, feature_values
            )
        elif method == 'variance_spike':
            # Would need variance data from Mode 3
            detected = []
        else:
            detected = []
        
        # Convert to ThresholdDetection objects
        for scale_idx, conf in detected:
            if conf >= min_confidence:
                threshold = ThresholdDetection(
                    scale_value=scales[scale_idx],
                    category=category,
                    threshold_type=f"{feature_name}_{method}",
                    feature_before=feature_values[scale_idx - 1] if scale_idx > 0 else 0.0,
                    feature_after=feature_values[scale_idx],
                    confidence=conf,
                    interpretation=f"{category} {feature_name} threshold at scale {scales[scale_idx]}",
                    metadata={'method': method, 'scales': scales, 'values': feature_values}
                )
                thresholds.append(threshold)
    
    return thresholds


def _detect_threshold_derivative(
    scales: List[float],
    values: List[float],
    sensitivity: float
) -> List[Tuple[int, float]]:
    """Detect thresholds using derivative method."""
    if len(values) < 2:
        return []
    
    # Compute discrete derivatives
    derivatives = []
    for i in range(1, len(values)):
        deriv = (values[i] - values[i-1]) / (scales[i] - scales[i-1])
        derivatives.append(abs(deriv))
    
    # Find peaks in derivative
    if not derivatives:
        return []
    
    threshold_value = sensitivity * np.mean(derivatives)
    
    detected = []
    for i, deriv in enumerate(derivatives):
        if deriv > threshold_value:
            confidence = min(1.0, deriv / (2 * threshold_value))
            detected.append((i + 1, confidence))  # i+1 because derivatives are offset
    
    return detected


def _detect_threshold_sign_change(
    scales: List[float],
    values: List[float]
) -> List[Tuple[int, float]]:
    """Detect thresholds using sign change method."""
    detected = []
    
    for i in range(1, len(values)):
        if np.sign(values[i]) != np.sign(values[i-1]) and values[i-1] != 0:
            # Sign changed
            confidence = 1.0  # High confidence for sign changes
            detected.append((i, confidence))
    
    return detected


######################################################################################
# UTILITY FUNCTIONS FOR VISUALIZATION
######################################################################################

def get_intervals_by_category(
    times: np.ndarray,
    process_types: List[str],
    target_category: str
) -> List[Tuple[float, float]]:
    """
    Extract time intervals where process type matches target category.
    
    Args:
        times: Array of time values
        process_types: List of process type strings for each time point
        target_category: Category to find intervals for
        
    Returns:
        List of (start_time, end_time) tuples
    """
    intervals = []
    in_interval = False
    start_time = None
    
    for i, ptype in enumerate(process_types):
        if target_category in ptype:
            if not in_interval:
                # Start new interval
                start_time = times[i]
                in_interval = True
        else:
            if in_interval:
                # End current interval
                intervals.append((start_time, times[i-1] if i > 0 else start_time))
                in_interval = False
    
    # Close final interval if still open
    if in_interval and start_time is not None:
        intervals.append((start_time, times[-1]))
    
    return intervals


def analyze_cone(
    S: np.ndarray,
    grid_max: Optional[float] = None,
    grid_res: int = 20,
    classify: bool = True,
    tol: float = 1e-6
) -> Dict:
    """
    Analyze the cone of feasible process vectors.
    
    Samples the process cone to find feasible vectors, computes nullspace,
    and optionally classifies each point.
    
    Args:
        S: Stoichiometric matrix (n_species x n_reactions)
        grid_max: Maximum value for grid sampling (auto-computed if None)
        grid_res: Resolution of grid sampling
        classify: Whether to classify each point
        tol: Tolerance for numerical operations
        
    Returns:
        Dictionary with:
            - feasible_points: List of feasible process vectors
            - nullspace_vectors: Basis for nullspace
            - classifications: List of classifications (if classify=True)
            - grid_max: Actual grid_max used
    """
    n_reactions = S.shape[1]
    
    # Auto-compute grid_max if not provided
    if grid_max is None:
        grid_max = 2.0  # Default reasonable value
    
    # Sample process vectors in positive orthant
    # Simple uniform sampling
    step = grid_max / grid_res
    grid_1d = np.arange(0, grid_max + step, step)
    
    feasible_points = []
    classifications = []
    
    # For computational efficiency, sample a subset rather than full grid
    # Use random sampling in the positive orthant
    n_samples = min(1000, grid_res ** min(n_reactions, 3))
    
    np.random.seed(42)  # For reproducibility
    for _ in range(n_samples):
        v = np.random.uniform(0, grid_max, n_reactions)
        
        # Check if non-negative (feasible in positive orthant)
        if np.all(v >= -tol):
            feasible_points.append(v)
            
            if classify:
                mode, _ = classify_process_mode(v, S, tol=tol)
                classifications.append(mode)
    
    # Compute nullspace of S
    try:
        _, _, Vt = np.linalg.svd(S)
        rank = np.linalg.matrix_rank(S, tol=tol)
        nullspace_vectors = Vt[rank:, :].T
    except:
        nullspace_vectors = np.array([])
    
    return {
        'feasible_points': feasible_points,
        'nullspace_vectors': nullspace_vectors,
        'classifications': classifications if classify else [],
        'grid_max': grid_max,
        'n_samples': len(feasible_points)
    }


def rescale_process_time_series(
    df: pd.DataFrame,
    window_size: int = 1
) -> pd.DataFrame:
    """
    Rescale/resample process time series by averaging over windows.
    
    Args:
        df: DataFrame with Time column and process data
        window_size: Size of rolling window for averaging
        
    Returns:
        Resampled DataFrame
    """
    if window_size <= 1:
        return df
    
    # Create windowed averages
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    numeric_cols = [col for col in numeric_cols if col != 'Time']
    
    if len(numeric_cols) == 0:
        return df
    
    # Compute rolling means
    result_data = []
    for i in range(0, len(df), window_size):
        window_data = df.iloc[i:i+window_size]
        
        row_data = {'Time': window_data['Time'].iloc[0]}
        for col in numeric_cols:
            row_data[col] = window_data[col].mean()
        
        result_data.append(row_data)
    
    return pd.DataFrame(result_data)


######################################################################################
# EXPORTS
######################################################################################

__all__ = [
    # Data structures
    'CategoryBehavior',
    'SemanticPattern',
    'ThresholdDetection',
    
    # Process classification
    'classify_process_mode',
    'is_cognitive_domain',
    
    # Category analysis
    'analyze_category_behavior',
    'track_category_sequence',
    'compute_category_statistics',
    
    # Pattern recognition
    'identify_semantic_patterns',
    
    # Analytical modes
    'analyze_pathway_variability',          # Mode 1
    'analyze_temporal_scale',               # Mode 2
    'analyze_full_multiscale_robustness',   # Mode 3
    
    # Threshold detection
    'detect_thresholds',
    
    # Utility functions for visualization
    'get_intervals_by_category',
    'analyze_cone',
    'rescale_process_time_series'
]