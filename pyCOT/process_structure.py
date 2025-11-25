"""
Process Structure Module - Mathematical Operations on Process Vectors

This module handles pure mathematical operations on process vectors without semantic
interpretation. It manages aggregation (summing processes) and disaggregation 
(decomposing processes into components), treating processes as mathematical objects 
in reaction space.

Core insight: Time evolution is process addition, so:
- Aggregation = simulating forward in time
- Disaggregation = reverse-engineering temporal construction

Part of the Challenge-Centered Semantic Pipeline (refactored architecture)
"""

import numpy as np
from typing import List, Dict, Tuple, Optional, Union
from dataclasses import dataclass, field
from collections import Counter


@dataclass
class ProcessAggregation:
    """
    Result of aggregating multiple processes over a time window.
    
    Attributes:
        aggregated_process: The summed process vector
        component_processes: List of original processes that were aggregated
        window_indices: Indices of processes in original sequence
        window_size: Number of processes aggregated
        metadata: Additional information (timestamps, etc.)
    """
    aggregated_process: np.ndarray
    component_processes: List[np.ndarray]
    window_indices: List[int]
    window_size: int
    metadata: Dict = field(default_factory=dict)
    
    def __repr__(self):
        return (f"ProcessAggregation(window_size={self.window_size}, "
                f"n_components={len(self.component_processes)}, "
                f"indices={self.window_indices[0]}:{self.window_indices[-1]})")
    def verify_sum(self, tol=1e-10):
        """Verify that components sum to original process."""
        reconstructed = sum(self.component_processes)
        difference = np.linalg.norm(self.aggregated_process - reconstructed)
        return difference < tol

@dataclass
class ProcessDisaggregation:
    """
    Result of disaggregating a process into components.
    
    Attributes:
        original_process: The process that was disaggregated
        components: List of component processes that sum to original
        n_components: Number of components
        strategy: Disaggregation strategy used
        metadata: Additional information
    """
    original_process: np.ndarray
    components: List[np.ndarray]
    n_components: int
    strategy: str
    metadata: Dict = field(default_factory=dict)
    
    def verify_sum(self, tol=1e-10):
        """Verify that components sum to original process."""
        reconstructed = sum(self.components)
        difference = np.linalg.norm(reconstructed - self.original_process)
        return difference < tol
    
    def __repr__(self):
        verified = "✓" if self.verify_sum() else "✗"
        return (f"ProcessDisaggregation(n_components={self.n_components}, "
                f"strategy='{self.strategy}', verified={verified})")


######################################################################################
# AGGREGATION FUNCTIONS
######################################################################################

def aggregate_processes(
    process_list: List[np.ndarray],
    indices: Optional[List[int]] = None,
    metadata: Optional[Dict] = None
) -> ProcessAggregation:
    """
    Aggregate (sum) a list of process vectors.
    
    This is the fundamental aggregation operation. Time evolution can be viewed
    as sequential aggregation: x[t+n] = x[t] + S·(v[t] + v[t+1] + ... + v[t+n-1])
    
    Args:
        process_list: List of process vectors to aggregate
        indices: Optional list of indices in original sequence
        metadata: Optional metadata dictionary
        
    Returns:
        ProcessAggregation object
        
    Example:
        >>> v1 = np.array([1.0, 2.0, 0.5])
        >>> v2 = np.array([0.5, 1.0, 1.5])
        >>> v3 = np.array([2.0, 0.5, 1.0])
        >>> agg = aggregate_processes([v1, v2, v3])
        >>> print(agg.aggregated_process)
        [3.5 3.5 3.0]
    """
    if not process_list:
        raise ValueError("process_list cannot be empty")
    
    # Verify all processes have same shape
    shapes = [v.shape for v in process_list]
    if len(set(shapes)) > 1:
        raise ValueError(f"All processes must have same shape, got {shapes}")
    
    aggregated = sum(process_list)
    
    if indices is None:
        indices = list(range(len(process_list)))
    
    if metadata is None:
        metadata = {}
    
    return ProcessAggregation(
        aggregated_process=aggregated,
        component_processes=process_list.copy(),
        window_indices=indices,
        window_size=len(process_list),
        metadata=metadata
    )


def rolling_window_aggregation(
    process_series: Union[List[np.ndarray], np.ndarray],
    window_size: int,
    stride: int = 1,
    include_metadata: bool = True
) -> List[ProcessAggregation]:
    """
    Apply rolling window aggregation to a sequence of processes.
    
    This creates overlapping (or non-overlapping if stride=window_size) temporal
    windows, aggregating processes within each window. Equivalent to changing
    temporal resolution - like changing video frame rate.
    
    Args:
        process_series: Sequence of process vectors (list or 2D array)
        window_size: Number of processes to aggregate in each window
        stride: Step size between windows (default 1 = maximum overlap)
        include_metadata: Whether to include window timing information
        
    Returns:
        List of ProcessAggregation objects, one per window
        
    Example:
        >>> # Create sequence of 10 processes
        >>> series = [np.random.rand(5) for _ in range(10)]
        >>> # Aggregate with window=3, stride=1 (overlapping)
        >>> aggregations = rolling_window_aggregation(series, window_size=3, stride=1)
        >>> len(aggregations)  # Should be 8 windows
        8
    """
    # Convert to list of arrays if needed
    if isinstance(process_series, np.ndarray):
        if process_series.ndim == 1:
            raise ValueError("process_series must be 2D array (n_steps x n_reactions)")
        process_series = [process_series[i] for i in range(len(process_series))]
    
    n_processes = len(process_series)
    
    if window_size < 1:
        raise ValueError(f"window_size must be >= 1, got {window_size}")
    
    if window_size > n_processes:
        raise ValueError(f"window_size ({window_size}) cannot exceed number of processes ({n_processes})")
    
    if stride < 1:
        raise ValueError(f"stride must be >= 1, got {stride}")
    
    aggregations = []
    
    # Compute number of windows
    n_windows = (n_processes - window_size) // stride + 1
    
    for i in range(n_windows):
        start_idx = i * stride
        end_idx = start_idx + window_size
        
        window_processes = process_series[start_idx:end_idx]
        window_indices = list(range(start_idx, end_idx))
        
        metadata = {}
        if include_metadata:
            metadata['window_id'] = i
            metadata['start_idx'] = start_idx
            metadata['end_idx'] = end_idx
            metadata['window_size'] = window_size
            metadata['stride'] = stride
        
        agg = aggregate_processes(window_processes, indices=window_indices, metadata=metadata)
        aggregations.append(agg)
    
    return aggregations


def multi_scale_aggregation(
    process_series: Union[List[np.ndarray], np.ndarray],
    window_sizes: List[int],
    stride: Optional[int] = None,
    include_metadata: bool = True
) -> Dict[int, List[ProcessAggregation]]:
    """
    Perform rolling window aggregation at multiple scales simultaneously.
    
    This enables multi-scale analysis by aggregating the same process sequence
    at different temporal resolutions.
    
    Args:
        process_series: Sequence of process vectors
        window_sizes: List of window sizes to apply
        stride: Step size (if None, uses stride=1 for all windows)
        include_metadata: Whether to include timing information
        
    Returns:
        Dictionary mapping window_size -> list of ProcessAggregations
        
    Example:
        >>> series = [np.random.rand(5) for _ in range(100)]
        >>> scales = [1, 10, 50, 100]
        >>> multi_agg = multi_scale_aggregation(series, window_sizes=scales)
        >>> for scale, aggs in multi_agg.items():
        ...     print(f"Scale {scale}: {len(aggs)} windows")
    """
    if stride is None:
        stride = 1
    
    results = {}
    
    for window_size in window_sizes:
        try:
            aggs = rolling_window_aggregation(
                process_series, 
                window_size=window_size, 
                stride=stride,
                include_metadata=include_metadata
            )
            results[window_size] = aggs
        except ValueError as e:
            # Skip window sizes that are too large
            results[window_size] = []
    
    return results


######################################################################################
# DISAGGREGATION FUNCTIONS
######################################################################################

def disaggregate_process(
    process: np.ndarray,
    n_components: int,
    strategy: str = 'mixed',
    seed: Optional[int] = None
) -> ProcessDisaggregation:
    """
    Randomly disaggregate a process into n components that sum to original.
    
    This generates a hypothesis about how the aggregate process might have been
    constructed from constituent sub-processes over time.
    
    Strategies:
    - 'coordinate': Assign each active reaction entirely to one component
    - 'value': Split each reaction's value across all components  
    - 'mixed': Randomly choose strategy per reaction
    
    Args:
        process: Process vector to disaggregate
        n_components: Number of components to create (must be >= 1)
        strategy: Disaggregation strategy
        seed: Random seed for reproducibility
        
    Returns:
        ProcessDisaggregation object
        
    Example:
        >>> v = np.array([3.0, 6.0, 1.5, 0.0, 4.5])
        >>> disagg = disaggregate_process(v, n_components=3, strategy='mixed')
        >>> # Verify components sum to original
        >>> assert disagg.verify_sum()
    """
    if seed is not None:
        np.random.seed(seed)
    
    if not isinstance(process, np.ndarray):
        raise TypeError("process must be numpy array")
    
    if n_components < 1:
        raise ValueError(f"n_components must be >= 1, got {n_components}")
    
    # Identify active reactions
    active_reactions = np.where(np.abs(process) > 1e-10)[0]
    k = len(active_reactions)
    
    if k == 0:
        # Zero process - return n zero components
        components = [np.zeros_like(process) for _ in range(n_components)]
        return ProcessDisaggregation(
            original_process=process,
            components=components,
            n_components=n_components,
            strategy=strategy,
            metadata={'n_active_reactions': 0}
        )
    
    if n_components > k and strategy == 'coordinate':
        raise ValueError(
            f"For strategy='coordinate', n_components ({n_components}) cannot exceed "
            f"number of active reactions ({k})"
        )
    
    # Initialize components
    components = [np.zeros_like(process) for _ in range(n_components)]
    
    # Disaggregate each active reaction
    for reaction_idx in active_reactions:
        rate = process[reaction_idx]
        
        if strategy == 'coordinate':
            # Assign entire reaction to one random component
            component_idx = np.random.randint(0, n_components)
            components[component_idx][reaction_idx] = rate
            
        elif strategy == 'value':
            # Split value across all components using Dirichlet distribution
            alphas = np.ones(n_components)
            fractions = np.random.dirichlet(alphas)
            for i in range(n_components):
                components[i][reaction_idx] = rate * fractions[i]
                
        elif strategy == 'mixed':
            # Randomly choose: coordinate or value split
            if np.random.random() < 0.5:
                # Coordinate assignment
                component_idx = np.random.randint(0, n_components)
                components[component_idx][reaction_idx] = rate
            else:
                # Value split
                alphas = np.ones(n_components)
                fractions = np.random.dirichlet(alphas)
                for i in range(n_components):
                    components[i][reaction_idx] = rate * fractions[i]
        else:
            raise ValueError(f"Unknown strategy: {strategy}")
    
    metadata = {
        'n_active_reactions': k,
        'seed': seed
    }
    
    return ProcessDisaggregation(
        original_process=process,
        components=components,
        n_components=n_components,
        strategy=strategy,
        metadata=metadata
    )


def disaggregation_random_sampling(
    process: np.ndarray,
    n_samples: int,
    n_components: int,
    strategy: str = 'mixed',
    seed: Optional[int] = None
) -> List[ProcessDisaggregation]:
    """
    Generate multiple random disaggregations of a process.
    
    This creates n_samples different hypotheses about how the process might have
    been constructed from n_components sub-processes. Essential for pathway
    variability analysis.
    
    Args:
        process: Process vector to disaggregate
        n_samples: Number of different disaggregations to generate (N1 in design)
        n_components: Number of components per disaggregation (N2 in design)
        strategy: Disaggregation strategy (applied to all samples)
        seed: Base random seed (incremented for each sample)
        
    Returns:
        List of ProcessDisaggregation objects
        
    Example:
        >>> v = np.array([5.0, 10.0, 3.0])
        >>> samples = disaggregation_random_sampling(v, n_samples=100, n_components=5)
        >>> # Analyze variability across samples
        >>> first_reactions = [s.components[0][0] for s in samples]
        >>> print(f"First reaction in first component: mean={np.mean(first_reactions):.2f}")
    """
    if n_samples < 1:
        raise ValueError(f"n_samples must be >= 1, got {n_samples}")
    
    disaggregations = []
    
    for i in range(n_samples):
        sample_seed = None if seed is None else seed + i
        
        disagg = disaggregate_process(
            process=process,
            n_components=n_components,
            strategy=strategy,
            seed=sample_seed
        )
        
        # Add sample metadata
        disagg.metadata['sample_id'] = i
        disagg.metadata['total_samples'] = n_samples
        
        disaggregations.append(disagg)
    
    return disaggregations


######################################################################################
# DECOMPOSITION FUNCTIONS (SPECIAL CASE: SINGLE-REACTION COMPONENTS)
######################################################################################

def decompose_to_single_reactions(
    process: np.ndarray,
    tol: float = 1e-10
) -> List[np.ndarray]:
    """
    Decompose process into single-reaction components.
    
    This is a special case of disaggregation where each component contains
    exactly one reaction. These are the atomic building blocks of any process.
    
    Args:
        process: Process vector
        tol: Tolerance for considering reaction active
        
    Returns:
        List of single-reaction process vectors
        
    Example:
        >>> v = np.array([2.0, 0.0, 1.5, 0.3, 0.0])
        >>> components = decompose_to_single_reactions(v)
        >>> len(components)  # Number of active reactions
        3
    """
    if not isinstance(process, np.ndarray):
        raise TypeError("process must be numpy array")
    
    components = []
    active_reactions = np.where(np.abs(process) > tol)[0]
    
    for reaction_idx in active_reactions:
        component = np.zeros_like(process)
        component[reaction_idx] = process[reaction_idx]
        components.append(component)
    
    return components


######################################################################################
# UTILITY FUNCTIONS
######################################################################################

def verify_aggregation(aggregation: ProcessAggregation, tol: float = 1e-10) -> bool:
    """
    Verify that aggregated process equals sum of components.
    
    Args:
        aggregation: ProcessAggregation object to verify
        tol: Numerical tolerance
        
    Returns:
        True if verification passes
    """
    reconstructed = sum(aggregation.component_processes)
    difference = np.linalg.norm(reconstructed - aggregation.aggregated_process)
    return difference < tol


def verify_disaggregation(disaggregation: ProcessDisaggregation, tol: float = 1e-10) -> bool:
    """
    Verify that disaggregation components sum to original.
    
    Args:
        disaggregation: ProcessDisaggregation object to verify
        tol: Numerical tolerance
        
    Returns:
        True if verification passes
    """
    return disaggregation.verify_sum(tol=tol)


def get_aggregation_statistics(aggregations: List[ProcessAggregation]) -> Dict:
    """
    Compute statistics across multiple aggregations.
    
    Args:
        aggregations: List of ProcessAggregation objects
        
    Returns:
        Dictionary with statistics
    """
    if not aggregations:
        return {'n_aggregations': 0}
    
    window_sizes = [agg.window_size for agg in aggregations]
    
    stats = {
        'n_aggregations': len(aggregations),
        'window_sizes': {
            'min': min(window_sizes),
            'max': max(window_sizes),
            'mean': np.mean(window_sizes),
            'unique': len(set(window_sizes))
        },
        'aggregated_processes_shape': aggregations[0].aggregated_process.shape
    }
    
    return stats


def get_disaggregation_statistics(disaggregations: List[ProcessDisaggregation]) -> Dict:
    """
    Compute statistics across multiple disaggregations.
    
    Args:
        disaggregations: List of ProcessDisaggregation objects
        
    Returns:
        Dictionary with statistics
    """
    if not disaggregations:
        return {'n_disaggregations': 0}
    
    n_components_list = [d.n_components for d in disaggregations]
    strategies = [d.strategy for d in disaggregations]
    
    # Verify all disaggregations
    verification_results = [d.verify_sum() for d in disaggregations]
    
    stats = {
        'n_disaggregations': len(disaggregations),
        'n_components': {
            'min': min(n_components_list),
            'max': max(n_components_list),
            'mean': np.mean(n_components_list),
            'unique': len(set(n_components_list))
        },
        'strategies': Counter(strategies),
        'verification': {
            'all_passed': all(verification_results),
            'pass_rate': sum(verification_results) / len(verification_results)
        }
    }
    
    return stats


######################################################################################
# CONE ANALYSIS FUNCTIONS (PROCESS SPACE GEOMETRY)
######################################################################################

def compute_nullspace_vectors(S: np.ndarray) -> List[np.ndarray]:
    """
    Compute basis vectors for the nullspace (kernel) of stoichiometric matrix S.
    
    The nullspace contains all processes v where Sv = 0 (stationary processes).
    These form the cone of steady-state processes.
    
    Args:
        S: Stoichiometric matrix (n_species x n_reactions)
        
    Returns:
        List of nullspace basis vectors (may be empty if nullspace is trivial)
        
    Example:
        >>> S = np.array([[1, -1, 0], [0, 1, -1]])
        >>> null_vecs = compute_nullspace_vectors(S)
        >>> # Verify: S @ v = 0 for each v in null_vecs
    """
    import sympy as sp
    
    if not isinstance(S, np.ndarray):
        raise TypeError("S must be numpy array")
    
    # Convert to sympy matrix for exact nullspace computation
    S_sympy = sp.Matrix(S)
    
    # Compute nullspace
    nullspace = S_sympy.nullspace()
    
    # Convert back to numpy arrays
    nullspace_vectors = []
    for vec in nullspace:
        # Convert sympy vector to numpy, handling rationals
        np_vec = np.array([float(x) for x in vec], dtype=float)
        nullspace_vectors.append(np_vec)
    
    return nullspace_vectors


def compute_feasible_region(
    S: np.ndarray,
    grid_max: Optional[float] = None,
    grid_res: int = 5,
    auto_scale: bool = True,
    tol: float = 1e-6
) -> np.ndarray:
    """
    Compute feasible process vectors in the region Sv >= 0 (self-maintaining cone).
    
    Samples a grid in process space and returns points where Sv >= 0, representing
    processes that don't deplete any species (cognitive domain candidates).
    
    Args:
        S: Stoichiometric matrix
        grid_max: Maximum value for grid (auto-computed if None)
        grid_res: Grid resolution (number of points per dimension)
        auto_scale: If True, scale grid_max based on nullspace
        tol: Tolerance for Sv >= 0 condition
        
    Returns:
        Array of feasible process vectors (n_feasible x n_reactions)
        
    Example:
        >>> feasible = compute_feasible_region(S, grid_max=2.0, grid_res=10)
        >>> # All returned vectors satisfy Sv >= 0
    """
    if not isinstance(S, np.ndarray):
        raise TypeError("S must be numpy array")
    
    n_reactions = S.shape[1]
    
    # Determine grid_max
    if grid_max is None:
        if auto_scale:
            nullspace_vecs = compute_nullspace_vectors(S)
            if nullspace_vecs:
                grid_max = np.max([np.max(np.abs(v)) for v in nullspace_vecs])
            else:
                grid_max = 1.0
        else:
            grid_max = 1.0
    
    # Generate grid
    grid_1d = np.linspace(0, grid_max, grid_res)
    
    # Create meshgrid for all reactions
    grids = np.meshgrid(*[grid_1d] * n_reactions, indexing='ij')
    
    # Flatten to get all combinations
    points = np.column_stack([g.ravel() for g in grids])
    
    # Filter: keep only points where Sv >= -tol (allowing small numerical errors)
    Sv_all = (S @ points.T).T
    is_feasible = np.all(Sv_all >= -tol, axis=1)
    
    feasible_points = points[is_feasible]
    
    return feasible_points


def classify_feasible_points(
    points: np.ndarray,
    S: np.ndarray,
    tol: float = 1e-3
) -> List[str]:
    """
    Classify feasible process points by their operational mode.
    
    This is a vectorized version that classifies multiple points efficiently.
    
    Args:
        points: Array of process vectors (n_points x n_reactions)
        S: Stoichiometric matrix
        tol: Tolerance for classification
        
    Returns:
        List of classification strings for each point
        
    Example:
        >>> feasible = compute_feasible_region(S, grid_res=10)
        >>> classifications = classify_feasible_points(feasible, S)
    """
    if not isinstance(points, np.ndarray):
        raise TypeError("points must be numpy array")
    
    if not isinstance(S, np.ndarray):
        raise TypeError("S must be numpy array")
    
    if points.ndim == 1:
        points = points.reshape(1, -1)
    
    classifications = []
    
    for point in points:
        # Use simplified classification for batch processing
        Sv = S @ point
        v_positive = np.all(point > tol)
        Sv_nonneg = np.all(Sv >= -tol)
        Sv_zero = np.all(np.abs(Sv) <= tol)
        
        if v_positive and Sv_nonneg:
            classification = "Cognitive Domain"
        elif Sv_zero:
            classification = "Stationary Mode"
        elif Sv_nonneg:
            classification = "Overproduction Mode"
        elif np.all(Sv <= tol):
            classification = "Problem"
        else:
            classification = "Challenge"
        
        classifications.append(classification)
    
    return classifications


def analyze_cone(
    S: np.ndarray,
    grid_max: Optional[float] = None,
    grid_res: int = 5,
    classify: bool = True,
    tol: float = 1e-3
) -> Dict:
    """
    Comprehensive analysis of the process cone structure.
    
    Computes nullspace, feasible region, and optionally classifies all feasible points.
    This provides complete geometric characterization of the process space.
    
    Args:
        S: Stoichiometric matrix
        grid_max: Maximum grid value
        grid_res: Grid resolution
        classify: Whether to classify feasible points
        tol: Tolerance
        
    Returns:
        Dictionary containing:
            - 'nullspace_vectors': Basis of steady-state processes
            - 'feasible_points': Points in self-maintaining cone
            - 'Sv_values': Stoichiometric effects for each point
            - 'grid_max': Grid maximum used
            - 'classifications': Point classifications (if classify=True)
            - 'classification_counts': Count by category (if classify=True)
            
    Example:
        >>> cone_data = analyze_cone(S, grid_max=3.0, grid_res=15, classify=True)
        >>> print(f"Nullspace dimension: {len(cone_data['nullspace_vectors'])}")
        >>> print(f"Feasible points: {cone_data['feasible_points'].shape[0]}")
        >>> print(f"Cognitive domain points: {cone_data['classification_counts']['Cognitive Domain']}")
    """
    if not isinstance(S, np.ndarray):
        raise TypeError("S must be numpy array")
    
    # Compute nullspace
    nullspace_vectors = compute_nullspace_vectors(S)
    
    # Compute feasible region
    feasible_points = compute_feasible_region(S, grid_max=grid_max, grid_res=grid_res, auto_scale=True)
    
    # Compute Sv for each point
    if feasible_points.shape[0] > 0:
        Sv_values = (S @ feasible_points.T).T
    else:
        Sv_values = np.array([])
    
    # Determine actual grid_max used
    if grid_max is None:
        if nullspace_vectors:
            grid_max = np.max([np.max(np.abs(v)) for v in nullspace_vectors])
        else:
            grid_max = 1.0
    
    # Build result dictionary
    result = {
        'nullspace_vectors': nullspace_vectors,
        'feasible_points': feasible_points,
        'Sv_values': Sv_values,
        'grid_max': grid_max
    }
    
    # Classify if requested
    if classify and feasible_points.shape[0] > 0:
        classifications = classify_feasible_points(feasible_points, S, tol=tol)
        result['classifications'] = classifications
        
        # Count classifications
        classification_counts = Counter(classifications)
        result['classification_counts'] = dict(classification_counts)
    
    return result


######################################################################################
# EXPORTS
######################################################################################

__all__ = [
    # Data structures
    'ProcessAggregation',
    'ProcessDisaggregation',
    
    # Aggregation
    'aggregate_processes',
    'rolling_window_aggregation',
    'multi_scale_aggregation',
    
    # Disaggregation
    'disaggregate_process',
    'disaggregation_random_sampling',
    'decompose_to_single_reactions',
    
    # Cone analysis (process space geometry)
    'compute_nullspace_vectors',
    'compute_feasible_region',
    'classify_feasible_points',
    'analyze_cone',
    
    # Utilities
    'verify_aggregation',
    'verify_disaggregation',
    'get_aggregation_statistics',
    'get_disaggregation_statistics'
]
