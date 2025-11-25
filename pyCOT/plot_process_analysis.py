"""
Plot Process Analysis Module - Visualization for Semantic Process Analysis

This module provides visualization functions for semantic process analysis results.
Pure visualization layer that consumes analysis outputs and produces interpretable
visual narratives.

Part of the Challenge-Centered Semantic Pipeline (refactored architecture)
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from typing import List, Dict, Optional, Union, Tuple
import warnings

# For heatmaps and advanced plots
try:
    import seaborn as sns
    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False
    warnings.warn("Seaborn not available, some visualizations will be limited")


######################################################################################
# COLOR SCHEMES AND STYLING
######################################################################################

# Standard category colors (can be overridden)
DEFAULT_CATEGORY_COLORS = {
    'conflict': '#d62728',    # Red
    'peace': '#2ca02c',       # Green
    'perception': '#ff7f0e',  # Orange
    'group A': '#1f77b4',     # Blue
    'group B': '#9467bd',     # Purple
}

# Balance state colors
BALANCE_STATE_COLORS = {
    'balanced': '#2ca02c',      # Green
    'overproduced': '#1f77b4',  # Blue
    'depleted': '#d62728'       # Red
}

# Pattern colors
PATTERN_COLORS = {
    'BALANCED_TRANSFORMATION': '#2ca02c',
    'REINFORCING_CYCLE': '#d62728',
    'CASCADE_COLLAPSE': '#ff7f0e',
    'STABILIZING_BALANCE': '#1f77b4'
}


def get_category_color(category: str, color_map: Optional[Dict] = None) -> str:
    """Get color for a category."""
    if color_map and category in color_map:
        return color_map[category]
    if category in DEFAULT_CATEGORY_COLORS:
        return DEFAULT_CATEGORY_COLORS[category]
    # Default fallback
    return '#808080'


######################################################################################
# CATEGORY BEHAVIOR VISUALIZATION
######################################################################################

def plot_category_dynamics_across_scales(
    scale_results: Dict,
    category: str,
    feature: str = 'net_effect',
    save_path: str = '.',
    show_fig: bool = False,
    ax: Optional[plt.Axes] = None,
    show_variance: bool = True,
    color: Optional[str] = None,
    label: Optional[str] = None
) -> plt.Axes:
    """
    Plot how a category's feature changes across aggregation scales.
    
    Creates line plot with optional error bands showing variance across samples.
    
    Args:
        scale_results: Results from Mode 2 or Mode 3 analysis
        category: Category to plot
        feature: Feature to plot ('net_effect', 'balance_quality', etc.)
        ax: Matplotlib axes (created if None)
        show_variance: Whether to show error bands (requires Mode 3 data)
        color: Line color (auto if None)
        label: Legend label
        
    Returns:
        Matplotlib axes
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))
    
    if color is None:
        color = get_category_color(category)
    
    if label is None:
        label = category
    
    # Extract data
    scales = sorted(scale_results['results_by_scale'].keys())
    means = []
    stds = []
    
    for scale in scales:
        scale_data = scale_results['results_by_scale'][scale]
        cat_stats = scale_data.get('category_statistics', {}).get(category, {})
        
        feature_dict = cat_stats.get(feature, {})
        mean_val = feature_dict.get('mean', 0.0)
        std_val = feature_dict.get('std', 0.0) if show_variance else 0.0
        
        means.append(mean_val)
        stds.append(std_val)
    
    # Plot
    ax.plot(scales, means, color=color, linewidth=2, label=label, marker='o')
    
    if show_variance and any(s > 0 for s in stds):
        means_arr = np.array(means)
        stds_arr = np.array(stds)
        ax.fill_between(scales, means_arr - stds_arr, means_arr + stds_arr,
                        color=color, alpha=0.2)
    
    ax.axhline(y=0, color='k', linestyle='--', linewidth=0.8, alpha=0.5)
    ax.set_xlabel('Aggregation Scale (window size)', fontsize=12)
    ax.set_ylabel(f'{feature.replace("_", " ").title()}', fontsize=12)
    ax.set_title(f'{category.title()} - {feature.replace("_", " ").title()} Across Scales', 
                 fontsize=14, fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    if save_path and ax is None:  # Only save if we created our own figure
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
 
    
    if show_fig:
        plt.show()
    return ax


def plot_all_categories_dynamics(
    scale_results: Dict,
    feature: str = 'net_effect',
    figsize: Tuple[int, int] = (14, 8),
    show_variance: bool = True,
    color_map: Optional[Dict] = None,
    save_path: Optional[str] = None,
    show_fig: bool = True
) -> plt.Figure:
    """
    Plot dynamics for all categories on same axes.
    
    Args:
        scale_results: Results from Mode 2 or Mode 3 analysis
        feature: Feature to plot
        figsize: Figure size
        show_variance: Whether to show error bands
        color_map: Custom color mapping
        save_path: Path to save figure
        show_fig: Whether to display figure
        
    Returns:
        Matplotlib figure
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    # Get all categories from first scale
    first_scale = list(scale_results['results_by_scale'].keys())[0]
    categories = list(scale_results['results_by_scale'][first_scale]['category_statistics'].keys())
    
    for category in categories:
        color = get_category_color(category, color_map)
        plot_category_dynamics_across_scales(
            scale_results, category, feature, ax=ax,
            show_variance=show_variance, color=color
        )
    
    plt.tight_layout()
    
    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved figure to {save_path}")
    
    if show_fig:
        plt.show()
    
    return fig


def plot_category_balance_heatmap(
    scale_results: Dict,
    figsize: Tuple[int, int] = (12, 6),
    cmap: str = 'RdYlGn',
    save_path: Optional[str] = None,
    show_fig: bool = True
) -> plt.Figure:
    """
    Create heatmap showing balance states across scales and categories.
    
    Args:
        scale_results: Results from Mode 2 or Mode 3 analysis
        figsize: Figure size
        cmap: Colormap name
        save_path: Path to save figure
        show_fig: Whether to display figure
        
    Returns:
        Matplotlib figure
    """
    # Extract data
    scales = sorted(scale_results['results_by_scale'].keys())
    first_scale = scales[0]
    categories = list(scale_results['results_by_scale'][first_scale]['category_statistics'].keys())
    
    # Build matrix: rows=categories, cols=scales
    data_matrix = np.zeros((len(categories), len(scales)))
    
    for i, category in enumerate(categories):
        for j, scale in enumerate(scales):
            scale_data = scale_results['results_by_scale'][scale]
            cat_stats = scale_data.get('category_statistics', {}).get(category, {})
            net_effect = cat_stats.get('net_effect', {}).get('mean', 0.0)
            data_matrix[i, j] = net_effect
    
    # Plot
    fig, ax = plt.subplots(figsize=figsize)
    
    if HAS_SEABORN:
        sns.heatmap(data_matrix, ax=ax, cmap=cmap, center=0,
                   xticklabels=scales, yticklabels=categories,
                   cbar_kws={'label': 'Net Effect'})
    else:
        im = ax.imshow(data_matrix, cmap=cmap, aspect='auto')
        ax.set_xticks(range(len(scales)))
        ax.set_xticklabels(scales)
        ax.set_yticks(range(len(categories)))
        ax.set_yticklabels(categories)
        plt.colorbar(im, ax=ax, label='Net Effect')
    
    ax.set_xlabel('Aggregation Scale', fontsize=12)
    ax.set_ylabel('Category', fontsize=12)
    ax.set_title('Category Balance Across Scales', fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    
    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved heatmap to {save_path}")
    
    if show_fig:
        plt.show()
    
    return fig


######################################################################################
# PATTERN VISUALIZATION
######################################################################################

def plot_pattern_frequency_heatmap(
    scale_results: Dict,
    figsize: Tuple[int, int] = (12, 6),
    cmap: str = 'YlOrRd',
    save_path: Optional[str] = None,
    show_fig: bool = True
) -> plt.Figure:
    """
    Create heatmap showing pattern frequencies across scales.
    
    Args:
        scale_results: Results from Mode 2 analysis
        figsize: Figure size
        cmap: Colormap
        save_path: Path to save
        show_fig: Whether to display
        
    Returns:
        Matplotlib figure
    """
    scales = sorted(scale_results['results_by_scale'].keys())
    
    # Collect all pattern types
    all_pattern_types = set()
    for scale in scales:
        scale_data = scale_results['results_by_scale'][scale]
        pattern_freq = scale_data.get('pattern_frequency', {})
        all_pattern_types.update(pattern_freq.keys())
    
    pattern_types = sorted(list(all_pattern_types))
    
    if not pattern_types:
        print("No patterns found to plot")
        return None
    
    # Build matrix
    data_matrix = np.zeros((len(pattern_types), len(scales)))
    
    for i, pattern_type in enumerate(pattern_types):
        for j, scale in enumerate(scales):
            scale_data = scale_results['results_by_scale'][scale]
            pattern_freq = scale_data.get('pattern_frequency', {})
            frequency = pattern_freq.get(pattern_type, 0.0)
            data_matrix[i, j] = frequency
    
    # Plot
    fig, ax = plt.subplots(figsize=figsize)
    
    if HAS_SEABORN:
        sns.heatmap(data_matrix, ax=ax, cmap=cmap,
                   xticklabels=scales, yticklabels=pattern_types,
                   cbar_kws={'label': 'Frequency'}, vmin=0, vmax=1,
                   annot=True, fmt='.2f')
    else:
        im = ax.imshow(data_matrix, cmap=cmap, aspect='auto', vmin=0, vmax=1)
        ax.set_xticks(range(len(scales)))
        ax.set_xticklabels(scales)
        ax.set_yticks(range(len(pattern_types)))
        ax.set_yticklabels(pattern_types)
        plt.colorbar(im, ax=ax, label='Frequency')
    
    ax.set_xlabel('Aggregation Scale', fontsize=12)
    ax.set_ylabel('Pattern Type', fontsize=12)
    ax.set_title('Pattern Frequency Across Scales', fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    
    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved pattern heatmap to {save_path}")
    
    if show_fig:
        plt.show()
    
    return fig


def plot_pattern_robustness(
    pathway_analysis: Dict,
    figsize: Tuple[int, int] = (10, 6),
    save_path: Optional[str] = None,
    show_fig: bool = True
) -> plt.Figure:
    """
    Plot pattern robustness from Mode 1 analysis.
    
    Shows how frequently each pattern appears across decomposition samples.
    
    Args:
        pathway_analysis: Results from analyze_pathway_variability (Mode 1)
        figsize: Figure size
        save_path: Save path
        show_fig: Whether to display
        
    Returns:
        Matplotlib figure
    """
    pattern_robustness = pathway_analysis.get('pattern_robustness', {})
    
    if not pattern_robustness:
        print("No pattern robustness data available")
        return None
    
    pattern_types = list(pattern_robustness.keys())
    frequencies = list(pattern_robustness.values())
    
    # Sort by frequency
    sorted_indices = np.argsort(frequencies)[::-1]
    pattern_types = [pattern_types[i] for i in sorted_indices]
    frequencies = [frequencies[i] for i in sorted_indices]
    
    # Plot
    fig, ax = plt.subplots(figsize=figsize)
    
    colors = [PATTERN_COLORS.get(pt, '#808080') for pt in pattern_types]
    bars = ax.barh(range(len(pattern_types)), frequencies, color=colors, alpha=0.7)
    
    ax.set_yticks(range(len(pattern_types)))
    ax.set_yticklabels(pattern_types)
    ax.set_xlabel('Frequency Across Decomposition Samples', fontsize=12)
    ax.set_ylabel('Pattern Type', fontsize=12)
    ax.set_title('Pattern Robustness (Mode 1)', fontsize=14, fontweight='bold')
    ax.set_xlim(0, 1)
    ax.grid(axis='x', alpha=0.3)
    
    # Add percentage labels
    for i, (bar, freq) in enumerate(zip(bars, frequencies)):
        ax.text(freq + 0.02, i, f'{freq:.1%}', va='center', fontsize=10)
    
    plt.tight_layout()
    
    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved robustness plot to {save_path}")
    
    if show_fig:
        plt.show()
    
    return fig


######################################################################################
# DECOMPOSITION ENSEMBLE VISUALIZATION
######################################################################################

def plot_decomposition_ensemble(
    pathway_analysis: Dict,
    category: str,
    feature: str = 'net_effect',
    figsize: Tuple[int, int] = (12, 6),
    alpha: float = 0.1,
    color: Optional[str] = None,
    save_path: Optional[str] = None,
    show_fig: bool = True
) -> plt.Figure:
    """
    Plot ensemble of decomposition pathways for a category.
    
    Shows how category evolves through decomposition steps, with multiple
    semi-transparent trajectories representing different samples.
    
    Args:
        pathway_analysis: Results from analyze_pathway_variability (Mode 1)
        category: Category to visualize
        feature: Feature to plot
        figsize: Figure size
        alpha: Transparency for individual trajectories
        color: Color (auto if None)
        save_path: Save path
        show_fig: Whether to display
        
    Returns:
        Matplotlib figure
    """
    if color is None:
        color = get_category_color(category)
    
    category_sequences = pathway_analysis.get('category_sequences', [])
    
    if not category_sequences:
        print("No decomposition sequences available")
        return None
    
    n_samples = len(category_sequences)
    n_components = len(category_sequences[0][category])
    
    # Extract trajectories
    trajectories = []
    for seq in category_sequences:
        cat_behaviors = seq[category]
        if feature == 'net_effect':
            values = [b.net_effect for b in cat_behaviors]
        elif feature == 'balance_quality':
            values = [b.balance_quality for b in cat_behaviors]
        else:
            values = [0] * len(cat_behaviors)
        trajectories.append(values)
    
    # Plot
    fig, ax = plt.subplots(figsize=figsize)
    
    steps = range(n_components)
    
    # Plot all trajectories with transparency
    for traj in trajectories:
        ax.plot(steps, traj, color=color, alpha=alpha, linewidth=1)
    
    # Plot mean trajectory
    mean_traj = np.mean(trajectories, axis=0)
    ax.plot(steps, mean_traj, color=color, linewidth=3, label='Mean pathway')
    
    # Plot std bands
    std_traj = np.std(trajectories, axis=0)
    ax.fill_between(steps, mean_traj - std_traj, mean_traj + std_traj,
                    color=color, alpha=0.3, label='±1 std')
    
    ax.axhline(y=0, color='k', linestyle='--', linewidth=0.8, alpha=0.5)
    ax.set_xlabel('Decomposition Step', fontsize=12)
    ax.set_ylabel(f'{feature.replace("_", " ").title()}', fontsize=12)
    ax.set_title(f'{category.title()} - Decomposition Ensemble (n={n_samples})',
                 fontsize=14, fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved ensemble plot to {save_path}")
    
    if show_fig:
        plt.show()
    
    return fig


######################################################################################
# THRESHOLD VISUALIZATION
######################################################################################

def plot_thresholds_on_dynamics(
    scale_results: Dict,
    thresholds: List,
    category: str,
    feature: str = 'net_effect',
    figsize: Tuple[int, int] = (12, 6),
    save_path: Optional[str] = None,
    show_fig: bool = True
) -> plt.Figure:
    """
    Plot category dynamics with detected thresholds marked.
    
    Args:
        scale_results: Results from Mode 2 or Mode 3
        thresholds: List of ThresholdDetection objects
        category: Category to visualize
        feature: Feature plotted
        figsize: Figure size
        save_path: Save path
        show_fig: Whether to display
        
    Returns:
        Matplotlib figure
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    # Plot dynamics
    plot_category_dynamics_across_scales(
        scale_results, category, feature, ax=ax, show_variance=True
    )
    
    # Mark thresholds
    category_thresholds = [t for t in thresholds if t.category == category]
    
    for threshold in category_thresholds:
        ax.axvline(x=threshold.scale_value, color='red', linestyle='--',
                  linewidth=2, alpha=0.7)
        
        # Add annotation
        y_pos = ax.get_ylim()[1] * 0.9
        ax.text(threshold.scale_value, y_pos, threshold.threshold_type,
               rotation=90, verticalalignment='top', fontsize=10,
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    ax.set_title(f'{category.title()} with Detected Thresholds', 
                 fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    
    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved threshold plot to {save_path}")
    
    if show_fig:
        plt.show()
    
    return fig


######################################################################################
# MODE COMPARISON VISUALIZATION
######################################################################################

def plot_three_mode_comparison(
    mode1_results: Dict,
    mode2_results: Dict,
    mode3_results: Dict,
    category: str,
    feature: str = 'net_effect',
    figsize: Tuple[int, int] = (18, 6),
    save_path: Optional[str] = None,
    show_fig: bool = True
) -> plt.Figure:
    """
    Create three-panel comparison showing all three analytical modes.
    
    Args:
        mode1_results: Results from Mode 1 (pathway variability)
        mode2_results: Results from Mode 2 (temporal scale)
        mode3_results: Results from Mode 3 (full robustness)
        category: Category to compare
        feature: Feature to analyze
        figsize: Figure size
        save_path: Save path
        show_fig: Whether to display
        
    Returns:
        Matplotlib figure
    """
    fig = plt.figure(figsize=figsize)
    gs = GridSpec(1, 3, figure=fig, wspace=0.3)
    
    # Panel 1: Mode 1 - Pathway ensemble
    ax1 = fig.add_subplot(gs[0, 0])
    category_sequences = mode1_results.get('category_sequences', [])
    if category_sequences:
        trajectories = []
        for seq in category_sequences:
            cat_behaviors = seq[category]
            values = [getattr(b, feature) for b in cat_behaviors]
            trajectories.append(values)
        
        color = get_category_color(category)
        steps = range(len(trajectories[0]))
        
        for traj in trajectories:
            ax1.plot(steps, traj, color=color, alpha=0.1, linewidth=1)
        
        mean_traj = np.mean(trajectories, axis=0)
        ax1.plot(steps, mean_traj, color=color, linewidth=3, label='Mean')
        
        ax1.axhline(y=0, color='k', linestyle='--', alpha=0.5)
        ax1.set_xlabel('Decomposition Step')
        ax1.set_ylabel(feature.replace('_', ' ').title())
        ax1.set_title('Mode 1: Pathway Variability')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
    
    # Panel 2: Mode 2 - Temporal scale
    ax2 = fig.add_subplot(gs[0, 1])
    plot_category_dynamics_across_scales(
        mode2_results, category, feature, ax=ax2, show_variance=False
    )
    ax2.set_title('Mode 2: Temporal Scale')
    
    # Panel 3: Mode 3 - Full robustness
    ax3 = fig.add_subplot(gs[0, 2])
    plot_category_dynamics_across_scales(
        mode3_results, category, feature, ax=ax3, show_variance=True
    )
    ax3.set_title('Mode 3: Multi-Scale Robustness')
    
    fig.suptitle(f'{category.title()} - Three-Mode Comparison', 
                 fontsize=16, fontweight='bold', y=1.02)
    
    plt.tight_layout()
    
    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved three-mode comparison to {save_path}")
    
    if show_fig:
        plt.show()
    
    return fig


######################################################################################
# CONE VISUALIZATION (PROCESS SPACE GEOMETRY)
######################################################################################

def plot_cone_and_region(
    S: np.ndarray,
    grid_max: Optional[float] = None,
    grid_res: int = 5,
    axis_names: Optional[List[str]] = None,
    show: bool = True,
    extra_vectors: Optional[List[np.ndarray]] = None,
    extra_vector_labels: Optional[List[str]] = None,
    extra_vector_colors: Optional[List[str]] = None,
    save_dir: str = "./visualizations/cone_projections_3D"
) -> Tuple[List[str], np.ndarray]:
    """
    Generate 3D projections of the process cone and feasible region.
    
    Visualizes the nullspace cone and the feasible region (Sv >= 0) in 3D projections.
    Points are colored by their classification (Cognitive Domain, Stationary, etc.).
    
    Args:
        S: Stoichiometric matrix
        grid_max: Maximum value for grid
        grid_res: Grid resolution
        axis_names: Names for reaction axes
        show: Whether to display plots
        extra_vectors: Additional vectors to plot
        extra_vector_labels: Labels for extra vectors
        extra_vector_colors: Colors for extra vectors
        save_dir: Directory to save plots
        
    Returns:
        Tuple of (list of saved file paths, feasible points array)
        
    Example:
        >>> from process_structure import analyze_cone
        >>> # S is your stoichiometric matrix
        >>> files, points = plot_cone_and_region(S, grid_max=3, grid_res=15)
    """
    from pyCOT.process_structure import analyze_cone
    import os
    from mpl_toolkits.mplot3d import Axes3D
    
    m, n = S.shape  # species x reactions
    
    if axis_names is None:
        axis_names = [f"v{i+1}" for i in range(n)]
    
    if extra_vector_labels is None and extra_vectors is not None:
        extra_vector_labels = [f"Extra {i+1}" for i in range(len(extra_vectors))]
    
    if extra_vector_colors is None and extra_vectors is not None:
        extra_vector_colors = ['orange', 'cyan', 'purple', 'brown', 'yellow', 'pink']
        extra_vector_colors = extra_vector_colors + ['gray'] * max(0, len(extra_vectors) - 6)
    
    os.makedirs(save_dir, exist_ok=True)
    
    # Adjust grid_max for extra vectors
    if extra_vectors is not None:
        extra_max = np.max([np.max(np.abs(v)) for v in extra_vectors])
        grid_max = max(grid_max or 0, extra_max)
    
    # Perform comprehensive cone analysis
    cone_data = analyze_cone(S, grid_max=grid_max, grid_res=grid_res, classify=True)
    
    null_vectors = cone_data['nullspace_vectors']
    points_pos = cone_data['feasible_points']
    classifications = cone_data.get('classifications', [])
    grid_max = cone_data['grid_max']
    
    print(f"Nullspace dimension: {len(null_vectors)}")
    print(f"Feasible points found: {points_pos.shape[0]}")
    if classifications:
        print(f"Classification breakdown: {cone_data.get('classification_counts', {})}")
    
    # Color mapping for classifications
    class_colors = {
        'Cognitive Domain': 'green',
        'Stationary Mode': 'blue',
        'Overproduction Mode': 'cyan',
        'Problem': 'red',
        'Challenge': 'orange'
    }
    
    # Assign colors to points
    point_colors = [class_colors.get(c, 'gray') for c in classifications]
    
    saved_files = []
    
    # Generate all 3D projections
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                fig = plt.figure(figsize=(10, 8))
                ax = fig.add_subplot(111, projection='3d')
                
                # Plot feasible points
                if points_pos.shape[0] > 0:
                    ax.scatter(points_pos[:, i], points_pos[:, j], points_pos[:, k],
                             c=point_colors, s=5, alpha=0.6)
                
                # Plot nullspace vectors
                for idx, null_vec in enumerate(null_vectors):
                    ax.quiver(0, 0, 0, 
                            null_vec[i], null_vec[j], null_vec[k],
                            color='red', arrow_length_ratio=0.1,
                            linewidth=2, label=f"Nullspace {idx+1}")
                
                # Plot extra vectors if provided
                if extra_vectors is not None:
                    for k_idx, vec in enumerate(extra_vectors):
                        ax.quiver(0, 0, 0,
                                vec[i], vec[j], vec[k],
                                color=extra_vector_colors[k_idx],
                                arrow_length_ratio=0.1,
                                label=extra_vector_labels[k_idx])
                
                ax.set_xlim([-grid_max*0.1, grid_max*1.1])
                ax.set_ylim([-grid_max*0.1, grid_max*1.1])
                ax.set_zlim([-grid_max*0.1, grid_max*1.1])
                
                ax.set_xlabel(axis_names[i])
                ax.set_ylabel(axis_names[j])
                ax.set_zlabel(axis_names[k])
                ax.set_title(f"Process Cone 3D ({axis_names[i]}, {axis_names[j]}, {axis_names[k]})")
                ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5))
                
                filename = os.path.join(save_dir, f"cone_projection_{i+1}_{j+1}_{k+1}.png")
                fig.savefig(filename, dpi=150, bbox_inches='tight')
                saved_files.append(filename)
                
                if show:
                    plt.show()
                else:
                    plt.close(fig)
    
    print(f"Saved {len(saved_files)} cone projections to: {save_dir}")
    
    return saved_files, points_pos


def plot_cone_2d_projections(
    S: np.ndarray,
    grid_max: Optional[float] = None,
    grid_res: int = 20,
    axis_names: Optional[List[str]] = None,
    figsize: Tuple[int, int] = (15, 10),
    save_path: Optional[str] = None,
    show_fig: bool = True
) -> plt.Figure:
    """
    Create 2D projections of process cone in a multi-panel figure.
    
    Args:
        S: Stoichiometric matrix
        grid_max: Maximum grid value
        grid_res: Grid resolution
        axis_names: Reaction axis names
        figsize: Figure size
        save_path: Path to save figure
        show_fig: Whether to display
        
    Returns:
        Matplotlib figure
    """
    from pyCOT.process_structure import analyze_cone
    
    n_reactions = S.shape[1]
    
    if axis_names is None:
        axis_names = [f"v{i+1}" for i in range(n_reactions)]
    
    # Analyze cone
    cone_data = analyze_cone(S, grid_max=grid_max, grid_res=grid_res, classify=True)
    points = cone_data['feasible_points']
    classifications = cone_data.get('classifications', [])
    
    if points.shape[0] == 0:
        print("No feasible points found")
        return None
    
    # Color mapping
    class_colors = {
        'Cognitive Domain': 'green',
        'Stationary Mode': 'blue',
        'Overproduction Mode': 'cyan',
        'Problem': 'red',
        'Challenge': 'orange'
    }
    point_colors = [class_colors.get(c, 'gray') for c in classifications]
    
    # Create subplots for all 2D projections
    n_pairs = n_reactions * (n_reactions - 1) // 2
    n_cols = min(3, n_pairs)
    n_rows = (n_pairs + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    axes = axes.flatten() if n_pairs > 1 else [axes]
    
    pair_idx = 0
    for i in range(n_reactions):
        for j in range(i+1, n_reactions):
            if pair_idx < len(axes):
                ax = axes[pair_idx]
                ax.scatter(points[:, i], points[:, j], 
                          c=point_colors, s=10, alpha=0.6)
                ax.set_xlabel(axis_names[i])
                ax.set_ylabel(axis_names[j])
                ax.set_title(f"{axis_names[i]} vs {axis_names[j]}")
                ax.grid(True, alpha=0.3)
                pair_idx += 1
    
    # Hide unused subplots
    for idx in range(pair_idx, len(axes)):
        axes[idx].set_visible(False)
    
    # Add legend
    legend_elements = [mpatches.Patch(facecolor=color, label=label) 
                      for label, color in class_colors.items()]
    fig.legend(handles=legend_elements, loc='upper right')
    
    plt.suptitle('Process Cone 2D Projections', fontsize=16, fontweight='bold')
    plt.tight_layout()
    
    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved 2D projections to {save_path}")
    
    if show_fig:
        plt.show()
    
    return fig


######################################################################################
# UTILITY FUNCTIONS
######################################################################################

def setup_publication_style():
    """Set matplotlib style for publication-quality figures."""
    plt.rcParams['font.size'] = 11
    plt.rcParams['axes.labelsize'] = 12
    plt.rcParams['axes.titlesize'] = 14
    plt.rcParams['xtick.labelsize'] = 10
    plt.rcParams['ytick.labelsize'] = 10
    plt.rcParams['legend.fontsize'] = 10
    plt.rcParams['figure.titlesize'] = 16
    plt.rcParams['figure.dpi'] = 100
    plt.rcParams['savefig.dpi'] = 300
    plt.rcParams['font.family'] = 'sans-serif'


######################################################################################
# EXPORTS
######################################################################################

__all__ = [
    # Category dynamics
    'plot_category_dynamics_across_scales',
    'plot_all_categories_dynamics',
    'plot_category_balance_heatmap',
    
    # Pattern visualization
    'plot_pattern_frequency_heatmap',
    'plot_pattern_robustness',
    
    # Decomposition
    'plot_decomposition_ensemble',
    
    # Thresholds
    'plot_thresholds_on_dynamics',
    
    # Mode comparison
    'plot_three_mode_comparison',
    
    # Cone visualization (process space geometry)
    'plot_cone_and_region',
    'plot_cone_2d_projections',
    
    # Utilities
    'setup_publication_style',
    'get_category_color'
]
