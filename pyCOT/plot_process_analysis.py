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
import pandas as pd
import os
import itertools
import math
from collections import Counter
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# Import from pyCOT modules
from pyCOT.process_analyzer import (
    classify_process_mode, 
    get_intervals_by_category,
    analyze_cone,
    rescale_process_time_series
)
from pyCOT.simulations import simulation

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
# # FunciÃ³n para clasificar un vector de proceso v en categorÃ­as base y extendidas
def plot_process_types_histogram(flux_vector, S, 
                                xlabel="Tipo de Proceso", ylabel="Frecuencia",
                                title="Histograma de Tipos de Proceso", 
                                excel_filename="classified_processes.xlsx",
                                filename="histograma.png", 
                                save_figure=True, 
                                ax=None, show_fig=False):
    """
    Clasifica, grafica histograma, guarda resultados en Excel y
    retorna tambiÃ©n las frecuencias (conteos) de cada tipo de proceso.
    """
    if isinstance(flux_vector, np.ndarray):
        flux_vector = pd.DataFrame(flux_vector, columns=[f"v{i+1}" for i in range(flux_vector.shape[1])])
        flux_vector.insert(0, "Time", range(len(flux_vector)))
    elif not isinstance(flux_vector, pd.DataFrame):
        raise TypeError("flux_vector must be a DataFrame or a NumPy ndarray.")    
    if not isinstance(S, np.ndarray):
        raise TypeError("S must be a NumPy array.")    

    # --- Classify ---
    flux_values = flux_vector.iloc[:, 1:] # remove Time column
    process_types = []
    for _, row in flux_values.iterrows(): # Iterate over each row of the DataFrame
        v = row.to_numpy()                # Convert the row to a NumPy array
        cat = classify_process_mode(v, S) # Classify the process mode
        cat_str = ",".join([cat[0]]) if cat else "None" # Join categories into a string
        process_types.append(cat_str)     # Append to the list
    
    # Calculate S*v 
    Sv_matrix = flux_values.apply(lambda v: S @ v.to_numpy(), axis=1)
    Sv_expanded = pd.DataFrame(Sv_matrix.tolist(), 
                               columns=[f"S*v_{i+1}" for i in range(S.shape[0])],
                               index=flux_vector.index)
    classified_df = pd.concat([flux_vector, Sv_expanded], axis=1)
    classified_df["Process_Type"] = process_types
    # print("Procesos clasificados:\n",classified_df)

    # --- Count ---
    process_counts = Counter(process_types) 
    category_order = color_map.keys()
    labels = [cat for cat in category_order if cat in process_counts]
    counts = [process_counts[cat] for cat in labels]
    colors = [color_map.get(label, "grey") for label in labels]

    # --- Plot ---
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))
    else:
        fig = ax.get_figure()

    bars = ax.bar(labels, counts, color=colors)
    for bar in bars:
        yval = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2.0, yval, int(yval),
                va='bottom', ha='center', fontweight='bold') 

    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.tick_params(axis="x", rotation=30)
    ax.grid(axis='y', linestyle='--', alpha=0.7)

    # --- Save to Excel ---
    out_dir = "./visualizations/process_classification"
    os.makedirs(out_dir, exist_ok=True)
    filepath_excel = os.path.join(out_dir, excel_filename)

    with pd.ExcelWriter(filepath_excel, engine="openpyxl") as writer:
        classified_df.to_excel(writer, sheet_name="All_Processes", index=False)
        for category, group in classified_df.groupby("Process_Type", sort=False):
            group.to_excel(writer, sheet_name=category[:31], index=False)
    
    print(f"Saved classified processes in: {filepath_excel}")

    # --- Guardar figura ---
    if save_figure:
        filepath_fig = os.path.join(out_dir, filename)
        fig.savefig(filepath_fig, dpi=300, bbox_inches="tight")
        print(f"Histogram saved in: {filepath_fig}")

    if show_fig:
        plt.show()

    # --- Retornar tambiÃ©n las frecuencias ---
    process_frequencies = dict(process_counts)  # convertir Counter â†’ dict

    return fig, ax, classified_df, process_frequencies 

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
                    color=color, alpha=0.3, label='Â±1 std')
    
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
# PROCESS TYPE CLASSIFICATION VISUALIZATION
######################################################################################

# Color map for process categories
PROCESS_TYPE_COLOR_MAP = {
    "Stationary Mode": "cyan", 
    "Overproduction Mode": "blue",
    "Cognitive Domain": "green",
    "Challenge": "orange",
    "Problem": "red",        
    "Other": "grey",
    "None": "black"
}


def plot_process_types_histogram(flux_vector, S, 
                                xlabel="Tipo de Proceso", ylabel="Frecuencia",
                                title="Histograma de Tipos de Proceso", 
                                excel_filename="classified_processes.xlsx",
                                filename="histograma.png", 
                                save_figure=True, 
                                ax=None, show_fig=False):
    """
    Classify processes, plot histogram, save results to Excel and
    return frequencies (counts) of each process type.
    
    Args:
        flux_vector: DataFrame or NumPy array with flux values
        S: Stoichiometric matrix (NumPy array)
        xlabel: X-axis label
        ylabel: Y-axis label
        title: Plot title
        excel_filename: Excel file name for classified processes
        filename: PNG file name for histogram
        save_figure: Whether to save the figure
        ax: Matplotlib axes (created if None)
        show_fig: Whether to display the figure
        
    Returns:
        fig, ax, classified_df, process_frequencies
    """
    if isinstance(flux_vector, np.ndarray):
        flux_vector = pd.DataFrame(flux_vector, columns=[f"v{i+1}" for i in range(flux_vector.shape[1])])
        flux_vector.insert(0, "Time", range(len(flux_vector)))
    elif not isinstance(flux_vector, pd.DataFrame):
        raise TypeError("flux_vector must be a DataFrame or a NumPy ndarray.")    
    if not isinstance(S, np.ndarray):
        raise TypeError("S must be a NumPy array.")    

    # --- Classify ---
    flux_values = flux_vector.iloc[:, 1:] # remove Time column
    process_types = []
    for _, row in flux_values.iterrows(): # Iterate over each row of the DataFrame
        v = row.to_numpy()                # Convert the row to a NumPy array
        cat = classify_process_mode(v, S) # Classify the process mode
        cat_str = ",".join([cat[0]]) if cat else "None" # Join categories into a string
        process_types.append(cat_str)     # Append to the list
    
    # Calculate S*v 
    Sv_matrix = flux_values.apply(lambda v: S @ v.to_numpy(), axis=1)
    Sv_expanded = pd.DataFrame(Sv_matrix.tolist(), 
                               columns=[f"S*v_{i+1}" for i in range(S.shape[0])],
                               index=flux_vector.index)
    classified_df = pd.concat([flux_vector, Sv_expanded], axis=1)
    classified_df["Process_Type"] = process_types

    # --- Count ---
    process_counts = Counter(process_types) 
    category_order = PROCESS_TYPE_COLOR_MAP.keys()
    labels = [cat for cat in category_order if cat in process_counts]
    counts = [process_counts[cat] for cat in labels]
    colors = [PROCESS_TYPE_COLOR_MAP.get(label, "grey") for label in labels]

    # --- Plot ---
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))
    else:
        fig = ax.get_figure()

    bars = ax.bar(labels, counts, color=colors)
    for bar in bars:
        yval = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2.0, yval, int(yval),
                va='bottom', ha='center', fontweight='bold') 

    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.tick_params(axis="x", rotation=30)
    ax.grid(axis='y', linestyle='--', alpha=0.7)

    # --- Save to Excel ---
    out_dir = "./visualizations/process_classification"
    os.makedirs(out_dir, exist_ok=True)
    filepath_excel = os.path.join(out_dir, excel_filename)

    with pd.ExcelWriter(filepath_excel, engine="openpyxl") as writer:
        classified_df.to_excel(writer, sheet_name="All_Processes", index=False)
        for category, group in classified_df.groupby("Process_Type", sort=False):
            group.to_excel(writer, sheet_name=category[:31], index=False)
    
    print(f"Saved classified processes in: {filepath_excel}")

    # --- Save figure ---
    if save_figure:
        filepath_fig = os.path.join(out_dir, filename)
        fig.savefig(filepath_fig, dpi=300, bbox_inches="tight")
        print(f"Histogram saved in: {filepath_fig}")

    if show_fig:
        plt.show()

    # --- Return frequencies ---
    process_frequencies = dict(process_counts)  # Convert Counter → dict

    return fig, ax, classified_df, process_frequencies


def plot_cone_and_region_3d(S,
                         grid_max=None, grid_res=5,
                         axis_names=None, show=True,
                         extra_vector=None, 
                         extra_vector_labels=None,
                         extra_vector_colors=None):   
    """
    Generate 3D projections of the cone defined by nullspace vectors
    and the feasible region Sv > 0 for a stoichiometric matrix S.
    Feasible points are classified with classify_process_mode() and colored
    according to their category.
    Also saves to Excel file the feasible points (v), S*v and
    Process_Type in a main sheet, and in separate sheets by category.
    
    NOTE: This function uses analyze_cone() from process_analysis for cone construction
    
    Args:
        S: Stoichiometric matrix
        grid_max: Maximum grid value
        grid_res: Grid resolution
        axis_names: Reaction axis names
        show: Whether to display plots
        extra_vector: Extra vectors to plot
        extra_vector_labels: Labels for extra vectors
        extra_vector_colors: Colors for extra vectors
        
    Returns:
        saved_files: List of saved file paths
        points_pos: Feasible points array
    """
    m = S.shape[0] # number of species
    n = S.shape[1] # number of reactions
    
    species_names = S.species # Species

    if axis_names is None:
        axis_names = [f"v{i+1}" for i in range(n)]
    
    if extra_vector_labels is None:
        extra_vector_labels = [f"Extra vector {i+1}" for i in range(len(extra_vector or []))]
    
    if extra_vector_colors is None:
        extra_vector_colors = ['orange', 'cyan', 'purple', 'brown', 'yellow', 'pink'] + \
                              ['gray'] * max(0, len(extra_vector or []) - 6)
    
    out_dir = "./visualizations/cone_projections_3D"
    os.makedirs(out_dir, exist_ok=True)

    # ==== REFACTORED: Use analyze_cone() from process_analysis ====
    # Adjust grid_max if extra vectors provided
    if extra_vector is not None:
        if isinstance(extra_vector, np.ndarray) and extra_vector.ndim == 1:
            extra_vector = [extra_vector]
        extra_max = np.max(np.abs(np.vstack(extra_vector))) if extra_vector else 0
        if grid_max is None:
            grid_max = extra_max
        else:
            grid_max = max(grid_max, extra_max)
    
    # Perform comprehensive cone analysis
    cone_data = analyze_cone(S, grid_max=grid_max, grid_res=grid_res, classify=True)
    
    null_vectors = cone_data['nullspace_vectors']
    points_pos = cone_data['feasible_points']
    classifications = cone_data.get('classifications', []) # may be empty
    grid_max = cone_data['grid_max']
    
    print(f"Null space vectors:\n {null_vectors}")
    
    # ---- Export classified processes to Excel ----
    if points_pos.shape[0] > 0:
        Sv = cone_data['Sv_values']
        df_points = pd.DataFrame(points_pos, columns=axis_names)
        df_Sv = pd.DataFrame(Sv, columns=[f"S*v=x{i+1}" for i in range(m)])
        df_class = pd.DataFrame({"Process_Type": classifications})

        # Combine all data
        df_all = pd.concat([df_points, df_Sv, df_class], axis=1)
        
        positive_species_list = []
        for idx, row in df_Sv.iterrows(): # Each row of Sv
            # Indices where Sv > 0
            pos_idx = np.where(row.values > 0)[0]
            # obtain corresponding species names
            if len(pos_idx) > 0:
                # Species with positive values
                pos_species = [species_names[i] for i in pos_idx]
            else:
                pos_species = []
            
            positive_species_list.append(pos_species)

        df_all["Abstractions"] = positive_species_list

        # Save to Excel with multiple sheets
        excel_path = os.path.join(out_dir, "classified_processes.xlsx")
        with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
            # Main sheet with all data
            df_all.to_excel(writer, sheet_name="All_Processes", index=False)
            
            # Separate sheets by category
            for category, group in df_all.groupby("Process_Type"):
                safe_name = str(category)[:31]  # Excel limits sheet names to 31 characters
                group.to_excel(writer, sheet_name=safe_name, index=False)
        
        print(f"Process vectors saved in: {excel_path}")
    else:
        print("No feasible points were found Sv>0.")

    # ---- 3. 3D projections ----
    saved_files = []
    for (i, j, k) in itertools.combinations(range(n), 3):
        fig = plt.figure(figsize=(7,6))
        ax = fig.add_subplot(111, projection='3d')

        # Cone (Sv=0)
        if len(null_vectors) == 0:
            ax.scatter([0], [0], [0], color='lightblue', s=50, label="Sv=0 (origen)")
        elif len(null_vectors) == 1:
            v1 = null_vectors[0]
            ax.quiver(0, 0, 0, *(v1[[i,j,k]]), color='black', linewidth=2, 
                      arrow_length_ratio=0.1, label="Stationary Mode 1")
        else:
            v1, v2 = null_vectors[:2] 
            tri = np.array([[0,0,0], v1[[i,j,k]], v2[[i,j,k]]])
            poly = Poly3DCollection([tri], alpha=0.4, facecolor='lightblue', 
                                   edgecolor='lightblue', label="Sv=0")
            ax.add_collection3d(poly)
            colors = ['black', 'lightblue', 'red', 'orange', 'purple', 'brown'] + \
                     ['purple'] * (len(null_vectors) - 2)
            for idx, v in enumerate(null_vectors):
                ax.quiver(0, 0, 0, *(v[[i,j,k]]), color=colors[idx], linewidth=2, 
                          arrow_length_ratio=0.1, label=f"Stationary Mode {idx+1}")

        # Feasible region Sv>0 classified
        if points_pos.shape[0] > 0:
            proj_pos = points_pos[:, [i, j, k]]
         
            categories = [c.split(",")[0] if c else "None" for c in classifications]

            for cat in sorted(set(categories)):
                mask = [c == cat for c in categories]
                ax.scatter(proj_pos[mask, 0],
                           proj_pos[mask, 1],
                           proj_pos[mask, 2],
                           alpha=0.6, s=20,
                           c=PROCESS_TYPE_COLOR_MAP.get(cat, "grey"),
                           label=cat)

        if extra_vector is not None:
            for k_idx, vec in enumerate(extra_vector):
                extra_proj = vec[[i, j, k]]
                ax.quiver(0, 0, 0, *extra_proj,
                          color=extra_vector_colors[k_idx % len(extra_vector_colors)],
                          linewidth=2, arrow_length_ratio=0.1,
                          label=extra_vector_labels[k_idx] if k_idx < len(extra_vector_labels) else f"Extra vector {k_idx+1}")

        ax.set_xlim([-grid_max*0.1, grid_max*1.1])
        ax.set_ylim([-grid_max*0.1, grid_max*1.1])
        ax.set_zlim([-grid_max*0.1, grid_max*1.1])

        ax.set_xlabel(axis_names[i])
        ax.set_ylabel(axis_names[j])
        ax.set_zlabel(axis_names[k])
        ax.set_title(f"Proyección 3D ({axis_names[i]}, {axis_names[j]}, {axis_names[k]})")
        ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5))

        filename = os.path.join(out_dir, f"Cone_and_region_{i+1}_{j+1}_{k+1}.png")
        fig.savefig(filename, dpi=150, bbox_inches='tight')
        saved_files.append(filename)
        if show:
            plt.show()
        else:
            plt.close(fig)

    print(f"Se guardaron {len(saved_files)} imágenes en: {out_dir}")
    return saved_files, points_pos


def plot_series_with_domain_intervals(time_series, flux_vector, S,
                                      title="Serie de Tiempo de Concentraciones",
                                      save_figure=False,
                                      ax=None,
                                      show_fig=False):
    """
    Plot time series with intervals for Cognitive Domain, Stationary Mode, and Problem.
    
    Args:
        time_series: DataFrame with Time column and concentration columns
        flux_vector: DataFrame with Time column and flux columns
        S: Stoichiometric matrix
        title: Plot title
        save_figure: Whether to save the figure
        ax: Matplotlib axes (created if None)
        show_fig: Whether to display the figure
        
    Returns:
        fig, ax
    """
    from pyCOT.plot_dynamics import plot_series_ode
    
    fig, ax = plot_series_ode(time_series, title=title, save_figure=save_figure, ax=ax)

    times = time_series["Time"].to_numpy()
    flux_values = flux_vector.iloc[:, 1:]
    process_types = [classify_process_mode(v.to_numpy(), S) for _, v in flux_values.iterrows()]

    # Using get_intervals_by_category from process_analysis

    # Masks
    is_cd = np.array(["Cognitive Domain" in cat for cat in process_types])
    is_sm = np.array(["Stationary Mode" in cat for cat in process_types])
    is_pb = np.array(["Problem" in cat for cat in process_types])

    # Intervals
    cd_intervals = get_intervals_by_category(times, process_types, "Cognitive Domain")
    sm_intervals = get_intervals_by_category(times, process_types, "Stationary Mode")
    pb_intervals = get_intervals_by_category(times, process_types, "Problem")

    # === Plot intervals ===
    for (t_start, t_end) in cd_intervals:
        ax.axvspan(t_start, t_end, color="green", alpha=0.2)
        ax.axvline(x=t_start, color="green", linestyle="--", alpha=0.8)
        ax.axvline(x=t_end, color="green", linestyle="--", alpha=0.8)

    for (t_start, t_end) in sm_intervals:
        ax.axvspan(t_start, t_end, color="cyan", alpha=0.15)
        ax.axvline(x=t_start, color="cyan", linestyle="--", alpha=0.7)
        ax.axvline(x=t_end, color="cyan", linestyle="--", alpha=0.7)

    for (t_start, t_end) in pb_intervals:
        ax.axvspan(t_start, t_end, color="red", alpha=0.25)
        ax.axvline(x=t_start, color="red", linestyle="--", alpha=0.8)
        ax.axvline(x=t_end, color="red", linestyle="--", alpha=0.8)

    # Legend
    ax.plot([], [], color="green", linestyle="--", label="Cognitive Domain")
    ax.plot([], [], color="cyan", linestyle="--", label="Stationary Mode")
    ax.plot([], [], color="red", linestyle="--", label="Problem")
    ax.legend(loc="upper right")

    if show_fig:
        plt.show()

    return fig, ax


def plot_flux_with_domain_intervals(flux_vector, S,
                                    title="Serie de Tiempo de Flujos",
                                    save_figure=False,
                                    ax=None,
                                    show_fig=False): 
    """
    Plot flux time series with intervals for Cognitive Domain, Stationary Mode, and Problem.
    
    Args:
        flux_vector: DataFrame with Time column and flux columns
        S: Stoichiometric matrix
        title: Plot title
        save_figure: Whether to save the figure
        ax: Matplotlib axes (created if None)
        show_fig: Whether to display the figure
        
    Returns:
        fig, ax
    """
    from pyCOT.plot_dynamics import plot_series_ode
    
    fig, ax = plot_series_ode(flux_vector, title=title, save_figure=save_figure, ax=ax)
    print("\nflux_vector_shape =", flux_vector.shape)
    
    flux_values = flux_vector.iloc[:, 1:]
    process_types = [classify_process_mode(v.to_numpy(), S) for _, v in flux_values.iterrows()]
    times = flux_vector["Time"].to_numpy()

    # Using get_intervals_by_category from process_analysis

    # Masks
    is_cd = np.array(["Cognitive Domain" in cat for cat in process_types])
    is_sm = np.array(["Stationary Mode" in cat for cat in process_types])
    is_pb = np.array(["Problem" in cat for cat in process_types])

    # Intervals
    cd_intervals = get_intervals_by_category(times, process_types, "Cognitive Domain")
    sm_intervals = get_intervals_by_category(times, process_types, "Stationary Mode")
    pb_intervals = get_intervals_by_category(times, process_types, "Problem")

    # Plot intervals
    for (t_start, t_end) in cd_intervals:
        ax.axvspan(t_start, t_end, color="green", alpha=0.2)
        ax.axvline(x=t_start, color="green", linestyle="--", alpha=0.8)
        ax.axvline(x=t_end, color="green", linestyle="--", alpha=0.8)

    for (t_start, t_end) in sm_intervals:
        ax.axvspan(t_start, t_end, color="cyan", alpha=0.15)
        ax.axvline(x=t_start, color="cyan", linestyle="--", alpha=0.7)
        ax.axvline(x=t_end, color="cyan", linestyle="--", alpha=0.7)

    for (t_start, t_end) in pb_intervals:
        ax.axvspan(t_start, t_end, color="red", alpha=0.25)
        ax.axvline(x=t_start, color="red", linestyle="--", alpha=0.8)
        ax.axvline(x=t_end, color="red", linestyle="--", alpha=0.8)

    # Statistics
    def print_stats(name, intervals):
        if len(intervals) == 0:
            print(f"No se detectaron intervalos de {name}.")
            return
        durations = [t2 - t1 for (t1, t2) in intervals]
        print(f"\n{name}: {len(intervals)} intervalos")
        print(f"  T_min = {min(durations):.4f}")
        print(f"  T_max = {max(durations):.4f}")
        print(f"  T_avg = {np.mean(durations):.4f}")

    print_stats("Cognitive Domain", cd_intervals)
    print_stats("Stationary Mode", sm_intervals)
    print_stats("Problem", pb_intervals)

    # Legend
    ax.plot([], [], color="green", linestyle="--", label="Cognitive Domain")
    ax.plot([], [], color="cyan", linestyle="--", label="Stationary Mode")
    ax.plot([], [], color="red", linestyle="--", label="Problem")
    ax.legend(loc="upper right")

    if show_fig:
        plt.show()

    return fig, ax


def histogram_flux_sum(S, flux1, flux2, 
                       title="Histograma de Tipos de Proceso (Suma Flux1 + Flux2)", 
                       filename="histograma_sum.png", csv_filename="sum_flux_data.csv",
                       excel_filename="sum_flux_data.xlsx", save_figure=True, show_fig=True, 
                       max_combinations=None):
    """
    Sum each vector from flux1 with each vector from flux2, classify resulting vectors,
    generate a histogram of process types and save results to CSV and Excel files.
    - CSV: Includes v_f1, v_f2, v_combined, S*v and Process_Type in automatic order.
    - Excel: Includes only v_combined, S*v and Process_Type, with sheets per Process_Type.

    Parameters:
    - S: Transformation matrix (numpy array).
    - flux1: DataFrame or NumPy array with Flux_r* columns and optionally Time (e.g., challenge vectors).
    - flux2: DataFrame or NumPy array with Flux_r* columns and optionally Time (e.g., cognitive control vectors).
    - title: Histogram title.
    - filename: File name to save the histogram.
    - csv_filename: CSV file name to save all data.
    - excel_filename: Excel file name to save v_combined, S*v and Process_Type.
    - save_figure: Boolean to save the figure.
    - show_fig: Boolean to display the figure.
    - max_combinations: Optional limit for the number of combinations to generate (default: None).

    Returns:
    - fig, ax, combined_df: Figure object, matplotlib axes and DataFrame with summed vectors.
    """
    # Convert flux1 and flux2 to DataFrame if they are NumPy arrays
    if isinstance(flux1, np.ndarray):
        flux1 = pd.DataFrame(flux1, columns=[f'Flux_r{i+1}' for i in range(flux1.shape[1])])
    if isinstance(flux2, np.ndarray):
        flux2 = pd.DataFrame(flux2, columns=[f'Flux_r{i+1}' for i in range(flux2.shape[1])])

    # Extract flux columns
    flux_columns = [col for col in flux1.columns if col.startswith('Flux_r')]
    if not all(col in flux2.columns for col in flux_columns):
        raise ValueError("flux1 and flux2 must have the same flux columns (Flux_r*).")

    print(f"Se encontraron {len(flux1)} vectores en flux1.")
    print(f"Se encontraron {len(flux2)} vectores en flux2.")

    # Verify that both DataFrames are not empty
    if flux1.empty or flux2.empty:
        raise ValueError("One or both DataFrames (flux1 or flux2) are empty.")

    # Check if 'Time' columns are present
    has_time_f1 = 'Time' in flux1.columns
    has_time_f2 = 'Time' in flux2.columns
    print(f"flux1 has 'Time' column: {has_time_f1}")
    print(f"flux2 has 'Time' column: {has_time_f2}")

    # Generate summed combinations of vectors
    flux_data = []
    combination_count = 0
    for idx_f2, row_f2 in flux2.iterrows():
        time_f2 = row_f2['Time'] if has_time_f2 else idx_f2
        v_f2 = row_f2[flux_columns].to_numpy()
        for idx_f1, row_f1 in flux1.iterrows():
            if max_combinations is not None and combination_count >= max_combinations:
                break
            time_f1 = row_f1['Time'] if has_time_f1 else idx_f1
            v_f1 = row_f1[flux_columns].to_numpy()
            # Sum the vectors
            v_combined = v_f2 + v_f1
            # Use time from flux2 or average if both have 'Time'
            time_value = time_f2 if not has_time_f1 else (time_f1 + time_f2) / 2 if has_time_f2 else time_f1
            # Create dictionary for the combined vector
            row_data = {}
            # Add v_combined
            for col, val in zip(flux_columns, v_combined):
                row_data[col] = val
            # Add v_f1
            for i, val in enumerate(v_f1):
                row_data[f'Flux1_r{i+1}'] = val
            # Add v_f2
            for i, val in enumerate(v_f2):
                row_data[f'Flux2_r{i+1}'] = val
            # Add Time if present
            if has_time_f1 or has_time_f2:
                row_data['Time'] = time_value
            flux_data.append(row_data)
            combination_count += 1
        if max_combinations is not None and combination_count >= max_combinations:
            break

    # Create DataFrame with combined vectors
    flux_vector = pd.DataFrame(flux_data)
    print(f"Total de vectores combinados generados: {len(flux_vector)}")

    # Classify processes in the combined flux_vector
    flux_values = flux_vector[flux_columns]  # Exclude Time column if exists
    process_types = []
    Sv_values = []
    for _, row in flux_values.iterrows():
        v = row.to_numpy()
        # Use classify_process logic
        Sv = S @ v
        Sv_values.append(Sv) 
        if np.all((0 < Sv) & (Sv <= 1e-8)):
            category = "Stationary Mode"
        elif np.all(Sv > 0) and np.any(Sv > 1e-8):
            category = "Cognitive Control"
        elif np.any(Sv < 0) and np.any(Sv > 0):
            category = "Challenge"
        elif np.any(Sv < 0):
            category = "Problem"
        else:
            category = "Other"
        process_types.append(category)

    # Calculate S*v
    Sv_expanded = pd.DataFrame(Sv_values, 
                               columns=[f"S*v_{i+1}" for i in range(S.shape[0])],
                               index=flux_vector.index)

    # Create combined DataFrame
    combined_df = pd.concat([flux_vector, Sv_expanded], axis=1)
    combined_df["Process_Type"] = process_types

    # Automatically determine column order for CSV
    csv_columns = []
    if has_time_f1 or has_time_f2:
        csv_columns.append('Time')
    csv_columns.extend([col for col in combined_df.columns if col.startswith('Flux1_r')])
    csv_columns.extend([col for col in combined_df.columns if col.startswith('Flux2_r')])
    csv_columns.extend([col for col in combined_df.columns if col.startswith('Flux_r')])
    csv_columns.extend([col for col in combined_df.columns if col.startswith('S*v_')])
    csv_columns.append('Process_Type')
    csv_available_columns = [col for col in csv_columns if col in combined_df.columns]
    csv_df = combined_df.reindex(columns=csv_available_columns)

    # Determine columns for Excel (v_combined, S*v, Process_Type)
    excel_columns = []
    if has_time_f1 or has_time_f2:
        excel_columns.append('Time')
    excel_columns.extend([col for col in combined_df.columns if col.startswith('Flux_r')])
    excel_columns.extend([col for col in combined_df.columns if col.startswith('S*v_')])
    excel_columns.append('Process_Type')
    excel_available_columns = [col for col in excel_columns if col in combined_df.columns]
    excel_df = combined_df.reindex(columns=excel_available_columns)

    # Count frequencies for histogram
    process_counts = Counter(process_types) 
    category_order = PROCESS_TYPE_COLOR_MAP.keys()
    labels = [cat for cat in category_order if cat in process_counts]
    counts = [process_counts[cat] for cat in labels]

    colors = [PROCESS_TYPE_COLOR_MAP.get(label, "grey") for label in labels]

    # Plot
    fig, ax = plt.subplots(figsize=(10, 6))
    bars = ax.bar(labels, counts, color=colors)
    for bar in bars:
        yval = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2.0, yval, int(yval),
                va='bottom', ha='center', fontweight='bold')

    ax.set_xlabel("Tipo de Proceso", fontsize=12)
    ax.set_ylabel("Frecuencia", fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.tick_params(axis="x", rotation=30)
    ax.grid(axis='y', linestyle='--', alpha=0.7)

    # Save to CSV
    out_dir = "./visualizations/process_classification/data_xslx"
    os.makedirs(out_dir, exist_ok=True)

    # Define output subfolder
    out_dir_sub = "./visualizations/process_classification/data_csv"
    os.makedirs(out_dir_sub, exist_ok=True)  # Create subfolder if it doesn't exist

    try:
        # Save main CSV file in subfolder
        filepath_csv = os.path.join(out_dir_sub, csv_filename)
        csv_df.to_csv(filepath_csv, index=False)
        print(f"Procesos clasificados guardados en CSV: {filepath_csv}")

        # Save separate CSV files by Process_Type in subfolder
        for category, group in csv_df.groupby("Process_Type"):
            # Ensure file name is valid (replace spaces and limit length)
            safe_category = category[:31].replace(' ', '_').replace('/', '_').replace('\\', '_')
            category_filepath = os.path.join(out_dir_sub, f"{csv_filename[:-4]}_{safe_category}.csv")
            group.to_csv(category_filepath, index=False)
            print(f"Procesos de tipo '{category}' guardados en CSV: {category_filepath}")

    except PermissionError as e:
        warnings.warn(f"Error de permisos al guardar el archivo CSV: {e}")
        print("No se pudo guardar el archivo CSV. Por favor, verifica los permisos de escritura en la carpeta.")
    except OSError as e:
        warnings.warn(f"Error del sistema al guardar el archivo CSV: {e}")
        print("No se pudo guardar el archivo CSV. Por favor, verifica el espacio en disco o la validez del nombre del archivo.")
    except Exception as e:
        warnings.warn(f"Error inesperado al guardar el archivo CSV: {e}")
        print("No se pudo guardar el archivo CSV. Ocurrió un error inesperado.")

    # Save to Excel
    filepath_excel = os.path.join(out_dir, excel_filename)
    with pd.ExcelWriter(filepath_excel, engine="openpyxl") as writer:
        excel_df.to_excel(writer, sheet_name="All_Processes", index=False)
        for category, group in excel_df.groupby("Process_Type"):
            group.to_excel(writer, sheet_name=category[:31], index=False)
    print(f"Procesos clasificados guardados en Excel: {filepath_excel}") 

    # Save figure
    out_dir_hist = "./visualizations/process_classification/histogram_flux_sum"
    os.makedirs(out_dir_hist, exist_ok=True)
    if save_figure:
        filepath_fig = os.path.join(out_dir_hist, filename)
        fig.savefig(filepath_fig, dpi=300, bbox_inches="tight")
        print(f"Histograma guardado en: {filepath_fig}")

    if show_fig:
        plt.show()

    return fig, ax, combined_df


def analyze_process_proportions_over_time(rn, S, rate_list, spec_vector, x0, t_span=(0, 200), n_steps=1001, n_cols = 4,
    window_sizes=[1, 2, 3, 4, 5], save_path="./visualizations/process_classification/"
):
    """
    Analyze the evolution of 'Cognitive Control' + 'Stationary Mode'
    in different simulation step windows, showing:
      - Subplots with histograms of process types.
      - Global plot of proportions.
      - Combined plots of the simulation with maximum total proportion.

    Parameters
    ----------
    rn : ReactionNetwork
        Reaction network object.
    S : np.ndarray
        Stoichiometric matrix of the system.
    rate_list : list
        List of reaction rates.
    spec_vector : list or np.ndarray
        Species vector.
    x0 : np.ndarray
        Initial conditions of the system.
    t_span : tuple, optional
        Time interval for the simulation (default (0, 200)).
    n_steps : int, optional
        Number of steps in the simulation (default 1001).
    n_cols : int, optional
        Number of columns in histogram subplots.
    window_sizes : list, optional
        List with window sizes (number of steps) to analyze.
    save_path : str, optional
        Base path to save results and figures.

    Returns
    -------
    results_df : pd.DataFrame
        DataFrame with proportions of each process type per window
        and the window corresponding to maximum 'Cognitive Domain'.
    """
    os.makedirs(save_path, exist_ok=True)
    results = []

    n_windows = len(window_sizes)
    n_rows = math.ceil(n_windows / n_cols)

    # ==================================================================
    # HISTOGRAM SUBPLOTS
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(12, 4 * n_rows), constrained_layout=True)
    axes = axes.flatten()

    # --- SIMULATION ---
    time_series, flux_series = simulation(
        rn, rate=rate_list, spec_vector=spec_vector, x0=x0,
        t_span=t_span, n_steps=n_steps
    )

    # --- ANALYZE EACH WINDOW SIZE ---
    for i, window_size in enumerate(window_sizes):

        # Special case: window size 1 → direct copy
        if window_size == 1:
            time_window = time_series.copy()
            flux_window = flux_series.copy()

        # Windows larger than 1 → moving sum
        else:
            # --- TIME SERIES ---
            time_numeric = time_series.drop(columns=['Time']).copy()
            time_rolling = (
                time_numeric
                .rolling(window=window_size, min_periods=window_size)
                .sum()
                .dropna()
                .reset_index(drop=True)
            )
            time_adjusted = time_series['Time'].iloc[window_size - 1:].reset_index(drop=True)
            time_window = pd.concat([time_adjusted, time_rolling], axis=1)

            # --- FLUX SERIES ---
            flux_numeric = flux_series.drop(columns=['Time']).copy()
            flux_rolling = (
                flux_numeric
                .rolling(window=window_size, min_periods=window_size)
                .sum()
                .dropna()
                .reset_index(drop=True)
            )
            flux_adjusted = flux_series['Time'].iloc[window_size - 1:].reset_index(drop=True)
            flux_window = pd.concat([flux_adjusted, flux_rolling], axis=1)

        # Classify and plot histograms of process types
        _, _, _, process_frequencies = plot_process_types_histogram(
            flux_window, S,
            title=f"Window={window_size}, n_steps={flux_window.shape[0]}",
            save_figure=False,
            ax=axes[i],
            show_fig=False
        ) 
        # Save all subplots in a single image 
        plt.savefig(os.path.join(save_path, "histogramas_subplots.png"), dpi=300, bbox_inches="tight") 

        # ==================================================================
        # Calculate proportions
        cognitive_control = process_frequencies.get("Cognitive Control", 0)
        stationary_mode = process_frequencies.get("Stationary Mode", 0)
        problem = process_frequencies.get("Problem", 0)
        total_cognitive_domain = cognitive_control + stationary_mode
        n=int(flux_window.shape[0])

        # Store results
        results.append({
            "Ventana": window_size,
            "Número de pasos": int(n),
            "CC": cognitive_control,
            "SM": stationary_mode,
            "PB": problem,
            "Total CD": total_cognitive_domain,
            "Proporción CC": cognitive_control / n,
            "Proporción SM": stationary_mode / n,
            "Proporción PB": problem / n,
            "Proporción Total": total_cognitive_domain / n
        })

    # Hide empty subplots
    for j in range(n_windows, n_rows * n_cols):
        axes[j].axis('off')
    plt.show()

    # ==================================================================
    # GLOBAL SUMMARY
    results_df = pd.DataFrame(results)
    max_row = results_df.select_dtypes(include=[np.number]).max()
    max_row["Ventana"] = "Máximo"
    results_df = pd.concat([results_df, pd.DataFrame([max_row])], ignore_index=True)

    idx_max = results_df["Proporción Total"].idxmax()
    max_row_total = results_df.loc[idx_max, "Proporción Total"]
    max_row_n_steps = results_df.loc[idx_max, "Número de pasos"]
    max_window_total = results_df.loc[idx_max, "Ventana"]

    print(f"\nMáximo Proporción Total Cognitive Domain: {max_row_total:.4f} (Ventana = {max_window_total}, n={int(max_row_n_steps)})\n")
    print("Resumen Global:")
    print(results_df)

    # ==================================================================
    # GLOBAL PLOT OF PROPORTIONS 
    results_numeric = results_df[results_df["Ventana"].apply(lambda x: isinstance(x, (int, float)))]
    plt.figure(figsize=(8, 5))
    plt.plot(results_numeric["Ventana"], 100 * results_numeric["Proporción PB"],
            marker='o', linestyle='-', color='red', label="Problem")    
    plt.plot(results_numeric["Ventana"], 100 * results_numeric["Proporción CC"],
            marker='o', linestyle='-', color='green', label="Cognitive Control")
    plt.plot(results_numeric["Ventana"], 100 * results_numeric["Proporción SM"],
            marker='o', linestyle='-', color='cyan', label="Stationary Mode")
    plt.plot(results_numeric["Ventana"], 100 * results_numeric["Proporción Total"],
            marker='o', linestyle='-', color='purple', label="Total Cognitive Domain")
    plt.plot(max_window_total, 100 * max_row_total, 'yo', label="Máximo Total", markersize=10)
    plt.xlabel("Tamaño de ventana (pasos)")
    plt.ylabel("Porcentaje (%)")
    plt.title("Cognitive Domain vs tamaño de ventana")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.savefig(os.path.join(save_path, "proportions_plots.png"), dpi=150, bbox_inches="tight")
    plt.show()

    # ==================================================================
    # FINAL SIMULATION WITH MAXIMUM COGNITIVE DOMAIN PROPORTION  
    window_size = max_window_total   

    # --- TIME SERIES ---
    time_numeric = time_series.drop(columns=['Time']).copy()
    time_rolling = (
        time_numeric
        .rolling(window=window_size, min_periods=window_size)
        .sum()
        .dropna()
        .reset_index(drop=True)
    )
    time_adjusted = time_series['Time'].iloc[window_size - 1:].reset_index(drop=True)
    time_window_max = pd.concat([time_adjusted, time_rolling], axis=1)

    # --- FLUX SERIES ---
    flux_numeric = flux_series.drop(columns=['Time']).copy()
    flux_rolling = (
        flux_numeric
        .rolling(window=window_size, min_periods=window_size)
        .sum()
        .dropna()
        .reset_index(drop=True)
    )
    flux_adjusted = flux_series['Time'].iloc[window_size - 1:].reset_index(drop=True)
    flux_window_max = pd.concat([flux_adjusted, flux_rolling], axis=1)

    # --- RESULTS ---
    print(f"Ventana máxima (window_size={window_size})")
    print("time_window_max =\n", time_window_max)
    print("flux_window_max =\n", flux_window_max)

    # ==================================================================
    # Combined plots: time series, fluxes and histogram
    fig, axes = plt.subplots(1, 3, figsize=(18, 4))
    plot_series_with_domain_intervals(time_window_max, flux_window_max, S,
                                      title=f"Serie de Tiempo - Ventana {max_window_total} pasos", save_figure=False, ax=axes[0])
    plot_flux_with_domain_intervals(flux_window_max, S,
                                    title=f"Flujos - Ventana {max_window_total} pasos", save_figure=False, ax=axes[1])
    plot_process_types_histogram(flux_window_max, S,
                                 title=f"Histograma de Tipos de Procesos (n={int(max_row_n_steps)})", save_figure=False, ax=axes[2])
    # General title
    fig.suptitle(f"Resumen de Análisis — Ventana Óptima de {max_window_total} pasos", fontsize=18)    
    plt.tight_layout()
    plt.savefig(os.path.join(save_path, "combined_plots.png"), dpi=150, bbox_inches="tight")
    plt.show()

    return results_df


def generate_rescaled_process_series(flux_vector, window_size=1):
    """
    Create rescaled (new timescale) process series from the simulation.
    
    Args:
        flux_vector: Flux vector array or DataFrame
        window_size: Window size for rescaling
        
    Returns:
        Rescaled process time series
    """
    df = pd.DataFrame(flux_vector)

    # If there's a time column already, use it
    if "Time" not in df.columns:
        df.insert(0, "Time", np.arange(len(df)))

    # Ensure Time is the first column (optional)
    df = df[["Time"] + [c for c in df.columns if c != "Time"]]

    return rescale_process_time_series(df, window_size=window_size)


def plot_process_classes_over_time(times, classifications):
    """
    Makes a single plot with:
    - x-axis = time
    - y-axis = process categories
    - red points = incomplete processes
    - blue points = complete processes
    
    Args:
        times: Array of time points
        classifications: List of classification tuples for each time point
    """
    # --- canonical order of categories for y-axis ---
    main_modes = [
        "Problem",
        "Challenge", 
        "Stationary Mode",
        "Overproduction Mode",
        "Other"
    ]

    # maps process classification to category
    def primary_class(cat_list):
        # first element is always the "main mode" in your system
        main = cat_list[0]
        if main not in main_modes:
            return "Other"
        return main

    # check complete vs incomplete
    def completeness(cat_list):
        return "Complete" if "Complete Process" in cat_list else "Incomplete"

    # Create the plot
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # Separate complete and incomplete processes
    complete_times = []
    complete_categories = []
    incomplete_times = []
    incomplete_categories = []
    
    for t, cat in zip(times, classifications):
        category = primary_class(cat)
        comp_status = completeness(cat)
        
        if comp_status == "Complete":
            complete_times.append(t)
            complete_categories.append(category)
        else:
            incomplete_times.append(t)
            incomplete_categories.append(category)
    
    # Map categories to y-axis positions
    category_to_y = {mode: i for i, mode in enumerate(main_modes)}
    
    # Plot incomplete processes (red)
    if incomplete_times:
        incomplete_y = [category_to_y[cat] for cat in incomplete_categories]
        ax.scatter(incomplete_times, incomplete_y, color='red', s=50, 
                  alpha=0.7, label='Incomplete', edgecolors='darkred', linewidth=0.5)
    
    # Plot complete processes (blue)  
    if complete_times:
        complete_y = [category_to_y[cat] for cat in complete_categories]
        ax.scatter(complete_times, complete_y, color='blue', s=50, 
                  alpha=0.7, label='Complete', edgecolors='darkblue', linewidth=0.5)
    
    # Format the plot
    ax.set_yticks(range(len(main_modes)))
    ax.set_yticklabels(main_modes)
    ax.set_ylabel('Process Type')
    ax.set_xlabel('Time')
    ax.set_title('Process Classification Over Time\n(Red = Incomplete, Blue = Complete)')
    
    # Add grid for better readability
    ax.grid(True, alpha=0.3, axis='both')
    ax.legend()
    
    # Set nicer y limits with some padding
    ax.set_ylim(-0.5, len(main_modes) - 0.5)
    
    # Format x-axis if there are time values
    if times.any():
        tmin, tmax = min(times), max(times)
        ax.set_xlim(tmin - 0.1*(tmax-tmin), tmax + 0.1*(tmax-tmin))
    
    plt.tight_layout()
    plt.show()


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
    
    # Process type classification visualization
    'PROCESS_TYPE_COLOR_MAP',
    'plot_process_types_histogram',
    'plot_cone_and_region_3d',
    'plot_series_with_domain_intervals',
    'plot_flux_with_domain_intervals',
    'histogram_flux_sum',
    'analyze_process_proportions_over_time',
    'generate_rescaled_process_series',
    'plot_process_classes_over_time',
    
    # Utilities
    'setup_publication_style',
    'get_category_color'
]