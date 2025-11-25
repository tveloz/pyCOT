#!/usr/bin/env python3
"""
Script 4: Multi-Scale Modes Testing with Visualizations

Tests the three analytical modes and multi-scale aggregation capabilities.
Replaces old script7.py with modular approach using new architecture.

Tests:
- Mode 1: Pathway Variability (one process, many decompositions)
- Mode 2: Temporal Scale (many processes, rolling windows)
- Mode 3: Full Multiscale Robustness (combined analysis)
- Aggregation/Disaggregation functions

Visualizations:
- Mode 1: Pathway variability plots, decomposition ensembles
- Mode 2: Category dynamics across scales, balance heatmaps
- Mode 3: Three-mode comparison, z-score heatmaps, robustness indicators
- Comparative plots: process proportion trends, multi-scale histograms
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import time

# Add project root to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pyCOT.io.functions import read_txt
from pyCOT.simulations import simulation
from pyCOT.semantic_partition import define_semantic_categories
from pyCOT.process_structure import (
    aggregate_processes,
    rolling_window_aggregation,
    multi_scale_aggregation,
    disaggregate_process,
    disaggregation_random_sampling,
)
from pyCOT.process_analyzer import (
    classify_process_mode,
    analyze_category_behavior,
    analyze_pathway_variability,
    analyze_temporal_scale,
    analyze_full_multiscale_robustness
)
from pyCOT.plot_process_analysis import (
    plot_category_dynamics_across_scales,
    plot_all_categories_dynamics,
    plot_category_balance_heatmap,
    plot_pattern_robustness,
    plot_decomposition_ensemble,
    plot_three_mode_comparison
)

# Create output directory
output_dir = "./outputs/script_multiscale_modes"
os.makedirs(output_dir, exist_ok=True)

print("=" * 80)
print("MULTI-SCALE MODES TEST SUITE")
print("=" * 80)

# ============================================================================
# SECTION 1: SETUP - LOAD NETWORK AND RUN SIMULATION
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 1: Setup and Simulation")
print("=" * 80)

# Load Lotka-Volterra for testing (simpler than conflict model)
# print("\nLoading Lotka-Volterra network...")
# file_path = 'Txt/Lotka_Volterra.txt'
# rn = read_txt(file_path)

# species_list = [s.name for s in rn.species()]
# S = rn.stoichiometry_matrix()
# n_reactions = len(rn.reactions())

# print(f"Species: {species_list}")
# print(f"Reactions: {n_reactions}")

# # Define semantic categories (peace vs predator)
# category_definitions = {
#     'peace': ['R'],
#     'predator': ['F']
# }
file_path = 'networks/Conflict_Theory/cause_driven_conflict_gov.txt'
rn = read_txt(file_path)

species_list = [s.name for s in rn.species()]
S = rn.stoichiometry_matrix()
n_reactions = len(rn.reactions())

print(f"\nLoaded model: {len(species_list)} species, {n_reactions} reactions")
print(f"Species: {species_list[:10]}..." if len(species_list) > 10 else f"Species: {species_list}")

# Create semantic partition
category_definitions = {
    'conflict': ['A_v', 'B_v', 'D_A', 'D_B', 'iD_A', 'iD_B'],
    'peace': ['A_p', 'B_p', 'G_A', 'G_B', 'iG_A', 'iG_B']
}

semantic_partition = define_semantic_categories(species_list, category_definitions)
print(f"\nSemantic partition: {semantic_partition.categories}")

# Run simulation
print("\nRunning simulation...")
#x0 = [2, 3]
rate_list = 'mak'
#spec_vector = [[.5], [.8], [0.1]]
tf = 100
n_steps = 1000

time_series, flux_vector = simulation(
    rn, rate=rate_list, 
#    spec_vector=spec_vector,x0=x0, 
    t_span=(0, tf), n_steps=n_steps+1
)

print(f"Simulation complete: {len(flux_vector)} time steps")

# Extract process vectors
process_series = [flux_vector.iloc[i, 1:].values for i in range(len(flux_vector))]
print(f"Process series extracted: {len(process_series)} processes")

# ============================================================================
# SECTION 2: AGGREGATION/DISAGGREGATION TESTS
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 2: Aggregation/Disaggregation Tests")
print("=" * 80)

# Test simple aggregation
print("\nTest 1: Simple Aggregation")
v1 = process_series[0]
v2 = process_series[1]
v3 = process_series[2]

agg_simple = aggregate_processes([v1, v2, v3])
print(f"  Aggregated process: {agg_simple.aggregated_process}")
print(f"  Sum verification: {agg_simple.verify_sum()}")

# Test rolling window aggregation
print("\nTest 2: Rolling Window Aggregation")
window_size = 10
stride = 5

rolling_aggs = rolling_window_aggregation(process_series, window_size, stride)
print(f"  Window size: {window_size}")
print(f"  Stride: {stride}")
print(f"  Number of windows: {len(rolling_aggs)}")
print(f"  First window verification: {rolling_aggs[0].verify_sum()}")

# Test multi-scale aggregation
print("\nTest 3: Multi-Scale Aggregation")
window_sizes = [1, 10, 50]

multi_scale_aggs = multi_scale_aggregation(process_series[:200], window_sizes)
print(f"  Window sizes: {window_sizes}")
for ws, aggs in multi_scale_aggs.items():
    print(f"    Window {ws}: {len(aggs)} aggregations")

# Test disaggregation
print("\nTest 4: Disaggregation Random Sampling")
v_aggregate = agg_simple.aggregated_process
n_samples = 5
n_components = 3

disagg_samples = disaggregation_random_sampling(
    v_aggregate, n_samples, n_components, strategy='mixed', seed=42
)

print(f"  Original process: {v_aggregate}")
print(f"  Generated {len(disagg_samples)} decomposition samples")
print(f"  Each with {n_components} components")
print(f"  Sample 1 verification: {disagg_samples[0].verify_sum()}")

# ============================================================================
# SECTION 3: MODE 1 - PATHWAY VARIABILITY
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 3: Mode 1 - Pathway Variability Analysis")
print("=" * 80)

# Select an aggregate process from simulation
print("\nAnalyzing pathway variability...")
t_start = time.time()

# Use aggregation of processes 100-120
v_aggregate_mode1 = aggregate_processes(process_series[100:120]).aggregated_process

mode1_results = analyze_pathway_variability(
    process=v_aggregate_mode1,
    S=S,
    semantic_partition=semantic_partition,
    n_samples=50,
    n_components=10,
    strategy='mixed',
    seed=42
)

t_mode1 = time.time() - t_start
print(f"Mode 1 analysis completed in {t_mode1:.2f}s")

print(f"\nCategory statistics:")
for cat, stats in mode1_results['category_statistics'].items():
    print(f"  {cat.upper()}:")
    print(f"    Mean net effect: {stats['net_effect']['mean']:.4f}")
    print(f"    Std net effect: {stats['net_effect']['std']:.4f}")
    print(f"    Pathway variance: {stats['pathway_variance']:.4f}")
    print(f"    Range: [{stats['net_effect']['min']:.4f}, {stats['net_effect']['max']:.4f}]")
    print(f"    Balance quality mean: {stats['balance_quality']['mean']:.4f}")
# ============================================================================
# SECTION 4: MODE 2 - TEMPORAL SCALE
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 4: Mode 2 - Temporal Scale Analysis")
print("=" * 80)

print("\nAnalyzing temporal scale...")
t_start = time.time()

# Use subset of process series for efficiency
process_subset = process_series[:500]
window_sizes = [1, 5, 10, 20, 50, 100]

mode2_results = analyze_temporal_scale(
    process_series=process_subset,
    S=S,
    semantic_partition=semantic_partition,
    window_sizes=window_sizes,
    stride=1
)

t_mode2 = time.time() - t_start
print(f"Mode 2 analysis completed in {t_mode2:.2f}s")

print(f"\nResults by scale:")
for ws in window_sizes:
    result = mode2_results['results_by_scale'][ws]
    print(f"\n  Window size {ws}:")
    print(f"    Number of aggregations: {result['aggregations']}")
    for cat, stats in result['category_statistics'].items():
        print(f"    Mean net effect: {stats['net_effect']['mean']:.4f}")

# ============================================================================
# SECTION 5: MODE 3 - FULL MULTISCALE ROBUSTNESS
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 5: Mode 3 - Full Multiscale Robustness Analysis")
print("=" * 80)

print("\nAnalyzing full multiscale robustness...")
print("WARNING: This is computationally intensive...")
t_start = time.time()

# Use even smaller subset for Mode 3
process_subset_mode3 = process_series[:200]
window_sizes_mode3 = [1, 10, 50]

mode3_results = analyze_full_multiscale_robustness(
    process_series=process_subset_mode3,
    S=S,
    semantic_partition=semantic_partition,
    window_sizes=window_sizes_mode3,
    n_decomp_samples=20,
    n_decomp_components=5,
    stride=5,
    seed=42
)

t_mode3 = time.time() - t_start
print(f"Mode 3 analysis completed in {t_mode3:.2f}s")

print(f"\nRobustness analysis results:")
for ws in window_sizes_mode3:
    results_at_scale = mode3_results['results_by_scale'][ws]
    print(f"\n  Window size {ws}:")
    print(f"    Type of results_at_scale: {type(results_at_scale)}")
    print(f"    Length: {len(results_at_scale) if hasattr(results_at_scale, '__len__') else 'N/A'}")
    
    # Try to access data based on actual type
    if isinstance(results_at_scale, dict):
        print(f"    Keys: {list(results_at_scale.keys())}")
        first_key = list(results_at_scale.keys())[0]
        result = results_at_scale[first_key]
    elif isinstance(results_at_scale, list):
        result = results_at_scale[0]
    else:
        print(f"    Unexpected type")
        continue
    
    print(f"    Type of result: {type(result)}")
    if isinstance(result, dict):
        print(f"    Result keys: {list(result.keys())}")
    elif isinstance(result, list):
        print(f"    Result is a list with {len(result)} items")
        if len(result) > 0:
            print(f"    Type of result[0]: {type(result[0])}")
            if isinstance(result[0], dict):
                print(f"    Result[0] keys: {list(result[0].keys())}")

# ============================================================================
# SECTION 6: VISUALIZATION - MODE 1 PATHWAY VARIABILITY
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 6: Mode 1 Visualizations")
print("=" * 80)

print("\nGenerating Mode 1 visualizations...")

# Pathway variance bar chart
fig, ax = plt.subplots(figsize=(10, 6))

categories = list(mode1_results['category_statistics'].keys())
variances = [mode1_results['category_statistics'][cat]['pathway_variance'] 
            for cat in categories]

colors = ['red' if cat == 'predator' else 'green' for cat in categories]
ax.bar(categories, variances, color=colors, alpha=0.7, edgecolor='black')

ax.set_xlabel('Category')
ax.set_ylabel('Pathway Variance')
ax.set_title('Mode 1: Pathway Variability by Category', fontsize=14, fontweight='bold')
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
save_path = os.path.join(output_dir, "mode1_pathway_variance.png")
plt.savefig(save_path, dpi=300, bbox_inches='tight')
print(f"Saved: {save_path}")
plt.close()

# Decomposition ensemble plot
print("Generating decomposition ensemble plot...")

fig = plot_decomposition_ensemble(
    mode1_results,
    category='peace',
    save_path=os.path.join(output_dir, "mode1_decomposition_ensemble_peace.png"),
    show_fig=False
)
print(f"Saved: {output_dir}/mode1_decomposition_ensemble_peace.png")

# ============================================================================
# SECTION 7: VISUALIZATION - MODE 2 TEMPORAL SCALE
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 7: Mode 2 Visualizations")
print("=" * 80)

print("\nGenerating Mode 2 visualizations...")

# Category dynamics across scales
fig = plot_category_dynamics_across_scales(
    mode2_results,
    category='peace',
    feature='net_effect',
    save_path=os.path.join(output_dir, "mode2_peace_dynamics.png"),
    show_fig=False
)
print(f"Saved: {output_dir}/mode2_peace_dynamics.png")

# All categories dynamics
fig = plot_all_categories_dynamics(
    mode2_results,
    feature='net_effect',
    save_path=os.path.join(output_dir, "mode2_all_categories.png"),
    show_fig=False
)
print(f"Saved: {output_dir}/mode2_all_categories.png")

# Balance heatmap
fig = plot_category_balance_heatmap(
    mode2_results,
    save_path=os.path.join(output_dir, "mode2_balance_heatmap.png"),
    show_fig=False
)
print(f"Saved: {output_dir}/mode2_balance_heatmap.png")

# ============================================================================
# SECTION 8: VISUALIZATION - MODE 3 ROBUSTNESS
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 8: Mode 3 Visualizations")
print("=" * 80)

print("\nGenerating Mode 3 visualizations...")

# Three-mode comparison
fig = plot_three_mode_comparison(
    mode1_results=mode1_results,
    mode2_results=mode2_results,
    mode3_results=mode3_results,
    category='peace',
    save_path=os.path.join(output_dir, "mode3_three_mode_comparison.png"),
    show_fig=False
)
print(f"Saved: {output_dir}/mode3_three_mode_comparison.png")

# Z-score heatmap
print("Generating z-score heatmap...")

fig, ax = plt.subplots(figsize=(10, 6))

# Extract z-scores - handle mode3_results structure
z_scores_data = {}

# Check if results_by_scale exists and has data
if 'results_by_scale' in mode3_results and mode3_results['results_by_scale']:
    for ws in window_sizes_mode3:
        if ws not in mode3_results['results_by_scale']:
            continue
            
        results_at_scale = mode3_results['results_by_scale'][ws]
        
        # Initialize scale z-scores
        scale_z_scores = {cat: [] for cat in semantic_partition.categories}
        
        # Handle different possible structures
        if isinstance(results_at_scale, list):
            for result in results_at_scale:
                if isinstance(result, dict) and 'temporal_vs_decomposition' in result:
                    comparison = result['temporal_vs_decomposition']
                    for cat in semantic_partition.categories:
                        if cat in comparison and 'z_score' in comparison[cat]:
                            scale_z_scores[cat].append(comparison[cat]['z_score'])
        
        # Store mean z-score for this scale
        for cat in semantic_partition.categories:
            if cat not in z_scores_data:
                z_scores_data[cat] = []
            if len(scale_z_scores[cat]) > 0:
                z_scores_data[cat].append(np.mean(scale_z_scores[cat]))
            else:
                z_scores_data[cat].append(0)

# Create heatmap only if we have data
if z_scores_data:
    categories_ordered = sorted(z_scores_data.keys())
    z_matrix = np.array([z_scores_data[cat] for cat in categories_ordered])
else:
    print("Warning: No z-score data available for heatmap")
    categories_ordered = []
    z_matrix = np.array([])

im = ax.imshow(z_matrix, aspect='auto', cmap='RdYlGn_r', vmin=-3, vmax=3)

ax.set_xticks(range(len(window_sizes_mode3)))
ax.set_xticklabels(window_sizes_mode3)
ax.set_yticks(range(len(categories_ordered)))
ax.set_yticklabels(categories_ordered)

ax.set_xlabel('Window Size')
ax.set_ylabel('Category')
ax.set_title('Mode 3: Z-Scores (Temporal vs Decomposition)', fontsize=14, fontweight='bold')

cbar = plt.colorbar(im, ax=ax)
cbar.set_label('Z-Score', rotation=270, labelpad=20)

plt.tight_layout()
save_path = os.path.join(output_dir, "mode3_zscore_heatmap.png")
plt.savefig(save_path, dpi=300, bbox_inches='tight')
print(f"Saved: {save_path}")
plt.close()

# ============================================================================
# SECTION 9: COMPARATIVE VISUALIZATION - REPLACING OLD SCRIPT7
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 9: Comparative Visualizations (replacing old script7.py)")
print("=" * 80)

print("\nGenerating process proportion trends...")

# Extract Cognitive Domain proportions at each scale from Mode 2
# Analyze proportion of Cognitive Domain-like processes at each scale
cognitive_proportions = []
for ws in window_sizes:
    result = mode2_results['results_by_scale'][ws]
    aggregations = result['aggregations']
    n_total = len(aggregations)
    n_cognitive = 0
    
    # This would require storing classification info in Mode 2
    # For now, approximate based on category behavior
    # A process is "cognitive domain-like" if all categories have net_effect >= 0
    
    for agg in aggregations[:10]:  # Sample first 10
        behavior = analyze_category_behavior(agg.aggregated_process, S, semantic_partition)
        if all(behavior[cat].net_effect >= 0 for cat in semantic_partition.categories):
            n_cognitive += 1
    
    proportion = n_cognitive / min(n_total, 10) if n_total > 0 else 0
    cognitive_proportions.append(proportion)  # MISSING - need to append!

ax.plot(window_sizes, cognitive_proportions, marker='o', linewidth=2,
       markersize=8, color='green', label='Cognitive Domain-like')

ax.set_xlabel('Window Size')
ax.set_ylabel('Proportion')
ax.set_title('Process Classification Proportions vs Window Size', fontsize=14, fontweight='bold')
ax.grid(True, alpha=0.3)
ax.legend()
ax.set_ylim([0, 1])

plt.tight_layout()
save_path = os.path.join(output_dir, "proportion_trends.png")
plt.savefig(save_path, dpi=300, bbox_inches='tight')
print(f"Saved: {save_path}")
plt.close()

# ============================================================================
# SECTION 10: PERFORMANCE BENCHMARKING
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 10: Performance Benchmarking")
print("=" * 80)

print("\nExecution times:")
print(f"  Mode 1 (50 samples, 10 components): {t_mode1:.2f}s")
print(f"  Mode 2 (500 processes, 6 scales): {t_mode2:.2f}s")
print(f"  Mode 3 (200 processes, 3 scales, 20 samples): {t_mode3:.2f}s")

# Create performance comparison plot
fig, ax = plt.subplots(figsize=(10, 6))

modes = ['Mode 1\nPathway\nVariability', 'Mode 2\nTemporal\nScale', 'Mode 3\nFull\nRobustness']
times = [t_mode1, t_mode2, t_mode3]
colors = ['steelblue', 'coral', 'mediumpurple']

bars = ax.bar(modes, times, color=colors, alpha=0.7, edgecolor='black')

# Add value labels on bars
for bar, time_val in zip(bars, times):
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2., height,
           f'{time_val:.2f}s',
           ha='center', va='bottom', fontweight='bold')

ax.set_ylabel('Execution Time (seconds)')
ax.set_title('Execution Time Comparison Across Modes', fontsize=14, fontweight='bold')
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
save_path = os.path.join(output_dir, "performance_comparison.png")
plt.savefig(save_path, dpi=300, bbox_inches='tight')
print(f"Saved: {save_path}")
plt.close()

# ============================================================================
# SECTION 11: SUMMARY
# ============================================================================
print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)

print(f"\nSimulation: {len(process_series)} time steps")
print(f"Semantic categories: {semantic_partition.categories}")

print(f"\nMode 1 (Pathway Variability):")
print(f"  Aggregate process analyzed with {mode1_results['n_samples']} decomposition samples")
print(f"  Each sample has {mode1_results['n_components']} components")

print(f"\nMode 2 (Temporal Scale):")
print(f"  Analyzed {len(process_subset)} processes")
print(f"  Window sizes: {window_sizes}")
print(f"  Total aggregations: {sum(len(mode2_results['results_by_scale'][ws]['aggregations']) for ws in window_sizes)}")

print(f"\nMode 3 (Full Multiscale Robustness):")
print(f"  Analyzed {len(process_subset_mode3)} processes")
print(f"  Window sizes: {window_sizes_mode3}")
print(f"  Decomposition samples per temporal aggregation: 20")

print(f"\nAll visualizations saved to: {output_dir}/")
print("  Mode 1:")
print("    - mode1_pathway_variance.png")
print("    - mode1_decomposition_ensemble_peace.png")
print("  Mode 2:")
print("    - mode2_peace_dynamics.png")
print("    - mode2_all_categories.png")
print("    - mode2_balance_heatmap.png")
print("  Mode 3:")
print("    - mode3_three_mode_comparison.png")
print("    - mode3_zscore_heatmap.png")
print("  Comparative:")
print("    - proportion_trends.png")
print("    - performance_comparison.png")

print("\n" + "=" * 80)
print("✅ SCRIPT COMPLETED SUCCESSFULLY")
print("=" * 80)
