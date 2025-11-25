#!/usr/bin/env python3
"""
Script 5: Threshold Detection Testing with Visualizations

Tests scale-dependent threshold detection capabilities (NEW functionality).

Tests:
- Threshold detection methods (derivative, sign change, variance spike, balance quality)
- Threshold types (category-level, system-level, robustness)
- Synthetic processes with designed thresholds
- Realistic processes with emergent thresholds
- Detection method comparison

Visualizations:
- Thresholds marked on category dynamics
- Before/after feature comparisons
- Detection method comparison plots
- Feature trajectory plots with threshold markers
- Confidence distribution histograms
- Multi-resolution state space
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Add project root to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pyCOT.io.functions import read_txt
from pyCOT.simulations import simulation
from pyCOT.semantic_partition import define_semantic_categories
from pyCOT.process_structure import multi_scale_aggregation
from pyCOT.process_analyzer import (
    analyze_category_behavior,
    analyze_temporal_scale,
    detect_thresholds
)
from pyCOT.plot_process_analysis import (
    plot_thresholds_on_dynamics,
    plot_category_dynamics_across_scales
)

# Create output directory
output_dir = "./outputs/script_threshold_detection"
os.makedirs(output_dir, exist_ok=True)

print("=" * 80)
print("THRESHOLD DETECTION TEST SUITE")
print("=" * 80)

# ============================================================================
# SECTION 1: SETUP - CREATE SYNTHETIC PROCESS WITH DESIGNED THRESHOLD
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 1: Synthetic Process with Designed Threshold")
print("=" * 80)

# Create a simple 2-species, 3-reaction network
print("\nCreating synthetic network...")

# Simple network: X <-> Y with production
# r1: -> X (production)
# r2: X -> Y (conversion)
# r3: Y -> (degradation)

S_synthetic = np.array([
    [1, -1, 0],   # X: produced by r1, consumed by r2
    [0, 1, -1]    # Y: produced by r2, consumed by r3
])

species_synthetic = ['X', 'Y']
n_reactions_syn = 3

print(f"Stoichiometric matrix:\n{S_synthetic}")

# Define semantic categories
category_def_syn = {
    'input': ['X'],
    'output': ['Y']
}

semantic_partition_syn = define_semantic_categories(species_synthetic, category_def_syn)

# Create process series with designed threshold at scale 50
print("\nCreating process series with threshold at scale ~50...")

n_steps_syn = 200
process_series_syn = []

for i in range(n_steps_syn):
    v = np.zeros(n_reactions_syn)
    
    # Design: Before scale 50, r1 dominates (X overproduced)
    # After scale 50, r2-r3 dominate (Y overproduced, X depleted)
    if i < 100:
        v[0] = 1.0 + 0.1 * np.random.random()  # r1: X production
        v[1] = 0.2 + 0.05 * np.random.random()  # r2: X->Y
        v[2] = 0.1 + 0.05 * np.random.random()  # r3: Y degradation
    else:
        v[0] = 0.2 + 0.05 * np.random.random()  # r1: X production (reduced)
        v[1] = 1.0 + 0.1 * np.random.random()  # r2: X->Y (increased)
        v[2] = 0.8 + 0.1 * np.random.random()  # r3: Y degradation
    
    process_series_syn.append(v)

print(f"Generated {len(process_series_syn)} synthetic processes")

# ============================================================================
# SECTION 2: MULTI-SCALE ANALYSIS OF SYNTHETIC PROCESS
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 2: Multi-Scale Analysis (Synthetic)")
print("=" * 80)

print("\nRunning temporal scale analysis...")

window_sizes_syn = [1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

mode2_synthetic = analyze_temporal_scale(
    process_series=process_series_syn,
    S=S_synthetic,
    semantic_partition=semantic_partition_syn,
    window_sizes=window_sizes_syn,
    stride=1
)

print(f"Analyzed at {len(window_sizes_syn)} scales")

# ============================================================================
# SECTION 3: THRESHOLD DETECTION - DERIVATIVE METHOD
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 3: Threshold Detection - Derivative Method")
print("=" * 80)

print("\nDetecting thresholds using derivative method...")

thresholds_derivative = detect_thresholds(
    mode2_synthetic,
    feature_name='net_effect',
    method='derivative',
    threshold_sensitivity=0.5
)

print(f"\nFound {len(thresholds_derivative)} thresholds (derivative method):")
for thresh in thresholds_derivative:
    print(f"  Scale {thresh.scale_value}: {thresh.category} - {thresh.threshold_type}")
    print(f"    Feature change: {thresh.feature_before:.4f} -> {thresh.feature_after:.4f}")
    print(f"    Confidence: {thresh.confidence:.3f}")
    print(f"    Interpretation: {thresh.interpretation}")

# ============================================================================
# SECTION 4: THRESHOLD DETECTION - SIGN CHANGE METHOD
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 4: Threshold Detection - Sign Change Method")
print("=" * 80)

print("\nDetecting thresholds using sign change method...")

thresholds_sign = detect_thresholds(
    mode2_synthetic,
    feature_name='net_effect',
    method='sign_change',
    threshold_sensitivity=0.5
)

print(f"\nFound {len(thresholds_sign)} thresholds (sign change method):")
for thresh in thresholds_sign:
    print(f"  Scale {thresh.scale_value}: {thresh.category} - {thresh.threshold_type}")
    print(f"    Sign change: {np.sign(thresh.feature_before):.0f} -> {np.sign(thresh.feature_after):.0f}")

# ============================================================================
# SECTION 5: THRESHOLD DETECTION - BALANCE QUALITY METHOD
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 5: Threshold Detection - Balance Quality Method")
print("=" * 80)

print("\nDetecting thresholds using balance quality method...")

thresholds_balance = detect_thresholds(
    mode2_synthetic,
    feature_name='balance_quality',
    method='derivative',
    threshold_sensitivity=0.3
)

print(f"\nFound {len(thresholds_balance)} thresholds (balance quality method):")
for thresh in thresholds_balance:
    print(f"  Scale {thresh.scale_value}: {thresh.category} - {thresh.threshold_type}")

# ============================================================================
# SECTION 6: REALISTIC PROCESS - LOTKA-VOLTERRA
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 6: Realistic Process (Lotka-Volterra)")
print("=" * 80)

print("\nLoading Lotka-Volterra and running simulation...")
rn_lv = read_txt('Txt/Lotka_Volterra.txt')
S_lv = rn_lv.stoichiometry_matrix()

x0_lv = [2, 3]
rate_list_lv = 'mak'
spec_vector_lv = [[.5], [.8], [0.1]]

time_series_lv, flux_vector_lv = simulation(
    rn_lv, rate=rate_list_lv, spec_vector=spec_vector_lv,
    x0=x0_lv, t_span=(0, 100), n_steps=1001
)

process_series_lv = [flux_vector_lv.iloc[i, 1:].values for i in range(len(flux_vector_lv))]

print(f"Generated {len(process_series_lv)} processes from simulation")

# Define semantic partition
category_def_lv = {
    'prey': ['x'],
    'predator': ['y']
}
semantic_partition_lv = define_semantic_categories(['x', 'y'], category_def_lv)

# Multi-scale analysis
print("\nRunning multi-scale analysis...")
window_sizes_lv = [1, 5, 10, 20, 50, 100, 200]

mode2_lv = analyze_temporal_scale(
    process_series=process_series_lv[:500],
    S=S_lv,
    semantic_partition=semantic_partition_lv,
    window_sizes=window_sizes_lv,
    stride=1
)

# Detect thresholds
print("\nDetecting thresholds in realistic process...")
thresholds_lv = detect_thresholds(
    mode2_lv,
    feature_name='net_effect',
    method='derivative',
    threshold_sensitivity=0.5
)

print(f"\nFound {len(thresholds_lv)} thresholds in Lotka-Volterra:")
for thresh in thresholds_lv:
    print(f"  Scale {thresh.scale_value}: {thresh.category} - {thresh.threshold_type}")

# ============================================================================
# SECTION 7: VISUALIZATION - THRESHOLDS ON DYNAMICS (SYNTHETIC)
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 7: Visualizations - Thresholds on Dynamics (Synthetic)")
print("=" * 80)

print("\nGenerating threshold visualization on category dynamics...")

# Plot category dynamics with thresholds
# Plot category dynamics with thresholds
result = plot_category_dynamics_across_scales(
    mode2_synthetic,
    category='input',
    feature='net_effect',
    show_fig=False
)

# Handle return value - could be Figure or Axes
if isinstance(result, plt.Figure):
    fig = result
    ax = fig.axes[0]
elif isinstance(result, plt.Axes):
    ax = result
    fig = ax.figure
else:
    # Fallback
    fig = plt.gcf()
    ax = plt.gca()

# Add threshold markers
for thresh in thresholds_derivative:
    if thresh.category == 'input':
        ax.axvline(x=thresh.scale_value, color='red', linestyle='--', 
                  linewidth=2, alpha=0.7, label=f'Threshold (deriv)')

for thresh in thresholds_sign:
    if thresh.category == 'input':
        ax.axvline(x=thresh.scale_value, color='orange', linestyle=':', 
                  linewidth=2, alpha=0.7, label=f'Threshold (sign)')

ax.legend()
plt.tight_layout()

save_path = os.path.join(output_dir, "synthetic_thresholds_input.png")
plt.savefig(save_path, dpi=300, bbox_inches='tight')
print(f"Saved: {save_path}")
plt.close()

# ============================================================================
# SECTION 8: VISUALIZATION - DETECTION METHOD COMPARISON
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 8: Detection Method Comparison")
print("=" * 80)

print("\nGenerating detection method comparison...")

fig, axes = plt.subplots(2, 1, figsize=(12, 10))

# Extract feature trajectories
scales = window_sizes_syn
categories = list(semantic_partition_syn.categories)

for cat_idx, cat in enumerate(categories):
    ax = axes[cat_idx]
    
    # Extract net effect trajectory
    net_effects = []
    for ws in scales:
        result = mode2_synthetic['results_by_scale'][ws]
        if cat in result['category_statistics']:
            # Access nested dictionary structure: ['net_effect']['mean']
            net_effects.append(result['category_statistics'][cat]['net_effect']['mean'])
        else:
            net_effects.append(0)
    
    # Plot trajectory
    ax.plot(scales, net_effects, 'b-o', linewidth=2, label='Net Effect', markersize=6)
    
    # Mark thresholds from different methods
    for thresh in thresholds_derivative:
        if thresh.category == cat:
            ax.axvline(x=thresh.scale_value, color='red', linestyle='--', 
                      linewidth=2, alpha=0.7)
            ax.plot(thresh.scale_value, thresh.feature_before, 'rs', markersize=10)
    
    for thresh in thresholds_sign:
        if thresh.category == cat:
            ax.axvline(x=thresh.scale_value, color='orange', linestyle=':', 
                      linewidth=2, alpha=0.7)
    
    ax.set_xlabel('Window Size (Scale)')
    ax.set_ylabel('Mean Net Effect')
    ax.set_title(f'Category: {cat.upper()} - Threshold Detection', 
                fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    
    # Add legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color='blue', linewidth=2, label='Net Effect'),
        Line2D([0], [0], color='red', linestyle='--', linewidth=2, label='Derivative Method'),
        Line2D([0], [0], color='orange', linestyle=':', linewidth=2, label='Sign Change Method')
    ]
    ax.legend(handles=legend_elements)

plt.suptitle('Threshold Detection: Method Comparison', fontsize=14, fontweight='bold')
plt.tight_layout()

save_path = os.path.join(output_dir, "method_comparison.png")
plt.savefig(save_path, dpi=300, bbox_inches='tight')
print(f"Saved: {save_path}")
plt.close()

# ============================================================================
# SECTION 9: VISUALIZATION - CONFIDENCE DISTRIBUTION
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 9: Confidence Distribution")
print("=" * 80)

print("\nGenerating confidence distribution plot...")

# Collect all thresholds
all_thresholds = thresholds_derivative + thresholds_sign + thresholds_balance

if len(all_thresholds) > 0:
    confidences = [t.confidence for t in all_thresholds]
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    ax.hist(confidences, bins=20, color='steelblue', alpha=0.7, edgecolor='black')
    ax.set_xlabel('Confidence Score')
    ax.set_ylabel('Frequency')
    ax.set_title('Threshold Detection Confidence Distribution', fontsize=14, fontweight='bold')
    ax.grid(axis='y', alpha=0.3)
    
    # Add mean line
    mean_conf = np.mean(confidences)
    ax.axvline(x=mean_conf, color='red', linestyle='--', linewidth=2, 
              label=f'Mean: {mean_conf:.3f}')
    ax.legend()
    
    plt.tight_layout()
    
    save_path = os.path.join(output_dir, "confidence_distribution.png")
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {save_path}")
    plt.close()
else:
    print("No thresholds detected, skipping confidence plot")

# ============================================================================
# SECTION 10: VISUALIZATION - THRESHOLD DENSITY
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 10: Threshold Density Across Scales")
print("=" * 80)

print("\nGenerating threshold density plot...")

# Count thresholds at each scale
threshold_counts = {ws: 0 for ws in window_sizes_syn}
for thresh in all_thresholds:
    if thresh.scale_value in threshold_counts:
        threshold_counts[thresh.scale_value] += 1

fig, ax = plt.subplots(figsize=(12, 6))

scales_sorted = sorted(threshold_counts.keys())
counts = [threshold_counts[s] for s in scales_sorted]

ax.bar(scales_sorted, counts, color='coral', alpha=0.7, edgecolor='black')
ax.set_xlabel('Window Size (Scale)')
ax.set_ylabel('Number of Thresholds')
ax.set_title('Threshold Density: How Many Thresholds at Each Scale', 
            fontsize=14, fontweight='bold')
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()

save_path = os.path.join(output_dir, "threshold_density.png")
plt.savefig(save_path, dpi=300, bbox_inches='tight')
print(f"Saved: {save_path}")
plt.close()

# ============================================================================
# SECTION 11: VISUALIZATION - BEFORE/AFTER COMPARISON
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 11: Before/After Feature Comparison")
print("=" * 80)

print("\nGenerating before/after comparison...")

if len(thresholds_derivative) > 0:
    fig, axes = plt.subplots(1, len(thresholds_derivative), 
                            figsize=(5*len(thresholds_derivative), 5))
    
    if len(thresholds_derivative) == 1:
        axes = [axes]
    
    for idx, thresh in enumerate(thresholds_derivative):
        ax = axes[idx]
        
        # Bar chart showing before vs after
        features = ['Before', 'After']
        values = [thresh.feature_before, thresh.feature_after]
        colors_bars = ['lightcoral' if v < 0 else 'lightgreen' for v in values]
        
        ax.bar(features, values, color=colors_bars, alpha=0.7, edgecolor='black')
        ax.axhline(y=0, color='black', linestyle='-', linewidth=1)
        ax.set_ylabel('Net Effect')
        ax.set_title(f'{thresh.category.upper()}\nScale {thresh.scale_value}', 
                    fontsize=11, fontweight='bold')
        ax.grid(axis='y', alpha=0.3)
        
        # Add value labels
        for i, (feat, val) in enumerate(zip(features, values)):
            ax.text(i, val, f'{val:.3f}', ha='center', 
                   va='bottom' if val > 0 else 'top', fontweight='bold')
    
    plt.suptitle('Threshold Feature Values: Before vs After', 
                fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    save_path = os.path.join(output_dir, "before_after_comparison.png")
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {save_path}")
    plt.close()

# ============================================================================
# SECTION 12: VISUALIZATION - REALISTIC PROCESS THRESHOLDS
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 12: Realistic Process Thresholds (Lotka-Volterra)")
print("=" * 80)

print("\nGenerating Lotka-Volterra threshold visualization...")

# Plot both categories
fig, axes = plt.subplots(2, 1, figsize=(12, 10))

for cat_idx, cat in enumerate(['prey', 'predator']):
    ax = axes[cat_idx]
    
    # Extract net effect trajectory
    net_effects = []
    for ws in window_sizes_lv:
        result = mode2_lv['results_by_scale'][ws]
        if cat in result['category_statistics']:
            # Access nested dictionary structure: ['net_effect']['mean']
            net_effects.append(result['category_statistics'][cat]['net_effect']['mean'])
        else:
            net_effects.append(0)
    
    # Plot trajectory
    ax.plot(window_sizes_lv, net_effects, 'b-o', linewidth=2, markersize=6)
    
    # Mark thresholds
    for thresh in thresholds_lv:
        if thresh.category == cat:
            ax.axvline(x=thresh.scale_value, color='red', linestyle='--', 
                      linewidth=2, alpha=0.7)
            ax.plot(thresh.scale_value, thresh.feature_before, 'rs', markersize=10)
            ax.text(thresh.scale_value, thresh.feature_before, 
                   f'  {thresh.threshold_type}', rotation=90, 
                   va='bottom', fontsize=9)
    
    ax.set_xlabel('Window Size (Scale)')
    ax.set_ylabel('Mean Net Effect')
    ax.set_title(f'{cat.upper()} Category', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5)

plt.suptitle('Lotka-Volterra: Detected Thresholds', fontsize=14, fontweight='bold')
plt.tight_layout()

save_path = os.path.join(output_dir, "lotka_volterra_thresholds.png")
plt.savefig(save_path, dpi=300, bbox_inches='tight')
print(f"Saved: {save_path}")
plt.close()

# ============================================================================
# SECTION 13: SUMMARY
# ============================================================================
print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)

print(f"\nSynthetic Process Analysis:")
print(f"  Process length: {len(process_series_syn)}")
print(f"  Window sizes tested: {len(window_sizes_syn)}")
print(f"  Thresholds detected (derivative): {len(thresholds_derivative)}")
print(f"  Thresholds detected (sign change): {len(thresholds_sign)}")
print(f"  Thresholds detected (balance quality): {len(thresholds_balance)}")
print(f"  Total unique thresholds: {len(all_thresholds)}")

print(f"\nRealistic Process Analysis (Lotka-Volterra):")
print(f"  Process length: {len(process_series_lv)}")
print(f"  Window sizes tested: {len(window_sizes_lv)}")
print(f"  Thresholds detected: {len(thresholds_lv)}")

print(f"\nThreshold Types Identified:")
threshold_types = set(t.threshold_type for t in all_thresholds)
for ttype in threshold_types:
    count = sum(1 for t in all_thresholds if t.threshold_type == ttype)
    print(f"  {ttype}: {count}")

print(f"\nDetection Method Performance:")
print(f"  Derivative method: {len(thresholds_derivative)} thresholds")
print(f"  Sign change method: {len(thresholds_sign)} thresholds")
print(f"  Balance quality method: {len(thresholds_balance)} thresholds")

if len(all_thresholds) > 0:
    print(f"\nConfidence Statistics:")
    confidences = [t.confidence for t in all_thresholds]
    print(f"  Mean confidence: {np.mean(confidences):.3f}")
    print(f"  Min confidence: {np.min(confidences):.3f}")
    print(f"  Max confidence: {np.max(confidences):.3f}")

print(f"\nAll visualizations saved to: {output_dir}/")
print("  - synthetic_thresholds_input.png")
print("  - method_comparison.png")
print("  - confidence_distribution.png")
print("  - threshold_density.png")
print("  - before_after_comparison.png")
print("  - lotka_volterra_thresholds.png")

print("\n" + "=" * 80)
print("✅ SCRIPT COMPLETED SUCCESSFULLY")
print("=" * 80)
