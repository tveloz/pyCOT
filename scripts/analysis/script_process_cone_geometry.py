#!/usr/bin/env python3
"""
Script 2: Cone Geometry Testing with Visualizations

Tests process space geometry and cone structure analysis.
Recovers functionality from old script8.py with updated imports and enhanced visualizations.

Tests:
- Nullspace computation and properties
- Feasible region computation (Sv >= 0)
- Linear combination reconstruction
- Classification of feasible points
- Multiple networks with different nullspace dimensions

Visualizations:
- 3D cone projections with classification coloring
- 2D cone projections (multi-panel)
- Nullspace basis visualization
- Extra vectors (challenges, problems, counteractions)
- Feasible region density plots
- Classification pie charts
- Reconstruction accuracy plots
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d import Axes3D

# Add project root to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pyCOT.io.functions import read_txt
from pyCOT.process_structure import (
    compute_nullspace_vectors, 
    compute_feasible_region,
    classify_feasible_points,
    analyze_cone
)
from pyCOT.plot_process_analysis import plot_cone_and_region, plot_cone_2d_projections

# Create output directory
output_dir = "./outputs/script_cone_geometry"
os.makedirs(output_dir, exist_ok=True)
os.makedirs(os.path.join(output_dir, "cone_3d_projections"), exist_ok=True)

print("=" * 80)
print("CONE GEOMETRY TEST SUITE")
print("=" * 80)

# ============================================================================
# SECTION 1: LOAD NETWORKS
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 1: Loading Networks")
print("=" * 80)

# Load multiple networks
networks = {}

print("\nLoading Lotka-Volterra (3 reactions, 1D nullspace)...")
networks['LV'] = {
    'rn': read_txt('Txt/Lotka_Volterra.txt'),
    'name': 'Lotka-Volterra'
}
networks['LV']['S'] = networks['LV']['rn'].stoichiometry_matrix()
print(f"  Shape: {networks['LV']['S'].shape}")

print("\nLoading Autopoietic (5 reactions, 2D nullspace)...")
networks['Auto'] = {
    'rn': read_txt('Txt/autopoietic.txt'),
    'name': 'Autopoietic'
}
networks['Auto']['S'] = networks['Auto']['rn'].stoichiometry_matrix()
print(f"  Shape: {networks['Auto']['S'].shape}")

# ============================================================================
# SECTION 2: NULLSPACE ANALYSIS
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 2: Nullspace Analysis")
print("=" * 80)

for key, net in networks.items():
    print(f"\n{net['name']} Network:")
    
    # Compute nullspace
    null_vectors = compute_nullspace_vectors(net['S'])
    net['null_vectors'] = null_vectors
    
    print(f"  Nullspace dimension: {len(null_vectors)}")
    
    for i, vec in enumerate(null_vectors):
        print(f"  Basis vector {i+1}: {vec}")
        
        # Verify Sv = 0
        Sv = net['S'] @ vec
        norm = np.linalg.norm(Sv)
        print(f"    ||S @ v_{i+1}||: {norm:.6e} (should be ~0)")
        assert norm < 1e-6, f"Nullspace verification failed for {net['name']}"
    
    print(f"  ✓ Nullspace verification passed")

# ============================================================================
# SECTION 3: LINEAR COMBINATION RECONSTRUCTION
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 3: Linear Combination Reconstruction")
print("=" * 80)

# Test on Autopoietic network
S = networks['Auto']['S']
null_vectors = networks['Auto']['null_vectors']

if len(null_vectors) > 0:
    N = np.column_stack(null_vectors)
    
    # Test case 1: Simple combination
    print("\nTest 1: Simple linear combination")
    v_test = np.array([0.5, 0.25, 0.5, 0.25, 0.25])
    print(f"  Original vector: {v_test}")
    
    # Solve N @ c = v
    c, residuals, rank, s = np.linalg.lstsq(N, v_test, rcond=None)
    v_reconstructed = N @ c
    
    print(f"  Coefficients: {c}")
    print(f"  Reconstructed: {v_reconstructed}")
    print(f"  Reconstruction error: {np.linalg.norm(v_test - v_reconstructed):.6e}")
    
    # Test case 2: From old script8.py
    print("\nTest 2: Original script8.py test case")
    v_original = np.array([1., 0., 1., 1., 0.]) + np.array([1., 1., 1., 0., 1.])
    print(f"  Original vector: {v_original}")
    
    c2, _, _, _ = np.linalg.lstsq(N, v_original, rcond=None)
    v_reconstructed2 = N @ c2
    
    print(f"  Coefficients: {c2}")
    print(f"  Reconstructed: {v_reconstructed2}")
    print(f"  Reconstruction error: {np.linalg.norm(v_original - v_reconstructed2):.6e}")

# ============================================================================
# SECTION 4: FEASIBLE REGION COMPUTATION
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 4: Feasible Region Computation")
print("=" * 80)

for key, net in networks.items():
    print(f"\n{net['name']} Network:")
    
    # Compute feasible region with different resolutions
    for grid_res in [5, 10]:
        print(f"  Grid resolution: {grid_res}")
        
        feasible_points = compute_feasible_region(
            net['S'], 
            grid_max=None, 
            grid_res=grid_res, 
            auto_scale=True
        )
        
        print(f"    Feasible points found: {feasible_points.shape[0]}")
        
        # Verify all points satisfy Sv >= 0
        if feasible_points.shape[0] > 0:
            Sv_all = (net['S'] @ feasible_points.T).T
            all_feasible = np.all(Sv_all >= -1e-6, axis=1)
            print(f"    All points satisfy Sv >= 0: {np.all(all_feasible)}")
            
            if not np.all(all_feasible):
                print(f"    WARNING: {np.sum(~all_feasible)} points violate feasibility")
        
        net[f'feasible_{grid_res}'] = feasible_points

# ============================================================================
# SECTION 5: COMPREHENSIVE CONE ANALYSIS
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 5: Comprehensive Cone Analysis")
print("=" * 80)

for key, net in networks.items():
    print(f"\n{net['name']} Network:")
    
    # Run comprehensive analysis
    cone_data = analyze_cone(net['S'], grid_max=None, grid_res=15, classify=True)
    net['cone_data'] = cone_data
    
    print(f"  Nullspace dimension: {len(cone_data['nullspace_vectors'])}")
    print(f"  Feasible points: {cone_data['feasible_points'].shape[0]}")
    print(f"  Grid max used: {cone_data['grid_max']:.3f}")
    
    if 'classification_counts' in cone_data:
        print(f"  Classification breakdown:")
        for cls, count in cone_data['classification_counts'].items():
            print(f"    {cls}: {count}")

# ============================================================================
# SECTION 6: VISUALIZATIONS - 3D CONE PROJECTIONS
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 6: 3D Cone Projections")
print("=" * 80)

# Plot Autopoietic cone
print("\nGenerating 3D cone projections for Autopoietic network...")
S_auto = networks['Auto']['S']

files_3d, points = plot_cone_and_region(
    S_auto,
    grid_max=None,
    grid_res=15,
    show=False,
    save_dir=os.path.join(output_dir, "cone_3d_projections")
)

print(f"  Generated {len(files_3d)} 3D projection files")

# ============================================================================
# SECTION 7: VISUALIZATIONS - 2D CONE PROJECTIONS
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 7: 2D Cone Projections")
print("=" * 80)

# Plot 2D projections for Autopoietic
print("\nGenerating 2D cone projections for Autopoietic network...")

fig_2d = plot_cone_2d_projections(
    S_auto,
    grid_max=None,
    grid_res=25,
    figsize=(15, 10),
    save_path=os.path.join(output_dir, "cone_2d_projections_autopoietic.png"),
    show_fig=False
)

print(f"  Saved: {output_dir}/cone_2d_projections_autopoietic.png")

# Plot 2D projections for Lotka-Volterra
print("\nGenerating 2D cone projections for Lotka-Volterra network...")

fig_2d_lv = plot_cone_2d_projections(
    networks['LV']['S'],
    grid_max=None,
    grid_res=25,
    figsize=(12, 4),
    save_path=os.path.join(output_dir, "cone_2d_projections_lotka_volterra.png"),
    show_fig=False
)

print(f"  Saved: {output_dir}/cone_2d_projections_lotka_volterra.png")

# ============================================================================
# SECTION 8: VISUALIZATIONS - NULLSPACE BASIS
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 8: Nullspace Basis Visualization")
print("=" * 80)

print("\nGenerating nullspace basis visualization...")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Autopoietic nullspace
ax1 = axes[0]
null_vecs_auto = networks['Auto']['null_vectors']
n_reactions = networks['Auto']['S'].shape[1]

if len(null_vecs_auto) > 0:
    for i, vec in enumerate(null_vecs_auto):
        ax1.plot(range(n_reactions), vec, marker='o', label=f'Basis {i+1}', linewidth=2)
    
    ax1.set_xlabel('Reaction Index')
    ax1.set_ylabel('Coefficient Value')
    ax1.set_title('Autopoietic Nullspace Basis', fontsize=12, fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.axhline(y=0, color='black', linestyle='-', linewidth=0.5)

# Lotka-Volterra nullspace
ax2 = axes[1]
null_vecs_lv = networks['LV']['null_vectors']
n_reactions_lv = networks['LV']['S'].shape[1]

if len(null_vecs_lv) > 0:
    for i, vec in enumerate(null_vecs_lv):
        ax2.plot(range(n_reactions_lv), vec, marker='o', label=f'Basis {i+1}', linewidth=2)
    
    ax2.set_xlabel('Reaction Index')
    ax2.set_ylabel('Coefficient Value')
    ax2.set_title('Lotka-Volterra Nullspace Basis', fontsize=12, fontweight='bold')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.axhline(y=0, color='black', linestyle='-', linewidth=0.5)

plt.suptitle('Nullspace Basis Vectors', fontsize=14, fontweight='bold')
plt.tight_layout()

save_path = os.path.join(output_dir, "nullspace_basis.png")
plt.savefig(save_path, dpi=300, bbox_inches='tight')
print(f"Saved: {save_path}")
plt.close()

# ============================================================================
# SECTION 9: VISUALIZATIONS - CLASSIFICATION PIE CHARTS
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 9: Classification Distribution")
print("=" * 80)

print("\nGenerating classification pie charts...")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

color_map = {
    'Cognitive Domain': '#2ecc71',
    'Stationary Mode': '#3498db',
    'Overproduction Mode': '#1abc9c',
    'Problem': '#e74c3c',
    'Challenge': '#f39c12'
}

# Autopoietic
ax1 = axes[0]
if 'classification_counts' in networks['Auto']['cone_data']:
    counts = networks['Auto']['cone_data']['classification_counts']
    labels = list(counts.keys())
    sizes = list(counts.values())
    colors = [color_map.get(label, 'gray') for label in labels]
    
    ax1.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90)
    ax1.set_title('Autopoietic Network', fontsize=12, fontweight='bold')

# Lotka-Volterra
ax2 = axes[1]
if 'classification_counts' in networks['LV']['cone_data']:
    counts = networks['LV']['cone_data']['classification_counts']
    labels = list(counts.keys())
    sizes = list(counts.values())
    colors = [color_map.get(label, 'gray') for label in labels]
    
    ax2.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90)
    ax2.set_title('Lotka-Volterra Network', fontsize=12, fontweight='bold')

plt.suptitle('Process Mode Distribution in Feasible Region', fontsize=14, fontweight='bold')
plt.tight_layout()

save_path = os.path.join(output_dir, "classification_pie_charts.png")
plt.savefig(save_path, dpi=300, bbox_inches='tight')
print(f"Saved: {save_path}")
plt.close()

# ============================================================================
# SECTION 10: EXTRA VECTORS (from original script8.py)
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 10: Extra Vectors Visualization")
print("=" * 80)

print("\nGenerating cone with extra vectors (challenges, problems, etc.)...")

# Define extra vectors for Autopoietic network (from original script8.py)
vc = np.array([2, 2, 1, 2, 1])  # Challenge
vp = np.array([0, 0, 1, 0, 1])  # Problem
v1 = np.array([2, 0, 2.5, 0, 0])  # Counteraction
v_combined = vc + v1  # Combined

extra_vectors = [vc, vp, v1, v_combined]
extra_labels = ["Challenge (vc)", "Problem (vp)", "Counteraction (v1)", "Combined (vc+v1)"]
extra_colors = ['orange', 'red', 'purple', 'brown']

# Verify stoichiometric effects
print("\nExtra vector stoichiometric effects:")
for label, vec in zip(extra_labels, extra_vectors):
    Sv = S_auto @ vec
    print(f"  {label}: Sv = {Sv}")

# Plot with extra vectors
files_extra, _ = plot_cone_and_region(
    S_auto,
    grid_max=None,
    grid_res=12,
    show=False,
    extra_vectors=extra_vectors,
    extra_vector_labels=extra_labels,
    extra_vector_colors=extra_colors,
    save_dir=os.path.join(output_dir, "cone_3d_with_extra_vectors")
)

print(f"  Generated {len(files_extra)} 3D projections with extra vectors")

# ============================================================================
# SECTION 11: RECONSTRUCTION ACCURACY VISUALIZATION
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 11: Reconstruction Accuracy")
print("=" * 80)

print("\nGenerating reconstruction accuracy plots...")

# Test multiple random vectors
np.random.seed(42)
n_test_vectors = 20
n_reactions = S_auto.shape[1]
N_auto = np.column_stack(networks['Auto']['null_vectors'])

test_vectors = []
reconstruction_errors = []

for i in range(n_test_vectors):
    v_test = np.random.uniform(0, 2, n_reactions)
    c, _, _, _ = np.linalg.lstsq(N_auto, v_test, rcond=None)
    v_reconstructed = N_auto @ c
    error = np.linalg.norm(v_test - v_reconstructed)
    
    test_vectors.append(v_test)
    reconstruction_errors.append(error)

# Plot reconstruction errors
fig, ax = plt.subplots(figsize=(10, 6))

ax.bar(range(n_test_vectors), reconstruction_errors, color='steelblue', alpha=0.7, edgecolor='black')
ax.set_xlabel('Test Vector Index')
ax.set_ylabel('Reconstruction Error ||v - v_reconstructed||')
ax.set_title('Nullspace Reconstruction Accuracy', fontsize=14, fontweight='bold')
ax.grid(axis='y', alpha=0.3)

# Add horizontal line at tolerance
ax.axhline(y=1e-6, color='red', linestyle='--', linewidth=2, label='Tolerance (1e-6)')
ax.legend()

plt.tight_layout()

save_path = os.path.join(output_dir, "reconstruction_accuracy.png")
plt.savefig(save_path, dpi=300, bbox_inches='tight')
print(f"Saved: {save_path}")
plt.close()

print(f"\nReconstruction statistics:")
print(f"  Mean error: {np.mean(reconstruction_errors):.6e}")
print(f"  Max error: {np.max(reconstruction_errors):.6e}")
print(f"  All errors < 1e-6: {np.all(np.array(reconstruction_errors) < 1e-6)}")

# ============================================================================
# SECTION 12: SUMMARY
# ============================================================================
print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)

print(f"\nNetworks analyzed: {len(networks)}")
for key, net in networks.items():
    print(f"\n{net['name']}:")
    print(f"  Nullspace dimension: {len(net['null_vectors'])}")
    if 'cone_data' in net:
        print(f"  Feasible points: {net['cone_data']['feasible_points'].shape[0]}")
        if 'classification_counts' in net['cone_data']:
            total = sum(net['cone_data']['classification_counts'].values())
            print(f"  Total classified: {total}")

print(f"\nAll visualizations saved to: {output_dir}/")
print("  - cone_3d_projections/ (directory)")
print("  - cone_3d_with_extra_vectors/ (directory)")
print("  - cone_2d_projections_autopoietic.png")
print("  - cone_2d_projections_lotka_volterra.png")
print("  - nullspace_basis.png")
print("  - classification_pie_charts.png")
print("  - reconstruction_accuracy.png")

print("\n" + "=" * 80)
print("✅ SCRIPT COMPLETED SUCCESSFULLY")
print("=" * 80)
