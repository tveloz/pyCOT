# Script: Spatial Dynamics Testing (Reaction-Diffusion on 2D Grid)
# Tests the new simulate_spatial_dynamics function with various models

# ========================================
# 1. LIBRARY LOADING AND CONFIGURATION
# ========================================
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from pyCOT.io.functions import read_txt
from pyCOT.simulations.spatial import simulate_spatial_dynamics
from pyCOT.plot_dynamics import (
    plot_diffusion_time_series_2D,
    plot_heatmaps_all_species_2D,
    animate_diffusion_heatmaps_all_species_2D
)

print("="*80)
print("SPATIAL DYNAMICS TEST SCRIPT")
print("Reaction-Diffusion Simulations on 2D Grids")
print("="*80)

# ========================================
# 2. MODEL SELECTION
# ========================================
# Choose one of the following models:

# --- Simple Models (Good for testing) ---
file_path = 'Txt/autopoietic.txt'  # Autopoietic network (default)
# file_path = 'networks/testing/Lotka_Volterra.txt'  # Predator-prey

# --- Complex Models ---
# file_path = 'networks/Conflict_Theory/conflict_toy_model0.txt'  # Conflict model
# file_path = 'Txt/BZ_cycle.txt'  # Belousov-Zhabotinsky
# file_path = 'Txt/SEIR.txt'  # Epidemic model

rn = read_txt(file_path)

# Get species and reactions
species = [s.name for s in rn.species()]
reactions = [r.name() for r in rn.reactions()]

print(f"\nLoaded model: {file_path}")
print(f"Species ({len(species)}): {species}")
print(f"Reactions ({len(reactions)}): {reactions}")

# ========================================
# 3. SPATIAL PARAMETERS CONFIGURATION
# ========================================
print("\n" + "="*80)
print("CONFIGURING SPATIAL PARAMETERS")
print("="*80)

# --- Grid Configuration ---
grid_shape = (5, 5)  # 5x5 grid (25 cells total)
# grid_shape = (10, 10)  # 10x10 grid for higher resolution
# grid_shape = (3, 3)   # 3x3 grid for quick testing

rows, cols = grid_shape
print(f"\nGrid shape: {rows} × {cols} = {rows*cols} cells")

# --- Kinetic Laws ---
# Option 1: Single kinetic for all reactions
rate = 'mak'  # Mass action kinetics

# Option 2: Mixed kinetics (uncomment if needed)
# rate = ['mak', 'mmk', 'hill', ...]  # One per reaction

print(f"Kinetic laws: {rate}")

# --- Diffusion Coefficients ---
# Option 1: Let system generate random diffusion coefficients
D_dict = None  # Auto-generates uniform [0.01, 0.2]

# Option 2: Specify diffusion coefficients manually
# Useful for testing different diffusion rates
# D_dict = {
#     species[0]: 0.1,  # Fast diffusion
#     species[1]: 0.05, # Medium diffusion
#     species[2]: 0.01, # Slow diffusion
# }

if D_dict is None:
    # Generate random diffusion coefficients
    np.random.seed(42)
    D_dict = {sp: np.round(np.random.uniform(0.01, 0.2), 3) for sp in species}

print(f"\nDiffusion coefficients:")
for sp, D in D_dict.items():
    print(f"  {sp}: D = {D}")

# --- Initial Spatial Distribution ---
# Option 1: Random uniform initial conditions
x0_dict = None  # Auto-generates uniform [0, 2.0]

# Option 2: Localized initial conditions (spot in center)
# x0_dict = {}
# for sp in species:
#     x0 = np.zeros(grid_shape)
#     center_i, center_j = rows // 2, cols // 2
#     x0[center_i, center_j] = 1.0  # Concentrated spot
#     x0_dict[sp] = x0

# Option 3: Gradient initial conditions
# x0_dict = {}
# for sp_idx, sp in enumerate(species):
#     x0 = np.zeros(grid_shape)
#     for i in range(rows):
#         for j in range(cols):
#             x0[i, j] = (i + j) / (rows + cols)  # Linear gradient
#     x0_dict[sp] = x0

if x0_dict is None:
    # Generate random initial conditions
    np.random.seed(42)
    x0_dict = {sp: np.round(np.random.uniform(0, 2.0, size=grid_shape), 2) 
               for sp in species}

print(f"\nInitial conditions: Random spatial distribution")
print(f"Sample ({species[0]}):")
print(x0_dict[species[0]])

# --- Reaction Parameters ---
# Option 1: Auto-generate random parameters
spec_vector = None

# Option 2: Manual parameters (example for simple model)
# For autopoietic with MAK kinetics, all reactions need [k] parameter
# spec_vector = [[0.5] for _ in range(len(reactions))]  # All k=0.5

# Option 3: Custom parameters per reaction
# spec_vector = [
#     [0.8],   # Reaction 1: k=0.8
#     [0.3],   # Reaction 2: k=0.3
#     [1.0],   # Reaction 3: k=1.0
#     # ... one per reaction
# ]

if spec_vector is None:
    print(f"\nReaction parameters: Auto-generated")
else:
    print(f"\nReaction parameters: {spec_vector}")

# --- Simulation Time Parameters ---
t_span = (0, 20)      # Time interval [0, 20]
n_steps = 100         # 100 time points

print(f"\nTime span: {t_span}")
print(f"Time steps: {n_steps}")

# ========================================
# 4. RUN SPATIAL SIMULATION
# ========================================
print("\n" + "="*80)
print("RUNNING SPATIAL SIMULATION")
print("="*80)

print("\nSimulating reaction-diffusion dynamics...")
print("This may take a moment for large grids...")

t_spatial, X_spatial, flux_spatial = simulate_spatial_dynamics(
    rn,
    rate=rate,
    grid_shape=grid_shape,
    D_dict=D_dict,
    x0_dict=x0_dict,
    spec_vector=spec_vector,
    t_span=t_span,
    n_steps=n_steps,
    method='RK45',   # Integration method
    rtol=1e-6,       # Relative tolerance
    atol=1e-8        # Absolute tolerance
)

print(f"✓ Simulation complete!")
print(f"  Time points: {len(t_spatial)}")
print(f"  Species tracked: {len(X_spatial)}")
print(f"  Reactions tracked: {len(flux_spatial)}")

# Check data shapes
print(f"\nData shapes:")
for sp in species[:3]:  # Show first 3 species
    print(f"  {sp}: {X_spatial[sp].shape} (time × rows × cols)")

# ========================================
# 5. VISUALIZATION: TIME SERIES
# ========================================
print("\n" + "="*80)
print("GENERATING VISUALIZATIONS")
print("="*80)

print("\n1. Plotting concentration time series for each grid cell...")

plot_diffusion_time_series_2D(
    time=t_spatial,
    concentration_data=X_spatial,
    grid_shape=grid_shape,
    xlabel='Time',
    ylabel='Concentration',
    main_title='Spatial Concentration Time Series (All Grid Cells)',
    legend_title='Species',
    cell_prefix='Cell',
    filename='spatial_concentration_timeseries.png',
    show=True
)

print("✓ Saved: spatial_concentration_timeseries.png")

print("\n2. Plotting reaction flux time series...")

plot_diffusion_time_series_2D(
    time=t_spatial,
    concentration_data=flux_spatial,
    grid_shape=grid_shape,
    xlabel='Time',
    ylabel='Reaction Flux',
    main_title='Spatial Reaction Flux Time Series',
    legend_title='Reactions',
    cell_prefix='Cell',
    save_path='spatial_flux_timeseries.png',
    show=True
)

print("✓ Saved: spatial_flux_timeseries.png")

# ========================================
# 6. VISUALIZATION: SPATIAL HEATMAPS
# ========================================
print("\n3. Creating spatial heatmap snapshots at different times...")

# Select time indices to visualize
time_indices = [0, n_steps//4, n_steps//2, 3*n_steps//4, n_steps-1]
time_labels = ['Initial', 'T=25%', 'T=50%', 'T=75%', 'Final']

print(f"   Time snapshots: {len(time_indices)} frames")
print(f"   Times: {[f'{t_spatial[i]:.2f}' for i in time_indices]}")

plot_heatmaps_all_species_2D(
    t_spatial,
    X_spatial,
    time_indices=time_indices,
    main_title="Spatial Evolution of Concentration Profiles",
    filename='spatial_heatmaps_snapshots.png',
    show=True
)

print("✓ Saved: spatial_heatmaps_snapshots.png")

# ========================================
# 7. VISUALIZATION: ANIMATED HEATMAPS
# ========================================
print("\n4. Creating animated heatmap (this may take a moment)...")

animate_diffusion_heatmaps_all_species_2D(
    t_spatial,
    X_spatial,
    main_title="Reaction-Diffusion Dynamics Animation",
    filename='spatial_dynamics_animation.gif',
    interval=100,  # milliseconds between frames
    show=True
)

print("✓ Saved: spatial_dynamics_animation.gif")

# ========================================
# 8. ANALYSIS: SPATIAL STATISTICS
# ========================================
print("\n" + "="*80)
print("SPATIAL STATISTICS SUMMARY")
print("="*80)

print("\nFinal spatial distribution statistics:")

for sp in species:
    final_state = X_spatial[sp][-1, :, :]  # Last time point
    
    print(f"\n{sp}:")
    print(f"  Mean concentration: {final_state.mean():.4f}")
    print(f"  Std deviation: {final_state.std():.4f}")
    print(f"  Min value: {final_state.min():.4f}")
    print(f"  Max value: {final_state.max():.4f}")
    print(f"  Spatial variance: {final_state.var():.4f}")

# ========================================
# 9. COMPARISON: INITIAL VS FINAL
# ========================================
print("\n" + "="*80)
print("INITIAL VS FINAL COMPARISON")
print("="*80)

import matplotlib.pyplot as plt

# Create comparison plot for first species
sp_to_plot = species[0]
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Initial state
im1 = axes[0].imshow(X_spatial[sp_to_plot][0, :, :], cmap='viridis', 
                     origin='lower', interpolation='nearest')
axes[0].set_title(f'{sp_to_plot} - Initial (t={t_spatial[0]:.2f})')
axes[0].set_xlabel('X')
axes[0].set_ylabel('Y')
plt.colorbar(im1, ax=axes[0], label='Concentration')

# Final state
im2 = axes[1].imshow(X_spatial[sp_to_plot][-1, :, :], cmap='viridis',
                     origin='lower', interpolation='nearest')
axes[1].set_title(f'{sp_to_plot} - Final (t={t_spatial[-1]:.2f})')
axes[1].set_xlabel('X')
axes[1].set_ylabel('Y')
plt.colorbar(im2, ax=axes[1], label='Concentration')

plt.tight_layout()
plt.savefig('spatial_comparison_initial_final.png', dpi=150)
plt.show()

print("✓ Saved: spatial_comparison_initial_final.png")

# ========================================
# 10. SUMMARY
# ========================================
print("\n" + "="*80)
print("SIMULATION COMPLETE")
print("="*80)

print(f"\nModel: {file_path}")
print(f"Grid: {grid_shape[0]} × {grid_shape[1]} = {grid_shape[0]*grid_shape[1]} cells")
print(f"Time: {t_span[0]} to {t_span[1]} ({n_steps} steps)")
print(f"Species: {len(species)}")
print(f"Reactions: {len(reactions)}")

print("\nGenerated files:")
print("  1. spatial_concentration_timeseries.png")
print("  2. spatial_flux_timeseries.png")
print("  3. spatial_heatmaps_snapshots.png")
print("  4. spatial_dynamics_animation.gif")
print("  5. spatial_comparison_initial_final.png")

print("\n" + "="*80)
print("All visualizations generated successfully!")
print("="*80)

# ========================================
# OPTIONAL: TEST DIFFERENT SCENARIOS
# ========================================

# Uncomment to test different scenarios:

# # Test 1: High diffusion vs low diffusion
# print("\n\nTEST: Comparing high vs low diffusion...")
# D_high = {sp: 0.5 for sp in species}
# D_low = {sp: 0.01 for sp in species}
# 
# t_high, X_high, _ = simulate_spatial_dynamics(rn, rate=rate, grid_shape=(5,5), 
#                                                D_dict=D_high, t_span=(0,10), n_steps=50)
# t_low, X_low, _ = simulate_spatial_dynamics(rn, rate=rate, grid_shape=(5,5),
#                                              D_dict=D_low, t_span=(0,10), n_steps=50)
# 
# # Compare final spatial variance
# sp_test = species[0]
# var_high = X_high[sp_test][-1].var()
# var_low = X_low[sp_test][-1].var()
# print(f"Spatial variance (high D): {var_high:.4f}")
# print(f"Spatial variance (low D): {var_low:.4f}")
# print(f"Diffusion smoothed by factor: {var_low/var_high:.2f}x")

# # Test 2: Localized vs distributed initial conditions
# print("\n\nTEST: Localized source spreading...")
# x0_localized = {}
# for sp in species:
#     x0 = np.zeros((5, 5))
#     x0[2, 2] = 10.0  # Single concentrated spot
#     x0_localized[sp] = x0
# 
# t_loc, X_loc, _ = simulate_spatial_dynamics(rn, rate=rate, grid_shape=(5,5),
#                                              x0_dict=x0_localized, t_span=(0,10), n_steps=50)
# 
# plot_heatmaps_all_species_2D(t_loc, X_loc, time_indices=[0, 10, 25, 49],
#                              main_title="Localized Source Spreading")