# Script 6: Spatial (PDE) and Metapopulation Simulations (UPDATED FOR NEW STRUCTURE)

# Import necessary libraries and modules
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from pyCOT.io.functions import read_txt
# UPDATED IMPORTS for new structure
from pyCOT.simulations.ode import simulation
from pyCOT.simulations.spatial import simulate_spatial_dynamics
from pyCOT.simulations.metapopulation import simulate_metapopulation_dynamics
from pyCOT.plot_dynamics import * 

# ################################################################################################# 
# # Example 1: Spatial (PDE) Reaction-Diffusion on 2D Grid
# ################################################################################################# 
# file_path = 'Txt/autopoietic.txt'  
# rn = read_txt(file_path)
#
# species = [s.name for s in rn.species()]
# print("Species:", species)
# reactions = [r.name() for r in rn.reactions()]
# print("Reactions:", reactions) 
#
# ##########################################################################
# # Spatial diffusion dynamics in 2D (UPDATED FUNCTION NAME)
# ##########################################################################
# rate = 'mak'
# grid_shape = (3, 3) 
#
# # Simulate spatial reaction-diffusion dynamics
# # NOTE: Function renamed from simulate_diffusion_dynamics_2D to simulate_spatial_dynamics
# t_spatial, X_spatial, flux_spatial = simulate_spatial_dynamics(
#     rn, 
#     rate=rate, 
#     grid_shape=grid_shape, 
#     t_span=(0, 10),
#     n_steps=100
# )
#
# # Visualization 
# plot_diffusion_time_series_2D(
#     time=t_spatial,
#     concentration_data=X_spatial, 
#     grid_shape=grid_shape, 
#     xlabel='Time', 
#     ylabel='Concentration', 
#     main_title='Spatial Concentration Time Series', 
#     legend_title='Species', 
#     cell_prefix='Cell'
# )
#
# plot_diffusion_time_series_2D(
#     time=t_spatial,
#     concentration_data=flux_spatial, 
#     grid_shape=grid_shape, 
#     xlabel='Time', 
#     ylabel='Flux', 
#     main_title='Spatial Flux Time Series', 
#     legend_title='Reactions', 
#     cell_prefix='Cell'
# )
#
# plot_heatmaps_all_species_2D(
#     t_spatial, 
#     X_spatial, 
#     time_indices=3, 
#     main_title="Evolution of Spatial Concentration Profiles"
# )
#
# animate_diffusion_heatmaps_all_species_2D(
#     t_spatial, 
#     X_spatial, 
#     main_title="Animation of Spatial Dynamics"
# )

# ##########################################################################
# # Example 2: Simple Metapopulation Dynamics
# ##########################################################################
# file_path = 'Txt/autopoietic.txt'  
# rn = read_txt(file_path)
#
# rate = 'mak'
# num_patches = 3  # Use num_patches instead of grid_shape
#
# # Simulate metapopulation dynamics
# t_MP, X_MP, flux_MP = simulate_metapopulation_dynamics(
#     rn, 
#     rate=rate, 
#     num_patches=num_patches, 
#     t_span=(0, 10),
#     n_steps=100
# )
#
# # Visualization  
# # Note: For vector-based metapopulation, we need vector-specific plots
# plot_diffusion_time_series_2D(
#     time=t_MP,
#     concentration_data=X_MP, 
#     grid_shape=(num_patches, 1),  # Treat as 1D array of patches
#     xlabel='Time', 
#     ylabel='Concentration', 
#     main_title='Metapopulation Concentration Time Series', 
#     legend_title='Species', 
#     cell_prefix='Patch'
# )
#
# plot_diffusion_time_series_2D(
#     time=t_MP,
#     concentration_data=flux_MP, 
#     grid_shape=(num_patches, 1), 
#     xlabel='Time', 
#     ylabel='Flux', 
#     main_title='Metapopulation Flux Time Series', 
#     legend_title='Reactions', 
#     cell_prefix='Patch'
# )

##########################################################################
# FULL EXAMPLE: Conflict Model with Metapopulation Dynamics
##########################################################################
print("="*80)
print("CONFLICT MODEL: METAPOPULATION DYNAMICS")
print("="*80)

file_path = 'Txt/cause_driven_conflict_gov.txt'
rn = read_txt(file_path)

species = [s.name for s in rn.species()]
print(f"\nSpecies ({len(species)}): {species}")
reactions = [r.name() for r in rn.reactions()]
print(f"Reactions ({len(reactions)}): {reactions}")

# === Kinetic rate specification ===
# Mix of mass action (mak) and Michaelis-Menten (mmk)
rates = [
    "mak",  # r1  A_p + G_A => A_p + G_A + iG_A
    "mak",  # r2  B_p + G_B => B_p + G_B + iG_B
    "mak",  # r3  A_v + D_A => A_v + D_A + iD_A
    "mak",  # r4  B_v + D_B => B_v + D_B + iD_B
    "mmk",  # r5  A_v + A_p + iD_A => 2A_v + iD_A
    "mmk",  # r6  B_v + B_p + iD_B => 2B_v + iD_B
    "mmk",  # r7  A_v + A_p + iG_A => 2A_p + iG_A
    "mmk",  # r8  B_v + B_p + iG_B => 2B_p + iG_B
    "mak",  # r9  A_v + G_B + iD_A => A_v + D_B
    "mak",  # r10 B_v + G_A + iD_B => B_v + D_A
    "mak",  # r11 A_p + D_A + iG_A => A_p + G_A
    "mak",  # r12 B_p + D_B => B_p + G_B
    "mmk",  # r13 2iG_A => (decay)
    "mmk",  # r14 2iG_B =>
    "mmk",  # r15 2iD_A =>
    "mmk",  # r16 2iD_B =>
    "mak",  # r17 B_p + D_B => B_p + D_B + iD_B
    "mak",  # r18 B_p + iD_B => B_v
    "mak",  # r19 A_p + D_A => A_v + D_A + iD_A
    "mak",  # r20 A_p + iD_A => A_v
    "mmk",  # r21 Gov_A + 2iD_A + D_A => Gov_A + G_A
    "mak",  # r22 Gov_A + G_A => Gov_A + G_A + iG_A
    "mmk",  # r23 Gov_B + 2iD_B + D_B => Gov_B + G_B
    "mak",  # r24 Gov_B + G_B => Gov_B + G_B + iG_B
]

# === Reaction parameters ===
spec_vector = [
    [0.01],          # r1
    [0.01],          # r2
    [0.01],          # r3
    [0.01],          # r4
    [0.1, 0.01],     # r5  [Vmax, Km]
    [0.1, 0.01],     # r6  [Vmax, Km]
    [0.1, 0.01],     # r7  [Vmax, Km]
    [0.1, 0.01],     # r8  [Vmax, Km]
    [0.01],          # r9
    [0.01],          # r10
    [0.01],          # r11
    [0.01],          # r12
    [0.001, 0.001],  # r13 [Vmax, Km]
    [0.001, 0.001],  # r14 [Vmax, Km]
    [0.001, 0.001],  # r15 [Vmax, Km]
    [0.001, 0.001],  # r16 [Vmax, Km]
    [0.01],          # r17
    [0],             # r18 (inactive)
    [0],             # r19 (inactive)
    [0],             # r20 (inactive)
    [0.0001, 0.0001], # r21 [Vmax, Km]
    [0.001],         # r22
    [0.0001, 0.0001], # r23 [Vmax, Km]
    [0.001],         # r24
]

# === Dispersal coefficients (set to 0 for isolated populations) ===
D_dict = {
    "A_p": 0.0, "A_v": 0.0,
    "B_p": 0.0, "B_v": 0.0,
    "G_A": 0.0, "G_B": 0.0,
    "D_A": 0.0, "D_B": 0.0,
    "iG_A": 0.0, "iG_B": 0.0,
    "iD_A": 0.0, "iD_B": 0.0,
    "Gov_A": 0.0, "Gov_B": 0.0,
}

# === Initial conditions for 2 populations ===
x0_dict = {
    "A_p":   np.array([0.5, 0.0]),
    "A_v":   np.array([0.05, 0.0]),
    "B_p":   np.array([0.55, 0.0]),
    "B_v":   np.array([0.0, 0.0]),
    "G_A":   np.array([0.3, 0.0]),
    "G_B":   np.array([0.32, 0.0]),
    "D_A":   np.array([0.2, 0.0]),
    "D_B":   np.array([0.22, 0.0]),
    "iG_A":  np.array([0.1, 0]),
    "iG_B":  np.array([0.12, 0.0]),
    "iD_A":  np.array([0.1, 0.0]),
    "iD_B":  np.array([0.12, 0.0]),
    "Gov_A": np.array([0.6, 0]),
    "Gov_B": np.array([0.6, 0.0]),
}

# === Connectivity matrix (isolated populations, no exchange) ===
connectivity_matrix = np.array([
    [1.0, 0.0],  # Population A: stays in place
    [0.0, 1.0],  # Population B: stays in place
])

num_patches = connectivity_matrix.shape[0]

print(f"\nSimulating {num_patches} populations with metapopulation dynamics...")

# === Run metapopulation simulation (UPDATED IMPORT) ===
t_MP, X_MP, flux_MP = simulate_metapopulation_dynamics(
    rn, 
    rate=rates,
    num_patches=num_patches,
    D_dict=D_dict, 
    x0_dict=x0_dict, 
    spec_vector=spec_vector, 
    t_span=(0, 100), 
    n_steps=200, 
    connectivity_matrix=connectivity_matrix
)

print("Simulation complete. Generating visualizations...")

# === Visualization ===
# Group species for clearer plots - SIMPLIFIED
species_pairs = [
    ('A_p', 'A_v'),  # Group A populations
    ('G_A', 'D_A'),  # Group A resources
]

# Alternative: plot individual species without pairs
# This avoids the subplot calculation issue
print("\nPlotting individual species time series...")

import matplotlib.pyplot as plt

fig, axes = plt.subplots(len(species_pairs), num_patches, 
                         figsize=(12, 8), sharex=True)

for row_idx, pair in enumerate(species_pairs):
    for patch_idx in range(num_patches):
        ax = axes[row_idx, patch_idx] if num_patches > 1 else axes[row_idx]
        
        for species_name in pair:
            ax.plot(t_MP, X_MP[species_name][:, patch_idx], 
                   label=species_name, linewidth=2)
        
        ax.set_ylabel('Concentration')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        if row_idx == 0:
            ax.set_title(f'Patch {patch_idx}')
        if row_idx == len(species_pairs) - 1:
            ax.set_xlabel('Time')

plt.tight_layout()
plt.savefig("metapopulation_concentration_conflict.png", dpi=150)
plt.show()

print("Concentration plot saved.")

# For flux, use simpler approach
print("\nPlotting reaction fluxes...")

flux_to_plot = ["r1", "r5", "r9", "r11"]  # Selected reactions

fig2, axes2 = plt.subplots(len(flux_to_plot), num_patches,
                           figsize=(12, 10), sharex=True)

for row_idx, reaction in enumerate(flux_to_plot):
    for patch_idx in range(num_patches):
        ax = axes2[row_idx, patch_idx] if num_patches > 1 else axes2[row_idx]
        
        ax.plot(t_MP, flux_MP[reaction][:, patch_idx],
               label=reaction, linewidth=2, color='red')
        
        ax.set_ylabel('Flux')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        if row_idx == 0:
            ax.set_title(f'Patch {patch_idx}')
        if row_idx == len(flux_to_plot) - 1:
            ax.set_xlabel('Time')

plt.tight_layout()
plt.savefig("metapopulation_flux_conflict.png", dpi=150)
plt.show()

print("Flux plot saved.")

print("\n" + "="*80)
print("SIMULATION COMPLETED SUCCESSFULLY")
print("="*80)
print("\nGenerated files:")
print("  - metapopulation_concentration_conflict.png")
print("  - metapopulation_flux_conflict.png")