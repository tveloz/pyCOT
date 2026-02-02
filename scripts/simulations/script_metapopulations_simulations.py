# Script 6: Diffusion and Metapopulations of Reaction Networks with pyCOT 

# Import necessary libraries and modules
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from pyCOT.io.functions import read_txt
from pyCOT.simulations.metapopulation import *
from pyCOT.visualization.plot_dynamics import * 

################################################################################################# 
# Ejemplo de uso  
file_path = 'Txt/Farm.txt'  
rn = read_txt(file_path)

species = [s.name for s in rn.species()]
print("Especies:", species)
reactions = [r.name() for r in rn.reactions()]
print("Reacciones:", reactions) 

# ##########################################################################
# # Diffusion dynamics in 2D
# ##########################################################################
# rate='mak'
# grid_shape = (3, 3) 

# # Simulación de la dinámica de difusión en 2D
# t_diff, X_diff, flux_diff = simulate_diffusion_dynamics_2D(rn, rate, grid_shape=grid_shape, t_span=(0,10)) 

# # Visualization 
# plot_diffusion_time_series_2D(time=t_diff,concentration_data= X_diff, grid_shape= grid_shape, xlabel='Tiempo', ylabel='Concentración', main_title='Serie de tiempo de las Concentraciones por Difusión', legend_title='Especies', cell_prefix='G')
# plot_diffusion_time_series_2D(time=t_diff,concentration_data=flux_diff, grid_shape=grid_shape, xlabel='Tiempo', ylabel='Concentración', main_title='Serie de tiempo del Flujo por Difusión', legend_title='Flujos', cell_prefix='G')

# plot_heatmaps_all_species_2D(t_diff, X_diff, time_indices=3, main_title="Evolución de los perfiles de concentración por Difusión")
# animate_diffusion_heatmaps_all_species_2D(t_diff, X_diff, main_title="Animación de la Dinámica de Difusión") 

# ##########################################################################
# # Metapopulation dynamics simulation 
# ##########################################################################
# rate = 'mak'
# grid_shape = (3, 3)
# t_MP, X_MP, flux_MP = simulate_metapopulation_dynamics(rn, rate, grid_shape = grid_shape, t_span=(0,10))

# # Visualization  
# plot_diffusion_time_series_2D(time=t_MP,concentration_data=X_MP, grid_shape=grid_shape, xlabel='Tiempo', ylabel='Concentración', main_title='Serie de tiempo de las Concentraciones de Metapoblaciones', legend_title='Especies', cell_prefix='G')
# plot_diffusion_time_series_2D(time=t_MP,concentration_data=flux_MP, grid_shape=grid_shape, xlabel='Tiempo', ylabel='Concentración', main_title='Serie de tiempo del Flujo de Metapoblaciones', legend_title='Flujos', cell_prefix='G')

# plot_heatmaps_all_species_2D(t_MP, X_MP, time_indices=3, main_title="Evolución de los perfiles de concentración de Metapoblaciones")
# animate_diffusion_heatmaps_all_species_2D(t_MP, X_MP, main_title="Animación de la Dinámica de Metapoblaciones")

##########################################################################
# EJEMPLO: Metapopulation dynamics with custom parameters for autopoietic
##########################################################################
##########################################################################
# EJEMPLO: Metapopulation dynamics with custom parameters for autopoietic (VECTOR-BASED)
##########################################################################
file_path = 'networks/Conflict_Theory/cause_driven_conflict_gov.txt'
rn = read_txt(file_path)
species = [s.name for s in rn.species()]
print("Especies:", species)
reactions = [r.name() for r in rn.reactions()]
print("Reacciones:", reactions)

num_patches = 2  # Número de poblaciones/parches

# Parámetros de reacciones (uno por cada reacción)
# spec_vector = [[0.3], [0.5], [1.0], [1.0], [1.0]] # Fig1a  
# spec_vector = [[1.], [0.5], [0.1], [1.0], [3.0]] # Fig1b
# spec_vector = [[0.2], [0.5], [1.5], [0.1], [1.5]] # Fig1c

# Diccionario de coeficientes de dispersión (uno por cada especie)
# D_dict = {"l": 0.1, "s1": 0.5, "s2": 0.9}

# Condiciones iniciales: ahora son VECTORES 1D de longitud num_patches para cada especie
# x0_dict = {
#     "l": np.array([1, 0, 0]),
#     "s1": np.array([0, 1, 0]),
#     "s2": np.array([0, 0, 1])
# }

# Matriz de conectividad de orden num_patches × num_patches
# Cada fila debe sumar 1 (representa la probabilidad de dispersión desde ese parche)
# connectivity_matrix = np.array([
#     [0.5, 0.0, 0.5],
#     [0.0, 1.0, 0.0],
#     [0.3, 0.3, 0.4],
# ]) 
import numpy as np

# === Reaction rate constants (spec_vector corresponds to r1–r24)
# Fast reactions: social contagion and emotional responses
# Medium: damage/reconstruction cycles
# Slow: institutional/government responses
# rates=[
#     "mak",
#     "mak",
#     "mak",
#     "mak",
#     "mmk",
#     "mmk",
#     "mmk",
#     "mmk",
#     "mak",
#     "mak",
#     "mak",
#     "mak",
#     "mmk",
#     "mmk",
#     "mmk",
#     "mmk",
#     "mak",
#     "mak",
#     "mak",
#     "mak",
#     "mmk",
#     "mak",
#     "mmk",
#     "mak"
# ]
spec_vector = [
    [0.01],  # r1  A_p + G_A => A_p + G_A + iG_A   (peace spreading idea of good)
    [0.01],  # r2  B_p + G_B => B_p + G_B + iG_B
    [0.01],  # r3  A_v + D_A => A_v + D_A + iD_A   (violence amplifies idea of damage)
    [0.01],  # r4  B_v + D_B => B_v + D_B + iD_B
    [0.1, 0.01],  # r5  A_v + A_p + iD_A => 2A_v + iD_A (violence recruitment)
    [0.1, 0.01],  # r6  B_v + B_p + iD_B => 2B_v + iD_B
    [0.1, 0.01],  # r7  A_v + A_p + iG_A => 2A_p + iG_A (peace reinforcement)
    [0.1, 0.01],  # r8  B_v + B_p + iG_B => 2B_p + iG_B
    [0.01],  # r9  A_v + G_B + iD_A => A_v + D_B   (A violence damages B)
    [0.01],  # r10 B_v + G_A + iD_B => B_v + D_A   (B violence damages A)
    [0.01],  # r11 A_p + D_A + iG_A => A_p + G_A   (reconstruction)
    [0.01],  # r12 B_p + D_B => B_p + G_B
    [0.001,0.001],  # r13 2iG_A => (decay of idea of good)
    [0.001,0.001],  # r14 2iG_B =>
    [0.001,0.001],  # r15 2iD_A => (decay of idea of damage)
    [0.001,0.001],  # r16 2iD_B =>
    [0.01],  # r17 B_p + D_B => B_p + D_B + iD_B (peaceful mind thinks of violence by witnessing damage)
    [0],  # r18:B_p + iD_B => B_v (radicalization)
    [0],  # r19: A_p + D_A => A_v + D_A + iD_A;  (peaceful mind thinks of violence by witnessing damage)
    [0],  # r20: A_p + iD_A => A_v (radicalization)
    [0.0001,0.0001],  # r21 Gov_A + 2iD_A + D_A => Gov_A + G_A (government reparation)
    [0.001],  # r22 Gov_A + G_A => Gov_A + G_A + iG_A (government spreading reconstruction)
    [0.0001,0.0001],  # r23 Gov_B + 2iD_B + D_B => Gov_B + G_B
    [0.001],  # r24 Gov_B + G_B => Gov_B + G_B + iG_B
]

# === Diffusion (dispersal) coefficients for each species
D_dict = {
    "A_p": 0.0, "A_v": 0.0,
    "B_p": 0.0, "B_v": 0.0,
    "G_A": 0.0, "G_B": 0.0,
    "D_A": 0.0, "D_B": 0.0,
    "iG_A": 0.0, "iG_B": 0.0,
    "iD_A": 0.0, "iD_B": 0.0,
    "Gov_A": 0.0, "Gov_B": 0.0,
}

# === Initial conditions (2 locations)
x0_dict = {
    "A_p":  np.array([0.5, 0.0]),
    "A_v":  np.array([0.05, 0.0]),
    "B_p":  np.array([0.55, 0.0]),
    "B_v":  np.array([0.0, 0.0]),
    "G_A":  np.array([0.3, 0.0]),
    "G_B":  np.array([0.32, 0.0]),
    "D_A":  np.array([0.2, 0.0]),
    "D_B":  np.array([0.22, 0.0]),
    "iG_A": np.array([0.1, 0]),
    "iG_B": np.array([0.12, 0.0]),
    "iD_A": np.array([0.1, 0.0]),
    "iD_B": np.array([0.12, 0.0]),
    "Gov_A": np.array([0.6, 0]),
    "Gov_B": np.array([0.6, 0.0]),
}

# === Connectivity matrix (A <-> B exchange)
connectivity_matrix = np.array([
    [1, 0],  # From A’s location
    [0, 1],  # From B’s location
])

num_patches = connectivity_matrix.shape[0]


# Simulation of metapopulation dynamics (vector-based)
t_MP, X_MP, flux_MP = simulate_metapopulation_dynamics(
    rn, 
    #rate=rates, 
    num_patches=num_patches,  # Usar num_patches en lugar de grid_shape
    D_dict=D_dict, 
    x0_dict=x0_dict, 
    spec_vector=spec_vector, 
    t_span=(0, 100), 
    n_steps=200, 
    connectivity_matrix=connectivity_matrix
) 

# Visualization
# Nota: Las funciones de visualización originales están diseñadas para datos en formato 2D (grid)
# Para visualizar datos vectoriales, necesitarás crear nuevas funciones de graficación o
# adaptar las existentes. A continuación se muestra un ejemplo simple usando matplotlib:
# print("species list")
# for sp in species:
#     print(sp)

import matplotlib.pyplot as plt

print(flux_MP["r7"])


species_to_plot_num = ['A_p', 'G_A', 'iG_A', 'B_p', 'G_B', 'iG_B']
species_to_plot_den = ['A_v', 'D_A', 'iD_A', 'B_v', 'D_B', 'iD_B']

# Graficar concentraciones para cada par de especies
fig, axes = plt.subplots(len(species_to_plot_num), 1, figsize=(10, 4 * len(species_to_plot_num)), sharex=True)
if len(species_to_plot_num) == 1:
    axes = [axes]

for idx, sp_up in enumerate(species_to_plot_num):
    sp_down = species_to_plot_den[idx]

    for patch_idx in range(num_patches):
        # Plot numerator species
        axes[idx].plot(
            t_MP,
            X_MP[sp_up][:, patch_idx],
            label=f'{sp_up} (Población {patch_idx+1})',
            linestyle='-'
        )

        # Plot denominator species
        axes[idx].plot(
            t_MP,
            X_MP[sp_down][:, patch_idx],
            label=f'{sp_down} (Población {patch_idx+1})',
            linestyle='--'
        )

    axes[idx].set_ylabel(f'{sp_up} y {sp_down}')
    axes[idx].legend()
    axes[idx].grid(True)

axes[-1].set_xlabel('Tiempo')
fig.suptitle('Series de Tiempo de Concentraciones de Metapoblaciones', fontsize=14)
plt.tight_layout(rect=[0, 0, 1, 0.97])
plt.savefig('metapopulation_concentration_timeseries_vector.png', dpi=150)
plt.show()


#Graficar flujos para cada reacción
# fig, axes = plt.subplots(len(reactions), 1, figsize=(10, 4*len(reactions)), sharex=True)
# if len(reactions) == 1:
#     axes = [axes]

# for idx, rxn in enumerate(reactions):
#     for patch_idx in range(num_patches):
#         axes[idx].plot(t_MP, flux_MP[rxn][:, patch_idx], label=f'Población {patch_idx+1}')
#     axes[idx].set_ylabel(f'{rxn}')
#     axes[idx].legend()
#     axes[idx].grid(True)

# axes[-1].set_xlabel('Tiempo')
# fig.suptitle('Series de Tiempo de Flujos de Metapoblaciones', fontsize=14)
# plt.tight_layout()
# plt.savefig('metapopulation_flux_timeseries_vector.png', dpi=150)
# plt.show()

print("\n✓ Simulación completada exitosamente")
print(f"✓ Forma de salida de concentraciones: {X_MP[species[0]].shape} = (n_steps, num_patches)")
print(f"✓ Forma de salida de flujos: {flux_MP[reactions[0]].shape} = (n_steps, num_patches)")

