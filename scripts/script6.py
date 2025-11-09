# Script 6: Diffusion and Metapopulations of Reaction Networks with pyCOT 

# Import necessary libraries and modules
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from pyCOT.io.functions import read_txt
from pyCOT.simulations import *
from pyCOT.plot_dynamics import * 

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

rate = 'mak'
num_patches = 2  # Número de poblaciones/parches

# Parámetros de reacciones (uno por cada reacción)
# spec_vector = [[0.3], [0.5], [1.0], [1.0], [1.0]] # Fig1a  
# spec_vector = [[1.], [0.5], [0.1], [1.0], [3.0]] # Fig1b
#spec_vector = [[0.2], [0.5], [1.5], [0.1], [1.5]] # Fig1c

# Diccionario de coeficientes de dispersión (uno por cada especie)
#D_dict = {"l": 0.1, "s1": 0.5, "s2": 0.9}

# Condiciones iniciales: ahora son VECTORES 1D de longitud num_patches para cada especie
# x0_dict = {
#     "l": np.array([1, 0, 0]),
#     "s1": np.array([0, 1, 0]),
#     "s2": np.array([0, 0, 1])
# }

# Matriz de conectividad de orden num_patches × num_patches
# Cada fila debe sumar 1 (representa la probabilidad de dispersión desde ese parche)
#connectivity_matrix = np.array([
#     [0.5, 0.0, 0.5],
#     [0.0, 1.0, 0.0],
#     [0.3, 0.3, 0.4],
# ]) 

# Simulation of metapopulation dynamics (vector-based)
t_MP, X_MP, flux_MP = simulate_metapopulation_dynamics(
    rn, 
    rate, 
    num_patches=num_patches,  # Usar num_patches en lugar de grid_shape
    #D_dict=D_dict, 
    #x0_dict=x0_dict, 
    #spec_vector=spec_vector, 
    t_span=(0, 100), 
    n_steps=200, 
    #connectivity_matrix=connectivity_matrix
) 

# Visualization
# Nota: Las funciones de visualización originales están diseñadas para datos en formato 2D (grid)
# Para visualizar datos vectoriales, necesitarás crear nuevas funciones de graficación o
# adaptar las existentes. A continuación se muestra un ejemplo simple usando matplotlib:

import matplotlib.pyplot as plt

# Graficar concentraciones para cada especie
fig, axes = plt.subplots(len(species), 1, figsize=(10, 4*len(species)), sharex=True)
if len(species) == 1:
    axes = [axes]

for idx, sp in enumerate(species):
    for patch_idx in range(num_patches):
        axes[idx].plot(t_MP, X_MP[sp][:, patch_idx], label=f'Población {patch_idx+1}')
    axes[idx].set_ylabel(f'{sp}')
    axes[idx].legend()
    axes[idx].grid(True)

axes[-1].set_xlabel('Tiempo')
fig.suptitle('Series de Tiempo de Concentraciones de Metapoblaciones', fontsize=14)
plt.tight_layout()
plt.savefig('metapopulation_concentration_timeseries_vector.png', dpi=150)
plt.show()

# Graficar flujos para cada reacción
fig, axes = plt.subplots(len(reactions), 1, figsize=(10, 4*len(reactions)), sharex=True)
if len(reactions) == 1:
    axes = [axes]

for idx, rxn in enumerate(reactions):
    for patch_idx in range(num_patches):
        axes[idx].plot(t_MP, flux_MP[rxn][:, patch_idx], label=f'Población {patch_idx+1}')
    axes[idx].set_ylabel(f'{rxn}')
    axes[idx].legend()
    axes[idx].grid(True)

axes[-1].set_xlabel('Tiempo')
fig.suptitle('Series de Tiempo de Flujos de Metapoblaciones', fontsize=14)
plt.tight_layout()
plt.savefig('metapopulation_flux_timeseries_vector.png', dpi=150)
plt.show()

print("\n✓ Simulación completada exitosamente")
print(f"✓ Forma de salida de concentraciones: {X_MP[species[0]].shape} = (n_steps, num_patches)")
print(f"✓ Forma de salida de flujos: {flux_MP[reactions[0]].shape} = (n_steps, num_patches)")

