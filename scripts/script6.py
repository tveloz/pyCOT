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
file_path = 'Txt/autopoietic.txt'
rn = read_txt(file_path)
species = [s.name for s in rn.species()]
print("Especies:", species)
reactions = [r.name() for r in rn.reactions()]
print("Reacciones:", reactions)

rate = 'mak'
grid_shape = (2, 2)
# spec_vector = [[0.3], [0.5], [1.0], [1.0], [1.0]] # Fig1a  # [[0.1]] * len(reactions) 
# spec_vector = [[1.], [0.5], [0.1], [1.0], [3.0]] # Fig1b
spec_vector = [[3.], [0.05], [6.0], [0.05], [2.0]] # Fig1c

D_dict = {"l": 0., "s1": 0.5, "s2": 0.9}
x0_dict = {
    "l": np.array([[0.5, 0.5], [0.5, 0.5]]),
    "s1": np.array([[0.5, 0.5], [0.5, 0.5]]),
    "s2": np.array([[0.5, 0.5], [0.5, 0.5]])
}

# Las filas deben sumar 1
connectivity_matrix = np.array([
    [0.5, 0, 0.5, 0.],
    [0., 1, 0., 0.],
    [0., 0., 0.5, 0.5],
    [0.2, 0., 0., 0.8]
]) 

# Simulation of metapopulation dynamics
t_MP, X_MP, flux_MP = simulate_metapopulation_dynamics(rn, rate, grid_shape=grid_shape, D_dict=D_dict, x0_dict=x0_dict, spec_vector=spec_vector, t_span=(0,100), n_steps=200, connectivity_matrix=connectivity_matrix) 

# Visualization  
plot_diffusion_time_series_2D(time=t_MP,concentration_data=X_MP, grid_shape=grid_shape, xlabel='Tiempo', ylabel='Concentración', main_title='Serie de tiempo de las Concentraciones de Metapoblaciones', legend_title='Especies', cell_prefix='G',filename='metapopulation_concentration_timeseries.png')
plot_diffusion_time_series_2D(time=t_MP,concentration_data=flux_MP, grid_shape=grid_shape, xlabel='Tiempo', ylabel='Concentración', main_title='Serie de tiempo del Flujo de Metapoblaciones', legend_title='Flujos', cell_prefix='G',filename='metapopulation_flux_timeseries.png')

plot_heatmaps_all_species_2D(t_MP, X_MP, time_indices=3, main_title="Evolución de los perfiles de concentración de Metapoblaciones",filename='metapopulation_concentration_heatmaps.png')
animate_diffusion_heatmaps_all_species_2D(t_MP, X_MP, main_title="Animación de la Dinámica de Metapoblaciones") 