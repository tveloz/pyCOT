# Script 6: Diffusion and Metapopulations of Reaction Networks with pyCOT 

# Import necessary libraries and modules
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from pyCOT.io.functions import read_txt
from pyCOT.simulations import *
from pyCOT.plot_dynamics_old import * 

################################################################################################# 
# Ejemplo de uso  
file_path = 'Txt/Farm.txt'  
file_path = 'networks/Conflict_Theory/cause_driven_conflict_gov.txt'

rn = read_txt(file_path)

species = [s.name for s in rn.species()]
print("Especies:", species)
reactions = [r.name() for r in rn.reactions()]
print("Reacciones:", reactions) 

# ##########################################################################
# # Diffusion dynamics in 2D
# ##########################################################################
rate='mak'
grid_shape = (2, 1) 

# Simulación de la dinámica de difusión en 2D
t_diff, X_diff, flux_diff = simulate_diffusion_dynamics_2D(rn, rate, grid_shape=grid_shape, t_span=(0,10)) 

# Visualization 
plot_diffusion_time_series_2D(time=t_diff,concentration_data= X_diff, grid_shape= grid_shape, xlabel='Tiempo', ylabel='Concentración', main_title='Serie de tiempo de las Concentraciones por Difusión', legend_title='Especies', cell_prefix='G')
plot_diffusion_time_series_2D(time=t_diff,concentration_data=flux_diff, grid_shape=grid_shape, xlabel='Tiempo', ylabel='Concentración', main_title='Serie de tiempo del Flujo por Difusión', legend_title='Flujos', cell_prefix='G')

plot_heatmaps_all_species_2D(t_diff, X_diff, time_indices=3, main_title="Evolución de los perfiles de concentración por Difusión")
animate_diffusion_heatmaps_all_species_2D(t_diff, X_diff, main_title="Animación de la Dinámica de Difusión") 


