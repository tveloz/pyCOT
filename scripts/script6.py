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
file_path = 'Txt/autopoietic.txt'
# file_path = 'Txt/Farm.txt'  
rn = read_txt(file_path)

species = [s.name for s in rn.species()]
print("Especies:", species)
reactions = [r.name() for r in rn.reactions()]
print("Reacciones:", reactions) 

# ##########################################################################
# # Diffusion dynamics in 2D
# ##########################################################################
grid_shape = (3, 3) 

D_dict = {sp: np.round(np.random.uniform(0.01, 0.06), 3) for sp in species} 
x0_dict = {sp: np.round(np.random.uniform(0, 2.0, size=grid_shape), 2) for sp in species}
print("Coeficientes de difusión generadas aleatoriamente:\nD_dict =", D_dict)
print("Condiciones iniciales generadas aleatoriamente:\nx0_dict =", x0_dict)

# Simulación de la dinámica de difusión en 2D
t, X = simulate_diffusion_dynamics_2D(rn, rate='mak', grid_shape=grid_shape,D_dict=D_dict, x0_dict=x0_dict, t_span=(0,5)) 
print("Tiempo de simulación:", t)

print("Concentraciones de especies en las celdas:")
print(X) 

# Dimensiones del vector de tiempo y de las concentraciones
print("Dimensión de t:", t.shape) 
print("Dimensión de X:")
for key, val in X.items():
    print(f"{key}: shape = {val.shape}")

# Visualization 
plot_series_diffusion_2D(t, X, grid_shape)

for sp in X.keys():
    plot_heatmaps_for_species_2D(X, sp, t, time_indices=np.linspace(0, len(t)-1, 3, dtype=int))
    animate_diffusion_heatmaps_for_species_2D(X, sp, t)

##########################################################################
# Metapopulation dynamics simulation
##########################################################################
rate = ['mak'] * len(reactions)
spec_vector = [[0.1]] * len(reactions)
grid_shape = (2, 2)
D_dict = {"l": 0., "s1": 0.1, "s2": 0.5}
x0_dict = {
    "l": np.array([[0.5, 0.5], [0.5, 0.5]]),
    "s1": np.array([[0.5, 0.5], [0.5, 0.5]]),
    "s2": np.array([[0.5, 0.5], [0.5, 0.5]])
}

# Filas sumar 1
prob_matrix = np.array([
    [0.5, 0, 0.5, 0.],
    [0., 1, 0., 0.],
    [0., 0., 0.5, 0.5],
    [0.2, 0., 0., 0.8]
])

# Simulation of metapopulation dynamics
# t, ts_MP, fv_MP = simulate_metapopulation_dynamics(rn, rate)
t, ts_MP, fv_MP = simulate_metapopulation_dynamics(rn, rate, grid_shape, D_dict, x0_dict, spec_vector, prob_matrix=prob_matrix,t_span=(0, 300), n_steps=600)

print("Tiempos:", t)
print("Series temporales (primeras filas):\n", ts_MP.head())
print("Flujos de reacción (primeras filas):\n", fv_MP.head())

X = reconstruct_tensor(ts_MP, species, grid_shape)

# Visualization 
plot_series_diffusion_2D(t, X, grid_shape, xlabel='Tiempo', ylabel='Concentración', title=' Serie de tiempo de Metapoblaciones')

for sp in X.keys():
    plot_heatmaps_for_species_2D(X, sp, t, time_indices=np.linspace(0, len(t)-1, 3, dtype=int))
    animate_diffusion_heatmaps_for_species_2D(X, sp, t)