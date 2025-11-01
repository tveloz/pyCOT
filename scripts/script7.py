# Script 8: Simulation with Process Classification and Flux Combination Histogram

# ========================================
# 1. LIBRARY LOADING AND CONFIGURATION
# ========================================
import os
import sys

# Add the project root directory to the PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pyCOT.io.functions import read_txt
from pyCOT.simulations import *
from pyCOT.plot_dynamics import *

# ========================================
# 2. CREATING THE REACTION_NETWORK OBJECT
# ========================================
file_path = 'Txt/autopoietic.txt'  
# file_path = 'Txt/Farm.txt'
file_path = 'networks/testing/Lotka_Volterra.txt'
rn = read_txt(file_path)

species = [specie.name for specie in rn.species()]
print("\nSpecies Set =",species)

reactions = [reaction.name() for reaction in rn.reactions()]
print("Reactions Set =",reactions,"\n")

S = rn.stoichiometry_matrix()

# ========================================      
# 3. SIMULATIONS
# ========================================
x0 = [2, 2, 2] 
rate_list = 'mak'  # Kinetics for all reactions

spec_vector = [[.7], [.5],[1],[1], [1]] # Sv=[-0.8,0.5,-1]#  Se estabiliza en 3 concentraciones distintas

time_series, flux_vector = simulation(
    rn, rate=rate_list, spec_vector=spec_vector, 
    x0=x0, 
    t_span=(0, 50), n_steps=1000+1 
)
# time_series, flux_vector = simulation(
#     rn, rate=rate_list,
#     t_span=(0, 250), n_steps=1000+1 
# )

# Gráfico combimado: Concentraciones, Flujos, Histograma para la configuración inicial
fig, axes = plt.subplots(1, 3, figsize=(18, 4))
plot_series_with_domain_intervals(time_series, flux_vector, S, title="Serie de Tiempo de Concentraciones", save_figure=False, ax=axes[0])
plot_flux_with_domain_intervals(flux_vector, S, title="Serie de Tiempo de Flujos", save_figure=False, ax=axes[1])
plot_process_types_histogram(flux_vector, S, title="Histograma de Tipos de Procesos", save_figure=False, ax=axes[2])
plt.tight_layout()
plt.savefig("./visualizations/process_classification/combined_plots.png", dpi=150, bbox_inches="tight")
plt.show()