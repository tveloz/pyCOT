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
# file_path = 'Txt/autopoietic_d.txt'  
# file_path = 'Txt/Farm.txt'

rn = read_txt(file_path)

species = [specie.name for specie in rn.species()]
print("\nSpecies Set =",species)

reactions = [reaction.name() for reaction in rn.reactions()]
print("Reactions Set =",reactions,"\n")

S = rn.stoichiometry_matrix()
print("Stoichiometry Matrix S =\n", S)
print("Shape of S =", S.shape)

# ========================================      
# 3. SIMULATIONS
# ========================================
x0 = [2, 2, 2] 
rate_list = 'mak'  # Kinetics for all reactions
spec_vector = [[.7], [.5], [1],[1], [1]]

# x0 = [2, 2, 2, 0.1] 
# rate_list = 'mak'  # Kinetics for all reactions
# spec_vector = [[.7], [.5],[1],[1], [1],[0.55],[0.55], [0.55]] # Sv=[-0.8,0.5,-1]#  Se estabiliza en 3 concentraciones distintas

tf=100
n=1000

time_series, flux_vector = simulation(rn, rate=rate_list, spec_vector=spec_vector, x0=x0, t_span=(0, tf), n_steps=n+1)

# # ========================================
# # 4. ANALYSIS: Process Type Proportions Over Time
# # ========================================
results_df = analyze_process_proportions_over_time(rn, S, rate_list, spec_vector, x0, t_span=(0, tf), n_steps=n+1,
    window_sizes=list(range(1,7,1)), save_path="./visualizations/process_classification/")