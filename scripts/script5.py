# Script 5: Exploring the dynamics and abstractions of Networks 

# ========================================
# 1. LIBRARY LOADING AND CONFIGURATION
# ========================================
import os
import sys

# Add the project root directory to the PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import pyCOT modules
from pyCOT.io.functions import read_txt 
from pyCOT.simulations import simulate_discrete_random
from pyCOT.plot_dynamics import plot_series_ode, plot_abstraction_size, plot_abstraction_sets, plot_abstraction_graph_movie_html
from pyCOT.abstractions import abstraction_ordinary

# ========================================
# 2. CREATING THE REACTION_NETWORK OBJECT
# ========================================
# # Network file path
file_path = 'Txt/Farm.txt'   
rn = read_txt(file_path)     

# ========================================
# 3. Discrete Random Simulation
# ========================================
# # Parameters  of the network
x01 = [14.0,  13.0,     5.0,   70.0,   5.0,    3.0,    3.0,     9.0,          10.0,        6.0,    10.0,    2.0,     12.0,     2.0]

ts, vs = simulate_discrete_random(rn, S=rn.stoichiometry_matrix(), x=x01, n_iter=31)

plot_series_ode(ts)

# ========================================
# 4. Abstracting the Time Series
# ======================================== 
abstract_time_series = abstraction_ordinary(ts, threshold=0.1)

plot_abstraction_size(abstract_time_series)

plot_abstraction_sets(abstract_time_series)

plot_abstraction_graph_movie_html(abstract_time_series, filename="abstraction_graph_movie_html_farm.html", interval=400, title="Abstraction Graph - Time")
