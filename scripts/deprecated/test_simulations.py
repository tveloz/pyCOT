import os
import sys
# Add the root directory to the PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pyCOT.simulations import *
from pyCOT.reaction_network import *
from pyCOT.closure_structure import *  
from pyCOT.file_manipulation import *   
from pyCOT.plot_dynamics import * 
from pyCOT.abstractions import *
from pyCOT.rn_types import StoichiometryMatrix

#########################################################################################
import numpy as np
import numpy as np
from scipy.integrate import odeint
 

############################################################################################ Load the reaction network
file_path = 'Txt/autopoietic.txt' 
# file_path = 'Txt/2019fig1.txt'
# file_path = 'Txt/2019fig2.txt'
# file_path = 'Txt/non_connected_example.txt' 
# file_path = 'Txt/PassiveUncomforableIndignated_problemsolution.txt'
# file_path = 'Txt/Farm.txt' 

# Generate the pyCOT object from the file containing the txt network
RN = load_pyCOT_from_file(file_path)
print(RN.SpStr)
print(RN.RnStr)

# List of kinetics for the reactions
rate_list = 'mak' 
# rate_list = ['mmk', 'mak', 'mmk', 'mmk', 'mak']
# rate_list = ['mmk', 'mak', 'mmk', 'mmk', 'mak', 'mak']
print("List of kinetics:")
print(rate_list)

# List of parameters for the reactions
spec_vector = [[60.0, 2750.0], [0.71], [270.0, 3550.0], [60.0, 2750.0], [0.89]]
# spec_vector = [[60.0, 2750.0], [0.71], [270.0, 3550.0], [60.0, 2750.0], [0.89], [0.9]]

# Ejecutar la simulación utilizando la cinética de Michaelis-Menten
time_series, flux_vector = simulations(RN, rate=rate_list)
# time_series, flux_vector = simulations(RN, rate=rate_list, spec_vector=spec_vector)

print("Time series of ODE:")
print(time_series)
plot_series_ode(time_series)

# print("Flux vector:")
# print(flux_vector)
# plot_series_ode(flux_vector)