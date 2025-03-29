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

# Load the reaction network
file_path = 'Txt/Langmuir_Hinshelwood.txt'  

# Generate the pyCOT object from the file containing the txt network
testRN = load_pyCOT_from_file(file_path)

#####################################################################################
# # Example 1.1: Time series of ODE with default parameters and initial conditions
###################################################################################### 
# # Solve EDO 
# time_series, flux = simulate_ode_mak(testRN)
# print("Time series of ODE:")
# print(time_series)

# # Plots the time series of ODE
# plot_series_ode(time_series)

###################################################################################### 
# # Example 1.2: Time series of ODE with parameters and initial conditions  
###################################################################################### 
# Initial values of species concentrations
x0=[1, 1, 0, 1, 0, 0]
print("Initial concentrations for species:")
print(x0)

# Time span
t_span = (0, 9)
print("Time span for the simulation:")
print(t_span)

# Steps for the simulation
n_steps = 200
print("Number of steps for the simulation:")
print(n_steps)

# Values of the reaction rate constants
k0=[3, 0.5, 2, 0.1, 1]
print("Rate constants for reactions:")
print(k0)

# Solve EDO
time_series, flux_series = simulate_ode_mak(testRN, x0, t_span, n_steps, k0) 
print("Time series of ODE:")
print(time_series)

# Plots the time series of ODE
plot_series_ode(time_series)