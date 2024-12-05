############################################################################## 
import os  # Imports the 'os' module to interact with the operating system
import sys  # Imports the 'sys' module to interact with the Python interpreter
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # Adds the project's root directory to PYTHONPATH to allow importing modules from that location

# # Imports from the pyCOT library   
from pyCOT.file_manipulation import load_pyCOT_from_file
from pyCOT.rn_simulations import * 

# # Example of using functions
file_path = 'Txt/autopoietic.txt'
# file_path = 'Txt/2019fig1.txt'
# file_path = 'Txt/2019fig2.txt'
# file_path = 'Txt/non_connected_example.txt'
# file_path = 'Txt/Farm.txt' 

# Generate the pyCOT object from the file containing the txt network
testRN=load_pyCOT_from_file(file_path)  
print(testRN.SpStr)

##############################################################################
# # Examples
##############################################################################
# Generates an initial state vector based on the number of species in the reaction network
x_inicial = generate_state_vector(len(testRN.SpStr))  
print(x_inicial)

# Calculates the stoichiometric matrix associated with the reaction network
S=stoichiometric_matrix(testRN) 
print(S)

# Simulates the time series for 20 iterations, using the stoichiometric matrix and the initial vector of states
time_series, flux_raction = iterate_state_vector(S, x_inicial, testRN.SpStr, n_iter=20)
print(time_series)
print(flux_raction)