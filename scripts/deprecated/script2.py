# # # Script 2: Simulation of Dynamics 
# This script shows how to simulate the dynamics of a reaction network (RN)
# with selected initial conditions and dynamic parameters. 
# The script is structured to: 
#    (a) load a reaction network, 
#    (b) select initial conditions of concentration and dynamic parameters for mak, 
#    (c) generate time series simulation,
#    (d) plot concentration and flux time series in all available ways
#    (e) plot abstractions in all possible ways 
 
# # Note: To comment a line of code, use Ctrl + C, and to uncomment it, use Ctrl + K.

#####################################################################
# Imports from the pyCOT library
#####################################################################
# Import Python standard library modules
import os  # Imports the 'os' module to interact with the operating system
import sys # Imports the 'sys' module to interact with the Python interpreter

# Add the root directory to the PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))) # Adds the project's root directory to PYTHONPATH to allow importing modules from that location
sys.stdout.reconfigure(encoding='utf-8')                                     # Set the standard output encoding to UTF-8 (pls commment out if using Windows)

# Imports from the pyCOT library 
from pyCOT.file_manipulation import load_pyCOT_from_file # Import only the load_pyCOT_from function of the file_manipulation module to load RN  
from pyCOT.rn_visualize import *                         # Imports all functions, classes and variables defined in the rn_visualize module

from pyCOT.simulations import *           # General simulation tools for dynamic systems modeled in pyCOT.
from pyCOT.abstractions import *          # Provides methods to compute abstractions, such as threshold-based state simplifications.
from pyCOT.plot_dynamics import *         # Functions for plotting dynamic behaviors, such as time series or abstracted data.
from pyCOT.rn_types import StoichiometryMatrix
from pyCOT.abstractions import abstraction_ordinary 

#####################################################################
# (a) Load a reaction network
#####################################################################
# # File path
#file_path = '../Txt/autopoietic.txt' 
# file_path = '../Txt/2019fig1.txt'
#file_path = '../Txt/2019fig2.txt'
# file_path = '../Txt/non_connected_example.txt' 
# file_path = '../Txt/PassiveUncomforableIndignated_problemsolution.txt'
file_path = '../Txt/Farm.txt' 
file_path = 'networks/testing/Riverland0.txt'
file_path = 'networks/testing/autopoietic.txt'

# # Loads the RN from the specified file
testRN = load_pyCOT_from_file(file_path)           # Creates an object called testRN, which represents the RN

#####################################################################
# (b) Select initial conditions of concentration and dynamic parameters for mak
#####################################################################
# Initial species concentration conditions
x0_values = generate_random_vector(len(testRN.SpStr), seed=11).round(2) # Generate random values for the initial concentrations of species
#x0_values = [0.6, 0.7, 0.3] # Specific initial concentrations for each species, need to modify based on current reaction network

# Dynamic parameters for mak
k0_values = generate_random_vector(len(testRN.RnStr), seed=71).round(2) # Generate random values for the specific reaction rate constants
# k0_values = [0.55, 0.71, 0.29, 0.51, 0.89] # Specific reaction rate constants for each reaction, need to modify based on number of reactions in network
print('x_0='+str(x0_values))
print('k_0='+str(k0_values))


#####################################################################
# (c) generate time series simulation
#####################################################################
#time_series, flux_vector = simulate_ode_mak(testRN, t_span = [0, 1], n_steps=10)   (for example 2019fig2.txt)
#time_series = simulate_ode_mak(testRN, x0=x0_values, k0=k0_values)
time_series, flux_vector = simulate_ode_mak(testRN, x0=x0_values, t_span = [0, 80], n_steps = 400, k0=k0_values)
#print("time_series:")
#print(time_series)
diff_series = time_series.shift(-1) - time_series[:-1]
#print("flux_vector:")
#print(flux_vector)

#####################################################################
# (d) plot concentration and flux time series in all available ways
#####################################################################
plot_series_ode(time_series)
diff_series.plot()
plt.show
#plot_series_ode(flux_vector)

#####################################################################
# (e) plot abstractions in all possible ways 
#####################################################################
# # # Define a fixed threshold or per-species thresholds
threshold = 0.01  # Fixed threshold
#threshold = {"l": 0.4, "s1": 0.5, "s2": 0.4}  # Per-species thresholds for autopoietic.txt

# # # Generates an abstraction with a threshold
# # abstract_time_series = abstraction_ordinary(time_series)          # Generates an abstraction with a default threshold of 0.5
abstract_time_series = abstraction_ordinary(time_series, threshold) # Generates an abstraction with a defined threshold for each species.
# print(abstract_time_series)

#############
# Plots of abstractions
#############
# # Plot of the number of abstractions
#plot_abstraction_size(abstract_time_series)

# # Plot the abstractions sorted by size
plot_abstraction_sets(abstract_time_series)

#############
# Histograms
#############
# # Generates a histogram of all species concentrations in the NR over time
#plot_join_concentration_histogram(time_series, bins=10)

# # Print reaction rates v(x,k) 
# # Generates histogram of reaction rates v(x,k)
#plot_reaction_rate_mak_histogram(flux_vector, bins=5)

# # Generates individual histograms for each species
#plot_species_histograms(time_series, species_names=testRN.SpStr)

# # Generates a combined histogram of the concentrations of all species.
#plot_combined_species_histogram(time_series, species_names=testRN.SpStr)

plot_combined_species_histogram(diff_series, species_names=testRN.SpStr)

#Obtener data frame a partir de diff que lleve cuenta del momento de la simulacion, de la abstraccion en operacion, del vector diff, y si el diff representa (challenge, problem, cognitive domain, cognitive control).
#Extender la construccion a barrido del dataframe que vaya sumando los diffs, y considere entonces intervalos de largo n, y agregue abstraccion inicial, abstraccion final, abtraccion mas repetida, y numero total de abstracciones.


