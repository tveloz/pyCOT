# Import necessary libraries and modules
import os
import sys
# Add the root directory to the PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import the necessary modules from pyCOT  
from pyCOT.io.functions import read_txt 
from pyCOT.plot_dynamics import plot_series_ode
from pyCOT.simulations import *

###################################################################################
# # # Load the reaction network from a file
file_path = 'Txt/autopoietic.txt' 
# file_path = 'Txt/2019fig1.txt'
# file_path = 'Txt/2019fig2.txt'
# file_path = 'Txt/non_connected_example.txt' 
# file_path = 'Txt/PassiveUncomforableIndignated_problemsolution.txt'
#file_path = 'Txt/Farm.txt' 
#file_path = 'networks/testing/Lake Chad Basin System_0.txt' 

rn = read_txt(file_path) # Creates the ReactionNetwork object 

species =   [specie.name for specie in rn.species()]
print(species)                  # Lista de nombres de las especies

reactions = [reaction.name() for reaction in rn.reactions()]
print(reactions)                # Lista de nombres de las reacciones
 
#########################################################
### Example 1: Simulation of the reaction network only with 'mak' 
#########################################################
x0=[0, 1, 0] # Initial conditions for the simulation

# rate_list = ['mmk','mak','mak','mak','mak']                               # List of kinetics for the reactions
# spec_vector=[[0.7,0.5], [0.5], [1.0], [1.0], [1.0]] # List of parameters for the reactions

rate_list = ['mmk','mmk','mmk','mmk','mmk']             # List of kinetics for the reactions
spec_vector = [[0.818, 17], [1.241, 19], [1.256, 0.13], [0.202, 0.7], [0.358, 13]]

time_series, flux_vector= simulation(rn, rate=rate_list, spec_vector=spec_vector, x0=x0, t_span=(0, 50),n_steps=500)

print("Time series of ODE:")
print(time_series)
plot_series_ode(time_series)

#########################################################
### Example 2: Simulation of the reaction network with 'mak' and 'mmk'
#########################################################
x0=[0, 1, 0] # Initial conditions for the simulation

rate_list = ['mmk', 'mmk', 'mmk', 'mmk', 'mmk']             # List of kinetics for the reactions
spec_vector = [[1.0, 0.3], [1.0, 0.4], [1.0, 0.3], [1.0, 0.4], [1.0, 0.4]] # List of parameters for the reactions
 
time_series, flux_vector= simulation(rn, rate=rate_list, spec_vector=spec_vector, x0=x0, t_span=(0, 50),n_steps=500)

print("Time series of ODE:")
print(time_series)
# plot_series_ode(time_series)

#########################################################
### Example 3: Register additional kinetics with randomly generated parameters
#########################################################
# rate_list = ['mak', 'mak', 'mmk', 'ping_pong', 'hill']  # List of kinetics for the reactions

# additional_laws = {'ping_pong': rate_ping_pong}
# time_series, flux_vector = simulation(rn, rate=rate_list, t_span=(0, 50), n_steps=500, additional_laws=additional_laws)

# print(build_reaction_dict(rn))

# print("Time series of ODE:")
# print(time_series)
# plot_series_ode(time_series)

<<<<<<< HEAD
#########################################################
### Example 4: Register additional kinetics with defined parameters 
#########################################################
#x0 = [0.1, 0.5, 0.2]  # Initial concentrations of species
#rate_list = ['mak', 'mak', 'mmk', 'ping_pong', 'hill']  # List of kinetics for the reactions
    # spec_vector = [
    #     [0.5],           # R1: mak
    #     [0.3],           # R2: mak
    #     [1.2, 0.4],      # R3: mmk
    #     [0.9, 0.5, 0.1], # R4: ping_pong
    #     [0.1, 0.9, 0.5]  # R5: hill
    # ]

additional_laws = {'mmk': rate_ping_pong}
time_series, flux_vector = simulation(rn,   
                          t_span=(0, 50), n_steps=500, 
                         additional_laws=additional_laws)
=======
# #########################################################
# ### Example 4: Register additional kinetics with defined parameters 
# #########################################################
# x0 = [0.1, 0.5, 0.2]  # Initial concentrations of species

# rate_list = ['mak', 'mak', 'mmk', 'ping_pong', 'hill']  # List of kinetics for the reactions
# spec_vector = [
#     [0.5],           # R1: mak
#     [0.3],           # R2: mak
#     [1.2, 0.4],      # R3: mmk
#     [0.9, 0.5, 0.1], # R4: ping_pong
#     [0.1, 0.9, 0.5]  # R5: hill
# ]
>>>>>>> abfe8eeb859921e23d6f38489786c1856857c8b6

# additional_laws = {'ping_pong': rate_ping_pong}
# time_series, flux_vector = simulation(rn, rate=rate_list, spec_vector=spec_vector, 
#                          x0=x0, t_span=(0, 50), n_steps=500, 
#                          additional_laws=additional_laws)

# print("Time series of ODE:")
# print(time_series)
# plot_series_ode(time_series)

# ###########################################################
# # Example 5: Simulation of the reaction network with discrete random simulation
# ###########################################################
# mm=universal_stoichiometric_matrix(rn)
# print("Universal Stoichiometric Matrix:\n",mm) # Imprime la representación de la matriz de estequiometría

# ts, vs =simulate_discrete_random(rn, S=rn.stoichiometry_matrix(), x=x0, n_iter=11)
# print("\nTime series of discrete random simulation:")
# print(ts)
# plot_series_ode(ts)

# print("\nFlux vector of discrete random simulation:")
# print(vs)
# plot_series_ode(vs)