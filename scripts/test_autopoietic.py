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
file_path = 'Txt/autopoietic.txt' 
# file_path = 'Txt/2019fig1.txt'
# file_path = 'Txt/2019fig2.txt' # No funcionó
# file_path = 'Txt/non_connected_example.txt' 
# file_path = 'Txt/PassiveUncomforableIndignated_problemsolution.txt'
# file_path = 'Txt/Farm.txt'
# file_path = 'Txt/RN_IN_04.txt' # No grafica la EDO, o depende de la escala 

# Generate the pyCOT object from the file containing the txt network
testRN = load_pyCOT_from_file(file_path)

#####################################################################################
# # Example 1.1: Time series of ODE with default parameters and initial conditions
###################################################################################### 
# # Solve EDO 
# time_series = simulate_ode_mak(testRN)
# print("Time series of ODE:")
# print(time_series)

# # Plots the time series of ODE
# plot_series_ode(time_series)

###################################################################################### 
# # Example 1.2: Time series of ODE with parameters and initial conditions with different seeds
###################################################################################### 
# Initial values of species concentrations
x0=[2, 1, 0]
print("Initial concentrations for species:")
print(x0)

# Time span
t_span = (0, 100)
print("Time span for the simulation:")
print(t_span)

# Steps for the simulation
n_steps = 500
print("Number of steps for the simulation:")
print(n_steps)

# Values of the reaction rate constants
k0=[1.0, 0.05, 0.3, 0.2, 0.1] 
print("Rate constants for reactions:")
print(k0)

# Solve EDO
time_series, flux = simulate_ode_mak(testRN, x0, t_span, n_steps, k0) 
print("Time series of ODE:")
print(time_series)

# Plots the time series of ODE
plot_series_ode(time_series)

######################################################################################
# # Example 1.3: PassiveUncomforableIndignated_problemsolution
######################################################################################
# # Load the reaction network
# file_path = 'Txt/PassiveUncomforableIndignated_problemsolution.txt' 
# testRN = load_pyCOT_from_file(file_path)
 
# # Define specific reaction rate constants
# k = {
#     testRN.RnStr[0]: 0.1,
#     testRN.RnStr[1]: 0.1,
#     testRN.RnStr[2]: 0.2,
#     testRN.RnStr[3]: 0.3,
#     testRN.RnStr[4]: 0.3, 
#     testRN.RnStr[5]: 0.3,
#     testRN.RnStr[6]: 0.3, 
#     testRN.RnStr[7]: 0.3, 
#     testRN.RnStr[8]: 0.15,
#     testRN.RnStr[9]: 0.15,                    
# }

# # Time span
# t_span = (0, 300)

# # Define initial concentrations for species
# y0 = [
#     0.3,  # Initial concentration for species 1: D
#     0.4,  # Initial concentration for species 2: U
#     0.3,  # Initial concentration for species 3: S
#     0.1,  # Initial concentration for species 4: s
#     0.9,  # Initial concentration for species 5: p
# ]

# ######################################################################################
# # Solve EDO
# time_series = simulate_ode_mak(testRN, k, t_span, x0,n_steps=1000)
# print(time_series)

# # Plots the time series of ODE
# plot_series_ode(time_series, xlabel="Time", ylabel="Concentration", title="Time Series of Concentrations")

########################################################################################
########################################################################################
########################################################################################
# # # Example 3: 
# file_path = 'Txt/autopoietic.txt'
# # file_path = 'Txt/2019fig1.txt'
# # file_path = 'Txt/2019fig2.txt'
# # file_path = 'Txt/non_connected_example.txt'
# file_path = 'Txt/Farm.txt' 

# # # Generate the pyCOT object from the file containing the txt network
# testRN=load_pyCOT_from_file(file_path)  
# print(testRN.SpStr)
 
# # # Generates an initial state vector based on the number of species in the reaction network
# x_inicial = generate_random_vector(len(testRN.SpStr))  
# print(x_inicial)

# # # Calculates the stoichiometric matrix associated with the reaction network
# matrix_data=universal_stoichiometric_matrix(testRN) 
# print(matrix_data)

# S=StoichiometryMatrix(matrix_data, species=testRN.SpStr, reactions=testRN.RnStr) 
# print(S)

# # # Simulates the time series for 20 iterations, using the stoichiometric matrix and the initial vector of states
# time_series, reaction_rate_maks = simulate_discrete_random(testRN, S=S, x=x_inicial, n_iter=300)
# print(time_series)
# # print(reaction_rate_maks)

# # # Define a fixed threshold or per-species thresholds
# threshold = 0.05  # Fixed threshold
# # threshold = {"A": 0.4, "B": 0.5, "C": 0.4}  # Per-species thresholds

# # # Calculate abstractions
# # abstract_time_series = abstraction_ordinary(time_series)
# abstract_time_series = abstraction_ordinary(time_series, threshold)
# print(abstract_time_series)

# Plots the time series of ODE
# plot_series_ode(time_series, xlabel="Time", ylabel="Concentration", title="Time Series of Concentrations")
# plot_abstraction_size(abstract_time_series, xlabel="Time", ylabel="Number of Species", title="Number of species per abstraction over time", marker='o', label="Abstraction Size")
# plot_abstraction_sets(abstract_time_series)

# # # Graphs
# # Static graph
# plot_static_abstraction_graph(abstract_time_series)  

# # Movie
# plot_abstraction_graph_movie(abstract_time_series, interval=1000)

# # Histograms
# plot_join_concentration_histofram(time_series, bins=30, alpha=0.7, color='skyblue', edgecolor='black', xlabel="Concentration", ylabel="Frequency", title="Histogram of Species Concentrations")
# plot_reaction_rate_mak_histogram(reaction_rate_maks, bins=10, color='skyblue', edgecolor='black', alpha=0.7, xlabel="Reaction Rate", ylabel="Frequency", title="Histogram of Reaction Rates")
# plot_species_histograms(time_series, species_names=testRN.SpStr, bins=10, alpha=0.7, color="skyblue", edgecolor="black", xlabel="Concentration", ylabel="Frequency", title_plot="Histogram of")
# plot_combined_species_histogram(time_series, species_names=testRN.SpStr, bins=10, alpha=0.7, edgecolor="black", xlabel="Concentration", ylabel="Frequency", title_plot="Combined Histogram")

########################################################################################
########################################################################################
########################################################################################
# # # # Example 4: Michaelis-Menten Kinetics
# # Define specific reaction rate constants
# k0 = {
#     testRN.RnStr[0]: 0.65,
#     testRN.RnStr[1]: 1.0,
#     testRN.RnStr[2]: 2.0,
#     testRN.RnStr[3]: 1.5,
#     testRN.RnStr[4]: 2.5,    
# }
 
# # Define specific reaction rate constants
# Vmax_dict = {
#     testRN.SpStr[0]: 1.0,
#     testRN.SpStr[1]: 1.2,
#     testRN.SpStr[2]: 0.8,    
# }

# # Michaelis constants
# Km_dict = {
#     testRN.SpStr[0]: 0.5,
#     testRN.SpStr[1]: 0.7,
#     testRN.SpStr[2]: 0.6,    
# }

# # Initial conditions
# x0 = [1.0, 2.0, 3.5]  # Initial concentrations of l, s1, s2

# # Time range for the simulation
# t_span = (0, 20)

# n_steps = 200
# #####################################################################################
# # Solve the ODEs
# # time_series, flux_vector  = simulate_ode_mmk(testRN) # Default values for all Reaction Network
# # time_series, flux_vector  = simulate_ode_mmk(testRN, x0) 
# # time_series, flux_vector  = simulate_ode_mmk(testRN, x0, t_span, n_steps, k0)
# time_series, flux_vector  = simulate_ode_mmk(testRN, x0, t_span, n_steps, k0, Vmax_dict, Km_dict)
# print("Time series of MMK:")
# print(time_series)

#####################################################################################
# Plot the time series of concentration and flux for the Michaelis-Menten system
plot_series_ode(time_series, title="Michaelis-Menten MMK: ODE")
# # plot_series_ode(flux_vector, title="Michaelis-Menten MMK: ODE")

# import numpy as np
# import matplotlib.pyplot as plt
# from scipy.integrate import odeint

# # Definimos los parámetros de la red autopoietic
# k1 = 1.0   # Producción de L (1/s)
# k2 = 0.5   # Conversión de L + S1 -> 2S1 (1/(M·s))
# k3 = 0.3   # Conversión de S1 -> S2 (1/s)
# k4 = 0.8   # Conversión de S2 + L -> S1 (1/(M·s))
# k5 = 0.6   # Conversión de S2 -> (1/s)
# L0 = 1.0   # Concentración inicial de L (M)
# S1_0 = 0.5 # Concentración inicial de S1 (M)
# S2_0 = 0.5 # Concentración inicial de S2 (M)

# # Definimos las ecuaciones diferenciales para la red autopoietic
# def autopoietic_network(y, t):
#     L, S1, S2 = y
    
#     dL_dt = k1 - k2 * L * S1
#     dS1_dt = k2 * L * S1 - k3 * S1 + k4 * S2 * L
#     dS2_dt = k3 * S1 - k5 * S2
    
#     return [dL_dt, dS1_dt, dS2_dt]

# # Condiciones iniciales
# y0 = [L0, S1_0, S2_0]

# # Vector de tiempo
# t = np.linspace(0, 50, 200)  # De 0 a 50 segundos, 200 puntos

# # Resolver la EDO
# sol = odeint(autopoietic_network, y0, t)

# # Extraer las soluciones
# L, S1, S2 = sol[:, 0], sol[:, 1], sol[:, 2]

# # Calcular velocidades de reacción
# v_S1 = k2 * L * S1  # Velocidad de formación de S1
# v_S2 = k3 * S1      # Velocidad de formación de S2

# # Determinar V_max y K_m para cada reacción
# V_max_S1 = k2 * max(L)  # Velocidad máxima de formación de S1
# K_m_S1 = k1 / k2        # Constante de Michaelis-Menten para S1

# V_max_S2 = k3 * max(S1)  # Velocidad máxima de formación de S2
# K_m_S2 = k3 / k4         # Constante de Michaelis-Menten para S2

# # Curvas de Michaelis-Menten
# vel_S1 = V_max_S1 * L / (K_m_S1 + L)
# vel_S2 = V_max_S2 * S1 / (K_m_S2 + S1)

# # Graficar la Curva de Michaelis-Menten para S1
# plt.figure(figsize=(8, 6))
# plt.plot(L, vel_S1, 'b', label='S1 vs v_S1')
# plt.xlabel('Concentración de L (M)')
# plt.ylabel('Velocidad de reacción v (M/s)')
# plt.title('Curva de Michaelis-Menten para S1')
# plt.legend()
# plt.grid()
# plt.show()

# # Graficar la Curva de Michaelis-Menten para S2
# plt.figure(figsize=(8, 6))
# plt.plot(S1, vel_S2, 'r', label='S1 vs v_S2')
# plt.xlabel('Concentración de S1 (M)')
# plt.ylabel('Velocidad de reacción v (M/s)')
# plt.title('Curva de Michaelis-Menten para S2')
# plt.legend()
# plt.grid()
# plt.show()


















########################################################################################
########################################################################################
########################################################################################
# # # Example 5: Metapolulation
# # Define the parameters  
# t_span=(0, 10)
# n_steps=200

# # Tasas de reacción
# k=np.array([1.1, 0.05, 0.2, 0.4, 0.5]) 

# # Coeficientes de difusión de las especies, que modula la tasa a la que cada especie se propaga en la grilla.
# D = np.array([0, 1, 0]) 

# Un valor mayor de D_i implica una difunsión más rápida y una mayor homogeneización espacial de la especie.
# Un valor menor de D_i implica una difunsión más lentamente, lo que retarda la homogenización espacial y puede mantener gradientes de concentración más marcados.

########################################################################################
# Ejemplo de malla 2x2
######################
# grid_shape=(2, 2)

# # Condiciones iniciales (x0) con forma (2, 2, 3)
# x0=np.array([
#     [[5, 8, 2], [1, 12, 1]],  # Fila inferior de la malla
#     [[2, 10, 3], [3, 5, 4]]   # Fila superior de la malla
# ])  

# # Tasas de intercambio entre nodos, en una grilla 2x2
# exchange_rates = np.array([
#     [[1.3, 0.3, 1.3], [2.5, 0.8, 0.3]],  
#     [[0.7, 1.1, 1.6], [0.5, 1.8, 0.9]]    
# ])  

# ########################################################################################
# # Ejemplo de malla 3x3
# ######################
# grid_shape = (3, 3) 

# # Condiciones iniciales (x0) con forma (3, 3, 3)
# x0 = np.array([
#     [[5, 8, 2], [1, 12, 1], [3, 4, 7]],
#     [[2, 10, 3], [3, 5, 4], [6, 2, 9]],
#     [[4, 7, 5], [8, 1, 2], [2, 3, 6]]
# ])  

# # Tasas de intercambio entre nodos, en una grilla 3x3
# exchange_rates = np.array([
#     [[0.7, 0.3, 1.3], [1.7, 0.8, 0.3], [1.2, 0.4, 0.9]],
#     [[0.1, 1.1, 1.6], [0.9, 1.8, 0.9], [0.5, 1.2, 1.0]],
#     [[0.8, 0.6, 1.4], [1.3, 1.0, 0.7], [1.0, 0.9, 1.2]]
# ]) 

# ########################################################################################
# # # Ejemplo de malla 5x5
# ######################
# grid_shape = (5, 5) 

# # Condiciones iniciales (x0) con forma (5, 5, 3) para 3 especies por nodo
# x0 = np.array([
#     [[5, 8, 2], [1, 12, 1], [3, 4, 7], [2, 9, 3], [4, 6, 5]],
#     [[2, 10, 3], [3, 5, 4], [6, 2, 9], [7, 8, 2], [1, 11, 4]],
#     [[4, 7, 5], [8, 1, 2], [2, 3, 6], [5, 9, 3], [3, 7, 1]],
#     [[3, 9, 4], [2, 6, 7], [4, 8, 5], [6, 3, 9], [2, 5, 8]],
#     [[7, 2, 6], [5, 4, 3], [3, 7, 1], [4, 9, 2], [8, 1, 5]]
# ])

# # Tasas de intercambio entre nodos para 5x5 (3 especies por nodo)
# exchange_rates = np.array([
#     [[0.7, 0.3, 1.3], [1.7, 0.8, 0.3], [1.2, 0.4, 0.9], [0.8, 0.6, 1.4], [1.0, 0.5, 1.1]],
#     [[0.1, 1.1, 1.6], [0.9, 1.8, 0.9], [0.5, 1.2, 1.0], [1.2, 0.7, 1.3], [0.7, 1.0, 1.5]],
#     [[0.8, 0.6, 1.4], [1.3, 1.0, 0.7], [1.0, 0.9, 1.2], [0.9, 1.1, 1.0], [1.4, 0.8, 0.6]],
#     [[1.2, 1.0, 0.5], [0.7, 1.2, 1.3], [1.1, 0.8, 1.0], [0.6, 0.9, 1.4], [1.0, 1.1, 0.8]],
#     [[0.9, 1.0, 1.2], [1.3, 0.7, 1.1], [0.8, 1.2, 0.6], [1.1, 0.9, 0.7], [1.0, 0.5, 1.3]]
# ])

########################################################################################
# # # Simulate the ODEs for the metapopulation model 
# modules = simulate_odes_metapop_mak(testRN, x0=x0, t_span=t_span, n_steps=n_steps, k=k, exchange_rates=exchange_rates, D=D, grid_shape=grid_shape) 
# # print("Modules:")
# # print(modules)

# # Plot the time series of the metapopulation model
# plot_series_MP(modules, testRN, t_span=t_span, n_steps=n_steps)

# # Plot the dynamics of the metapopulation model in 2D
# plot_all_species_2d(results=modules, species_names=testRN.SpStr)

# # Plot the dynamics of the metapopulation model in 3D
# plot_all_species_3d(results=modules, species_names=testRN.SpStr)

#######################################
# # Definir un umbral para la abstracción
# threshold_value = 0.5  # Puede ser un diccionario para umbrales específicos por especie

# # Calcular la abstracción para cada módulo (parche)
# # abstracted_results = abstract_metapopulation(modules, threshold=threshold_value)
# abstracted_results = apply_abstraction_to_simulation(modules, species_names=testRN.SpStr, threshold=threshold_value, t_span=t_span, title="Abstraction")

# # Mostrar la abstracción de todos los parches
# for patch, abstraction in abstracted_results.items():
#     print(f"Abstracción del {patch}:")
#     print(abstraction, "\n")

########################################################################################
# # Example 6: Michaelis-Menten Kinetics in the metapopulation model
# modules = simulate_odes_metapop_mmk(testRN, t_span=t_span, n_steps=n_steps, grid_shape=grid_shape)

# plot_series_MP(modules, testRN, t_span=t_span, n_steps=n_steps,title_plot="Metapopulation Michaelis-Menten Kinetics")

########################################################################################
# # # Example 7: PDE simulation
# # simulation_data=simulate_pde_rd(testRN, t_span=t_span, n_steps=n_steps,grid_shape=grid_shape)
# simulation_data=simulate_pde_rd(testRN, x0=x0, t_span=t_span, n_steps=n_steps, k=k, exchange_rates=exchange_rates, D=D,grid_shape=grid_shape)
# # print("Simulation data:")
# # print(simulation_data)

# # Plot the time series of the PDE simulation
# plot_series_PDE(simulation_data, species_names=testRN.SpStr, t_span=t_span)  

# # Animate the time series of the PDE simulation
# animate_series_PDE(simulation_data, species_names=testRN.SpStr, t_span=t_span)