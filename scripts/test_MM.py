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
file_path = 'Txt/Michaelis_Menten.txt'  

# Generate the pyCOT object from the file containing the txt network
testRN = load_pyCOT_from_file(file_path)

#####################################################################################
# # Example 1.1: Time series of ODE with default parameters and initial conditions
###################################################################################### 
# # Solve EDO 
# time_series, flux_series = simulate_ode_mak(testRN)
# print("Time series of ODE:")
# print(time_series)

# # Plots the time series of ODE
# plot_series_ode(time_series)

###################################################################################### 
# # Example 1.2: Time series of ODE with parameters and initial conditions with different seeds
###################################################################################### 
# # Initial values of species concentrations
# x0=[0.4, 0.03, 0.79, 1]
# print("Initial concentrations for species:")
# print(x0)

# # Time span
# t_span = (0, 25)
# print("Time span for the simulation:")
# print(t_span)

# # Steps for the simulation
# n_steps = 500
# print("Number of steps for the simulation:")
# print(n_steps)

# # Values of the reaction rate constants
# k0=[0.2036, 0.485, 0.1478]
# print("Rate constants for reactions:")
# print(k0)

# # Solve EDO
# time_series, flux_series = simulate_ode_mak(testRN, x0, t_span, n_steps, k0) 
# print("Time series of ODE:")
# print(time_series)

# # Plots the time series of ODE
# plot_series_ode(time_series)

######################################################################################
# # # Example 1.3: Michaelis-Menten Kinetics
######################################################################################
# Initial conditions for species
x0 = [1.0, 10.0, 0.0, 0.0] # E, S, C (ES), P

# Time range for the simulation
t_span = (0, 50)

# Number of steps for the simulation
n_steps = 200

# Values of the reaction rate constants for reactions
k0 = {1.0, 0.5, 0.3}
 
# Maximum reaction rates for species

Vmax_dict = {1.0, 25.0, 0.5, 5.0}

# Michaelis constants for species
Km_dict = {1.5, 1.9, 0.5, 0.8} 
#####################################################################################
# Solve the ODEs
# time_series, flux_vector  = simulate_ode_mmk(testRN) # Default values for all Reaction Network
# time_series, flux_vector  = simulate_ode_mmk(testRN, x0) 
# time_series, flux_vector  = simulate_ode_mmk(testRN, x0, t_span, n_steps, k0)
time_series, flux_vector  = simulate_ode_mmk(testRN, x0, t_span, n_steps, k0, Vmax_dict, Km_dict)
print("Time series of MMK:")
print(time_series)

# Plot the time series of concentration and flux for the Michaelis-Menten system
plot_series_ode(time_series, title="Michaelis-Menten MMK: ODE")










# import numpy as np
# import matplotlib.pyplot as plt
# from scipy.integrate import odeint

# # Definimos los parámetros
# k1 = 1.0   # Constante de unión (1/(M·s))
# k_1 = 0.5  # Constante de disociación (1/s)
# k2 = 0.3   # Constante de conversión a producto (1/s)
# E_total = 1.0  # Concentración total de la enzima (M)
# E0 = E_total # Concentración inicial de enzima (M)
# S0 = 10.0    # Concentración inicial de sustrato (M)
# C0 = 0.0     # No hay complejo inicialmente
# P0 = 0.0     # No hay producto inicialmente

# # Definimos las ecuaciones diferenciales
# def michaelis_menten_system(y, t):
#     E, S, C, P = y
#     E = E_total - C  # Enzima libre
    
#     dE_dt = -k1 * E * S + (k_1 + k2) * C
#     dS_dt = -k1 * E * S + k_1 * C
#     dC_dt =  k1 * E * S - (k_1 + k2) * C
#     dP_dt =  k2 * C  # Velocidad de formación del producto
    
#     return [dE_dt, dS_dt, dC_dt, dP_dt]

# # Condiciones iniciales
# y0 = [E0, S0, C0, P0]

# # Vector de tiempo
# t = np.linspace(0, 50, 200)  # De 0 a 50 segundos, 200 puntos

# # Resolver la EDO
# sol = odeint(michaelis_menten_system, y0, t)

# # Extraer las soluciones
# E, S, C, P = sol[:, 0], sol[:, 1], sol[:, 2] , sol[:, 3]

# V_max = k2 * E0  # Velocidad máxima de la reacción
# print("V_max:")
# print(V_max)

# Km = (k_1 + k2) / k1  # Constante de Michaelis-Menten
# print("Km:")
# print(Km)

# vel=V_max*S/(Km+S)
 

# # Graficar la velocidad de formación del producto
# plt.figure(figsize=(8, 6))
# plt.plot(t, E, 'c', label='E')
# plt.plot(t, S, 'r', label='S')
# plt.plot(t, C, 'k', label='C')
# plt.plot(t, P, 'b', label='P')
# plt.xlabel('Tiempo (s)')
# plt.ylabel('Contentration (M)')
# plt.title('Cinética de Michaelis-Menten')
# plt.legend()
# plt.grid()
# plt.show()

# plt.figure(figsize=(8, 6)) 
# plt.plot(S, vel, 'b', label='S vs vel')
# plt.xlabel('Concentración de Sustrato [S] (M)')
# plt.ylabel('Velocidad de reacción v (M/s)')
# plt.title('Curva de Michaelis-Menten')
# plt.legend()
# plt.grid()
# plt.show()