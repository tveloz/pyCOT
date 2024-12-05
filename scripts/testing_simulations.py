import os
import sys
# Add the root directory to the PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pyCOT.simulations import *
from pyCOT.reaction_network import *
from pyCOT.closure_structure import * 
from pyCOT.file_manipulation import *   
from pyCOT.plot_dynamics import *

#####################################################################################
# Example 1: autopoietic 
######################################################################################
# Load the reaction network
file_path = 'Txt/autopoietic.txt' 
testRN = load_pyCOT_from_file(file_path)
 
# Define specific reaction rate constants
k = {
    testRN.RnStr[0]: 0.5,
    testRN.RnStr[1]: 1.0,
    testRN.RnStr[2]: 2.0,
    testRN.RnStr[3]: 1.5,
    testRN.RnStr[4]: 2.5,    
} 

# Time span
t_span = (0, 40)

# Define initial concentrations for species
y0 = [
    0.8,  # Initial concentration for species 1: l
    0.3,  # Initial concentration for species 2: s1
    0.5,  # Initial concentration for species 3: s2
]

######################################################################################
# Solve EDO
sol=solve_ode(testRN, k, t_span, y0,n_steps=100)

#CONTRUIR grafico con plot_dynamics: plot_series_ode(time_series)

# Plot the results
plt.figure(figsize=(10, 6))
for i, species in enumerate(testRN.SpStr):
    plt.plot(sol.t, sol.y[i], label=species)
plt.xlabel('Time')
plt.ylabel('Concentration')
# plt.title('Mass action kinetics MAK: ODE')
plt.legend()
plt.grid()
plt.show() 

######################################################################################
for i, species in enumerate(testRN.SpStr):
    data_concentrations = {
        "Time": sol.t,
        i: sol.y[i]
    }
time_series = pd.DataFrame(data_concentrations) 
 


# #####################################################################################
# # Example 2: PassiveUncomforableIndignated_problemsolution
# ######################################################################################
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
# sol=solve_ode(testRN, k, t_span, y0)

# # Plot the results
# plt.figure(figsize=(10, 6))
# for i, species in enumerate(testRN.SpStr):
#     plt.plot(sol.t, sol.y[i], label=species)
# plt.xlabel('Time')
# plt.ylabel('Concentration')
# # plt.title('Mass action kinetics MAK: ODE')
# plt.legend()
# plt.grid()
# plt.show() 

