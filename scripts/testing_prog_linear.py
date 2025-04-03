#####################################################################
# Imports from the pyCOT library
#####################################################################
# Import Python standard library modules 
import os  # Imports the 'os' module to interact with the operating system
import sys # Imports the 'sys' module to interact with the Python interpreter

# Add the root directory to the PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))) # Adds the project's root directory to PYTHONPATH to allow importing modules from that location
sys.stdout.reconfigure(encoding='utf-8')                                     # Set the standard output encoding to UTF-8

# Imports from the pyCOT library 
#from pyCOT.file_manipulation import load_pyCOT_from_file # Import only the load_pyCOT_from function of the file_manipulation module to load RN  
from pyCOT.io.functions import read_txt # Import the load_pyCOT_from_file function from the io.functions module to load RN
from pyCOT.closure_structure import *
#from pyCOT.simulations import *            # General simulation tools for dynamic systems modeled in pyCOT.
from pyCOT.rn_types import StoichiometryMatrix

import numpy as np # Import the numpy library as np
from scipy.optimize import linprog # Import the linprog function from the scipy.optimize module

#####################################################################
# (a) Load a reaction network
#####################################################################
# # File path
file_path = 'networks/autopoietic.txt' 
file_path = 'Txt/2019fig1.txt'        # No solution was found
# file_path = 'Txt/2019fig2.txt' 
# file_path = 'Txt/non_connected_example.txt' 
# file_path = 'Txt/PassiveUncomforableIndignated_problemsolution.txt'
#file_path = 'networks/FarmVariants/Farm.txt' 

# (a) Loads the RN from the specified file
testRN = read_txt(file_path)           # Creates an object called testRN, which represents the RN
# print(testRN.RnMprod)
# print(testRN.RnMsupp) 
# print(testRN.SpStr)

#####################################################################
# # # Example
#####################################################################
# Stoichiometric matrix  
# S = universal_stoichiometric_matrix(testRN)  


# # # Nombres de las filas (índices)
# columns = testRN.RnStr
# index = testRN.SpStr
# # index = ['water', 'grass', 'cows', 'infr', 'milk', 'dung', 'worms', 'fertilizer', 
# #          'chickens', 'eggs', 'grain', 'straw', 'money', 'farmer']
# S = pd.DataFrame(S, columns=testRN.RnStr, index=index)  # Creates a DataFrame object with the stoichiometric matrix data

# # # Intercambiar las filas 'money' y 'farmer'
# S.loc[['farmer', 'money']] = S.loc[['money', 'farmer']].values
# print(S) # Muestra el DataFrame en la consola

# # v1 = np.array([2, 1, 2, 1, 1])  # Process Vector of autopoietic
# # # v1 = [11.,  3.,  3.,  1.,  2.,  5.,  1.,  2.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.]  # Process Vector of Farm

# v1 = [5.3, 1.4, 1.4, 0.5, 1.0, 2.4, 0.5, 1.0, 0.5, 1.0, 1.9, 7.7, 0.5, 0.5, 0.4, 0.5, 0.5] # S.v1 < 0 (Falló)
# # v1= [14.0, 4.0, 3.0, 1.0, 5.0, 7.0, 1.0, 2.0, 1.0, 5.0, 5.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0] # S.v2>=0
# # print(v1) 

# # Comprobation: S.v>=0
# xv10 = S @ v1
# print("x_v1 = S.v1 =",xv10)

#####################################################################
print("checking self-maintainance of the whole network")
X=testRN.species() # Get the names of all species in the reaction network
res = is_self_maintaining(testRN,X)
s_orgs=reactive_semi_orgs(testRN)

print("checking self-maintainance of each semi_self")
for x in s_orgs:
    print("set X="+str(x))
    res = is_self_maintaining(testRN,x)
    print(res)

# Calculate S.v_opt
# x_v = S @ v_opt
# print("x_v=S.v_opt:")
# for i, val in enumerate(x_v):
#     print(f"x_v[{i+1}] = {val} {'✅ Yes' if val >= 0 else '❌ No'}") 

#####################################################################  


# #####################################################################
# # Ejemplo: Página 13 de Paper 1
# #####################################################################
# v1 = np.array([2, 1, 2, 1, 1])
# v2 = np.array([3, 1, 2, 1, 1]) 
# v3 = np.array([3, 2, 2, 1, 1]) 
# v4 = np.array([3, 2, 3, 1, 1]) 
# v5 = np.array([4, 3, 3, 1, 1])

# B= np.transpose([v1, v2, v3, v4, v5])
# print("Matriz B:\n",B)
# ref_matrix_B = row_echelon_form(B)
# print("Matriz B escalonada:\n",ref_matrix_B)

# #####################################################################
# # Ejemplo 9 (Kolman-Hill, Álgebra Lineal, 8va edición, p. 311)
# #####################################################################
# A= np.array([[1, -1, 1,0,0,0], [0,-1,0,1,0,0], [1,-1,0,0,1,0], [0,0,0,0,0,1]])
# ref_matrix_A = row_echelon_form(A)
# print("Matriz A escalonada:\n",ref_matrix_A)
