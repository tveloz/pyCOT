# Script 8: Cone Polytopes in Reaction Networks
##########################################
# Ejemplo de uso
# ========================================
# 1. LIBRARY LOADING AND CONFIGURATION
# ========================================
import os
import sys # Add the project root directory to the PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pyCOT.io.functions import read_txt
from pyCOT.simulations import *
from pyCOT.plot_dynamics import *

# ========================================
# 2. CREATING THE REACTION_NETWORK OBJECT
# ========================================
# file_path = 'Txt/Lotka_Volterra.txt' # 3 reacciones, 1 rayo, 1 proyección en 3D
file_path = 'Txt/autopoietic.txt'    # 5 reacciones, 2 rayos, 10 proyecciones en 3D
# file_path = 'TxT/2019fig1.txt'       # 6 reacciones, 3 rayos, 20 proyecciones en 3D
# file_path = 'TxT/2019fig2.txt'       # 6 reacciones, 3 rayos, 20 proyecciones en 3D
# file_path = 'Txt/autopoietic_d.txt'  # 8 reacciones, 4 rayos, 56 proyecciones en 3D
# 9 reacciones, 5 rayos, 84 proyecciones en 3D

# Con redes con más de 9 reacciones, el cálculo del cono puede ser muy lento y consumir mucha memoria.
# file_path = 'Txt/PassiveUncomforableIndignated_problemsolution.txt'    # NO CORRE: 10 reacciones, 6 rayos, 120 proyecciones en 3D
# file_path = 'Txt/Scenario1_baseline_only_reactions.txt'                # NO CORRE: 17 reacciones, 9 rayos, 680 proyecciones en 3D

rn = read_txt(file_path)

species = [specie.name for specie in rn.species()]
print("\nSpecies Set =",species)

reactions = [reaction.name() for reaction in rn.reactions()]
print("Reactions Set =",reactions,"\n")

S = rn.stoichiometry_matrix() 
print(S.species)
print(S.reactions)

#########################################################################  
# # Ejemplo 1: Cone
# #######################################################################
plot_cone_and_region(S, show=True) 
# print("Null vectors of S (S*v=0):\n",null_vectors)

null_vectors=compute_nullspace_vectors(S)
N = np.column_stack(null_vectors)

# Ejemplo de cálculo de combinación lineal N c = v
import numpy as np
np.random.seed(42)
# v=np.random.uniform(low=0.0, high=10.0, size=(S.shape[1])).round(1) 

# v = np.array([1., 0., 1., 1., 0.]) + np.array([1., 1., 1., 0., 1.])
v =[0.5, 0.25, 0.5, 0.25, 0.25]
print("Vector original v:", v) 

# Resolver N c = v
c, residuals, rank, s = np.linalg.lstsq(N, v, rcond=None)

print("Coeficientes de combinación lineal:", c)
print("Residual ||N c - v||:", np.linalg.norm(N @ c - v))
v_reconstructed = N @ c #c@N.T
print("v reconstruido:", np.round(v_reconstructed,2))



# ################################################################################
# # Ejemplo 2: Vectores extra en la gráfica del cono  de AUTOPOIETIC
# ################################################################################
# x0 = np.array([2, 2,2])
# vc = np.array([2, 2, 1, 2, 1]) # Es un desafío para l y s2, porque se consumen (v_i<0)
# print(S@vc,": Challenge") # [-2.  3. -2.]

# vp = np.array([0, 0, 1, 0, 1]) # Es una problema para s1, porque se consume (v_i<0) y no produce ninguna especie 
# print(S@vp,": Problem") # [ 0. -1.  0.] 

# v1 = np.array([2, 0, 2.5, 0, 0]) # Es una Contramedida al desafío vc (v(0) = (2, 2, 1, 2, 1) problema, v(1) = (2, 0, 2.5, 0, 0) solución)
# print(S@v1,": Challenge")    # [0.  0.5 0.5]

# v=vc+v1 
# print(S@(v),": Counteraction")    # [0.  0.5 0.5]

# # Definir los vectores extra y sus etiquetas personalizadas
# extra_vectors = [vc, vp, v1, v]
# extra_labels = ["Vector desafío (vc)", "Vector problema (vp)", "Vector Contramedida (v1)", "Vector combinado (v=vc+v1)"]
# extra_colors = ['orange', 'cyan', 'purple', 'brown'] 

# print("Vectores extra =", extra_vectors)

# # # Llamar a la función con etiquetas y colores personalizados 
# plot_cone_and_region(S, show=False, extra_vector=extra_vectors, extra_vector_labels=extra_labels, extra_vector_colors=extra_colors)