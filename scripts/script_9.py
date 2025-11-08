# Script 9: Process Classification Tests for Reaction Network
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
file_path = 'Txt/autopoietic.txt'    # 2 rayos, 10 proyecciones en 3D 
rn = read_txt(file_path)

species = [specie.name for specie in rn.species()]
print("\nSpecies Set =",species)

reactions = [reaction.name() for reaction in rn.reactions()]
print("Reactions Set =",reactions,"\n")

S = rn.stoichiometry_matrix()
print("Stoichiometry Matrix S:\n", S.shape)

#########################################################################
# 3. PROCESS CLASSIFICATION TESTS
#########################################################################
# --- Función auxiliar para imprimir resultados de forma clara ---
def test_and_print(process_name, v, S, v_prev=None): 
    print(process_name)
    print(f"  Vector de Proceso (v): {v}")
    sv_result = S @ v
    print(f"  Cambios Netos (S*v) : {sv_result}") 
    classification = classify_process_mode(v, S)
    print(f"classified as {classification}")
    is_cog = is_cognitive_domain(v, S)
    print(f" Is in Cognitive Domain? {is_cog}")
    if v_prev is not None:
        print(f" Disturbance:  {v_prev}")
        sv_prev_result = S @ v_prev
        print(f"  Disturbance effect (S*v_prev): {sv_prev_result}")
        print(f"  Disturbance and response effect (S*(v + v_prev)): {S @ (v + v_prev)}")
        disturbance_classification = classify_response_to_disturbance(v, S, v_prev, verbose=True)
        print(f"disturbance classified as {disturbance_classification}")
    return
        

# --- Casos de prueba basados en el documento ---

# 0. Proceso v0: Un modo estacionario.
v0 = np.array([2, 1, 2, 1, 1])
test_and_print("v0 ∈ (Stationary Mode)", v0, S)

#1. Proceso v1: Un modo estacionario perfecto[cite: 249].
v1 = np.array([4, 2, 1, 0, 0])
test_and_print("v1 ∈ (Feasible)", v1, S)    

# 2. Proceso v2: Una solución que sobreproduce s1[cite: 251].
v2 = np.array([3, 1, 2, 1, 1])
test_and_print("v2 ∈ (Cognitive Domain)", v2, S)

# 3. Proceso vc3 de la Eq. (9)[cite_start]: Un desafío que consume s1[cite: 281].
vc3 = np.array([0, 0, 1, 0, 0])
test_and_print("vc3 ∈ (Challenge)", vc3, S)

# 4. Proceso que solo consume s2 (r5): Un problema[cite: 281].
vp4 = np.array([0, 0, 0, 0, 1])
test_and_print("vp4 ∈ (Problem)", vp4, S) # [ 0  0 -1] -> vp4 es un ['Problem']

#5. Proceso vc5 como contramedida (Counteraction) del desafío vc5.
vc5 = np.array([2, 2, 1, 2, 1])
v1 = np.array([2, 0, 2.5, 0, 0]) # COMPENSACIÖN
# test_and_print("vc5 (Challenge)", vc5, S) # [-2  3 -2] -> vc5 es un ['Challenge'] 
test_and_print("v5 (Counteraction)", vc5, S, v_prev=v1) # [0.  0.5 0.5] ->  v1 es un ['Counteraction']

#6. Proceso v6 (Solution) controlando el problema vp6. 
v6 = np.array([4, 2, 1, 0, 0]) # Feasible np.array([5, 1, 2, 2, 1]) # [ 2  1 -1] DESAFIO
vp6 = np.array([0, 0, 2, 2, 1]) # [-2  0 -1] Problema # vp6 = np.array([0, 0, 1, 1, 0]) # [-1  0  0] Problema # 
test_and_print("v6 ∈ (Solución)", v6, S, v_prev=vp6) 

# 7. Proceso v7 (Cognitive Control) controlando el problema v7. 
v7 = np.array([4, 1, 2, 1, 1]) # [2 0 0] Cognitive domain
vp7 = np.array([0, 0, 1, 1, 0]) # [-1  0  0] Problema
test_and_print("v7 ∈ (Cognitive Control)", v7, S, v_prev=vp7)