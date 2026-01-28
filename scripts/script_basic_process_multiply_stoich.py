#!/usr/bin/env python3
"""
Simple script to load a reaction network and multiply a process vector 
by the stoichiometric matrix.
"""

# Add the root directory to the PYTHONPATH
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from pyCOT.io.functions import read_txt

# Load reaction network
file_path = 'networks/Conflict_Theory/Resource_Community_Insurgency_Loops_model1.txt'
rn = read_txt(file_path)

# Get stoichiometric matrix
S = rn.stoichiometry_matrix()

# Choose process vector coefficients (one for each reaction)
# Here we choose arbitrary positive values for a few reactions
n_reactions = S.shape[1]
# Create a process vector with NON-NEGATIVE random numbers
process_vector = np.abs(np.random.randn(n_reactions))
# process_vector = np.zeros(n_reactions)
# process_vector[0] = 1.0   # r1
# process_vector[5] = 2.0   # r6
# process_vector[10] = 0.5  # r11

# Multiply stoichiometric matrix by process vector
result = S @ process_vector

# Print results
print("Process vector (reaction rates):")
print(process_vector)
print(f"\nShape: {process_vector.shape}")

print("\n" + "="*50)
print("Result (species production rates):")
print(result)
print(f"\nShape: {result.shape}")

print("\n" + "="*50)
print("Species names:")
for i, species in enumerate(rn._species_map.keys()):
    print(f"{i}: {species} -> {result[i]:.2f}")