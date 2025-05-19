######################################################################################
import os
import sys
# Add the root directory to the PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import pandas as pd

from pyCOT.file_manipulation import load_pyCOT_from_file
from pyCOT.abstractions import abstraction_ordinary
from pyCOT.rn_simulations import *

#######################################################################################
# Example 1
#######################################################################################
# Create a time series
data = {
    "Time": [0, 1, 2, 3, 4, 5],
    "A": [0.5, 0.8, 0.3, 0.1, 0.0, 0.4],
    "B": [0.0, 0.2, 0.4, 0.6, 0.8, 0.5],
    "C": [0.9, 0.7, 0.1, 0.0, 0.0, 0.2]
}
time_series = pd.DataFrame(data)

# Define a fixed threshold or per-species thresholds
threshold = 0.2  # Fixed threshold
# threshold = {"A": 0.4, "B": 0.5, "C": 0.2}  # Per-species thresholds

# Calculate the abstraction
# abstract_time_series = abstraction_ordinary(time_series)
abstract_time_series = abstraction_ordinary(time_series, threshold)
print(abstract_time_series)

#######################################################################################
# # Example 2

# file_path = 'Txt/autopoietic.txt'
# # file_path = 'Txt/2019fig1.txt'
# # file_path = 'Txt/2019fig2.txt'
# # file_path = 'Txt/non_connected_example.txt'
# # file_path = 'Txt/Farm.txt'

# testRN = load_pyCOT_from_file(file_path)

# # Generate and simulate the initial state vector
# x_initial = generate_state_vector(len(testRN.SpStr))  # Should take the number of rows (species) in S
# print(x_initial)

# S = stoichiometric_matrix(testRN)
# print(S)

# time_series2,flux_vector = iterate_state_vector(S, x_initial, testRN.SpStr, n_iter=10)
# print(time_series2)

# # Define a fixed threshold or per-species thresholds
# # threshold2 = 0.2  # Fixed threshold
# threshold2 = {"l": 0.4, "s1": 0.5, "s2": 0.4}  # Per-species thresholds

# # Calculate abstractions
# # abstract_time_series = abstraction_ordinary(time_series)
# abstract_time_series2 = abstraction_ordinary(time_series2, threshold2)

# print(abstract_time_series2)