#!/usr/bin/env python3
"""
Test script for the Gillespie (SSA) simulator.

Adapted to the DataFrame-based API: `gillespie` now mirrors
`pyCOT.simulations.ode.simulation` and returns (time_series_df, flux_vector_df)
instead of a GillespieResult object.

Main changes with respect to the old script:
  - `rate_constants={"R1": 0.7, ...}`  ->  `rate='mak'` + `spec_vector=[[0.7], ...]`
  - `t_max=10`                         ->  `t_span=(0, 10)`
  - `result.times` / `result.populations` -> columns of `ts` (Time, species...)
  - `result.plot()`                    ->  plotted here from the DataFrame

Note on ordering: `spec_vector` is indexed by the order of
`stoichiometry_matrix().reactions`, NOT by the order in which you write a dict.
The helper `spec_from_dict` below builds it safely from reaction names.
"""

import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt

# -- Path setup ---------------------------------------------------------------
root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(root_dir, 'src'))  # Add src/ to PYTHONPATH for pyCOT

from pyCOT.io.functions import read_txt, print_reaction_network
from pyCOT.visualization.rn_visualize import *
from pyCOT.visualization.plot_dynamics import plot_series_ode
# Adjust this import to wherever you placed the module
from pyCOT.simulations.stochastic import gillespie

# -----------------------------------------------------------------------------
# Network
# -----------------------------------------------------------------------------

file_path = 'data/Examples_tests/autopoietic.txt'  # Change this path as needed
# file_path = 'data/Examples_tests/Farm.txt'
# file_path = 'data/Examples_tests/2019fig2.txt'
# file_path = 'data/Ecological_models/AMF.txt'
# file_path = 'data/Ecological_models/AMF_final.txt'

rn = read_txt(file_path)


species = [s.name for s in rn.species()]
print("Species:", species)

reactions = [reaction.name() for reaction in rn.reactions()]
print("Reactions:", reactions)

# -----------------------------------------------------------------------------
# Simulation
# -----------------------------------------------------------------------------
# # # Example 1: Automatic initial condition and parameters  
# ts, flux = gillespie(
#     rn,
#     rate='mak', 
#     t_span=(0, 100),
#     n_steps=100+1,
#     seed=42,
#     verbose=False,      # Print with True: x0, spec_vector, ODEs, velocity expressions
# )                       

# # Example 2: Explicit initial condition and parameters
x0 = {"l": [10], "s1": [15], "s2": [5]}  
spec_vector=[[0.5], [0.6], [0.2], [0.1], [0.1]] 
ts, flux = gillespie(
    rn,
    rate='mak',
    spec_vector=spec_vector,  
    x0=x0,
    t_span=(0, 100),
    n_steps=101,
    seed=32 
)                      







# -----------------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------------
# Figure 1: species populations over time
print("\nTime series:\n", ts) 
print(f"\nTime series final at t={ts['Time'].iloc[-1]:.4f}:\n", ts[species].iloc[-1].to_dict())
plot_series_ode(ts, xlabel="Time", ylabel="Number of molecules", 
                title=f"Gillespie SSA - {Path(file_path).stem}",
                filename=f"gillespie_{Path(file_path).stem}.png",show_fig=True)

# # Figure 2: flux vector over time
# print("\nFlux vector:\n", flux) 
# print(f"\nFlux vector final at t={flux['Time'].iloc[-1]:.4f}:\n", flux.iloc[-1].drop('Time').to_dict())
# plot_series_ode(flux, xlabel="Time", ylabel="Reaction firings", title=f"Gillespie SSA - {Path(file_path).stem} (flux vector)", filename=f"gillespie_flux_{Path(file_path).stem}.png", show_fig=True)    