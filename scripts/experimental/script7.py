#!/usr/bin/env python3
"""
Script 7 (REFACTORED): Simulation with Process Classification and Flux Combination Histogram

This script demonstrates:
1. Loading a reaction network from a text file
2. Running ODE simulations
3. Analyzing process type proportions over multiple time scales
4. Visualizing cognitive domain dynamics

UPDATES:
- Now uses refactored plot_dynamics.py and process_analysis.py
- Cleaner separation between analysis and visualization
- All classification logic uses process_analysis module
"""

# ========================================
# 1. LIBRARY LOADING AND CONFIGURATION
# ========================================
import os
import sys

# Add the project root directory to the PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pyCOT.io.functions import read_txt
from pyCOT.simulations.ode import *
from pyCOT.visualization.plot_dynamics import *
from pyCOT.visualization.plot_process_analysis import *      # Updated: uses refactored version
from pyCOT.analysis.process_analyzer import *   # New: explicit import of analysis functions

# ========================================
# 2. CREATING THE REACTION_NETWORK OBJECT
# ========================================
print("="*60)
print("LOADING REACTION NETWORK")
print("="*60)

# Choose one of the following reaction networks:
# file_path = 'Txt/autopoietic.txt'  
# file_path = 'Txt/autopoietic_d.txt'  
file_path = 'Txt/Lotka_Volterra.txt'
# file_path = 'Txt/Farm.txt'

rn = read_txt(file_path)

species = [specie.name for specie in rn.species()]
print(f"\n📊 Species Set = {species}")

reactions = [reaction.name() for reaction in rn.reactions()]
print(f"⚗️  Reactions Set = {reactions}\n")

S = rn.stoichiometry_matrix()
print("Stoichiometry Matrix S =")
print(S)
print(f"Shape of S = {S.shape}")
print()

# ========================================      
# 3. INITIAL CONDITIONS & PARAMETERS
# ========================================
print("="*60)
print("SIMULATION PARAMETERS")
print("="*60)

# ---- Autopoietic system ----
# x0 = [2, 2, 2] 
# rate_list = 'mak'  # Mass-action kinetics for all reactions
# spec_vector = [[.7], [.5], [1], [1], [1]]

# ---- Autopoietic_d system ----
# x0 = [2, 2, 2, 0.1] 
# rate_list = 'mak'
# spec_vector = [[.7], [.5], [1], [1], [1], [0.55], [0.55], [0.55]]
# # Note: Sv=[-0.8, 0.5, -1], stabilizes at 3 different concentrations

# ---- Lotka-Volterra system ----
x0 = [2, 3] 
rate_list = 'mak'  # Mass-action kinetics
spec_vector = [[.5], [.8], [0.1]]

# ---- Simulation time parameters ----
tf = 200        # Final time
n = 2000        # Number of steps

print(f"Initial conditions: x0 = {x0}")
print(f"Time span: t ∈ [0, {tf}]")
print(f"Number of steps: {n+1}")
print(f"Kinetics: {rate_list}")
print()

# ========================================
# 4. RUN SIMULATION
# ========================================
print("="*60)
print("RUNNING SIMULATION")
print("="*60)

time_series, flux_vector = simulation(
    rn, 
    rate=rate_list, 
    spec_vector=spec_vector, 
    x0=x0, 
    t_span=(0, tf), 
    n_steps=n+1
)

print("✓ Simulation completed successfully!")
print(f"\nTime series shape: {time_series.shape}")
print(f"Flux vector shape: {flux_vector.shape}")
print("\nTime Series (first 5 rows):")
print(time_series) #.head())
print("\nFlux Vector (first 5 rows):")
print(flux_vector) #.head())
print()

# ========================================
# 5. BASIC VISUALIZATION
# ========================================
print("="*60)
print("BASIC VISUALIZATION")
print("="*60)

print("Generating time series plot...")
plot_series_ode(
    time_series, 
    filename="time_series_plot.png",
    show_fig=True
)
print("✓ Time series plot generated")
print()


# ========================================
# 6. MULTI-SCALE ANALYSIS: Process Type Proportions Over Time
# ========================================
print("="*60)
print("MULTI-SCALE TEMPORAL ANALYSIS")
print("="*60)
print("Analyzing process proportions across different time scales...")
print("This will generate:")
print("  1. Histograms for each window size")
print("  2. Global proportion trends plot")
print("  3. Combined plots for optimal window")
print()

# Define window sizes to analyze
# These represent different temporal aggregation scales
window_sizes = list(range(1, 601, 50))  # [1, 51, 101, 151, ..., 551]
print(f"Window sizes to analyze: {len(window_sizes)} scales")
print(f"  Min window: {min(window_sizes)} steps")
print(f"  Max window: {max(window_sizes)} steps")
print()

# Run the comprehensive analysis
results_df = analyze_process_proportions_over_time(
    rn=rn,
    S=S,
    rate_list=rate_list,
    spec_vector=spec_vector,
    x0=x0,
    t_span=(0, tf),
    n_steps=n+1,
    window_sizes=window_sizes,
    save_path="./visualizations/process_classification/"
)

print("\n" + "="*60)
print("ANALYSIS RESULTS")
print("="*60)
print(results_df)
print()

# ========================================
# 7. SUMMARY AND INSIGHTS
# ========================================
print("="*60)
print("SUMMARY")
print("="*60)

# Find optimal window
max_idx = results_df["Proporción Total"].idxmax()
optimal_window = results_df.loc[max_idx, "Ventana"]
max_proportion = results_df.loc[max_idx, "Proporción Total"]
optimal_steps = results_df.loc[max_idx, "Número de pasos"]

print(f"\n🎯 Optimal Window Size: {optimal_window} steps")
print(f"   Maximum Cognitive Domain proportion: {100*max_proportion:.2f}%")
print(f"   Effective time steps: {int(optimal_steps)}")
print()

# Show breakdown at optimal window
if optimal_window != "Máximo":
    optimal_row = results_df.loc[max_idx]
    print(f"Breakdown at optimal window ({optimal_window} steps):")
    print(f"  Cognitive Domain: {100*optimal_row['Proporción CC']:.2f}%")
    print(f"  Stationary Mode:   {100*optimal_row['Proporción SM']:.2f}%")
    print(f"  Problem:           {100*optimal_row['Proporción PB']:.2f}%")
    print()

print("📁 All visualizations saved to: ./visualizations/process_classification/")
print("   - histogramas_subplots.png    (All window histograms)")
print("   - proportions_plots.png       (Proportion trends)")
print("   - combined_plots.png          (Optimal window analysis)")
print()

print("="*60)
print("✅ SCRIPT COMPLETED SUCCESSFULLY")
print("="*60)

