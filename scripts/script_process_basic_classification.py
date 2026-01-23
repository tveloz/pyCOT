#!/usr/bin/env python3
"""
Script 1: Basic Process Classification Testing with Visualizations

Tests fundamental process classification capabilities with the refactored architecture.
Recovers functionality from old script_9.py with updated imports and adds new visualizations.

Tests:
- Stationary Mode, Cognitive Domain, Challenge, Problem
- Counteraction, Solution, Cognitive Control
- Disturbance-response analysis
- Multi-network comparison
- Edge cases

Visualizations:
- Classification summary plots
- Stoichiometric effect visualizations
- Disturbance response comparisons
- Multi-network classification comparison
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Add project root to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pyCOT.io.functions import read_txt
from pyCOT.simulations.ode import simulation
from pyCOT.process_analyzer import classify_process_mode, is_cognitive_domain, classify_response_to_disturbance
from pyCOT.semantic_partition import *
# Create output directory
output_dir = "./outputs/script_basic_process_classification"
os.makedirs(output_dir, exist_ok=True)

print("=" * 80)
print("BASIC PROCESS CLASSIFICATION TEST SUITE")
print("=" * 80)

# ============================================================================
# SECTION 1: LOAD AUTOPOIETIC NETWORK (PRIMARY TEST NETWORK)
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 1: Loading Autopoietic Network")
print("=" * 80)

file_path = 'Txt/autopoietic.txt'
rn = read_txt(file_path)

species = [specie.name for specie in rn.species()]
reactions = [reaction.name() for reaction in rn.reactions()]
S = rn.stoichiometry_matrix()

print(f"\nSpecies: {species}")
print(f"Reactions: {reactions}")
print(f"Stoichiometry matrix shape: {S.shape}")

# ============================================================================
# SECTION 2: TEST ALL CLASSIFICATION CATEGORIES
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 2: Testing All Classification Categories")
print("=" * 80)

def test_and_record(process_name, v, S, v_prev=None):
    """Test a process and record results for visualization."""
    result = {
        'name': process_name,
        'v': v,
        'Sv': S @ v,
        'mode': classify_process_mode(v, S),
        'is_cognitive': is_cognitive_domain(v, S)
    }
    
    if v_prev is not None:
        result['v_prev'] = v_prev
        result['Sv_prev'] = S @ v_prev
        result['Sv_combined'] = S @ (v + v_prev)
        result['disturbance_class'] = classify_response_to_disturbance(v, S, v_prev, verbose=False)
    
    # Print results
    print(f"\n{process_name}:")
    print(f"  Vector: {v}")
    print(f"  Stoichiometric effect (Sv): {result['Sv']}")
    print(f"  Classification: {result['mode']}")
    print(f"  In Cognitive Domain: {result['is_cognitive']}")
    
    if v_prev is not None:
        print(f"  Disturbance (v_prev): {v_prev}")
        print(f"  Disturbance effect (Sv_prev): {result['Sv_prev']}")
        print(f"  Combined effect (S(v+v_prev)): {result['Sv_combined']}")
        print(f"  Disturbance classification: {result['disturbance_class']}")
    
    return result

# Define test cases (from original script_9.py)
test_cases = []

# Test 1: Stationary Mode
v0 = np.array([2, 1, 2, 1, 1])
test_cases.append(test_and_record("v0 - Stationary Mode", v0, S))

# Test 2: Feasible (close to stationary)
v1 = np.array([4, 2, 1, 0, 0])
test_cases.append(test_and_record("v1 - Feasible", v1, S))

# Test 3: Cognitive Domain
v2 = np.array([3, 1, 2, 1, 1])
test_cases.append(test_and_record("v2 - Cognitive Domain", v2, S))

# Test 4: Challenge
vc3 = np.array([0, 0, 1, 0, 0])
test_cases.append(test_and_record("vc3 - Challenge", vc3, S))

# Test 5: Problem
vp4 = np.array([0, 0, 0, 0, 1])
test_cases.append(test_and_record("vp4 - Problem", vp4, S))

# Test 6: Counteraction
vc5 = np.array([2, 2, 1, 2, 1])
v1_counter = np.array([2, 0, 2.5, 0, 0])
test_cases.append(test_and_record("v5 - Counteraction", vc5, S, v_prev=v1_counter))

# Test 7: Solution
v6 = np.array([4, 2, 1, 0, 0])
vp6 = np.array([0, 0, 2, 2, 1])
test_cases.append(test_and_record("v6 - Solution", v6, S, v_prev=vp6))

# Test 8: Cognitive Control
v7 = np.array([4, 1, 2, 1, 1])
vp7 = np.array([0, 0, 1, 1, 0])
test_cases.append(test_and_record("v7 - Cognitive Control", v7, S, v_prev=vp7))

# ============================================================================
# SECTION 3: EDGE CASES
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 3: Testing Edge Cases")
print("=" * 80)

edge_cases = []

# Zero vector
v_zero = np.zeros(S.shape[1])
edge_cases.append(test_and_record("Zero vector", v_zero, S))

# All reactions equal
v_equal = np.ones(S.shape[1])
edge_cases.append(test_and_record("All reactions equal", v_equal, S))

# Single reaction
v_single = np.zeros(S.shape[1])
v_single[0] = 1.0
edge_cases.append(test_and_record("Single reaction (r1)", v_single, S))

# ============================================================================
# SECTION 4: VISUALIZATIONS
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 4: Generating Visualizations")
print("=" * 80)

# --- Visualization 1: Classification Summary ---
print("\nGenerating classification summary plot...")

fig, ax = plt.subplots(figsize=(12, 6))

# Extract primary classification for each test case
labels = [tc['name'].split(' - ')[0] for tc in test_cases]
primary_classes = [tc['mode'][0] if isinstance(tc['mode'], (list,tuple)) else tc['mode'] for tc in test_cases]

# Color mapping
color_map = {
    'Stationary Mode': 'blue',
    'Cognitive Domain': 'green',
    'Challenge': 'orange',
    'Problem': 'red',
    'Overproduction Mode': 'cyan',
    'Other': 'gray'
}

colors = [color_map.get(cls, 'gray') for cls in primary_classes]

# Create bar plot
bars = ax.bar(range(len(labels)), [1]*len(labels), color=colors, alpha=0.7, edgecolor='black')

# Add labels
ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels, rotation=45, ha='right')
ax.set_ylabel('Process Count')
ax.set_title('Process Classification Summary (Autopoietic Network)', fontsize=14, fontweight='bold')

# Add legend
legend_elements = [mpatches.Patch(facecolor=color, label=label) 
                  for label, color in color_map.items()]
ax.legend(handles=legend_elements, loc='upper right')

ax.set_ylim([0, 1.5])
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
save_path = os.path.join(output_dir, "classification_summary.png")
plt.savefig(save_path, dpi=300, bbox_inches='tight')
print(f"Saved: {save_path}")
plt.close()

# --- Visualization 2: Stoichiometric Effects ---
print("\nGenerating stoichiometric effects plot...")

fig, axes = plt.subplots(2, 4, figsize=(16, 8))
axes = axes.flatten()

for idx, tc in enumerate(test_cases[:8]):  # First 8 test cases
    ax = axes[idx]
    
    Sv = tc['Sv']
    species_names = species
    
    colors_bars = ['red' if x < 0 else 'green' if x > 0 else 'gray' for x in Sv]
    ax.bar(range(len(Sv)), Sv, color=colors_bars, alpha=0.7, edgecolor='black')
    
    ax.set_xticks(range(len(species_names)))
    ax.set_xticklabels(species_names)
    ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    ax.set_ylabel('Net Change (Sv)')
    ax.set_title(tc['name'].split(' - ')[1] if ' - ' in tc['name'] else tc['name'], 
                fontsize=10, fontweight='bold')
    ax.grid(axis='y', alpha=0.3)

plt.suptitle('Stoichiometric Effects for Each Process Type', fontsize=16, fontweight='bold')
plt.tight_layout()

save_path = os.path.join(output_dir, "stoichiometric_effects.png")
plt.savefig(save_path, dpi=300, bbox_inches='tight')
print(f"Saved: {save_path}")
plt.close()

# --- Visualization 3: Disturbance Response ---
print("\nGenerating disturbance response plot...")

# Filter test cases with disturbances
disturbance_cases = [tc for tc in test_cases if 'v_prev' in tc]

if disturbance_cases:
    fig, axes = plt.subplots(len(disturbance_cases), 3, figsize=(15, 4*len(disturbance_cases)))
    
    if len(disturbance_cases) == 1:
        axes = axes.reshape(1, -1)
    
    for idx, tc in enumerate(disturbance_cases):
        # Plot 1: Disturbance effect
        ax1 = axes[idx, 0]
        Sv_prev = tc['Sv_prev']
        colors_bars = ['red' if x < 0 else 'green' if x > 0 else 'gray' for x in Sv_prev]
        ax1.bar(range(len(Sv_prev)), Sv_prev, color=colors_bars, alpha=0.7, edgecolor='black')
        ax1.set_xticks(range(len(species)))
        ax1.set_xticklabels(species)
        ax1.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
        ax1.set_ylabel('Net Change')
        ax1.set_title(f'{tc["name"]} - Disturbance', fontsize=10, fontweight='bold')
        ax1.grid(axis='y', alpha=0.3)
        
        # Plot 2: Response effect
        ax2 = axes[idx, 1]
        Sv = tc['Sv']
        colors_bars = ['red' if x < 0 else 'green' if x > 0 else 'gray' for x in Sv]
        ax2.bar(range(len(Sv)), Sv, color=colors_bars, alpha=0.7, edgecolor='black')
        ax2.set_xticks(range(len(species)))
        ax2.set_xticklabels(species)
        ax2.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
        ax2.set_ylabel('Net Change')
        ax2.set_title(f'{tc["name"]} - Response', fontsize=10, fontweight='bold')
        ax2.grid(axis='y', alpha=0.3)
        
        # Plot 3: Combined effect
        ax3 = axes[idx, 2]
        Sv_combined = tc['Sv_combined']
        colors_bars = ['red' if x < 0 else 'green' if x > 0 else 'gray' for x in Sv_combined]
        ax3.bar(range(len(Sv_combined)), Sv_combined, color=colors_bars, alpha=0.7, edgecolor='black')
        ax3.set_xticks(range(len(species)))
        ax3.set_xticklabels(species)
        ax3.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
        ax3.set_ylabel('Net Change')
        ax3.set_title(f'{tc["name"]} - Combined', fontsize=10, fontweight='bold')
        ax3.grid(axis='y', alpha=0.3)
    
    plt.suptitle('Disturbance-Response Analysis', fontsize=16, fontweight='bold')
    plt.tight_layout()
    
    save_path = os.path.join(output_dir, "disturbance_response.png")
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {save_path}")
    plt.close()

# --- Visualization 4: Edge Cases ---
print("\nGenerating edge cases plot...")

fig, axes = plt.subplots(1, len(edge_cases), figsize=(15, 4))

for idx, tc in enumerate(edge_cases):
    ax = axes[idx]
    
    Sv = tc['Sv']
    colors_bars = ['red' if x < 0 else 'green' if x > 0 else 'gray' for x in Sv]
    ax.bar(range(len(Sv)), Sv, color=colors_bars, alpha=0.7, edgecolor='black')
    
    ax.set_xticks(range(len(species)))
    ax.set_xticklabels(species)
    ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    ax.set_ylabel('Net Change (Sv)')
    ax.set_title(tc['name'], fontsize=10, fontweight='bold')
    ax.grid(axis='y', alpha=0.3)
    
    # Add classification as text
    mode_text = tc['mode'][0] if isinstance(tc['mode'], (list, tuple)) else tc['mode']
    ax.text(0.5, 0.95, mode_text, transform=ax.transAxes,
           ha='center', va='top', fontsize=8, 
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.suptitle('Edge Case Analysis', fontsize=16, fontweight='bold')
plt.tight_layout()

save_path = os.path.join(output_dir, "edge_cases.png")
plt.savefig(save_path, dpi=300, bbox_inches='tight')
print(f"Saved: {save_path}")
plt.close()

# ============================================================================
# SECTION 5: MULTI-NETWORK COMPARISON
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 5: Multi-Network Comparison")
print("=" * 80)

# Test on Lotka-Volterra
print("\nTesting on Lotka-Volterra network...")
rn_lv = read_txt('Txt/Lotka_Volterra.txt')
S_lv = rn_lv.stoichiometry_matrix()

# Create comparable test processes
v_lv_stationary = np.array([1, 1, 0])  # Balanced
v_lv_challenge = np.array([1, 0, 0])   # Produces prey, challenge for predator
v_lv_problem = np.array([0, 0, 1])     # Death only

lv_results = [
    {'name': 'LV Stationary', 'mode': classify_process_mode(v_lv_stationary, S_lv)},
    {'name': 'LV Challenge', 'mode': classify_process_mode(v_lv_challenge, S_lv)},
    {'name': 'LV Problem', 'mode': classify_process_mode(v_lv_problem, S_lv)}
]

print("\nLotka-Volterra classifications:")
for result in lv_results:
    print(f"  {result['name']}: {result['mode']}")

# Visualization: Multi-network comparison
print("\nGenerating multi-network comparison plot...")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Autopoietic distribution
auto_classes = [tc['mode'][0] if isinstance(tc['mode'], (list, tuple)) else tc['mode'] 
                for tc in test_cases]
from collections import Counter
auto_counts = Counter(auto_classes)

ax1.bar(auto_counts.keys(), auto_counts.values(), color='steelblue', alpha=0.7, edgecolor='black')
ax1.set_ylabel('Count')
ax1.set_title('Autopoietic Network', fontsize=12, fontweight='bold')
ax1.tick_params(axis='x', rotation=45)
ax1.grid(axis='y', alpha=0.3)

# Lotka-Volterra distribution
lv_classes = [r['mode'][0] if isinstance(r['mode'], (list, tuple)) else r['mode'] for r in lv_results]
lv_counts = Counter(lv_classes)

ax2.bar(lv_counts.keys(), lv_counts.values(), color='coral', alpha=0.7, edgecolor='black')
ax2.set_ylabel('Count')
ax2.set_title('Lotka-Volterra Network', fontsize=12, fontweight='bold')
ax2.tick_params(axis='x', rotation=45)
ax2.grid(axis='y', alpha=0.3)

plt.suptitle('Process Classification Distribution Across Networks', fontsize=14, fontweight='bold')
plt.tight_layout()

save_path = os.path.join(output_dir, "multi_network_comparison.png")
plt.savefig(save_path, dpi=300, bbox_inches='tight')
print(f"Saved: {save_path}")
plt.close()

# ============================================================================
# SECTION 6: SUMMARY
# ============================================================================
print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)

print(f"\nTotal test cases: {len(test_cases) + len(edge_cases)}")
print(f"  Standard cases: {len(test_cases)}")
print(f"  Edge cases: {len(edge_cases)}")
print(f"  Disturbance-response cases: {len(disturbance_cases)}")

print(f"\nClassification breakdown (Autopoietic):")
for mode, count in auto_counts.items():
    print(f"  {mode}: {count}")

print(f"\nAll visualizations saved to: {output_dir}/")
print("  - classification_summary.png")
print("  - stoichiometric_effects.png")
print("  - disturbance_response.png")
print("  - edge_cases.png")
print("  - multi_network_comparison.png")

print("\n" + "=" * 80)
print("✅ SCRIPT COMPLETED SUCCESSFULLY")
print("=" * 80)
