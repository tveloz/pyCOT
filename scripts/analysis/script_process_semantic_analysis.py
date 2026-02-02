#!/usr/bin/env python3
"""
Script 3: Semantic Analysis Testing with Visualizations

Tests semantic interpretation and category-based analysis.
Recovers functionality from script_semantic_process.py with updated function names
and comprehensive category visualizations.

Tests:
- Category behavior analysis
- Single vs multi-reaction processes
- Pattern recognition (balanced transformation, reinforcing cycles, etc.)
- Observer-dependent classification
- Decomposition into single reactions

Visualizations:
- Category net effect bar charts
- Balance state indicators
- Pattern detection summaries
- Observer-dependent classification comparisons
- Decomposition breakdowns
- Semantic effect flow diagrams
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Add project root to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pyCOT.io.functions import read_txt
from pyCOT.core.semantic_partition import define_semantic_categories
from pyCOT.analysis.process_structure import decompose_to_single_reactions
from pyCOT.analysis.process_analyzer import (
    classify_process_mode,
    analyze_category_behavior,
    track_category_sequence,
    compute_category_statistics,
    identify_semantic_patterns
)

# Create output directory
output_dir = "./outputs/script_semantic_analysis"
os.makedirs(output_dir, exist_ok=True)

print("=" * 80)
print("SEMANTIC ANALYSIS TEST SUITE")
print("=" * 80)

# ============================================================================
# SECTION 1: LOAD CONFLICT MODEL AND DEFINE SEMANTICS
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 1: Loading Conflict Model")
print("=" * 80)

file_path = 'networks/Conflict_Theory/cause_driven_conflict_gov.txt'
rn = read_txt(file_path)

species_list = [s.name for s in rn.species()]
S = rn.stoichiometry_matrix()
n_reactions = len(rn.reactions())

print(f"\nLoaded model: {len(species_list)} species, {n_reactions} reactions")
print(f"Species: {species_list[:10]}..." if len(species_list) > 10 else f"Species: {species_list}")

# Create semantic partition
category_definitions = {
    'conflict': ['A_v', 'B_v', 'D_A', 'D_B', 'iD_A', 'iD_B'],
    'peace': ['A_p', 'B_p', 'G_A', 'G_B', 'iG_A', 'iG_B']
}

semantic_partition = define_semantic_categories(species_list, category_definitions)
print(f"\nSemantic partition created: {len(semantic_partition.categories)} categories")
for cat in semantic_partition.categories:
    indices = semantic_partition.category_indices[cat]
    print(f"  {cat}: {len(indices)} species")

# ============================================================================
# SECTION 2: SINGLE REACTION ANALYSIS
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 2: Single Reaction Analysis")
print("=" * 80)

# Test r7
print("\nAnalyzing reaction r7:")
print(f"  {rn.get_reaction('r7')}")

v_r7 = np.zeros(n_reactions)
v_r7[6] = 1.0  # r7 is index 6

# Analyze category behavior
behavior_r7 = analyze_category_behavior(v_r7, S, semantic_partition)

print(f"\nCategory behavior for r7:")
for cat, behavior in behavior_r7.items():
    print(f"\n  {cat.upper()}:")
    print(f"    Net effect: {behavior.net_effect:.4f}")
    print(f"    Production: {behavior.total_production:.4f}")
    print(f"    Consumption: {behavior.total_consumption:.4f}")
    print(f"    Balance state: {behavior.balance_state}")
    print(f"    Balance quality: {behavior.balance_quality:.4f}")

# Identify patterns
patterns_r7 = identify_semantic_patterns(behavior_r7)
print(f"\nPatterns identified in r7:")
for pattern in patterns_r7:
    print(f"  {pattern.pattern_type}: {pattern.description}")
    print(f"    Strength: {pattern.strength:.3f}")

# ============================================================================
# SECTION 3: MULTI-REACTION PROCESS ANALYSIS
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 3: Multi-Reaction Process Analysis")
print("=" * 80)

# Create a process with multiple reactions
v_multi = np.zeros(n_reactions)
v_multi[6] = 0.8   # r7
v_multi[10] = 1.2  # r11
v_multi[15] = 0.5  # r16

print(f"\nMulti-reaction process with {np.count_nonzero(v_multi)} active reactions")
print(f"  Active reactions: {np.where(v_multi > 0)[0] + 1}")

# Analyze category behavior
behavior_multi = analyze_category_behavior(v_multi, S, semantic_partition)

print(f"\nCategory behavior for multi-reaction process:")
for cat, behavior in behavior_multi.items():
    print(f"\n  {cat.upper()}:")
    print(f"    Net effect: {behavior.net_effect:.4f}")
    print(f"    Production: {behavior.total_production:.4f}")
    print(f"    Consumption: {behavior.total_consumption:.4f}")
    print(f"    Balance state: {behavior.balance_state}")

# Identify patterns
patterns_multi = identify_semantic_patterns(behavior_multi)
print(f"\nPatterns identified:")
for pattern in patterns_multi:
    print(f"  {pattern.pattern_type}: {pattern.description}")

# ============================================================================
# SECTION 4: PROCESS DECOMPOSITION
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 4: Process Decomposition")
print("=" * 80)

# Decompose multi-reaction process
print("\nDecomposing multi-reaction process into single reactions...")
components = decompose_to_single_reactions(v_multi)

# Debug: Check structure of first component
if len(components) > 0:
    first_comp = components[0]
    print(f"\nDEBUG - Component structure:")
    print(f"  Type: {type(first_comp)}")
    if isinstance(first_comp, dict):
        print(f"  Keys: {first_comp.keys()}")
    else:
        print(f"  Attributes: {dir(first_comp)}")

print(f"\nDecomposed into {len(components)} single reactions:")
for i, comp in enumerate(components):
    # Handle both dict and object attribute access
    if isinstance(comp, dict):
        reaction_idx = comp['reaction_idx']
        rate = comp['rate']
    else:
        # Assume it's an object with attributes
        reaction_idx = getattr(comp, 'reaction_idx', i)
        rate = getattr(comp, 'rate', 1.0)
    
    print(f"  {i+1}. Reaction r{reaction_idx+1}, rate={rate:.3f}")

# Verify decomposition
Sv_original = S @ v_multi

# Handle different component structures for reconstruction
if isinstance(components[0], dict):
    Sv_reconstructed = sum(comp['Sv_component'] for comp in components)
else:
    # If components are ProcessDisaggregation or similar objects
    # Reconstruct from individual process vectors
    Sv_reconstructed = np.zeros_like(Sv_original)
    for comp in components:
        if hasattr(comp, 'process'):
            Sv_reconstructed += S @ comp.process
        elif hasattr(comp, 'component'):
            Sv_reconstructed += S @ comp.component
        else:
            # Fallback: reconstruct from reaction_idx if available
            v_comp = np.zeros_like(v_multi)
            if hasattr(comp, 'reaction_idx'):
                v_comp[comp.reaction_idx] = getattr(comp, 'rate', 1.0)
                Sv_reconstructed += S @ v_comp

difference = np.linalg.norm(Sv_original - Sv_reconstructed)

print(f"\nDecomposition verification:")
print(f"  Original Sv norm: {np.linalg.norm(Sv_original):.6f}")
print(f"  Reconstructed Sv norm: {np.linalg.norm(Sv_reconstructed):.6f}")
print(f"  Difference: {difference:.6e}")
print(f"  Valid: {difference < 1e-6}")

# ============================================================================
# SECTION 5: OBSERVER-DEPENDENT CLASSIFICATION
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 5: Observer-Dependent Classification")
print("=" * 80)

# Create a test process
np.random.seed(42)
v_test = np.random.uniform(0, 1, n_reactions)

print(f"\nTest process with {np.count_nonzero(v_test > 0.01)} significant reactions")

# Get species indices for each category
conflict_indices = semantic_partition.category_indices['conflict']
peace_indices = semantic_partition.category_indices['peace']

# Create restricted stoichiometric matrices for each observer
S_conflict = S[conflict_indices, :]  # Only conflict species rows
S_peace = S[peace_indices, :]        # Only peace species rows

print(f"\nMatrix dimensions:")
print(f"  Full S: {S.shape}")
print(f"  Conflict view S: {S_conflict.shape} ({len(conflict_indices)} conflict species)")
print(f"  Peace view S: {S_peace.shape} ({len(peace_indices)} peace species)")

# Classify from different perspectives
# Each observer only sees effects on their species of interest
mode_full = classify_process_mode(v_test, S)
mode_conflict = classify_process_mode(v_test, S_conflict)
mode_peace = classify_process_mode(v_test, S_peace)

print(f"\nClassification from different perspectives:")
print(f"  Full network: {mode_full[0] if isinstance(mode_full, (list, tuple)) else mode_full}")
print(f"  Conflict observer (only sees conflict species): {mode_conflict[0] if isinstance(mode_conflict, (list, tuple)) else mode_conflict}")
print(f"  Peace observer (only sees peace species): {mode_peace[0] if isinstance(mode_peace, (list, tuple)) else mode_peace}")

# Show Sv effects for each subset
Sv_test = S @ v_test
Sv_conflict = S_conflict @ v_test
Sv_peace = S_peace @ v_test

print(f"\nStoichiometric effects visible to each observer:")
print(f"  Conflict observer sees Sv range: [{Sv_conflict.min():.3f}, {Sv_conflict.max():.3f}]")
print(f"  Peace observer sees Sv range: [{Sv_peace.min():.3f}, {Sv_peace.max():.3f}]")
print(f"  Full network Sv range: [{Sv_test.min():.3f}, {Sv_test.max():.3f}]")

# Analyze if observers would reach different conclusions
print(f"\nObserver agreement:")
if mode_full == mode_conflict == mode_peace:
    print("  All observers agree on classification")
else:
    print("  Observers DISAGREE - same process classified differently!")
    print("  This demonstrates observer-dependent interpretation")

# ============================================================================
# SECTION 6: VISUALIZATION - CATEGORY NET EFFECTS
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 6: Category Net Effects Visualization")
print("=" * 80)

print("\nGenerating category net effects plots...")

# Compare single vs multi-reaction
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Single reaction (r7)
categories = list(behavior_r7.keys())
net_effects_r7 = [behavior_r7[cat].net_effect for cat in categories]
production_r7 = [behavior_r7[cat].total_production for cat in categories]
consumption_r7 = [behavior_r7[cat].total_consumption for cat in categories]

x = np.arange(len(categories))
width = 0.25

ax1.bar(x - width, production_r7, width, label='Production', color='green', alpha=0.7)
ax1.bar(x, consumption_r7, width, label='Consumption', color='red', alpha=0.7)
ax1.bar(x + width, net_effects_r7, width, label='Net Effect', color='blue', alpha=0.7)

ax1.set_xlabel('Category')
ax1.set_ylabel('Effect Magnitude')
ax1.set_title('Single Reaction (r7)', fontsize=12, fontweight='bold')
ax1.set_xticks(x)
ax1.set_xticklabels(categories)
ax1.legend()
ax1.grid(axis='y', alpha=0.3)
ax1.axhline(y=0, color='black', linestyle='-', linewidth=0.5)

# Multi-reaction process
net_effects_multi = [behavior_multi[cat].net_effect for cat in categories]
production_multi = [behavior_multi[cat].total_production for cat in categories]
consumption_multi = [behavior_multi[cat].total_consumption for cat in categories]

ax2.bar(x - width, production_multi, width, label='Production', color='green', alpha=0.7)
ax2.bar(x, consumption_multi, width, label='Consumption', color='red', alpha=0.7)
ax2.bar(x + width, net_effects_multi, width, label='Net Effect', color='blue', alpha=0.7)

ax2.set_xlabel('Category')
ax2.set_ylabel('Effect Magnitude')
ax2.set_title('Multi-Reaction Process', fontsize=12, fontweight='bold')
ax2.set_xticks(x)
ax2.set_xticklabels(categories)
ax2.legend()
ax2.grid(axis='y', alpha=0.3)
ax2.axhline(y=0, color='black', linestyle='-', linewidth=0.5)

plt.suptitle('Category Production/Consumption/Net Effects', fontsize=14, fontweight='bold')
plt.tight_layout()

save_path = os.path.join(output_dir, "category_net_effects.png")
plt.savefig(save_path, dpi=300, bbox_inches='tight')
print(f"Saved: {save_path}")
plt.close()

# ============================================================================
# SECTION 7: VISUALIZATION - BALANCE STATES
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 7: Balance State Visualization")
print("=" * 80)

print("\nGenerating balance state visualization...")

# Create balance state visualization for multiple processes
test_processes = {
    'r7 (single)': v_r7,
    'Multi-reaction': v_multi,
    'Random': v_test
}

balance_data = {}
for name, v in test_processes.items():
    behavior = analyze_category_behavior(v, S, semantic_partition)
    balance_data[name] = {cat: behavior[cat].balance_state for cat in categories}

# Create heatmap-style visualization
fig, ax = plt.subplots(figsize=(10, 6))

# Balance state to numeric mapping
state_map = {'balanced': 0, 'overproduced': 1, 'depleted': -1}
color_map_balance = {0: 'green', 1: 'blue', -1: 'red'}

process_names = list(balance_data.keys())
n_processes = len(process_names)
n_cats = len(categories)

for i, process in enumerate(process_names):
    for j, cat in enumerate(categories):
        state = balance_data[process][cat]
        numeric_state = state_map[state]
        color = color_map_balance[numeric_state]
        
        rect = mpatches.Rectangle((j, i), 1, 1, 
                                  facecolor=color, edgecolor='black', alpha=0.7)
        ax.add_patch(rect)
        
        # Add text label
        ax.text(j + 0.5, i + 0.5, state.upper()[:3], 
               ha='center', va='center', fontsize=9, fontweight='bold')

ax.set_xlim(0, n_cats)
ax.set_ylim(0, n_processes)
ax.set_xticks(np.arange(n_cats) + 0.5)
ax.set_yticks(np.arange(n_processes) + 0.5)
ax.set_xticklabels(categories)
ax.set_yticklabels(process_names)
ax.set_xlabel('Category')
ax.set_ylabel('Process')
ax.set_title('Balance States: Balanced (Green) / Overproduced (Blue) / Depleted (Red)', 
            fontsize=12, fontweight='bold')

plt.tight_layout()

save_path = os.path.join(output_dir, "balance_states.png")
plt.savefig(save_path, dpi=300, bbox_inches='tight')
print(f"Saved: {save_path}")
plt.close()

# ============================================================================
# SECTION 8: VISUALIZATION - PATTERN DETECTION
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 8: Pattern Detection Visualization")
print("=" * 80)

print("\nGenerating pattern detection visualization...")

# Analyze patterns for all test processes
pattern_data = {}
for name, v in test_processes.items():
    behavior = analyze_category_behavior(v, S, semantic_partition)
    patterns = identify_semantic_patterns(behavior)
    pattern_data[name] = patterns

# Create pattern summary plot
fig, ax = plt.subplots(figsize=(12, 6))

# Collect all unique pattern types
all_pattern_types = set()
for patterns in pattern_data.values():
    for pattern in patterns:
        all_pattern_types.add(pattern.pattern_type)

pattern_types = sorted(list(all_pattern_types))
n_patterns = len(pattern_types)

# Create matrix showing which patterns appear in which processes
x_pos = np.arange(n_patterns)
bar_width = 0.25

for i, (process_name, patterns) in enumerate(pattern_data.items()):
    pattern_strengths = []
    for ptype in pattern_types:
        # Find pattern strength
        strength = 0
        for pattern in patterns:
            if pattern.pattern_type == ptype:
                strength = pattern.strength
                break
        pattern_strengths.append(strength)
    
    offset = (i - len(pattern_data)/2 + 0.5) * bar_width
    ax.bar(x_pos + offset, pattern_strengths, bar_width, 
          label=process_name, alpha=0.7, edgecolor='black')

ax.set_xlabel('Pattern Type')
ax.set_ylabel('Pattern Strength')
ax.set_title('Semantic Pattern Detection Across Processes', fontsize=14, fontweight='bold')
ax.set_xticks(x_pos)
ax.set_xticklabels(pattern_types, rotation=45, ha='right')
ax.legend()
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()

save_path = os.path.join(output_dir, "pattern_detection.png")
plt.savefig(save_path, dpi=300, bbox_inches='tight')
print(f"Saved: {save_path}")
plt.close()

# ============================================================================
# SECTION 9: VISUALIZATION - OBSERVER DEPENDENCE
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 9: Observer-Dependent Classification Visualization")
print("=" * 80)

print("\nGenerating observer-dependent classification plot...")

# Test multiple random processes
np.random.seed(123)
n_test = 10
observer_results = []

for i in range(n_test):
    v = np.random.uniform(0, 1, n_reactions)
    
    # Each observer sees different species (restricted stoichiometric matrix)
    mode_full = classify_process_mode(v, S)
    mode_conflict = classify_process_mode(v, S_conflict)  # Only conflict species
    mode_peace = classify_process_mode(v, S_peace)        # Only peace species
    
    # Extract first element if tuple/list
    if isinstance(mode_full, (list, tuple)):
        mode_full = mode_full[0]
    if isinstance(mode_conflict, (list, tuple)):
        mode_conflict = mode_conflict[0]
    if isinstance(mode_peace, (list, tuple)):
        mode_peace = mode_peace[0]
    
    observer_results.append({
        'full': mode_full,
        'conflict': mode_conflict,
        'peace': mode_peace
    })

# Count disagreements
disagreements = sum(1 for r in observer_results 
                   if not (r['full'] == r['conflict'] == r['peace']))
print(f"\nOut of {n_test} test processes:")
print(f"  {disagreements} showed observer-dependent classification")
print(f"  {n_test - disagreements} had agreement across all observers")

# Create comparison plot
fig, ax = plt.subplots(figsize=(14, 8))

# Mode to color mapping
mode_colors = {
    'Cognitive Domain': 'green',
    'Stationary Mode': 'blue',
    'Challenge': 'orange',
    'Problem': 'red',
    'Overproduction Mode': 'cyan',
    'Counteraction': 'purple',
    'Solution': 'lime',
    'Other': 'gray'
}

perspectives = ['Full Network', 'Conflict View', 'Peace View']
n_perspectives = len(perspectives)

for i, perspective_key in enumerate(['full', 'conflict', 'peace']):
    for j, result in enumerate(observer_results):
        mode = result[perspective_key]
        color = mode_colors.get(mode, 'gray')
        
        rect = mpatches.Rectangle((i, j), 1, 1, 
                                  facecolor=color, edgecolor='black', alpha=0.7)
        ax.add_patch(rect)
        
        # Add text label (abbreviated)
        mode_abbrev = ''.join([word[0] for word in mode.split()])
        ax.text(i + 0.5, j + 0.5, mode_abbrev, 
               ha='center', va='center', fontsize=8, fontweight='bold')

ax.set_xlim(0, n_perspectives)
ax.set_ylim(0, n_test)
ax.set_xticks(np.arange(n_perspectives) + 0.5)
ax.set_yticks(np.arange(n_test) + 0.5)
ax.set_xticklabels(perspectives)
ax.set_yticklabels([f'Process {i+1}' for i in range(n_test)])
ax.set_xlabel('Observer Perspective')
ax.set_ylabel('Test Process')
ax.set_title('Observer-Dependent Process Classification\n(Same process, different perspectives)', 
            fontsize=14, fontweight='bold')

# Add legend
legend_elements = [mpatches.Patch(facecolor=color, label=mode) 
                  for mode, color in mode_colors.items()]
ax.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5))

plt.tight_layout()

save_path = os.path.join(output_dir, "observer_dependence.png")
plt.savefig(save_path, dpi=300, bbox_inches='tight')
print(f"Saved: {save_path}")
plt.close()

# ============================================================================
# SECTION 10: VISUALIZATION - DECOMPOSITION BREAKDOWN
# ============================================================================
print("\n" + "=" * 80)
print("SECTION 10: Decomposition Breakdown Visualization")
print("=" * 80)

print("\nGenerating decomposition breakdown plot...")

# Analyze each component's contribution to categories
component_effects = []
for comp in components:
    # Handle different data structures
    if isinstance(comp, dict):
        reaction_idx = comp['reaction_idx']
        rate = comp['rate']
    else:
        reaction_idx = getattr(comp, 'reaction_idx', 0)
        rate = getattr(comp, 'rate', 1.0)
    
    v_comp = np.zeros(n_reactions)
    v_comp[reaction_idx] = rate
    
    behavior = analyze_category_behavior(v_comp, S, semantic_partition)
    effects = {cat: behavior[cat].net_effect for cat in categories}
    effects['reaction'] = f"r{reaction_idx+1}"
    component_effects.append(effects)

# Create stacked bar chart
fig, ax = plt.subplots(figsize=(12, 6))

# Prepare data for stacking
n_components = len(component_effects)
category_data = {cat: [] for cat in categories}

for comp_effect in component_effects:
    for cat in categories:
        category_data[cat].append(comp_effect[cat])

# Create stacked bars
x = np.arange(n_components)
bottom = np.zeros(n_components)

category_colors = {'conflict': 'red', 'peace': 'green'}

for cat in categories:
    values = category_data[cat]
    ax.bar(x, values, bottom=bottom, label=cat.capitalize(), 
          color=category_colors.get(cat, 'gray'), alpha=0.7, edgecolor='black')
    bottom += np.array(values)

# Add original process total as line
total_behavior = analyze_category_behavior(v_multi, S, semantic_partition)
for i, cat in enumerate(categories):
    total = total_behavior[cat].net_effect
    ax.axhline(y=total, color=category_colors.get(cat, 'gray'), 
              linestyle='--', linewidth=2, alpha=0.5)

ax.set_xlabel('Component Reaction')
ax.set_ylabel('Net Effect Contribution')
ax.set_title('Decomposition: Single Reaction Contributions to Category Effects', 
            fontsize=14, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels([comp['reaction'] for comp in component_effects])
ax.legend()
ax.grid(axis='y', alpha=0.3)
ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5)

plt.tight_layout()

save_path = os.path.join(output_dir, "decomposition_breakdown.png")
plt.savefig(save_path, dpi=300, bbox_inches='tight')
print(f"Saved: {save_path}")
plt.close()

# ============================================================================
# SECTION 11: SUMMARY
# ============================================================================
print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)

print(f"\nProcesses analyzed: {len(test_processes)}")
print(f"Categories defined: {len(categories)}")
print(f"Patterns identified across all processes:")
for ptype in pattern_types:
    count = sum(1 for patterns in pattern_data.values() 
               if any(p.pattern_type == ptype for p in patterns))
    print(f"  {ptype}: {count} processes")

print(f"\nObserver-dependent classification tests: {n_test}")
print(f"Decomposition verified: Difference = {difference:.6e}")

print(f"\nAll visualizations saved to: {output_dir}/")
print("  - category_net_effects.png")
print("  - balance_states.png")
print("  - pattern_detection.png")
print("  - observer_dependence.png")
print("  - decomposition_breakdown.png")

print("\n" + "=" * 80)
print("✅ SCRIPT COMPLETED SUCCESSFULLY")
print("=" * 80)
