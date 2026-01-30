#!/usr/bin/env python3
"""
Intervention Multiscale Semantic Analysis

Analyzes how intervention strategies operate across temporal scales using
the semantic multiscale architecture. Compares baseline vs intervention
to determine if strategies work through:
- Short-term accumulation (immediate effects that stack)
- Long-term adaptation (system reorganization over time)
- Mixed profiles

Usage:
1. Load strategy from CSV ranking file
2. OR manually set intervention parameters
3. Compare with baseline across multiple temporal scales
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.gridspec import GridSpec
from pyCOT.simulations.ode import simulation
from pyCOT.kinetics.deterministic_advanced import rate_saturated, rate_threshold_memory, rate_cosine
from pyCOT.io.functions import read_txt
from pyCOT.semantic_partition import define_semantic_categories
from pyCOT.process_analyzer import analyze_temporal_scale

# ========================================
# CONFIGURATION
# ========================================

# Model loading
FILE_PATH = 'networks/Conflict_Theory/Resource_Scarcity_Toy_Model2.txt'

# Semantic categories for conflict model
CATEGORY_DEFINITIONS = {
    'peace': ['SR', 'R', 'E', 'T'],      # Strong regions, Resources, Economy, Trust
    'conflict': ['DT', 'V', 'WR']        # Detached, Violence, Weak regions
}

# Temporal scales to analyze (in timesteps)
WINDOW_SIZES = [1, 5, 10, 20, 40]

# Simulation parameters
SIMULATION_TIME = 200
N_STEPS = 400
INITIAL_CONDITION = [1.0, 0.0, 0.0, 0.5, 1.5, 0.0, 0.0]  # [SR, R, E, WR, DT, T, V]

# ========================================
# BASELINE PARAMETERS (REFERENCE)
# ========================================

BASELINE_PARAMS = {
    'Vmax_production': 0.5,
    'Km_production': 0.5,
    'E_threshold': 0.1,
    'T_threshold_weak_strong': 1.5,
    'T_threshold_detached_weak': 1.5,
    'k_recovery': 0.1,
    'k_degradation': 0.1,
    'detachment_threshold': 1.5,
    'k_detachment_base': 5.0,
    'Vmax_trust': 0.1,
    'Km_trust_economy': 0.2,
    'k_trust_destruction': 0.01,
    'k_trust_decay': 0.01,
    'k_violence_decay': 0.01,
    'k_violence_weak': 0.01,
    'k_violence_detached_weak': 0.05,
    'k_violence_detached': 0.1,
    'R_amplitude': 0.5,
    'R_frequency': 1.0,
    'k_resource_depletion': 0.02,
    'k_economic_decay': 0.02
}

# ========================================
# MANUAL INTERVENTION EXAMPLE
# ========================================
# Modify these values to test custom interventions
# Set to None to use baseline, or provide modified value

MANUAL_INTERVENTION = {
    # Example: Boost trust generation and reduce violence
    'Vmax_trust': 0.15,              # +50% trust generation
    'k_violence_weak': 0.005,        # -50% violence from weak populations
    'k_trust_decay': 0.005,          # -50% trust decay
    # Leave others as None to use baseline
}

# ========================================
# HELPER FUNCTIONS
# ========================================

def build_spec_vector(params):
    """Build spec_vector from parameter dictionary"""
    return [
        [params['Vmax_production'], params['Km_production']],
        [params['E_threshold'], params['k_recovery']],
        [params['E_threshold'], params['k_recovery']],
        [params['T_threshold_weak_strong'], params['k_recovery']], 
        [params['T_threshold_detached_weak'], params['k_recovery']],
        [params['k_degradation']],
        [params['detachment_threshold'], params['k_detachment_base']],
        [params['Vmax_trust'], params['Km_trust_economy']],
        [params['k_trust_destruction']],
        [params['k_trust_decay']],
        [params['k_violence_weak']],
        [params['k_violence_detached_weak']],
        [params['k_violence_detached']],
        [params['k_violence_decay']],
        [params['R_amplitude'], params['R_frequency']],
        [params['k_resource_depletion']],
        [params['k_economic_decay']]
    ]

def load_strategy_from_csv(csv_path, strategy_index=0, strategy_type='single'):
    """
    Load intervention strategy from CSV ranking file
    
    Args:
        csv_path: Path to intervention_rankings.csv or intervention_combinations.csv
        strategy_index: Which ranked strategy to load (0 = best)
        strategy_type: 'single' or 'combination'
    
    Returns:
        params_dict: Modified parameters
        strategy_info: Dict with strategy metadata
    """
    df = pd.read_csv(csv_path)
    
    if strategy_index >= len(df):
        raise ValueError(f"Strategy index {strategy_index} out of range (max: {len(df)-1})")
    
    strategy = df.iloc[strategy_index]
    params = BASELINE_PARAMS.copy()
    
    if strategy_type == 'single':
        # Single intervention format: "parameter ↑/↓ magnitude%"
        param_name = strategy['Parameter']
        magnitude = strategy['Magnitude (%)'] / 100.0
        direction = strategy['Direction']
        
        if direction == '↑':
            params[param_name] = BASELINE_PARAMS[param_name] * (1 + magnitude)
        else:  # '↓'
            params[param_name] = BASELINE_PARAMS[param_name] * (1 - magnitude)
        
        strategy_info = {
            'name': strategy['Intervention'],
            'type': 'single',
            'success_rate': strategy['Success Rate (%)'],
            'cost': strategy['Cost'],
            'cost_effectiveness': strategy['Cost-Effectiveness']
        }
    
    else:  # combination
        # Parse combination format: "param1↑X% + param2↓Y% + ..."
        strategy_str = strategy['Strategy']
        params_list = strategy['Parameters'].split(', ')
        magnitudes_list = [float(m.strip('%')) / 100 for m in strategy['Magnitudes'].split(', ')]
        
        # Parse directions from strategy string
        for param_name, magnitude in zip(params_list, magnitudes_list):
            # Find direction from strategy string
            if param_name in strategy_str:
                if f"{param_name}↑" in strategy_str:
                    params[param_name] = BASELINE_PARAMS[param_name] * (1 + magnitude)
                else:
                    params[param_name] = BASELINE_PARAMS[param_name] * (1 - magnitude)
        
        strategy_info = {
            'name': strategy_str[:60] + '...' if len(strategy_str) > 60 else strategy_str,
            'type': 'combination',
            'success_rate': strategy['Success Rate (%)'],
            'cost': strategy['Total Cost'],
            'cost_effectiveness': strategy['Cost-Effectiveness']
        }
    
    return params, strategy_info

def apply_manual_intervention(manual_dict):
    """Apply manual intervention to baseline parameters"""
    params = BASELINE_PARAMS.copy()
    modified_params = []
    
    for key, value in manual_dict.items():
        if value is not None:
            params[key] = value
            pct_change = ((value / BASELINE_PARAMS[key]) - 1) * 100
            direction = "↑" if pct_change > 0 else "↓"
            modified_params.append(f"{key}{direction}{abs(pct_change):.1f}%")
    
    strategy_info = {
        'name': 'Manual: ' + ', '.join(modified_params),
        'type': 'manual',
        'success_rate': None,
        'cost': None,
        'cost_effectiveness': None
    }
    
    return params, strategy_info

def run_simulation_with_params(params, rn, rate_list, additional_laws):
    """Run simulation with given parameters"""
    spec_vector = build_spec_vector(params)
    
    ts, fv = simulation(
        rn,
        rate=rate_list,
        spec_vector=spec_vector,
        x0=INITIAL_CONDITION,
        t_span=(0, SIMULATION_TIME),
        n_steps=N_STEPS,
        additional_laws=additional_laws,
        verbose=False
    )
    
    # Extract process series
    process_series = [fv.iloc[i, 1:].values for i in range(len(fv))]
    
    return ts, fv, process_series

def classify_profile(scale_results, category):
    """
    Classify intervention profile for a category
    
    Returns: 'accumulative', 'adaptive', 'mixed', or 'neutral'
    """
    window_sizes = scale_results['scales']
    stats = [scale_results['results_by_scale'][ws]['category_statistics'][category] 
             for ws in window_sizes]
    
    net_effects = [s['net_effect']['mean'] for s in stats]
    
    # Calculate trend
    if len(net_effects) < 3:
        return 'neutral'
    
    # Linear regression to see if effect strengthens with scale
    x = np.log10(window_sizes)  # Log scale for window sizes
    y = net_effects
    slope = np.polyfit(x, y, 1)[0]
    
    # Check variance
    variance = np.std(net_effects)
    
    # Classification logic
    if variance < 0.05:  # Very stable across scales
        return 'accumulative'
    elif abs(slope) > 0.1:  # Strong trend with scale
        return 'adaptive' if slope != 0 else 'accumulative'
    else:
        return 'mixed'

# ========================================
# VISUALIZATION FUNCTIONS
# ========================================

def plot_category_comparison(baseline_results, intervention_results, strategy_info, 
                            category_name, save_path):
    """
    Plot category behavior comparison across scales
    
    Shows baseline vs intervention for one category
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    window_sizes = baseline_results['scales']
    
    # Extract data
    baseline_stats = [baseline_results['results_by_scale'][ws]['category_statistics'][category_name] 
                      for ws in window_sizes]
    intervention_stats = [intervention_results['results_by_scale'][ws]['category_statistics'][category_name] 
                         for ws in window_sizes]
    
    baseline_net = [s['net_effect']['mean'] for s in baseline_stats]
    baseline_std = [s['net_effect']['std'] for s in baseline_stats]
    
    intervention_net = [s['net_effect']['mean'] for s in intervention_stats]
    intervention_std = [s['net_effect']['std'] for s in intervention_stats]
    
    # Plot 1: Net effect across scales
    ax1.plot(window_sizes, baseline_net, 'o-', linewidth=2, markersize=8, 
             label='Baseline', color='gray', alpha=0.7)
    ax1.fill_between(window_sizes, 
                     np.array(baseline_net) - np.array(baseline_std),
                     np.array(baseline_net) + np.array(baseline_std),
                     alpha=0.2, color='gray')
    
    color = 'green' if category_name == 'peace' else 'red'
    ax1.plot(window_sizes, intervention_net, 'o-', linewidth=2, markersize=8,
             label='Intervention', color=color)
    ax1.fill_between(window_sizes,
                     np.array(intervention_net) - np.array(intervention_std),
                     np.array(intervention_net) + np.array(intervention_std),
                     alpha=0.2, color=color)
    
    ax1.axhline(y=0, color='black', linestyle='--', alpha=0.3)
    ax1.set_xscale('log')
    ax1.set_xlabel('Window Size (timesteps)', fontsize=11)
    ax1.set_ylabel('Mean Net Effect', fontsize=11)
    ax1.set_title(f'{category_name.upper()} - Net Effect Across Scales', 
                  fontsize=12, fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Difference (intervention - baseline)
    difference = np.array(intervention_net) - np.array(baseline_net)
    
    colors_diff = ['green' if d > 0 else 'red' for d in difference]
    if category_name == 'conflict':
        colors_diff = ['green' if d < 0 else 'red' for d in difference]  # Invert for conflict
    
    ax2.bar(range(len(window_sizes)), difference, color=colors_diff, alpha=0.6, edgecolor='black')
    ax2.axhline(y=0, color='black', linestyle='-', linewidth=1)
    ax2.set_xticks(range(len(window_sizes)))
    ax2.set_xticklabels(window_sizes)
    ax2.set_xlabel('Window Size (timesteps)', fontsize=11)
    ax2.set_ylabel('Δ Net Effect (Intervention - Baseline)', fontsize=11)
    ax2.set_title(f'{category_name.upper()} - Intervention Impact', 
                  fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3, axis='y')
    
    plt.suptitle(f'Strategy: {strategy_info["name"]}\n{category_name.upper()} Category Analysis',
                 fontsize=13, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {save_path}")
    plt.close()

def plot_comprehensive_comparison(baseline_results, intervention_results, strategy_info, save_path):
    """
    Comprehensive comparison plot showing all categories
    """
    fig = plt.figure(figsize=(16, 10))
    gs = GridSpec(3, 2, figure=fig, hspace=0.35, wspace=0.3)
    
    categories = list(baseline_results['results_by_scale'][WINDOW_SIZES[0]]['category_statistics'].keys())
    window_sizes = baseline_results['scales']
    
    # Determine which category is which
    peace_cat = 'peace' if 'peace' in categories else categories[0]
    conflict_cat = 'conflict' if 'conflict' in categories else categories[1]
    
    # Top row: Peace category
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    plot_category_on_axes(ax1, ax2, baseline_results, intervention_results, 
                         peace_cat, window_sizes, 'peace')
    
    # Middle row: Conflict category
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])
    plot_category_on_axes(ax3, ax4, baseline_results, intervention_results,
                         conflict_cat, window_sizes, 'conflict')
    
    # Bottom row: Balance comparison
    ax5 = fig.add_subplot(gs[2, :])
    plot_balance_comparison(ax5, baseline_results, intervention_results,
                           peace_cat, conflict_cat, window_sizes)
    
    # Title
    title_text = f'Multiscale Semantic Analysis\nStrategy: {strategy_info["name"]}'
    if strategy_info['success_rate'] is not None:
        title_text += f'\nSuccess Rate: {strategy_info["success_rate"]:.1f}% | Cost: ${strategy_info["cost"]:.0f}'
    
    plt.suptitle(title_text, fontsize=14, fontweight='bold', y=0.98)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {save_path}")
    plt.close()

def plot_category_on_axes(ax1, ax2, baseline_results, intervention_results, 
                         category_name, window_sizes, category_type):
    """Helper to plot category on given axes"""
    baseline_stats = [baseline_results['results_by_scale'][ws]['category_statistics'][category_name] 
                      for ws in window_sizes]
    intervention_stats = [intervention_results['results_by_scale'][ws]['category_statistics'][category_name] 
                         for ws in window_sizes]
    
    baseline_net = [s['net_effect']['mean'] for s in baseline_stats]
    baseline_std = [s['net_effect']['std'] for s in baseline_stats]
    intervention_net = [s['net_effect']['mean'] for s in intervention_stats]
    intervention_std = [s['net_effect']['std'] for s in intervention_stats]
    
    # Net effect plot
    ax1.plot(window_sizes, baseline_net, 'o-', linewidth=2, markersize=6,
             label='Baseline', color='gray', alpha=0.7)
    ax1.fill_between(window_sizes, 
                     np.array(baseline_net) - np.array(baseline_std),
                     np.array(baseline_net) + np.array(baseline_std),
                     alpha=0.2, color='gray')
    
    color = 'green' if category_type == 'peace' else 'red'
    ax1.plot(window_sizes, intervention_net, 'o-', linewidth=2, markersize=6,
             label='Intervention', color=color)
    ax1.fill_between(window_sizes,
                     np.array(intervention_net) - np.array(intervention_std),
                     np.array(intervention_net) + np.array(intervention_std),
                     alpha=0.2, color=color)
    
    ax1.axhline(y=0, color='black', linestyle='--', alpha=0.3)
    ax1.set_xscale('log')
    ax1.set_xlabel('Window Size', fontsize=10)
    ax1.set_ylabel('Net Effect', fontsize=10)
    ax1.set_title(f'{category_name.upper()} - Net Effect Across Scales', fontsize=11, fontweight='bold')
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3)
    
    # Difference plot
    difference = np.array(intervention_net) - np.array(baseline_net)
    colors_diff = ['green' if d > 0 else 'red' for d in difference]
    if category_type == 'conflict':
        colors_diff = ['green' if d < 0 else 'red' for d in difference]
    
    ax2.bar(range(len(window_sizes)), difference, color=colors_diff, alpha=0.6, edgecolor='black')
    ax2.axhline(y=0, color='black', linestyle='-', linewidth=1)
    ax2.set_xticks(range(len(window_sizes)))
    ax2.set_xticklabels(window_sizes, fontsize=9)
    ax2.set_xlabel('Window Size', fontsize=10)
    ax2.set_ylabel('Δ Net Effect', fontsize=10)
    ax2.set_title(f'{category_name.upper()} - Impact at Each Scale', fontsize=11, fontweight='bold')
    ax2.grid(True, alpha=0.3, axis='y')

def plot_balance_comparison(ax, baseline_results, intervention_results,
                           peace_cat, conflict_cat, window_sizes):
    """Plot peace/conflict balance across scales"""
    baseline_peace = [baseline_results['results_by_scale'][ws]['category_statistics'][peace_cat]['net_effect']['mean'] 
                      for ws in window_sizes]
    baseline_conflict = [baseline_results['results_by_scale'][ws]['category_statistics'][conflict_cat]['net_effect']['mean'] 
                        for ws in window_sizes]
    
    intervention_peace = [intervention_results['results_by_scale'][ws]['category_statistics'][peace_cat]['net_effect']['mean'] 
                         for ws in window_sizes]
    intervention_conflict = [intervention_results['results_by_scale'][ws]['category_statistics'][conflict_cat]['net_effect']['mean'] 
                            for ws in window_sizes]
    
    baseline_balance = np.array(baseline_peace) - np.array(baseline_conflict)
    intervention_balance = np.array(intervention_peace) - np.array(intervention_conflict)
    
    x = np.arange(len(window_sizes))
    width = 0.35
    
    ax.bar(x - width/2, baseline_balance, width, label='Baseline', color='gray', alpha=0.7)
    ax.bar(x + width/2, intervention_balance, width, label='Intervention', color='steelblue', alpha=0.7)
    
    ax.axhline(y=0, color='black', linestyle='-', linewidth=1)
    ax.set_xticks(x)
    ax.set_xticklabels(window_sizes)
    ax.set_xlabel('Window Size (timesteps)', fontsize=11)
    ax.set_ylabel('Balance (Peace - Conflict)', fontsize=11)
    ax.set_title('Overall Balance Across Scales', fontsize=12, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, axis='y')

def print_profile_analysis(baseline_results, intervention_results, categories):
    """Print profile classification for each category"""
    print("\n" + "="*80)
    print("PROFILE ANALYSIS")
    print("="*80)
    
    for category in categories:
        baseline_profile = classify_profile(baseline_results, category)
        intervention_profile = classify_profile(intervention_results, category)
        
        print(f"\n{category.upper()}:")
        print(f"  Baseline Profile:     {baseline_profile.upper()}")
        print(f"  Intervention Profile: {intervention_profile.upper()}")
        
        # Interpretation
        if intervention_profile == 'accumulative':
            print(f"  → Immediate effects that remain constant across scales")
        elif intervention_profile == 'adaptive':
            print(f"  → Effects strengthen/weaken with longer temporal windows")
        elif intervention_profile == 'mixed':
            print(f"  → Combination of immediate and scale-dependent effects")

# ========================================
# MAIN EXECUTION
# ========================================

def main():
    print("="*80)
    print("INTERVENTION MULTISCALE SEMANTIC ANALYSIS")
    print("="*80)
    
    # ========================================
    # SETUP
    # ========================================
    print("\n1. Loading model...")
    rn = read_txt(FILE_PATH)
    species_list = [s.name for s in rn.species()]
    S = rn.stoichiometry_matrix()
    
    semantic_partition = define_semantic_categories(species_list, CATEGORY_DEFINITIONS)
    print(f"   Species: {species_list}")
    print(f"   Semantic categories: {semantic_partition.categories}")
    
    additional_laws = {
        'saturated': rate_saturated,
        'threshold_memory': rate_threshold_memory,
        'cosine': rate_cosine
    }
    
    rate_list = [
        'saturated', 'threshold_memory', 'threshold_memory', 'threshold_memory',
        'threshold_memory', 'mak', 'threshold_memory', 'saturated',
        'mak', 'mak', 'mak', 'mak', 'mak', 'mak', 'cosine', 'mak', 'mak'
    ]
    
    # ========================================
    # CHOOSE INTERVENTION MODE
    # ========================================
    print("\n2. Intervention Selection:")
    print("-"*80)
    print("Choose mode:")
    print("  1. Load from CSV ranking file")
    print("  2. Use manual intervention (modify MANUAL_INTERVENTION in script)")
    print("  3. Use baseline only (no intervention)")
    
    mode = input("\nEnter choice (1/2/3) [default=2]: ").strip() or "2"
    
    if mode == "1":
        csv_path = input("Enter CSV path [intervention_rankings.csv]: ").strip() or "intervention_rankings.csv"
        strategy_type = input("Strategy type (single/combination) [single]: ").strip() or "single"
        strategy_index = int(input("Strategy index (0=best) [0]: ").strip() or "0")
        
        intervention_params, strategy_info = load_strategy_from_csv(
            csv_path, strategy_index, strategy_type
        )
    
    elif mode == "2":
        intervention_params, strategy_info = apply_manual_intervention(MANUAL_INTERVENTION)
    
    else:  # mode == "3"
        intervention_params = BASELINE_PARAMS.copy()
        strategy_info = {
            'name': 'Baseline (no intervention)',
            'type': 'baseline',
            'success_rate': None,
            'cost': None,
            'cost_effectiveness': None
        }
    
    print(f"\n   Selected: {strategy_info['name']}")
    if strategy_info['success_rate'] is not None:
        print(f"   Success Rate: {strategy_info['success_rate']:.1f}%")
        print(f"   Cost: ${strategy_info['cost']:.0f}")
    
    # ========================================
    # RUN SIMULATIONS
    # ========================================
    print("\n3. Running simulations...")
    
    print("   Running baseline simulation...")
    ts_baseline, fv_baseline, process_baseline = run_simulation_with_params(
        BASELINE_PARAMS, rn, rate_list, additional_laws
    )
    
    print("   Running intervention simulation...")
    ts_intervention, fv_intervention, process_intervention = run_simulation_with_params(
        intervention_params, rn, rate_list, additional_laws
    )
    
    # ========================================
    # MULTISCALE ANALYSIS
    # ========================================
    print("\n4. Performing multiscale analysis...")
    print(f"   Window sizes: {WINDOW_SIZES}")
    
    print("   Analyzing baseline...")
    baseline_results = analyze_temporal_scale(
        process_series=process_baseline,
        S=S,
        semantic_partition=semantic_partition,
        window_sizes=WINDOW_SIZES,
        stride=1
    )
    
    print("   Analyzing intervention...")
    intervention_results = analyze_temporal_scale(
        process_series=process_intervention,
        S=S,
        semantic_partition=semantic_partition,
        window_sizes=WINDOW_SIZES,
        stride=1
    )
    
    # ========================================
    # PROFILE ANALYSIS
    # ========================================
    categories = list(baseline_results['results_by_scale'][WINDOW_SIZES[0]]['category_statistics'].keys())
    print_profile_analysis(baseline_results, intervention_results, categories)
    
    # ========================================
    # NUMERICAL SUMMARY
    # ========================================
    print("\n" + "="*80)
    print("NUMERICAL SUMMARY")
    print("="*80)
    
    for category in categories:
        print(f"\n{category.upper()}:")
        print(f"  {'Window':<10} {'Baseline Net':<15} {'Intervention Net':<18} {'Difference':<12} {'% Change'}")
        print(f"  {'-'*75}")
        
        for ws in WINDOW_SIZES:
            b_stat = baseline_results['results_by_scale'][ws]['category_statistics'][category]
            i_stat = intervention_results['results_by_scale'][ws]['category_statistics'][category]
            
            b_net = b_stat['net_effect']['mean']
            i_net = i_stat['net_effect']['mean']
            diff = i_net - b_net
            
            if abs(b_net) > 1e-6:
                pct_change = (diff / abs(b_net)) * 100
            else:
                pct_change = float('inf') if diff != 0 else 0
            
            pct_str = f"{pct_change:+.1f}%" if abs(pct_change) < 999 else "N/A"
            
            print(f"  {ws:<10} {b_net:<+15.4f} {i_net:<+18.4f} {diff:<+12.4f} {pct_str}")
    
    # ========================================
    # VISUALIZATION
    # ========================================
    print("\n5. Generating visualizations...")
    
    output_dir = "./outputs/intervention_multiscale"
    os.makedirs(output_dir, exist_ok=True)
    
    # Comprehensive comparison
    plot_comprehensive_comparison(
        baseline_results, intervention_results, strategy_info,
        os.path.join(output_dir, "comprehensive_comparison.png")
    )
    
    # Individual category plots
    for category in categories:
        plot_category_comparison(
            baseline_results, intervention_results, strategy_info, category,
            os.path.join(output_dir, f"category_{category}_analysis.png")
        )
    
    # ========================================
    # TIME SERIES COMPARISON
    # ========================================
    print("\n6. Generating time series comparison...")
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Get species indices for each category
    peace_species = semantic_partition.category_to_indices['peace']
    conflict_species = semantic_partition.category_to_indices['conflict']
    
    time = ts_baseline['Time'].values
    
    # Peace species
    ax = axes[0, 0]
    for idx in peace_species:
        species_name = species_list[idx]
        ax.plot(time, ts_baseline.iloc[:, idx+1], '--', alpha=0.5, label=f'{species_name} (base)')
        ax.plot(time, ts_intervention.iloc[:, idx+1], '-', linewidth=2, label=f'{species_name} (int)')
    ax.set_xlabel('Time', fontsize=10)
    ax.set_ylabel('Concentration', fontsize=10)
    ax.set_title('Peace Species Dynamics', fontsize=11, fontweight='bold')
    ax.legend(fontsize=8, ncol=2)
    ax.grid(True, alpha=0.3)
    
    # Conflict species
    ax = axes[0, 1]
    for idx in conflict_species:
        species_name = species_list[idx]
        ax.plot(time, ts_baseline.iloc[:, idx+1], '--', alpha=0.5, label=f'{species_name} (base)')
        ax.plot(time, ts_intervention.iloc[:, idx+1], '-', linewidth=2, label=f'{species_name} (int)')
    ax.set_xlabel('Time', fontsize=10)
    ax.set_ylabel('Concentration', fontsize=10)
    ax.set_title('Conflict Species Dynamics', fontsize=11, fontweight='bold')
    ax.legend(fontsize=8, ncol=2)
    ax.grid(True, alpha=0.3)
    
    # Category aggregates
    ax = axes[1, 0]
    peace_baseline = ts_baseline.iloc[:, [i+1 for i in peace_species]].sum(axis=1)
    peace_intervention = ts_intervention.iloc[:, [i+1 for i in peace_species]].sum(axis=1)
    ax.plot(time, peace_baseline, '--', linewidth=2, label='Baseline', color='gray')
    ax.plot(time, peace_intervention, '-', linewidth=2, label='Intervention', color='green')
    ax.set_xlabel('Time', fontsize=10)
    ax.set_ylabel('Total Peace', fontsize=10)
    ax.set_title('Peace Category Aggregate', fontsize=11, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    
    ax = axes[1, 1]
    conflict_baseline = ts_baseline.iloc[:, [i+1 for i in conflict_species]].sum(axis=1)
    conflict_intervention = ts_intervention.iloc[:, [i+1 for i in conflict_species]].sum(axis=1)
    ax.plot(time, conflict_baseline, '--', linewidth=2, label='Baseline', color='gray')
    ax.plot(time, conflict_intervention, '-', linewidth=2, label='Intervention', color='red')
    ax.set_xlabel('Time', fontsize=10)
    ax.set_ylabel('Total Conflict', fontsize=10)
    ax.set_title('Conflict Category Aggregate', fontsize=11, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    
    plt.suptitle(f'Time Series Comparison\n{strategy_info["name"]}', 
                 fontsize=13, fontweight='bold', y=0.98)
    plt.tight_layout()
    save_path = os.path.join(output_dir, "time_series_comparison.png")
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {save_path}")
    plt.close()
    
    print("\n" + "="*80)
    print("ANALYSIS COMPLETE")
    print("="*80)
    print(f"\nAll outputs saved to: {output_dir}/")
    print("\nGenerated files:")
    print("  1. comprehensive_comparison.png - Full multiscale comparison")
    print("  2. category_peace_analysis.png - Peace category detailed analysis")
    print("  3. category_conflict_analysis.png - Conflict category detailed analysis")
    print("  4. time_series_comparison.png - Dynamic trajectories comparison")

if __name__ == "__main__":
    main()