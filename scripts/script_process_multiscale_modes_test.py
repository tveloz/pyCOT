#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Intervention Multiscale Semantic Analysis

Analyzes intervention strategies using multiscale process classification.
Two modes:
1. TOP3 MODE: Analyzes top 3 strategies from intervention_combinations.csv
2. CUSTOM MODE: Analyzes custom intervention parameters

For each strategy, produces:
- Time series dynamics plot (trajectories of all variables)
- Process distribution plots for each time window (1, 5, 10, 20)
"""

import sys
import os

# Set UTF-8 encoding for output
if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pyCOT.simulations.ode import simulation
from pyCOT.kinetics.deterministic_advanced import rate_saturated, rate_threshold_memory, rate_cosine
from pyCOT.io.functions import read_txt
from pyCOT.semantic_partition import define_semantic_categories
from pyCOT.process_analyzer import analyze_temporal_scale

# ========================================
# MODE SWITCH - CHANGE THIS TO SWITCH MODES
# ========================================
# MODE options:
#   'top3'             - Analyze top 3 strategies from CSV (ranks 1, 2, 3)
#   'list'             - Analyze specific strategies by index (set STRATEGY_LIST below)
#   'custom'           - Analyze custom intervention parameters (set CUSTOM_INTERVENTION below)
#   'net_value_success' - Analyze all 100% success strategies, plot cost vs net value,
#                        generate detailed outputs for top 10 by (peace-conflict)/cost
MODE = 'net_value_success'

# Strategy list (only used when MODE = 'list')
# List the indices of strategies you want to analyze from the CSV
# Examples:
#   [0, 1, 2] - Analyze top 3 strategies (same as 'top3' mode)
#   [0, 5, 10] - Analyze strategies at positions 0, 5, and 10
#   [15, 20, 25] - Compare strategies further down the ranking
# Note: Index 0 = Rank 1 (first strategy in CSV), Index 1 = Rank 2, etc.
STRATEGY_LIST = [0, 24, 110, 274]  # Change this to select different strategies

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
WINDOW_SIZES = [1, 5, 10, 20]

# Simulation parameters
SIMULATION_TIME = 200
N_STEPS = 400
INITIAL_CONDITION = [1.0, 0.0, 0.0, 0.5, 1.5, 0.0, 0.0]  # [SR, R, E, WR, DT, T, V]

# CSV path for top strategies
CSV_PATH = 'intervention_combinations.csv'

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
# CUSTOM INTERVENTION PARAMETERS
# (Only used when MODE = 'custom')
# ========================================
# Modify these values to test custom interventions
CUSTOM_INTERVENTION = {
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

def load_strategies_from_csv(csv_path, strategy_indices):
    """
    Load specific strategies from intervention combinations CSV by index

    Args:
        csv_path: Path to CSV file
        strategy_indices: List of indices or single integer for number of top strategies

    Returns:
        List of (params_dict, strategy_info) tuples
    """
    df = pd.read_csv(csv_path)

    # Get unique strategies (there are duplicates for different budgets)
    # Keep only the first occurrence of each unique strategy
    df_unique = df.drop_duplicates(subset=['Strategy'])

    # Handle both list of indices and integer for top-N
    if isinstance(strategy_indices, int):
        # If integer, take top N strategies
        indices_to_load = list(range(min(strategy_indices, len(df_unique))))
    else:
        # If list, use those specific indices
        indices_to_load = strategy_indices

    # Check if all indices are valid
    max_index = len(df_unique) - 1
    invalid_indices = [i for i in indices_to_load if i > max_index]
    if invalid_indices:
        print(f"Warning: Indices {invalid_indices} are out of range (max: {max_index}). Skipping them.")
        indices_to_load = [i for i in indices_to_load if i <= max_index]

    if not indices_to_load:
        raise ValueError("No valid strategy indices provided!")

    strategies = []
    for idx in indices_to_load:
        strategy = df_unique.iloc[idx]
        params = BASELINE_PARAMS.copy()

        # Parse combination format
        params_list = strategy['Parameters'].split(', ')
        magnitudes_list = [float(m.strip('%')) / 100 for m in strategy['Magnitudes'].split(', ')]
        strategy_str = strategy['Strategy']

        # Apply each parameter modification
        for param_name, magnitude in zip(params_list, magnitudes_list):
            if f"{param_name}↑" in strategy_str:
                params[param_name] = BASELINE_PARAMS[param_name] * (1 + magnitude)
            else:  # ↓
                params[param_name] = BASELINE_PARAMS[param_name] * (1 - magnitude)

        strategy_info = {
            'name': strategy_str,
            'success_rate': strategy['Success Rate (%)'],
            'cost': strategy['Total Cost'],
            'cost_effectiveness': strategy['Cost-Effectiveness'],
            'csv_index': idx,  # The actual index in the CSV (0-based)
            'rank': idx + 1    # Human-readable rank (1-based)
        }

        strategies.append((params, strategy_info))

    return strategies

def load_successful_strategies(csv_path, success_threshold=100.0):
    """
    Load all strategies with 100% success rate from CSV

    Args:
        csv_path: Path to CSV file
        success_threshold: Minimum success rate (default 100%)

    Returns:
        List of (params_dict, strategy_info) tuples for successful strategies
    """
    df = pd.read_csv(csv_path)

    # Filter for successful strategies
    df_success = df[df['Success Rate (%)'] >= success_threshold]

    # Get unique strategies
    df_unique = df_success.drop_duplicates(subset=['Strategy'])

    print(f"   Found {len(df_unique)} unique strategies with {success_threshold}% success rate")

    strategies = []
    for idx in range(len(df_unique)):
        strategy = df_unique.iloc[idx]
        params = BASELINE_PARAMS.copy()

        # Parse combination format
        params_list = strategy['Parameters'].split(', ')
        magnitudes_list = [float(m.strip('%')) / 100 for m in strategy['Magnitudes'].split(', ')]
        strategy_str = strategy['Strategy']

        # Apply each parameter modification
        for param_name, magnitude in zip(params_list, magnitudes_list):
            if f"{param_name}↑" in strategy_str:
                params[param_name] = BASELINE_PARAMS[param_name] * (1 + magnitude)
            else:  # ↓
                params[param_name] = BASELINE_PARAMS[param_name] * (1 - magnitude)

        strategy_info = {
            'name': strategy_str,
            'success_rate': strategy['Success Rate (%)'],
            'cost': strategy['Total Cost'],
            'cost_effectiveness': strategy['Cost-Effectiveness'],
            'csv_index': idx,
            'rank': idx + 1
        }

        strategies.append((params, strategy_info))

    return strategies

def apply_custom_intervention(custom_dict):
    """Apply custom intervention to baseline parameters"""
    params = BASELINE_PARAMS.copy()
    modified_params = []

    for key, value in custom_dict.items():
        if value is not None:
            params[key] = value
            pct_change = ((value / BASELINE_PARAMS[key]) - 1) * 100
            direction = "↑" if pct_change > 0 else "↓"
            modified_params.append(f"{key}{direction}{abs(pct_change):.1f}%")

    strategy_info = {
        'name': 'Custom: ' + ', '.join(modified_params),
        'success_rate': None,
        'cost': None,
        'cost_effectiveness': None,
        'rank': 'Custom'
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

# ========================================
# VISUALIZATION FUNCTIONS
# ========================================

def plot_dynamics(ts, species_list, semantic_partition, strategy_info, save_path):
    """
    Plot time series dynamics of all variables
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Get species indices for each category
    peace_species = semantic_partition.category_indices['peace']
    conflict_species = semantic_partition.category_indices['conflict']

    time = ts['Time'].values

    # Peace species
    ax = axes[0]
    for idx in peace_species:
        species_name = species_list[idx]
        ax.plot(time, ts.iloc[:, idx+1], linewidth=2, label=species_name)
    ax.set_xlabel('Time', fontsize=11)
    ax.set_ylabel('Concentration', fontsize=11)
    ax.set_title('Peace Variables', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # Conflict species
    ax = axes[1]
    for idx in conflict_species:
        species_name = species_list[idx]
        ax.plot(time, ts.iloc[:, idx+1], linewidth=2, label=species_name)
    ax.set_xlabel('Time', fontsize=11)
    ax.set_ylabel('Concentration', fontsize=11)
    ax.set_title('Conflict Variables', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # Overall title
    rank_str = f"Rank #{strategy_info['rank']}: " if strategy_info['rank'] != 'Custom' else ""
    plt.suptitle(f'{rank_str}{strategy_info["name"]}\nDynamics Trajectories',
                 fontsize=13, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {save_path}")
    plt.close()

def plot_process_distributions(scale_results, strategy_info, save_path):
    """
    Plot process distributions for each category across all time windows

    Creates one plot per time window showing distribution of balance states
    (overproduced, depleted, balanced) for each category (peace, conflict)
    """
    window_sizes = scale_results['scales']
    categories = list(scale_results['results_by_scale'][window_sizes[0]]['category_statistics'].keys())

    # Create figure with subplots: one row per window size
    fig, axes = plt.subplots(len(window_sizes), 1, figsize=(12, 4*len(window_sizes)))
    if len(window_sizes) == 1:
        axes = [axes]

    for i, ws in enumerate(window_sizes):
        ax = axes[i]

        # Get process type distributions for this window size
        result = scale_results['results_by_scale'][ws]

        # Balance states to plot
        balance_states = ['overproduced', 'depleted', 'balanced']

        # Collect data for each category
        data = {}
        for category in categories:
            cat_stats = result['category_statistics'][category]
            balance_dist = cat_stats.get('balance_state_distribution', {})
            data[category] = [balance_dist.get(state, 0) for state in balance_states]

        # Bar plot
        x = np.arange(len(balance_states))
        width = 0.35

        colors = {'peace': 'green', 'conflict': 'red'}

        for j, category in enumerate(categories):
            offset = (j - 0.5) * width
            ax.bar(x + offset, data[category], width,
                   label=category.upper(),
                   color=colors.get(category, 'blue'),
                   alpha=0.7)

        ax.set_xlabel('Balance State', fontsize=11)
        ax.set_ylabel('Count', fontsize=11)
        ax.set_title(f'Window Size = {ws} timesteps', fontsize=12, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(balance_states, rotation=45, ha='right')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3, axis='y')

    # Overall title
    rank_str = f"Rank #{strategy_info['rank']}: " if strategy_info['rank'] != 'Custom' else ""
    plt.suptitle(f'{rank_str}{strategy_info["name"]}\nBalance State Distributions Across Time Windows',
                 fontsize=13, fontweight='bold', y=0.995)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {save_path}")
    plt.close()

def plot_phase_space(ts, species_list, strategy_info, save_path):
    """
    Plot phase space showing resilience dynamics
    - SR vs WR (Strong Regions vs Weak Regions)
    - V vs T (Violence vs Trust)

    Args:
        ts: Time series dataframe
        species_list: List of species names
        strategy_info: Strategy metadata
        save_path: Where to save the plot
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Get indices for relevant species
    sr_idx = species_list.index('SR')
    wr_idx = species_list.index('WR')
    v_idx = species_list.index('V')
    t_idx = species_list.index('T')

    # Get time values for coloring
    time = ts['Time'].values

    # Plot 1: SR vs WR (Regional Resilience)
    ax = axes[0]
    sr_vals = ts.iloc[:, sr_idx + 1].values
    wr_vals = ts.iloc[:, wr_idx + 1].values

    scatter = ax.scatter(sr_vals, wr_vals, c=time, cmap='viridis',
                        s=30, alpha=0.6, edgecolors='black', linewidth=0.3)

    # Add trajectory direction with arrows
    step = max(1, len(time) // 10)  # Show ~10 arrows
    for i in range(0, len(time) - step, step):
        ax.annotate('', xy=(sr_vals[i+step], wr_vals[i+step]),
                   xytext=(sr_vals[i], wr_vals[i]),
                   arrowprops=dict(arrowstyle='->', color='gray', alpha=0.5, lw=1))

    # Mark start and end
    ax.plot(sr_vals[0], wr_vals[0], 'go', markersize=10, label='Start', zorder=5)
    ax.plot(sr_vals[-1], wr_vals[-1], 'r*', markersize=15, label='End', zorder=5)

    ax.set_xlabel('SR (Strong Regions)', fontsize=11)
    ax.set_ylabel('WR (Weak Regions)', fontsize=11)
    ax.set_title('Regional Resilience Phase Space', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    plt.colorbar(scatter, ax=ax, label='Time')

    # Plot 2: V vs T (Violence-Trust Dynamics)
    ax = axes[1]
    v_vals = ts.iloc[:, v_idx + 1].values
    t_vals = ts.iloc[:, t_idx + 1].values

    scatter = ax.scatter(t_vals, v_vals, c=time, cmap='viridis',
                        s=30, alpha=0.6, edgecolors='black', linewidth=0.3)

    # Add trajectory direction with arrows
    for i in range(0, len(time) - step, step):
        ax.annotate('', xy=(t_vals[i+step], v_vals[i+step]),
                   xytext=(t_vals[i], v_vals[i]),
                   arrowprops=dict(arrowstyle='->', color='gray', alpha=0.5, lw=1))

    # Mark start and end
    ax.plot(t_vals[0], v_vals[0], 'go', markersize=10, label='Start', zorder=5)
    ax.plot(t_vals[-1], v_vals[-1], 'r*', markersize=15, label='End', zorder=5)

    ax.set_xlabel('T (Trust)', fontsize=11)
    ax.set_ylabel('V (Violence)', fontsize=11)
    ax.set_title('Violence-Trust Phase Space', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    plt.colorbar(scatter, ax=ax, label='Time')

    # Overall title
    rank_str = f"Rank #{strategy_info['rank']}: " if strategy_info['rank'] != 'Custom' else ""
    plt.suptitle(f'{rank_str}{strategy_info["name"]}\nResilience Phase Space Dynamics',
                 fontsize=13, fontweight='bold', y=1.00)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {save_path}")
    plt.close()

def plot_cost_vs_net_value(results_data, save_path):
    """
    Plot cost vs net value for all successful strategies
    Uses different colors for each window size to show overlaps clearly

    Args:
        results_data: List of dicts with keys: 'cost', 'net_values_by_window', 'strategy_name', 'metric'
        save_path: Where to save the plot
    """
    window_sizes = WINDOW_SIZES

    fig, ax = plt.subplots(1, 1, figsize=(12, 8))

    # Color map for different window sizes
    colors = plt.cm.viridis(np.linspace(0, 1, len(window_sizes)))

    # Plot each window size with a different color
    for i, ws in enumerate(window_sizes):
        # Extract data for this window
        costs = [r['cost'] for r in results_data]
        net_values = [r['net_values_by_window'][ws] for r in results_data]

        # Scatter plot
        ax.scatter(costs, net_values, alpha=0.7, s=60, c=[colors[i]],
                  label=f'Window = {ws}', edgecolors='black', linewidth=0.5)

    ax.set_xlabel('Cost ($)', fontsize=13)
    ax.set_ylabel('Net Value (Peace - Conflict)', fontsize=13)
    ax.set_title('Cost vs Net Value for All Successful Strategies\nColors Represent Different Time Windows',
                 fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color='red', linestyle='--', linewidth=1.5, alpha=0.6, label='Neutral (Net=0)')
    ax.legend(fontsize=11, loc='best', framealpha=0.9)

    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {save_path}")
    plt.close()

# ========================================
# MAIN EXECUTION
# ========================================

def main():
    print("="*80)
    print("INTERVENTION MULTISCALE SEMANTIC ANALYSIS")
    print("="*80)
    print(f"\nMODE: {MODE.upper()}")
    print("-"*80)

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
    # LOAD STRATEGIES
    # ========================================
    print("\n2. Loading strategies...")

    if MODE == 'top3':
        strategies = load_strategies_from_csv(CSV_PATH, 3)
        print(f"\n   Loaded top 3 strategies from {CSV_PATH}:")
        for params, info in strategies:
            print(f"   - Rank #{info['rank']}: {info['name']}")
            print(f"     Success Rate: {info['success_rate']:.1f}%, Cost: ${info['cost']:.0f}")

    elif MODE == 'list':
        strategies = load_strategies_from_csv(CSV_PATH, STRATEGY_LIST)
        print(f"\n   Loaded {len(strategies)} strategies from {CSV_PATH}:")
        for params, info in strategies:
            print(f"   - Index {info['csv_index']} (Rank #{info['rank']}): {info['name']}")
            print(f"     Success Rate: {info['success_rate']:.1f}%, Cost: ${info['cost']:.0f}")

    elif MODE == 'net_value_success':
        all_strategies = load_successful_strategies(CSV_PATH, success_threshold=100.0)
        strategies = all_strategies  # Will process all, then rank by metric
        print(f"\n   Will analyze all {len(strategies)} successful strategies")
        print("   Then rank by metric: (Peace Net Effect - Conflict Net Effect) / Cost")

    else:  # custom mode
        params, info = apply_custom_intervention(CUSTOM_INTERVENTION)
        strategies = [(params, info)]
        print(f"\n   Using custom intervention:")
        print(f"   - {info['name']}")

    # ========================================
    # PROCESS EACH STRATEGY
    # ========================================

    output_dir = "./outputs/intervention_multiscale"
    os.makedirs(output_dir, exist_ok=True)

    # Special handling for net_value_success mode
    if MODE == 'net_value_success':
        results_data = []

        for idx, (intervention_params, strategy_info) in enumerate(strategies):
            print(f"\n[{idx+1}/{len(strategies)}] Processing: {strategy_info['name'][:60]}...")

            # Run simulation (quietly)
            ts, _, process_series = run_simulation_with_params(
                intervention_params, rn, rate_list, additional_laws
            )

            # Multiscale analysis
            scale_results = analyze_temporal_scale(
                process_series=process_series,
                S=S,
                semantic_partition=semantic_partition,
                window_sizes=WINDOW_SIZES,
                stride=1
            )

            # Extract net effects for peace and conflict at each window
            net_values_by_window = {}
            for ws in WINDOW_SIZES:
                peace_net = scale_results['results_by_scale'][ws]['category_statistics']['peace']['net_effect']['mean']
                conflict_net = scale_results['results_by_scale'][ws]['category_statistics']['conflict']['net_effect']['mean']
                # Since conflict is negative, subtracting it adds the absolute value
                net_values_by_window[ws] = peace_net - conflict_net

            # Compute metric: average net value across windows divided by cost
            avg_net_value = np.mean(list(net_values_by_window.values()))
            metric = avg_net_value / strategy_info['cost'] if strategy_info['cost'] > 0 else 0

            results_data.append({
                'params': intervention_params,
                'strategy_info': strategy_info,
                'scale_results': scale_results,
                'ts': ts,
                'net_values_by_window': net_values_by_window,
                'avg_net_value': avg_net_value,
                'cost': strategy_info['cost'],
                'metric': metric
            })

        # Sort by metric (descending)
        results_data.sort(key=lambda x: x['metric'], reverse=True)

        # Print top 10
        print("\n" + "="*80)
        print("TOP 10 STRATEGIES BY (Peace - Conflict) / Cost")
        print("="*80)
        for i in range(min(10, len(results_data))):
            r = results_data[i]
            print(f"\n{i+1}. {r['strategy_info']['name']}")
            print(f"   Cost: ${r['cost']:.0f}")
            print(f"   Avg Net Value: {r['avg_net_value']:+.4f}")
            print(f"   Metric: {r['metric']:.6f}")
            print(f"   Net Values by Window: {', '.join([f'{ws}={v:+.4f}' for ws, v in r['net_values_by_window'].items()])}")

        # Generate cost vs net value plot
        print("\n3. Generating cost vs net value plot...")
        plot_cost_vs_net_value(results_data, os.path.join(output_dir, "cost_vs_net_value.png"))

        # Generate detailed plots for ALL successful strategies
        print(f"\n4. Generating detailed plots for all {len(results_data)} successful strategies...")
        for i in range(len(results_data)):
            r = results_data[i]
            print(f"   [{i+1}/{len(results_data)}] Generating plots for: {r['strategy_info']['name'][:60]}...")

            # Use strategy index for filename (sorted by metric)
            base_name = f"strategy_{i+1:03d}_metric"

            # Dynamics
            dynamics_path = os.path.join(output_dir, f"{base_name}_dynamics.png")
            plot_dynamics(r['ts'], species_list, semantic_partition, r['strategy_info'], dynamics_path)

            # Process distributions
            processes_path = os.path.join(output_dir, f"{base_name}_process_distributions.png")
            plot_process_distributions(r['scale_results'], r['strategy_info'], processes_path)

            # Phase space
            phase_space_path = os.path.join(output_dir, f"{base_name}_phase_space.png")
            plot_phase_space(r['ts'], species_list, r['strategy_info'], phase_space_path)

    else:
        # Original mode logic
        for intervention_params, strategy_info in strategies:
            print(f"\n{'='*80}")
            print(f"PROCESSING: {strategy_info['name']}")
            print(f"{'='*80}")

            # Run simulation
            print("\n3. Running simulation...")
            ts, _, process_series = run_simulation_with_params(
                intervention_params, rn, rate_list, additional_laws
            )

            # Multiscale analysis
            print(f"\n4. Performing multiscale analysis (windows: {WINDOW_SIZES})...")
            scale_results = analyze_temporal_scale(
                process_series=process_series,
                S=S,
                semantic_partition=semantic_partition,
                window_sizes=WINDOW_SIZES,
                stride=1
            )

            # Print summary
            print("\n5. Process Analysis Summary:")
            categories = list(scale_results['results_by_scale'][WINDOW_SIZES[0]]['category_statistics'].keys())

            for category in categories:
                print(f"\n   {category.upper()}:")
                print(f"   {'Window':<10} {'Overprod':<12} {'Depleted':<12} {'Balanced':<12} {'Net Effect':<12}")
                print(f"   {'-'*60}")

                for ws in WINDOW_SIZES:
                    stats = scale_results['results_by_scale'][ws]['category_statistics'][category]
                    balance_dist = stats.get('balance_state_distribution', {})
                    overprod = balance_dist.get('overproduced', 0)
                    depleted = balance_dist.get('depleted', 0)
                    balanced = balance_dist.get('balanced', 0)
                    net_effect = stats['net_effect']['mean']
                    print(f"   {ws:<10} {overprod:<12} {depleted:<12} {balanced:<12} {net_effect:<+12.4f}")

            # Generate plots
            print("\n6. Generating plots...")

            # File naming
            if MODE == 'custom':
                base_name = "custom"
            else:  # 'top3' or 'list' mode
                base_name = f"rank{strategy_info['rank']}"

            # Plot 1: Dynamics
            dynamics_path = os.path.join(output_dir, f"{base_name}_dynamics.png")
            plot_dynamics(ts, species_list, semantic_partition, strategy_info, dynamics_path)

            # Plot 2: Process distributions
            processes_path = os.path.join(output_dir, f"{base_name}_process_distributions.png")
            plot_process_distributions(scale_results, strategy_info, processes_path)

    print("\n" + "="*80)
    print("ANALYSIS COMPLETE")
    print("="*80)
    print(f"\nAll outputs saved to: {output_dir}/")
    print("\nGenerated files per strategy:")
    print("  1. [rank]_dynamics.png - Time series trajectories of all variables")
    print("  2. [rank]_process_distributions.png - Balance states (overproduced/depleted/balanced)")
    print("                                         for each time window (1, 5, 10, 20)")
    print("\nHow to interpret balance states:")
    print("  - Overproduced: Category has net positive production")
    print("  - Depleted: Category has net negative production (consumption > production)")
    print("  - Balanced: Category production and consumption are nearly equal")

if __name__ == "__main__":
    main()
