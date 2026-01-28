
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pyCOT.simulations.ode import simulation
from pyCOT.kinetics.deterministic_advanced import rate_saturated, rate_threshold_memory, rate_cosine
from pyCOT.io.functions import read_txt

# ========================================
# LOAD MODEL
# ========================================
file_path = 'networks/Conflict_Theory/Resource_Scarcity_Toy_Model2.txt'
rn = read_txt(file_path)

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
# BASELINE PARAMETERS
# ========================================
baseline_params = {
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
# COST PER % CHANGE (INTERVENTION ECONOMICS)
# ========================================
# Cost in arbitrary units per 1% change
# Higher cost = harder/more expensive to change

intervention_costs = {
    # EXPENSIVE: Infrastructure & Production (100-200 per %)
    'Vmax_production': 150,      # Building productive capacity - very expensive
    'R_amplitude': 120,           # Environmental restoration - very expensive
    
    # HIGH: Institutional capacity (80-120 per %)
    'Vmax_trust': 100,            # Building trust institutions - expensive
    'k_recovery': 90,             # Recovery programs - expensive
    
    # MODERATE-HIGH: Security & Prevention (50-80 per %)
    'k_violence_weak': 70,        # Community conflict prevention
    'k_violence_detached_weak': 75,   # Targeted violence prevention
    'k_violence_detached': 80,    # High-risk group interventions
    'k_detachment_base': 65,      # Preventing disengagement
    
    # MODERATE: Social programs (30-60 per %)
    'k_trust_destruction': 50,    # Protecting trust from violence
    'k_trust_decay': 45,          # Maintaining social cohesion
    'k_degradation': 55,          # Preventing degradation
    'k_violence_decay': 40,       # Violence resolution mechanisms
    
    # MODERATE-LOW: Policy & Thresholds (20-40 per %)
    'E_threshold': 35,            # Reducing barriers to economic recovery
    'T_threshold_weak_strong': 30,     # Easier trust-based strengthening
    'T_threshold_detached_weak': 38,   # Reducing detachment barriers
    'detachment_threshold': 32,   # Raising detachment barriers
    
    # LOW: Efficiency & Decay (10-25 per %)
    'k_resource_depletion': 25,   # Resource management efficiency
    'k_economic_decay': 22,       # Economic sustainability
    'Km_production': 20,          # Production efficiency
    'Km_trust_economy': 18,       # Trust-building efficiency
}

# ========================================
# INTERVENTION LOGIC (DIRECTION OF CHANGE)
# ========================================
# +1 means INCREASE parameter for peace
# -1 means DECREASE parameter for peace

intervention_direction = {
    # INCREASE for peace
    'Vmax_production': +1,        # More production capacity
    'k_recovery': +1,             # Faster recovery
    'Vmax_trust': +1,             # More trust generation
    'R_amplitude': +1,            # More resources
    'k_violence_decay': +1,       # Faster violence dissipation
    'detachment_threshold': +1,   # Harder to become detached
    
    # DECREASE for peace
    'k_degradation': -1,          # Slower degradation
    'k_detachment_base': -1,      # Slower detachment
    'k_trust_destruction': -1,    # Less trust destruction
    'k_trust_decay': -1,          # Slower trust decay
    'k_violence_weak': -1,        # Less violence from weak populations
    'k_violence_detached_weak': -1,   # Less detached-weak violence
    'k_violence_detached': -1,    # Less detached violence
    'E_threshold': -1,            # Easier economic recovery
    'T_threshold_weak_strong': -1,     # Easier trust strengthening
    'T_threshold_detached_weak': -1,   # Easier detachment recovery
    'k_resource_depletion': -1,   # Slower resource loss
    'k_economic_decay': -1,       # Slower economic decay
    'Km_production': -1,          # Lower Km means easier to reach max production
    'Km_trust_economy': -1,       # Lower Km means easier to generate trust
}

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

def check_peace(ts_data, fv_data, last_n_steps=50):
    """Check if last N steps show peaceful state

    ts_data: pandas DataFrame with Time column + species columns
    fv_data: pandas DataFrame with Time column + species columns
    """
    # Peace criteria: low violence, high SR, low DT
    # x0 = [SR, R, E, WR, DT, T, V]
    # ts_data has columns: ['Time', 'SR', 'R', 'E', 'WR', 'DT', 'T', 'V']

    # Get last N steps
    last_steps = ts_data.tail(last_n_steps)

    # Calculate averages using column names
    avg_V = last_steps['V'].mean()
    avg_SR = last_steps['SR'].mean()
    avg_DT = last_steps['DT'].mean()

    # Peace = V < 0.5 AND SR > 0.5 AND DT < 1.0
    is_peace = (avg_V < 0.5) and (avg_SR > 1.0) and (avg_DT < 0.5)
    return is_peace

# ========================================
# MONTE CARLO SIMULATION
# ========================================
def run_intervention_montecarlo(n_trials=500, seed=42, verbose=False):
    """
    Run Monte Carlo simulation of interventions

    Returns:
    - costs: array of intervention costs
    - successes: array of success (1) or failure (0)
    """
    np.random.seed(seed)

    # Initial condition
    x0 = [1.0, 0.0, 0.0, 0.5, 1.5, 0.0, 0.0]  # [SR, R, E, WR, DT, T, V]

    costs = []
    successes = []

    # Get intervene-able parameters (those with costs defined)
    param_names = list(intervention_costs.keys())

    # Debugging counters
    skipped_budget = 0
    failed_sim = 0

    for trial in range(n_trials):
        # Random budget between 500 and 5000 units (increased upper limit)
        budget = 10000

        # Random number of parameters to intervene on (1-3)
        n_params = np.random.randint(1, 4)

        # Random parameter selection
        selected_params = np.random.choice(param_names, size=n_params, replace=False)

        # Random intervention magnitude (5%, 10%, 15%, or 20%)
        magnitude = np.random.choice([5,10,15,20,25,30,35,40,45,50])  # percent change
        
        # Calculate cost of this intervention
        total_cost = sum(intervention_costs[p] * magnitude for p in selected_params)
        
        # Only proceed if within budget
        if total_cost > budget:
            skipped_budget += 1
            if verbose and trial < 10:
                print(f"Trial {trial}: Cost ${total_cost:.0f} exceeds budget ${budget:.0f}")
            continue

        # Apply interventions
        modified_params = baseline_params.copy()
        for param in selected_params:
            direction = intervention_direction[param]
            pct_change = magnitude / 100.0
            if direction > 0:
                # Increase parameter
                modified_params[param] = baseline_params[param] * (1 + pct_change)
            else:
                # Decrease parameter
                modified_params[param] = baseline_params[param] * (1 - pct_change)

        # Build spec_vector
        spec_vector = build_spec_vector(modified_params)

        # Run simulation
        try:
            ts, fv = simulation(
                rn,
                rate=rate_list,
                spec_vector=spec_vector,
                x0=x0,
                t_span=(0, 200),
                n_steps=400,
                additional_laws=additional_laws,
                verbose=False
            )
            # Debug: print shape on first successful run
            if len(costs) == 0 and verbose:
                print(f"\nDEBUG: fv.shape = {fv.shape}, ts.shape = {ts.shape}")
                print(f"Expected fv shape: (n_steps={400}, n_species={len(x0)})")

            # Check peace
            peace = check_peace(ts, fv, last_n_steps=50)

            costs.append(total_cost)
            successes.append(1 if peace else 0)

        except Exception as e:
            failed_sim += 1
            if verbose:
                print(f"Trial {trial} failed: {e}")
            continue

    # Print summary
    print(f"\n  Budget rejections: {skipped_budget}/{n_trials}")
    print(f"  Simulation failures: {failed_sim}/{n_trials}")
    print(f"  Successful runs: {len(costs)}/{n_trials}")

    return np.array(costs), np.array(successes)

# ========================================
# SYSTEMATIC INTERVENTION EVALUATION
# ========================================
def evaluate_single_intervention(param_name, magnitude_pct):
    """
    Evaluate a single intervention strategy (deterministic, single run)

    Returns: success (0 or 1), cost
    """
    x0 = [1.0, 0.0, 0.0, 0.5, 1.5, 0.0, 0.0]

    cost = intervention_costs[param_name] * magnitude_pct
    direction = intervention_direction[param_name]
    pct_change = magnitude_pct / 100.0

    modified_params = baseline_params.copy()

    if direction > 0:
        modified_params[param_name] = baseline_params[param_name] * (1 + pct_change)
    else:
        modified_params[param_name] = baseline_params[param_name] * (1 - pct_change)

    spec_vector = build_spec_vector(modified_params)

    try:
        ts, fv = simulation(
            rn,
            rate=rate_list,
            spec_vector=spec_vector,
            x0=x0,
            t_span=(0, 200),
            n_steps=400,
            additional_laws=additional_laws,
            verbose=False
        )

        peace = check_peace(ts, fv, last_n_steps=50)
        return (1 if peace else 0), cost

    except Exception:
        return 0, cost

def evaluate_combination_intervention(param_names, magnitudes_pct):
    """
    Evaluate a combination intervention strategy (deterministic, single run)

    param_names: list of parameter names
    magnitudes_pct: list of magnitudes (same length as param_names)

    Returns: success (0 or 1), total_cost, cost_effectiveness
    """
    x0 = [1.0, 0.0, 0.0, 0.5, 1.5, 0.0, 0.0]

    # Calculate total cost
    total_cost = sum(intervention_costs[p] * m for p, m in zip(param_names, magnitudes_pct))

    modified_params = baseline_params.copy()

    for param_name, magnitude_pct in zip(param_names, magnitudes_pct):
        direction = intervention_direction[param_name]
        pct_change = magnitude_pct / 100.0

        if direction > 0:
            modified_params[param_name] = baseline_params[param_name] * (1 + pct_change)
        else:
            modified_params[param_name] = baseline_params[param_name] * (1 - pct_change)

    spec_vector = build_spec_vector(modified_params)

    try:
        ts, fv = simulation(
            rn,
            rate=rate_list,
            spec_vector=spec_vector,
            x0=x0,
            t_span=(0, 200),
            n_steps=400,
            additional_laws=additional_laws,
            verbose=False
        )

        peace = check_peace(ts, fv, last_n_steps=50)
        success = 1 if peace else 0
        cost_effectiveness = success / total_cost if total_cost > 0 else 0

        return success, total_cost, cost_effectiveness

    except Exception:
        return 0, total_cost, 0.0

def rank_interventions(magnitudes=[10, 20, 30]):
    """
    Systematically evaluate and rank all single-parameter interventions

    Returns: DataFrame with rankings sorted by cost-effectiveness
    """
    print("\n" + "=" * 80)
    print("SYSTEMATIC INTERVENTION RANKING")
    print("=" * 80)

    results = []
    param_names = list(intervention_costs.keys())

    total_evals = len(param_names) * len(magnitudes)
    eval_count = 0

    for param_name in param_names:
        for magnitude in magnitudes:
            eval_count += 1
            print(f"Evaluating {eval_count}/{total_evals}: {param_name} @ {magnitude}%...", end='\r')

            success, cost = evaluate_single_intervention(param_name, magnitude)

            cost_effectiveness = success / cost if cost > 0 else 0

            direction = "↑" if intervention_direction[param_name] > 0 else "↓"

            results.append({
                'Parameter': param_name,
                'Direction': direction,
                'Magnitude (%)': magnitude,
                'Success Rate (%)': success * 100,
                'Cost': cost,
                'Cost-Effectiveness': cost_effectiveness * 1000,  # per 1000 units
                'Intervention': f"{param_name} {direction}{magnitude}%"
            })

    print(" " * 100, end='\r')  # Clear progress line

    df = pd.DataFrame(results)
    df_sorted = df.sort_values('Cost-Effectiveness', ascending=False)

    return df_sorted

def find_best_combinations(budget=100):
    """
    Find best combination strategies within a budget
    """
    print("\n" + "=" * 80)
    print(f"FINDING BEST COMBINATION STRATEGIES (Budget: ${budget})")
    print("=" * 80)

    param_names = list(intervention_costs.keys())
    results = []

    # Test 2-parameter combinations
    print("\nTesting 2-parameter combinations...")
    combo_count = 0
    for i, p1 in enumerate(param_names):
        for p2 in param_names[i+1:]:
            combo_count += 1
            # Test different magnitude combinations
            for m1 in [10, 20, 30]:
                for m2 in [10, 20, 30]:
                    cost = intervention_costs[p1] * m1 + intervention_costs[p2] * m2
                    if cost <= budget:
                        success, total_cost, ce = evaluate_combination_intervention(
                            [p1, p2], [m1, m2]
                        )

                        d1 = "↑" if intervention_direction[p1] > 0 else "↓"
                        d2 = "↑" if intervention_direction[p2] > 0 else "↓"

                        results.append({
                            'Strategy': f"{p1}{d1}{m1}% + {p2}{d2}{m2}%",
                            'N_Params': 2,
                            'Success Rate (%)': success * 100,
                            'Total Cost': total_cost,
                            'Cost-Effectiveness': ce * 1000,
                            'Parameters': f"{p1}, {p2}",
                            'Magnitudes': f"{m1}%, {m2}%"
                        })
            print(f"  Tested {combo_count} combinations...", end='\r')

    print(" " * 100, end='\r')

    # Test 3-parameter combinations (limited sampling)
    print("\nTesting 3-parameter combinations (sampling)...")
    np.random.seed(42)
    for _ in range(100):  # Sample 100 random 3-param combinations
        params = np.random.choice(param_names, size=3, replace=False)
        mags = np.random.choice([10, 20, 30], size=3)

        cost = sum(intervention_costs[p] * m for p, m in zip(params, mags))
        if cost <= budget:
            success, total_cost, ce = evaluate_combination_intervention(
                params.tolist(), mags.tolist()
            )

            strategy_str = " + ".join([
                f"{p}{'↑' if intervention_direction[p] > 0 else '↓'}{m}%"
                for p, m in zip(params, mags)
            ])

            results.append({
                'Strategy': strategy_str,
                'N_Params': 3,
                'Success Rate (%)': success * 100,
                'Total Cost': total_cost,
                'Cost-Effectiveness': ce * 1000,
                'Parameters': ", ".join(params),
                'Magnitudes': ", ".join([f"{m}%" for m in mags])
            })

    df = pd.DataFrame(results)
    df_sorted = df.sort_values('Cost-Effectiveness', ascending=False)

    return df_sorted

# ========================================
# RUN SIMULATION & PLOT
# ========================================
print("=" * 80)
print("INTERVENTION COST-BENEFIT ANALYSIS")
print("=" * 80)
print("\nINTERVENTION COSTS (per 1% change):")
print("-" * 80)
for param, cost in sorted(intervention_costs.items(), key=lambda x: -x[1]):
    direction = "INCREASE" if intervention_direction[param] > 0 else "DECREASE"
    print(f"{param:30s} | {direction:8s} | ${cost:6.1f} per 1%")

# ========================================
# PART 1: RANK SINGLE INTERVENTIONS
# ========================================
rankings_df = rank_interventions(magnitudes=[10, 20, 30, 40, 50])

print("\n" + "=" * 80)
print("TOP 20 MOST COST-EFFECTIVE SINGLE INTERVENTIONS")
print("=" * 80)
print(rankings_df.head(20).to_string(index=False))

print("\n" + "=" * 80)
print("BOTTOM 10 LEAST COST-EFFECTIVE SINGLE INTERVENTIONS")
print("=" * 80)
print(rankings_df.tail(10).to_string(index=False))

# Save full rankings
rankings_df.to_csv('intervention_rankings.csv', index=False)
print("\nFull rankings saved to: intervention_rankings.csv")

# ========================================
# PART 2: FIND BEST COMBINATION STRATEGIES
# ========================================
# Test multiple budget levels
budgets = [2000, 4000, 6000]
all_combos = []

for budget in budgets:
    combos_df = find_best_combinations(budget=budget)
    combos_df['Budget'] = budget
    all_combos.append(combos_df)

combined_df = pd.concat(all_combos, ignore_index=True)
combined_df_sorted = combined_df.sort_values('Cost-Effectiveness', ascending=False)

print("\n" + "=" * 80)
print("TOP 20 COMBINATION STRATEGIES (ALL BUDGETS)")
print("=" * 80)
print(combined_df_sorted.head(20).to_string(index=False))

# Save combination results
combined_df_sorted.to_csv('intervention_combinations.csv', index=False)
print("\nFull combination rankings saved to: intervention_combinations.csv")

# ========================================
# PART 3: POLICY MAKER RECOMMENDATIONS
# ========================================
print("\n" + "=" * 80)
print("POLICY MAKER RECOMMENDATIONS")
print("=" * 80)

# Best overall strategy
best = combined_df_sorted.iloc[0]
print(f"\n1. MOST COST-EFFECTIVE OVERALL:")
print(f"   Strategy: {best['Strategy']}")
print(f"   Success Rate: {best['Success Rate (%)']:.1f}%")
print(f"   Total Cost: ${best['Total Cost']:.0f}")
print(f"   Cost-Effectiveness: {best['Cost-Effectiveness']:.3f} per $1000")

# Best by budget category
print("\n2. BEST STRATEGY BY BUDGET:")
for budget in budgets:
    budget_strategies = combined_df_sorted[combined_df_sorted['Budget'] == budget]
    if len(budget_strategies) > 0:
        best_in_budget = budget_strategies.iloc[0]
        print(f"\n   Budget ${budget}:")
        print(f"   - {best_in_budget['Strategy']}")
        print(f"   - Success: {best_in_budget['Success Rate (%)']:.1f}% | Cost: ${best_in_budget['Total Cost']:.0f}")

# Best single interventions by category
print("\n3. BEST SINGLE INTERVENTIONS BY CATEGORY:")

# Infrastructure interventions
infra = rankings_df[rankings_df['Parameter'].isin(['Vmax_production', 'R_amplitude'])].head(3)
print("\n   Infrastructure & Production:")
for _, row in infra.iterrows():
    print(f"   - {row['Intervention']}: {row['Success Rate (%)']:.1f}% success, ${row['Cost']:.0f}")

# Social programs
social = rankings_df[rankings_df['Parameter'].isin(['k_trust_destruction', 'k_trust_decay',
                                                     'Vmax_trust', 'k_violence_decay'])].head(3)
print("\n   Social Programs:")
for _, row in social.iterrows():
    print(f"   - {row['Intervention']}: {row['Success Rate (%)']:.1f}% success, ${row['Cost']:.0f}")

# Policy changes
policy = rankings_df[rankings_df['Parameter'].isin(['E_threshold', 'T_threshold_weak_strong',
                                                     'T_threshold_detached_weak', 'detachment_threshold'])].head(3)
print("\n   Policy Changes:")
for _, row in policy.iterrows():
    print(f"   - {row['Intervention']}: {row['Success Rate (%)']:.1f}% success, ${row['Cost']:.0f}")

# ========================================
# PART 4: VISUALIZATIONS
# ========================================
print("\n" + "=" * 80)
print("GENERATING VISUALIZATIONS")
print("=" * 80)

fig = plt.figure(figsize=(18, 12))
gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)

# Plot 1: Top 15 Single Interventions
ax1 = fig.add_subplot(gs[0, :2])
top15 = rankings_df.head(15)
colors = plt.cm.RdYlGn(top15['Success Rate (%)'] / 100)
bars = ax1.barh(range(len(top15)), top15['Cost-Effectiveness'], color=colors, edgecolor='black')
ax1.set_yticks(range(len(top15)))
ax1.set_yticklabels(top15['Intervention'], fontsize=9)
ax1.set_xlabel('Cost-Effectiveness (success rate per $1000)', fontsize=11)
ax1.set_title('Top 15 Single Interventions by Cost-Effectiveness', fontsize=13, fontweight='bold')
ax1.invert_yaxis()
ax1.grid(True, alpha=0.3, axis='x')

# Plot 2: Success Rate vs Cost (scatter)
ax2 = fig.add_subplot(gs[0, 2])
scatter = ax2.scatter(rankings_df['Cost'], rankings_df['Success Rate (%)'],
                     c=rankings_df['Cost-Effectiveness'], cmap='viridis',
                     s=100, alpha=0.6, edgecolors='black')
ax2.set_xlabel('Cost ($)', fontsize=10)
ax2.set_ylabel('Success Rate (%)', fontsize=10)
ax2.set_title('Cost vs Success Rate', fontsize=11, fontweight='bold')
ax2.grid(True, alpha=0.3)
plt.colorbar(scatter, ax=ax2, label='Cost-Effectiveness')

# Plot 3: Top 10 Combinations
ax3 = fig.add_subplot(gs[1, :2])
top10_combos = combined_df_sorted.head(10)
colors = plt.cm.RdYlGn(top10_combos['Success Rate (%)'] / 100)
ax3.barh(range(len(top10_combos)), top10_combos['Cost-Effectiveness'], color=colors, edgecolor='black')
ax3.set_yticks(range(len(top10_combos)))
labels = [s[:50] + '...' if len(s) > 50 else s for s in top10_combos['Strategy']]
ax3.set_yticklabels(labels, fontsize=8)
ax3.set_xlabel('Cost-Effectiveness (success rate per $1000)', fontsize=11)
ax3.set_title('Top 10 Combination Strategies', fontsize=13, fontweight='bold')
ax3.invert_yaxis()
ax3.grid(True, alpha=0.3, axis='x')

# Plot 4: Success Rate by Budget Level
ax4 = fig.add_subplot(gs[1, 2])
budget_analysis = combined_df_sorted.groupby('Budget').agg({
    'Success Rate (%)': 'mean',
    'Cost-Effectiveness': 'mean'
}).reset_index()
ax4.plot(budget_analysis['Budget'], budget_analysis['Success Rate (%)'],
         marker='o', linewidth=2, markersize=8, color='steelblue')
ax4.set_xlabel('Budget ($)', fontsize=10)
ax4.set_ylabel('Avg Success Rate (%)', fontsize=10)
ax4.set_title('Success Rate by Budget', fontsize=11, fontweight='bold')
ax4.grid(True, alpha=0.3)

# Plot 5: Intervention Type Distribution
ax5 = fig.add_subplot(gs[2, 0])
direction_counts = rankings_df['Direction'].value_counts()
ax5.pie(direction_counts.values, labels=['Decrease', 'Increase'], autopct='%1.1f%%',
        colors=['lightcoral', 'lightgreen'], startangle=90)
ax5.set_title('Intervention Direction\nDistribution', fontsize=11, fontweight='bold')

# Plot 6: Parameter Sensitivity (magnitude effect)
ax6 = fig.add_subplot(gs[2, 1])
# Get top 5 parameters and plot how their effectiveness changes with magnitude
top_params = rankings_df.groupby('Parameter')['Cost-Effectiveness'].max().nlargest(5).index
for param in top_params:
    param_data = rankings_df[rankings_df['Parameter'] == param].sort_values('Magnitude (%)')
    ax6.plot(param_data['Magnitude (%)'], param_data['Success Rate (%)'],
             marker='o', label=param[:15], linewidth=2)
ax6.set_xlabel('Magnitude (%)', fontsize=10)
ax6.set_ylabel('Success Rate (%)', fontsize=10)
ax6.set_title('Top Parameters: Magnitude Effect', fontsize=11, fontweight='bold')
ax6.legend(fontsize=8)
ax6.grid(True, alpha=0.3)

# Plot 7: Cost Distribution
ax7 = fig.add_subplot(gs[2, 2])
ax7.hist(rankings_df['Cost'], bins=30, color='steelblue', alpha=0.7, edgecolor='black')
ax7.set_xlabel('Intervention Cost ($)', fontsize=10)
ax7.set_ylabel('Frequency', fontsize=10)
ax7.set_title('Cost Distribution', fontsize=11, fontweight='bold')
ax7.grid(True, alpha=0.3, axis='y')

plt.savefig('intervention_rankings_complete.png', dpi=300, bbox_inches='tight')
print("Visualization saved to: intervention_rankings_complete.png")

# ========================================
# PART 5: INTERACTIVE BUDGET TOOL
# ========================================
print("\n" + "=" * 80)
print("INTERACTIVE BUDGET ANALYSIS")
print("=" * 80)

def analyze_budget(budget, top_n=5):
    """Show best strategies for a given budget"""
    feasible = combined_df_sorted[combined_df_sorted['Total Cost'] <= budget]

    print(f"\nBudget: ${budget}")
    print(f"Feasible strategies: {len(feasible)}")

    if len(feasible) > 0:
        print(f"\nTop {top_n} strategies:")
        for i, (_, row) in enumerate(feasible.head(top_n).iterrows(), 1):
            print(f"\n{i}. {row['Strategy']}")
            print(f"   Success Rate: {row['Success Rate (%)']:.1f}%")
            print(f"   Cost: ${row['Total Cost']:.0f} (${budget - row['Total Cost']:.0f} remaining)")
            print(f"   Cost-Effectiveness: {row['Cost-Effectiveness']:.3f}")
        return feasible
    else:
        print("No strategies found within budget!")
        return pd.DataFrame()

# Example budget analyses
# for example_budget in [2000, 5000, 10000]:
#     analyze_budget(example_budget, top_n=3)

print("\n" + "=" * 80)
print("ANALYSIS COMPLETE")
print("=" * 80)
print("\nGenerated files:")
print("  1. intervention_rankings.csv - All single interventions ranked")
print("  2. intervention_combinations.csv - All combination strategies ranked")
print("  3. intervention_rankings_complete.png - Comprehensive visualization")

plt.show()