"""
Diagnostic script to verify reaction rates at steady state.
This script runs a single security strategy simulation and examines the final state in detail.
"""
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..')))

import numpy as np
from pathlib import Path

# Import the simulator from the main script
exec(open("script_two_budget_adaptive.py").read(), globals())

def diagnostic_run():
    """Run diagnostic test of steady state mechanism."""

    print("="*80)
    print("STEADY STATE DIAGNOSTIC TEST")
    print("="*80)
    print()

    # Model file path
    model_file = Path(__file__).parent.parent.parent.parent / "networks" / "Conflict_Theory" / "Resource_Community_Insurgency_Loops_model3.txt"
    if not model_file.exists():
        print(f"ERROR: Model file not found: {model_file}")
        return

    # Create simulator with medium conflict scenario
    sim = ConflictSimulator(str(model_file))
    sim.set_government_strategy("security", quiet=True)
    sim.set_armed_groups_strategy("balanced", quiet=True)
    sim.set_climate(0.01)

    print(f"Running simulation with:")
    print(f"  Government Strategy: security")
    print(f"  Armed Groups Strategy: balanced")
    print(f"  Climate parameter: 0.01")
    print(f"  Duration: 10 years")
    print()

    # Run simulation
    sim.simulate(years=10, dt=0.01)

    print("="*80)
    print("FINAL STATE ANALYSIS")
    print("="*80)
    print()

    # Get final state
    final_state = sim.state_history[-1]
    print("Final State Values:")
    print(f"  AG_SL:  {final_state['AG_SL']:.6f}")
    print(f"  AG_RL:  {final_state['AG_RL']:.6f}")
    print(f"  V:      {final_state['V']:.6f}")
    print(f"  T:      {final_state['T']:.6f}")
    print(f"  Gov:    {final_state['Gov']:.6f}")
    print(f"  E:      {final_state['E']:.6f}")
    print(f"  WR_SL:  {final_state['WR_SL']:.6f}")
    print()

    # Compute final allocations
    E_alloc, Gov_alloc, V_alloc = sim.compute_allocations(final_state)

    print("Final Budget Allocations:")
    print(f"  E available:   {final_state['E']:.6f}")
    print(f"  Gov available: {final_state['Gov']:.6f}")
    print(f"  V available:   {final_state['V']:.6f}")
    print()
    print("  E_alloc:")
    for rxn, val in sorted(E_alloc.items()):
        print(f"    {rxn}: {val:.6f}")
    print("  Gov_alloc:")
    for rxn, val in sorted(Gov_alloc.items()):
        print(f"    {rxn}: {val:.6f}")
    print("  V_alloc:")
    for rxn, val in sorted(V_alloc.items()):
        print(f"    {rxn}: {val:.6f}")
    print()

    # Check kinetic constants for critical reactions
    print("="*80)
    print("KINETIC CONSTANTS FOR DECAY REACTIONS")
    print("="*80)
    print(f"  k25 (trust destruction): {sim.k['r25']}")
    print(f"  k26 (governance decay):  {sim.k['r26']}")
    print(f"  k31 (violence decay):    {sim.k['r31']}")
    print()

    # Manually compute rates for critical reactions
    print("="*80)
    print("MANUAL RATE COMPUTATION FOR CRITICAL REACTIONS")
    print("="*80)
    print()

    reactions_to_check = [
        ('r25', 'Trust destruction: V + T => ;'),
        ('r26', 'Governance decay: Gov => ;'),
        ('r27', 'Violence from desperation: WR_SL => WR_SL + V'),
        ('r28', 'Violence from AG attacks: AG_SL + SR_RL => ...'),
        ('r29', 'Violence from Gov-AG conflict: AG_SL + Gov => ...'),
        ('r30', 'Violence from displacement: AG_RL + Gov => ...'),
        ('r31', 'Violence decay: V => ;'),
        ('r22', 'Trust building from economy: SR_RL + E => SR_RL + T'),
        ('r23', 'Trust building from gov: Gov + SR_RL => ... + T'),
    ]

    for rxn, description in reactions_to_check:
        rate, consumption = sim.compute_reaction_rate(rxn, final_state, E_alloc, Gov_alloc, V_alloc)
        print(f"{rxn}: {description}")
        print(f"  Rate: {rate:.8f}")
        print(f"  Consumption: E={consumption['E']:.8f}, Gov={consumption['Gov']:.8f}, V={consumption['V']:.8f}")
        print()

    # Check rates from history
    print("="*80)
    print("RATES FROM HISTORY (Last 5 timesteps)")
    print("="*80)
    print()

    for rxn, description in reactions_to_check:
        print(f"{rxn}: {description}")
        last_5_rates = [sim.rate_history[i][rxn] for i in range(-5, 0)]
        print(f"  Rates: {[f'{r:.8f}' for r in last_5_rates]}")
        print()

    # Analyze V dynamics in detail
    print("="*80)
    print("VIOLENCE DYNAMICS BALANCE")
    print("="*80)
    print()

    producers = ['r27', 'r28', 'r29', 'r30']
    consumers = ['r15', 'r16', 'r25', 'r31']

    print("Violence Production:")
    total_production = 0
    for rxn in producers:
        rate = sim.rate_history[-1][rxn]
        total_production += rate
        print(f"  {rxn}: {rate:.8f}")
    print(f"  TOTAL: {total_production:.8f}")
    print()

    print("Violence Consumption:")
    total_consumption = 0
    for rxn in consumers:
        rate = sim.rate_history[-1][rxn]
        total_consumption += rate
        print(f"  {rxn}: {rate:.8f}")
    print(f"  TOTAL: {total_consumption:.8f}")
    print()

    print(f"Net V change rate: {total_production - total_consumption:.8f}")
    print()

    # Check if r25 implementation is correct
    print("="*80)
    print("r25 IMPLEMENTATION CHECK")
    print("="*80)
    print()

    k25 = sim.k['r25']
    V_val = final_state['V']
    T_val = final_state['T']
    V_alloc_r25 = V_alloc.get('r25', 0)

    print("r25: V + T => ; (Trust destruction by violence)")
    print()
    print("EXPECTED (pure mass action):")
    print(f"  rate = k25 × V × T")
    print(f"  rate = {k25} × {V_val:.6f} × {T_val:.6f}")
    print(f"  rate = {k25 * V_val * T_val:.8f}")
    print()
    print("CURRENT IMPLEMENTATION (budget-controlled):")
    print(f"  rate = k25 × V_alloc['r25'] × T")
    print(f"  rate = {k25} × {V_alloc_r25:.6f} × {T_val:.6f}")
    print(f"  rate = {k25 * V_alloc_r25 * T_val:.8f}")
    print()
    print("ACTUAL COMPUTED RATE:")
    rate_r25, _ = sim.compute_reaction_rate('r25', final_state, E_alloc, Gov_alloc, V_alloc)
    print(f"  rate = {rate_r25:.8f}")
    print()

    if abs(rate_r25 - (k25 * V_alloc_r25 * T_val)) < 1e-6:
        print("[BUG CONFIRMED] r25 is using budget allocation instead of pure mass action!")
        print("  This makes r25 inactive when V_alloc['r25'] is small,")
        print("  even though violence and trust both exist in the system.")
    else:
        print("[UNEXPECTED] Rate doesn't match either formula - investigate further.")
    print()

    # Check stability of variables
    print("="*80)
    print("VARIABLE STABILITY CHECK (Last 50 timesteps)")
    print("="*80)
    print()

    variables = ['V', 'T', 'Gov', 'E', 'AG_SL', 'AG_RL']
    for var in variables:
        values = [sim.state_history[i][var] for i in range(-50, 0)]
        mean_val = np.mean(values)
        std_val = np.std(values)
        rel_std = std_val / mean_val if mean_val != 0 else std_val

        print(f"{var}:")
        print(f"  Mean: {mean_val:.6f}")
        print(f"  Std:  {std_val:.6f}")
        print(f"  Rel Std: {rel_std:.6f} ({rel_std*100:.4f}%)")
        print(f"  Status: {'STABLE' if rel_std < 0.01 else 'VARYING'}")
        print()

if __name__ == "__main__":
    diagnostic_run()
