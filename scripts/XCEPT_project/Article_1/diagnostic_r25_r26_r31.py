"""
Enhanced diagnostic script to investigate r25 bug fix and r26/r31 inactivity.

This script:
1. Verifies r25 fix (now pure mass action)
2. Investigates why r26 and r31 become inactive in security strategy
3. Tracks the cascade of reactions when AG_SL -> 0
"""
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..')))

import numpy as np
from pathlib import Path

def run_enhanced_diagnostics():
    """Run enhanced diagnostics with the fixed r25 implementation."""

    print("="*80)
    print("ENHANCED DIAGNOSTIC: r25 FIX VERIFICATION & r26/r31 INVESTIGATION")
    print("="*80)
    print()

    # Model file path
    project_root = Path(__file__).parent.parent.parent.parent
    model_file = project_root / "networks" / "Conflict_Theory" / "Resource_Community_Insurgency_Loops_model3.txt"

    if not model_file.exists():
        print(f"ERROR: Model file not found: {model_file}")
        return

    # Import the ConflictSimulator class from the main script
    import importlib.util
    spec = importlib.util.spec_from_file_location("simulator",
                                                   Path(__file__).parent / "script_two_budget_adaptive.py")
    simulator_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(simulator_module)
    ConflictSimulator = simulator_module.ConflictSimulator

    # Run simulation with security strategy
    print("Setting up simulation...")
    print(f"  Model: {model_file.name}")
    print(f"  Strategy: security")
    print(f"  Duration: 10 years")
    print(f"  Timestep: 0.01 years")
    print()

    sim = ConflictSimulator(str(model_file))
    sim.set_government_strategy("security", quiet=True)
    sim.set_armed_groups_strategy("balanced", quiet=True)

    # Get initial state (medium conflict)
    initial_state = {
        "SR_RL": 30.0, "SR_SL": 20.0, "WR_RL": 25.0, "WR_SL": 50.0,
        "AG_RL": 0.8, "AG_SL": 2.5, "RL": 40.0, "SL": 60.0,
        "E": 80.0, "T": 40.0, "V": 30.0, "Gov": 20.0,
    }
    sim.state = initial_state
    sim.state_history = [initial_state.copy()]
    sim.time_history = [0.0]
    sim.rate_history = []
    sim.allocation_history = []
    sim.consumption_history = []
    sim.generation_history = []

    print("Running simulation...")
    dt = 0.01
    total_time = 10.0
    steps = int(total_time / dt)

    # Track when AG_SL crosses thresholds
    ag_sl_threshold_time = None
    ag_sl_zero_time = None

    # Track detailed diagnostics at key timepoints
    diagnostic_timepoints = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
    diagnostics_log = []

    for step in range(steps):
        # Check AG_SL thresholds
        current_ag_sl = sim.state.get("AG_SL", 0)
        if ag_sl_threshold_time is None and current_ag_sl < 0.5:
            ag_sl_threshold_time = sim.time
            print(f"  AG_SL dropped below 0.5 at t={sim.time:.3f}")
        if ag_sl_zero_time is None and current_ag_sl < 0.001:
            ag_sl_zero_time = sim.time
            print(f"  AG_SL ~0 at t={sim.time:.3f}")

        # Capture diagnostics at key timepoints
        if any(abs(sim.time - tp) < dt/2 for tp in diagnostic_timepoints):
            # Compute process vector
            process_vector, allocation_info, total_consumption, total_generation = \
                sim.compute_process_vector(sim.state, verbose=False)

            # Get decay reaction diagnostics
            decay_diag = sim.diagnose_decay_reactions(sim.state, process_vector, scale_factor=1.0)

            diagnostics_log.append({
                "time": sim.time,
                "state": sim.state.copy(),
                "decay_diagnostics": decay_diag,
                "process_vector": process_vector.copy(),
            })

        # Apply timestep
        sim.apply_timestep(dt, verbose=False)

    print(f"  Simulation complete!")
    print()

    # ========================================================================
    # ANALYSIS: r25 FIX VERIFICATION
    # ========================================================================
    print("="*80)
    print("PART 1: r25 FIX VERIFICATION")
    print("="*80)
    print()

    print("r25: V + T => ; (Trust destruction by violence)")
    print("  Expected behavior: Pure mass action (rate = k25 × V × T)")
    print("  Previous bug: Was budget-controlled (rate = k25 × V_alloc['r25'] × T)")
    print()

    for diag_entry in diagnostics_log:
        t = diag_entry["time"]
        decay_diag = diag_entry["decay_diagnostics"]
        r25_diag = decay_diag["r25"]

        if r25_diag["match"]:
            status = "[OK]"
        else:
            status = "[MISMATCH]"

        print(f"  t={t:5.1f}: V={r25_diag['V']:6.2f}, T={r25_diag['T']:6.2f}, "
              f"expected={r25_diag['expected_rate']:8.4f}, "
              f"actual={r25_diag['actual_rate']:8.4f} {status}")

    print()

    # ========================================================================
    # ANALYSIS: r26 and r31 INVESTIGATION
    # ========================================================================
    print("="*80)
    print("PART 2: r26 AND r31 INVESTIGATION")
    print("="*80)
    print()

    print("r26: Gov => ; (Governance decay)")
    print("  Expected: Always active when Gov > 0 (rate = k26 × Gov)")
    print()

    for diag_entry in diagnostics_log:
        t = diag_entry["time"]
        decay_diag = diag_entry["decay_diagnostics"]
        r26_diag = decay_diag["r26"]

        if r26_diag["match"]:
            status = "[OK]"
        else:
            status = "[MISMATCH]"

        is_active = r26_diag["actual_rate"] > 1e-10
        activity = "ACTIVE  " if is_active else "INACTIVE"

        print(f"  t={t:5.1f}: Gov={r26_diag['Gov']:6.2f}, "
              f"expected={r26_diag['expected_rate']:8.4f}, "
              f"actual={r26_diag['actual_rate']:8.4f} {status} {activity}")

    print()
    print("r31: V => ; (Violence decay)")
    print("  Expected: Always active when V > 0 (rate = k31 × V)")
    print()

    for diag_entry in diagnostics_log:
        t = diag_entry["time"]
        decay_diag = diag_entry["decay_diagnostics"]
        r31_diag = decay_diag["r31"]

        if r31_diag["match"]:
            status = "[OK]"
        else:
            status = "[MISMATCH]"

        is_active = r31_diag["actual_rate"] > 1e-10
        activity = "ACTIVE  " if is_active else "INACTIVE"

        print(f"  t={t:5.1f}: V={r31_diag['V']:6.2f}, "
              f"expected={r31_diag['expected_rate']:8.4f}, "
              f"actual={r31_diag['actual_rate']:8.4f} {status} {activity}")

    print()

    # ========================================================================
    # ANALYSIS: WHEN DOES INACTIVITY START?
    # ========================================================================
    print("="*80)
    print("PART 3: TIMELINE OF INACTIVITY")
    print("="*80)
    print()

    print("Tracking when reactions become inactive (rate < 1e-10):")
    print()

    # Check each reaction for inactivity transition
    reactions_to_track = ["r25", "r26", "r31"]

    for rxn in reactions_to_track:
        print(f"{rxn}:")
        was_active = True
        for diag_entry in diagnostics_log:
            t = diag_entry["time"]
            decay_diag = diag_entry["decay_diagnostics"]
            rxn_diag = decay_diag[rxn]

            is_active = rxn_diag["actual_rate"] > 1e-10

            if was_active and not is_active:
                print(f"  Became INACTIVE at t={t:.1f}")
                was_active = False
            elif not was_active and is_active:
                print(f"  Became ACTIVE again at t={t:.1f}")
                was_active = True

        if was_active:
            print(f"  Remained ACTIVE throughout simulation")
        else:
            print(f"  Remained INACTIVE after transition")
        print()

    # ========================================================================
    # ANALYSIS: REACTION CASCADE WHEN AG_SL -> 0
    # ========================================================================
    print("="*80)
    print("PART 4: REACTION CASCADE ANALYSIS")
    print("="*80)
    print()

    print("When AG_SL -> 0, these reactions should become inactive:")
    print("  r17: AG_SL + RL => ...")
    print("  r21: E + Gov + AG_SL => ...")
    print("  r28: AG_SL + SR_RL => ... + V  (violence generation)")
    print("  r29: AG_SL + Gov => ... + V  (violence generation)")
    print("  r32: WR_SL + AG_SL + V => 2AG_SL  (recruitment)")
    print("  r33: AG_SL + WR_RL => ...")
    print("  r34: AG_SL + SR_RL => ...")
    print("  r35: AG_SL + E => ...")
    print()

    ag_dependent_reactions = ["r17", "r21", "r28", "r29", "r32", "r33", "r34", "r35"]

    print("Checking AG_SL-dependent reactions at key timepoints:")
    print()

    for diag_entry in diagnostics_log:
        t = diag_entry["time"]
        state = diag_entry["state"]
        pv = diag_entry["process_vector"]

        ag_sl = state.get("AG_SL", 0)

        # Count active/inactive
        active_count = 0
        inactive_count = 0

        for rxn in ag_dependent_reactions:
            idx = sim.REACTIONS.index(rxn)
            rate = pv[idx]
            if rate > 1e-10:
                active_count += 1
            else:
                inactive_count += 1

        print(f"  t={t:5.1f}: AG_SL={ag_sl:6.3f}, "
              f"Active={active_count}, Inactive={inactive_count}")

    print()

    # ========================================================================
    # FINAL SUMMARY
    # ========================================================================
    print("="*80)
    print("SUMMARY")
    print("="*80)
    print()

    final_diag = diagnostics_log[-1]
    final_decay = final_diag["decay_diagnostics"]

    print("r25 Fix Status:")
    if final_decay["r25"]["match"]:
        print("  [OK] r25 now uses pure mass action correctly")
    else:
        print("  [ERROR] r25 still has implementation issues")
    print()

    print("r26 Status:")
    if final_decay["r26"]["actual_rate"] > 1e-10:
        print("  [OK] r26 remains active")
    else:
        print(f"  [ISSUE] r26 is INACTIVE despite Gov={final_decay['r26']['Gov']:.2f}")
        print(f"         Expected rate: {final_decay['r26']['expected_rate']:.6f}")
        print(f"         Actual rate: {final_decay['r26']['actual_rate']:.6f}")
    print()

    print("r31 Status:")
    if final_decay["r31"]["actual_rate"] > 1e-10:
        print("  [OK] r31 remains active")
    else:
        print(f"  [ISSUE] r31 is INACTIVE despite V={final_decay['r31']['V']:.2f}")
        print(f"         Expected rate: {final_decay['r31']['expected_rate']:.6f}")
        print(f"         Actual rate: {final_decay['r31']['actual_rate']:.6f}")
    print()

    print("="*80)
    print("END OF DIAGNOSTIC")
    print("="*80)


if __name__ == "__main__":
    run_enhanced_diagnostics()
