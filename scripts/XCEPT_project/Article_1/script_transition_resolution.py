"""
CONFLICT RESOLUTION ANALYSIS - Lake Chad Basin Model
=====================================================
Based on the Three Interlocking Loops Framework:
1. Resource-Scarcity Loop: Climate -> Degradation -> Resilience Loss
2. Scarcity-Conflict Loop: Scarcity -> Violence -> Trust Erosion -> Institutional Vacuum
3. Insurgency-Grievance Loop: Grievances -> Recruitment -> Violence -> State Withdrawal

Key insight: Single-focus interventions fail. Multi-pronged approaches required.

This analysis compares:
- Mono-view strategies (governance-only, economic-only, etc.) - EXPECTED TO FAIL
- Multi-pronged optimization - EXPECTED TO SUCCEED (but costly)
"""

import pulp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Any
from dataclasses import dataclass, field
import warnings
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))

try:
    from pyCOT.io.functions import read_txt
    HAVE_COT = True
except ImportError:
    print("ERROR: pyCOT not found.")
    sys.exit(1)

warnings.filterwarnings('ignore')

# ==================== LOAD REACTION NETWORK ====================

def load_reaction_network(file_path: str):
    """Load the reaction network"""
    if not os.path.exists(file_path):
        print(f"ERROR: File not found: {file_path}")
        return None, None, None, None
    try:
        rn = read_txt(file_path)
        S = rn.stoichiometry_matrix()
        species_names = list(rn._species_map.keys())
        reaction_names = [reaction.name() for reaction in rn.reactions()]
        print(f"Loaded: {len(species_names)} species, {len(reaction_names)} reactions")
        return S, species_names, reaction_names, rn
    except Exception as e:
        print(f"ERROR: {e}")
        return None, None, None, None

FILE_PATH = 'networks/Conflict_Theory/Resource_Community_Insurgency_Loops_model1.txt'

S, SPECIES_NAMES, REACTION_NAMES, RN = load_reaction_network(FILE_PATH)
if S is None:
    sys.exit(1)

RXN_IDX = {name: idx for idx, name in enumerate(REACTION_NAMES)}

# ==================== CONSTANTS ====================

POP_TOTAL = 40.00  # million

# Realistic AG levels
AG_HIGH = 0.02 * POP_TOTAL    # 0.8M (2%)
AG_MIDDLE = 0.01 * POP_TOTAL  # 0.4M (1%)
AG_LOW = 0.005 * POP_TOTAL    # 0.2M (0.5%)

CLIMATES = [0.1, 0.25, 0.5, 0.75]
BUDGETS = [200, 400, 800]
AG_TARGETS = [0.05, 0.1, 0.2]

CONFLICT_SPECIES = ['AG_RL', 'AG_SL', 'V', 'WR_RL', 'WR_SL']
WEIGHTS = {'AG_RL': 10.0, 'AG_SL': 10.0, 'V': 5.0, 'WR_RL': 1.0, 'WR_SL': 1.0}

# ==================== INTERVENTION STRATEGIES ====================

@dataclass
class InterventionStrategy:
    """
    Defines which reaction types can be activated by an intervention.

    Mono-view strategies focus on one aspect and typically fail.
    Multi-pronged strategies address all loops simultaneously.
    """
    name: str
    description: str

    # Which reaction categories are enabled (0-1 scale, 0=blocked, 1=full)
    governance_enabled: float = 1.0      # r9, r10, r21, r23 (institution building)
    economic_enabled: float = 1.0        # r7, r8, r16, r17 (economic production, resilience)
    trust_enabled: float = 1.0           # r21, r22, r24 (trust building, violence reduction)
    demobilization_enabled: float = 1.0  # r19 (DDR programs)
    security_enabled: float = 1.0        # Violence suppression capacity

    # Coupling parameters (from policy context)
    alpha: float = 0.4   # Trust efficiency
    beta: float = 0.5    # DDR-governance coupling
    gamma: float = 0.5   # DDR-economy coupling
    delta: float = 0.6   # Institution-economy coupling
    epsilon: float = 1.0 # External intervention limit

    def __str__(self):
        return f"{self.name}"


# Define intervention strategies
STRATEGIES = {
    # === MONO-VIEW STRATEGIES (expected to fail or be very costly) ===

    'governance_only': InterventionStrategy(
        name='Governance Focus',
        description='Military/institutional approach - build state capacity only',
        governance_enabled=1.0,
        economic_enabled=0.2,      # Minimal economic investment
        trust_enabled=0.2,         # Minimal trust building
        demobilization_enabled=0.3,
        security_enabled=1.0,
        alpha=0.2, beta=0.3, gamma=0.2, delta=0.3, epsilon=2.0
    ),

    'economic_only': InterventionStrategy(
        name='Economic Development',
        description='Development-first approach - jobs and livelihoods',
        governance_enabled=0.2,
        economic_enabled=1.0,
        trust_enabled=0.3,
        demobilization_enabled=0.3,
        security_enabled=0.2,
        alpha=0.2, beta=0.2, gamma=0.6, delta=0.8, epsilon=0.5
    ),

    'trust_only': InterventionStrategy(
        name='Social Reconciliation',
        description='Community-based peace building and dialogue',
        governance_enabled=0.3,
        economic_enabled=0.3,
        trust_enabled=1.0,
        demobilization_enabled=0.2,
        security_enabled=0.2,
        alpha=0.5, beta=0.2, gamma=0.2, delta=0.3, epsilon=0.5
    ),

    'security_only': InterventionStrategy(
        name='Military Response',
        description='Counter-insurgency operations only',
        governance_enabled=0.2,
        economic_enabled=0.1,
        trust_enabled=0.1,
        demobilization_enabled=0.1,
        security_enabled=1.0,
        alpha=0.1, beta=0.1, gamma=0.1, delta=0.2, epsilon=2.0
    ),

    # === MULTI-PRONGED STRATEGIES ===

    'multipronged_weak': InterventionStrategy(
        name='Multi-Pronged (Weak Context)',
        description='Balanced approach in fragile institutional environment',
        governance_enabled=0.8,
        economic_enabled=0.8,
        trust_enabled=0.8,
        demobilization_enabled=0.8,
        security_enabled=0.8,
        alpha=0.3, beta=0.4, gamma=0.4, delta=0.5, epsilon=1.5
    ),

    'multipronged_strong': InterventionStrategy(
        name='Multi-Pronged (Strong Context)',
        description='Balanced approach with strong local institutions',
        governance_enabled=1.0,
        economic_enabled=1.0,
        trust_enabled=1.0,
        demobilization_enabled=1.0,
        security_enabled=1.0,
        alpha=0.5, beta=0.6, gamma=0.6, delta=0.8, epsilon=1.0
    ),

    'multipronged_optimal': InterventionStrategy(
        name='Multi-Pronged (Optimal)',
        description='Fully integrated approach addressing all three loops',
        governance_enabled=1.0,
        economic_enabled=1.0,
        trust_enabled=1.0,
        demobilization_enabled=1.0,
        security_enabled=1.0,
        alpha=0.6, beta=0.7, gamma=0.7, delta=0.9, epsilon=1.0
    ),
}

# ==================== CASE DEFINITIONS ====================

def define_cases():
    """
    Define initial conflict states.
    Population conserved: SR_RL + SR_SL + WR_RL + WR_SL + AG_RL + AG_SL = 40M
    """
    return {
        'N1_High': {
            'SR_RL': 8.0, 'WR_RL': 6.0, 'AG_RL': 0.3,
            'SR_SL': 10.0, 'WR_SL': 15.2, 'AG_SL': 0.5,
            'V': 3.0, 'Gov': 0.5, 'E': 1.0, 'T': 0.2,
            'RL': 0.0, 'SL': 0.0
        },
        'N2_Mid': {
            'SR_RL': 10.0, 'WR_RL': 5.0, 'AG_RL': 0.15,
            'SR_SL': 12.0, 'WR_SL': 12.6, 'AG_SL': 0.25,
            'V': 2.0, 'Gov': 0.8, 'E': 1.5, 'T': 0.4,
            'RL': 0.0, 'SL': 0.0
        },
        'N3_Low': {
            'SR_RL': 12.0, 'WR_RL': 4.0, 'AG_RL': 0.08,
            'SR_SL': 14.0, 'WR_SL': 9.8, 'AG_SL': 0.12,
            'V': 1.0, 'Gov': 1.0, 'E': 2.0, 'T': 0.6,
            'RL': 0.0, 'SL': 0.0
        },
        'M1_High': {
            'SR_RL': 7.0, 'WR_RL': 7.0, 'AG_RL': 0.35,
            'SR_SL': 9.0, 'WR_SL': 16.2, 'AG_SL': 0.45,
            'V': 3.0, 'Gov': 0.6, 'E': 1.5, 'T': 0.3,
            'RL': 20.0, 'SL': 20.0
        },
        'M2_Mid': {
            'SR_RL': 10.0, 'WR_RL': 5.5, 'AG_RL': 0.18,
            'SR_SL': 11.0, 'WR_SL': 13.1, 'AG_SL': 0.22,
            'V': 2.0, 'Gov': 0.9, 'E': 2.0, 'T': 0.5,
            'RL': 15.0, 'SL': 25.0
        },
        'M3_Low': {
            'SR_RL': 12.0, 'WR_RL': 4.5, 'AG_RL': 0.08,
            'SR_SL': 13.0, 'WR_SL': 10.3, 'AG_SL': 0.12,
            'V': 1.2, 'Gov': 1.1, 'E': 2.5, 'T': 0.7,
            'RL': 25.0, 'SL': 15.0
        },
    }


# ==================== LP SOLVER WITH THREE LOOPS MODEL ====================

def solve_with_loops_model(
    case_data: Dict[str, float],
    climate: float,
    budget: float,
    target_ag: float,
    strategy: InterventionStrategy,
    S: np.ndarray,
    species_names: List[str]
) -> Dict[str, Any]:
    """
    LP solver implementing the Three Interlocking Loops framework.

    Key constraint: AG reduction REQUIRES coordinated intervention across:
    - Governance (state capacity for DDR)
    - Economics (livelihood alternatives)
    - Trust (community acceptance)
    - Security (prevent backsliding)

    Mono-view strategies fail because they cannot achieve the minimum
    required levels across all dimensions simultaneously.
    """

    prob = pulp.LpProblem("Conflict_Resolution_Loops", pulp.LpMinimize)

    # Setup
    species_to_idx = {name: idx for idx, name in enumerate(species_names)}
    present_species = [s for s in case_data.keys() if s in species_to_idx]

    if not present_species:
        return {'optimal': False, 'error': 'Species mismatch'}

    n_reactions = S.shape[1]

    # Decision variables
    flux_vars = [pulp.LpVariable(f"v_{j}", lowBound=0) for j in range(n_reactions)]
    target_vars = {sp: pulp.LpVariable(f"x_{sp}", lowBound=0) for sp in present_species}
    delta_vars = {sp: pulp.LpVariable(f"d_{sp}", lowBound=None) for sp in present_species}

    # ===== STOICHIOMETRIC CONSTRAINT =====
    for species in present_species:
        i = species_to_idx[species]
        reaction_effect = pulp.lpSum([S[i, j] * flux_vars[j] for j in range(n_reactions)])
        prob += delta_vars[species] == reaction_effect
        prob += target_vars[species] == case_data[species] + delta_vars[species]

    # ===== REACTION INDICES =====
    idx_r7 = RXN_IDX.get('r7', 6)     # SR economic production
    idx_r8 = RXN_IDX.get('r8', 7)     # WR economic production
    idx_r9 = RXN_IDX.get('r9', 8)     # Economy -> Governance
    idx_r10 = RXN_IDX.get('r10', 9)   # External governance
    idx_r16 = RXN_IDX.get('r16', 15)  # Resilience RL
    idx_r17 = RXN_IDX.get('r17', 16)  # Resilience SL
    idx_r19 = RXN_IDX.get('r19', 18)  # Demobilization: E + Gov + AG_SL => WR_SL
    idx_r21 = RXN_IDX.get('r21', 20)  # Gov + SR_RL => Gov + SR_RL + T
    idx_r22 = RXN_IDX.get('r22', 21)  # Gov + WR_RL => SR_RL
    idx_r24 = RXN_IDX.get('r24', 23)  # 2Gov =>
    idx_r30 = RXN_IDX.get('r30', 29)  # WR_SL + AG_SL + V => 2AG_SL (recruitment)
    idx_r33 = RXN_IDX.get('r33', 32)  # AG_SL + E => AG_RL

    # ===== INITIAL CONDITIONS =====
    initial_E = case_data.get('E', 1.0)
    initial_V = case_data.get('V', 1.0)
    initial_T = case_data.get('T', 0.1)
    initial_Gov = case_data.get('Gov', 0.5)
    initial_WR = case_data.get('WR_RL', 0) + case_data.get('WR_SL', 0)
    initial_AG = case_data.get('AG_RL', 0) + case_data.get('AG_SL', 0)

    # ===== LOOP FACTORS (computed from initial state) =====

    # Loop 1: Resource-Scarcity - how stressed is the system?
    scarcity_factor = min(1.0, (initial_WR / POP_TOTAL) / (initial_E + 0.1))

    # Loop 2: Violence erosion - how much does violence destroy capacity?
    violence_factor = initial_V / (initial_V + 1.0)

    # Loop 3: Grievance-recruitment - how attractive is joining AG?
    grievance_level = min(1.0, max(0.0,
        (1.0 - initial_Gov / 2.0) * (1.0 - initial_E / 3.0) * (initial_WR / 20.0)))

    # ===== STRATEGY CAPACITY LIMITS =====
    # Mono-view strategies have HARD LIMITS on non-focus areas

    # Economic capacity (for economic-focused and multi-pronged)
    econ_capacity = climate * strategy.economic_enabled
    prob += flux_vars[idx_r7] <= 5.0 * econ_capacity
    prob += flux_vars[idx_r8] <= 5.0 * econ_capacity
    prob += flux_vars[idx_r16] <= 5.0 * econ_capacity
    prob += flux_vars[idx_r17] <= 5.0 * econ_capacity

    # Governance capacity (for governance-focused and multi-pronged)
    gov_capacity = climate * strategy.governance_enabled
    prob += flux_vars[idx_r9] <= 5.0 * gov_capacity
    prob += flux_vars[idx_r10] <= strategy.epsilon * strategy.governance_enabled

    # Trust capacity (for trust-focused and multi-pronged)
    trust_capacity = climate * strategy.trust_enabled
    prob += flux_vars[idx_r21] <= 5.0 * trust_capacity
    prob += flux_vars[idx_r22] <= 5.0 * trust_capacity

    # Demobilization capacity
    demob_capacity = climate * strategy.demobilization_enabled
    prob += flux_vars[idx_r19] <= 5.0 * demob_capacity

    # ===== LOOP 1: RESOURCE-SCARCITY FEEDBACK =====
    # Scarcity reduces economic productivity (harder to generate E)
    scarcity_penalty = max(0.2, 1.0 - 0.5 * scarcity_factor)
    prob += flux_vars[idx_r7] + flux_vars[idx_r8] <= 8.0 * scarcity_penalty * econ_capacity

    # ===== LOOP 2: VIOLENCE DESTRUCTION FEEDBACK =====
    # Violence destroys trust-building capacity
    violence_trust_penalty = max(0.2, 1.0 - 0.6 * violence_factor)
    prob += flux_vars[idx_r21] <= 3.0 * violence_trust_penalty * trust_capacity

    # Violence destroys governance-building capacity
    violence_gov_penalty = max(0.2, 1.0 - 0.5 * violence_factor)
    prob += flux_vars[idx_r9] <= 3.0 * violence_gov_penalty * gov_capacity

    # ===== LOOP 3: GRIEVANCE-RECRUITMENT FEEDBACK =====
    # High grievances make demobilization MUCH harder
    grievance_penalty = max(0.15, 1.0 - 0.7 * grievance_level)

    # ===== CRITICAL COUPLING CONSTRAINTS =====
    # These are what make mono-view strategies FAIL

    # Define activity variables for coupling
    economic_activity = flux_vars[idx_r7] + flux_vars[idx_r8]
    governance_activity = flux_vars[idx_r9] + flux_vars[idx_r10]
    trust_activity = flux_vars[idx_r21] + flux_vars[idx_r22]
    demobilization = flux_vars[idx_r19]

    # C1: Demobilization requires BOTH governance AND economics
    # You can't just throw money at it (economic_only fails)
    # You can't just build institutions (governance_only fails)
    prob += demobilization <= strategy.beta * governance_activity * grievance_penalty
    prob += demobilization <= strategy.gamma * economic_activity * grievance_penalty

    # C2: Institution building requires economic foundation
    # Governance-only fails without economy
    prob += flux_vars[idx_r9] <= strategy.delta * economic_activity

    # C3: Trust building requires governance presence
    # Trust-only fails without governance
    prob += trust_activity <= 2.0 * (initial_Gov + governance_activity) * strategy.trust_enabled

    # C4: Violence reduction requires trust
    # Security-only fails because you can't sustain peace without trust
    prob += flux_vars[idx_r24] <= strategy.alpha * (initial_T + trust_activity)

    # C5: Economic development requires security
    # Economic-only fails in high violence environments
    violence_econ_penalty = max(0.3, 1.0 - 0.4 * violence_factor)
    prob += economic_activity <= 10.0 * violence_econ_penalty * strategy.economic_enabled

    # ===== MINIMUM COORDINATED INTERVENTION FOR AG REDUCTION =====
    # This is the key constraint that makes single-focus strategies fail
    # To reduce AG by X, you need minimum levels across ALL dimensions

    ag_reduction_requested = initial_AG - target_ag  # How much AG we want to eliminate

    if ag_reduction_requested > 0.01:  # If significant AG reduction is required
        # Minimum requirements scale with how much AG we want to eliminate
        min_factor = ag_reduction_requested / initial_AG if initial_AG > 0 else 0.5

        # You need MINIMUM governance activity (can't do DDR without state)
        prob += governance_activity >= 0.3 * min_factor * strategy.governance_enabled

        # You need MINIMUM economic activity (can't reintegrate without jobs)
        prob += economic_activity >= 0.4 * min_factor * strategy.economic_enabled

        # You need MINIMUM trust activity (communities must accept ex-combatants)
        prob += trust_activity >= 0.2 * min_factor * strategy.trust_enabled

        # You need MINIMUM security commitment (prevent recruitment backsliding)
        security_commitment = flux_vars[idx_r24] + flux_vars[idx_r22]
        prob += security_commitment >= 0.15 * min_factor * strategy.security_enabled

    # ===== BUDGET CONSTRAINT =====
    flux_costs = []
    for j in range(n_reactions):
        base_cost = 10.0 / climate
        flux_costs.append(base_cost * flux_vars[j])
    prob += pulp.lpSum(flux_costs) <= budget

    # ===== POPULATION CONSERVATION =====
    pop_species = ['SR_RL', 'WR_RL', 'AG_RL', 'SR_SL', 'WR_SL', 'AG_SL']
    present_pop = [s for s in pop_species if s in target_vars]
    if present_pop:
        prob += pulp.lpSum([target_vars[s] for s in present_pop]) == POP_TOTAL

    # ===== ORGANIZATION CONSTRAINTS =====
    scenario = 'migration' if case_data.get('RL', 0) > 0.1 else 'no_migration'
    if scenario == 'migration':
        if 'RL' in target_vars:
            prob += target_vars['RL'] >= 0.01
        if 'SL' in target_vars:
            prob += target_vars['SL'] >= 0.01
    else:
        if 'RL' in target_vars:
            prob += target_vars['RL'] <= 0.001
        if 'SL' in target_vars:
            prob += target_vars['SL'] <= 0.001

    # ===== AG TARGET =====
    if 'AG_RL' in target_vars and 'AG_SL' in target_vars:
        prob += target_vars['AG_RL'] + target_vars['AG_SL'] <= target_ag

    # ===== MINIMUM GOVERNANCE =====
    if 'Gov' in target_vars:
        prob += target_vars['Gov'] >= 0.1

    # ===== OBJECTIVE: MINIMIZE CONFLICT =====
    conflict_terms = []
    for species in CONFLICT_SPECIES:
        if species in target_vars:
            weight = WEIGHTS.get(species, 1.0)
            conflict_terms.append(weight * target_vars[species])
    prob += pulp.lpSum(conflict_terms)

    # ===== SOLVE =====
    solver = pulp.PULP_CBC_CMD(msg=False, timeLimit=30)
    status = prob.solve(solver)

    # ===== RESULTS =====
    results = {
        'status': pulp.LpStatus[status],
        'optimal': status == pulp.LpStatusOptimal,
        'objective': pulp.value(prob.objective) if status == pulp.LpStatusOptimal else None,
        'strategy': strategy.name,
        'scarcity_factor': scarcity_factor,
        'violence_factor': violence_factor,
        'grievance_level': grievance_level,
    }

    if results['optimal']:
        results['targets'] = {s: pulp.value(target_vars[s]) for s in present_species}
        results['deltas'] = {s: pulp.value(delta_vars[s]) for s in present_species}

        ag_initial = case_data.get('AG_RL', 0) + case_data.get('AG_SL', 0)
        ag_final = (pulp.value(target_vars.get('AG_RL')) or 0) + (pulp.value(target_vars.get('AG_SL')) or 0)
        results['ag_initial'] = ag_initial
        results['ag_final'] = ag_final
        results['ag_reduction'] = ag_initial - ag_final
        results['ag_reduction_pct'] = (ag_initial - ag_final) / ag_initial * 100 if ag_initial > 0 else 0

        # Budget used
        flux_cost = sum((10.0 / climate) * (pulp.value(flux_vars[j]) or 0) for j in range(n_reactions))
        results['budget_used'] = min(flux_cost, budget)

        # Key activities
        results['economic_activity'] = (pulp.value(flux_vars[idx_r7]) or 0) + (pulp.value(flux_vars[idx_r8]) or 0)
        results['governance_activity'] = (pulp.value(flux_vars[idx_r9]) or 0) + (pulp.value(flux_vars[idx_r10]) or 0)
        results['trust_activity'] = (pulp.value(flux_vars[idx_r21]) or 0) + (pulp.value(flux_vars[idx_r22]) or 0)
        results['demobilization_activity'] = pulp.value(flux_vars[idx_r19]) or 0

    return results


# ==================== ANALYSIS ====================

def run_analysis():
    """Run comprehensive strategy comparison"""

    print("\n" + "="*70)
    print("CONFLICT RESOLUTION ANALYSIS - Three Loops Framework")
    print("="*70)
    print("\nComparing mono-view vs multi-pronged intervention strategies")

    cases = define_cases()
    print(f"\nCases: {len(cases)}")

    print(f"\nStrategies: {len(STRATEGIES)}")
    mono_strategies = [s for s in STRATEGIES.keys() if 'only' in s]
    multi_strategies = [s for s in STRATEGIES.keys() if 'multipronged' in s]
    print(f"  Mono-view: {mono_strategies}")
    print(f"  Multi-pronged: {multi_strategies}")

    all_results = []
    total = len(cases) * len(STRATEGIES) * len(CLIMATES) * len(BUDGETS) * len(AG_TARGETS)
    current = 0

    print(f"\nTotal LP problems: {total}")
    print("Running...")

    for case_name, case_data in cases.items():
        scenario = 'migration' if 'M' in case_name else 'no_migration'
        ag_level = 'high' if 'High' in case_name else ('mid' if 'Mid' in case_name else 'low')

        for strat_name, strategy in STRATEGIES.items():
            strategy_type = 'mono' if 'only' in strat_name else 'multi'

            for climate in CLIMATES:
                for budget in BUDGETS:
                    for target_ag in AG_TARGETS:
                        current += 1
                        if current % 100 == 0:
                            print(f"  {current}/{total}")

                        try:
                            result = solve_with_loops_model(
                                case_data, climate, budget, target_ag,
                                strategy, S, SPECIES_NAMES
                            )

                            record = {
                                'case': case_name,
                                'scenario': scenario,
                                'ag_level': ag_level,
                                'strategy': strat_name,
                                'strategy_name': strategy.name,
                                'strategy_type': strategy_type,
                                'climate': climate,
                                'budget': budget,
                                'target_ag': target_ag,
                                'feasible': result['optimal'],
                                'objective': result.get('objective'),
                                'ag_initial': result.get('ag_initial', 0),
                                'ag_final': result.get('ag_final', 0),
                                'ag_reduction': result.get('ag_reduction', 0),
                                'ag_reduction_pct': result.get('ag_reduction_pct', 0),
                                'budget_used': result.get('budget_used', 0),
                                'scarcity_factor': result.get('scarcity_factor', 0),
                                'violence_factor': result.get('violence_factor', 0),
                                'grievance_level': result.get('grievance_level', 0),
                                'economic_activity': result.get('economic_activity', 0),
                                'governance_activity': result.get('governance_activity', 0),
                                'trust_activity': result.get('trust_activity', 0),
                                'demobilization_activity': result.get('demobilization_activity', 0),
                            }
                            all_results.append(record)

                        except Exception as e:
                            all_results.append({
                                'case': case_name,
                                'scenario': scenario,
                                'strategy': strat_name,
                                'strategy_type': strategy_type,
                                'climate': climate,
                                'budget': budget,
                                'target_ag': target_ag,
                                'feasible': False,
                                'error': str(e)[:50]
                            })

    df = pd.DataFrame(all_results)

    print("\n" + "="*70)
    print("RESULTS SUMMARY")
    print("="*70)

    print(f"\nTotal runs: {len(df)}")
    print(f"Feasible: {df['feasible'].sum()} ({df['feasible'].mean()*100:.1f}%)")

    return df


def print_strategy_comparison(df: pd.DataFrame):
    """Print comparison of mono vs multi-pronged strategies"""

    print("\n" + "="*70)
    print("STRATEGY COMPARISON: Mono-View vs Multi-Pronged")
    print("="*70)

    # Overall feasibility
    print("\n1. FEASIBILITY BY STRATEGY TYPE:")
    print("-"*50)
    for stype in ['mono', 'multi']:
        data = df[df['strategy_type'] == stype]
        rate = data['feasible'].mean() * 100
        print(f"  {stype:12s}: {rate:5.1f}% feasible")

    # By individual strategy
    print("\n2. FEASIBILITY BY STRATEGY:")
    print("-"*50)
    strat_feas = df.groupby('strategy_name')['feasible'].mean().sort_values() * 100
    for strat, rate in strat_feas.items():
        print(f"  {strat:30s}: {rate:5.1f}%")

    feasible = df[df['feasible']].copy()
    if len(feasible) == 0:
        print("\nNo feasible solutions found!")
        return

    # Cost effectiveness
    print("\n3. COST-EFFECTIVENESS (feasible solutions only):")
    print("-"*50)
    feasible['cost_eff'] = feasible['ag_reduction'] / feasible['budget_used'].clip(lower=1)

    for stype in ['mono', 'multi']:
        data = feasible[feasible['strategy_type'] == stype]
        if len(data) > 0:
            avg_cost = data['budget_used'].mean()
            avg_red = data['ag_reduction_pct'].mean()
            print(f"  {stype:12s}: avg budget={avg_cost:.0f}, avg reduction={avg_red:.1f}%")

    # Climate sensitivity
    print("\n4. CLIMATE SENSITIVITY:")
    print("-"*50)
    for climate in sorted(df['climate'].unique()):
        mono_rate = df[(df['strategy_type'] == 'mono') & (df['climate'] == climate)]['feasible'].mean() * 100
        multi_rate = df[(df['strategy_type'] == 'multi') & (df['climate'] == climate)]['feasible'].mean() * 100
        print(f"  climate={climate}: mono={mono_rate:5.1f}%, multi={multi_rate:5.1f}%")

    # High conflict cases
    print("\n5. HIGH CONFLICT CASES (AG=2%):")
    print("-"*50)
    high = df[df['ag_level'] == 'high']
    for stype in ['mono', 'multi']:
        data = high[high['strategy_type'] == stype]
        rate = data['feasible'].mean() * 100
        print(f"  {stype:12s}: {rate:5.1f}% feasible")

    # Activity balance in successful multi-pronged
    print("\n6. ACTIVITY BALANCE (multi-pronged, feasible):")
    print("-"*50)
    multi_feas = feasible[feasible['strategy_type'] == 'multi']
    if len(multi_feas) > 0:
        for activity in ['economic_activity', 'governance_activity', 'trust_activity', 'demobilization_activity']:
            avg = multi_feas[activity].mean()
            print(f"  {activity:25s}: {avg:.2f}")


def create_visualizations(df: pd.DataFrame):
    """Create policy comparison visualizations"""

    feasible = df[df['feasible']].copy()

    fig = plt.figure(figsize=(18, 14))

    # 1. Feasibility comparison
    ax1 = plt.subplot(2, 3, 1)
    strat_feas = df.groupby(['strategy_name', 'strategy_type'])['feasible'].mean() * 100
    strat_feas = strat_feas.reset_index()
    colors = ['#e74c3c' if t == 'mono' else '#27ae60' for t in strat_feas['strategy_type']]
    bars = ax1.barh(range(len(strat_feas)), strat_feas['feasible'], color=colors)
    ax1.set_yticks(range(len(strat_feas)))
    ax1.set_yticklabels(strat_feas['strategy_name'], fontsize=9)
    ax1.set_xlabel('Feasibility (%)')
    ax1.set_title('Peace Achievability by Strategy\n(Red=Mono-view, Green=Multi-pronged)', fontweight='bold')
    ax1.axvline(50, color='gray', linestyle='--', alpha=0.5)

    # 2. Climate impact by strategy type
    ax2 = plt.subplot(2, 3, 2)
    for stype, color, marker in [('mono', '#e74c3c', 'o'), ('multi', '#27ae60', 's')]:
        data = df[df['strategy_type'] == stype]
        climate_feas = data.groupby('climate')['feasible'].mean() * 100
        ax2.plot(climate_feas.index, climate_feas.values, f'{marker}-',
                color=color, linewidth=2, markersize=8, label=stype.title())
    ax2.set_xlabel('Climate Favorability')
    ax2.set_ylabel('Feasibility (%)')
    ax2.set_title('Climate Sensitivity by Strategy Type', fontweight='bold')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    # 3. Budget effectiveness
    ax3 = plt.subplot(2, 3, 3)
    if len(feasible) > 0:
        for stype, color in [('mono', '#e74c3c'), ('multi', '#27ae60')]:
            data = feasible[feasible['strategy_type'] == stype]
            if len(data) > 0:
                budget_red = data.groupby('budget')['ag_reduction_pct'].mean()
                ax3.plot(budget_red.index, budget_red.values, 'o-',
                        color=color, linewidth=2, markersize=8, label=stype.title())
    ax3.set_xlabel('Budget')
    ax3.set_ylabel('AG Reduction (%)')
    ax3.set_title('Budget Effectiveness', fontweight='bold')
    ax3.legend()
    ax3.grid(True, alpha=0.3)

    # 4. Activity balance (multi-pronged)
    ax4 = plt.subplot(2, 3, 4)
    if len(feasible) > 0:
        multi_feas = feasible[feasible['strategy_type'] == 'multi']
        if len(multi_feas) > 0:
            activities = ['economic_activity', 'governance_activity', 'trust_activity', 'demobilization_activity']
            act_names = ['Economic', 'Governance', 'Trust', 'Demobilization']
            act_means = [multi_feas[a].mean() for a in activities]

            colors_act = ['#3498db', '#9b59b6', '#f39c12', '#1abc9c']
            ax4.bar(act_names, act_means, color=colors_act)
            ax4.set_ylabel('Average Activity Level')
            ax4.set_title('Activity Balance in Multi-Pronged Strategy', fontweight='bold')
            plt.xticks(rotation=45, ha='right')

    # 5. Conflict level impact
    ax5 = plt.subplot(2, 3, 5)
    for ag_level in ['high', 'mid', 'low']:
        data = df[df['ag_level'] == ag_level]
        mono_rate = data[data['strategy_type'] == 'mono']['feasible'].mean() * 100
        multi_rate = data[data['strategy_type'] == 'multi']['feasible'].mean() * 100
        x = ['mono', 'multi']
        ax5.plot(x, [mono_rate, multi_rate], 'o-', linewidth=2, markersize=10, label=f'AG={ag_level}')
    ax5.set_ylabel('Feasibility (%)')
    ax5.set_title('Strategy Effectiveness by Conflict Level', fontweight='bold')
    ax5.legend()
    ax5.grid(True, alpha=0.3)

    # 6. Key finding summary
    ax6 = plt.subplot(2, 3, 6)
    ax6.axis('off')

    mono_feas = df[df['strategy_type'] == 'mono']['feasible'].mean() * 100
    multi_feas_rate = df[df['strategy_type'] == 'multi']['feasible'].mean() * 100

    summary = [
        "KEY FINDINGS",
        "-" * 40,
        "",
        f"Mono-view strategies: {mono_feas:.1f}% feasible",
        f"Multi-pronged strategies: {multi_feas_rate:.1f}% feasible",
        "",
        "Mono-view approaches fail because:",
        "  - Single-focus cannot break all loops",
        "  - Violence destroys other investments",
        "  - Grievances regenerate quickly",
        "",
        "Multi-pronged succeeds because:",
        "  - Addresses all three loops simultaneously",
        "  - Economic + Governance + Trust + Security",
        "  - Breaks reinforcing feedback cycles",
    ]

    ax6.text(0.05, 0.95, '\n'.join(summary), fontsize=11, va='top', family='monospace',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    plt.tight_layout()
    plt.savefig('strategy_comparison_results.png', dpi=150, bbox_inches='tight')
    print("\nVisualization saved to 'strategy_comparison_results.png'")
    plt.show()


# ==================== MAIN ====================

if __name__ == "__main__":

    print("\n" + "="*70)
    print("LAKE CHAD BASIN CONFLICT MODEL")
    print("Three Interlocking Loops Analysis")
    print("="*70)

    # Run analysis
    results_df = run_analysis()

    # Save results
    results_df.to_csv('strategy_comparison_results.csv', index=False)
    print(f"\nResults saved to 'strategy_comparison_results.csv'")

    # Print comparison
    print_strategy_comparison(results_df)

    # Visualize
    create_visualizations(results_df)

    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)
