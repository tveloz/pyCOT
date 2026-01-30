"""
CALIBRATED CONFLICT MODEL - Lake Chad Basin
============================================
Based on real-world data and literature:

UNITS (per year):
- Population: 100,000 people per unit
- Budget: $1 million USD per unit
- Economy (E): $100 million GDP-equivalent per unit
- Violence (V): 1,000 conflict deaths per unit
- Governance (Gov): 0-1 scale (0=absent, 1=full state presence)
- Trust (T): 0-1 scale (0=no trust, 1=high social cohesion)
- Land (RL, SL): 100,000 hectares per unit

TIME FRAME: 1 year per optimization step

CALIBRATION SOURCES:
- UN OCHA Lake Chad Basin appeals (~$2B/year requested)
- Boko Haram estimated strength: 15,000-25,000 fighters
- DDR costs: $3,000-$10,000 per ex-combatant (World Bank, UNDP)
- Regional population: ~17 million in conflict-affected areas
- Displacement: 2-3 million people
"""

import pulp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Any
from dataclasses import dataclass
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))

from pyCOT.io.functions import read_txt

# ==================== CALIBRATION CONSTANTS ====================

# Population scale: 100,000 people per unit
POP_UNIT = 100_000  # people per unit
POP_TOTAL_UNITS = 170  # 17 million people total

# Budget scale: $1 million per unit
BUDGET_UNIT = 1_000_000  # USD per unit

# Real-world DDR costs (USD per person)
DDR_COST_PER_PERSON = 5_000  # $5,000 per ex-combatant (mid-range estimate)

# Real-world program costs ($ million per year for meaningful impact)
COST_DDR_PROGRAM = 50       # $50M/year for DDR program processing ~10,000
COST_GOVERNANCE = 100       # $100M/year for governance capacity building
COST_ECONOMIC = 30          # $30M/year for livelihood programs
COST_TRUST = 20             # $20M/year for reconciliation programs
COST_SECURITY = 150         # $150M/year for security operations

# AG levels as fraction of population (literature-based)
AG_SEVERE = 0.0015    # 0.15% = 25,500 fighters (Boko Haram peak)
AG_MODERATE = 0.0006  # 0.06% = 10,200 fighters
AG_LOW = 0.0002       # 0.02% = 3,400 fighters

# ==================== LOAD NETWORK ====================

FILE_PATH = 'networks/Conflict_Theory/Resource_Community_Insurgency_Loops_model1.txt'

def load_network():
    rn = read_txt(FILE_PATH)
    S = rn.stoichiometry_matrix()
    species = list(rn._species_map.keys())
    reactions = [r.name() for r in rn.reactions()]
    return S, species, reactions

S, SPECIES, REACTIONS = load_network()
RXN_IDX = {name: idx for idx, name in enumerate(REACTIONS)}
SP_IDX = {name: idx for idx, name in enumerate(SPECIES)}

# ==================== REACTION COST STRUCTURE ====================

def get_reaction_costs():
    """
    Assign realistic costs to each reaction type.
    Costs in budget units ($1M) per unit of flux.

    Flux interpretation:
    - Population reactions: flux=1 means 100,000 people affected
    - For DDR: flux=0.1 means 10,000 demobilized
    """

    costs = {}

    # Default cost (moderate)
    for rxn in REACTIONS:
        costs[rxn] = 10.0  # $10M per flux unit

    # DDR reactions - EXPENSIVE
    # r19: E + Gov + AG_SL => WR_SL (demobilization)
    # flux=0.1 = 10,000 people, cost = 10,000 * $5,000 = $50M
    costs['r19'] = 500.0  # $500M per 0.1 flux unit = $50M per 0.01 (1,000 people)

    # Governance building - EXPENSIVE
    # r9: 2E => Gov (build institutions from economy)
    # r10: => Gov (external governance injection)
    costs['r9'] = 100.0   # $100M per unit of governance built
    costs['r10'] = 150.0  # External intervention more expensive

    # Economic production - MODERATE
    # r7: SR_RL => SR_SL + E
    # r8: WR_RL => WR_SL + E
    costs['r7'] = 20.0    # $20M per unit of economic activity
    costs['r8'] = 30.0    # WR less efficient

    # Resilience building - MODERATE
    # r16: E + WR_RL => SR_RL
    # r17: E + WR_SL => SR_SL
    costs['r16'] = 40.0   # $40M to convert 100,000 WR to SR
    costs['r17'] = 40.0

    # Trust building - CHEAPER (community programs)
    # r20: SR_RL + 2E => SR_RL + 2E + T
    # r21: Gov + SR_RL => Gov + SR_RL + T
    costs['r20'] = 15.0   # $15M per trust unit
    costs['r21'] = 20.0

    # Violence-related - VERY EXPENSIVE (military operations)
    # r24: 2Gov => (governance decay from violence)
    # r29: 2V => (violence decay)
    costs['r24'] = 5.0    # Natural decay, low cost
    costs['r29'] = 5.0

    # AG recruitment/violence - AG BUDGET
    # r25-r28: Violence generation
    # r30: WR_SL + AG_SL + V => 2AG_SL (recruitment)
    # r31-r34: AG expansion
    for rxn in ['r25', 'r26', 'r27', 'r28', 'r30', 'r31', 'r32', 'r33', 'r34']:
        costs[rxn] = 30.0  # AG operations cost

    # Land/migration - LOW COST (natural processes)
    for rxn in ['r1', 'r1p', 'r1q', 'r2', 'r3', 'r4', 'r5', 'r6',
                'r11', 'r12', 'r13', 'r14', 'r15']:
        costs[rxn] = 2.0  # Mostly uncontrolled

    return costs

REACTION_COSTS = get_reaction_costs()

# ==================== CASE DEFINITIONS (CALIBRATED) ====================

def define_calibrated_cases():
    """
    Realistic scenarios based on Lake Chad Basin data.
    All population values in units of 100,000 people.
    Total population = 170 units (17 million)
    """

    return {
        # Severe: Peak Boko Haram (2014-2015 levels)
        # 25,000 fighters, 3M displaced, state nearly absent
        'Severe': {
            # Population distribution (total = 170 units = 17M people)
            'SR_RL': 30.0,   # 3M strong-resilient on resource land
            'SR_SL': 40.0,   # 4M strong-resilient on subsistence land
            'WR_RL': 25.0,   # 2.5M weak-resilient on resource land
            'WR_SL': 50.0,   # 5M weak-resilient on subsistence land (includes displaced)
            'AG_RL': 0.10,   # 10,000 AG on resource land
            'AG_SL': 0.15,   # 15,000 AG on subsistence land
            # Total AG = 0.25 units = 25,000 fighters (0.15% of population)

            # Other state variables
            'V': 8.0,        # 8,000 deaths/year (peak violence)
            'Gov': 0.15,     # 15% state presence (near collapse)
            'E': 10.0,       # $1B GDP (severely depressed)
            'T': 0.05,       # 5% trust (destroyed social fabric)
            'RL': 50.0,      # 5M hectares resource land
            'SL': 100.0,     # 10M hectares subsistence land
        },

        # Moderate: Current situation (~2020-2023)
        # 10,000 fighters, 2M displaced, fragile state
        'Moderate': {
            'SR_RL': 45.0,   # 4.5M
            'SR_SL': 55.0,   # 5.5M
            'WR_RL': 20.0,   # 2M
            'WR_SL': 35.0,   # 3.5M
            'AG_RL': 0.05,   # 5,000
            'AG_SL': 0.07,   # 7,000
            # Total AG = 0.12 units = 12,000 fighters

            'V': 3.0,        # 3,000 deaths/year
            'Gov': 0.35,     # 35% state presence
            'E': 18.0,       # $1.8B GDP
            'T': 0.20,       # 20% trust
            'RL': 60.0,
            'SL': 90.0,
        },

        # Stabilizing: Post-conflict recovery target
        # 3,000 fighters, minimal displacement
        'Stabilizing': {
            'SR_RL': 55.0,   # 5.5M
            'SR_SL': 65.0,   # 6.5M
            'WR_RL': 15.0,   # 1.5M
            'WR_SL': 20.0,   # 2M
            'AG_RL': 0.02,   # 2,000
            'AG_SL': 0.02,   # 2,000
            # Total AG = 0.04 units = 4,000 fighters

            'V': 1.0,        # 1,000 deaths/year
            'Gov': 0.55,     # 55% state presence
            'E': 25.0,       # $2.5B GDP
            'T': 0.40,       # 40% trust
            'RL': 70.0,
            'SL': 80.0,
        },
    }

# ==================== PEACE STRATEGIES ====================

@dataclass
class PeaceStrategy:
    name: str
    description: str

    # Budget allocation weights (must sum to 1.0)
    governance_weight: float = 0.25
    economic_weight: float = 0.25
    trust_weight: float = 0.25
    ddr_weight: float = 0.25

    # Coupling efficiency parameters
    alpha: float = 0.4   # Trust -> violence reduction efficiency
    beta: float = 0.5    # Gov -> DDR efficiency
    gamma: float = 0.5   # Econ -> DDR efficiency
    delta: float = 0.6   # Econ -> Gov coupling

PEACE_STRATEGIES = {
    'governance_only': PeaceStrategy(
        'Governance Only',
        'State-building focus: institutions, security forces',
        governance_weight=0.70, economic_weight=0.15,
        trust_weight=0.05, ddr_weight=0.10,
        alpha=0.2, beta=0.3, gamma=0.2, delta=0.3
    ),

    'economic_only': PeaceStrategy(
        'Economic Only',
        'Development focus: livelihoods, infrastructure',
        governance_weight=0.10, economic_weight=0.70,
        trust_weight=0.10, ddr_weight=0.10,
        alpha=0.2, beta=0.2, gamma=0.5, delta=0.7
    ),

    'security_only': PeaceStrategy(
        'Security Only',
        'Military focus: counter-insurgency, DDR',
        governance_weight=0.15, economic_weight=0.10,
        trust_weight=0.05, ddr_weight=0.70,
        alpha=0.1, beta=0.2, gamma=0.1, delta=0.2
    ),

    'balanced': PeaceStrategy(
        'Balanced',
        'Multi-pronged: equal investment across all dimensions',
        governance_weight=0.25, economic_weight=0.30,
        trust_weight=0.20, ddr_weight=0.25,
        alpha=0.5, beta=0.6, gamma=0.6, delta=0.7
    ),

    'integrated': PeaceStrategy(
        'Integrated',
        'Coordinated approach with strong coupling',
        governance_weight=0.25, economic_weight=0.35,
        trust_weight=0.25, ddr_weight=0.15,
        alpha=0.6, beta=0.7, gamma=0.7, delta=0.8
    ),
}

# ==================== AG STRATEGIES ====================

@dataclass
class AGStrategy:
    name: str
    recruitment_priority: float = 0.4
    violence_priority: float = 0.3
    expansion_priority: float = 0.3

AG_STRATEGIES = {
    'aggressive': AGStrategy('Aggressive', 0.3, 0.5, 0.2),
    'recruitment': AGStrategy('Recruitment Focus', 0.6, 0.2, 0.2),
    'territorial': AGStrategy('Territorial', 0.2, 0.3, 0.5),
    'defensive': AGStrategy('Defensive', 0.4, 0.4, 0.2),
}

# ==================== SOLVER ====================

def solve_calibrated_game(
    case_data: Dict[str, float],
    peace_budget: float,       # $million per year
    ag_budget: float,          # $million per year
    climate: float,            # 0-1 (climate favorability)
    peace_strategy: PeaceStrategy,
    ag_strategy: AGStrategy,
    target_ag_reduction: float  # Fraction of AG to demobilize (0-1)
) -> Dict[str, Any]:
    """
    Solve the calibrated two-budget conflict game.

    Parameters:
    - peace_budget: Government/international investment ($M/year)
    - ag_budget: Armed group resources ($M/year)
    - climate: Environmental favorability (affects all costs)
    - target_ag_reduction: What fraction of current AG should be reduced
    """

    prob = pulp.LpProblem("Calibrated_Conflict", pulp.LpMinimize)

    n_reactions = S.shape[1]

    # Decision variables
    flux = [pulp.LpVariable(f"v_{j}", lowBound=0) for j in range(n_reactions)]
    target = {sp: pulp.LpVariable(f"x_{sp}", lowBound=0) for sp in SPECIES}
    delta = {sp: pulp.LpVariable(f"d_{sp}", lowBound=None) for sp in SPECIES}

    # ===== STOICHIOMETRIC CONSTRAINTS =====
    for sp in SPECIES:
        i = SP_IDX[sp]
        reaction_effect = pulp.lpSum([S[i, j] * flux[j] for j in range(n_reactions)])
        prob += delta[sp] == reaction_effect
        initial = case_data.get(sp, 0)
        prob += target[sp] == initial + delta[sp]

    # ===== INITIAL STATE =====
    initial_AG = case_data.get('AG_RL', 0) + case_data.get('AG_SL', 0)
    initial_V = case_data.get('V', 1.0)
    initial_Gov = case_data.get('Gov', 0.1)
    initial_E = case_data.get('E', 10.0)
    initial_T = case_data.get('T', 0.1)
    initial_WR = case_data.get('WR_RL', 0) + case_data.get('WR_SL', 0)

    # Target AG level
    target_ag = initial_AG * (1 - target_ag_reduction)

    # ===== CONDITION FACTORS =====
    violence_factor = initial_V / (initial_V + 5.0)  # Saturates around V=5
    scarcity_factor = min(1.0, initial_WR / (POP_TOTAL_UNITS * 0.5))
    grievance_factor = min(1.0, (1 - initial_Gov) * (1 - initial_E/30) * scarcity_factor)

    # ===== CLIMATE COST MULTIPLIER =====
    # Bad climate makes everything more expensive
    climate_multiplier = 1.0 / max(0.2, climate)

    # ===== PEACE BUDGET ALLOCATION =====

    # Categorize reactions
    gov_reactions = ['r9', 'r10', 'r21', 'r22']
    econ_reactions = ['r7', 'r8', 'r16', 'r17']
    trust_reactions = ['r20', 'r21']
    ddr_reactions = ['r19']

    # Compute allocated budgets
    gov_budget = peace_budget * peace_strategy.governance_weight
    econ_budget = peace_budget * peace_strategy.economic_weight
    trust_budget = peace_budget * peace_strategy.trust_weight
    ddr_budget = peace_budget * peace_strategy.ddr_weight

    # Budget constraints by category
    gov_cost = pulp.lpSum([
        REACTION_COSTS.get(REACTIONS[j], 10) * climate_multiplier * flux[j]
        for j in range(n_reactions) if REACTIONS[j] in gov_reactions
    ])
    prob += gov_cost <= gov_budget

    econ_cost = pulp.lpSum([
        REACTION_COSTS.get(REACTIONS[j], 10) * climate_multiplier * flux[j]
        for j in range(n_reactions) if REACTIONS[j] in econ_reactions
    ])
    prob += econ_cost <= econ_budget

    trust_cost = pulp.lpSum([
        REACTION_COSTS.get(REACTIONS[j], 10) * climate_multiplier * flux[j]
        for j in range(n_reactions) if REACTIONS[j] in trust_reactions
    ])
    prob += trust_cost <= trust_budget

    ddr_cost = pulp.lpSum([
        REACTION_COSTS.get(REACTIONS[j], 10) * climate_multiplier * flux[j]
        for j in range(n_reactions) if REACTIONS[j] in ddr_reactions
    ])
    prob += ddr_cost <= ddr_budget

    # ===== AG BUDGET CONSTRAINT =====
    ag_reactions = ['r25', 'r26', 'r27', 'r28', 'r30', 'r31', 'r32', 'r33', 'r34']

    ag_cost = pulp.lpSum([
        REACTION_COSTS.get(REACTIONS[j], 10) * climate_multiplier * flux[j]
        for j in range(n_reactions) if REACTIONS[j] in ag_reactions
    ])
    prob += ag_cost <= ag_budget

    # AG allocation by priority
    recruitment_budget = ag_budget * ag_strategy.recruitment_priority
    violence_budget = ag_budget * ag_strategy.violence_priority

    # Recruitment constraint
    if 'r30' in RXN_IDX:
        recruitment_cost = REACTION_COSTS['r30'] * climate_multiplier * flux[RXN_IDX['r30']]
        prob += recruitment_cost <= recruitment_budget

        # Grievances boost recruitment effectiveness
        recruitment_boost = 1.0 + grievance_factor
        prob += flux[RXN_IDX['r30']] <= 0.5 * recruitment_boost  # Max recruitment rate

    # Violence constraint
    for rxn in ['r25', 'r26', 'r27', 'r28']:
        if rxn in RXN_IDX:
            prob += flux[RXN_IDX[rxn]] <= violence_budget / (REACTION_COSTS[rxn] * climate_multiplier)

    # ===== COUPLING CONSTRAINTS =====

    # Get activity levels
    gov_activity = pulp.lpSum([flux[RXN_IDX[r]] for r in gov_reactions if r in RXN_IDX])
    econ_activity = pulp.lpSum([flux[RXN_IDX[r]] for r in econ_reactions if r in RXN_IDX])
    trust_activity = pulp.lpSum([flux[RXN_IDX[r]] for r in trust_reactions if r in RXN_IDX])
    ddr_activity = pulp.lpSum([flux[RXN_IDX[r]] for r in ddr_reactions if r in RXN_IDX])

    # C1: DDR requires governance capacity
    # Can't demobilize without state to manage program
    prob += ddr_activity <= peace_strategy.beta * (initial_Gov + gov_activity * 0.1)

    # C2: DDR requires economic alternatives
    # Ex-combatants need jobs
    prob += ddr_activity <= peace_strategy.gamma * (initial_E + econ_activity) / 10.0

    # C3: Governance building needs economic foundation
    prob += gov_activity <= peace_strategy.delta * (initial_E + econ_activity) / 5.0

    # C4: Trust building needs governance presence
    prob += trust_activity <= 2.0 * (initial_Gov + gov_activity * 0.1)

    # C5: Violence destroys trust-building
    violence_penalty = max(0.1, 1.0 - 0.15 * initial_V)
    prob += trust_activity <= trust_activity * violence_penalty

    # C6: Violence destroys economic activity
    econ_penalty = max(0.2, 1.0 - 0.1 * initial_V)
    prob += econ_activity <= (econ_budget / 20.0) * econ_penalty

    # ===== MINIMUM REQUIREMENTS FOR AG REDUCTION =====

    ag_to_demobilize = initial_AG * target_ag_reduction

    if ag_to_demobilize > 0.01:
        # DDR cost: $5,000 per person = $50M per 0.1 units (10,000 people)
        # To demobilize ag_to_demobilize units, need flux of that amount through r19
        required_ddr_flux = ag_to_demobilize

        # Must have minimum governance to run DDR
        min_gov_for_ddr = 0.2 * ag_to_demobilize
        prob += initial_Gov + gov_activity * 0.1 >= min_gov_for_ddr

        # Must have minimum economy for reintegration
        min_econ_for_ddr = ag_to_demobilize * 5  # $500M econ per 0.1 AG
        prob += initial_E + econ_activity >= min_econ_for_ddr

        # Must have minimum trust for community acceptance
        min_trust_for_ddr = 0.1 * ag_to_demobilize
        prob += initial_T + trust_activity * 0.1 >= min_trust_for_ddr

    # ===== POPULATION CONSERVATION =====
    pop_species = ['SR_RL', 'WR_RL', 'AG_RL', 'SR_SL', 'WR_SL', 'AG_SL']
    prob += pulp.lpSum([target[sp] for sp in pop_species]) == POP_TOTAL_UNITS

    # ===== AG TARGET =====
    prob += target['AG_RL'] + target['AG_SL'] <= target_ag

    # ===== MINIMUM GOVERNANCE =====
    prob += target['Gov'] >= 0.05

    # ===== OBJECTIVE: MINIMIZE CONFLICT COST =====
    # Weighted combination of remaining conflict indicators
    conflict_score = (
        100 * (target['AG_RL'] + target['AG_SL']) +  # AG very bad
        20 * target['V'] +                            # Violence bad
        5 * (target['WR_RL'] + target['WR_SL']) +    # Vulnerability bad
        -10 * target['Gov'] +                         # Governance good
        -5 * target['T'] +                            # Trust good
        -2 * target['E']                              # Economy good
    )
    prob += conflict_score

    # ===== SOLVE =====
    solver = pulp.PULP_CBC_CMD(msg=False, timeLimit=30)
    status = prob.solve(solver)

    # ===== RESULTS =====
    results = {
        'status': pulp.LpStatus[status],
        'optimal': status == pulp.LpStatusOptimal,
        'peace_strategy': peace_strategy.name,
        'ag_strategy': ag_strategy.name,
        'violence_factor': violence_factor,
        'grievance_factor': grievance_factor,
    }

    if results['optimal']:
        results['objective'] = pulp.value(prob.objective)

        ag_final = pulp.value(target['AG_RL']) + pulp.value(target['AG_SL'])
        results['ag_initial'] = initial_AG
        results['ag_final'] = ag_final
        results['ag_initial_people'] = int(initial_AG * POP_UNIT)
        results['ag_final_people'] = int(ag_final * POP_UNIT)
        results['ag_reduction_people'] = int((initial_AG - ag_final) * POP_UNIT)
        results['ag_reduction_pct'] = (initial_AG - ag_final) / initial_AG * 100 if initial_AG > 0 else 0

        # Budget utilization
        results['gov_budget_used'] = pulp.value(gov_cost) if gov_cost else 0
        results['econ_budget_used'] = pulp.value(econ_cost) if econ_cost else 0
        results['trust_budget_used'] = pulp.value(trust_cost) if trust_cost else 0
        results['ddr_budget_used'] = pulp.value(ddr_cost) if ddr_cost else 0

        # Cost per demobilized (if any demobilization)
        if results['ag_reduction_people'] > 0:
            total_spent = (results['gov_budget_used'] + results['econ_budget_used'] +
                          results['trust_budget_used'] + results['ddr_budget_used'])
            results['cost_per_demobilized'] = total_spent * 1_000_000 / results['ag_reduction_people']
        else:
            results['cost_per_demobilized'] = 0

    return results

# ==================== ANALYSIS ====================

def run_calibrated_analysis():
    """Run analysis with calibrated parameters"""

    print("\n" + "="*70)
    print("CALIBRATED CONFLICT MODEL - Lake Chad Basin")
    print("="*70)
    print("\nUnit interpretations:")
    print(f"  Population unit: {POP_UNIT:,} people")
    print(f"  Budget unit: ${BUDGET_UNIT:,}")
    print(f"  Total population: {POP_TOTAL_UNITS} units = {POP_TOTAL_UNITS * POP_UNIT / 1e6:.1f}M people")

    cases = define_calibrated_cases()

    # Print case summaries
    print("\nConflict scenarios:")
    for name, data in cases.items():
        ag = data['AG_RL'] + data['AG_SL']
        print(f"  {name}: AG={ag:.2f} units ({int(ag * POP_UNIT):,} fighters), V={data['V']:.0f}, Gov={data['Gov']:.0%}")

    # Realistic budget levels ($M/year)
    peace_budgets = [500, 1000, 2000, 3000]  # $500M to $3B
    ag_budgets = [50, 100, 200, 500]          # $50M to $500M
    climates = [0.3, 0.5, 0.7, 0.9]
    ag_reduction_targets = [0.25, 0.50, 0.75]  # 25%, 50%, 75% reduction

    print(f"\nPeace budgets: ${peace_budgets} million/year")
    print(f"AG budgets: ${ag_budgets} million/year")
    print(f"AG reduction targets: {[f'{t:.0%}' for t in ag_reduction_targets]}")

    all_results = []

    total = (len(cases) * len(peace_budgets) * len(ag_budgets) *
             len(climates) * len(PEACE_STRATEGIES) * len(AG_STRATEGIES) *
             len(ag_reduction_targets))

    print(f"\nTotal scenarios: {total}")
    print("Running...")

    count = 0
    for case_name, case_data in cases.items():
        for peace_budget in peace_budgets:
            for ag_budget in ag_budgets:
                for climate in climates:
                    for ps_name, peace_strat in PEACE_STRATEGIES.items():
                        for ag_name, ag_strat in AG_STRATEGIES.items():
                            for target_red in ag_reduction_targets:
                                count += 1
                                if count % 1000 == 0:
                                    print(f"  {count}/{total}")

                                try:
                                    result = solve_calibrated_game(
                                        case_data, peace_budget, ag_budget,
                                        climate, peace_strat, ag_strat, target_red
                                    )

                                    record = {
                                        'case': case_name,
                                        'peace_budget': peace_budget,
                                        'ag_budget': ag_budget,
                                        'budget_ratio': peace_budget / ag_budget,
                                        'climate': climate,
                                        'peace_strategy': ps_name,
                                        'ag_strategy': ag_name,
                                        'target_reduction': target_red,
                                        'feasible': result['optimal'],
                                        'strategy_type': 'mono' if 'only' in ps_name.lower() else 'multi',
                                        **{k: v for k, v in result.items()
                                           if k not in ['optimal', 'status']}
                                    }
                                    all_results.append(record)

                                except Exception as e:
                                    all_results.append({
                                        'case': case_name,
                                        'peace_budget': peace_budget,
                                        'ag_budget': ag_budget,
                                        'climate': climate,
                                        'peace_strategy': ps_name,
                                        'ag_strategy': ag_name,
                                        'target_reduction': target_red,
                                        'feasible': False,
                                        'strategy_type': 'mono' if 'only' in ps_name.lower() else 'multi',
                                        'error': str(e)[:50]
                                    })

    df = pd.DataFrame(all_results)
    return df

def print_calibrated_results(df: pd.DataFrame):
    """Print results with real-world interpretations"""

    print("\n" + "="*70)
    print("RESULTS SUMMARY (Real-World Interpretation)")
    print("="*70)

    print(f"\nTotal scenarios: {len(df)}")
    print(f"Feasible: {df['feasible'].sum()} ({df['feasible'].mean()*100:.1f}%)")

    print("\n1. FEASIBILITY BY STRATEGY:")
    print("-"*50)
    for strat in sorted(df['peace_strategy'].unique()):
        data = df[df['peace_strategy'] == strat]
        rate = data['feasible'].mean() * 100
        stype = 'MONO ' if 'only' in strat.lower() else 'MULTI'
        print(f"  [{stype}] {strat:20s}: {rate:5.1f}%")

    print("\n2. FEASIBILITY BY CONFLICT SEVERITY:")
    print("-"*50)
    for case in df['case'].unique():
        mono = df[(df['case'] == case) & (df['strategy_type'] == 'mono')]
        multi = df[(df['case'] == case) & (df['strategy_type'] == 'multi')]
        print(f"  {case:12s}: mono={mono['feasible'].mean()*100:5.1f}%, multi={multi['feasible'].mean()*100:5.1f}%")

    print("\n3. FEASIBILITY BY BUDGET RATIO (Peace/AG):")
    print("-"*50)
    for ratio in sorted(df['budget_ratio'].unique())[:6]:
        mono = df[(df['budget_ratio'] == ratio) & (df['strategy_type'] == 'mono')]
        multi = df[(df['budget_ratio'] == ratio) & (df['strategy_type'] == 'multi')]
        if len(mono) > 0 and len(multi) > 0:
            print(f"  ratio={ratio:5.1f}x: mono={mono['feasible'].mean()*100:5.1f}%, multi={multi['feasible'].mean()*100:5.1f}%")

    print("\n4. FEASIBILITY BY REDUCTION TARGET:")
    print("-"*50)
    for target in sorted(df['target_reduction'].unique()):
        mono = df[(df['target_reduction'] == target) & (df['strategy_type'] == 'mono')]
        multi = df[(df['target_reduction'] == target) & (df['strategy_type'] == 'multi')]
        print(f"  {target:.0%} reduction: mono={mono['feasible'].mean()*100:5.1f}%, multi={multi['feasible'].mean()*100:5.1f}%")

    feasible = df[df['feasible']].copy()
    if len(feasible) > 0 and 'cost_per_demobilized' in feasible.columns:
        print("\n5. COST PER DEMOBILIZED (feasible only):")
        print("-"*50)
        valid_costs = feasible[feasible['cost_per_demobilized'] > 0]['cost_per_demobilized']
        if len(valid_costs) > 0:
            print(f"  Average: ${valid_costs.mean():,.0f} per ex-combatant")
            print(f"  Range: ${valid_costs.min():,.0f} - ${valid_costs.max():,.0f}")
            print(f"  (Literature estimate: $3,000-$10,000 per person)")

def create_calibrated_visualizations(df: pd.DataFrame):
    """Create visualizations with real-world labels"""

    fig = plt.figure(figsize=(18, 14))

    # 1. Feasibility by strategy
    ax1 = plt.subplot(2, 3, 1)
    strat_feas = df.groupby('peace_strategy')['feasible'].mean() * 100
    colors = ['#e74c3c' if 'only' in s.lower() else '#27ae60' for s in strat_feas.index]
    ax1.barh(range(len(strat_feas)), strat_feas.values, color=colors)
    ax1.set_yticks(range(len(strat_feas)))
    ax1.set_yticklabels(strat_feas.index, fontsize=9)
    ax1.set_xlabel('Peace Achievability (%)')
    ax1.set_title('Strategy Effectiveness\n(Red=Single-Focus, Green=Multi-Pronged)', fontweight='bold')
    ax1.axvline(50, color='gray', linestyle='--', alpha=0.5)

    # 2. Budget ratio impact
    ax2 = plt.subplot(2, 3, 2)
    for stype, color in [('mono', '#e74c3c'), ('multi', '#27ae60')]:
        data = df[df['strategy_type'] == stype]
        ratio_feas = data.groupby('budget_ratio')['feasible'].mean() * 100
        ax2.plot(ratio_feas.index, ratio_feas.values, 'o-',
                color=color, linewidth=2, markersize=8, label=stype.title())
    ax2.set_xlabel('Peace Budget / AG Budget')
    ax2.set_ylabel('Peace Achievability (%)')
    ax2.set_title('Budget Ratio Impact\n(How much must peace outspend violence?)', fontweight='bold')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_xscale('log')

    # 3. Reduction target impact
    ax3 = plt.subplot(2, 3, 3)
    for stype, color in [('mono', '#e74c3c'), ('multi', '#27ae60')]:
        data = df[df['strategy_type'] == stype]
        target_feas = data.groupby('target_reduction')['feasible'].mean() * 100
        ax3.plot([f'{t:.0%}' for t in target_feas.index], target_feas.values, 'o-',
                color=color, linewidth=2, markersize=10, label=stype.title())
    ax3.set_xlabel('AG Reduction Target')
    ax3.set_ylabel('Peace Achievability (%)')
    ax3.set_title('Ambition vs Success\n(More reduction = harder to achieve)', fontweight='bold')
    ax3.legend()
    ax3.grid(True, alpha=0.3)

    # 4. Conflict severity heatmap
    ax4 = plt.subplot(2, 3, 4)
    pivot = df.pivot_table(values='feasible', index='case',
                          columns='peace_strategy', aggfunc='mean') * 100
    sns.heatmap(pivot, annot=True, fmt='.0f', cmap='RdYlGn', ax=ax4,
                vmin=0, vmax=100, cbar_kws={'label': 'Feasibility %'})
    ax4.set_title('Severity x Strategy', fontweight='bold')

    # 5. Cost effectiveness (if available)
    ax5 = plt.subplot(2, 3, 5)
    feasible = df[df['feasible']].copy()
    if 'cost_per_demobilized' in feasible.columns:
        valid = feasible[feasible['cost_per_demobilized'] > 0]
        if len(valid) > 0:
            for stype, color in [('mono', '#e74c3c'), ('multi', '#27ae60')]:
                data = valid[valid['strategy_type'] == stype]
                if len(data) > 0:
                    costs = data.groupby('peace_budget')['cost_per_demobilized'].mean() / 1000
                    ax5.plot(costs.index, costs.values, 'o-', color=color,
                            linewidth=2, markersize=8, label=stype.title())
            ax5.axhline(5, color='blue', linestyle='--', alpha=0.5, label='Literature estimate ($5k)')
            ax5.set_xlabel('Peace Budget ($M/year)')
            ax5.set_ylabel('Cost per Demobilized ($thousands)')
            ax5.set_title('DDR Cost Efficiency', fontweight='bold')
            ax5.legend()
            ax5.grid(True, alpha=0.3)

    # 6. Key findings
    ax6 = plt.subplot(2, 3, 6)
    ax6.axis('off')

    mono_feas = df[df['strategy_type'] == 'mono']['feasible'].mean() * 100
    multi_feas = df[df['strategy_type'] == 'multi']['feasible'].mean() * 100

    severe = df[df['case'] == 'Severe']
    severe_mono = severe[severe['strategy_type'] == 'mono']['feasible'].mean() * 100
    severe_multi = severe[severe['strategy_type'] == 'multi']['feasible'].mean() * 100

    summary = [
        "CALIBRATED MODEL FINDINGS",
        "-" * 40,
        f"Population: {POP_TOTAL_UNITS * POP_UNIT / 1e6:.0f}M people",
        f"AG range: 4,000-25,000 fighters",
        "",
        f"Overall feasibility:",
        f"  Single-focus: {mono_feas:.1f}%",
        f"  Multi-pronged: {multi_feas:.1f}%",
        "",
        f"Severe conflict (25k fighters):",
        f"  Single-focus: {severe_mono:.1f}%",
        f"  Multi-pronged: {severe_multi:.1f}%",
        "",
        "Key insight:",
        "  Multi-pronged approaches required",
        "  for sustained peace, especially",
        "  in severe conflict scenarios.",
    ]

    ax6.text(0.05, 0.95, '\n'.join(summary), fontsize=10, va='top', family='monospace',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    plt.tight_layout()
    plt.savefig('calibrated_conflict_results.png', dpi=150, bbox_inches='tight')
    print("\nVisualization saved to 'calibrated_conflict_results.png'")
    plt.show()

# ==================== MAIN ====================

if __name__ == "__main__":

    print("\n" + "="*70)
    print("LAKE CHAD BASIN CONFLICT MODEL")
    print("Calibrated with Real-World Data")
    print("="*70)

    # Run analysis
    results_df = run_calibrated_analysis()

    # Save results
    results_df.to_csv('calibrated_conflict_results.csv', index=False)
    print(f"\nResults saved to 'calibrated_conflict_results.csv'")

    # Print summary
    print_calibrated_results(results_df)

    # Visualize
    create_calibrated_visualizations(results_df)

    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)
