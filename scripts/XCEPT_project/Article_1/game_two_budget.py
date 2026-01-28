"""
CONFLICT RESOLUTION GAME
A gamified version of the two-budget conflict model where you (Government)
try to achieve peace with minimal budget against an adaptive Armed Group.
"""

import pulp
import numpy as np
import sys
import os
from typing import Dict, List, Any, Tuple
from dataclasses import dataclass
import json
from datetime import datetime

# Try to load the reaction network from file
try:
    sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))
    from pyCOT.io.functions import read_txt
    
    # Load reaction network
    project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
    FILE_PATH = os.path.join(project_root, 'networks/Conflict_Theory/Resource_Community_Insurgency_Loops_model1.txt')
    rn = read_txt(FILE_PATH)
    S = rn.stoichiometry_matrix()
    SPECIES = list(rn._species_map.keys())
    REACTIONS = [r.name() for r in rn.reactions()]
    
    print(f"[Game] Loaded reaction network with {len(SPECIES)} species and {len(REACTIONS)} reactions")
    
except ImportError as e:
    print(f"[Game] Warning: Could not load pyCOT: {e}")
    print("[Game] Using fallback stoichiometry")
    # Fallback to simplified network
    REACTIONS = [
        'r1', 'r2', 'r3', 'r4', 'r5', 'r6', 'r7', 'r8', 'r9', 'r10',
        'r11', 'r12', 'r13', 'r14', 'r15', 'r16', 'r17', 'r18', 'r19', 'r20',
        'r21', 'r22', 'r23', 'r24', 'r25', 'r26', 'r27', 'r28', 'r29', 'r30',
        'r31', 'r32', 'r33', 'r34', 'r35', 'r36'
    ]
    SPECIES = [
        'SR_RL', 'SR_SL', 'WR_RL', 'WR_SL', 'AG_RL', 'AG_SL',
        'RL', 'SL', 'E', 'T', 'V', 'Gov'
    ]
    # Create placeholder stoichiometry matrix
    S = np.zeros((len(SPECIES), len(REACTIONS)))
    
    # Set up some basic stoichiometry for gameplay
    # This is simplified - in reality you'd parse the reaction network
    reaction_definitions = {
        'r1': {'SL': -1, 'RL': 1},
        'r2': {'SR_SL': -1, 'SR_RL': 1},
        'r8': {'Gov': -1, 'SR_SL': -1, 'SR_RL': 1},
        'r9': {'SR_RL': -1, 'SR_SL': 1, 'E': 1},
        'r11': {'E': -2, 'Gov': 1},
        'r12': {'Gov': 1},
        'r21': {'E': -1, 'Gov': -1, 'AG_SL': -1, 'WR_SL': 1},
        'r29': {'AG_SL': -1, 'Gov': -1, 'WR_SL': 1, 'V': 1},
    }
    
    for rxn_name, stoich in reaction_definitions.items():
        if rxn_name in REACTIONS:
            rxn_idx = REACTIONS.index(rxn_name)
            for sp, coeff in stoich.items():
                if sp in SPECIES:
                    sp_idx = SPECIES.index(sp)
                    S[sp_idx, rxn_idx] = coeff

# ==================== GAME PARAMETERS ====================

@dataclass
class GameConfig:
    """Configurable game parameters"""
    # Budget settings
    INITIAL_PEACE_BUDGET: float = 100.0
    INITIAL_AG_BUDGET: float = 50.0
    BUDGET_INCREMENT_PER_ROUND: float = 5.0  # Both sides get more budget each round
    
    # Climate settings
    CLIMATE_MIN: float = 0.3  # Harsh
    CLIMATE_MAX: float = 0.7  # Favorable
    DEFAULT_CLIMATE: float = 0.5
    
    # Victory/Defeat thresholds
    PEACE_VICTORY_THRESHOLD: float = 10  # Peace index to win
    PEACE_DEFEAT_THRESHOLD: float = 1   # Peace index to lose
    MAX_ROUNDS: int = 10
    
    # Scoring weights
    SCORE_BUDGET_WEIGHT: float = 0.4
    SCORE_CLIMATE_WEIGHT: float = 0.3
    SCORE_ROUNDS_WEIGHT: float = 0.3
    
    # Game difficulty
    AG_STRATEGY: str = "adaptive"  # Options: "adaptive", "aggressive", "defensive"
    
    def __post_init__(self):
        self.save_file = "game_save.json"

# ==================== REACTION CLASSIFICATION ====================

# Policy reactions from the network
GOV_POLICY_RXNS = ['r8', 'r11', 'r12', 'r21', 'r23', 'r24', 'r29']
AG_POLICY_RXNS = ['r17', 'r28', 'r30', 'r33', 'r34', 'r35', 'r36']
FIXED_RXNS = [r for r in REACTIONS if r not in GOV_POLICY_RXNS + AG_POLICY_RXNS]

RXN_IDX = {name: idx for idx, name in enumerate(REACTIONS)}
SP_IDX = {name: idx for idx, name in enumerate(SPECIES)}

# Define strategies dictionary as a global variable
STRATEGIES = {
    'Multi-Pronged': GOV_POLICY_RXNS,
    'Security-Only': ['r29'],
    'DDR-Only': ['r21'],
    'Economic-Only': ['r11', 'r24'],
    'Governance-Only': ['r12', 'r23'],
    'Land-Only': ['r8'],
}

# ==================== REACTION SEMANTICS ====================

# Human-readable names for each reaction
REACTION_NAMES = {
    'r1': 'Land Recovery (natural)',
    'r2': 'SR move to restored land',
    'r3': 'WR move to restored land',
    'r4': 'Land Degradation',
    'r5': 'SR crowding on scarce land',
    'r6': 'WR crowding on scarce land',
    'r7': 'AG crowding on scarce land',
    'r8': 'Gov-assisted land restoration',
    'r9': 'SR economic production',
    'r10': 'WR economic production',
    'r11': 'Economy -> Governance',
    'r12': 'Base governance presence',
    'r13': 'SR migration (scarcity-driven)',
    'r14': 'WR migration (scarcity-driven)',
    'r15': 'SR displacement (violence)',
    'r16': 'WR displacement (violence)',
    'r17': 'AG movement to resources',
    'r18': 'WR uplift to SR (economy)',
    'r19': 'WR uplift to SR (economy)',
    'r20': 'SR stress -> WR',
    'r21': 'DDR: AG -> WR (E+Gov required)',
    'r22': 'Trust generation',
    'r23': 'Gov trust-building',
    'r24': 'Gov direct uplift',
    'r25': 'Violence destroys trust',
    'r26': 'Governance decay',
    'r27': 'Grievance -> Violence',
    'r28': 'AG predation on SR',
    'r29': 'Counter-insurgency (Gov)',
    'r30': 'Counter-insurgency (AG wins)',
    'r31': 'Violence decay',
    'r32': 'AG recruitment from WR',
    'r33': 'AG extraction from WR',
    'r34': 'AG extraction from SR',
    'r35': 'AG economic predation',
    'r36': 'AG expels government'
}

# Classify reactions by their role in peace dynamics
PEACE_LOOPS = {
    'Land_Restoration': {
        'reactions': ['r1', 'r2', 'r3', 'r8'],
        'description': 'Natural and policy-driven land recovery',
        'icon': '🌱'
    },
    'Economic_Development': {
        'reactions': ['r9', 'r10', 'r11'],
        'description': 'Economic production and Gov building',
        'icon': '💰'
    },
    'Governance_Building': {
        'reactions': ['r12', 'r23', 'r24'],
        'description': 'Institutional presence and trust-building',
        'icon': '🏛️'
    },
    'DDR_Reintegration': {
        'reactions': ['r21'],
        'description': 'Demobilization: AG -> WR (requires E, Gov)',
        'icon': '🕊️'
    },
    'Counter_Insurgency': {
        'reactions': ['r29'],
        'description': 'State-led AG reduction (creates violence)',
        'icon': '⚔️'
    },
    'Resilience_Uplift': {
        'reactions': ['r18', 'r19', 'r22'],
        'description': 'WR -> SR transitions and trust generation',
        'icon': '💪'
    }
}

# Classify reactions by their role in conflict dynamics
CONFLICT_LOOPS = {
    'Recruitment': {
        'reactions': ['r32'],
        'description': 'WR + AG + V -> 2AG (grievance-driven)',
        'icon': '📢'
    },
    'Predation': {
        'reactions': ['r28', 'r33', 'r34', 'r35'],
        'description': 'AG extraction from populations and economy',
        'icon': '💀'
    },
    'Violence_Generation': {
        'reactions': ['r27', 'r30'],
        'description': 'Grievance -> Violence, AG attacks',
        'icon': '🔥'
    },
    'State_Erosion': {
        'reactions': ['r36', 'r26'],
        'description': 'AG expels Gov, natural Gov decay',
        'icon': '🏚️'
    },
    'Displacement': {
        'reactions': ['r15', 'r16'],
        'description': 'Violence-driven population displacement',
        'icon': '🏃'
    },
    'Land_Degradation': {
        'reactions': ['r4', 'r5', 'r6', 'r7'],
        'description': 'Resource depletion and crowding',
        'icon': '🏜️'
    }
}

# ==================== MODEL PARAMETERS ====================

@dataclass
class ModelParams:
    """Game model parameters"""
    # Base rates
    k_land_recovery: float = 0.3
    k_land_degradation: float = 0.2
    base_productivity: float = 0.8
    wr_productivity_factor: float = 0.5
    base_migration_rate: float = 0.15
    base_displacement_rate: float = 0.25
    base_uplift_rate: float = 0.2
    base_stress_decay: float = 0.1
    base_trust_gen: float = 0.4
    base_trust_destruction: float = 0.3
    base_violence_decay: float = 0.15
    base_gov_decay: float = 0.05
    base_recruitment: float = 0.15
    base_extraction: float = 0.2
    
    # Thresholds
    scarcity_threshold: float = 1.5
    econ_threshold: float = 1.0
    gov_threshold: float = 0.3
    gov_uplift_threshold: float = 0.4
    
    # Cost parameters
    C_gov_land: float = 8.0
    C_gov_build: float = 10.0
    C_base_gov: float = 2.0
    C_ddr: float = 12.0
    C_trust_gov: float = 6.0
    C_gov_uplift: float = 8.0
    C_counter_gov: float = 15.0
    
    C_ag_move: float = 3.0
    C_ag_predation: float = 5.0
    C_counter_ag: float = 12.0
    C_ag_extract_wr: float = 2.5
    C_ag_extract_sr: float = 3.0
    C_ag_econ_pred: float = 4.0
    C_ag_expel: float = 6.0
    
    epsilon: float = 0.01
    
    # Derived
    @property
    def V_max_econ(self):
        return 2.0 * self.base_productivity
    
    @property
    def K_m_econ(self):
        return 5.0
    
    @property
    def V_max_trust(self):
        return 2.0 * self.base_trust_gen
    
    @property
    def K_m_trust(self):
        return 2.0
    
    @property
    def V_max_grievance(self):
        return 1.5
    
    @property
    def K_m_grievance(self):
        return 10.0
    
    @property
    def V_max_recruit(self):
        return 2.0 * self.base_recruitment
    
    @property
    def K_m_recruit(self):
        return 5.0

# ==================== GAME STATE ====================

class GameState:
    """Tracks the complete game state"""
    
    def __init__(self, config: GameConfig):
        self.config = config
        self.round = 0
        self.peace_budget = config.INITIAL_PEACE_BUDGET
        self.ag_budget = config.INITIAL_AG_BUDGET
        self.climate = config.DEFAULT_CLIMATE
        self.score = 0
        self.history = []
        
        # Initial system state
        self.state = {
            'SR_RL': 10.0, 'SR_SL': 12.0,
            'WR_RL': 6.0, 'WR_SL': 11.0,
            'AG_RL': 0.08, 'AG_SL': 0.08,
            'RL': 20.0, 'SL': 20.0,
            'E': 1.0, 'T': 0.15, 'V': 2.5, 'Gov': 0.3
        }
        
        # Player choices
        self.player_strategy = "Multi-Pronged"
        self.player_budget_allocated = 0
        
    def calculate_peace_index(self) -> float:
        """Calculate peace index (0-inf) based on current state"""
        violence_trust_norm = 10*(self.state['T']+self.state['Gov']+self.state['E'])/(self.state['V']+1e-5)
        population_norm = 0.1 * (self.state['SR_RL'] + self.state['WR_RL']) / (self.state['AG_RL'] + self.state['AG_SL'] + 1e-5)
        print("peace_index calculation:")
        print(f"  violence_trust_norm: {violence_trust_norm}")
        print(f"  population_norm: {population_norm}")
        print("minimum of the two used as peace_index")
        peace_index = min(violence_trust_norm,population_norm)
        return peace_index
    
    def get_game_status(self) -> str:
        """Check if game is won, lost, or ongoing"""
        peace_index = self.calculate_peace_index()
        
        if peace_index >= self.config.PEACE_VICTORY_THRESHOLD:
            return "VICTORY"
        elif peace_index <= self.config.PEACE_DEFEAT_THRESHOLD:
            return "DEFEAT"
        elif self.round >= self.config.MAX_ROUNDS:
            return "TIMEOUT"
        else:
            return "ONGOING"
    
    def calculate_score(self) -> float:
        """Calculate final score if game is won"""
        if self.get_game_status() != "VICTORY":
            return 0
        
        # Base score
        base_score = 1000
        
        # Budget efficiency (lower is better)
        total_budget_used = sum(h['player_budget'] for h in self.history)
        max_possible_budget = self.config.INITIAL_PEACE_BUDGET + \
                             self.config.BUDGET_INCREMENT_PER_ROUND * self.round
        
        budget_efficiency = 1.0 - (total_budget_used / max_possible_budget)
        
        # Climate penalty (harsher climate gives more points)
        avg_climate = np.mean([h['climate'] for h in self.history])
        climate_factor = 1.0 + (self.config.CLIMATE_MIN / avg_climate) * 0.5
        
        # Round efficiency (fewer rounds is better, but balanced with cheapness)
        round_factor = 1.0 + (self.config.MAX_ROUNDS - self.round) / self.config.MAX_ROUNDS
        
        final_score = base_score * budget_efficiency * climate_factor * round_factor
        
        return final_score
    
    def save_game(self):
        """Save game state to file"""
        save_data = {
            'round': self.round,
            'peace_budget': self.peace_budget,
            'ag_budget': self.ag_budget,
            'climate': self.climate,
            'state': self.state,
            'player_strategy': self.player_strategy,
            'history': self.history,
            'config': self.config.__dict__,
            'timestamp': datetime.now().isoformat()
        }
        
        with open(self.config.save_file, 'w') as f:
            json.dump(save_data, f, indent=2)
        
        print(f"[Game] Game saved to {self.config.save_file}")
    
    def load_game(self):
        """Load game state from file"""
        try:
            with open(self.config.save_file, 'r') as f:
                save_data = json.load(f)
            
            self.round = save_data['round']
            self.peace_budget = save_data['peace_budget']
            self.ag_budget = save_data['ag_budget']
            self.climate = save_data['climate']
            self.state = save_data['state']
            self.player_strategy = save_data['player_strategy']
            self.history = save_data['history']
            
            print(f"[Game] Game loaded from {self.config.save_file}")
            return True
        except FileNotFoundError:
            print("[Game] No saved game found")
            return False

# ==================== GAME MECHANICS ====================

class GameEngine:
    """Handles the game mechanics and LP solving"""
    
    def __init__(self, config: GameConfig):
        self.config = config
        self.params = ModelParams()
        self.game_state = GameState(config)
        self.strategies = STRATEGIES  # Add strategies as instance attribute
    
    def calculate_fixed_rates(self, state: Dict[str, float], climate: float) -> Dict[str, float]:
        """Calculate rates for fixed reactions"""
        rates = {}
        
        # Extract state variables
        SR_RL = state['SR_RL']; SR_SL = state['SR_SL']
        WR_RL = state['WR_RL']; WR_SL = state['WR_SL']
        AG_RL = state['AG_RL']; AG_SL = state['AG_SL']
        RL = state['RL']; SL = state['SL']
        E = state['E']; T = state['T']; V = state['V']; Gov = state['Gov']
        
        total_pop = SR_RL + SR_SL + WR_RL + WR_SL + AG_RL + AG_SL
        
        # Simplified rate calculations for gameplay
        # Land dynamics
        rates['r1'] = self.params.k_land_recovery * climate * SL
        rates['r2'] = self.params.k_land_recovery * SR_SL * RL * 0.5
        rates['r3'] = self.params.k_land_recovery * WR_SL * RL * 0.5
        
        pop_pressure = total_pop / (RL + self.params.epsilon)
        rates['r4'] = self.params.k_land_degradation * min(1.0, pop_pressure/5) * RL
        
        # Economic production
        rates['r9'] = self.params.V_max_econ * SR_RL / (self.params.K_m_econ + SR_RL + self.params.epsilon)
        rates['r10'] = self.params.wr_productivity_factor * self.params.V_max_econ * WR_RL / \
                      (self.params.K_m_econ + WR_RL + self.params.epsilon)
        
        # Migration
        scarcity_ratio = SL / (RL + self.params.epsilon)
        if scarcity_ratio > self.params.scarcity_threshold:
            rates['r13'] = self.params.base_migration_rate * SR_SL * RL
            rates['r14'] = self.params.base_migration_rate * WR_SL * RL
        else:
            rates['r13'] = rates['r14'] = 0.0
        
        # Violence effects
        violence_factor = V / (V + 1.0)
        rates['r15'] = self.params.base_displacement_rate * violence_factor * SR_RL * SL
        rates['r16'] = self.params.base_displacement_rate * violence_factor * WR_RL * SL
        
        # Resilience
        if E > self.params.econ_threshold:
            rates['r18'] = self.params.base_uplift_rate * E * WR_RL
            rates['r19'] = self.params.base_uplift_rate * E * WR_SL
        else:
            rates['r18'] = rates['r19'] = 0.0
        
        stress = 0.3 * (1 - E) + 0.4 * violence_factor + 0.3 * (SL/(RL+SL+self.params.epsilon))
        rates['r20'] = self.params.base_stress_decay * stress * SR_SL
        
        # Trust & violence
        rates['r22'] = self.params.V_max_trust * SR_RL * E / (self.params.K_m_trust + E + self.params.epsilon)
        rates['r25'] = self.params.base_trust_destruction * V * T
        rates['r26'] = self.params.base_gov_decay * Gov * Gov
        rates['r27'] = self.params.V_max_grievance * WR_SL * WR_SL / \
                      (self.params.K_m_grievance + WR_SL * WR_SL + self.params.epsilon)
        rates['r31'] = self.params.base_violence_decay * V * V
        rates['r32'] = self.params.V_max_recruit * AG_SL * WR_SL * V / \
                      (self.params.K_m_recruit + WR_SL + self.params.epsilon)
        
        return rates
    
    def calculate_policy_costs(self, state: Dict[str, float]) -> Dict[str, float]:
        """Calculate cost per unit for policy reactions"""
        costs = {}
        eps = self.params.epsilon
        
        E = state['E']; Gov = state['Gov']; V = state['V']; T = state['T']
        RL = state['RL']; SR_RL = state['SR_RL']; WR_RL = state['WR_RL']
        AG_RL = state['AG_RL']; AG_SL = state['AG_SL']
        
        # Government costs
        costs['r8'] = self.params.C_gov_land / (Gov + eps)
        costs['r11'] = self.params.C_gov_build / (E + eps)
        costs['r12'] = self.params.C_base_gov
        costs['r21'] = self.params.C_ddr/(E + eps) + self.params.C_ddr/(Gov + eps)
        costs['r23'] = self.params.C_trust_gov / (Gov + eps)
        costs['r24'] = self.params.C_gov_uplift / (Gov + eps)
        
        violence_factor = V / (V + 1.0)
        gov_capacity = Gov * (1 - violence_factor)
        costs['r29'] = self.params.C_counter_gov / (gov_capacity + eps)
        
        # AG costs
        costs['r17'] = self.params.C_ag_move / (RL + eps)
        costs['r28'] = self.params.C_ag_predation / (SR_RL * (1 - Gov) + eps)
        costs['r30'] = self.params.C_counter_ag / (AG_SL * (1 - Gov) + eps)
        costs['r33'] = self.params.C_ag_extract_wr / (WR_RL + eps)
        costs['r34'] = self.params.C_ag_extract_sr / (SR_RL * (1 - T) + eps)
        costs['r35'] = self.params.C_ag_econ_pred / (E + eps)
        costs['r36'] = self.params.C_ag_expel / (AG_RL * (1 - Gov) + eps)
        
        return costs
    
    def get_strategy_reactions(self, strategy: str) -> List[str]:
        """Get reactions allowed for a given strategy"""
        return self.strategies.get(strategy, GOV_POLICY_RXNS)
    
    def solve_player_turn(self, player_budget: float, strategy: str) -> Dict[str, Any]:
        """Solve LP for player's turn with relaxed constraints and fallback"""

        # Try first with normal constraints, then relax if infeasible
        for attempt, (pop_tolerance, min_gov) in enumerate([
            (0.01, 0.05),   # Attempt 1: Tight constraints
            (0.1, 0.02),    # Attempt 2: Relaxed population, lower min gov
            (0.5, 0.0),     # Attempt 3: Very relaxed
        ]):
            result = self._solve_player_lp(player_budget, strategy, pop_tolerance, min_gov)
            if result['optimal']:
                if attempt > 0:
                    result['relaxed'] = True
                    result['relaxation_level'] = attempt
                return result

        # All attempts failed - return diagnostic info
        return self._create_failure_result(strategy, player_budget)

    def _solve_player_lp(self, player_budget: float, strategy: str,
                         pop_tolerance: float, min_gov: float) -> Dict[str, Any]:
        """Internal LP solver with configurable constraint tightness"""
        prob = pulp.LpProblem("Player_Turn", pulp.LpMinimize)

        n_reactions = len(REACTIONS)
        flux_vars = [pulp.LpVariable(f"v_{i}", lowBound=0) for i in range(n_reactions)]
        target_vars = {sp: pulp.LpVariable(f"t_{sp}", lowBound=0) for sp in SPECIES}

        # Calculate rates and costs
        fixed_rates = self.calculate_fixed_rates(self.game_state.state, self.game_state.climate)
        policy_costs = self.calculate_policy_costs(self.game_state.state)

        # Stoichiometric constraints
        for i, sp in enumerate(SPECIES):
            total_change = pulp.lpSum([S[i, j] * flux_vars[j] for j in range(n_reactions)])
            prob += target_vars[sp] == self.game_state.state.get(sp, 0) + total_change

        # Fixed reaction bounds - use upper bounds only to allow flexibility
        TIME_SCALE = 0.1
        for rxn in FIXED_RXNS:
            if rxn in RXN_IDX:
                idx = RXN_IDX[rxn]
                rate = fixed_rates.get(rxn, 0) * TIME_SCALE
                if rate > 0:
                    prob += flux_vars[idx] <= rate

        # Strategy constraint
        allowed_rxns = self.get_strategy_reactions(strategy)
        for rxn in GOV_POLICY_RXNS:
            if rxn not in allowed_rxns and rxn in RXN_IDX:
                prob += flux_vars[RXN_IDX[rxn]] == 0

        # Budget constraint
        gov_cost = pulp.lpSum([
            policy_costs.get(rxn, 0) * flux_vars[RXN_IDX[rxn]]
            for rxn in allowed_rxns if rxn in RXN_IDX
        ])
        prob += gov_cost <= player_budget

        # Population conservation with TOLERANCE (not strict equality)
        pop_species = ['SR_RL', 'SR_SL', 'WR_RL', 'WR_SL', 'AG_RL', 'AG_SL']
        initial_pop = sum(self.game_state.state.get(sp, 0) for sp in pop_species)
        final_pop = pulp.lpSum([target_vars[sp] for sp in pop_species])

        # Allow population to vary within tolerance
        prob += final_pop >= initial_pop * (1 - pop_tolerance)
        prob += final_pop <= initial_pop * (1 + pop_tolerance)

        # Minimum governance (configurable)
        if min_gov > 0:
            prob += target_vars['Gov'] >= min_gov

        # Non-negative species (already set via lowBound=0 in target_vars)

        # Objective: minimize conflict AND maximize peace
        # Penalize conflict indicators (positive weights = minimize)
        conflict_weights = {'AG_RL': 10.0, 'AG_SL': 10.0, 'V': 5.0, 'WR_SL': 2.0}
        conflict_score = pulp.lpSum([
            conflict_weights.get(sp, 0) * target_vars[sp] for sp in SPECIES
        ])

        # Reward peace indicators (negative weights = maximize)
        peace_weights = {'Gov': -4.0, 'T': -3.0, 'E': -2.0, 'SR_RL': -1.0, 'SR_SL': -0.5}
        peace_score = pulp.lpSum([
            peace_weights.get(sp, 0) * target_vars[sp] for sp in SPECIES
        ])

        # Small reward for budget utilization (encourages action)
        budget_utilization = -0.01 * gov_cost

        prob += conflict_score + peace_score + budget_utilization

        # Solve
        solver = pulp.PULP_CBC_CMD(msg=False, timeLimit=10)
        status = prob.solve(solver)

        results = {
            'optimal': status == pulp.LpStatusOptimal,
            'strategy': strategy,
            'budget_used': pulp.value(gov_cost) if status == pulp.LpStatusOptimal else 0,
            'initial_state': self.game_state.state.copy(),
            'allowed_rxns': allowed_rxns,
            'lp_status': pulp.LpStatus[status]
        }

        # Track ALL reaction fluxes for process analysis
        results['all_fluxes'] = {}
        for rxn in REACTIONS:
            if rxn in RXN_IDX:
                val = pulp.value(flux_vars[RXN_IDX[rxn]])
                results['all_fluxes'][rxn] = val if val is not None else 0.0

        if results['optimal']:
            results['final_state'] = {sp: pulp.value(target_vars[sp]) for sp in SPECIES}

            # Track Gov policy allocations
            results['gov_allocations'] = {}
            for rxn in allowed_rxns:
                if rxn in RXN_IDX:
                    val = pulp.value(flux_vars[RXN_IDX[rxn]]) or 0
                    results['gov_allocations'][rxn] = val

            # Compute state changes
            results['state_changes'] = {}
            for sp in SPECIES:
                initial = self.game_state.state.get(sp, 0)
                final = results['final_state'].get(sp, 0)
                results['state_changes'][sp] = final - initial

            # Compute loop activities
            results['peace_loops'] = {}
            for loop_name, loop_info in PEACE_LOOPS.items():
                total_flux = sum(results['all_fluxes'].get(r, 0) for r in loop_info['reactions'])
                if total_flux > 0.001:
                    results['peace_loops'][loop_name] = {
                        'flux': total_flux,
                        'reactions': {r: results['all_fluxes'].get(r, 0)
                                     for r in loop_info['reactions']
                                     if results['all_fluxes'].get(r, 0) > 0.001},
                        'description': loop_info['description'],
                        'icon': loop_info['icon']
                    }

            results['conflict_loops'] = {}
            for loop_name, loop_info in CONFLICT_LOOPS.items():
                total_flux = sum(results['all_fluxes'].get(r, 0) for r in loop_info['reactions'])
                if total_flux > 0.001:
                    results['conflict_loops'][loop_name] = {
                        'flux': total_flux,
                        'reactions': {r: results['all_fluxes'].get(r, 0)
                                     for r in loop_info['reactions']
                                     if results['all_fluxes'].get(r, 0) > 0.001},
                        'description': loop_info['description'],
                        'icon': loop_info['icon']
                    }

            # Violence reduction metric
            initial_v = self.game_state.state.get('V', 0)
            final_v = results['final_state'].get('V', 0)
            results['violence_reduction'] = ((initial_v - final_v) / initial_v * 100) if initial_v > 0 else 0

        return results

    def _create_failure_result(self, strategy: str, _player_budget: float) -> Dict[str, Any]:
        """Create diagnostic result when all LP attempts fail"""
        state = self.game_state.state

        # Analyze why it might be failing
        diagnostics = []

        # Check for depleted resources
        for sp in ['E', 'Gov', 'T', 'RL']:
            if state.get(sp, 0) < 0.1:
                diagnostics.append(f"{sp} critically low ({state.get(sp, 0):.3f})")

        # Check for population imbalance
        pop_species = ['SR_RL', 'SR_SL', 'WR_RL', 'WR_SL', 'AG_RL', 'AG_SL']
        for sp in pop_species:
            if state.get(sp, 0) < 0.01:
                diagnostics.append(f"{sp} depleted")

        # Check for high violence
        if state.get('V', 0) > 3.0:
            diagnostics.append(f"Violence very high ({state.get('V', 0):.2f})")

        return {
            'optimal': False,
            'strategy': strategy,
            'budget_used': 0,
            'initial_state': state.copy(),
            'allowed_rxns': self.get_strategy_reactions(strategy),
            'all_fluxes': {rxn: 0.0 for rxn in REACTIONS if rxn in RXN_IDX},
            'lp_status': 'Infeasible',
            'diagnostics': diagnostics,
            'failure_reason': "System state has evolved to an unrecoverable configuration. " +
                            "; ".join(diagnostics) if diagnostics else "Constraint conflicts detected."
        }
    
    def ag_strategy_decision(self) -> Dict[str, Any]:
        """Determine AG's strategy based on game state"""
        state = self.game_state.state
        
        if self.config.AG_STRATEGY == "aggressive":
            # AG focuses on attacks when strong
            if state['AG_RL'] + state['AG_SL'] > 0.1:
                return {'priority': 'attack', 'budget_focus': 0.7}
            else:
                return {'priority': 'recruit', 'budget_focus': 0.5}
        
        elif self.config.AG_STRATEGY == "defensive":
            # AG focuses on survival when weak
            if state['Gov'] > 0.5:
                return {'priority': 'hide', 'budget_focus': 0.3}
            else:
                return {'priority': 'expand', 'budget_focus': 0.6}
        
        else:  # adaptive
            # AG adapts to situation
            violence_factor = state['V'] / (state['V'] + 1)
            gov_factor = state['Gov']
            
            if violence_factor > 0.6 and gov_factor < 0.4:
                return {'priority': 'attack', 'budget_focus': 0.8}
            elif state['AG_RL'] + state['AG_SL'] < 0.05:
                return {'priority': 'recruit', 'budget_focus': 0.9}
            elif gov_factor > 0.6:
                return {'priority': 'hide', 'budget_focus': 0.4}
            else:
                return {'priority': 'expand', 'budget_focus': 0.6}
    
    def play_round(self, player_budget: float, strategy: str, climate: float):
        """Play one round of the game"""
        print(f"\n{'='*60}")
        print(f"ROUND {self.game_state.round + 1}")
        print('='*60)
        
        # Update game state
        self.game_state.round += 1
        self.game_state.climate = climate
        self.game_state.player_strategy = strategy
        self.game_state.player_budget_allocated = player_budget
        
        print(f"Peace Budget: {self.game_state.peace_budget:.1f}")
        print(f"AG Budget: {self.game_state.ag_budget:.1f}")
        print(f"Climate: {climate:.2f} ({'Harsh' if climate < 0.4 else 'Favorable' if climate > 0.6 else 'Normal'})")
        print(f"Peace Index: {self.game_state.calculate_peace_index():.3f}")
        
        # Player's turn
        print(f"\n[Your Turn] Strategy: {strategy}, Budget: {player_budget:.1f}")
        player_result = self.solve_player_turn(player_budget, strategy)

        if not player_result['optimal']:
            print("⚠️  Your strategy had difficulties!")
            if 'diagnostics' in player_result:
                for diag in player_result['diagnostics'][:3]:
                    print(f"   • {diag}")
            # Still apply some effects
            player_result['final_state'] = self.game_state.state.copy()
        else:
            msg = f"✓ Strategy successful! Budget used: {player_result['budget_used']:.1f}"
            if player_result.get('relaxed'):
                msg += " (constraints relaxed)"
            print(msg)
        
        # Update state with player's actions
        intermediate_state = player_result.get('final_state', self.game_state.state.copy())
        
        # AG's turn (simplified - AG spends its budget)
        print(f"\n[AG's Turn] Budget: {self.game_state.ag_budget:.1f}")
        ag_plan = self.ag_strategy_decision()
        print(f"AG priority: {ag_plan['priority']}")
        
        # Apply AG effects (simplified simulation)
        ag_effect = self.simulate_ag_effects(intermediate_state, self.game_state.ag_budget, ag_plan)
        final_state = ag_effect['final_state']
        
        # Update game state
        old_state = self.game_state.state.copy()
        self.game_state.state = final_state
        
        # Update budgets for next round
        self.game_state.peace_budget -= player_result.get('budget_used', 0)
        self.game_state.peace_budget += self.config.BUDGET_INCREMENT_PER_ROUND
        self.game_state.ag_budget += self.config.BUDGET_INCREMENT_PER_ROUND
        
        # Record history with full process information
        history_entry = {
            'round': self.game_state.round,
            'player_budget': player_result.get('budget_used', 0),
            'player_strategy': strategy,
            'ag_budget': self.game_state.ag_budget,
            'ag_priority': ag_plan['priority'],
            'climate': climate,
            'peace_index': self.game_state.calculate_peace_index(),
            'state_changes': {sp: final_state[sp] - old_state[sp] for sp in SPECIES},
            # Process information from player
            'player_result': player_result,
            'player_gov_allocations': player_result.get('gov_allocations', {}),
            'player_peace_loops': player_result.get('peace_loops', {}),
            'player_conflict_loops': player_result.get('conflict_loops', {}),
            'player_all_fluxes': player_result.get('all_fluxes', {}),
            'violence_reduction': player_result.get('violence_reduction', 0),
            # Process information from AG
            'ag_result': ag_effect,
            'ag_allocations': ag_effect.get('ag_allocations', {}),
            'ag_conflict_loops': ag_effect.get('conflict_loops', {}),
            'ag_all_fluxes': ag_effect.get('all_fluxes', {})
        }
        self.game_state.history.append(history_entry)

        return history_entry
    
    def solve_ag_turn(self, state: Dict[str, float], ag_budget: float, ag_plan: Dict) -> Dict[str, Any]:
        """Solve LP for AG's turn with relaxed constraints"""
        # Try with progressively relaxed constraints
        for pop_tolerance in [0.01, 0.1, 0.5]:
            result = self._solve_ag_lp(state, ag_budget, ag_plan, pop_tolerance)
            if result['optimal']:
                return result

        # All attempts failed - return fallback
        return {
            'optimal': False,
            'ag_plan': ag_plan,
            'budget_used': 0,
            'initial_state': state.copy(),
            'final_state': state.copy(),
            'all_fluxes': {rxn: 0.0 for rxn in REACTIONS if rxn in RXN_IDX},
            'ag_allocations': {},
            'state_changes': {sp: 0 for sp in SPECIES},
            'conflict_loops': {}
        }

    def _solve_ag_lp(self, state: Dict[str, float], ag_budget: float,
                     ag_plan: Dict, pop_tolerance: float) -> Dict[str, Any]:
        """Internal AG LP solver with configurable constraint tightness"""
        prob = pulp.LpProblem("AG_Turn", pulp.LpMaximize)  # AG maximizes conflict

        n_reactions = len(REACTIONS)
        flux_vars = [pulp.LpVariable(f"v_{i}", lowBound=0) for i in range(n_reactions)]
        target_vars = {sp: pulp.LpVariable(f"t_{sp}", lowBound=0) for sp in SPECIES}

        # Calculate rates and costs
        fixed_rates = self.calculate_fixed_rates(state, self.game_state.climate)
        policy_costs = self.calculate_policy_costs(state)

        # Stoichiometric constraints
        for i, sp in enumerate(SPECIES):
            total_change = pulp.lpSum([S[i, j] * flux_vars[j] for j in range(n_reactions)])
            prob += target_vars[sp] == state.get(sp, 0) + total_change

        # Fixed reaction bounds
        TIME_SCALE = 0.1
        for rxn in FIXED_RXNS:
            if rxn in RXN_IDX:
                idx = RXN_IDX[rxn]
                rate = fixed_rates.get(rxn, 0) * TIME_SCALE
                if rate > 0:
                    prob += flux_vars[idx] <= rate

        # AG can only use AG policy reactions
        for rxn in GOV_POLICY_RXNS:
            if rxn in RXN_IDX:
                prob += flux_vars[RXN_IDX[rxn]] == 0

        # Budget constraint for AG
        ag_cost = pulp.lpSum([
            policy_costs.get(rxn, 0) * flux_vars[RXN_IDX[rxn]]
            for rxn in AG_POLICY_RXNS if rxn in RXN_IDX
        ])
        prob += ag_cost <= ag_budget * ag_plan['budget_focus']

        # Population conservation with TOLERANCE
        pop_species = ['SR_RL', 'SR_SL', 'WR_RL', 'WR_SL', 'AG_RL', 'AG_SL']
        initial_pop = sum(state.get(sp, 0) for sp in pop_species)
        final_pop = pulp.lpSum([target_vars[sp] for sp in pop_species])
        prob += final_pop >= initial_pop * (1 - pop_tolerance)
        prob += final_pop <= initial_pop * (1 + pop_tolerance)

        # AG objective depends on priority
        # Base conflict objectives + budget utilization incentive
        budget_bonus = 0.01 * ag_cost  # Incentivize using budget

        if ag_plan['priority'] == 'attack':
            # Maximize violence and gov reduction
            objective = (3.0 * target_vars['V']
                        - 2.0 * target_vars['Gov']
                        - target_vars['T']
                        - 0.5 * target_vars['E']
                        + budget_bonus)
        elif ag_plan['priority'] == 'recruit':
            # Maximize AG population, destabilize WR
            objective = (5.0 * (target_vars['AG_RL'] + target_vars['AG_SL'])
                        + 0.5 * target_vars['V']
                        - target_vars['WR_SL']  # Consume WR into AG
                        + budget_bonus)
        elif ag_plan['priority'] == 'expand':
            # Maximize AG_RL and extract resources
            objective = (3.0 * target_vars['AG_RL']
                        + target_vars['AG_SL']
                        - 0.5 * target_vars['E']
                        - 0.3 * target_vars['Gov']
                        + budget_bonus)
        else:  # hide
            # Minimize losses, survive
            objective = (target_vars['AG_RL'] + target_vars['AG_SL']
                        - 0.05 * target_vars['V']
                        + 0.5 * budget_bonus)  # Less incentive to act when hiding

        prob += objective

        # Solve
        solver = pulp.PULP_CBC_CMD(msg=False, timeLimit=10)
        status = prob.solve(solver)

        results = {
            'optimal': status == pulp.LpStatusOptimal,
            'ag_plan': ag_plan,
            'budget_used': pulp.value(ag_cost) if status == pulp.LpStatusOptimal else 0,
            'initial_state': state.copy()
        }

        # Track ALL reaction fluxes
        results['all_fluxes'] = {}
        for rxn in REACTIONS:
            if rxn in RXN_IDX:
                val = pulp.value(flux_vars[RXN_IDX[rxn]])
                results['all_fluxes'][rxn] = val if val is not None else 0.0

        if results['optimal']:
            results['final_state'] = {sp: pulp.value(target_vars[sp]) for sp in SPECIES}

            # Track AG allocations
            results['ag_allocations'] = {}
            for rxn in AG_POLICY_RXNS:
                if rxn in RXN_IDX:
                    val = pulp.value(flux_vars[RXN_IDX[rxn]]) or 0
                    results['ag_allocations'][rxn] = val

            # Compute state changes
            results['state_changes'] = {}
            for sp in SPECIES:
                initial = state.get(sp, 0)
                final = results['final_state'].get(sp, 0)
                results['state_changes'][sp] = final - initial

            # Compute conflict loop activities
            results['conflict_loops'] = {}
            for loop_name, loop_info in CONFLICT_LOOPS.items():
                total_flux = sum(results['all_fluxes'].get(r, 0) for r in loop_info['reactions'])
                if total_flux > 0.001:
                    results['conflict_loops'][loop_name] = {
                        'flux': total_flux,
                        'reactions': {r: results['all_fluxes'].get(r, 0)
                                     for r in loop_info['reactions']
                                     if results['all_fluxes'].get(r, 0) > 0.001},
                        'description': loop_info['description'],
                        'icon': loop_info['icon']
                    }

        return results

    def simulate_ag_effects(self, state: Dict[str, float], ag_budget: float, ag_plan: Dict) -> Dict[str, Any]:
        """Wrapper that uses LP solver for AG turn"""
        return self.solve_ag_turn(state, ag_budget, ag_plan)
    
    def display_state(self):
        """Display current game state"""
        print(f"\n{'='*60}")
        print("CURRENT GAME STATE")
        print('='*60)
        
        state = self.game_state.state
        peace_index = self.game_state.calculate_peace_index()
        
        print(f"\nPeace Index: {peace_index:.3f} ", end="")
        if peace_index >= 0.7:
            print("✅ (Good)")
        elif peace_index >= 0.4:
            print("⚠️  (Moderate)")
        else:
            print("❌ (Critical)")
        
        print(f"\nKey Indicators:")
        print(f"  Violence (V):       {state['V']:.2f}")
        print(f"  Armed Groups:       {state['AG_RL'] + state['AG_SL']:.2f}")
        print(f"  Governance:         {state['Gov']:.2f}")
        print(f"  Economy:            {state['E']:.2f}")
        print(f"  Trust:              {state['T']:.2f}")
        print(f"  Secure Population:  {state['SR_RL'] + state['SR_SL']:.1f}")
        print(f"  Vulnerable Pop:     {state['WR_RL'] + state['WR_SL']:.1f}")
        print(f"  Resource Land:      {state['RL']:.1f}")
        print(f"  Scarce Land:        {state['SL']:.1f}")
        
        print(f"\nResources:")
        print(f"  Your Budget:        {self.game_state.peace_budget:.1f}")
        print(f"  AG Budget:          {self.game_state.ag_budget:.1f}")
        print(f"  Round:              {self.game_state.round}/{self.config.MAX_ROUNDS}")
        print(f"  Climate:            {self.game_state.climate:.2f}")

# ==================== GAME UI ====================

class ConflictGame:
    """Main game interface"""
    
    def __init__(self):
        self.config = GameConfig()
        self.engine = GameEngine(self.config)
        self.running = True
        
    def show_help(self):
        """Display help information"""
        print(f"\n{'='*60}")
        print("CONFLICT RESOLUTION GAME - HELP")
        print('='*60)

        print("\n🎯 OBJECTIVE:")
        print("  Achieve peace (Peace Index >= {}) while minimizing budget use.".format(
            self.config.PEACE_VICTORY_THRESHOLD))
        print("  Avoid collapse (Peace Index <= {}).".format(
            self.config.PEACE_DEFEAT_THRESHOLD))

        print("\n💰 SCORING (if you win):")
        print("  - Higher score for using LESS budget")
        print("  - Higher score in HARSHER climate")
        print("  - Higher score for FEWER rounds")

        print("\n⚔️  STRATEGIES & THEIR REACTIONS:")
        for strategy_name, reactions in self.engine.strategies.items():
            print(f"\n  {strategy_name}:")
            for rxn in reactions:
                name = REACTION_NAMES.get(rxn, rxn)
                print(f"    {rxn}: {name}")

        print("\n🔄 KEY PEACE LOOPS:")
        for loop_name, loop_info in PEACE_LOOPS.items():
            print(f"  {loop_info['icon']} {loop_name}: {loop_info['description']}")

        print("\n💀 KEY CONFLICT LOOPS (AG uses these):")
        for loop_name, loop_info in CONFLICT_LOOPS.items():
            print(f"  {loop_info['icon']} {loop_name}: {loop_info['description']}")

        print("\n🌦️  CLIMATE EFFECTS:")
        print(f"  - Range: {self.config.CLIMATE_MIN} (Harsh) to {self.config.CLIMATE_MAX} (Favorable)")
        print("  - Harsher climate makes all reactions more costly")
        print("  - Favorable climate helps land restoration")

        print("\n📊 PROCESS VISIBILITY:")
        print("  After each round you will see:")
        print("  - Which reactions YOU activated (and their flux rates)")
        print("  - Which peace loops were strengthened")
        print("  - Which reactions the AG activated")
        print("  - Which conflict loops fired")
        print("  - Net state changes (green = good, red = bad)")

        print("\n🎮 COMMANDS:")
        print("  play     - Play next round")
        print("  state    - Show current state")
        print("  loops    - Show active loop dynamics and potentials")
        print("  history  - Show game history with process details")
        print("  help     - Show this help")
        print("  save     - Save game")
        print("  load     - Load game")
        print("  quit     - Quit game")
    
    def get_player_input(self) -> Tuple[str, float, float]:
        """Get player's decisions for next round"""
        print(f"\n{'='*60}")
        print("YOUR DECISIONS FOR NEXT ROUND")
        print('='*60)
        
        # Strategy selection
        strategies = list(self.engine.strategies.keys())
        print("\nAvailable Strategies:")
        for i, strat in enumerate(strategies, 1):
            print(f"  {i}. {strat}")
        
        while True:
            try:
                choice = input(f"\nChoose strategy (1-{len(strategies)} or 'back'): ").strip()
                if choice.lower() == 'back':
                    return None, None, None
                
                if choice.isdigit():
                    idx = int(choice) - 1
                    if 0 <= idx < len(strategies):
                        strategy = strategies[idx]
                        break
                
                print(f"Please enter a number between 1 and {len(strategies)} or 'back'")
            except (ValueError, IndexError):
                print("Invalid choice")
        
        # Budget allocation
        max_budget = self.engine.game_state.peace_budget
        print(f"\nAvailable Budget: {max_budget:.1f}")
        while True:
            try:
                budget_input = input(f"Budget to allocate (0-{max_budget:.1f}): ").strip()
                if budget_input.lower() == 'back':
                    return None, None, None
                
                budget = float(budget_input)
                if 0 <= budget <= max_budget:
                    break
                else:
                    print(f"Budget must be between 0 and {max_budget:.1f}")
            except ValueError:
                print("Please enter a valid number")
        
        # Climate setting
        print(f"\nClimate Settings: {self.config.CLIMATE_MIN} (Harsh) to {self.config.CLIMATE_MAX} (Favorable)")
        print(f"Current: {self.engine.game_state.climate:.2f}")
        while True:
            try:
                climate_input = input(f"Set climate ({self.config.CLIMATE_MIN}-{self.config.CLIMATE_MAX}): ").strip()
                if climate_input.lower() == 'back':
                    return None, None, None
                
                climate = float(climate_input)
                if self.config.CLIMATE_MIN <= climate <= self.config.CLIMATE_MAX:
                    break
                else:
                    print(f"Climate must be between {self.config.CLIMATE_MIN} and {self.config.CLIMATE_MAX}")
            except ValueError:
                print("Please enter a valid number")
        
        return strategy, budget, climate
    
    def show_round_results(self, round_data: Dict):
        """Display results of the round with full process analysis"""
        print(f"\n{'='*60}")
        print(f"ROUND {round_data['round']} RESULTS")
        print('='*60)

        # Overall status
        print(f"\nPeace Index: {round_data['peace_index']:.3f} ", end="")
        if round_data['peace_index'] >= 0.7:
            print("✅ (Progress toward peace!)")
        elif round_data['peace_index'] <= 0.3:
            print("🚨 (Conflict worsening!)")
        else:
            print("⚠️  (Contested)")

        print(f"Strategy: {round_data['player_strategy']}")
        print(f"Climate: {round_data['climate']:.2f}")

        # Violence reduction
        v_red = round_data.get('violence_reduction', 0)
        if v_red != 0:
            print(f"Violence Change: {v_red:+.1f}%")

        # ===== YOUR ACTIONS (PEACE SIDE) =====
        print(f"\n{'─'*60}")
        print("🕊️  YOUR PEACE ACTIONS")
        print(f"{'─'*60}")
        print(f"Budget Used: {round_data['player_budget']:.1f}")

        # Show if constraints were relaxed
        player_result = round_data.get('player_result', {})
        if player_result.get('relaxed'):
            print(f"  ⚠️  (Constraints relaxed to find solution)")

        # Show diagnostics if strategy failed
        if not player_result.get('optimal', True):
            print(f"\n  ❌ Strategy encountered issues:")
            for diag in player_result.get('diagnostics', [])[:3]:
                print(f"     • {diag}")
            failure_reason = player_result.get('failure_reason', '')
            if failure_reason:
                print(f"     Reason: {failure_reason[:60]}...")

        # Show government policy allocations (reactions used)
        gov_alloc = round_data.get('player_gov_allocations', {})
        active_gov = {r: v for r, v in gov_alloc.items() if v > 0.001}

        if active_gov:
            print(f"\nReactions Activated:")
            for rxn, flux in sorted(active_gov.items(), key=lambda x: -x[1]):
                name = REACTION_NAMES.get(rxn, rxn)
                print(f"  {rxn}: {flux:6.3f} - {name}")
        else:
            print(f"\n  (No policy reactions activated)")

        # Show active peace loops
        peace_loops = round_data.get('player_peace_loops', {})
        if peace_loops:
            print(f"\nActive Peace Loops:")
            for loop_name, loop_data in peace_loops.items():
                icon = loop_data.get('icon', '•')
                flux = loop_data.get('flux', 0)
                desc = loop_data.get('description', '')
                print(f"  {icon} {loop_name}: {flux:.3f}")
                print(f"      {desc}")

        # ===== AG ACTIONS (CONFLICT SIDE) =====
        print(f"\n{'─'*60}")
        print("💀 ARMED GROUP ACTIONS")
        print(f"{'─'*60}")
        print(f"AG Priority: {round_data.get('ag_priority', 'unknown')}")

        # Show AG policy allocations
        ag_alloc = round_data.get('ag_allocations', {})
        active_ag = {r: v for r, v in ag_alloc.items() if v > 0.001}

        if active_ag:
            print(f"\nReactions Activated:")
            for rxn, flux in sorted(active_ag.items(), key=lambda x: -x[1]):
                name = REACTION_NAMES.get(rxn, rxn)
                print(f"  {rxn}: {flux:6.3f} - {name}")
        else:
            print(f"\n  (AG could not activate reactions)")

        # Show active conflict loops
        conflict_loops = round_data.get('ag_conflict_loops', {})
        if conflict_loops:
            print(f"\nActive Conflict Loops:")
            for loop_name, loop_data in conflict_loops.items():
                icon = loop_data.get('icon', '•')
                flux = loop_data.get('flux', 0)
                desc = loop_data.get('description', '')
                print(f"  {icon} {loop_name}: {flux:.3f}")
                print(f"      {desc}")

        # ===== NET STATE CHANGES =====
        print(f"\n{'─'*60}")
        print("📊 NET STATE CHANGES")
        print(f"{'─'*60}")

        changes = round_data.get('state_changes', {})
        # Order by importance
        priority_species = ['AG_RL', 'AG_SL', 'V', 'Gov', 'T', 'E', 'WR_SL', 'SR_RL']
        significant = {k: v for k, v in changes.items() if abs(v) > 0.005}

        for sp in priority_species:
            if sp in significant:
                change = significant[sp]
                direction = "↑" if change > 0 else "↓"
                # Color code: green for peace-positive, red for conflict-positive
                if sp in ['AG_RL', 'AG_SL', 'V', 'WR_SL']:
                    status = "✅" if change < 0 else "🔺"  # Less is better
                else:
                    status = "✅" if change > 0 else "🔻"  # More is better
                print(f"  {status} {sp:8s}: {direction} {abs(change):.3f}")

        # Show other changes
        other_changes = {k: v for k, v in significant.items() if k not in priority_species}
        for sp, change in other_changes.items():
            direction = "↑" if change > 0 else "↓"
            print(f"     {sp:8s}: {direction} {abs(change):.3f}")
    
    def check_game_over(self):
        """Check if game is over and handle end game"""
        status = self.engine.game_state.get_game_status()
        
        if status == "ONGOING":
            return False
        
        print(f"\n{'='*60}")
        print("GAME OVER")
        print('='*60)
        
        final_peace = self.engine.game_state.calculate_peace_index()
        
        if status == "VICTORY":
            score = self.engine.game_state.calculate_score()
            print(f"\n🎉 VICTORY! 🎉")
            print(f"Final Peace Index: {final_peace:.3f}")
            print(f"Rounds: {self.engine.game_state.round}")
            print(f"\n🏆 FINAL SCORE: {score:.0f}")
            print(f"\nScore Breakdown:")
            total_budget_used = sum(h['player_budget'] for h in self.engine.game_state.history)
            max_possible = self.config.INITIAL_PEACE_BUDGET + self.config.BUDGET_INCREMENT_PER_ROUND * self.engine.game_state.round
            print(f"  • Budget Efficiency: {1 - total_budget_used/max_possible:.1%}")
            avg_climate = np.mean([h['climate'] for h in self.engine.game_state.history])
            print(f"  • Climate Challenge: {avg_climate:.2f} avg")
            print(f"  • Round Efficiency: {self.engine.game_state.round}/{self.config.MAX_ROUNDS}")
        
        elif status == "DEFEAT":
            print(f"\n💥 DEFEAT! 💥")
            print(f"Peace collapsed below threshold: {final_peace:.3f}")
            print(f"\nThe conflict spiraled out of control.")
            print(f"Armed groups strengthened: {self.engine.game_state.state['AG_RL'] + self.engine.game_state.state['AG_SL']:.2f}")
            print(f"Violence level: {self.engine.game_state.state['V']:.2f}")
        
        elif status == "TIMEOUT":
            print(f"\n⏰ TIME'S UP! ⏰")
            print(f"Maximum rounds ({self.config.MAX_ROUNDS}) reached.")
            print(f"Final Peace Index: {final_peace:.3f}")
            if final_peace > 0.5:
                print("Partial success, but not complete peace.")
            else:
                print("Conflict persists.")
        
        # Show final statistics
        print(f"\n{'='*60}")
        print("FINAL STATISTICS")
        print('='*60)
        
        history = self.engine.game_state.history
        if history:
            total_budget = sum(h['player_budget'] for h in history)
            avg_climate = np.mean([h['climate'] for h in history])
            strategies_used = set(h['player_strategy'] for h in history)
            
            print(f"\nTotal Budget Used: {total_budget:.1f}")
            print(f"Average Climate: {avg_climate:.2f}")
            print(f"Strategies Used: {', '.join(strategies_used)}")
            print(f"Peace Index Change: {history[0]['peace_index']:.3f} → {final_peace:.3f}")
        
        return True
    
    def run(self):
        """Main game loop"""
        print(f"\n{'='*60}")
        print("CONFLICT RESOLUTION GAME")
        print("Achieve peace with minimal resources")
        print('='*60)
        
        # Check for saved game
        if os.path.exists(self.config.save_file):
            load = input("Load saved game? (y/n): ").lower().strip()
            if load == 'y':
                if not self.engine.game_state.load_game():
                    print("Starting new game...")
        
        # Show initial state
        self.engine.display_state()
        
        # Main game loop
        while self.running:
            print(f"\n{'='*60}")
            print("MAIN MENU")
            print('='*60)
            print("Commands: play, state, loops, history, help, save, load, quit")
            
            command = input("\nCommand: ").lower().strip()
            
            if command == 'play':
                # Check if game is already over
                if self.check_game_over():
                    play_again = input("\nPlay again? (y/n): ").lower().strip()
                    if play_again == 'y':
                        self.engine = GameEngine(self.config)
                        self.engine.display_state()
                    else:
                        self.running = False
                    continue
                
                # Get player input
                strategy, budget, climate = self.get_player_input()
                if strategy is None:  # User chose 'back'
                    continue
                
                # Play round
                round_data = self.engine.play_round(budget, strategy, climate)
                self.show_round_results(round_data)
                
                # Check if game over
                if self.check_game_over():
                    play_again = input("\nPlay again? (y/n): ").lower().strip()
                    if play_again == 'y':
                        self.engine = GameEngine(self.config)
                    else:
                        self.running = False
            
            elif command == 'state':
                self.engine.display_state()
            
            elif command == 'help':
                self.show_help()

            elif command == 'loops':
                self.show_loop_dynamics()

            elif command == 'history':
                self.show_history()

            elif command == 'save':
                self.engine.game_state.save_game()

            elif command == 'load':
                if self.engine.game_state.load_game():
                    self.engine.display_state()

            elif command == 'quit':
                save = input("Save game before quitting? (y/n): ").lower().strip()
                if save == 'y':
                    self.engine.game_state.save_game()
                print("Thanks for playing!")
                self.running = False

            else:
                print("Unknown command. Type 'help' for commands.")

    def show_loop_dynamics(self):
        """Show the current loop dynamics based on state"""
        print(f"\n{'='*60}")
        print("CURRENT LOOP DYNAMICS")
        print('='*60)

        state = self.engine.game_state.state
        history = self.engine.game_state.history

        # Calculate potential loop strengths based on current state
        print("\n🕊️  PEACE LOOP POTENTIALS:")
        print("  (Higher = more effective if activated)")

        # Land Restoration potential
        land_potential = state['SL'] * self.engine.game_state.climate
        print(f"  🌱 Land_Restoration: {land_potential:.2f} (SL={state['SL']:.1f}, climate={self.engine.game_state.climate:.2f})")

        # Economic Development potential
        econ_potential = (state['SR_RL'] + state['WR_RL']) / 10.0 * state.get('E', 1)
        print(f"  💰 Economic_Development: {econ_potential:.2f} (Workers={state['SR_RL']+state['WR_RL']:.1f})")

        # Governance Building potential
        gov_potential = state['Gov'] * state['E']
        print(f"  🏛️  Governance_Building: {gov_potential:.2f} (Gov={state['Gov']:.2f}, E={state['E']:.2f})")

        # DDR potential
        ddr_potential = state['Gov'] * state['E'] / (state['AG_SL'] + state['AG_RL'] + 0.1)
        print(f"  🕊️  DDR_Reintegration: {ddr_potential:.2f} (requires Gov+E)")

        print("\n💀 CONFLICT LOOP POTENTIALS:")
        print("  (Higher = more dangerous)")

        # Recruitment potential
        ag_total = state['AG_SL'] + state['AG_RL']
        wr_pool = state['WR_SL'] + state['WR_RL']
        recruitment_pot = ag_total * wr_pool * state['V'] / (5.0 + wr_pool)
        print(f"  📢 Recruitment: {recruitment_pot:.3f} (AG={ag_total:.2f}, WR={wr_pool:.1f}, V={state['V']:.2f})")

        # Predation potential
        predation_pot = ag_total * (state['SR_RL'] + state['E']) * (1 - state['Gov'])
        print(f"  💀 Predation: {predation_pot:.2f} (AG strength vs weak governance)")

        # Violence generation
        grievance = state['WR_SL'] * state['WR_SL'] / (10 + state['WR_SL'] * state['WR_SL'])
        print(f"  🔥 Violence_Generation: {grievance:.2f} (grievance from WR_SL)")

        # State erosion
        erosion_pot = ag_total * (1 - state['Gov']) * state['V']
        print(f"  🏚️  State_Erosion: {erosion_pot:.2f} (AG vs Gov)")

        # Show last round's actual loop activities if available
        if history:
            last_round = history[-1]
            print(f"\n📊 LAST ROUND ACTUAL ACTIVITY:")

            peace = last_round.get('player_peace_loops', {})
            if peace:
                print("  Peace loops activated:")
                for loop, data in peace.items():
                    print(f"    {data.get('icon', '•')} {loop}: {data.get('flux', 0):.3f}")

            conflict = last_round.get('ag_conflict_loops', {})
            if conflict:
                print("  Conflict loops activated:")
                for loop, data in conflict.items():
                    print(f"    {data.get('icon', '•')} {loop}: {data.get('flux', 0):.3f}")

    def show_history(self):
        """Show game history with process information"""
        print(f"\n{'='*60}")
        print("GAME HISTORY")
        print('='*60)

        history = self.engine.game_state.history
        if not history:
            print("\nNo rounds played yet.")
            return

        for entry in history:
            print(f"\n--- Round {entry['round']} ---")
            print(f"  Strategy: {entry['player_strategy']}, Budget: {entry['player_budget']:.1f}")
            print(f"  AG Priority: {entry.get('ag_priority', '?')}")
            print(f"  Peace Index: {entry['peace_index']:.3f}")

            # Show key reactions used
            gov_alloc = entry.get('player_gov_allocations', {})
            active = [f"{r}({v:.2f})" for r, v in gov_alloc.items() if v > 0.01]
            if active:
                print(f"  Your reactions: {', '.join(active[:3])}")

            ag_alloc = entry.get('ag_allocations', {})
            ag_active = [f"{r}({v:.2f})" for r, v in ag_alloc.items() if v > 0.01]
            if ag_active:
                print(f"  AG reactions: {', '.join(ag_active[:3])}")

# ==================== MAIN ====================

if __name__ == "__main__":
    # Check for pulp installation
    try:
        import pulp
    except ImportError:
        print("Error: pulp library not installed.")
        print("Install it with: pip install pulp")
        sys.exit(1)
    
    # Run the game
    game = ConflictGame()
    game.run()