#!/usr/bin/env python3
"""
Conflict Dynamics Simulator with Budget-Based Resource Allocation

Features:
1. Manual consumption tracking for E, Gov, V (not via stoichiometry)
2. Geometric mean for multi-resource reactions (r21)
3. Properly scaled kinetic constants by reaction order
4. Rate limiter for conservation
5. Batch simulation across conflict levels, climate parameters, and strategies
6. Comprehensive comparison plots
"""

import os
import sys
import numpy as np
from math import sqrt
from typing import Dict, List, Tuple, Any, Optional
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# Add the root directory to the PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))
from pyCOT.io.functions import read_txt

# ============================================================================
# SCENARIO DEFINITIONS
# ============================================================================

CONFLICT_SCENARIOS = {
    "low": {
        "description": "Low Conflict - Stable region with minor tensions",
        "initial_state": {
            "SR_RL": 40.0,   # More strong resilient
            "SR_SL": 25.0,
            "WR_RL": 30.0,
            "WR_SL": 35.0,   # Less weak on stressed land
            "AG_RL": 0.5,    # Few armed groups
            "AG_SL": 2.0,
            "RL": 50.0,      # More restored land
            "SL": 50.0,
            "E": 100.0,      # Strong economy
            "T": 60.0,       # High trust
            "V": 10.0,       # Low violence
            "Gov": 30.0,     # Strong governance
        }
    },
    "medium": {
        "description": "Medium Conflict - Regional instability with active tensions",
        "initial_state": {
            "SR_RL": 30.0,
            "SR_SL": 20.0,
            "WR_RL": 25.0,
            "WR_SL": 50.0,
            "AG_RL": 0.8,
            "AG_SL": 2.5,
            "RL": 40.0,
            "SL": 60.0,
            "E": 80.0,
            "T": 40.0,
            "V": 30.0,
            "Gov": 20.0,
        }
    },
    "severe": {
        "description": "Severe Conflict - Active insurgency with widespread violence",
        "initial_state": {
            "SR_RL": 20.0,   # Few strong resilient
            "SR_SL": 15.0,
            "WR_RL": 20.0,
            "WR_SL": 65.0,   # Many weak on stressed land
            "AG_RL": 1.2,    # More armed groups
            "AG_SL": 3.0,
            "RL": 30.0,      # Less restored land
            "SL": 70.0,
            "E": 50.0,       # Weak economy
            "T": 20.0,       # Low trust
            "V": 50.0,       # High violence
            "Gov": 10.0,     # Weak governance
        }
    }
}

CLIMATE_PARAMETERS = [0.01, 0.02, 0.03, 0.04]  # From slow to fast climate dynamics

GOVERNMENT_STRATEGIES = ["development", "security", "balanced"]
ARMED_GROUP_STRATEGIES = ["recruitment", "displacement", "balanced"]

class ConflictSimulator:
    """Simple simulator for conflict dynamics with budget allocation."""

    def __init__(self, model_file: str):
        """Initialize simulator with reaction network."""
        self.rn = read_txt(model_file)
        self.S = self.rn.stoichiometry_matrix()
        self.SPECIES = list(self.rn._species_map.keys())
        self.REACTIONS = [f"r{i+1}" for i in range(self.S.shape[1])]

        # Find indices for E, Gov, V in species list (for manual updates)
        self.E_idx = self.SPECIES.index("E") if "E" in self.SPECIES else None
        self.Gov_idx = self.SPECIES.index("Gov") if "Gov" in self.SPECIES else None
        self.V_idx = self.SPECIES.index("V") if "V" in self.SPECIES else None

        # Map reactions to their descriptions
        self.reaction_descriptions = {
            # Land and Resource Dynamics
            "r1": "SL -> RL (climate regeneration)",
            "r2": "SR_SL -> SR_RL (climate regeneration)",
            "r3": "WR_SL -> WR_RL (climate regeneration)",
            "r4": "RL -> SL (climate degradation)",
            "r5": "SR_RL -> SR_SL (climate degradation)",
            "r6": "WR_RL -> WR_SL (climate degradation)",
            "r7": "AG_RL -> AG_SL (climate degradation)",
            "r8": "Gov + SR_SL -> SR_RL (gov land restoration)",

            # Economic Production
            "r9": "SR_RL -> SR_SL + E (economic production)",
            "r10": "WR_RL -> WR_SL + E (economic production)",
            "r11": "E -> Gov (economy to governance)",
            "r12": "-> Gov (baseline governance)",

            # Population Migration
            "r13": "SR_SL + RL -> SR_RL + SL (migration strong)",
            "r14": "WR_SL + RL -> WR_RL + SL (migration weak)",
            "r15": "SR_RL + SL + V -> RL + SR_SL (violent displacement strong)",
            "r16": "WR_RL + SL + V -> RL + WR_SL (violent displacement weak)",
            "r17": "AG_SL + RL -> AG_RL + SL (AG migration)",

            # Resilience Transitions
            "r18": "E + WR_RL -> SR_RL (economic resilience)",
            "r19": "E + WR_SL -> SR_SL (economic resilience)",
            "r20": "SR_SL -> WR_SL (resilience degradation)",
            "r21": "E + Gov + AG_SL -> WR_SL (reintegration)",

            # Trust and Social Cohesion
            "r22": "SR_RL + E -> SR_RL + T (trust from prosperity)",
            "r23": "Gov + SR_RL -> Gov + SR_RL + T (gov trust building)",
            "r24": "Gov + WR_RL -> SR_RL (gov outreach)",
            "r25": "V + T -> (trust destruction)",
            "r26": "Gov -> (gov corruption)",

            # Violence Generation
            "r27": "WR_SL -> WR_SL + V (desperation violence)",
            "r28": "AG_SL + SR_RL -> AG_RL + WR_SL + V (AG attacks)",
            "r29": "AG_SL + Gov -> WR_SL + V (gov-AG conflict)",
            "r30": "AG_RL + Gov -> AG_SL + V (gov displacement ops)",
            "r31": "V -> (violence decay)",

            # Armed Group Dynamics
            "r32": "WR_SL + AG_SL + V -> 2AG_SL (AG recruitment)",
            "r33": "AG_SL + WR_RL -> AG_RL + WR_SL (AG expansion weak)",
            "r34": "AG_SL + SR_RL -> AG_RL + WR_SL (AG expansion strong)",
            "r35": "AG_SL + E -> AG_RL (AG corruption)",
            "r36": "AG_RL + Gov -> AG_RL (AG institutional capture)",
        }

        # Define which reactions CONSUME allocated resources
        # These will be tracked manually, not via stoichiometry
        self.E_consuming_reactions = ["r11", "r18", "r19", "r21", "r35"]  # r35 uses total E
        self.Gov_consuming_reactions = ["r8", "r21", "r24", "r29", "r30", "r36"]  # r36 uses total Gov
        self.V_consuming_reactions = ["r15", "r16", "r25", "r31"]  # r31 uses total V

        # Define which reactions use ALLOCATED resources (strategy-dependent)
        self.E_allocated_reactions = ["r11", "r18", "r19", "r21", "r22"]
        self.Gov_allocated_reactions = ["r8", "r21", "r23", "r24", "r29", "r30"]
        self.V_allocated_reactions = ["r15", "r16", "r25", "r32"]

        # Set default strategies
        self.set_government_strategy("balanced")
        self.set_armed_groups_strategy("balanced")

        # Set kinetic parameters (per year)
        self.set_kinetic_parameters()

        # Initialize state
        self.state = self.get_initial_state()
        self.time = 0.0

        # History tracking
        self.time_history = [0.0]
        self.state_history = [self.state.copy()]
        self.rate_history = []
        self.allocation_history = []
        self.consumption_history = []

        print(f"Model loaded: {len(self.SPECIES)} species, {len(self.REACTIONS)} reactions")

    def reset(self, initial_state: Optional[Dict[str, float]] = None):
        """Reset simulator to initial state for new simulation run."""
        if initial_state is not None:
            self.state = initial_state.copy()
        else:
            self.state = self.get_initial_state()
        self.time = 0.0
        self.time_history = [0.0]
        self.state_history = [self.state.copy()]
        self.rate_history = []
        self.allocation_history = []
        self.consumption_history = []

    def set_conflict_scenario(self, scenario: str):
        """Set initial state based on conflict scenario."""
        if scenario not in CONFLICT_SCENARIOS:
            raise ValueError(f"Unknown scenario: {scenario}. Choose from: {list(CONFLICT_SCENARIOS.keys())}")
        self.current_scenario = scenario
        self.reset(CONFLICT_SCENARIOS[scenario]["initial_state"])

    def set_climate_parameter(self, climate_k: float):
        """Set climate-related kinetic parameters (r1-r7).

        Higher values = faster climate dynamics (both regeneration and degradation).
        """
        self.climate_k = climate_k
        # Climate regeneration (positive effects)
        self.k["r1"] = climate_k      # SL -> RL
        self.k["r2"] = climate_k      # SR_SL -> SR_RL
        self.k["r3"] = climate_k      # WR_SL -> WR_RL
        # Climate degradation (negative effects)
        self.k["r4"] = climate_k      # RL -> SL
        self.k["r5"] = climate_k      # SR_RL -> SR_SL
        self.k["r6"] = climate_k      # WR_RL -> WR_SL
        self.k["r7"] = climate_k * 2  # AG_RL -> AG_SL (AG more affected)

    def get_metrics(self) -> Dict[str, Any]:
        """Extract key metrics from simulation history."""
        if len(self.state_history) < 2:
            return {}

        # Time series
        times = np.array(self.time_history)

        # Extract species time series
        E_series = np.array([s["E"] for s in self.state_history])
        T_series = np.array([s["T"] for s in self.state_history])
        Gov_series = np.array([s["Gov"] for s in self.state_history])
        V_series = np.array([s["V"] for s in self.state_history])

        SR_total = np.array([s["SR_RL"] + s["SR_SL"] for s in self.state_history])
        WR_total = np.array([s["WR_RL"] + s["WR_SL"] for s in self.state_history])
        AG_total = np.array([s["AG_RL"] + s["AG_SL"] for s in self.state_history])

        RL_series = np.array([s["RL"] for s in self.state_history])

        # Cumulative budget consumption
        if self.consumption_history:
            E_consumed = np.cumsum([c["E"] * 0.1 for c in self.consumption_history])
            Gov_consumed = np.cumsum([c["Gov"] * 0.1 for c in self.consumption_history])
            V_consumed = np.cumsum([c["V"] * 0.1 for c in self.consumption_history])
            # Pad to match time_history length
            E_consumed = np.concatenate([[0], E_consumed])
            Gov_consumed = np.concatenate([[0], Gov_consumed])
            V_consumed = np.concatenate([[0], V_consumed])
        else:
            E_consumed = np.zeros_like(times)
            Gov_consumed = np.zeros_like(times)
            V_consumed = np.zeros_like(times)

        return {
            "times": times,
            "E": E_series,
            "T": T_series,
            "Gov": Gov_series,
            "V": V_series,
            "SR_total": SR_total,
            "WR_total": WR_total,
            "AG_total": AG_total,
            "RL": RL_series,
            "E_consumed": E_consumed,
            "Gov_consumed": Gov_consumed,
            "V_consumed": V_consumed,
            # Final values
            "final_E": E_series[-1],
            "final_T": T_series[-1],
            "final_Gov": Gov_series[-1],
            "final_V": V_series[-1],
            "final_SR": SR_total[-1],
            "final_AG": AG_total[-1],
            "final_RL": RL_series[-1],
            "total_E_consumed": E_consumed[-1],
            "total_Gov_consumed": Gov_consumed[-1],
            "total_V_consumed": V_consumed[-1],
        }

    def set_kinetic_parameters(self):
        """Set kinetic parameters (per year timescale).

        IMPORTANT: Rate constants must be scaled by reaction order!
        For mass action: rate = k * [A] * [B] * ...

        With concentrations ~10-100, we need:
        - 1st order (A ->): k ~ 0.01-0.1
        - 2nd order (A + B ->): k ~ 0.001-0.01  (divides by ~100)
        - 3rd order (A + B + C ->): k ~ 0.0001-0.001 (divides by ~10000)

        This ensures rates stay O(1) and don't violate conservation.
        """
        # Base rates (per year) - PROPERLY SCALED BY ORDER
        self.k = {
            # Land and Resource Dynamics (1st order)
            "r1": 0.02,   # SL -> RL (1st order)
            "r2": 0.02,   # SR_SL -> SR_RL (1st order)
            "r3": 0.02,   # WR_SL -> WR_RL (1st order)
            "r4": 0.02,   # RL -> SL (1st order)
            "r5": 0.02,   # SR_RL -> SR_SL (1st order)
            "r6": 0.02,   # WR_RL -> WR_SL (1st order)
            "r7": 0.05,   # AG_RL -> AG_SL (1st order, higher for AG scale)
            "r8": 0.002,  # Gov + SR_SL -> SR_RL (2nd order)

            # Economic Production (1st order)
            "r9": 0.05,   # SR_RL -> SR_SL + E (1st order)
            "r10": 0.03,  # WR_RL -> WR_SL + E (1st order)
            "r11": 0.02,  # E -> Gov (1st order)
            "r12": 1.0,   # -> Gov (0th order, constant input)

            # Population Migration (2nd order)
            "r13": 0.001, # SR_SL + RL -> SR_RL + SL (2nd order)
            "r14": 0.001, # WR_SL + RL -> WR_RL + SL (2nd order)
            "r15": 0.0001,# SR_RL + SL + V -> RL + SR_SL (3rd order!)
            "r16": 0.0001,# WR_RL + SL + V -> RL + WR_SL (3rd order!)
            "r17": 0.005, # AG_SL + RL -> AG_RL + SL (2nd order, higher for AG)

            # Resilience Transitions
            "r18": 0.001, # E + WR_RL -> SR_RL (2nd order)
            "r19": 0.001, # E + WR_SL -> SR_SL (2nd order)
            "r20": 0.02,  # SR_SL -> WR_SL (1st order)
            "r21": 0.001, # E + Gov + AG_SL -> WR_SL (3rd order! - reintegration)

            # Trust and Social Cohesion
            "r22": 0.001, # SR_RL + E -> SR_RL + T (2nd order)
            "r23": 0.002, # Gov + SR_RL -> Gov + SR_RL + T (2nd order)
            "r24": 0.002, # Gov + WR_RL -> SR_RL (2nd order)
            "r25": 0.002, # V + T -> (2nd order)
            "r26": 0.02,  # Gov -> (1st order, corruption)

            # Violence Generation
            "r27": 0.01,  # WR_SL -> WR_SL + V (1st order)
            "r28": 0.005, # AG_SL + SR_RL -> AG_RL + WR_SL + V (2nd order)
            "r29": 0.01,  # AG_SL + Gov -> WR_SL + V (2nd order)
            "r30": 0.02,  # AG_RL + Gov -> AG_SL + V (2nd order, higher for AG)
            "r31": 0.1,   # V -> (1st order, decay)

            # Armed Group Dynamics
            "r32": 0.0002,# WR_SL + AG_SL + V -> 2AG_SL (3rd order! recruitment)
            "r33": 0.005, # AG_SL + WR_RL -> AG_RL + WR_SL (2nd order)
            "r34": 0.005, # AG_SL + SR_RL -> AG_RL + WR_SL (2nd order)
            "r35": 0.002, # AG_SL + E -> AG_RL (2nd order, corruption)
            "r36": 0.005, # AG_RL + Gov -> AG_RL (2nd order, capture)
        }

    def get_initial_state(self) -> Dict[str, float]:
        """Return initial state for simulation."""
        return {
            "SR_RL": 30.0,   # 3 million strong resilient on restored land
            "SR_SL": 20.0,   # 2 million strong resilient on stressed land
            "WR_RL": 25.0,   # 2.5 million weak resilient on restored land
            "WR_SL": 50.0,   # 5 million weak resilient on stressed land
            "AG_RL": 1.0,    # 0.1 million armed group on restored land
            "AG_SL": 5.0,    # 0.5 million armed group on stressed land
            "RL": 40.0,      # 400,000 km2 restored land
            "SL": 60.0,      # 600,000 km2 stressed land
            "E": 80.0,       # 80 million USD economic capital
            "T": 40.0,       # Trust index
            "V": 30.0,       # Violence index
            "Gov": 20.0,     # Governance index
        }

    def set_government_strategy(self, strategy: str, quiet: bool = False):
        """Set government allocation strategy for E and Gov."""
        if strategy == "development":
            # Focus on development and trust building
            self.E_weights = {"r11": 0.2, "r18": 0.3, "r19": 0.2, "r21": 0.1, "r22": 0.2}
            self.Gov_weights = {"r8": 0.3, "r21": 0.2, "r23": 0.3, "r24": 0.2, "r29": 0.0, "r30": 0.0}
        elif strategy == "security":
            # Focus on security operations
            self.E_weights = {"r11": 0.1, "r18": 0.1, "r19": 0.1, "r21": 0.4, "r22": 0.3}
            self.Gov_weights = {"r8": 0.1, "r21": 0.2, "r23": 0.1, "r24": 0.1, "r29": 0.25, "r30": 0.25}
        elif strategy == "balanced":
            # Balanced approach
            self.E_weights = {"r11": 0.2, "r18": 0.2, "r19": 0.2, "r21": 0.2, "r22": 0.2}
            self.Gov_weights = {"r8": 0.2, "r21": 0.2, "r23": 0.2, "r24": 0.2, "r29": 0.1, "r30": 0.1}
        else:
            raise ValueError(f"Unknown government strategy: {strategy}")

        self.gov_strategy = strategy
        if not quiet:
            print(f"Government strategy set to: {strategy}")

    def set_armed_groups_strategy(self, strategy: str, quiet: bool = False):
        """Set armed groups allocation strategy for V."""
        if strategy == "recruitment":
            # Focus on recruitment and expansion
            self.V_weights = {"r15": 0.1, "r16": 0.1, "r25": 0.2, "r32": 0.6}
        elif strategy == "displacement":
            # Focus on displacing population
            self.V_weights = {"r15": 0.4, "r16": 0.4, "r25": 0.1, "r32": 0.1}
        elif strategy == "balanced":
            # Balanced approach
            self.V_weights = {"r15": 0.25, "r16": 0.25, "r25": 0.25, "r32": 0.25}
        else:
            raise ValueError(f"Unknown armed groups strategy: {strategy}")

        self.ag_strategy = strategy
        if not quiet:
            print(f"Armed groups strategy set to: {strategy}")

    def compute_allocations(self, state: Dict[str, float]) -> Tuple[Dict[str, float], Dict[str, float], Dict[str, float]]:
        """Compute allocations for E, Gov, and V based on current stocks and strategies."""
        # Get total available resources
        total_E = state.get("E", 0)
        total_Gov = state.get("Gov", 0)
        total_V = state.get("V", 0)

        # Compute allocations
        E_alloc = {rxn: weight * total_E for rxn, weight in self.E_weights.items()}
        Gov_alloc = {rxn: weight * total_Gov for rxn, weight in self.Gov_weights.items()}
        V_alloc = {rxn: weight * total_V for rxn, weight in self.V_weights.items()}

        return E_alloc, Gov_alloc, V_alloc

    def compute_reaction_rate(self, reaction: str, state: Dict[str, float],
                            E_alloc: Dict[str, float], Gov_alloc: Dict[str, float],
                            V_alloc: Dict[str, float]) -> Tuple[float, Dict[str, float]]:
        """Compute reaction rate using mass action with allocated budgets.

        Returns:
            rate: The reaction rate
            consumption: Dict with E, Gov, V consumption for this reaction
        """
        k = self.k[reaction]
        rate = k
        consumption = {"E": 0.0, "Gov": 0.0, "V": 0.0}

        # Determine which reactants this reaction uses
        if reaction == "r1":  # SL -> RL
            rate *= state.get("SL", 0)
        elif reaction == "r2":  # SR_SL -> SR_RL
            rate *= state.get("SR_SL", 0)
        elif reaction == "r3":  # WR_SL -> WR_RL
            rate *= state.get("WR_SL", 0)
        elif reaction == "r4":  # RL -> SL
            rate *= state.get("RL", 0)
        elif reaction == "r5":  # SR_RL -> SR_SL
            rate *= state.get("SR_RL", 0)
        elif reaction == "r6":  # WR_RL -> WR_SL
            rate *= state.get("WR_RL", 0)
        elif reaction == "r7":  # AG_RL -> AG_SL
            rate *= state.get("AG_RL", 0)
        elif reaction == "r8":  # Gov + SR_SL -> SR_RL
            gov_used = Gov_alloc.get("r8", 0)
            rate *= gov_used * state.get("SR_SL", 0)
            consumption["Gov"] = rate  # Consume Gov proportional to rate
        elif reaction == "r9":  # SR_RL -> SR_SL + E (produces E, handled by stoich)
            rate *= state.get("SR_RL", 0)
        elif reaction == "r10":  # WR_RL -> WR_SL + E (produces E, handled by stoich)
            rate *= state.get("WR_RL", 0)
        elif reaction == "r11":  # E -> Gov
            e_used = E_alloc.get("r11", 0)
            rate *= e_used
            consumption["E"] = rate  # Consume E proportional to rate
        elif reaction == "r12":  # -> Gov (produces Gov, no consumption)
            pass  # Constant rate
        elif reaction == "r13":  # SR_SL + RL -> SR_RL + SL
            rate *= state.get("SR_SL", 0) * state.get("RL", 0)
        elif reaction == "r14":  # WR_SL + RL -> WR_RL + SL
            rate *= state.get("WR_SL", 0) * state.get("RL", 0)
        elif reaction == "r15":  # SR_RL + SL + V -> RL + SR_SL
            v_used = V_alloc.get("r15", 0)
            rate *= state.get("SR_RL", 0) * state.get("SL", 0) * v_used
            consumption["V"] = rate
        elif reaction == "r16":  # WR_RL + SL + V -> RL + WR_SL
            v_used = V_alloc.get("r16", 0)
            rate *= state.get("WR_RL", 0) * state.get("SL", 0) * v_used
            consumption["V"] = rate
        elif reaction == "r17":  # AG_SL + RL -> AG_RL + SL
            rate *= state.get("AG_SL", 0) * state.get("RL", 0)
        elif reaction == "r18":  # E + WR_RL -> SR_RL
            e_used = E_alloc.get("r18", 0)
            rate *= e_used * state.get("WR_RL", 0)
            consumption["E"] = rate
        elif reaction == "r19":  # E + WR_SL -> SR_SL
            e_used = E_alloc.get("r19", 0)
            rate *= e_used * state.get("WR_SL", 0)
            consumption["E"] = rate
        elif reaction == "r20":  # SR_SL -> WR_SL
            rate *= state.get("SR_SL", 0)
        elif reaction == "r21":  # E + Gov + AG_SL -> WR_SL (reintegration)
            # FIX 2: Use geometric mean for multi-resource reaction
            e_alloc = E_alloc.get("r21", 0)
            gov_alloc = Gov_alloc.get("r21", 0)
            if e_alloc > 0 and gov_alloc > 0:
                combined_budget = sqrt(e_alloc * gov_alloc)
            else:
                combined_budget = 0
            rate *= combined_budget * state.get("AG_SL", 0)
            # Both E and Gov are consumed proportionally
            if combined_budget > 0:
                consumption["E"] = rate * sqrt(e_alloc / gov_alloc) if gov_alloc > 0 else rate
                consumption["Gov"] = rate * sqrt(gov_alloc / e_alloc) if e_alloc > 0 else rate
        elif reaction == "r22":  # SR_RL + E -> SR_RL + T (E is catalyst, not consumed)
            rate *= state.get("SR_RL", 0) * E_alloc.get("r22", 0)
            # E is catalyst here, not consumed (produces T)
        elif reaction == "r23":  # Gov + SR_RL -> Gov + SR_RL + T (Gov is catalyst)
            rate *= Gov_alloc.get("r23", 0) * state.get("SR_RL", 0)
            # Gov is catalyst here, not consumed (produces T)
        elif reaction == "r24":  # Gov + WR_RL -> SR_RL
            gov_used = Gov_alloc.get("r24", 0)
            rate *= gov_used * state.get("WR_RL", 0)
            consumption["Gov"] = rate
        elif reaction == "r25":  # V + T -> (destroys both)
            v_used = V_alloc.get("r25", 0)
            rate *= v_used * state.get("T", 0)
            consumption["V"] = rate
        elif reaction == "r26":  # Gov -> (corruption decay, handled by stoich)
            rate *= state.get("Gov", 0)
        elif reaction == "r27":  # WR_SL -> WR_SL + V (produces V)
            rate *= state.get("WR_SL", 0)
        elif reaction == "r28":  # AG_SL + SR_RL -> AG_RL + WR_SL + V
            rate *= state.get("AG_SL", 0) * state.get("SR_RL", 0)
        elif reaction == "r29":  # AG_SL + Gov -> WR_SL + V
            gov_used = Gov_alloc.get("r29", 0)
            rate *= state.get("AG_SL", 0) * gov_used
            consumption["Gov"] = rate
        elif reaction == "r30":  # AG_RL + Gov -> AG_SL + V
            gov_used = Gov_alloc.get("r30", 0)
            rate *= state.get("AG_RL", 0) * gov_used
            consumption["Gov"] = rate
        elif reaction == "r31":  # V -> (violence decay, handled by stoich)
            rate *= state.get("V", 0)
        elif reaction == "r32":  # WR_SL + AG_SL + V -> 2AG_SL
            v_used = V_alloc.get("r32", 0)
            rate *= state.get("WR_SL", 0) * state.get("AG_SL", 0) * v_used
            consumption["V"] = rate
        elif reaction == "r33":  # AG_SL + WR_RL -> AG_RL + WR_SL
            rate *= state.get("AG_SL", 0) * state.get("WR_RL", 0)
        elif reaction == "r34":  # AG_SL + SR_RL -> AG_RL + WR_SL
            rate *= state.get("AG_SL", 0) * state.get("SR_RL", 0)
        elif reaction == "r35":  # AG_SL + E -> AG_RL (uses total E, corruption)
            rate *= state.get("AG_SL", 0) * state.get("E", 0)
            consumption["E"] = rate
        elif reaction == "r36":  # AG_RL + Gov -> AG_RL (uses total Gov, capture)
            rate *= state.get("AG_RL", 0) * state.get("Gov", 0)
            consumption["Gov"] = rate

        return max(0.0, rate), consumption

    def compute_process_vector(self, state: Dict[str, float], verbose: bool = False) -> Tuple[np.ndarray, Dict[str, Any], Dict[str, float]]:
        """Compute process vector with current allocations.

        Returns:
            process_vector: Rates for all reactions
            allocation_info: Info about allocations for display
            total_consumption: Total consumption of E, Gov, V
        """
        # Compute allocations
        E_alloc, Gov_alloc, V_alloc = self.compute_allocations(state)

        # Compute rates for all reactions
        process_vector = np.zeros(len(self.REACTIONS))
        rates = {}
        total_consumption = {"E": 0.0, "Gov": 0.0, "V": 0.0}
        consumption_by_rxn = {}

        for i, reaction in enumerate(self.REACTIONS):
            rate, consumption = self.compute_reaction_rate(reaction, state, E_alloc, Gov_alloc, V_alloc)
            process_vector[i] = rate
            rates[reaction] = rate
            consumption_by_rxn[reaction] = consumption

            # Accumulate total consumption
            for resource in ["E", "Gov", "V"]:
                total_consumption[resource] += consumption[resource]

        # Store allocation info for display
        allocation_info = {
            "E": {"total": state.get("E", 0), "allocations": E_alloc},
            "Gov": {"total": state.get("Gov", 0), "allocations": Gov_alloc},
            "V": {"total": state.get("V", 0), "allocations": V_alloc},
            "consumption_by_rxn": consumption_by_rxn,
        }

        if verbose:
            print("\n" + "="*60)
            print("ALLOCATIONS:")
            print("="*60)

            print(f"\nECONOMIC BUDGET (E): {state.get('E', 0):.1f} units")
            for rxn, alloc in E_alloc.items():
                desc = self.reaction_descriptions.get(rxn, "")
                print(f"  {rxn}: {alloc:.2f} units ({desc})")

            print(f"\nGOVERNANCE BUDGET (Gov): {state.get('Gov', 0):.1f} units")
            for rxn, alloc in Gov_alloc.items():
                desc = self.reaction_descriptions.get(rxn, "")
                print(f"  {rxn}: {alloc:.2f} units ({desc})")

            print(f"\nVIOLENCE BUDGET (V): {state.get('V', 0):.1f} units")
            for rxn, alloc in V_alloc.items():
                desc = self.reaction_descriptions.get(rxn, "")
                print(f"  {rxn}: {alloc:.2f} units ({desc})")

            print("\n" + "="*60)
            print("REACTION RATES:")
            print("="*60)

            # Group reactions by category
            categories = {
                "Land & Resources (r1-r8)": [f"r{i}" for i in range(1, 9)],
                "Economic & Gov (r9-r12)": [f"r{i}" for i in range(9, 13)],
                "Migration (r13-r17)": [f"r{i}" for i in range(13, 18)],
                "Resilience (r18-r21)": [f"r{i}" for i in range(18, 22)],
                "Trust (r22-r26)": [f"r{i}" for i in range(22, 27)],
                "Violence (r27-r31)": [f"r{i}" for i in range(27, 32)],
                "Armed Groups (r32-r36)": [f"r{i}" for i in range(32, 37)],
            }

            for category, reactions in categories.items():
                print(f"\n{category}:")
                for rxn in reactions:
                    rate = rates.get(rxn, 0)
                    if rate > 1e-6:  # Only show non-negligible rates
                        desc = self.reaction_descriptions.get(rxn, "")
                        cons = consumption_by_rxn.get(rxn, {})
                        cons_str = ""
                        if any(v > 0 for v in cons.values()):
                            cons_parts = [f"{k}:{v:.2f}" for k, v in cons.items() if v > 0]
                            cons_str = f" [consumes {', '.join(cons_parts)}]"
                        print(f"  {rxn}: {rate:.4f} ({desc}){cons_str}")

            print(f"\nTOTAL CONSUMPTION THIS STEP:")
            for resource, amount in total_consumption.items():
                print(f"  {resource}: {amount:.4f}")

        return process_vector, allocation_info, total_consumption

    def apply_timestep(self, dt: float = 0.1, verbose: bool = True):
        """Apply one timestep of the simulation.

        FIX 1: Manual consumption tracking for E, Gov, V.
        FIX 4: Rate limiter to prevent over-consumption (conservation).
        """
        # Store current state for history
        old_state = self.state.copy()

        # Compute process vector with allocations and consumption
        process_vector, allocation_info, total_consumption = self.compute_process_vector(self.state, verbose)

        # FIX 4: Apply rate limiter - scale down all rates if any species would go negative
        # This preserves relative rates while ensuring conservation
        delta_raw = self.S @ process_vector * dt
        scale_factor = 1.0
        for i, species in enumerate(self.SPECIES):
            current = self.state.get(species, 0)
            if delta_raw[i] < 0 and current > 0:  # Species being consumed
                max_consumption = current * 0.5  # Don't consume more than 50% per step
                if abs(delta_raw[i]) > max_consumption:
                    species_scale = max_consumption / abs(delta_raw[i])
                    scale_factor = min(scale_factor, species_scale)

        if scale_factor < 1.0:
            if verbose:
                print(f"  [Rate limiter applied: scale={scale_factor:.3f}]")
            process_vector = process_vector * scale_factor
            total_consumption = {k: v * scale_factor for k, v in total_consumption.items()}

        # Compute state changes via stoichiometry for non-budget species
        delta = self.S @ process_vector * dt

        # Update state
        new_state = {}
        for i, species in enumerate(self.SPECIES):
            new_value = self.state.get(species, 0) + delta[i]
            new_state[species] = max(0.0, new_value)

        # FIX 1: Override E, Gov, V with manual consumption tracking
        # The stoichiometry handles production (r9, r10 produce E; r11, r12 produce Gov; r27 produces V)
        # But we need to subtract our manually tracked consumption

        # For E: stoichiometry adds production, we subtract manual consumption
        # Note: Some reactions have E in stoichiometry as consumed - we zero those out
        # and use our manual tracking instead
        e_production = 0.0
        for i, rxn in enumerate(self.REACTIONS):
            if rxn in ["r9", "r10"]:  # E-producing reactions
                e_production += process_vector[i] * dt
        new_state["E"] = max(0.0, old_state["E"] + e_production - total_consumption["E"] * dt)

        # For Gov: similar approach
        gov_production = 0.0
        for i, rxn in enumerate(self.REACTIONS):
            if rxn in ["r11", "r12"]:  # Gov-producing reactions
                gov_production += process_vector[i] * dt
        # Also subtract corruption decay (r26)
        r26_idx = self.REACTIONS.index("r26")
        gov_decay = process_vector[r26_idx] * dt
        new_state["Gov"] = max(0.0, old_state["Gov"] + gov_production - gov_decay - total_consumption["Gov"] * dt)

        # For V: similar approach
        v_production = 0.0
        for i, rxn in enumerate(self.REACTIONS):
            if rxn in ["r27", "r28", "r29", "r30"]:  # V-producing reactions
                v_production += process_vector[i] * dt
        # Also subtract natural decay (r31)
        r31_idx = self.REACTIONS.index("r31")
        v_decay = process_vector[r31_idx] * dt
        new_state["V"] = max(0.0, old_state["V"] + v_production - v_decay - total_consumption["V"] * dt)

        # Update time and state
        self.time += dt
        self.state = new_state

        # Store history
        self.time_history.append(self.time)
        self.state_history.append(self.state.copy())
        self.rate_history.append({rxn: process_vector[i] for i, rxn in enumerate(self.REACTIONS)})
        self.allocation_history.append(allocation_info)
        self.consumption_history.append(total_consumption)

        if verbose:
            print(f"\nTIME: {self.time:.1f} years")
            print("="*60)
            print("STATE CHANGES:")
            print("="*60)

            print(f"{'Species':<8} {'Old':>8} {'New':>8} {'Delta':>8}")
            print("-" * 36)

            for species in self.SPECIES:
                old_val = old_state.get(species, 0)
                new_val = self.state.get(species, 0)
                delta_val = new_val - old_val
                print(f"{species:<8} {old_val:>8.2f} {new_val:>8.2f} {delta_val:>+8.2f}")

            # Show budget changes with consumption breakdown
            print(f"\nBUDGET TRACKING (manual):")
            for resource in ["E", "Gov", "V"]:
                old = old_state.get(resource, 0)
                new = self.state.get(resource, 0)
                consumed = total_consumption[resource] * dt
                print(f"  {resource}: {old:.2f} -> {new:.2f} (consumed: {consumed:.2f})")

            # CONSERVATION CHECKS
            print(f"\nCONSERVATION CHECKS:")
            # Land conservation: RL + SL = const
            old_land = old_state.get("RL", 0) + old_state.get("SL", 0)
            new_land = self.state.get("RL", 0) + self.state.get("SL", 0)
            print(f"  Land (RL+SL): {old_land:.2f} -> {new_land:.2f} (delta: {new_land-old_land:+.4f})")

            # Population check (SR + WR + AG should be roughly conserved,
            # though AG can recruit from WR so it's not strict)
            old_pop = (old_state.get("SR_RL", 0) + old_state.get("SR_SL", 0) +
                       old_state.get("WR_RL", 0) + old_state.get("WR_SL", 0) +
                       old_state.get("AG_RL", 0) + old_state.get("AG_SL", 0))
            new_pop = (self.state.get("SR_RL", 0) + self.state.get("SR_SL", 0) +
                       self.state.get("WR_RL", 0) + self.state.get("WR_SL", 0) +
                       self.state.get("AG_RL", 0) + self.state.get("AG_SL", 0))
            print(f"  Total people (SR+WR+AG): {old_pop:.2f} -> {new_pop:.2f} (delta: {new_pop-old_pop:+.4f})")

    def run_simulation(self, total_time: float = 10.0, dt: float = 0.1, verbose: bool = True):
        """Run simulation for specified time."""
        steps = int(total_time / dt)

        if verbose:
            print(f"Running simulation for {total_time} years ({steps} steps)")
            print(f"Initial state:")
            for species, value in self.state_history[0].items():
                print(f"  {species}: {value:.1f}")

        for step in range(steps):
            if verbose and step == 0:
                print(f"\nStep {step+1}/{steps}:")
                self.apply_timestep(dt, verbose=True)
            else:
                self.apply_timestep(dt, verbose=False)

        if verbose:
            print(f"\nSimulation complete. Final state at t={self.time:.1f} years:")
            for species, value in self.state.items():
                print(f"  {species}: {value:.1f}")

    def plot_results(self):
        """Plot time series of species and reaction rates."""
        fig, axes = plt.subplots(4, 3, figsize=(15, 12))
        axes = axes.flatten()

        # Plot species concentrations
        print("\nPlotting species concentrations...")
        for idx, species in enumerate(self.SPECIES):
            if idx < len(axes):
                values = [state[species] for state in self.state_history]
                axes[idx].plot(self.time_history, values, 'b-', linewidth=2)
                axes[idx].set_title(species)
                axes[idx].set_xlabel("Time (years)")
                axes[idx].set_ylabel("Concentration")
                axes[idx].grid(True, alpha=0.3)

        plt.suptitle("Species Concentrations Over Time", fontsize=16)
        plt.tight_layout()
        plt.show()

        # Plot reaction rates for selected reactions
        fig2, axes2 = plt.subplots(3, 4, figsize=(16, 10))
        axes2 = axes2.flatten()

        # Select key reactions to plot
        key_reactions = ["r8", "r11", "r15", "r16", "r21", "r23", "r25", "r28", "r29", "r30", "r32", "r35"]

        print("\nPlotting key reaction rates...")
        for idx, reaction in enumerate(key_reactions):
            if idx < len(axes2):
                rates = [rates_dict.get(reaction, 0) for rates_dict in self.rate_history]
                # Pad with initial rate of 0
                if len(rates) < len(self.time_history):
                    rates = [0] + rates

                axes2[idx].plot(self.time_history[:len(rates)], rates, 'r-', linewidth=2)
                axes2[idx].set_title(f"{reaction}: {self.reaction_descriptions.get(reaction, '')[:30]}...")
                axes2[idx].set_xlabel("Time (years)")
                axes2[idx].set_ylabel("Rate")
                axes2[idx].grid(True, alpha=0.3)

        plt.suptitle("Key Reaction Rates Over Time", fontsize=16)
        plt.tight_layout()
        plt.show()

        # Plot budget allocations over time
        fig3, axes3 = plt.subplots(1, 3, figsize=(15, 5))

        resources = ["E", "Gov", "V"]
        titles = ["Economic Budget Allocation", "Governance Budget Allocation", "Violence Budget Allocation"]

        for i, (resource, title) in enumerate(zip(resources, titles)):
            # Get unique reactions for this resource
            reactions = []
            if resource == "E":
                reactions = list(self.E_weights.keys())
            elif resource == "Gov":
                reactions = list(self.Gov_weights.keys())
            elif resource == "V":
                reactions = list(self.V_weights.keys())

            # Plot allocations over time
            for rxn in reactions:
                allocations = []
                for alloc_dict in self.allocation_history:
                    alloc = alloc_dict[resource]["allocations"].get(rxn, 0)
                    allocations.append(alloc)

                # Pad with initial allocation
                if len(allocations) < len(self.time_history):
                    allocations = [0] + allocations

                axes3[i].plot(self.time_history[:len(allocations)], allocations, label=f"{rxn}", linewidth=2)

            axes3[i].set_title(title)
            axes3[i].set_xlabel("Time (years)")
            axes3[i].set_ylabel("Allocated Amount")
            axes3[i].legend(loc='best')
            axes3[i].grid(True, alpha=0.3)

        plt.suptitle("Budget Allocations Over Time", fontsize=16)
        plt.tight_layout()
        plt.show()

# ============================================================================
# BATCH SIMULATION FUNCTIONS
# ============================================================================

def run_batch_simulation(model_file: str, total_time: float = 10.0, dt: float = 0.1):
    """Run comprehensive batch simulation across all parameter combinations.

    Returns dict with results for each combination:
    results[scenario][climate_k][gov_strategy][ag_strategy] = metrics
    """
    print("="*70)
    print("BATCH CONFLICT DYNAMICS SIMULATION")
    print("="*70)
    print(f"\nSimulation parameters:")
    print(f"  Total time: {total_time} years")
    print(f"  Timestep: {dt} years")
    print(f"  Conflict scenarios: {list(CONFLICT_SCENARIOS.keys())}")
    print(f"  Climate parameters: {CLIMATE_PARAMETERS}")
    print(f"  Government strategies: {GOVERNMENT_STRATEGIES}")
    print(f"  Armed group strategies: {ARMED_GROUP_STRATEGIES}")

    total_runs = (len(CONFLICT_SCENARIOS) * len(CLIMATE_PARAMETERS) *
                  len(GOVERNMENT_STRATEGIES) * len(ARMED_GROUP_STRATEGIES))
    print(f"\nTotal simulation runs: {total_runs}")
    print("="*70)

    # Initialize simulator once
    simulator = ConflictSimulator(model_file)

    # Store results
    results = {}
    run_count = 0

    for scenario in CONFLICT_SCENARIOS:
        results[scenario] = {}
        print(f"\n>>> Scenario: {scenario.upper()} - {CONFLICT_SCENARIOS[scenario]['description']}")

        for climate_k in CLIMATE_PARAMETERS:
            results[scenario][climate_k] = {}

            for gov_strategy in GOVERNMENT_STRATEGIES:
                results[scenario][climate_k][gov_strategy] = {}

                for ag_strategy in ARMED_GROUP_STRATEGIES:
                    run_count += 1
                    print(f"  Run {run_count}/{total_runs}: climate={climate_k}, "
                          f"gov={gov_strategy}, ag={ag_strategy}", end="")

                    # Configure simulation
                    simulator.set_conflict_scenario(scenario)
                    simulator.set_climate_parameter(climate_k)
                    simulator.set_government_strategy(gov_strategy, quiet=True)
                    simulator.set_armed_groups_strategy(ag_strategy, quiet=True)

                    # Run simulation
                    simulator.run_simulation(total_time=total_time, dt=dt, verbose=False)

                    # Extract metrics
                    metrics = simulator.get_metrics()
                    metrics["scenario"] = scenario
                    metrics["climate_k"] = climate_k
                    metrics["gov_strategy"] = gov_strategy
                    metrics["ag_strategy"] = ag_strategy

                    results[scenario][climate_k][gov_strategy][ag_strategy] = metrics
                    print(f" -> AG_final={metrics['final_AG']:.1f}, E_final={metrics['final_E']:.1f}")

    print("\n" + "="*70)
    print("BATCH SIMULATION COMPLETE")
    print("="*70)

    return results


def plot_comparison_results(results: Dict, save_dir: str = None):
    """Generate comprehensive comparison plots from batch simulation results."""

    if save_dir:
        os.makedirs(save_dir, exist_ok=True)

    # ========================================================================
    # PLOT 1: Budget Spent by Strategy (for each scenario and climate)
    # ========================================================================
    print("\nGenerating Plot 1: Budget Spent Comparison...")

    fig, axes = plt.subplots(3, 4, figsize=(20, 15))
    fig.suptitle("Total Budget Spent (E + Gov) by Strategy Combination", fontsize=16)

    for row, scenario in enumerate(CONFLICT_SCENARIOS.keys()):
        for col, climate_k in enumerate(CLIMATE_PARAMETERS):
            ax = axes[row, col]

            # Collect data for this scenario/climate
            x_labels = []
            e_consumed = []
            gov_consumed = []

            for gov_strat in GOVERNMENT_STRATEGIES:
                for ag_strat in ARMED_GROUP_STRATEGIES:
                    metrics = results[scenario][climate_k][gov_strat][ag_strat]
                    x_labels.append(f"{gov_strat[:3]}-{ag_strat[:3]}")
                    e_consumed.append(metrics["total_E_consumed"])
                    gov_consumed.append(metrics["total_Gov_consumed"])

            x = np.arange(len(x_labels))
            width = 0.35

            ax.bar(x - width/2, e_consumed, width, label='E consumed', color='gold')
            ax.bar(x + width/2, gov_consumed, width, label='Gov consumed', color='steelblue')

            ax.set_title(f"{scenario.title()} | climate={climate_k}")
            ax.set_xticks(x)
            ax.set_xticklabels(x_labels, rotation=45, ha='right', fontsize=8)
            ax.set_ylabel("Total Consumed")
            if row == 0 and col == 0:
                ax.legend(loc='upper right', fontsize=8)

    plt.tight_layout()
    if save_dir:
        plt.savefig(os.path.join(save_dir, "01_budget_spent.png"), dpi=150)
    plt.show()

    # ========================================================================
    # PLOT 2: Prosperity (E and T) Final Values
    # ========================================================================
    print("Generating Plot 2: Prosperity (E and T) Comparison...")

    fig, axes = plt.subplots(3, 4, figsize=(20, 15))
    fig.suptitle("Final Prosperity: Economy (E) and Trust (T)", fontsize=16)

    for row, scenario in enumerate(CONFLICT_SCENARIOS.keys()):
        for col, climate_k in enumerate(CLIMATE_PARAMETERS):
            ax = axes[row, col]

            x_labels = []
            final_e = []
            final_t = []

            for gov_strat in GOVERNMENT_STRATEGIES:
                for ag_strat in ARMED_GROUP_STRATEGIES:
                    metrics = results[scenario][climate_k][gov_strat][ag_strat]
                    x_labels.append(f"{gov_strat[:3]}-{ag_strat[:3]}")
                    final_e.append(metrics["final_E"])
                    final_t.append(metrics["final_T"])

            x = np.arange(len(x_labels))
            width = 0.35

            ax.bar(x - width/2, final_e, width, label='E (Economy)', color='green')
            ax.bar(x + width/2, final_t, width, label='T (Trust)', color='teal')

            ax.set_title(f"{scenario.title()} | climate={climate_k}")
            ax.set_xticks(x)
            ax.set_xticklabels(x_labels, rotation=45, ha='right', fontsize=8)
            ax.set_ylabel("Final Value")
            if row == 0 and col == 0:
                ax.legend(loc='upper right', fontsize=8)

    plt.tight_layout()
    if save_dir:
        plt.savefig(os.path.join(save_dir, "02_prosperity.png"), dpi=150)
    plt.show()

    # ========================================================================
    # PLOT 3: Governance Final Values
    # ========================================================================
    print("Generating Plot 3: Governance Comparison...")

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle("Final Governance by Scenario and Strategy", fontsize=16)

    for idx, scenario in enumerate(CONFLICT_SCENARIOS.keys()):
        ax = axes[idx]

        # Group by climate parameter
        climate_groups = {k: [] for k in CLIMATE_PARAMETERS}
        labels = []

        for gov_strat in GOVERNMENT_STRATEGIES:
            for ag_strat in ARMED_GROUP_STRATEGIES:
                label = f"{gov_strat[:3]}-{ag_strat[:3]}"
                if label not in labels:
                    labels.append(label)

                for climate_k in CLIMATE_PARAMETERS:
                    metrics = results[scenario][climate_k][gov_strat][ag_strat]
                    climate_groups[climate_k].append(metrics["final_Gov"])

        x = np.arange(len(labels))
        width = 0.2

        for i, climate_k in enumerate(CLIMATE_PARAMETERS):
            offset = (i - len(CLIMATE_PARAMETERS)/2 + 0.5) * width
            ax.bar(x + offset, climate_groups[climate_k], width,
                   label=f'climate={climate_k}', alpha=0.8)

        ax.set_title(f"{scenario.title()} Conflict")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=9)
        ax.set_ylabel("Final Governance")
        ax.legend(loc='upper right', fontsize=8)

    plt.tight_layout()
    if save_dir:
        plt.savefig(os.path.join(save_dir, "03_governance.png"), dpi=150)
    plt.show()

    # ========================================================================
    # PLOT 4: Strong Resilience (SR_total)
    # ========================================================================
    print("Generating Plot 4: Community Strong Resilience Comparison...")

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle("Final Strong Resilient Population (SR_RL + SR_SL)", fontsize=16)

    for idx, scenario in enumerate(CONFLICT_SCENARIOS.keys()):
        ax = axes[idx]

        climate_groups = {k: [] for k in CLIMATE_PARAMETERS}
        labels = []

        for gov_strat in GOVERNMENT_STRATEGIES:
            for ag_strat in ARMED_GROUP_STRATEGIES:
                label = f"{gov_strat[:3]}-{ag_strat[:3]}"
                if label not in labels:
                    labels.append(label)

                for climate_k in CLIMATE_PARAMETERS:
                    metrics = results[scenario][climate_k][gov_strat][ag_strat]
                    climate_groups[climate_k].append(metrics["final_SR"])

        x = np.arange(len(labels))
        width = 0.2

        for i, climate_k in enumerate(CLIMATE_PARAMETERS):
            offset = (i - len(CLIMATE_PARAMETERS)/2 + 0.5) * width
            ax.bar(x + offset, climate_groups[climate_k], width,
                   label=f'climate={climate_k}', alpha=0.8)

        ax.set_title(f"{scenario.title()} Conflict")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=9)
        ax.set_ylabel("Final SR Population")
        ax.legend(loc='upper right', fontsize=8)

    plt.tight_layout()
    if save_dir:
        plt.savefig(os.path.join(save_dir, "04_strong_resilience.png"), dpi=150)
    plt.show()

    # ========================================================================
    # PLOT 5: Armed Group Sustainability (AG_total)
    # ========================================================================
    print("Generating Plot 5: Armed Group Sustainability Comparison...")

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle("Final Armed Group Population (AG_RL + AG_SL)", fontsize=16)

    for idx, scenario in enumerate(CONFLICT_SCENARIOS.keys()):
        ax = axes[idx]

        climate_groups = {k: [] for k in CLIMATE_PARAMETERS}
        labels = []

        for gov_strat in GOVERNMENT_STRATEGIES:
            for ag_strat in ARMED_GROUP_STRATEGIES:
                label = f"{gov_strat[:3]}-{ag_strat[:3]}"
                if label not in labels:
                    labels.append(label)

                for climate_k in CLIMATE_PARAMETERS:
                    metrics = results[scenario][climate_k][gov_strat][ag_strat]
                    climate_groups[climate_k].append(metrics["final_AG"])

        x = np.arange(len(labels))
        width = 0.2

        for i, climate_k in enumerate(CLIMATE_PARAMETERS):
            offset = (i - len(CLIMATE_PARAMETERS)/2 + 0.5) * width
            ax.bar(x + offset, climate_groups[climate_k], width,
                   label=f'climate={climate_k}', alpha=0.8)

        ax.set_title(f"{scenario.title()} Conflict")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=9)
        ax.set_ylabel("Final AG Population")
        ax.legend(loc='upper right', fontsize=8)

        # Add initial value reference line
        initial_ag = (CONFLICT_SCENARIOS[scenario]["initial_state"]["AG_RL"] +
                      CONFLICT_SCENARIOS[scenario]["initial_state"]["AG_SL"])
        ax.axhline(y=initial_ag, color='red', linestyle='--', linewidth=2,
                   label=f'Initial: {initial_ag:.1f}')

    plt.tight_layout()
    if save_dir:
        plt.savefig(os.path.join(save_dir, "05_armed_groups.png"), dpi=150)
    plt.show()

    # ========================================================================
    # PLOT 6: Time Series Comparison (selected scenarios)
    # ========================================================================
    print("Generating Plot 6: Time Series Comparison...")

    # Pick medium conflict, middle climate for detailed time series
    scenario = "medium"
    climate_k = 0.02  # Use actual climate parameter from CLIMATE_PARAMETERS

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle(f"Time Series: {scenario.title()} Conflict, Climate={climate_k}", fontsize=16)

    metrics_to_plot = [
        ("E", "Economy (E)", "green"),
        ("T", "Trust (T)", "teal"),
        ("Gov", "Governance", "blue"),
        ("SR_total", "Strong Resilient", "darkgreen"),
        ("AG_total", "Armed Groups", "red"),
        ("V", "Violence", "darkred"),
    ]

    for ax_idx, (metric_key, metric_label, _) in enumerate(metrics_to_plot):
        ax = axes.flatten()[ax_idx]

        for gov_strat in GOVERNMENT_STRATEGIES:
            for ag_strat in ARMED_GROUP_STRATEGIES:
                metrics = results[scenario][climate_k][gov_strat][ag_strat]
                label = f"{gov_strat[:3]}-{ag_strat[:3]}"

                linestyle = '-' if gov_strat == "development" else ('--' if gov_strat == "security" else ':')
                marker = 'o' if ag_strat == "recruitment" else ('s' if ag_strat == "displacement" else '^')

                ax.plot(metrics["times"], metrics[metric_key],
                        label=label, linestyle=linestyle, marker=marker,
                        markevery=20, markersize=4, alpha=0.7)

        ax.set_title(metric_label)
        ax.set_xlabel("Time (years)")
        ax.set_ylabel("Value")
        ax.grid(True, alpha=0.3)
        if ax_idx == 0:
            ax.legend(loc='best', fontsize=7, ncol=3)

    plt.tight_layout()
    if save_dir:
        plt.savefig(os.path.join(save_dir, "06_time_series.png"), dpi=150)
    plt.show()

    # ========================================================================
    # PLOT 7: Heatmap Summary
    # ========================================================================
    print("Generating Plot 7: Strategy Effectiveness Heatmap...")

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.suptitle("Strategy Effectiveness: Final AG Population (lower is better for government)", fontsize=14)

    # First pass: compute all matrices and find global min/max for consistent color scale
    matrices = {}
    global_min = float('inf')
    global_max = float('-inf')

    for scenario in CONFLICT_SCENARIOS.keys():
        matrix = np.zeros((len(GOVERNMENT_STRATEGIES), len(ARMED_GROUP_STRATEGIES)))
        for i, gov_strat in enumerate(GOVERNMENT_STRATEGIES):
            for j, ag_strat in enumerate(ARMED_GROUP_STRATEGIES):
                values = [results[scenario][ck][gov_strat][ag_strat]["final_AG"]
                          for ck in CLIMATE_PARAMETERS]
                matrix[i, j] = np.mean(values)
        matrices[scenario] = matrix
        global_min = min(global_min, matrix.min())
        global_max = max(global_max, matrix.max())

    # Second pass: plot with consistent color scale
    for idx, scenario in enumerate(CONFLICT_SCENARIOS.keys()):
        ax = axes[idx]
        matrix = matrices[scenario]

        im = ax.imshow(matrix, cmap='RdYlGn_r', aspect='auto',
                       vmin=global_min, vmax=global_max)
        ax.set_xticks(range(len(ARMED_GROUP_STRATEGIES)))
        ax.set_yticks(range(len(GOVERNMENT_STRATEGIES)))
        ax.set_xticklabels([s[:4] for s in ARMED_GROUP_STRATEGIES])
        ax.set_yticklabels([s[:4] for s in GOVERNMENT_STRATEGIES])
        ax.set_xlabel("Armed Group Strategy")
        ax.set_ylabel("Government Strategy")
        ax.set_title(f"{scenario.title()} Conflict")

        # Add value annotations
        for i in range(len(GOVERNMENT_STRATEGIES)):
            for j in range(len(ARMED_GROUP_STRATEGIES)):
                ax.text(j, i, f"{matrix[i, j]:.1f}",
                        ha="center", va="center", color="black", fontsize=10)

    # Single colorbar for all subplots
    fig.colorbar(im, ax=axes, label="Final AG Population", shrink=0.8)

    plt.tight_layout()
    if save_dir:
        plt.savefig(os.path.join(save_dir, "07_heatmap.png"), dpi=150)
    plt.show()

    print("\nAll plots generated!")


def main():
    """Main entry point for batch simulation."""
    print("="*70)
    print("CONFLICT DYNAMICS BATCH SIMULATOR")
    print("="*70)
    print("\nThis will run simulations across:")
    print(f"  - 3 conflict scenarios: {list(CONFLICT_SCENARIOS.keys())}")
    print(f"  - 4 climate parameters: {CLIMATE_PARAMETERS}")
    print(f"  - 3 government strategies: {GOVERNMENT_STRATEGIES}")
    print(f"  - 3 armed group strategies: {ARMED_GROUP_STRATEGIES}")
    print(f"\nTotal: {3*4*3*3} = 108 simulations, each for 10 years")

    # Load the model
    project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
    model_file = os.path.join(project_root, 'networks/Conflict_Theory/Resource_Community_Insurgency_Loops_model1.txt')

    # Run batch simulation
    results = run_batch_simulation(model_file, total_time=10.0, dt=0.1)

    # Generate plots
    save_dir = "visualizations/conflict_simulation_results"
    plot_comparison_results(results, save_dir=save_dir)

    # Print summary table
    print("\n" + "="*70)
    print("SUMMARY: Best Government Strategy by Scenario")
    print("="*70)

    for scenario in CONFLICT_SCENARIOS.keys():
        print(f"\n{scenario.upper()} CONFLICT:")
        print(f"  Initial AG: {CONFLICT_SCENARIOS[scenario]['initial_state']['AG_RL'] + CONFLICT_SCENARIOS[scenario]['initial_state']['AG_SL']:.1f}")

        best_result = None
        best_ag = float('inf')

        for climate_k in CLIMATE_PARAMETERS:
            for gov_strat in GOVERNMENT_STRATEGIES:
                for ag_strat in ARMED_GROUP_STRATEGIES:
                    metrics = results[scenario][climate_k][gov_strat][ag_strat]
                    if metrics["final_AG"] < best_ag:
                        best_ag = metrics["final_AG"]
                        best_result = (climate_k, gov_strat, ag_strat, metrics)

        if best_result:
            ck, gs, ags, m = best_result
            print(f"  Best outcome (lowest final AG):")
            print(f"    Climate: {ck}, Gov: {gs}, AG: {ags}")
            print(f"    Final AG: {m['final_AG']:.2f}, Final E: {m['final_E']:.2f}, Final Gov: {m['final_Gov']:.2f}")

    print("\n" + "="*70)
    print(f"Results saved to: {save_dir}/")
    print("="*70)


if __name__ == "__main__":
    main()
