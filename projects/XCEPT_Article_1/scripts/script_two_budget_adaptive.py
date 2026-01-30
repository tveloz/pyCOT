#!/usr/bin/env python3
"""
Conflict Dynamics Simulator with Adaptive Budget Allocation

Enhanced Features:
1. Manual consumption tracking for E, Gov, V (not via stoichiometry)
2. Geometric mean for multi-resource reactions (r21)
3. Properly scaled kinetic constants by reaction order
4. Rate limiter for conservation
5. Flexible climate parameters (up to 10)
6. Budget renewal rate parameter (r12)
7. Baseline values shown in comparison plots
8. Adaptive government strategy that optimizes multiple objectives
9. Budget generation and consumption tracking
"""

import os
import sys
import numpy as np
from math import sqrt
from typing import Dict, List, Tuple, Any, Optional
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend - save plots without displaying
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# Add the root directory to the PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))
from pyCOT.io.functions import read_txt
from pyCOT.semantic_partition import define_semantic_categories, SemanticPartition
from pyCOT.process_analyzer import (classify_process_mode, is_cognitive_domain,
                                     analyze_category_behavior, track_category_sequence,
                                     compute_category_statistics)
from pyCOT.process_structure import rolling_window_aggregation

# ============================================================================
# SCENARIO DEFINITIONS
# ============================================================================

CONFLICT_SCENARIOS = {
    "low": {
        "description": "Low Conflict - Stable region with minor tensions",
        "initial_state": {
            "SR_RL": 15.0,   # More strong resilient
            "SR_SL": 12.0,
            "WR_RL": 10.0,
            "WR_SL": 8.0,   # Less weak on stressed land
            "AG_RL": 0.5,    # Few armed groups
            "AG_SL": 1.5,
            "RL": 20.0,      # More restored land
            "SL": 100.0,
            "E": 100.0,      # Strong economy
            "T": 60.0,       # High trust
            "V": 10.0,       # Low violence
            "Gov": 30.0,     # Strong governance
        }
    },
    "medium": {
        "description": "Medium Conflict - Regional instability with active tensions",
        "initial_state": {
            "SR_RL": 12.0,
            "SR_SL": 8.0,
            "WR_RL": 15.0,
            "WR_SL": 10.0,
            "AG_RL": 0.7,
            "AG_SL": 1.7,
            "RL": 20.0,
            "SL": 100.0,
            "E": 80.0,
            "T": 40.0,
            "V": 30.0,
            "Gov": 20.0,
        }
    },
    "severe": {
        "description": "Severe Conflict - Active insurgency with widespread violence",
        "initial_state": {
            "SR_RL": 8.0,   # Few strong resilient
            "SR_SL": 12.0,
            "WR_RL": 10.0,
            "WR_SL": 15.0,   # Many weak on stressed land
            "AG_RL": 0.9,    # More armed groups
            "AG_SL": 1.9,
            "RL": 20.0,      # Less restored land
            "SL": 100.0,
            "E": 50.0,       # Weak economy
            "T": 20.0,       # Low trust
            "V": 50.0,       # High violence
            "Gov": 10.0,     # Weak governance
        }
    }
}

# Flexible parameters - modify these as needed (max 10 climate parameters)
CLIMATE_PARAMETERS = [0.01, 0.1, 1]
BUDGET_RENEWAL_RATES = [0.05, 0.4, 1]  # r12 governance baseline rate

GOVERNMENT_STRATEGIES = ["development", "security", "balanced", "adaptive"]
ARMED_GROUP_STRATEGIES = ["recruitment", "displacement", "balanced"]

class ConflictSimulator:
    """Simulator for conflict dynamics with adaptive budget allocation."""

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
            "r9": "SR_RL -> SR_SL + E (economic production)",
            "r10": "WR_RL -> WR_SL + E (economic production)",
            "r11": "E -> Gov (economy to governance)",
            "r12": "-> Gov (baseline governance)",
            "r13": "SR_SL + RL -> SR_RL + SL (migration strong)",
            "r14": "WR_SL + RL -> WR_RL + SL (migration weak)",
            "r15": "SR_RL + SL + V -> RL + SR_SL (violent displacement strong)",
            "r16": "WR_RL + SL + V -> RL + WR_SL (violent displacement weak)",
            "r17": "AG_SL + RL -> AG_RL + SL (AG migration)",
            "r18": "E + WR_RL -> SR_RL (economic resilience)",
            "r19": "E + WR_SL -> SR_SL (economic resilience)",
            "r20": "SR_SL -> WR_SL (resilience degradation)",
            "r21": "E + Gov + AG_SL -> WR_SL (reintegration)",
            "r22": "SR_RL + E -> SR_RL + T (trust from prosperity)",
            "r23": "Gov + SR_RL -> Gov + SR_RL + T (gov trust building)",
            "r24": "Gov + WR_RL -> SR_RL (gov outreach)",
            "r25": "V + T -> (trust destruction)",
            "r26": "Gov -> (gov corruption)",
            "r27": "WR_SL -> WR_SL + V (desperation violence)",
            "r28": "AG_SL + SR_RL -> AG_RL + WR_SL + V (AG attacks)",
            "r29": "AG_SL + Gov -> WR_SL + V (gov-AG conflict)",
            "r30": "AG_RL + Gov -> AG_SL + V (gov displacement ops)",
            "r31": "V -> (violence decay)",
            "r32": "WR_SL + AG_SL + V -> 2AG_SL (AG recruitment)",
            "r33": "AG_SL + WR_RL -> AG_RL + WR_SL (AG expansion weak)",
            "r34": "AG_SL + SR_RL -> AG_RL + WR_SL (AG expansion strong)",
            "r35": "AG_SL + E -> AG_RL (AG corruption)",
            "r36": "AG_RL + Gov -> AG_RL (AG institutional capture)",
        }

        # Define which reactions CONSUME allocated resources (for manual tracking)
        # Note: r25 (V+T->; trust destruction) and r31 (V->; violence decay) consume V
        # through stoichiometry, NOT manual tracking - they are pure mass action reactions
        self.E_consuming_reactions = ["r11", "r18", "r19", "r21", "r35"]
        self.Gov_consuming_reactions = ["r8", "r21", "r24", "r29", "r30", "r36"]
        self.V_consuming_reactions = ["r15", "r16", "r32"]  # Removed r25, r31 (stoichiometric only)

        # Define which reactions use ALLOCATED resources (strategy-dependent)
        self.E_allocated_reactions = ["r11", "r18", "r19", "r21", "r22"]
        self.Gov_allocated_reactions = ["r8", "r21", "r23", "r24", "r29", "r30"]
        self.V_allocated_reactions = ["r15", "r16", "r32"]  # Removed r25 (pure mass action)

        # Set default strategies
        self.set_government_strategy("balanced")
        self.set_armed_groups_strategy("balanced")

        # Set kinetic parameters
        self.set_kinetic_parameters()

        # Seasonal variation parameters
        self.enable_seasonal_variation = True  # Toggle for seasonal effects
        self.seasonal_amplitude = 0.3  # 30% variation (±30% from baseline)
        self.seasonal_period = 1.0  # One year period
        self.seasonal_phase = 0.0  # Phase shift (0 = max at t=0)

        # Initialize state
        self.state = self.get_initial_state()
        self.time = 0.0

        # History tracking
        self.time_history = [0.0]
        self.state_history = [self.state.copy()]
        self.rate_history = []
        self.allocation_history = []
        self.consumption_history = []
        self.generation_history = []

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
        self.generation_history = []

    def get_seasonal_modifier(self, reaction: str) -> float:
        """
        Compute seasonal modifier for resource recovery reactions.

        Returns a rate multiplier that varies sinusoidally over the year,
        representing seasonal variations in resource renewal/recovery.

        Formula: (A/2) * (1 + cos(2π * t / period + phase))

        This formula ranges from [0, A] where A is the seasonal amplitude:
        - When cos = 1: rate = A (maximum renewal)
        - When cos = -1: rate = 0 (no renewal)
        - Average over year: rate = A/2

        Args:
            reaction: Reaction name

        Returns:
            Seasonal rate multiplier in range [0, A]
        """
        if not self.enable_seasonal_variation:
            return self.seasonal_amplitude  # Return constant amplitude if disabled

        # Reactions affected by seasonal variation
        # All resource renewal/recovery reactions (r1-r7):
        # r1: SL => RL (land restoration)
        # r2: SR_SL => SR_RL (strong resilient population resource renewal)
        # r3: WR_SL => WR_RL (weak resilient population resource renewal)
        # r4: RL => SL (land degradation)
        # r5: SR_RL => SR_SL (strong resilient degradation)
        # r6: WR_RL => WR_SL (weak resilient degradation)
        # r7: AG_RL => AG_SL (armed groups degradation)
        seasonal_reactions = ["r1", "r2", "r3", "r4", "r5", "r6", "r7"]

        if reaction not in seasonal_reactions:
            return 1.0  # Non-seasonal reactions use constant rate

        # Compute sinusoidal seasonal factor
        # cos(2π*t/T + φ) varies from -1 to +1
        # (A/2) * (1 + cos(...)) varies from 0 to A
        import math
        seasonal_factor = (self.seasonal_amplitude / 2.0) * (
            1.0 + math.cos(2.0 * math.pi * self.time / self.seasonal_period + self.seasonal_phase)
        )

        return seasonal_factor

    def set_conflict_scenario(self, scenario: str):
        """Set initial state based on conflict scenario."""
        if scenario not in CONFLICT_SCENARIOS:
            raise ValueError(f"Unknown scenario: {scenario}. Choose from: {list(CONFLICT_SCENARIOS.keys())}")
        self.current_scenario = scenario
        self.reset(CONFLICT_SCENARIOS[scenario]["initial_state"])

    def set_climate_parameter(self, climate_amplitude: float):
        """
        Set seasonal amplitude for resource renewal reactions (r1-r7).

        The climate parameter now controls the amplitude A in the seasonal formula:
        rate = (A/2) * (1 + cos(2π*t/T))

        This formula produces rates that vary from 0 to A over the year, with
        an annual average of A/2.

        Args:
            climate_amplitude: The maximum rate for resource renewal (A in formula).
                             For example, A=0.04 means rates vary from 0 to 0.04.
        """
        self.climate_k = climate_amplitude  # Store for reference
        self.seasonal_amplitude = climate_amplitude  # Set the seasonal amplitude

    def set_budget_renewal_rate(self, renewal_rate: float):
        """Set budget renewal rate (r12 - baseline governance production)."""
        self.budget_renewal_rate = renewal_rate
        self.k["r12"] = renewal_rate

    def get_metrics(self) -> Dict[str, Any]:
        """Extract key metrics from simulation history."""
        if len(self.state_history) < 2:
            return {}

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

        # Budget tracking - both consumption and generation
        if self.consumption_history:
            E_consumed_per_step = [c["E"] * 0.1 for c in self.consumption_history]
            Gov_consumed_per_step = [c["Gov"] * 0.1 for c in self.consumption_history]
            E_consumed = np.cumsum(E_consumed_per_step)
            Gov_consumed = np.cumsum(Gov_consumed_per_step)
            E_consumed = np.concatenate([[0], E_consumed])
            Gov_consumed = np.concatenate([[0], Gov_consumed])
        else:
            E_consumed = np.zeros_like(times)
            Gov_consumed = np.zeros_like(times)

        if self.generation_history:
            E_generated_per_step = [g["E"] * 0.1 for g in self.generation_history]
            Gov_generated_per_step = [g["Gov"] * 0.1 for g in self.generation_history]
            E_generated = np.cumsum(E_generated_per_step)
            Gov_generated = np.cumsum(Gov_generated_per_step)
            E_generated = np.concatenate([[0], E_generated])
            Gov_generated = np.concatenate([[0], Gov_generated])
        else:
            E_generated = np.zeros_like(times)
            Gov_generated = np.zeros_like(times)

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
            "E_generated": E_generated,
            "Gov_generated": Gov_generated,
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
            "total_E_generated": E_generated[-1],
            "total_Gov_generated": Gov_generated[-1],
            "total_budget_consumed": E_consumed[-1] + Gov_consumed[-1],
            "total_budget_generated": E_generated[-1] + Gov_generated[-1],
            "net_budget": (E_generated[-1] - E_consumed[-1]) + (Gov_generated[-1] - Gov_consumed[-1]),
        }

    def set_kinetic_parameters(self):
        """
        Set kinetic parameters (per year timescale).

        Note: r1-r7 are resource renewal reactions controlled by seasonal variation.
        Their base rate is set to 1.0 and the actual rate is determined by the
        seasonal modifier which uses CLIMATE_PARAMETERS as amplitude.
        """
        self.k = {
            # Resource renewal/degradation reactions (r1-r7) - base rate = 1.0
            # Actual rate = seasonal_modifier * base_rate, where seasonal_modifier = (A/2)*(1+cos(...))
            "r1": 1.0, "r2": 1.0, "r3": 1.0, "r4": 1.0, "r5": 1.0,
            "r6": 1.0, "r7": 1.0,
            # Other reactions maintain their original kinetic constants
            "r8": 0.002, "r9": 0.05, "r10": 0.03,
            "r11": 0.02, "r12": 1.0, "r13": 0.001, "r14": 0.001, "r15": 0.0001,
            "r16": 0.0001, "r17": 0.005, "r18": 0.001, "r19": 0.001, "r20": 0.02,
            "r21": 0.001, "r22": 0.001, "r23": 0.002, "r24": 0.002, "r25": 0.002,
            "r26": 0.02, "r27": 0.01, "r28": 0.005, "r29": 0.01, "r30": 0.02,
            "r31": 0.1, "r32": 0.0002, "r33": 0.005, "r34": 0.005, "r35": 0.002,
            "r36": 0.005,
        }

    def get_initial_state(self) -> Dict[str, float]:
        """Return initial state for simulation."""
        return {
            "SR_RL": 30.0, "SR_SL": 20.0, "WR_RL": 25.0, "WR_SL": 50.0,
            "AG_RL": 0.8, "AG_SL": 2.5, "RL": 0.0, "SL": 0.0,
            "E": 80.0, "T": 40.0, "V": 30.0, "Gov": 20.0,
        }

    def set_government_strategy(self, strategy: str, quiet: bool = False):
        """Set government allocation strategy for E and Gov."""
        if strategy == "development":
            self.E_weights = {"r11": 0.2, "r18": 0.3, "r19": 0.2, "r21": 0.1, "r22": 0.2}
            self.Gov_weights = {"r8": 0.3, "r21": 0.2, "r23": 0.3, "r24": 0.2, "r29": 0.0, "r30": 0.0}
        elif strategy == "security":
            self.E_weights = {"r11": 0.1, "r18": 0.1, "r19": 0.1, "r21": 0.4, "r22": 0.3}
            self.Gov_weights = {"r8": 0.1, "r21": 0.2, "r23": 0.1, "r24": 0.1, "r29": 0.25, "r30": 0.25}
        elif strategy == "balanced":
            self.E_weights = {"r11": 0.2, "r18": 0.2, "r19": 0.2, "r21": 0.2, "r22": 0.2}
            self.Gov_weights = {"r8": 0.2, "r21": 0.2, "r23": 0.2, "r24": 0.2, "r29": 0.1, "r30": 0.1}
        elif strategy == "adaptive":
            # Adaptive strategy computes weights dynamically
            self.E_weights = None
            self.Gov_weights = None
        else:
            raise ValueError(f"Unknown government strategy: {strategy}")

        self.gov_strategy = strategy
        if not quiet:
            print(f"Government strategy set to: {strategy}")

    def compute_adaptive_weights(self, state: Dict[str, float]) -> Tuple[Dict[str, float], Dict[str, float]]:
        """
        Enhanced adaptive strategy with multi-objective optimization.
        
        Goals (in order of priority):
        1. Crisis Management: Reduce immediate threats (high V, high AG)
        2. Stability Building: Increase SR, reduce WR
        3. Trust Development: Increase T
        4. Economic Sustainability: Maintain/Grow E
        5. Budget Efficiency: Optimize resource allocation
        
        Returns dynamically computed E and Gov allocation weights.
        """
        # =========================================================================
        # 1. STATE ASSESSMENT & NORMALIZATION
        # =========================================================================
        # Get current state with safe defaults
        AG_total = state.get("AG_RL", 0) + state.get("AG_SL", 0)
        V = state.get("V", 0)
        SR_total = state.get("SR_RL", 0) + state.get("SR_SL", 0)
        WR_total = state.get("WR_RL", 0) + state.get("WR_SL", 0)
        T = state.get("T", 0)
        E = state.get("E", 0)
        Gov = state.get("Gov", 0)
        RL = state.get("RL", 0)
        
        # Normalized threat indicators [0-1] where 1 = severe
        ag_threat = min(1.0, AG_total / 5.0)  # Threshold at 5 AG units
        v_threat = min(1.0, V / 50.0)  # Threshold at 50 violence
        economic_pressure = max(0.0, 1.0 - (E / 100.0))  # Inverse: low E = high pressure
        trust_deficit = max(0.0, 1.0 - (T / 80.0))  # Inverse: low T = high deficit
        instability = (ag_threat + v_threat) / 2.0
        
        # Resilience metrics
        resilience_ratio = SR_total / (WR_total + 0.1)  # SR dominance
        resilience_target = min(1.0, resilience_ratio / 2.0)  # Target SR = 2x WR
        
        # =========================================================================
        # 2. STRATEGIC MODE DETERMINATION
        # =========================================================================
        # Determine which mode we're in based on threat levels
        if instability > 0.7:
            mode = "CRISIS"  # Immediate security threat
            mode_factor = {"security": 0.6, "reintegration": 0.3, "development": 0.1}
        elif instability > 0.4:
            mode = "STABILIZATION"  # Mixed threats
            mode_factor = {"security": 0.4, "reintegration": 0.3, "development": 0.3}
        elif economic_pressure > 0.6:
            mode = "ECONOMIC_RECOVERY"  # Economic crisis
            mode_factor = {"security": 0.2, "reintegration": 0.3, "development": 0.5}
        elif trust_deficit > 0.5:
            mode = "TRUST_BUILDING"  # Low trust
            mode_factor = {"security": 0.1, "reintegration": 0.2, "development": 0.7}
        else:
            mode = "SUSTAINABLE_GROWTH"  # Stable state
            mode_factor = {"security": 0.1, "reintegration": 0.3, "development": 0.6}
        
        # =========================================================================
        # 3. DYNAMIC WEIGHT CALCULATION
        # =========================================================================
        # Base allocation factors
        security_factor = mode_factor["security"]
        reintegration_factor = mode_factor["reintegration"]
        development_factor = mode_factor["development"]
        
        # Efficiency adjustment: boost when resources are high (baseline at 1.0)
        resource_efficiency = min(1.0, (E + Gov) / 150.0)  # Normalize to typical capacity
        efficiency_multiplier = 1.0 + (0.3 * resource_efficiency)  # 1.0-1.3 range
        
        # =========================================================================
        # 4. ECONOMIC (E) ALLOCATION WEIGHTS
        # =========================================================================
        # Reaction purposes:
        # r11: E → Gov conversion (institutional capacity)
        # r18: WR_RL → SR_RL (economic uplift on restored land)
        # r19: WR_SL → SR_SL (economic uplift on stressed land)
        # r21: AG_SL → WR_SL (reintegration via economic opportunity)
        # r22: SR_RL + E → T (trust from prosperity)
        
        E_weights = {
            "r11": 0.15 * efficiency_multiplier,  # Base governance support

            # Development-focused (additive: baseline + modulation bonus)
            "r18": (0.15 + 0.15 * development_factor) * efficiency_multiplier,  # WR_RL → SR_RL uplift
            "r19": (0.15 + 0.15 * development_factor) * efficiency_multiplier,  # WR_SL → SR_SL uplift
            "r22": (0.10 + 0.15 * development_factor) * efficiency_multiplier,  # SR+E → T trust building

            # Reintegration-focused (baseline + threat-modulated bonus)
            "r21": (0.15 + 0.20 * reintegration_factor * ag_threat) * efficiency_multiplier,
        }
        
        # =========================================================================
        # 5. GOVERNANCE (Gov) ALLOCATION WEIGHTS
        # =========================================================================
        # Reaction purposes:
        # r8: Gov + SR_SL → SR_RL (land restoration)
        # r21: E + Gov + AG_SL → WR_SL (reintegration - shared with E)
        # r23: Gov + SR_RL → T (trust building with strong communities)
        # r24: Gov + WR_RL → SR_RL (outreach to weak on restored land)
        # r29: Gov + AG_SL → WR_SL + V (direct conflict with AG on stressed land)
        # r30: Gov + AG_RL → AG_SL + V (displacement operations)
        
        Gov_weights = {
            # Development-focused (baseline + modulation)
            "r8": (0.10 + 0.15 * development_factor) * min(1.0, RL / 40.0) * efficiency_multiplier,  # SR_SL → SR_RL restoration
            "r23": (0.10 + 0.15 * development_factor) * efficiency_multiplier,  # Gov + SR_RL → T trust building
            "r24": (0.10 + 0.15 * development_factor) * min(1.0, WR_total / 50.0) * efficiency_multiplier,  # Gov + WR_RL → SR_RL outreach

            # Reintegration-focused (baseline + threat-modulated bonus)
            "r21": (0.15 + 0.20 * reintegration_factor * ag_threat) * efficiency_multiplier,

            # Security-focused (stronger baseline + crisis bonus)
            "r29": (0.05 + 0.25 * security_factor) * min(1.0, ag_threat + v_threat) * efficiency_multiplier,  # Gov + AG_SL → WR_SL
            "r30": (0.05 + 0.25 * security_factor) * min(1.0, ag_threat + v_threat) * efficiency_multiplier,  # Gov + AG_RL → AG_SL
        }
        
        # =========================================================================
        # 6. NORMALIZATION & VALIDATION
        # =========================================================================
        # Ensure weights sum to reasonable values (1.0-1.3 based on efficiency)
        target_sum = efficiency_multiplier
        
        # Normalize E weights
        E_sum = sum(E_weights.values())
        if E_sum > 0:
            E_weights = {k: (v / E_sum) * target_sum for k, v in E_weights.items()}
        
        # Normalize Gov weights  
        Gov_sum = sum(Gov_weights.values())
        if Gov_sum > 0:
            Gov_weights = {k: (v / Gov_sum) * target_sum for k, v in Gov_weights.items()}
        
        # Debug output (optional)
        if hasattr(self, 'debug') and self.debug:
            print(f"\n[ADAPTIVE MODE: {mode}]")
            print(f"  Threats: AG={ag_threat:.2f}, V={v_threat:.2f}, Instability={instability:.2f}")
            print(f"  Factors: Security={security_factor:.2f}, Reintegration={reintegration_factor:.2f}, Development={development_factor:.2f}")
            print(f"  Efficiency: {efficiency_multiplier:.2f}")
            print(f"  E weights: {E_weights}")
            print(f"  Gov weights: {Gov_weights}")
        
        return E_weights, Gov_weights

    def set_armed_groups_strategy(self, strategy: str, quiet: bool = False):
        """
        Set armed groups allocation strategy for V.

        Note: r25 (V+T->; trust destruction) removed - it's pure mass action, not strategic.
        Weights renormalized after r25 removal.
        """
        if strategy == "recruitment":
            # Focus on recruitment (r32), minimal displacement
            self.V_weights = {"r15": 0.125, "r16": 0.125, "r32": 0.75}
        elif strategy == "displacement":
            # Focus on displacement (r15, r16), minimal recruitment
            self.V_weights = {"r15": 0.444, "r16": 0.444, "r32": 0.111}
        elif strategy == "balanced":
            # Equal allocation across all strategic uses
            self.V_weights = {"r15": 0.333, "r16": 0.333, "r32": 0.333}
        else:
            raise ValueError(f"Unknown armed groups strategy: {strategy}")

        self.ag_strategy = strategy
        if not quiet:
            print(f"Armed groups strategy set to: {strategy}")

    def compute_allocations(self, state: Dict[str, float]) -> Tuple[Dict[str, float], Dict[str, float], Dict[str, float]]:
        """Compute allocations for E, Gov, and V based on current stocks and strategies."""
        total_E = state.get("E", 0)
        total_Gov = state.get("Gov", 0)
        total_V = state.get("V", 0)

        # For adaptive strategy, compute weights dynamically
        if self.gov_strategy == "adaptive":
            E_weights, Gov_weights = self.compute_adaptive_weights(state)
        else:
            E_weights = self.E_weights
            Gov_weights = self.Gov_weights

        # Compute allocations
        E_alloc = {rxn: weight * total_E for rxn, weight in E_weights.items()}
        Gov_alloc = {rxn: weight * total_Gov for rxn, weight in Gov_weights.items()}
        V_alloc = {rxn: weight * total_V for rxn, weight in self.V_weights.items()}

        return E_alloc, Gov_alloc, V_alloc

    def compute_reaction_rate(self, reaction: str, state: Dict[str, float],
                            E_alloc: Dict[str, float], Gov_alloc: Dict[str, float],
                            V_alloc: Dict[str, float]) -> Tuple[float, Dict[str, float]]:
        """Compute reaction rate using mass action with allocated budgets."""
        k = self.k[reaction]
        rate = k
        consumption = {"E": 0.0, "Gov": 0.0, "V": 0.0}

        if reaction == "r1":
            rate *= state.get("SL", 0)
        elif reaction == "r2":
            rate *= state.get("SR_SL", 0)
        elif reaction == "r3":
            rate *= state.get("WR_SL", 0)
        elif reaction == "r4":
            rate *= state.get("RL", 0)
        elif reaction == "r5":
            rate *= state.get("SR_RL", 0)
        elif reaction == "r6":
            rate *= state.get("WR_RL", 0)
        elif reaction == "r7":
            rate *= state.get("AG_RL", 0)
        elif reaction == "r8":
            gov_used = Gov_alloc.get("r8", 0)
            rate *= gov_used * state.get("SR_SL", 0)
            consumption["Gov"] = rate
        elif reaction == "r9":
            rate *= state.get("SR_RL", 0)
        elif reaction == "r10":
            rate *= state.get("WR_RL", 0)
        elif reaction == "r11":
            e_used = E_alloc.get("r11", 0)
            rate *= e_used
            consumption["E"] = rate
        elif reaction == "r12":
            pass
        elif reaction == "r13":
            rate *= state.get("SR_SL", 0) * state.get("RL", 0)
        elif reaction == "r14":
            rate *= state.get("WR_SL", 0) * state.get("RL", 0)
        elif reaction == "r15":
            v_used = V_alloc.get("r15", 0)
            rate *= state.get("SR_RL", 0) * state.get("SL", 0) * v_used
            consumption["V"] = rate
        elif reaction == "r16":
            v_used = V_alloc.get("r16", 0)
            rate *= state.get("WR_RL", 0) * state.get("SL", 0) * v_used
            consumption["V"] = rate
        elif reaction == "r17":
            rate *= state.get("AG_SL", 0) * state.get("RL", 0)
        elif reaction == "r18":
            e_used = E_alloc.get("r18", 0)
            rate *= e_used * state.get("WR_RL", 0)
            consumption["E"] = rate
        elif reaction == "r19":
            e_used = E_alloc.get("r19", 0)
            rate *= e_used * state.get("WR_SL", 0)
            consumption["E"] = rate
        elif reaction == "r20":
            rate *= state.get("SR_SL", 0)
        elif reaction == "r21":
            e_alloc = E_alloc.get("r21", 0)
            gov_alloc = Gov_alloc.get("r21", 0)
            if e_alloc > 0 and gov_alloc > 0:
                combined_budget = sqrt(e_alloc * gov_alloc)
            else:
                combined_budget = 0
            rate *= combined_budget * state.get("AG_SL", 0)
            if combined_budget > 0:
                consumption["E"] = rate * sqrt(e_alloc / gov_alloc) if gov_alloc > 0 else rate
                consumption["Gov"] = rate * sqrt(gov_alloc / e_alloc) if e_alloc > 0 else rate
        elif reaction == "r22":
            rate *= state.get("SR_RL", 0) * E_alloc.get("r22", 0)
        elif reaction == "r23":
            rate *= Gov_alloc.get("r23", 0) * state.get("SR_RL", 0)
        elif reaction == "r24":
            gov_used = Gov_alloc.get("r24", 0)
            rate *= gov_used * state.get("WR_RL", 0)
            consumption["Gov"] = rate
        elif reaction == "r25":
            # r25: V + T => ; (Trust destruction by violence - pure mass action)
            # This is NOT budget-controlled - violence automatically destroys trust
            rate *= state.get("V", 0) * state.get("T", 0)
            # V consumption handled by stoichiometry, not manual tracking
        elif reaction == "r26":
            rate *= state.get("Gov", 0)
        elif reaction == "r27":
            rate *= state.get("WR_SL", 0)
        elif reaction == "r28":
            rate *= state.get("AG_SL", 0) * state.get("SR_RL", 0)
        elif reaction == "r29":
            gov_used = Gov_alloc.get("r29", 0)
            rate *= state.get("AG_SL", 0) * gov_used
            consumption["Gov"] = rate
        elif reaction == "r30":
            gov_used = Gov_alloc.get("r30", 0)
            rate *= state.get("AG_RL", 0) * gov_used
            consumption["Gov"] = rate
        elif reaction == "r31":
            rate *= state.get("V", 0)
        elif reaction == "r32":
            v_used = V_alloc.get("r32", 0)
            rate *= state.get("WR_SL", 0) * state.get("AG_SL", 0) * v_used
            consumption["V"] = rate
        elif reaction == "r33":
            rate *= state.get("AG_SL", 0) * state.get("WR_RL", 0)
        elif reaction == "r34":
            rate *= state.get("AG_SL", 0) * state.get("SR_RL", 0)
        elif reaction == "r35":
            rate *= state.get("AG_SL", 0) * state.get("E", 0)
            consumption["E"] = rate
        elif reaction == "r36":
            rate *= state.get("AG_RL", 0) * state.get("Gov", 0)
            consumption["Gov"] = rate

        # Apply seasonal variation to resource recovery reactions
        seasonal_modifier = self.get_seasonal_modifier(reaction)
        rate *= seasonal_modifier

        # Update consumption proportionally with seasonal modifier
        for resource in consumption:
            consumption[resource] *= seasonal_modifier

        return max(0.0, rate), consumption

    def diagnose_decay_reactions(self, state: Dict[str, float], process_vector: np.ndarray,
                                 scale_factor: float = 1.0) -> Dict[str, Any]:
        """
        Diagnostic function to investigate why r25, r26, r31 might become inactive.

        Returns detailed information about decay reaction rates and their dependencies.
        """
        diagnostics = {}

        # r25: V + T => ; (trust destruction - FIXED to pure mass action)
        r25_idx = self.REACTIONS.index("r25")
        V = state.get("V", 0)
        T = state.get("T", 0)
        k25 = self.k["r25"]
        expected_r25_rate = k25 * V * T
        actual_r25_rate = process_vector[r25_idx]
        diagnostics["r25"] = {
            "V": V,
            "T": T,
            "k": k25,
            "expected_rate": expected_r25_rate,
            "actual_rate": actual_r25_rate,
            "scaled_rate": actual_r25_rate * scale_factor,
            "match": abs(expected_r25_rate - actual_r25_rate) < 1e-10,
            "description": "Trust destruction by violence (pure mass action)"
        }

        # r26: Gov => ; (governance decay - should be pure mass action)
        r26_idx = self.REACTIONS.index("r26")
        Gov = state.get("Gov", 0)
        k26 = self.k["r26"]
        expected_r26_rate = k26 * Gov
        actual_r26_rate = process_vector[r26_idx]
        diagnostics["r26"] = {
            "Gov": Gov,
            "k": k26,
            "expected_rate": expected_r26_rate,
            "actual_rate": actual_r26_rate,
            "scaled_rate": actual_r26_rate * scale_factor,
            "match": abs(expected_r26_rate - actual_r26_rate) < 1e-10,
            "description": "Governance decay (pure mass action)"
        }

        # r31: V => ; (violence decay - should be pure mass action)
        r31_idx = self.REACTIONS.index("r31")
        k31 = self.k["r31"]
        expected_r31_rate = k31 * V
        actual_r31_rate = process_vector[r31_idx]
        diagnostics["r31"] = {
            "V": V,
            "k": k31,
            "expected_rate": expected_r31_rate,
            "actual_rate": actual_r31_rate,
            "scaled_rate": actual_r31_rate * scale_factor,
            "match": abs(expected_r31_rate - actual_r31_rate) < 1e-10,
            "description": "Violence decay (pure mass action)"
        }

        return diagnostics

    def compute_process_vector(self, state: Dict[str, float], verbose: bool = False) -> Tuple[np.ndarray, Dict[str, Any], Dict[str, float], Dict[str, float]]:
        """Compute process vector with current allocations."""
        E_alloc, Gov_alloc, V_alloc = self.compute_allocations(state)

        process_vector = np.zeros(len(self.REACTIONS))
        rates = {}
        total_consumption = {"E": 0.0, "Gov": 0.0, "V": 0.0}
        total_generation = {"E": 0.0, "Gov": 0.0}
        consumption_by_rxn = {}

        for i, reaction in enumerate(self.REACTIONS):
            rate, consumption = self.compute_reaction_rate(reaction, state, E_alloc, Gov_alloc, V_alloc)
            process_vector[i] = rate
            rates[reaction] = rate
            consumption_by_rxn[reaction] = consumption

            # Accumulate consumption
            for resource in ["E", "Gov", "V"]:
                total_consumption[resource] += consumption[resource]

            # Track generation
            if reaction in ["r9", "r10"]:  # E production
                total_generation["E"] += rate
            if reaction in ["r11", "r12"]:  # Gov production
                total_generation["Gov"] += rate

        allocation_info = {
            "E": {"total": state.get("E", 0), "allocations": E_alloc},
            "Gov": {"total": state.get("Gov", 0), "allocations": Gov_alloc},
            "V": {"total": state.get("V", 0), "allocations": V_alloc},
            "consumption_by_rxn": consumption_by_rxn,
        }

        return process_vector, allocation_info, total_consumption, total_generation

    def apply_timestep(self, dt: float = 0.1, verbose: bool = False):
        """Apply one timestep of the simulation."""
        old_state = self.state.copy()

        # Compute process vector with allocations, consumption, and generation
        process_vector, allocation_info, total_consumption, total_generation = self.compute_process_vector(self.state, verbose)

        # Rate limiter
        delta_raw = self.S @ process_vector * dt
        scale_factor = 1.0
        for i, species in enumerate(self.SPECIES):
            current = self.state.get(species, 0)
            if delta_raw[i] < 0 and current > 0:
                max_consumption = current * 0.5
                if abs(delta_raw[i]) > max_consumption:
                    species_scale = max_consumption / abs(delta_raw[i])
                    scale_factor = min(scale_factor, species_scale)

        if scale_factor < 1.0:
            process_vector = process_vector * scale_factor
            total_consumption = {k: v * scale_factor for k, v in total_consumption.items()}
            total_generation = {k: v * scale_factor for k, v in total_generation.items()}

        # Compute state changes via stoichiometry
        delta = self.S @ process_vector * dt

        # Update state
        new_state = {}
        for i, species in enumerate(self.SPECIES):
            new_value = self.state.get(species, 0) + delta[i]
            new_state[species] = max(0.0, new_value)

        # Manual tracking for E, Gov, V
        e_production = sum(process_vector[i] * dt for i, rxn in enumerate(self.REACTIONS) if rxn in ["r9", "r10"])
        new_state["E"] = max(0.0, old_state["E"] + e_production - total_consumption["E"] * dt)

        gov_production = sum(process_vector[i] * dt for i, rxn in enumerate(self.REACTIONS) if rxn in ["r11", "r12"])
        r26_idx = self.REACTIONS.index("r26")
        gov_decay = process_vector[r26_idx] * dt
        new_state["Gov"] = max(0.0, old_state["Gov"] + gov_production - gov_decay - total_consumption["Gov"] * dt)

        v_production = sum(process_vector[i] * dt for i, rxn in enumerate(self.REACTIONS) if rxn in ["r27", "r28", "r29", "r30"])
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
        self.generation_history.append(total_generation)

    def run_simulation(self, total_time: float = 20.0, dt: float = 0.1, verbose: bool = False):
        """Run simulation for specified time."""
        steps = int(total_time / dt)
        for step in range(steps):
            self.apply_timestep(dt, verbose=False)

# ============================================================================
# BATCH SIMULATION FUNCTIONS
# ============================================================================

def run_batch_simulation(model_file: str, climate_params: List[float],
                        budget_renewal_rates: List[float],
                        total_time: float = 20.0, dt: float = 0.1):
    """Run comprehensive batch simulation."""
    print("="*70)
    print("BATCH CONFLICT DYNAMICS SIMULATION - ADAPTIVE")
    print("="*70)
    print(f"\nSimulation parameters:")
    print(f"  Total time: {total_time} years")
    print(f"  Timestep: {dt} years")
    print(f"  Conflict scenarios: {list(CONFLICT_SCENARIOS.keys())}")
    print(f"  Climate parameters: {climate_params}")
    print(f"  Budget renewal rates: {budget_renewal_rates}")
    print(f"  Government strategies: {GOVERNMENT_STRATEGIES}")
    print(f"  Armed group strategies: {ARMED_GROUP_STRATEGIES}")

    total_runs = (len(CONFLICT_SCENARIOS) * len(climate_params) * len(budget_renewal_rates) *
                  len(GOVERNMENT_STRATEGIES) * len(ARMED_GROUP_STRATEGIES))
    print(f"\nTotal simulation runs: {total_runs}")
    print("="*70)

    simulator = ConflictSimulator(model_file)
    results = {}
    run_count = 0

    for scenario in CONFLICT_SCENARIOS:
        results[scenario] = {}
        print(f"\n>>> Scenario: {scenario.upper()} - {CONFLICT_SCENARIOS[scenario]['description']}")

        for climate_k in climate_params:
            results[scenario][climate_k] = {}

            for renewal_rate in budget_renewal_rates:
                results[scenario][climate_k][renewal_rate] = {}

                for gov_strategy in GOVERNMENT_STRATEGIES:
                    results[scenario][climate_k][renewal_rate][gov_strategy] = {}

                    for ag_strategy in ARMED_GROUP_STRATEGIES:
                        run_count += 1
                        print(f"  Run {run_count}/{total_runs}: climate={climate_k}, "
                              f"renewal={renewal_rate}, gov={gov_strategy}, ag={ag_strategy}", end="")

                        # Configure simulation
                        simulator.set_conflict_scenario(scenario)
                        simulator.set_climate_parameter(climate_k)
                        simulator.set_budget_renewal_rate(renewal_rate)
                        simulator.set_government_strategy(gov_strategy, quiet=True)
                        simulator.set_armed_groups_strategy(ag_strategy, quiet=True)

                        # Run simulation
                        simulator.run_simulation(total_time=total_time, dt=dt, verbose=False)

                        # Extract metrics
                        metrics = simulator.get_metrics()
                        metrics["scenario"] = scenario
                        metrics["climate_k"] = climate_k
                        metrics["renewal_rate"] = renewal_rate
                        metrics["gov_strategy"] = gov_strategy
                        metrics["ag_strategy"] = ag_strategy

                        # Store simulator history data for process analysis
                        # Create a lightweight copy with only essential data
                        import copy
                        sim_copy = type('SimulatorData', (), {})()
                        sim_copy.rate_history = copy.deepcopy(simulator.rate_history)
                        sim_copy.time_history = copy.deepcopy(simulator.time_history)
                        sim_copy.state_history = copy.deepcopy(simulator.state_history)
                        sim_copy.S = simulator.S
                        sim_copy.REACTIONS = simulator.REACTIONS

                        metrics["simulator"] = sim_copy
                        results[scenario][climate_k][renewal_rate][gov_strategy][ag_strategy] = metrics
                        print(f" -> AG={metrics['final_AG']:.1f}, Budget_net={metrics['net_budget']:.1f}")

    print("\n" + "="*70)
    print("BATCH SIMULATION COMPLETE")
    print("="*70)

    return results


# ============================================================================
# PROCESS-BASED TEMPORAL ANALYSIS FUNCTIONS
# ============================================================================

def setup_semantic_partition(species_list: List[str]) -> SemanticPartition:
    """
    Create semantic partition for conflict dynamics.

    Categories:
    - peace: Species associated with stability and prosperity (SR, E, T, WR_RL, RL)
    - conflict: Species associated with instability and violence (AG, V, WR_SL, SL)

    Args:
        species_list: List of species names from the reaction network

    Returns:
        SemanticPartition object
    """
    # Define peace and conflict categories based on species semantics
    category_definitions = {
        'peace': ['SR_RL', 'SR_SL', 'E', 'T', 'WR_RL', 'RL'],  # Prosperity and stability
        'conflict': ['AG_RL', 'AG_SL', 'V', 'WR_SL', 'SL']     # Violence and instability
    }

    # Filter to only include species that exist in the model
    filtered_definitions = {}
    for category, species in category_definitions.items():
        filtered_species = [sp for sp in species if sp in species_list]
        if filtered_species:  # Only add category if it has at least one species
            filtered_definitions[category] = filtered_species

    # Ensure both categories exist, even if empty
    if 'peace' not in filtered_definitions:
        filtered_definitions['peace'] = []
    if 'conflict' not in filtered_definitions:
        filtered_definitions['conflict'] = []

    # Prevent empty categories from causing errors
    if not filtered_definitions['peace'] or not filtered_definitions['conflict']:
        raise ValueError(f"Empty semantic categories! Peace: {filtered_definitions['peace']}, Conflict: {filtered_definitions['conflict']}, Available species: {species_list}")

    return define_semantic_categories(species_list, filtered_definitions)


def analyze_process_evolution(sim: 'ConflictSimulator',
                              semantic_partition: SemanticPartition,
                              time_scales: Optional[List[int]] = None) -> Dict:
    """
    Analyze how processes evolve over time across semantic categories and time scales.

    Args:
        sim: ConflictSimulator object after simulation
        semantic_partition: Semantic partition defining peace/conflict categories
        time_scales: List of window sizes for multi-scale analysis (default: [1, 5, 10, 20])

    Returns:
        Dictionary with analysis results
    """
    if time_scales is None:
        # Expanded temporal scales: from single timestep to half the simulation
        n_steps = len(sim.rate_history)
        time_scales = [1, 3, 5, 10, 15, 20, 30, 40, 50]  # More granular temporal resolutions
        # Filter to only include scales smaller than available data
        time_scales = [s for s in time_scales if s < n_steps]

    # Convert rate_history to process vectors (rates are already process vectors)
    process_sequence = []
    for rates in sim.rate_history:
        # Extract process vector from reaction rates
        process_vector = np.array([rates[rxn] for rxn in sim.REACTIONS])
        process_sequence.append(process_vector)

    if not process_sequence:
        return {"error": "No process history available"}

    # Get stoichiometric matrix
    S = sim.S

    results = {
        'time_points': sim.time_history,
        'n_timesteps': len(process_sequence),
        'semantic_categories': semantic_partition.categories,
        'time_scales': time_scales,
        'single_step_analysis': {},
        'multiscale_analysis': {}
    }

    # ========================================================================
    # 1. SINGLE-STEP PROCESS MODE CLASSIFICATION
    # ========================================================================
    # Track how each timestep's process is classified for each category
    peace_modes = []
    conflict_modes = []
    overall_modes = []
    peace_cognitive = []
    conflict_cognitive = []

    for process in process_sequence:
        # Overall classification
        mode, completeness = classify_process_mode(process, S)
        overall_modes.append(mode)

        # Peace category classification
        peace_indices = semantic_partition.category_indices['peace']
        peace_mode, _ = classify_process_mode(process, S, species_subset=peace_indices)
        peace_modes.append(peace_mode)
        peace_cognitive.append(is_cognitive_domain(process, S, species_subset=peace_indices))

        # Conflict category classification
        conflict_indices = semantic_partition.category_indices['conflict']
        conflict_mode, _ = classify_process_mode(process, S, species_subset=conflict_indices)
        conflict_modes.append(conflict_mode)
        conflict_cognitive.append(is_cognitive_domain(process, S, species_subset=conflict_indices))

    results['single_step_analysis'] = {
        'overall_modes': overall_modes,
        'peace_modes': peace_modes,
        'conflict_modes': conflict_modes,
        'peace_in_cognitive_domain': peace_cognitive,
        'conflict_in_cognitive_domain': conflict_cognitive,
        'mode_distribution': {
            'overall': dict(zip(*np.unique(overall_modes, return_counts=True))),
            'peace': dict(zip(*np.unique(peace_modes, return_counts=True))),
            'conflict': dict(zip(*np.unique(conflict_modes, return_counts=True)))
        },
        'cognitive_domain_fraction': {
            'peace': np.mean(peace_cognitive),
            'conflict': np.mean(conflict_cognitive)
        }
    }

    # ========================================================================
    # 2. CATEGORY BEHAVIOR TRACKING
    # ========================================================================
    # Track production/consumption for each category over time
    category_sequences = track_category_sequence(process_sequence, S, semantic_partition)

    category_stats = {}
    for category, behaviors in category_sequences.items():
        category_stats[category] = compute_category_statistics(behaviors)

    results['category_behavior'] = {
        'sequences': category_sequences,
        'statistics': category_stats
    }

    # ========================================================================
    # 3. MULTI-SCALE TEMPORAL ANALYSIS
    # ========================================================================
    # Aggregate processes at different time scales and analyze resulting patterns
    multiscale_results = {}

    for window_size in time_scales:
        if window_size >= len(process_sequence):
            continue  # Skip windows larger than available data

        # Apply rolling window aggregation
        aggregations = rolling_window_aggregation(
            process_sequence,
            window_size=window_size,
            stride=max(1, window_size // 2)  # 50% overlap
        )

        # Analyze aggregated processes
        scale_peace_modes = []
        scale_conflict_modes = []
        scale_overall_modes = []

        for agg in aggregations:
            aggregated_v = agg.aggregated_process

            # Classify aggregated process
            overall_mode, _ = classify_process_mode(aggregated_v, S)
            scale_overall_modes.append(overall_mode)

            peace_mode, _ = classify_process_mode(aggregated_v, S,
                                                  species_subset=semantic_partition.category_indices['peace'])
            scale_peace_modes.append(peace_mode)

            conflict_mode, _ = classify_process_mode(aggregated_v, S,
                                                     species_subset=semantic_partition.category_indices['conflict'])
            scale_conflict_modes.append(conflict_mode)

        multiscale_results[window_size] = {
            'n_windows': len(aggregations),
            'overall_modes': scale_overall_modes,
            'peace_modes': scale_peace_modes,
            'conflict_modes': scale_conflict_modes,
            'mode_distribution': {
                'overall': dict(zip(*np.unique(scale_overall_modes, return_counts=True))) if scale_overall_modes else {},
                'peace': dict(zip(*np.unique(scale_peace_modes, return_counts=True))) if scale_peace_modes else {},
                'conflict': dict(zip(*np.unique(scale_conflict_modes, return_counts=True))) if scale_conflict_modes else {}
            }
        }

    results['multiscale_analysis'] = multiscale_results

    return results


def analyze_steady_state_mechanism(sim: 'ConflictSimulator', variable: str,
                                   time_window: int = 20) -> Dict:
    """
    Analyze why a variable reaches steady state by examining reaction balance.

    Args:
        sim: ConflictSimulator with completed simulation
        variable: Variable to analyze ('V', 'T', 'Gov', 'E')
        time_window: Number of final timesteps to analyze

    Returns:
        Dictionary with production/consumption analysis
    """
    # Map variables to their producing and consuming reactions
    reaction_effects = {
        'V': {
            'producers': ['r27', 'r28', 'r29', 'r30'],
            'consumers': ['r15', 'r16', 'r25', 'r31'],
            'descriptions': {
                'r27': 'WR_SL desperation violence',
                'r28': 'AG attacks on SR',
                'r29': 'Gov-AG conflict (stressed)',
                'r30': 'Gov displacement ops',
                'r15': 'Violence in displacement',
                'r16': 'Violence in displacement',
                'r25': 'Trust destruction',
                'r31': 'Violence decay'
            }
        },
        'T': {
            'producers': ['r22', 'r23'],
            'consumers': ['r25'],
            'descriptions': {
                'r22': 'Trust from prosperity',
                'r23': 'Gov trust building',
                'r25': 'Violence destroys trust'
            }
        },
        'Gov': {
            'producers': ['r11', 'r12'],
            'consumers': ['r8', 'r21', 'r24', 'r29', 'r30', 'r26', 'r36'],
            'descriptions': {
                'r11': 'E → Gov conversion',
                'r12': 'Baseline Gov renewal',
                'r8': 'Land restoration',
                'r21': 'Reintegration',
                'r24': 'Outreach to weak',
                'r29': 'Conflict with AG_SL',
                'r30': 'Displacement of AG_RL',
                'r26': 'Gov corruption/decay',
                'r36': 'AG institutional capture'
            }
        },
        'E': {
            'producers': ['r9', 'r10'],
            'consumers': ['r11', 'r18', 'r19', 'r21', 'r22', 'r35'],
            'descriptions': {
                'r9': 'SR_RL production',
                'r10': 'WR_RL production',
                'r11': 'E → Gov conversion',
                'r18': 'Economic uplift (RL)',
                'r19': 'Economic uplift (SL)',
                'r21': 'Reintegration',
                'r22': 'Trust building',
                'r35': 'AG corruption'
            }
        }
    }

    if variable not in reaction_effects:
        return {"error": f"Unknown variable: {variable}"}

    effects = reaction_effects[variable]

    # Analyze final timesteps
    n_steps = len(sim.rate_history)
    start_idx = max(0, n_steps - time_window)

    # Compute average rates for producers and consumers
    producer_rates = {rxn: [] for rxn in effects['producers']}
    consumer_rates = {rxn: [] for rxn in effects['consumers']}

    for i in range(start_idx, n_steps):
        rates = sim.rate_history[i]
        for rxn in effects['producers']:
            producer_rates[rxn].append(rates.get(rxn, 0))
        for rxn in effects['consumers']:
            consumer_rates[rxn].append(rates.get(rxn, 0))

    # Average over time window
    avg_producers = {rxn: np.mean(rates_list) for rxn, rates_list in producer_rates.items()}
    avg_consumers = {rxn: np.mean(rates_list) for rxn, rates_list in consumer_rates.items()}

    total_production = sum(avg_producers.values())
    total_consumption = sum(avg_consumers.values())
    net_rate = total_production - total_consumption

    # Check if variable is at steady state
    state_history = sim.state_history
    var_values = [state[variable] for state in state_history[start_idx:]]
    var_std = np.std(var_values)
    var_mean = np.mean(var_values)
    is_steady_state = var_std < 0.01 * abs(var_mean) if var_mean != 0 else var_std < 0.01

    # Identify inactive reactions (rate < threshold)
    inactive_threshold = 1e-6
    inactive_producers = [rxn for rxn, rate in avg_producers.items() if rate < inactive_threshold]
    inactive_consumers = [rxn for rxn, rate in avg_consumers.items() if rate < inactive_threshold]

    return {
        'variable': variable,
        'is_steady_state': is_steady_state,
        'mean_value': var_mean,
        'std_value': var_std,
        'total_production': total_production,
        'total_consumption': total_consumption,
        'net_rate': net_rate,
        'balance_quality': abs(net_rate) / max(total_production, total_consumption, 1e-10),
        'producer_rates': avg_producers,
        'consumer_rates': avg_consumers,
        'inactive_producers': inactive_producers,
        'inactive_consumers': inactive_consumers,
        'active_producers': [rxn for rxn in effects['producers'] if rxn not in inactive_producers],
        'active_consumers': [rxn for rxn in effects['consumers'] if rxn not in inactive_consumers],
        'descriptions': effects['descriptions']
    }


def batch_process_analysis(results: Dict, model_file: str,
                           climate_params: List[float],
                           budget_renewal_rates: List[float]) -> Dict:
    """
    Perform process-based analysis on all batch simulation results.

    Args:
        results: Results dictionary from run_batch_simulation
        model_file: Path to model file
        climate_params: List of climate parameters
        budget_renewal_rates: List of budget renewal rates

    Returns:
        Dictionary with aggregated process analysis results
    """
    print("\n" + "="*70)
    print("PROCESS-BASED TEMPORAL ANALYSIS")
    print("="*70)

    # Load model to get species list
    rn = read_txt(model_file)
    species_objects = rn.species()
    species_list = [sp.name for sp in species_objects]  # Extract names from Species objects

    # Setup semantic partition
    semantic_partition = setup_semantic_partition(species_list)
    print(f"\nSemantic Categories:")
    print(f"  Peace: {semantic_partition.category_to_species['peace']}")
    print(f"  Conflict: {semantic_partition.category_to_species['conflict']}")

    # Analyze a representative simulation (medium conflict, middle climate, middle renewal, adaptive strategy)
    scenario = "medium"
    climate_k = climate_params[len(climate_params)//2]
    renewal_rate = budget_renewal_rates[len(budget_renewal_rates)//2]

    print(f"\nAnalyzing representative simulation:")
    print(f"  Scenario: {scenario}, Climate: {climate_k}, Renewal: {renewal_rate}")

    process_analyses = {}
    for gov_strat in GOVERNMENT_STRATEGIES:
        ag_strat = "balanced"  # Use balanced AG strategy for consistency
        sim_data = results[scenario][climate_k][renewal_rate][gov_strat][ag_strat]

        if 'simulator' in sim_data:
            sim = sim_data['simulator']
            analysis = analyze_process_evolution(sim, semantic_partition)
            process_analyses[gov_strat] = analysis

            print(f"\n  {gov_strat.upper()} Strategy:")
            print(f"    Peace cognitive domain: {analysis['single_step_analysis']['cognitive_domain_fraction']['peace']:.1%}")
            print(f"    Conflict cognitive domain: {analysis['single_step_analysis']['cognitive_domain_fraction']['conflict']:.1%}")
            print(f"    Peace mode distribution: {analysis['single_step_analysis']['mode_distribution']['peace']}")
            print(f"    Conflict mode distribution: {analysis['single_step_analysis']['mode_distribution']['conflict']}")

    print("\nProcess analysis complete!")

    # ========================================================================
    # STEADY STATE MECHANISM ANALYSIS
    # ========================================================================
    print("\n" + "="*70)
    print("STEADY STATE MECHANISM ANALYSIS")
    print("="*70)

    steady_state_analyses = {}
    variables_to_check = ['V', 'T', 'Gov']

    for gov_strat in GOVERNMENT_STRATEGIES:
        ag_strat = "balanced"
        sim_data = results[scenario][climate_k][renewal_rate][gov_strat][ag_strat]

        if 'simulator' not in sim_data:
            continue

        sim = sim_data['simulator']
        steady_state_analyses[gov_strat] = {}

        print(f"\n{gov_strat.upper()} Strategy:")
        for var in variables_to_check:
            analysis = analyze_steady_state_mechanism(sim, var)
            steady_state_analyses[gov_strat][var] = analysis

            status = "STEADY STATE" if analysis['is_steady_state'] else "DYNAMIC"
            print(f"\n  {var} ({status}):")
            print(f"    Mean: {analysis['mean_value']:.3f}, Std: {analysis['std_value']:.3f}")
            print(f"    Production: {analysis['total_production']:.4f}")
            print(f"    Consumption: {analysis['total_consumption']:.4f}")
            print(f"    Net Rate: {analysis['net_rate']:.4f} (balance: {analysis['balance_quality']:.1%})")

            if analysis['is_steady_state']:
                print(f"    ✓ BALANCED: Production ≈ Consumption")
                if analysis['active_producers']:
                    print(f"    Active producers: {', '.join(analysis['active_producers'])}")
                if analysis['active_consumers']:
                    print(f"    Active consumers: {', '.join(analysis['active_consumers'])}")
                if analysis['inactive_producers']:
                    print(f"    Inactive producers: {', '.join(analysis['inactive_producers'])}")
                if analysis['inactive_consumers']:
                    print(f"    Inactive consumers: {', '.join(analysis['inactive_consumers'])}")
            else:
                print(f"    → Net {'production' if analysis['net_rate'] > 0 else 'consumption'}")

    # Add to process_analyses for visualization
    for gov_strat in steady_state_analyses:
        if gov_strat in process_analyses:
            process_analyses[gov_strat]['steady_state_analysis'] = steady_state_analyses[gov_strat]

    return process_analyses


def plot_comparison_results(results: Dict, climate_params: List[float],
                           budget_renewal_rates: List[float], save_dir: str = None):
    """Generate comprehensive comparison plots with flexible parameters."""

    if save_dir:
        os.makedirs(save_dir, exist_ok=True)

    n_climates = len(climate_params)
    n_renewals = len(budget_renewal_rates)

    # ========================================================================
    # PLOT 1: Budget Spent vs Generated (for each strategy)
    # ========================================================================
    print("\nGenerating Plot 1: Budget Spent vs Generated Comparison...")

    # Use middle scenario and middle renewal rate for comparison
    scenario = "medium"
    renewal_rate = budget_renewal_rates[len(budget_renewal_rates)//2]

    fig, axes = plt.subplots(1, n_climates, figsize=(5*n_climates, 5))
    if n_climates == 1:
        axes = [axes]
    fig.suptitle(f"Total Budget: Generated vs Consumed ({scenario.title()}, renewal={renewal_rate})", fontsize=14)

    for col, climate_k in enumerate(climate_params):
        ax = axes[col]

        x_labels = []
        budget_generated = []
        budget_consumed = []

        for gov_strat in GOVERNMENT_STRATEGIES:
            for ag_strat in ARMED_GROUP_STRATEGIES:
                metrics = results[scenario][climate_k][renewal_rate][gov_strat][ag_strat]
                x_labels.append(f"{gov_strat[:4]}-{ag_strat[:4]}")
                budget_generated.append(metrics["total_budget_generated"])
                budget_consumed.append(metrics["total_budget_consumed"])

        x = np.arange(len(x_labels))
        width = 0.35

        ax.bar(x - width/2, budget_generated, width, label='Generated (E+Gov)', color='lightgreen', alpha=0.8)
        ax.bar(x + width/2, budget_consumed, width, label='Consumed (E+Gov)', color='coral', alpha=0.8)

        # Show baseline (initial E + Gov)
        initial_budget = (CONFLICT_SCENARIOS[scenario]["initial_state"]["E"] +
                         CONFLICT_SCENARIOS[scenario]["initial_state"]["Gov"])
        ax.axhline(y=initial_budget, color='gray', linestyle='--', linewidth=1.5,
                   label=f'Initial: {initial_budget:.0f}')

        ax.set_title(f"Climate={climate_k}")
        ax.set_xticks(x)
        ax.set_xticklabels(x_labels, rotation=45, ha='right', fontsize=8)
        ax.set_ylabel("Total Budget (E+Gov)")
        if col == 0:
            ax.legend(loc='upper left', fontsize=8)
        ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    if save_dir:
        plt.savefig(os.path.join(save_dir, "01_budget_comparison.png"), dpi=150)
    plt.close()

    # ========================================================================
    # PLOT 2: Net Budget Efficiency by Renewal Rate
    # ========================================================================
    print("Generating Plot 2: Net Budget Efficiency...")

    fig, axes = plt.subplots(len(CONFLICT_SCENARIOS), 1, figsize=(14, 4*len(CONFLICT_SCENARIOS)))
    if len(CONFLICT_SCENARIOS) == 1:
        axes = [axes]
    fig.suptitle("Net Budget (Generated - Consumed) by Renewal Rate", fontsize=14)

    for row, scenario in enumerate(CONFLICT_SCENARIOS.keys()):
        ax = axes[row]

        # Use middle climate
        climate_k = climate_params[len(climate_params)//2]

        for gov_strat in GOVERNMENT_STRATEGIES:
            net_budgets = []
            for renewal_rate in budget_renewal_rates:
                # Average across AG strategies
                values = [results[scenario][climate_k][renewal_rate][gov_strat][ag_strat]["net_budget"]
                         for ag_strat in ARMED_GROUP_STRATEGIES]
                net_budgets.append(np.mean(values))

            marker = 'o' if gov_strat == "development" else ('s' if gov_strat == "security" else ('^' if gov_strat == "balanced" else '*'))
            linewidth = 3 if gov_strat == "adaptive" else 2
            ax.plot(budget_renewal_rates, net_budgets, marker=marker, markersize=10 if gov_strat == "adaptive" else 8,
                   label=gov_strat, linewidth=linewidth)

        ax.axhline(y=0, color='black', linestyle='--', linewidth=1, alpha=0.5)
        ax.set_title(f"{scenario.title()} Conflict (climate={climate_k})")
        ax.set_xlabel("Budget Renewal Rate (r12)")
        ax.set_ylabel("Net Budget (Generated - Consumed)")
        ax.legend(loc='best')
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    if save_dir:
        plt.savefig(os.path.join(save_dir, "02_net_budget_efficiency.png"), dpi=150)
    plt.close()

    # ========================================================================
    # PLOT 3: Prosperity (E and T) with Baselines
    # ========================================================================
    print("Generating Plot 3: Prosperity (E and T) with Baselines...")

    scenario = "medium"
    renewal_rate = budget_renewal_rates[len(budget_renewal_rates)//2]

    fig, axes = plt.subplots(2, n_climates, figsize=(5*n_climates, 10))
    if n_climates == 1:
        axes = axes.reshape(-1, 1)
    fig.suptitle(f"Final Prosperity: Economy & Trust ({scenario.title()}, renewal={renewal_rate})", fontsize=14)

    # Get baselines
    baseline_E = CONFLICT_SCENARIOS[scenario]["initial_state"]["E"]
    baseline_T = CONFLICT_SCENARIOS[scenario]["initial_state"]["T"]

    for col, climate_k in enumerate(climate_params):
        # Economy subplot
        ax_e = axes[0, col]
        x_labels = []
        final_e = []

        for gov_strat in GOVERNMENT_STRATEGIES:
            for ag_strat in ARMED_GROUP_STRATEGIES:
                metrics = results[scenario][climate_k][renewal_rate][gov_strat][ag_strat]
                x_labels.append(f"{gov_strat[:4]}-{ag_strat[:4]}")
                final_e.append(metrics["final_E"])

        x = np.arange(len(x_labels))
        colors = ['darkgreen' if 'adap' in label else 'green' for label in x_labels]
        ax_e.bar(x, final_e, color=colors, alpha=0.7)
        ax_e.axhline(y=baseline_E, color='red', linestyle='--', linewidth=2,
                    label=f'Baseline: {baseline_E:.0f}')
        ax_e.set_title(f"Economy (E) - Climate={climate_k}")
        ax_e.set_xticks(x)
        ax_e.set_xticklabels(x_labels, rotation=45, ha='right', fontsize=7)
        ax_e.set_ylabel("Final E")
        ax_e.legend(loc='upper right', fontsize=8)
        ax_e.grid(True, alpha=0.3, axis='y')

        # Trust subplot
        ax_t = axes[1, col]
        final_t = []

        for gov_strat in GOVERNMENT_STRATEGIES:
            for ag_strat in ARMED_GROUP_STRATEGIES:
                metrics = results[scenario][climate_k][renewal_rate][gov_strat][ag_strat]
                final_t.append(metrics["final_T"])

        colors = ['darkcyan' if 'adap' in label else 'teal' for label in x_labels]
        ax_t.bar(x, final_t, color=colors, alpha=0.7)
        ax_t.axhline(y=baseline_T, color='red', linestyle='--', linewidth=2,
                    label=f'Baseline: {baseline_T:.0f}')
        ax_t.set_title(f"Trust (T) - Climate={climate_k}")
        ax_t.set_xticks(x)
        ax_t.set_xticklabels(x_labels, rotation=45, ha='right', fontsize=7)
        ax_t.set_ylabel("Final T")
        ax_t.legend(loc='upper right', fontsize=8)
        ax_t.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    if save_dir:
        plt.savefig(os.path.join(save_dir, "03_prosperity_baselines.png"), dpi=150)
    plt.close()

    # ========================================================================
    # PLOT 4: Strong Resilience with Baseline
    # ========================================================================
    print("Generating Plot 4: Strong Resilience with Baseline...")

    fig, axes = plt.subplots(1, len(CONFLICT_SCENARIOS), figsize=(6*len(CONFLICT_SCENARIOS), 5))
    if len(CONFLICT_SCENARIOS) == 1:
        axes = [axes]
    fig.suptitle(f"Final Strong Resilient Population (SR) - renewal={renewal_rate}", fontsize=14)

    climate_k = climate_params[len(climate_params)//2]

    for idx, scenario in enumerate(CONFLICT_SCENARIOS.keys()):
        ax = axes[idx]

        baseline_SR = (CONFLICT_SCENARIOS[scenario]["initial_state"]["SR_RL"] +
                      CONFLICT_SCENARIOS[scenario]["initial_state"]["SR_SL"])

        x_labels = []
        final_sr = []

        for gov_strat in GOVERNMENT_STRATEGIES:
            for ag_strat in ARMED_GROUP_STRATEGIES:
                metrics = results[scenario][climate_k][renewal_rate][gov_strat][ag_strat]
                x_labels.append(f"{gov_strat[:4]}-{ag_strat[:4]}")
                final_sr.append(metrics["final_SR"])

        x = np.arange(len(x_labels))

        # Color adaptive strategy differently
        colors = ['darkgreen' if 'adap' in label else 'steelblue' for label in x_labels]
        ax.bar(x, final_sr, color=colors, alpha=0.7)
        ax.axhline(y=baseline_SR, color='red', linestyle='--', linewidth=2,
                  label=f'Baseline: {baseline_SR:.1f}')

        ax.set_title(f"{scenario.title()} Conflict (climate={climate_k})")
        ax.set_xticks(x)
        ax.set_xticklabels(x_labels, rotation=45, ha='right', fontsize=8)
        ax.set_ylabel("Final SR Population")
        ax.legend(loc='best', fontsize=8)
        ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    if save_dir:
        plt.savefig(os.path.join(save_dir, "04_strong_resilience_baseline.png"), dpi=150)
    plt.close()

    # ========================================================================
    # PLOT 5: Armed Group Population with Baseline
    # ========================================================================
    print("Generating Plot 5: Armed Group Population with Baseline...")

    fig, axes = plt.subplots(1, len(CONFLICT_SCENARIOS), figsize=(6*len(CONFLICT_SCENARIOS), 5))
    if len(CONFLICT_SCENARIOS) == 1:
        axes = [axes]
    fig.suptitle(f"Final Armed Group Population (AG) - renewal={renewal_rate}", fontsize=14)

    for idx, scenario in enumerate(CONFLICT_SCENARIOS.keys()):
        ax = axes[idx]

        baseline_AG = (CONFLICT_SCENARIOS[scenario]["initial_state"]["AG_RL"] +
                      CONFLICT_SCENARIOS[scenario]["initial_state"]["AG_SL"])

        x_labels = []
        final_ag = []

        for gov_strat in GOVERNMENT_STRATEGIES:
            for ag_strat in ARMED_GROUP_STRATEGIES:
                metrics = results[scenario][climate_k][renewal_rate][gov_strat][ag_strat]
                x_labels.append(f"{gov_strat[:4]}-{ag_strat[:4]}")
                final_ag.append(metrics["final_AG"])

        x = np.arange(len(x_labels))

        # Color adaptive strategy differently
        colors = ['darkgreen' if 'adap' in label else 'coral' for label in x_labels]
        ax.bar(x, final_ag, color=colors, alpha=0.7)
        ax.axhline(y=baseline_AG, color='red', linestyle='--', linewidth=2,
                  label=f'Baseline: {baseline_AG:.1f}')

        ax.set_title(f"{scenario.title()} Conflict (climate={climate_k})")
        ax.set_xticks(x)
        ax.set_xticklabels(x_labels, rotation=45, ha='right', fontsize=8)
        ax.set_ylabel("Final AG Population")
        ax.legend(loc='best', fontsize=8)
        ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    if save_dir:
        plt.savefig(os.path.join(save_dir, "05_armed_groups_baseline.png"), dpi=150)
    plt.close()

    # ========================================================================
    # PLOT 6: Adaptive Strategy Performance Comparison
    # ========================================================================
    print("Generating Plot 6: Adaptive Strategy Performance...")

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle("Adaptive Strategy Performance vs Other Strategies", fontsize=14)

    scenario = "medium"
    climate_k = climate_params[len(climate_params)//2]
    renewal_rate = budget_renewal_rates[len(budget_renewal_rates)//2]

    metrics_to_compare = [
        ("final_AG", "Final AG Population (lower is better)", axes[0, 0], True),
        ("final_SR", "Final SR Population (higher is better)", axes[0, 1], False),
        ("final_T", "Final Trust (higher is better)", axes[1, 0], False),
        ("total_budget_consumed", "Total Budget Consumed (lower is better)", axes[1, 1], True),
    ]

    for metric_key, title, ax, lower_better in metrics_to_compare:
        values_by_strategy = {}

        for gov_strat in GOVERNMENT_STRATEGIES:
            values = []
            for ag_strat in ARMED_GROUP_STRATEGIES:
                metrics = results[scenario][climate_k][renewal_rate][gov_strat][ag_strat]
                values.append(metrics[metric_key])
            values_by_strategy[gov_strat] = np.mean(values)

        strategies = list(values_by_strategy.keys())
        values = list(values_by_strategy.values())
        colors = ['darkgreen' if s == 'adaptive' else 'steelblue' for s in strategies]

        x = np.arange(len(strategies))
        bars = ax.bar(x, values, color=colors, alpha=0.7)

        # Highlight adaptive
        for i, s in enumerate(strategies):
            if s == 'adaptive':
                bars[i].set_edgecolor('gold')
                bars[i].set_linewidth(3)

        ax.set_title(title)
        ax.set_xticks(x)
        ax.set_xticklabels(strategies)
        ax.set_ylabel(metric_key.replace('_', ' ').title())
        ax.grid(True, alpha=0.3, axis='y')

        # Mark best performer
        if lower_better:
            best_idx = np.argmin(values)
        else:
            best_idx = np.argmax(values)
        ax.scatter(best_idx, values[best_idx], s=300, marker='*', color='gold',
                  edgecolors='black', linewidths=2, zorder=10, label='Best')
        ax.legend()

    plt.tight_layout()
    if save_dir:
        plt.savefig(os.path.join(save_dir, "06_adaptive_performance.png"), dpi=150)
    plt.close()

    # ========================================================================
    # PLOT 7: Heatmap - Strategy Effectiveness
    # ========================================================================
    print("Generating Plot 7: Strategy Effectiveness Heatmap...")

    fig, axes = plt.subplots(1, len(CONFLICT_SCENARIOS), figsize=(6*len(CONFLICT_SCENARIOS), 6))
    if len(CONFLICT_SCENARIOS) == 1:
        axes = [axes]
    fig.suptitle("Strategy Effectiveness: Final AG Population (lower is better)", fontsize=14)

    # Average across climate and renewal parameters
    matrices = {}
    global_min = float('inf')
    global_max = float('-inf')

    for scenario in CONFLICT_SCENARIOS.keys():
        matrix = np.zeros((len(GOVERNMENT_STRATEGIES), len(ARMED_GROUP_STRATEGIES)))
        for i, gov_strat in enumerate(GOVERNMENT_STRATEGIES):
            for j, ag_strat in enumerate(ARMED_GROUP_STRATEGIES):
                values = []
                for ck in climate_params:
                    for rr in budget_renewal_rates:
                        values.append(results[scenario][ck][rr][gov_strat][ag_strat]["final_AG"])
                matrix[i, j] = np.mean(values)
        matrices[scenario] = matrix
        global_min = min(global_min, matrix.min())
        global_max = max(global_max, matrix.max())

    # Plot with consistent color scale
    for idx, scenario in enumerate(CONFLICT_SCENARIOS.keys()):
        ax = axes[idx]
        matrix = matrices[scenario]

        im = ax.imshow(matrix, cmap='RdYlGn_r', aspect='auto',
                       vmin=global_min, vmax=global_max)
        ax.set_xticks(range(len(ARMED_GROUP_STRATEGIES)))
        ax.set_yticks(range(len(GOVERNMENT_STRATEGIES)))
        ax.set_xticklabels([s[:6] for s in ARMED_GROUP_STRATEGIES])
        ax.set_yticklabels([s[:6] for s in GOVERNMENT_STRATEGIES])
        ax.set_xlabel("Armed Group Strategy")
        ax.set_ylabel("Government Strategy")
        ax.set_title(f"{scenario.title()} Conflict")

        # Annotate values
        for i in range(len(GOVERNMENT_STRATEGIES)):
            for j in range(len(ARMED_GROUP_STRATEGIES)):
                color = "white" if matrix[i, j] > (global_min + global_max)/2 else "black"
                ax.text(j, i, f"{matrix[i, j]:.1f}",
                        ha="center", va="center", color=color, fontsize=9, fontweight='bold' if i == 3 else 'normal')

    fig.colorbar(im, ax=axes, label="Final AG Population", shrink=0.8)

    plt.tight_layout()
    if save_dir:
        plt.savefig(os.path.join(save_dir, "07_heatmap_effectiveness.png"), dpi=150)
    plt.close()

    print("\nAll plots generated!")


def generate_comprehensive_report(results: Dict, climate_params: List[float],
                                  budget_renewal_rates: List[float], save_dir: str = None):
    """Generate comprehensive analysis report with ALL budget scenarios."""

    n_renewals = len(budget_renewal_rates)

    if save_dir:
        os.makedirs(save_dir, exist_ok=True)

    report_lines = []
    report_lines.append("="*80)
    report_lines.append("COMPREHENSIVE CONFLICT DYNAMICS ANALYSIS REPORT")
    report_lines.append("Adaptive Strategy Performance Across All Scenarios")
    report_lines.append("="*80)
    report_lines.append("")

    # ========================================================================
    # PLOT A: Budget Renewal Rate Comparison - ALL Scenarios
    # ========================================================================
    print("\nGenerating Plot A: Budget Impact Across ALL Renewal Rates...")

    fig, axes = plt.subplots(len(CONFLICT_SCENARIOS), n_renewals,
                            figsize=(6*n_renewals, 5*len(CONFLICT_SCENARIOS)))
    if len(CONFLICT_SCENARIOS) == 1:
        axes = axes.reshape(1, -1)
    if n_renewals == 1:
        axes = axes.reshape(-1, 1)

    fig.suptitle("Strategy Performance Across ALL Budget Renewal Rates", fontsize=16)

    for row, scenario in enumerate(CONFLICT_SCENARIOS.keys()):
        for col, renewal_rate in enumerate(budget_renewal_rates):
            ax = axes[row, col]

            # Average across all climates and AG strategies for each gov strategy
            gov_scores = {}
            for gov_strat in GOVERNMENT_STRATEGIES:
                ag_values = []
                for climate_k in climate_params:
                    for ag_strat in ARMED_GROUP_STRATEGIES:
                        m = results[scenario][climate_k][renewal_rate][gov_strat][ag_strat]
                        # Composite score: lower AG, higher SR, higher T, lower budget
                        score = -m['final_AG'] + m['final_SR']*0.5 + m['final_T']*0.3 - m['total_budget_consumed']*0.01
                        ag_values.append(score)
                gov_scores[gov_strat] = np.mean(ag_values)

            strategies = list(gov_scores.keys())
            scores = list(gov_scores.values())
            colors = ['darkgreen' if s == 'adaptive' else 'steelblue' for s in strategies]

            bars = ax.bar(range(len(strategies)), scores, color=colors, alpha=0.7)

            # Highlight best
            best_idx = np.argmax(scores)
            bars[best_idx].set_edgecolor('gold')
            bars[best_idx].set_linewidth(3)

            ax.set_title(f"{scenario.title()}\nRenewal={renewal_rate:.2f}")
            ax.set_xticks(range(len(strategies)))
            ax.set_xticklabels(strategies, rotation=45, ha='right', fontsize=8)
            ax.set_ylabel("Composite Score")
            ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    if save_dir:
        plt.savefig(os.path.join(save_dir, "0A_budget_renewal_comparison_all.png"), dpi=150)
    plt.close()

    # ========================================================================
    # PLOT B: Best Strategy Per Budget Renewal Rate
    # ========================================================================
    print("Generating Plot B: Best Strategy Rankings...")

    fig, axes = plt.subplots(1, len(CONFLICT_SCENARIOS), figsize=(6*len(CONFLICT_SCENARIOS), 6))
    if len(CONFLICT_SCENARIOS) == 1:
        axes = [axes]
    fig.suptitle("Best Government Strategy by Budget Renewal Rate", fontsize=14)

    report_lines.append("\n" + "="*80)
    report_lines.append("1. BEST STRATEGY BY BUDGET RENEWAL RATE")
    report_lines.append("="*80)

    for idx, scenario in enumerate(CONFLICT_SCENARIOS.keys()):
        ax = axes[idx]

        report_lines.append(f"\n{scenario.upper()} CONFLICT:")
        report_lines.append("-" * 40)

        # For each renewal rate, find best strategy
        best_strategies = []
        best_scores = []

        for renewal_rate in budget_renewal_rates:
            gov_performance = {}
            for gov_strat in GOVERNMENT_STRATEGIES:
                metrics_list = []
                for climate_k in climate_params:
                    for ag_strat in ARMED_GROUP_STRATEGIES:
                        m = results[scenario][climate_k][renewal_rate][gov_strat][ag_strat]
                        metrics_list.append(m)

                # Compute average performance metrics
                avg_ag = np.mean([m['final_AG'] for m in metrics_list])
                avg_sr = np.mean([m['final_SR'] for m in metrics_list])
                avg_t = np.mean([m['final_T'] for m in metrics_list])
                avg_budget = np.mean([m['total_budget_consumed'] for m in metrics_list])

                # Composite score (normalized)
                score = -avg_ag/5.0 + avg_sr/50.0 + avg_t/50.0 - avg_budget/1000.0
                gov_performance[gov_strat] = {
                    'score': score,
                    'ag': avg_ag,
                    'sr': avg_sr,
                    't': avg_t,
                    'budget': avg_budget
                }

            # Find best
            best_strat = max(gov_performance.keys(), key=lambda k: gov_performance[k]['score'])
            best_strategies.append(best_strat)
            best_scores.append(gov_performance[best_strat]['score'])

            # Report
            report_lines.append(f"\nRenewal Rate: {renewal_rate:.2f}")
            for strat, perf in sorted(gov_performance.items(), key=lambda x: x[1]['score'], reverse=True):
                marker = " [*BEST*]" if strat == best_strat else ""
                report_lines.append(f"  {strat:12s}: Score={perf['score']:6.2f}, AG={perf['ag']:5.2f}, "
                                  f"SR={perf['sr']:5.1f}, T={perf['t']:5.1f}, Budget={perf['budget']:6.1f}{marker}")

        # Plot best strategy evolution
        x = range(len(budget_renewal_rates))
        colors_map = {'development': 'blue', 'security': 'red', 'balanced': 'gray', 'adaptive': 'green'}
        colors = [colors_map.get(s, 'black') for s in best_strategies]

        ax.bar(x, best_scores, color=colors, alpha=0.7)
        ax.set_title(f"{scenario.title()} Conflict")
        ax.set_xticks(x)
        ax.set_xticklabels([f"{r:.2f}" for r in budget_renewal_rates])
        ax.set_xlabel("Budget Renewal Rate")
        ax.set_ylabel("Best Strategy Score")
        ax.grid(True, alpha=0.3, axis='y')

        # Add legend
        for strat, color in colors_map.items():
            if strat in best_strategies:
                ax.bar([], [], color=color, label=strat, alpha=0.7)
        ax.legend(loc='best', fontsize=8)

    plt.tight_layout()
    if save_dir:
        plt.savefig(os.path.join(save_dir, "0B_best_strategy_by_renewal.png"), dpi=150)
    plt.close()

    # ========================================================================
    # PLOT C: Detailed Performance Metrics Grid
    # ========================================================================
    print("Generating Plot C: Detailed Performance Metrics...")

    n_renewals_actual = len(budget_renewal_rates)
    fig, axes = plt.subplots(4, n_renewals_actual, figsize=(5*n_renewals_actual, 16))
    if n_renewals_actual == 1:
        axes = axes.reshape(-1, 1)

    fig.suptitle("Detailed Performance Metrics Across Budget Renewal Rates", fontsize=16)

    metrics_to_plot = [
        ('final_AG', 'Final AG Population', 'RdYlGn', True),
        ('final_SR', 'Final SR Population', 'RdYlGn_r', False),
        ('final_T', 'Final Trust', 'RdYlGn_r', False),
        ('total_budget_consumed', 'Total Budget Consumed', 'RdYlGn', True)
    ]

    for row, (metric_key, title, cmap, lower_better) in enumerate(metrics_to_plot):
        for col, renewal_rate in enumerate(budget_renewal_rates):
            ax = axes[row, col]

            # Average across all scenarios and climates
            strategy_values = {}
            for gov_strat in GOVERNMENT_STRATEGIES:
                values = []
                for scenario in CONFLICT_SCENARIOS.keys():
                    for climate_k in climate_params:
                        for ag_strat in ARMED_GROUP_STRATEGIES:
                            m = results[scenario][climate_k][renewal_rate][gov_strat][ag_strat]
                            values.append(m[metric_key])
                strategy_values[gov_strat] = np.mean(values)

            strategies = list(strategy_values.keys())
            values = list(strategy_values.values())
            colors = ['darkgreen' if s == 'adaptive' else 'steelblue' for s in strategies]

            bars = ax.bar(range(len(strategies)), values, color=colors, alpha=0.7)

            # Highlight best
            if lower_better:
                best_idx = np.argmin(values)
            else:
                best_idx = np.argmax(values)
            bars[best_idx].set_edgecolor('gold')
            bars[best_idx].set_linewidth(3)

            if row == 0:
                ax.set_title(f"Renewal={renewal_rate:.2f}")
            if col == 0:
                ax.set_ylabel(title)
            ax.set_xticks(range(len(strategies)))
            ax.set_xticklabels(strategies, rotation=45, ha='right', fontsize=8)
            ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    if save_dir:
        plt.savefig(os.path.join(save_dir, "0C_detailed_metrics_grid.png"), dpi=150)
    plt.close()

    # ========================================================================
    # PLOT D: Overall Winner Analysis
    # ========================================================================
    print("Generating Plot D: Overall Winner Analysis...")

    report_lines.append("\n" + "="*80)
    report_lines.append("2. OVERALL STRATEGY PERFORMANCE (ALL SCENARIOS)")
    report_lines.append("="*80)

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    fig.suptitle("Overall Strategy Performance: Adaptive vs Others", fontsize=16)

    # Compute overall statistics
    overall_stats = {}
    for gov_strat in GOVERNMENT_STRATEGIES:
        all_metrics = []
        for scenario in CONFLICT_SCENARIOS.keys():
            for climate_k in climate_params:
                for renewal_rate in budget_renewal_rates:
                    for ag_strat in ARMED_GROUP_STRATEGIES:
                        all_metrics.append(results[scenario][climate_k][renewal_rate][gov_strat][ag_strat])

        overall_stats[gov_strat] = {
            'ag_mean': np.mean([m['final_AG'] for m in all_metrics]),
            'ag_std': np.std([m['final_AG'] for m in all_metrics]),
            'sr_mean': np.mean([m['final_SR'] for m in all_metrics]),
            'sr_std': np.std([m['final_SR'] for m in all_metrics]),
            't_mean': np.mean([m['final_T'] for m in all_metrics]),
            't_std': np.std([m['final_T'] for m in all_metrics]),
            'budget_mean': np.mean([m['total_budget_consumed'] for m in all_metrics]),
            'budget_std': np.std([m['total_budget_consumed'] for m in all_metrics]),
        }

    report_lines.append("\nAveraged across ALL scenarios, climates, renewals, and AG strategies:")
    report_lines.append("-" * 80)
    report_lines.append(f"{'Strategy':<15} {'AG (mean±std)':<20} {'SR (mean±std)':<20} "
                       f"{'T (mean±std)':<20} {'Budget (mean±std)':<20}")
    report_lines.append("-" * 80)

    for strat, stats in overall_stats.items():
        marker = " [*]" if strat == "adaptive" else ""
        report_lines.append(
            f"{strat:<15}{marker} "
            f"{stats['ag_mean']:5.2f}±{stats['ag_std']:4.2f}       "
            f"{stats['sr_mean']:5.1f}±{stats['sr_std']:4.1f}       "
            f"{stats['t_mean']:5.1f}±{stats['t_std']:4.1f}       "
            f"{stats['budget_mean']:6.1f}±{stats['budget_std']:5.1f}"
        )

    # Plot 1: AG population comparison
    ax = axes[0, 0]
    strategies = list(overall_stats.keys())
    ag_means = [overall_stats[s]['ag_mean'] for s in strategies]
    ag_stds = [overall_stats[s]['ag_std'] for s in strategies]
    colors = ['darkgreen' if s == 'adaptive' else 'steelblue' for s in strategies]

    ax.bar(range(len(strategies)), ag_means, yerr=ag_stds, color=colors, alpha=0.7, capsize=5)
    ax.set_title("Armed Group Population (lower is better)")
    ax.set_xticks(range(len(strategies)))
    ax.set_xticklabels(strategies)
    ax.set_ylabel("Final AG ± std")
    ax.grid(True, alpha=0.3, axis='y')

    # Plot 2: SR population comparison
    ax = axes[0, 1]
    sr_means = [overall_stats[s]['sr_mean'] for s in strategies]
    sr_stds = [overall_stats[s]['sr_std'] for s in strategies]

    ax.bar(range(len(strategies)), sr_means, yerr=sr_stds, color=colors, alpha=0.7, capsize=5)
    ax.set_title("Strong Resilient Population (higher is better)")
    ax.set_xticks(range(len(strategies)))
    ax.set_xticklabels(strategies)
    ax.set_ylabel("Final SR ± std")
    ax.grid(True, alpha=0.3, axis='y')

    # Plot 3: Trust comparison
    ax = axes[1, 0]
    t_means = [overall_stats[s]['t_mean'] for s in strategies]
    t_stds = [overall_stats[s]['t_std'] for s in strategies]

    ax.bar(range(len(strategies)), t_means, yerr=t_stds, color=colors, alpha=0.7, capsize=5)
    ax.set_title("Trust (higher is better)")
    ax.set_xticks(range(len(strategies)))
    ax.set_xticklabels(strategies)
    ax.set_ylabel("Final T ± std")
    ax.grid(True, alpha=0.3, axis='y')

    # Plot 4: Budget comparison
    ax = axes[1, 1]
    budget_means = [overall_stats[s]['budget_mean'] for s in strategies]
    budget_stds = [overall_stats[s]['budget_std'] for s in strategies]

    ax.bar(range(len(strategies)), budget_means, yerr=budget_stds, color=colors, alpha=0.7, capsize=5)
    ax.set_title("Budget Consumed (lower is better)")
    ax.set_xticks(range(len(strategies)))
    ax.set_xticklabels(strategies)
    ax.set_ylabel("Total Budget ± std")
    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    if save_dir:
        plt.savefig(os.path.join(save_dir, "0D_overall_winner_analysis.png"), dpi=150)
    plt.close()

    # ========================================================================
    # PLOT E: Win Rate Analysis
    # ========================================================================
    print("Generating Plot E: Win Rate Analysis...")

    report_lines.append("\n" + "="*80)
    report_lines.append("3. WIN RATE ANALYSIS")
    report_lines.append("="*80)

    # Count wins for each strategy
    win_counts = {s: 0 for s in GOVERNMENT_STRATEGIES}
    total_comparisons = 0

    for scenario in CONFLICT_SCENARIOS.keys():
        for climate_k in climate_params:
            for renewal_rate in budget_renewal_rates:
                # For each combination, find best strategy (avg over AG strategies)
                strategy_scores = {}
                for gov_strat in GOVERNMENT_STRATEGIES:
                    ag_values = []
                    for ag_strat in ARMED_GROUP_STRATEGIES:
                        m = results[scenario][climate_k][renewal_rate][gov_strat][ag_strat]
                        # Composite score
                        score = -m['final_AG'] + m['final_SR']*0.5 + m['final_T']*0.3 - m['total_budget_consumed']*0.01
                        ag_values.append(score)
                    strategy_scores[gov_strat] = np.mean(ag_values)

                best_strategy = max(strategy_scores.keys(), key=lambda k: strategy_scores[k])
                win_counts[best_strategy] += 1
                total_comparisons += 1

    report_lines.append(f"\nTotal scenario combinations tested: {total_comparisons}")
    report_lines.append(f"(3 scenarios × {len(climate_params)} climates × {len(budget_renewal_rates)} renewals)")
    report_lines.append("\nWin counts (best strategy per combination):")
    report_lines.append("-" * 40)

    fig, ax = plt.subplots(1, 1, figsize=(10, 6))

    strategies = list(win_counts.keys())
    wins = [win_counts[s] for s in strategies]
    win_rates = [w / total_comparisons * 100 for w in wins]
    colors = ['darkgreen' if s == 'adaptive' else 'steelblue' for s in strategies]

    bars = ax.bar(range(len(strategies)), win_rates, color=colors, alpha=0.7)

    # Add win count labels
    for i, (w, wr) in enumerate(zip(wins, win_rates)):
        ax.text(i, wr + 1, f"{w}/{total_comparisons}\n({wr:.1f}%)",
               ha='center', va='bottom', fontweight='bold')
        report_lines.append(f"  {strategies[i]:<15}: {w:3d} wins ({wr:5.1f}%)")

    # Highlight winner
    best_idx = np.argmax(wins)
    bars[best_idx].set_edgecolor('gold')
    bars[best_idx].set_linewidth(4)

    ax.set_title("Win Rate: Best Strategy Across All Scenarios", fontsize=14)
    ax.set_xticks(range(len(strategies)))
    ax.set_xticklabels(strategies, fontsize=12)
    ax.set_ylabel("Win Rate (%)", fontsize=12)
    ax.set_ylim(0, max(win_rates) * 1.15)
    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    if save_dir:
        plt.savefig(os.path.join(save_dir, "0E_win_rate_analysis.png"), dpi=150)
    plt.close()

    # ========================================================================
    # PLOT F: Time Series - Strategy Comparison (ALL BUDGET RATES)
    # ========================================================================
    print("Generating Plot F: Time Series Across ALL Budget Renewal Rates...")

    # Select medium scenario and middle climate for detailed time series
    scenario = "medium"
    climate_k = climate_params[len(climate_params)//2]

    # Create subplot grid: n_renewals columns x 6 rows (one per variable)
    fig, axes = plt.subplots(6, n_renewals, figsize=(7*n_renewals, 24))
    if n_renewals == 1:
        axes = axes.reshape(-1, 1)

    fig.suptitle(f"Time Series Evolution: {scenario.title()} Conflict, Climate={climate_k}", fontsize=16)

    variables_to_plot = [
        ('AG_total', 'Armed Groups (AG)', 'red', True),
        ('SR_total', 'Strong Resilient (SR)', 'green', False),
        ('E', 'Economy (E)', 'gold', False),
        ('T', 'Trust (T)', 'teal', False),
        ('Gov', 'Governance (Gov)', 'blue', False),
        ('V', 'Violence (V)', 'darkred', True),
    ]

    # Get baseline values for this scenario
    scenario_init = CONFLICT_SCENARIOS[scenario]["initial_state"]
    baselines = {
        'AG_total': scenario_init["AG_RL"] + scenario_init["AG_SL"],
        'SR_total': scenario_init["SR_RL"] + scenario_init["SR_SL"],
        'E': scenario_init["E"],
        'T': scenario_init["T"],
        'Gov': scenario_init["Gov"],
        'V': scenario_init["V"],
    }

    for row, (var_key, var_label, var_color, lower_better) in enumerate(variables_to_plot):
        for col, renewal_rate in enumerate(budget_renewal_rates):
            ax = axes[row, col]

            # Plot time series for each government strategy (average over AG strategies)
            for gov_strat in GOVERNMENT_STRATEGIES:
                # Average over AG strategies
                all_times = None
                all_series = []

                for ag_strat in ARMED_GROUP_STRATEGIES:
                    m = results[scenario][climate_k][renewal_rate][gov_strat][ag_strat]
                    if all_times is None:
                        all_times = m["times"]
                    all_series.append(m[var_key])

                # Compute mean across AG strategies
                mean_series = np.mean(all_series, axis=0)

                # Line style and color
                if gov_strat == "adaptive":
                    linestyle = '-'
                    linewidth = 3
                    color = 'darkgreen'
                    zorder = 10
                elif gov_strat == "development":
                    linestyle = '--'
                    linewidth = 2
                    color = 'blue'
                    zorder = 5
                elif gov_strat == "security":
                    linestyle = ':'
                    linewidth = 2
                    color = 'red'
                    zorder = 5
                else:  # balanced
                    linestyle = '-.'
                    linewidth = 2
                    color = 'gray'
                    zorder = 5

                ax.plot(all_times, mean_series, linestyle=linestyle, linewidth=linewidth,
                       color=color, label=gov_strat, alpha=0.8, zorder=zorder)

            # Add baseline reference
            baseline_val = baselines[var_key]
            ax.axhline(y=baseline_val, color='black', linestyle='--', linewidth=1,
                      alpha=0.4, label=f'Initial: {baseline_val:.1f}')

            # Formatting
            if row == 0:
                ax.set_title(f"Renewal={renewal_rate:.2f}", fontsize=11)
            if col == 0:
                ax.set_ylabel(var_label, fontsize=10)
            if row == len(variables_to_plot) - 1:
                ax.set_xlabel("Time (years)", fontsize=9)

            ax.grid(True, alpha=0.2)

            # Legend only on first plot
            if row == 0 and col == 0:
                ax.legend(loc='best', fontsize=8, ncol=2)

    plt.tight_layout()
    if save_dir:
        plt.savefig(os.path.join(save_dir, "0F_timeseries_all_budgets.png"), dpi=150)
    plt.close()

    # ========================================================================
    # PLOT G: Time Series - Scenario Comparison (ALL SCENARIOS)
    # ========================================================================
    print("Generating Plot G: Time Series Across ALL Conflict Scenarios...")

    # Use middle climate and middle renewal rate
    climate_k = climate_params[len(climate_params)//2]
    renewal_rate = budget_renewal_rates[len(budget_renewal_rates)//2]

    # Create subplot grid: 3 scenarios x 6 variables
    fig, axes = plt.subplots(6, len(CONFLICT_SCENARIOS), figsize=(7*len(CONFLICT_SCENARIOS), 24))
    if len(CONFLICT_SCENARIOS) == 1:
        axes = axes.reshape(-1, 1)

    fig.suptitle(f"Time Series Evolution Across Scenarios (climate={climate_k}, renewal={renewal_rate})", fontsize=16)

    for col, scenario in enumerate(CONFLICT_SCENARIOS.keys()):
        scenario_init = CONFLICT_SCENARIOS[scenario]["initial_state"]
        baselines_scenario = {
            'AG_total': scenario_init["AG_RL"] + scenario_init["AG_SL"],
            'SR_total': scenario_init["SR_RL"] + scenario_init["SR_SL"],
            'E': scenario_init["E"],
            'T': scenario_init["T"],
            'Gov': scenario_init["Gov"],
            'V': scenario_init["V"],
        }

        for row, (var_key, var_label, var_color, lower_better) in enumerate(variables_to_plot):
            ax = axes[row, col]

            # Plot time series for each government strategy
            for gov_strat in GOVERNMENT_STRATEGIES:
                # Average over AG strategies
                all_times = None
                all_series = []

                for ag_strat in ARMED_GROUP_STRATEGIES:
                    m = results[scenario][climate_k][renewal_rate][gov_strat][ag_strat]
                    if all_times is None:
                        all_times = m["times"]
                    all_series.append(m[var_key])

                mean_series = np.mean(all_series, axis=0)

                # Line style and color
                if gov_strat == "adaptive":
                    linestyle = '-'
                    linewidth = 3
                    color = 'darkgreen'
                    zorder = 10
                elif gov_strat == "development":
                    linestyle = '--'
                    linewidth = 2
                    color = 'blue'
                    zorder = 5
                elif gov_strat == "security":
                    linestyle = ':'
                    linewidth = 2
                    color = 'red'
                    zorder = 5
                else:  # balanced
                    linestyle = '-.'
                    linewidth = 2
                    color = 'gray'
                    zorder = 5

                ax.plot(all_times, mean_series, linestyle=linestyle, linewidth=linewidth,
                       color=color, label=gov_strat, alpha=0.8, zorder=zorder)

            # Add baseline reference
            baseline_val = baselines_scenario[var_key]
            ax.axhline(y=baseline_val, color='black', linestyle='--', linewidth=1,
                      alpha=0.4, label=f'Initial: {baseline_val:.1f}')

            # Formatting
            if row == 0:
                ax.set_title(f"{scenario.title()} Conflict", fontsize=11)
            if col == 0:
                ax.set_ylabel(var_label, fontsize=10)
            if row == len(variables_to_plot) - 1:
                ax.set_xlabel("Time (years)", fontsize=9)

            ax.grid(True, alpha=0.2)

            # Legend only on first plot
            if row == 0 and col == 0:
                ax.legend(loc='best', fontsize=8, ncol=2)

    plt.tight_layout()
    if save_dir:
        plt.savefig(os.path.join(save_dir, "0G_timeseries_all_scenarios.png"), dpi=150)
    plt.close()

    # ========================================================================
    # PLOT H: Adaptive Strategy Mode Tracking (if debug enabled)
    # ========================================================================
    print("Generating Plot H: Budget Consumed vs Generated Over Time...")

    # Show budget evolution over time for adaptive vs others
    scenario = "medium"
    climate_k = climate_params[len(climate_params)//2]
    renewal_rate = budget_renewal_rates[len(budget_renewal_rates)//2]

    fig, axes = plt.subplots(2, 1, figsize=(14, 10))
    fig.suptitle(f"Budget Evolution: {scenario.title()}, climate={climate_k}, renewal={renewal_rate}", fontsize=14)

    # Plot 1: Cumulative budget consumed
    ax = axes[0]
    for gov_strat in GOVERNMENT_STRATEGIES:
        # Average over AG strategies
        all_times = None
        e_consumed_list = []
        gov_consumed_list = []

        for ag_strat in ARMED_GROUP_STRATEGIES:
            m = results[scenario][climate_k][renewal_rate][gov_strat][ag_strat]
            if all_times is None:
                all_times = m["times"]
            e_consumed_list.append(m["E_consumed"])
            gov_consumed_list.append(m["Gov_consumed"])

        # Mean across AG strategies
        total_consumed = np.mean(e_consumed_list, axis=0) + np.mean(gov_consumed_list, axis=0)

        # Line style
        if gov_strat == "adaptive":
            linestyle = '-'
            linewidth = 3
            color = 'darkgreen'
        elif gov_strat == "development":
            linestyle = '--'
            linewidth = 2
            color = 'blue'
        elif gov_strat == "security":
            linestyle = ':'
            linewidth = 2
            color = 'red'
        else:
            linestyle = '-.'
            linewidth = 2
            color = 'gray'

        ax.plot(all_times, total_consumed, linestyle=linestyle, linewidth=linewidth,
               color=color, label=gov_strat, alpha=0.8)

    ax.set_ylabel("Cumulative Budget Consumed (E+Gov)", fontsize=11)
    ax.set_xlabel("Time (years)", fontsize=11)
    ax.legend(loc='best', fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_title("Total Budget Consumed Over Time")

    # Plot 2: Cumulative budget generated
    ax = axes[1]
    for gov_strat in GOVERNMENT_STRATEGIES:
        all_times = None
        e_generated_list = []
        gov_generated_list = []

        for ag_strat in ARMED_GROUP_STRATEGIES:
            m = results[scenario][climate_k][renewal_rate][gov_strat][ag_strat]
            if all_times is None:
                all_times = m["times"]
            e_generated_list.append(m["E_generated"])
            gov_generated_list.append(m["Gov_generated"])

        total_generated = np.mean(e_generated_list, axis=0) + np.mean(gov_generated_list, axis=0)

        # Line style
        if gov_strat == "adaptive":
            linestyle = '-'
            linewidth = 3
            color = 'darkgreen'
        elif gov_strat == "development":
            linestyle = '--'
            linewidth = 2
            color = 'blue'
        elif gov_strat == "security":
            linestyle = ':'
            linewidth = 2
            color = 'red'
        else:
            linestyle = '-.'
            linewidth = 2
            color = 'gray'

        ax.plot(all_times, total_generated, linestyle=linestyle, linewidth=linewidth,
               color=color, label=gov_strat, alpha=0.8)

    ax.set_ylabel("Cumulative Budget Generated (E+Gov)", fontsize=11)
    ax.set_xlabel("Time (years)", fontsize=11)
    ax.legend(loc='best', fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_title("Total Budget Generated Over Time")

    plt.tight_layout()
    if save_dir:
        plt.savefig(os.path.join(save_dir, "0H_budget_evolution_timeseries.png"), dpi=150)
    plt.close()

    # ========================================================================
    # Save Report
    # ========================================================================
    report_lines.append("\n" + "="*80)
    report_lines.append("END OF REPORT")
    report_lines.append("="*80)

    if save_dir:
        report_path = os.path.join(save_dir, "COMPREHENSIVE_REPORT.txt")
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write('\n'.join(report_lines))
        print(f"\n[OK] Comprehensive report saved to: {report_path}")

    # Print report to console
    print("\n" + '\n'.join(report_lines))

    return '\n'.join(report_lines)


def plot_process_analysis(process_analyses: Dict, save_dir: str = None):
    """
    Generate visualizations for process-based temporal analysis.

    Args:
        process_analyses: Dictionary with process analysis results from batch_process_analysis
        save_dir: Directory to save plots
    """
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)

    print("\n" + "="*70)
    print("GENERATING PROCESS ANALYSIS VISUALIZATIONS")
    print("="*70)

    # ========================================================================
    # PLOT I: Process Mode Evolution Over Time
    # ========================================================================
    print("\nGenerating Plot I: Process Mode Evolution...")

    fig, axes = plt.subplots(2, len(process_analyses), figsize=(6*len(process_analyses), 10))
    if len(process_analyses) == 1:
        axes = axes.reshape(-1, 1)

    fig.suptitle("Process Mode Evolution: Peace vs Conflict Categories", fontsize=16)

    mode_colors = {
        'Stationary Mode': 'green',
        'Overproduction Mode': 'blue',
        'Challenge Mode': 'orange',
        'Problem Mode': 'red'
    }

    for col, (gov_strat, analysis) in enumerate(process_analyses.items()):
        # Peace category modes over time
        ax = axes[0, col]
        peace_modes = analysis['single_step_analysis']['peace_modes']
        time_points_raw = analysis['time_points']

        # Ensure time_points and modes have same length
        # (time_points might have extra element for initial state)
        time_points = time_points_raw[:len(peace_modes)]

        # Create mode transition visualization
        mode_indices = [list(mode_colors.keys()).index(m) for m in peace_modes]
        ax.scatter(time_points, mode_indices, c=[mode_colors[m] for m in peace_modes],
                  alpha=0.6, s=20)

        ax.set_ylabel("Process Mode (Peace)")
        ax.set_yticks(range(len(mode_colors)))
        ax.set_yticklabels(list(mode_colors.keys()), fontsize=8)
        ax.set_title(f"{gov_strat.title()} Strategy", fontsize=12)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, time_points[-1] if time_points else 1)

        # Conflict category modes over time
        ax = axes[1, col]
        conflict_modes = analysis['single_step_analysis']['conflict_modes']

        mode_indices = [list(mode_colors.keys()).index(m) for m in conflict_modes]
        ax.scatter(time_points, mode_indices, c=[mode_colors[m] for m in conflict_modes],
                  alpha=0.6, s=20)

        ax.set_ylabel("Process Mode (Conflict)")
        ax.set_xlabel("Time (years)")
        ax.set_yticks(range(len(mode_colors)))
        ax.set_yticklabels(list(mode_colors.keys()), fontsize=8)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, time_points[-1])

    plt.tight_layout()
    if save_dir:
        plt.savefig(os.path.join(save_dir, "0I_process_mode_evolution.png"), dpi=150)
    plt.close()

    # ========================================================================
    # PLOT J: Multi-Scale Process Mode Distribution
    # ========================================================================
    print("Generating Plot J: Multi-Scale Process Mode Distribution...")

    fig, axes = plt.subplots(2, len(process_analyses), figsize=(6*len(process_analyses), 10))
    if len(process_analyses) == 1:
        axes = axes.reshape(-1, 1)

    fig.suptitle("Multi-Scale Process Mode Distribution", fontsize=16)

    for col, (gov_strat, analysis) in enumerate(process_analyses.items()):
        multiscale = analysis['multiscale_analysis']

        # Peace category multi-scale
        ax = axes[0, col]
        scales = sorted(multiscale.keys())
        mode_fractions = {mode: [] for mode in mode_colors.keys()}

        for scale in scales:
            scale_data = multiscale[scale]
            peace_dist = scale_data['mode_distribution']['peace']
            total = sum(peace_dist.values())
            if total > 0:
                for mode in mode_colors.keys():
                    mode_fractions[mode].append(peace_dist.get(mode, 0) / total)
            else:
                for mode in mode_colors.keys():
                    mode_fractions[mode].append(0)

        # Stacked bar chart
        bottom = np.zeros(len(scales))
        for mode, fractions in mode_fractions.items():
            ax.bar(range(len(scales)), fractions, bottom=bottom,
                  label=mode, color=mode_colors[mode], alpha=0.8)
            bottom += fractions

        ax.set_xlabel("Time Window Size")
        ax.set_ylabel("Mode Fraction (Peace)")
        ax.set_title(f"{gov_strat.title()} Strategy", fontsize=12)
        ax.set_xticks(range(len(scales)))
        # Rotate labels if too many scales
        rotation = 45 if len(scales) > 6 else 0
        ax.set_xticklabels([f"{s}" for s in scales], rotation=rotation, ha='right' if rotation else 'center', fontsize=8)
        if col == 0:
            ax.legend(loc='best', fontsize=7)
        ax.grid(True, alpha=0.3, axis='y')

        # Conflict category multi-scale
        ax = axes[1, col]
        mode_fractions = {mode: [] for mode in mode_colors.keys()}

        for scale in scales:
            scale_data = multiscale[scale]
            conflict_dist = scale_data['mode_distribution']['conflict']
            total = sum(conflict_dist.values())
            if total > 0:
                for mode in mode_colors.keys():
                    mode_fractions[mode].append(conflict_dist.get(mode, 0) / total)
            else:
                for mode in mode_colors.keys():
                    mode_fractions[mode].append(0)

        # Stacked bar chart
        bottom = np.zeros(len(scales))
        for mode, fractions in mode_fractions.items():
            ax.bar(range(len(scales)), fractions, bottom=bottom,
                  label=mode, color=mode_colors[mode], alpha=0.8)
            bottom += fractions

        ax.set_xlabel("Time Window Size")
        ax.set_ylabel("Mode Fraction (Conflict)")
        ax.set_xticks(range(len(scales)))
        # Rotate labels if too many scales
        rotation = 45 if len(scales) > 6 else 0
        ax.set_xticklabels([f"{s}" for s in scales], rotation=rotation, ha='right' if rotation else 'center', fontsize=8)
        ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    if save_dir:
        plt.savefig(os.path.join(save_dir, "0J_multiscale_mode_distribution.png"), dpi=150)
    plt.close()

    # ========================================================================
    # PLOT K: Cognitive Domain Fraction
    # ========================================================================
    print("Generating Plot K: Cognitive Domain Fraction...")

    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    fig.suptitle("Fraction of Time in Cognitive Domain (Self-Maintaining Processes)", fontsize=14)

    strategies = list(process_analyses.keys())
    peace_fractions = [analysis['single_step_analysis']['cognitive_domain_fraction']['peace']
                      for analysis in process_analyses.values()]
    conflict_fractions = [analysis['single_step_analysis']['cognitive_domain_fraction']['conflict']
                         for analysis in process_analyses.values()]

    x = np.arange(len(strategies))
    width = 0.35

    bars1 = ax.bar(x - width/2, peace_fractions, width, label='Peace Category',
                   color='skyblue', alpha=0.8)
    bars2 = ax.bar(x + width/2, conflict_fractions, width, label='Conflict Category',
                   color='coral', alpha=0.8)

    # Add value labels on bars
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{height:.1%}',
                   ha='center', va='bottom', fontsize=9)

    ax.set_ylabel("Fraction of Time in Cognitive Domain", fontsize=12)
    ax.set_xlabel("Government Strategy", fontsize=12)
    ax.set_xticks(x)
    ax.set_xticklabels(strategies, fontsize=11)
    ax.legend(fontsize=10)
    ax.set_ylim(0, 1.1)
    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    if save_dir:
        plt.savefig(os.path.join(save_dir, "0K_cognitive_domain_fraction.png"), dpi=150)
    plt.close()

    # ========================================================================
    # PLOT L: Seasonal Resource Recovery Effects
    # ========================================================================
    print("Generating Plot L: Seasonal Resource Recovery Effects...")

    # Get time series data from one representative analysis
    representative_analysis = list(process_analyses.values())[0]
    time_points = representative_analysis['time_points']

    # Compute seasonal modifier over time
    import math
    seasonal_amplitude = 0.3
    seasonal_period = 1.0
    # Create seasonal modifiers for each climate parameter
    seasonal_modifiers_by_climate = {}
    for climate_amplitude in CLIMATE_PARAMETERS:
        seasonal_modifiers = [
            (climate_amplitude / 2.0) * (1.0 + math.cos(2.0 * math.pi * t / seasonal_period))
            for t in time_points
        ]
        seasonal_modifiers_by_climate[climate_amplitude] = seasonal_modifiers

    fig, axes = plt.subplots(2, 1, figsize=(14, 10))
    fig.suptitle("Seasonal Variation Effects on Resource Recovery", fontsize=16)

    # Plot 1: Seasonal modifier over time (use a representative climate parameter)
    ax = axes[0]
    representative_climate = CLIMATE_PARAMETERS[len(CLIMATE_PARAMETERS)//2]  # Use middle value
    seasonal_modifiers = seasonal_modifiers_by_climate[representative_climate]
    ax.plot(time_points, seasonal_modifiers, linewidth=2, color='darkgreen', label='Seasonal Modifier')
    ax.axhline(y=0, color='black', linestyle='-', linewidth=1, alpha=0.3)
    ax.set_ylabel("Resource Recovery Rate", fontsize=12)
    ax.set_xlabel("Time (years)", fontsize=12)
    ax.set_title(f"Seasonal Variation in Resource Recovery (Climate={representative_climate:.3f})", fontsize=12)
    ax.legend(loc='upper right', fontsize=10)
    ax.grid(True, alpha=0.3)

    # Plot 2: Compare strategies with seasonal effects
    ax = axes[1]

    # Extract land restoration (RL) and economy (E) for each strategy
    for gov_strat, analysis in process_analyses.items():
        # We'll need to extract state history from simulator data
        # For now, show conceptual seasonal impact on different metrics
        # (This would require storing state_history in process_analyses)
        pass

    # Placeholder: Show seasonal impact on key metrics conceptually
    ax.text(0.5, 0.5, "Seasonal effects on resource recovery:\n\n" +
            "• Land restoration (r1, r2): ±30% variation\n" +
            "• Economic production (r9, r10): ±30% variation\n" +
            "• Peak recovery in summer (t=0, 1, 2,...)\n" +
            "• Minimum recovery in winter (t=0.5, 1.5, 2.5,...)\n\n" +
            "This periodic variation represents seasonal climate impacts\n" +
            "on agricultural productivity and land use dynamics.",
            transform=ax.transAxes, fontsize=11, verticalalignment='center',
            horizontalalignment='center', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    plt.tight_layout()
    if save_dir:
        plt.savefig(os.path.join(save_dir, "0L_seasonal_variation.png"), dpi=150)
    plt.close()

    # ========================================================================
    # PLOT M: Steady State Mechanism (Production vs Consumption Balance)
    # ========================================================================
    print("Generating Plot M: Steady State Mechanism Analysis...")

    # Check if steady state analysis data is available
    has_steady_state_data = any('steady_state_analysis' in analysis
                                for analysis in process_analyses.values())

    if has_steady_state_data:
        variables_analyzed = ['V', 'T', 'Gov']
        n_vars = len(variables_analyzed)
        n_strats = len(process_analyses)

        fig, axes = plt.subplots(n_vars, n_strats, figsize=(5*n_strats, 4*n_vars))
        if n_strats == 1:
            axes = axes.reshape(-1, 1)

        fig.suptitle("Production vs Consumption Balance (Steady State Mechanism)", fontsize=16)

        for col, (gov_strat, analysis) in enumerate(process_analyses.items()):
            if 'steady_state_analysis' not in analysis:
                continue

            for row, var in enumerate(variables_analyzed):
                ax = axes[row, col]
                ss_analysis = analysis['steady_state_analysis'][var]

                # Extract production and consumption data
                producers = ss_analysis['producer_rates']
                consumers = ss_analysis['consumer_rates']

                # Create data for stacked bar chart
                prod_labels = list(producers.keys())
                prod_values = [producers[k] for k in prod_labels]
                cons_labels = list(consumers.keys())
                cons_values = [consumers[k] for k in cons_labels]

                # Plot production (positive) and consumption (negative)
                x_pos = 0
                y_offset_prod = 0
                y_offset_cons = 0

                # Production bars (stacked upward)
                for i, (label, value) in enumerate(zip(prod_labels, prod_values)):
                    if value > 1e-6:  # Only show active reactions
                        color = f'C{i % 10}'
                        ax.bar(x_pos, value, bottom=y_offset_prod, color=color,
                              alpha=0.7, label=f"{label} (prod)", width=0.35)
                        y_offset_prod += value

                # Consumption bars (stacked downward)
                for i, (label, value) in enumerate(zip(cons_labels, cons_values)):
                    if value > 1e-6:  # Only show active reactions
                        color = f'C{(i+len(prod_labels)) % 10}'
                        ax.bar(x_pos + 0.4, -value, bottom=-y_offset_cons, color=color,
                              alpha=0.7, label=f"{label} (cons)", width=0.35)
                        y_offset_cons += value

                # Add net rate line
                net_rate = ss_analysis['net_rate']
                ax.axhline(y=0, color='black', linestyle='-', linewidth=2, alpha=0.5)
                ax.axhline(y=net_rate, color='red', linestyle='--', linewidth=2,
                          label=f"Net: {net_rate:.4f}")

                # Formatting
                ax.set_ylabel(f"{var} Rate", fontsize=10)
                if row == 0:
                    ax.set_title(f"{gov_strat.title()}\n{'[STEADY STATE]' if ss_analysis['is_steady_state'] else '[DYNAMIC]'}",
                                fontsize=11, fontweight='bold' if ss_analysis['is_steady_state'] else 'normal')
                ax.set_xticks([x_pos, x_pos + 0.4])
                ax.set_xticklabels(['Production', 'Consumption'], fontsize=9)
                ax.legend(loc='best', fontsize=6, ncol=2)
                ax.grid(True, alpha=0.3, axis='y')

                # Highlight if steady state
                if ss_analysis['is_steady_state']:
                    ax.set_facecolor('#f0f0f0')

        plt.tight_layout()
        if save_dir:
            plt.savefig(os.path.join(save_dir, "0M_steady_state_mechanism.png"), dpi=150)
        plt.close()

    print("\nProcess analysis visualizations complete!")


def generate_comprehensive_individual_plots(results: Dict, climate_params: List[float],
                                            budget_renewal_rates: List[float],
                                            save_dir: str = None):
    """
    Generate individual plots for EVERY combination of scenario, climate, and renewal rate.
    This ensures all 45 combinations (3 scenarios × 5 climates × 3 renewals) are visualized.
    """
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
        # Create subdirectory for individual combination plots
        combo_dir = os.path.join(save_dir, "individual_combinations")
        os.makedirs(combo_dir, exist_ok=True)
    else:
        return

    print("\n" + "="*70)
    print("GENERATING COMPREHENSIVE INDIVIDUAL PLOTS FOR ALL COMBINATIONS")
    print("="*70)

    total_combos = len(CONFLICT_SCENARIOS) * len(climate_params) * len(budget_renewal_rates)
    combo_count = 0

    for scenario in CONFLICT_SCENARIOS.keys():
        for climate_k in climate_params:
            for renewal_rate in budget_renewal_rates:
                combo_count += 1
                print(f"  Generating plot {combo_count}/{total_combos}: "
                      f"{scenario}, climate={climate_k:.3f}, renewal={renewal_rate:.2f}")

                # Create figure with 3x3 grid showing key metrics for all strategy combinations
                fig, axes = plt.subplots(3, 3, figsize=(20, 15))
                fig.suptitle(f"{scenario.title()} Conflict | Climate={climate_k:.3f} | "
                           f"Budget Renewal={renewal_rate:.2f}", fontsize=16, fontweight='bold')

                # Variables to plot
                variables = [
                    ('final_AG', 'Final Armed Groups (AG)', 'red', True),
                    ('final_SR', 'Final Strong Resilient (SR)', 'green', False),
                    ('final_E', 'Final Economy (E)', 'gold', False),
                    ('final_T', 'Final Trust (T)', 'teal', False),
                    ('final_Gov', 'Final Governance', 'blue', False),
                    ('final_V', 'Final Violence (V)', 'darkred', True),
                    ('total_budget_consumed', 'Total Budget Consumed', 'purple', True),
                    ('total_budget_generated', 'Total Budget Generated', 'cyan', False),
                    ('net_budget', 'Net Budget (Gen - Cons)', 'navy', False),
                ]

                for idx, (metric_key, ylabel, color, lower_better) in enumerate(variables):
                    row = idx // 3
                    col = idx % 3
                    ax = axes[row, col]

                    # Collect data for all strategy combinations
                    x_labels = []
                    values = []

                    for gov_strat in GOVERNMENT_STRATEGIES:
                        for ag_strat in ARMED_GROUP_STRATEGIES:
                            try:
                                metrics = results[scenario][climate_k][renewal_rate][gov_strat][ag_strat]
                                x_labels.append(f"{gov_strat[:3]}-{ag_strat[:3]}")
                                values.append(metrics[metric_key])
                            except KeyError:
                                continue

                    if values:
                        x = np.arange(len(x_labels))
                        bars = ax.bar(x, values, color=color, alpha=0.7)

                        # Highlight adaptive strategy
                        for i, label in enumerate(x_labels):
                            if 'ada' in label.lower():
                                bars[i].set_color('darkgreen')
                                bars[i].set_alpha(0.9)
                                bars[i].set_edgecolor('black')
                                bars[i].set_linewidth(2)

                        # Mark best/worst
                        if lower_better:
                            best_idx = np.argmin(values)
                            worst_idx = np.argmax(values)
                        else:
                            best_idx = np.argmax(values)
                            worst_idx = np.argmin(values)

                        bars[best_idx].set_edgecolor('lime')
                        bars[best_idx].set_linewidth(3)
                        bars[worst_idx].set_edgecolor('red')
                        bars[worst_idx].set_linewidth(3)

                        ax.set_xticks(x)
                        ax.set_xticklabels(x_labels, rotation=45, ha='right', fontsize=8)
                        ax.set_ylabel(ylabel, fontsize=10)
                        ax.grid(True, alpha=0.3, axis='y')
                        ax.set_title(ylabel, fontsize=11, fontweight='bold')

                        # Add value labels on bars
                        for i, (bar, val) in enumerate(zip(bars, values)):
                            height = bar.get_height()
                            ax.text(bar.get_x() + bar.get_width()/2., height,
                                  f'{val:.1f}',
                                  ha='center', va='bottom', fontsize=7)

                plt.tight_layout()

                # Save with descriptive filename
                filename = f"combo_{scenario}_{climate_k:.3f}_{renewal_rate:.2f}.png"
                filepath = os.path.join(combo_dir, filename)
                plt.savefig(filepath, dpi=150, bbox_inches='tight')
                plt.close()

    print(f"\n[OK] Generated {combo_count} comprehensive individual plots in: {combo_dir}/")
    print()


def generate_comprehensive_time_series_plots(results: Dict, climate_params: List[float],
                                             budget_renewal_rates: List[float],
                                             save_dir: str = None):
    """
    Generate time series plots for EVERY combination showing temporal evolution.
    This creates detailed time series for each of the 45 combinations.
    """
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
        # Create subdirectory for time series plots
        ts_dir = os.path.join(save_dir, "time_series_all")
        os.makedirs(ts_dir, exist_ok=True)
    else:
        return

    print("\n" + "="*70)
    print("GENERATING COMPREHENSIVE TIME SERIES FOR ALL COMBINATIONS")
    print("="*70)

    total_combos = len(CONFLICT_SCENARIOS) * len(climate_params) * len(budget_renewal_rates)
    combo_count = 0

    for scenario in CONFLICT_SCENARIOS.keys():
        for climate_k in climate_params:
            for renewal_rate in budget_renewal_rates:
                combo_count += 1
                print(f"  Generating time series {combo_count}/{total_combos}: "
                      f"{scenario}, climate={climate_k:.3f}, renewal={renewal_rate:.2f}")

                # Create figure with 2x3 grid for key variables
                fig, axes = plt.subplots(2, 3, figsize=(20, 12))
                fig.suptitle(f"Time Series: {scenario.title()} | Climate={climate_k:.3f} | "
                           f"Renewal={renewal_rate:.2f}", fontsize=16, fontweight='bold')

                variables_to_plot = [
                    ('AG_total', 'Armed Groups (AG)', 'red'),
                    ('SR_total', 'Strong Resilient (SR)', 'green'),
                    ('E', 'Economy (E)', 'gold'),
                    ('T', 'Trust (T)', 'teal'),
                    ('V', 'Violence (V)', 'darkred'),
                    ('Gov', 'Governance', 'blue'),
                ]

                for idx, (var_key, ylabel, color) in enumerate(variables_to_plot):
                    row = idx // 3
                    col = idx % 3
                    ax = axes[row, col]

                    # Plot all government strategies
                    for gov_strat in GOVERNMENT_STRATEGIES:
                        # Average over AG strategies for cleaner visualization
                        all_series = []
                        all_times = None

                        for ag_strat in ARMED_GROUP_STRATEGIES:
                            try:
                                m = results[scenario][climate_k][renewal_rate][gov_strat][ag_strat]
                                if all_times is None:
                                    all_times = m["times"]
                                all_series.append(m[var_key])
                            except KeyError:
                                continue

                        if all_series:
                            mean_series = np.mean(all_series, axis=0)
                            linestyle = '-' if gov_strat != 'adaptive' else '-'
                            linewidth = 3 if gov_strat == 'adaptive' else 2
                            alpha = 1.0 if gov_strat == 'adaptive' else 0.7

                            ax.plot(all_times, mean_series, label=gov_strat,
                                  linestyle=linestyle, linewidth=linewidth, alpha=alpha)

                    # Add baseline
                    if var_key in ['AG_total', 'SR_total', 'E', 'T', 'V', 'Gov']:
                        scenario_init = CONFLICT_SCENARIOS[scenario]["initial_state"]
                        if var_key == 'AG_total':
                            baseline = scenario_init["AG_RL"] + scenario_init["AG_SL"]
                        elif var_key == 'SR_total':
                            baseline = scenario_init["SR_RL"] + scenario_init["SR_SL"]
                        else:
                            baseline = scenario_init.get(var_key[0] if len(var_key) == 1 else var_key, 0)

                        ax.axhline(y=baseline, color='black', linestyle='--', linewidth=1.5,
                                 alpha=0.5, label=f'Initial: {baseline:.1f}')

                    ax.set_xlabel("Time (years)", fontsize=10)
                    ax.set_ylabel(ylabel, fontsize=10)
                    ax.set_title(ylabel, fontsize=11, fontweight='bold')
                    ax.legend(loc='best', fontsize=8)
                    ax.grid(True, alpha=0.3)

                plt.tight_layout()

                # Save with descriptive filename
                filename = f"ts_{scenario}_{climate_k:.3f}_{renewal_rate:.2f}.png"
                filepath = os.path.join(ts_dir, filename)
                plt.savefig(filepath, dpi=150, bbox_inches='tight')
                plt.close()

    print(f"\n[OK] Generated {combo_count} comprehensive time series plots in: {ts_dir}/")
    print()


def main():
    """Main entry point for batch simulation."""
    print("="*70)
    print("CONFLICT DYNAMICS BATCH SIMULATOR - ADAPTIVE")
    print("="*70)

    # Use first 10 climate parameters (or fewer if list is shorter)
    climate_params = CLIMATE_PARAMETERS[:min(10, len(CLIMATE_PARAMETERS))]

    print("\nThis will run simulations across:")
    print(f"  - {len(CONFLICT_SCENARIOS)} conflict scenarios: {list(CONFLICT_SCENARIOS.keys())}")
    print(f"  - {len(climate_params)} climate parameters: {climate_params}")
    print(f"  - {len(BUDGET_RENEWAL_RATES)} budget renewal rates: {BUDGET_RENEWAL_RATES}")
    print(f"  - {len(GOVERNMENT_STRATEGIES)} government strategies: {GOVERNMENT_STRATEGIES}")
    print(f"  - {len(ARMED_GROUP_STRATEGIES)} armed group strategies: {ARMED_GROUP_STRATEGIES}")

    total_sims = (len(CONFLICT_SCENARIOS) * len(climate_params) * len(BUDGET_RENEWAL_RATES) *
                  len(GOVERNMENT_STRATEGIES) * len(ARMED_GROUP_STRATEGIES))
    print(f"\nTotal: {total_sims} simulations, each for 20 years")

    # Load the model
    project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
    model_file = os.path.join(project_root, 'networks/Conflict_Theory/Resource_Community_Insurgency_Loops_model3.txt')

    # Run batch simulation
    results = run_batch_simulation(model_file, climate_params, BUDGET_RENEWAL_RATES,
                                   total_time=20.0, dt=0.1)

    # Generate comprehensive report FIRST (most important)
    save_dir = "visualizations/conflict_simulation_results_adaptive"
    print("\n" + "="*70)
    print("GENERATING COMPREHENSIVE REPORT WITH ALL BUDGET SCENARIOS")
    print("="*70)
    generate_comprehensive_report(results, climate_params, BUDGET_RENEWAL_RATES, save_dir=save_dir)

    # Generate standard comparison plots
    print("\n" + "="*70)
    print("GENERATING ADDITIONAL COMPARISON PLOTS")
    print("="*70)
    plot_comparison_results(results, climate_params, BUDGET_RENEWAL_RATES, save_dir=save_dir)

    # Perform process-based temporal analysis
    process_analyses = batch_process_analysis(results, model_file, climate_params, BUDGET_RENEWAL_RATES)

    # Generate process analysis visualizations
    plot_process_analysis(process_analyses, save_dir=save_dir)

    # Generate comprehensive individual plots for ALL combinations
    generate_comprehensive_individual_plots(results, climate_params, BUDGET_RENEWAL_RATES, save_dir=save_dir)

    # Generate comprehensive time series plots for ALL combinations
    generate_comprehensive_time_series_plots(results, climate_params, BUDGET_RENEWAL_RATES, save_dir=save_dir)

    # Print summary
    print("\n" + "="*70)
    print("SUMMARY: Adaptive Strategy Performance")
    print("="*70)

    for scenario in CONFLICT_SCENARIOS.keys():
        print(f"\n{scenario.upper()} CONFLICT:")

        # Compare adaptive vs other strategies
        climate_k = climate_params[len(climate_params)//2]
        renewal_rate = BUDGET_RENEWAL_RATES[len(BUDGET_RENEWAL_RATES)//2]

        print(f"  Parameters: climate={climate_k}, renewal={renewal_rate}")

        for gov_strat in GOVERNMENT_STRATEGIES:
            ag_values = []
            budget_values = []
            sr_values = []
            t_values = []
            for ag_strat in ARMED_GROUP_STRATEGIES:
                m = results[scenario][climate_k][renewal_rate][gov_strat][ag_strat]
                ag_values.append(m['final_AG'])
                budget_values.append(m['total_budget_consumed'])
                sr_values.append(m['final_SR'])
                t_values.append(m['final_T'])

            avg_ag = np.mean(ag_values)
            avg_budget = np.mean(budget_values)
            avg_sr = np.mean(sr_values)
            avg_t = np.mean(t_values)
            marker = " [*]" if gov_strat == "adaptive" else ""
            print(f"    {gov_strat:12s}{marker}: AG={avg_ag:5.2f}, Budget={avg_budget:6.1f}, SR={avg_sr:5.1f}, T={avg_t:5.1f}")

    print("\n" + "="*70)
    print(f"Results saved to: {save_dir}/")
    print("="*70)


if __name__ == "__main__":
    main()
