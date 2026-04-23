#!/usr/bin/env python3
"""
aps_lib_stochastic.py  —  APS Stochastic Model Library (paper-specific layer)
==============================================================================

Encodes the Lake Chad conflict model (Conflict_Model_Article_2_v1.txt) as
data structures consumed by pyCOT's BudgetAllocationSimulator engine.

Nothing here knows about ODE integration or strategy parameters.
The only job of this module is:
  1. Define REACTION_DEFS   — topology + consumed/enabling/net_change/seasonal
  2. Define DEFAULT_KINETICS — k and K values for each reaction
  3. Define REGIME_IC_RANGES — initial condition ranges for regimes R1-R6
  4. Expose make_simulator() — returns a ready-to-run BudgetAllocationSimulator

Kinetic parameter guide
-----------------------
k_r : base rate constant [events / year / (units of consumed species)]
      Controls how frequently reaction r fires per unit of its reactants.
      Tune relative to other reactions of the same "timescale group".

K   : saturation threshold for a specific species in reaction r.
      When [X] >> K_X the reaction saturates (capacity ceiling).
      When [X] << K_X the reaction is approximately linear in [X].
      Reasonable starting value ≈ expected [X] at regime boundary.

dt  : time step.  Must satisfy  sum_r p(X_i -> r) <= 1  for all species.
      Default 0.05 years (~18 days).  Reduce if warnings appear.

climate : float in [0,1].  0=very dry, 0.5=neutral, 1=very wet.
          Scales amplitude of seasonal reactions (r1-r4).
"""

import os
import sys
import warnings
from typing import Dict, Tuple, Optional

# ── pyCOT engine ─────────────────────────────────────────────────────────────
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT  = os.path.normpath(os.path.join(_SCRIPT_DIR, '..', '..', '..'))
sys.path.insert(0, os.path.join(_REPO_ROOT, 'src'))

from pyCOT.kinetics.stochastic_events import BudgetAllocationSimulator

# ─────────────────────────────────────────────────────────────────────────────
# SPECIES
# ─────────────────────────────────────────────────────────────────────────────
SPECIES = [
    'SR_RL', 'SR_SL',       # Resilient population on good / degraded land
    'WR_RL', 'WR_SL',       # Vulnerable population on good / degraded land
    'AG_RL', 'AG_SL',       # Armed group on good / degraded land
    'RL',    'SL',           # Resourceful / Scarce land (area units)
    'E_peace', 'E_conflict', # Peace / conflict economy stocks
    'T',     'V',    'Gov',  # Trust, Violence, Governance
]

# ─────────────────────────────────────────────────────────────────────────────
# REACTION DEFINITIONS
# ─────────────────────────────────────────────────────────────────────────────
# Each entry: {
#   'consumed' : {species: stoich}   draws budget; enters propensity
#   'enabling' : [species, ...]      catalyst; enters propensity only (no budget)
#                                    repeat species for higher power
#   'net_change': {species: delta}   applied as delta * n_firings to state
#   'seasonal'  : None | 'recovery' | 'degradation'
# }
#
# Rule: a species with net_change == 0 (true catalyst) goes in 'enabling'.
#       a species with net_change > 0 (net produced)   also goes in 'enabling'
#       if it appears on the reactant side (e.g. Gov in r10).
#       Only species with net_change < 0 appear in 'consumed'.

REACTION_DEFS: Dict[str, dict] = {

    # ── Land and Resource Dynamics (r1–r6) ───────────────────────────────────

    'r1': {   # SL => RL  — seasonal land recovery (Natural)
        'consumed'  : {'SL': 1},
        'enabling'  : [],
        'net_change': {'SL': -1, 'RL': +1},
        'seasonal'  : 'recovery',
    },
    'r2': {   # SR_SL => SR_RL  — seasonal recovery of resilient-land pop (Natural)
        'consumed'  : {'SR_SL': 1},
        'enabling'  : [],
        'net_change': {'SR_SL': -1, 'SR_RL': +1},
        'seasonal'  : 'recovery',
    },
    'r3': {   # WR_SL => WR_RL  — seasonal recovery of vulnerable-land pop (Natural)
        'consumed'  : {'WR_SL': 1},
        'enabling'  : [],
        'net_change': {'WR_SL': -1, 'WR_RL': +1},
        'seasonal'  : 'recovery',
    },
    'r4': {   # RL => SL  — seasonal land degradation (Natural, opposite phase)
        'consumed'  : {'RL': 1},
        'enabling'  : [],
        'net_change': {'RL': -1, 'SL': +1},
        'seasonal'  : 'degradation',
    },
    'r5': {   # AG_RL => AG_SL   — diminishing extraction returns (Natural)
        # K[AG_RL]: as more AG crowd territory per-unit extraction falls
        'consumed'  : {'AG_RL': 1},
        'enabling'  : [],
        'net_change': {'AG_RL': -1, 'AG_SL': +1},
        'seasonal'  : 'degradation',
    },
    'r6': {   # Gov + SR_SL => SR_RL  — governance-driven rehabilitation (Development)
        # K[Gov]: governance capacity has a ceiling
        'consumed'  : {'Gov': 1, 'SR_SL': 1},
        'enabling'  : [],
        'net_change': {'Gov': -1, 'SR_SL': -1, 'SR_RL': +1},
        'seasonal'  : None,
    },

    # ── Economic Production and Governance (r7–r10) ──────────────────────────

    'r7': {   # SR_RL => SR_SL + E_peace  — economic productivity SR (Natural)
        'consumed'  : {'SR_RL': 1},
        'enabling'  : [],
        'net_change': {'SR_RL': -1, 'SR_SL': +1, 'E_peace': +1},
        'seasonal'  : None,
    },
    'r8': {   # WR_RL => WR_SL + E_peace  — economic productivity WR (Natural, weaker)
        'consumed'  : {'WR_RL': 1},
        'enabling'  : [],
        'net_change': {'WR_RL': -1, 'WR_SL': +1, 'E_peace': +1},
        'seasonal'  : None,
    },
    'r9': {   # WR_RL => WR_SL  — lazy behavior / soil degradation (Natural)
        'consumed'  : {'WR_RL': 1},
        'enabling'  : [],
        'net_change': {'WR_RL': -1, 'WR_SL': +1},
        'seasonal'  : None,
    },
    'r10': {  # Gov + E_peace => 2Gov  — taxation sustains governance (Development)
        # Gov is net-produced (+1): it acts as an enabler (does not draw budget)
        # E_peace is consumed (-1)
        'consumed'  : {'E_peace': 1},
        'enabling'  : ['Gov'],
        'net_change': {'E_peace': -1, 'Gov': +1},
        'seasonal'  : None,
    },

    # ── Population Migration (r12–r16) ───────────────────────────────────────

    'r12': {  # SR_SL + RL => SR_RL + SL  — SR seeks better soil (Natural)
        # K[RL]: land availability saturates migration benefit
        'consumed'  : {'SR_SL': 1, 'RL': 1},
        'enabling'  : [],
        'net_change': {'SR_SL': -1, 'RL': -1, 'SR_RL': +1, 'SL': +1},
        'seasonal'  : None,
    },
    'r13': {  # WR_SL + RL => WR_RL + SL  — WR seeks better soil, less able (Natural)
        'consumed'  : {'WR_SL': 1, 'RL': 1},
        'enabling'  : [],
        'net_change': {'WR_SL': -1, 'RL': -1, 'WR_RL': +1, 'SL': +1},
        'seasonal'  : None,
    },
    'r14': {  # SR_RL + AG_SL + V => AG_RL + SR_SL  — violence displaces SR (Extraction)
        # K[V]: moderate violence already triggers full displacement
        'consumed'  : {'SR_RL': 1, 'AG_SL': 1, 'V': 1},
        'enabling'  : [],
        'net_change': {'SR_RL': -1, 'AG_SL': -1, 'V': -1,
                       'AG_RL': +1, 'SR_SL': +1},
        'seasonal'  : None,
    },
    'r15': {  # WR_RL + AG_SL + V => AG_RL + WR_SL  — violence displaces WR (Extraction)
        'consumed'  : {'WR_RL': 1, 'AG_SL': 1, 'V': 1},
        'enabling'  : [],
        'net_change': {'WR_RL': -1, 'AG_SL': -1, 'V': -1,
                       'AG_RL': +1, 'WR_SL': +1},
        'seasonal'  : None,
    },
    'r16': {  # AG_SL + RL => AG_RL + SL  — AG occupies resourceful land (Coerce)
        'consumed'  : {'AG_SL': 1, 'RL': 1},
        'enabling'  : [],
        'net_change': {'AG_SL': -1, 'RL': -1, 'AG_RL': +1, 'SL': +1},
        'seasonal'  : None,
    },

    # ── Resilience Transitions (r17–r20) ─────────────────────────────────────

    'r17': {  # E_peace + WR_RL => SR_RL  — investment upgrades WR on good land (Development)
        # K[WR_RL]: upgrade capacity saturates (bottleneck in absorptive capacity)
        'consumed'  : {'E_peace': 1, 'WR_RL': 1},
        'enabling'  : [],
        'net_change': {'E_peace': -1, 'WR_RL': -1, 'SR_RL': +1},
        'seasonal'  : None,
    },
    'r18': {  # E_peace + WR_SL => SR_SL  — investment upgrades WR on bad land (Development)
        # Less effective than r17 (lower k)
        'consumed'  : {'E_peace': 1, 'WR_SL': 1},
        'enabling'  : [],
        'net_change': {'E_peace': -1, 'WR_SL': -1, 'SR_SL': +1},
        'seasonal'  : None,
    },
    'r19': {  # SR_SL => WR_SL  — resilience degradation without support (Natural)
        'consumed'  : {'SR_SL': 1},
        'enabling'  : [],
        'net_change': {'SR_SL': -1, 'WR_SL': +1},
        'seasonal'  : None,
    },
    'r20': {  # E_peace + Gov + AG_SL => WR_RL  — DDR reintegration (Defense)
        # K[AG_SL]: DDR throughput bottleneck; intake capacity saturates
        'consumed'  : {'E_peace': 1, 'Gov': 1, 'AG_SL': 1},
        'enabling'  : [],
        'net_change': {'E_peace': -1, 'Gov': -1, 'AG_SL': -1, 'WR_RL': +1},
        'seasonal'  : None,
    },

    # ── Trust and Social Cohesion (r21–r25) ──────────────────────────────────

    'r21': {  # SR_RL + Gov => SR_RL + Gov + 2T  — gov presence generates trust (Development)
        # Both SR_RL and Gov are catalysts (net_change=0); purely enabling
        # K values dampen trust at high amounts to avoid runaway production
        'consumed'  : {},
        'enabling'  : ['SR_RL', 'Gov'],
        'net_change': {'T': +2},
        'seasonal'  : None,
    },
    'r22': {  # Gov + WR_RL + T => SR_RL  — gov+trust improves WR resilience (Development)
        # K[T]: trust needs to cross a threshold to drive this transition
        'consumed'  : {'Gov': 1, 'WR_RL': 1, 'T': 1},
        'enabling'  : [],
        'net_change': {'Gov': -1, 'WR_RL': -1, 'T': -1, 'SR_RL': +1},
        'seasonal'  : None,
    },
    'r23': {  # V + T => (nothing)  — violence and trust annihilate (Natural)
        'consumed'  : {'V': 1, 'T': 1},
        'enabling'  : [],
        'net_change': {'V': -1, 'T': -1},
        'seasonal'  : None,
    },
    'r24': {  # Gov => (nothing)  — governance requires constant renewal (Natural)
        'consumed'  : {'Gov': 1},
        'enabling'  : [],
        'net_change': {'Gov': -1},
        'seasonal'  : None,
    },
    'r25': {  # T => (nothing)  — trust requires constant renewal (Natural)
        'consumed'  : {'T': 1},
        'enabling'  : [],
        'net_change': {'T': -1},
        'seasonal'  : None,
    },

    # ── Violence Generation (r26–r30) ─────────────────────────────────────────

    'r26': {  # 2WR_SL => 2WR_SL + V  — grievance-driven violence (Natural, catalytic)
        # WR_SL is a catalyst (2 in, 2 out); listed twice for squared effect
        # K[WR_SL] saturates each factor to avoid explosive growth
        'consumed'  : {},
        'enabling'  : ['WR_SL', 'WR_SL'],
        'net_change': {'V': +1},
        'seasonal'  : None,
    },
    'r27': {  # AG_SL + SR_RL + E_conflict => AG_RL + WR_SL + V  — AG destroys resilience (Coerce)
        # K[E_conflict]: conflict economy funds this; saturates at campaign capacity
        'consumed'  : {'AG_SL': 1, 'SR_RL': 1, 'E_conflict': 1},
        'enabling'  : [],
        'net_change': {'AG_SL': -1, 'SR_RL': -1, 'E_conflict': -1,
                       'AG_RL': +1, 'WR_SL': +1, 'V': +1},
        'seasonal'  : None,
    },
    'r28': {  # AG_SL + Gov => WR_SL + V  — governance turns AG scouts into WR (Defense)
        # K[Gov]: governance capacity saturates counter-insurgency throughput
        'consumed'  : {'AG_SL': 1, 'Gov': 1},
        'enabling'  : [],
        'net_change': {'AG_SL': -1, 'Gov': -1, 'WR_SL': +1, 'V': +1},
        'seasonal'  : None,
    },
    'r29': {  # Gov + AG_RL => AG_RL + V  — AG fighters defeat governance (Coerce)
        # AG_RL is catalyst (unchanged fighters); K[AG_RL] saturates combat power
        'consumed'  : {'Gov': 1},
        'enabling'  : ['AG_RL'],
        'net_change': {'Gov': -1, 'V': +1},
        'seasonal'  : None,
    },
    'r30': {  # V => (nothing)  — violence dissipates (Natural)
        'consumed'  : {'V': 1},
        'enabling'  : [],
        'net_change': {'V': -1},
        'seasonal'  : None,
    },

    # ── Armed Group Dynamics (r31–r34) ────────────────────────────────────────

    'r31': {  # WR_SL + AG_RL + E_conflict => AG_RL + AG_SL  — AG recruits WR (Coerce)
        # AG_RL is catalyst (recruiter pool); K[AG_RL] caps recruitment capacity
        # E_conflict funds it; K[E_conflict] saturates funding effect
        'consumed'  : {'WR_SL': 1, 'E_conflict': 1},
        'enabling'  : ['AG_RL'],
        'net_change': {'WR_SL': -1, 'E_conflict': -1, 'AG_SL': +1},
        'seasonal'  : None,
    },
    'r32': {  # AG_RL + WR_RL + E_peace => AG_RL + WR_SL + E_conflict  — looting WR (Extract)
        # AG_RL catalyst (looting capacity per fighter); K[AG_RL] caps it
        # WR_RL and E_peace are consumed (land access + extracted wealth)
        'consumed'  : {'WR_RL': 1, 'E_peace': 1},
        'enabling'  : ['AG_RL'],
        'net_change': {'WR_RL': -1, 'E_peace': -1, 'WR_SL': +1, 'E_conflict': +1},
        'seasonal'  : None,
    },
    'r33': {  # AG_RL + SR_RL + E_peace => AG_RL + WR_SL + E_conflict  — looting SR (Extract)
        # Same structure as r32 but less effective (lower k; SR more resilient)
        'consumed'  : {'SR_RL': 1, 'E_peace': 1},
        'enabling'  : ['AG_RL'],
        'net_change': {'SR_RL': -1, 'E_peace': -1, 'WR_SL': +1, 'E_conflict': +1},
        'seasonal'  : None,
    },
    'r34': {  # 2AG_RL + Gov => 2AG_RL + E_conflict  — corruption (Extract)
        # AG_RL catalyst listed twice for squared collusion effect
        # Gov is consumed (corrupted resource); K[Gov] saturates extraction
        'consumed'  : {'Gov': 1},
        'enabling'  : ['AG_RL', 'AG_RL'],
        'net_change': {'Gov': -1, 'E_conflict': +1},
        'seasonal'  : None,
    },
}

# ─────────────────────────────────────────────────────────────────────────────
# DEFAULT KINETICS
# ─────────────────────────────────────────────────────────────────────────────
# Starting values designed for dt=0.05 years and regime-1 typical amounts.
# Verify no allocation-sum > 1 warnings appear; tune from here.
#
# Group timescale rationale (annual rate at typical amounts):
#   Fast  (>> 1/yr) : land exchange, economic production, decay
#   Medium (~1/yr)  : migration, investment upgrades, trust dynamics
#   Slow  (<< 1/yr) : DDR, AG recruitment, governance corruption
#
# K values set to ~50% of expected [X] at regime 3 (mid-conflict) so the
# saturating form kicks in near the transition zone.

DEFAULT_KINETICS: Dict[str, dict] = {

    # ── Land dynamics ────────────────────────────────────────────────────────
    # Seasonal amplitude is separately scaled by climate; k here is the peak rate.
    'r1' : {'k': 0.40},                              # SL -> RL, seasonal recovery
    'r2' : {'k': 0.30},                              # SR_SL -> SR_RL
    'r3' : {'k': 0.20},                              # WR_SL -> WR_RL (weaker)
    'r4' : {'k': 0.25},                              # RL -> SL, seasonal degradation
    'r5' : {'k': 0.15, 'K': {'AG_RL': 3.0}},        # AG extraction; saturates in territory
    'r6' : {'k': 0.04, 'K': {'Gov': 20.0}},          # governance rehab; gov capacity ceiling

    # ── Economic production & governance ────────────────────────────────────
    'r7' : {'k': 1},                              # SR_RL -> productivity (fast)
    'r8' : {'k': 0.8},                              # WR_RL -> productivity (slower)
    'r9' : {'k': 0.05},                              # WR_RL lazy degradation
    'r10': {'k': 0.10, 'K': {'Gov': 25.0}},          # taxation; gov enabling saturates

    # ── Migration ────────────────────────────────────────────────────────────
    'r12': {'k': 0.03, 'K': {'RL': 60.0}},           # SR land-seeking; land availability sat.
    'r13': {'k': 0.05, 'K': {'RL': 60.0}},           # WR land-seeking (less able)
    'r14': {'k': 0.05, 'K': {'V': 10.0}},             # violence displaces SR; V saturates
    'r15': {'k': 0.05, 'K': {'V': 10.0}},             # violence displaces WR (easier)
    'r16': {'k': 0.05, 'K': {'RL': 60.0}},           # AG occupies available land

    # ── Resilience transitions ────────────────────────────────────────────────
    'r17': {'k': 0.04, 'K': {'WR_RL': 15.0}},        # invest upgrades WR on good land
    'r18': {'k': 0.02, 'K': {'WR_SL': 15.0}},        # invest upgrades WR on bad land (harder)
    'r19': {'k': 0.05},                              # SR_SL -> WR_SL, resilience decay
    'r20': {'k': 0.01, 'K': {'AG_SL': 20.0}},        # DDR; intake capacity bottleneck (slow, high threshold)

    # ── Trust and social cohesion ────────────────────────────────────────────
    # r21 is catalytic: k must be calibrated so T_eq ≈ 10 in regime 1.
    # With [SR_RL]=40, [Gov]=43, K_SR=30, K_Gov=30:
    #   lambda_21 ≈ k21 * (40/70) * (43/73) ≈ 0.335 * k21
    #   T_eq = 2*lambda_21 / (k25 + k23*[V]) ≈ 2*0.335*k21 / 0.12 ≈ 5.6*k21
    #   Target T_eq=10 → k21 ≈ 1.8
    'r21': {'k': 1, 'K': {'SR_RL': 30.0, 'Gov': 30.0}},  # catalytic trust gen.
    'r22': {'k': 0.03, 'K': {'T': 5.0}},             # gov+trust improves resilience
    'r23': {'k': 0.10},                              # V + T annihilation
    'r24': {'k': 0.12},                              # Gov decay
    'r25': {'k': 0.10},                              # T decay

    # ── Violence generation ──────────────────────────────────────────────────
    # r26 is catalytic with squared WR_SL effect.
    # With [WR_SL]=11 (R1), K_WR=20: lambda_26 = k26*(11/31)^2 ≈ 0.126*k26
    # Target V equilibrium from r26 alone ≈ 0.25 → k30*[V]=k26*0.126
    #   k26=0.25/0.126/k30 — with k30=0.20 → k26 ≈ 10  (catalytic, so this is ok)
    'r26': {'k': 1.00, 'K': {'WR_SL': 20.0}},        # grievance violence (catalytic sq.)
    'r27': {'k': 0.05, 'K': {'E_conflict': 10.0}},    # AG destroys resilience
    'r28': {'k': 0.1, 'K': {'Gov': 50.0}},          # counter-insurgency
    'r29': {'k': 0.1, 'K': {'AG_RL': 50.0}},         # AG defeats governance
    'r30': {'k': 0.20},                              # V decay

    # ── Armed group dynamics ─────────────────────────────────────────────────
    'r31': {'k': 0.05, 'K': {'AG_RL': 2.0, 'E_conflict': 5.0}},  # recruitment
    'r32': {'k': 0.05, 'K': {'AG_RL': 29.0}},         # AG loots WR
    'r33': {'k': 0.05, 'K': {'AG_RL': 20.0}},         # AG loots SR (harder)
    'r34': {'k': 0.05, 'K': {'AG_RL': 30.0, 'Gov': 25.0}},        # corruption (sq. AG)
}

# ─────────────────────────────────────────────────────────────────────────────
# REGIME DEFINITIONS
# ─────────────────────────────────────────────────────────────────────────────

REGIME_NAMES: Dict[int, str] = {
    1: 'I: Functional Stability',
    2: 'II: Fragile Governance',
    3: 'III: Resource Competition',
    4: 'IV: Insurgency Periphery',
    5: 'V: Insurgency Stronghold',
    6: 'VI: State Fragmentation',
}

REGIME_COLORS: Dict[int, str] = {
    1: '#27ae60',
    2: '#f1c40f',
    3: '#e67e22',
    4: '#e74c3c',
    5: '#8e44ad',
    6: '#2c3e50',
}

# IC ranges [lo, hi] per species per regime.
# Population units: tens-of-thousands of persons.
# E_peace, E_conflict, T, V, Gov: dimensionless state indicators.
REGIME_IC_RANGES: Dict[int, Dict[str, Tuple[float, float]]] = {
    1: {  # Functional Stability — high peace economy, minimal conflict
        'SR_RL': (32.14, 50.07), 'SR_SL': (8.93,  17.86),
        'WR_RL': (14.29, 21.43), 'WR_SL': (7.14,  14.29),
        'AG_RL': (0.00,   0.36), 'AG_SL': (0.00,   0.89),
        'RL':    (100.0, 100.0), 'SL':    (100.0, 100.0),
        'E_peace': (45.0, 55.0), 'E_conflict': (0.0,   3.0),
        'T':    (8.0,   15.0),  'V':    (0.0,   0.5),   'Gov': (35.0, 50.0),
    },
    2: {  # Fragile Governance
        'SR_RL': (20.00, 30.00), 'SR_SL': (8.33,  15.00),
        'WR_RL': (16.67, 26.67), 'WR_SL': (10.00, 20.00),
        'AG_RL': (0.00,   1.11), 'AG_SL': (1.11,   3.33),
        'RL':    (90.0,  90.0),  'SL':    (110.0, 110.0),
        'E_peace': (35.0, 45.0), 'E_conflict': (5.0,  11.0),
        'T':    (4.0,   10.0),  'V':    (0.5,   1.0),   'Gov': (25.0, 35.0),
    },
    3: {  # Resource Competition
        'SR_RL': (10.00, 20.00), 'SR_SL': (5.71,  11.43),
        'WR_RL': (11.43, 20.00), 'WR_SL': (21.43, 35.71),
        'AG_RL': (0.00,   0.43), 'AG_SL': (0.00,   1.43),
        'RL':    (80.0,  80.0),  'SL':    (120.0, 120.0),
        'E_peace': (18.0, 26.0), 'E_conflict': (11.0, 17.0),
        'T':    (2.0,    8.0),  'V':    (1.0,   1.5),   'Gov': (18.0, 25.0),
    },
    4: {  # Insurgency Periphery
        'SR_RL': (7.14,  14.29), 'SR_SL': (4.29,  10.00),
        'WR_RL': (8.57,  14.29), 'WR_SL': (21.43, 35.71),
        'AG_RL': (0.71,   2.14), 'AG_SL': (2.14,   4.29),
        'RL':    (70.0,  70.0),  'SL':    (130.0, 130.0),
        'E_peace': (16.0, 24.0), 'E_conflict': (12.0, 18.0),
        'T':    (1.0,    5.0),  'V':    (1.5,   2.0),   'Gov': (12.0, 18.0),
    },
    5: {  # Insurgency Stronghold
        'SR_RL': (2.86,  10.00), 'SR_SL': (2.86,   7.14),
        'WR_RL': (4.29,  10.00), 'WR_SL': (28.57, 45.71),
        'AG_RL': (2.14,   4.29), 'AG_SL': (2.86,   5.71),
        'RL':    (60.0,  60.0),  'SL':    (140.0, 140.0),
        'E_peace': (7.0,  13.0), 'E_conflict': (17.0, 23.0),
        'T':    (0.0,    3.0),  'V':    (2.0,   3.0),   'Gov': (8.0,  12.0),
    },
    6: {  # State Fragmentation
        'SR_RL': (4.29,  11.43), 'SR_SL': (2.86,   7.14),
        'WR_RL': (4.29,  10.00), 'WR_SL': (28.57, 42.86),
        'AG_RL': (2.86,   5.71), 'AG_SL': (2.86,   6.43),
        'RL':    (50.0,  50.0),  'SL':    (150.0, 150.0),
        'E_peace': (1.0,   7.0), 'E_conflict': (30.0, 40.0),
        'T':    (0.0,    3.0),  'V':    (3.0,   6.0),   'Gov': (0.0,   8.0),
    },
}

# ─────────────────────────────────────────────────────────────────────────────
# UTILITY FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────

def regime_mean_ic(regime: int) -> Dict[str, float]:
    """Return mean of IC ranges for a given regime as a state dict."""
    ranges = REGIME_IC_RANGES[regime]
    return {sp: (lo + hi) / 2.0 for sp, (lo, hi) in ranges.items()}


def regime_sample_ic(regime: int, rng=None) -> Dict[str, float]:
    """Return a uniform random sample from IC ranges for a given regime."""
    import numpy as np
    rng = rng or np.random.default_rng()
    ranges = REGIME_IC_RANGES[regime]
    return {sp: float(rng.uniform(lo, hi)) for sp, (lo, hi) in ranges.items()}


# ── OLD classifier (raw-species, two-ratio heuristic) ────────────────────────
# Kept for reference and backward-compatibility.  Not used by default.
#
# _RATIO_IDEAL       = 0.03
# _RATIO_WORRYING    = 0.09
# _RATIO_PROBLEMATIC = 0.20
# _CHAOS_MILD        = 0.50
# _CHAOS_HEAVY       = 0.70
#
# def classify_regime_raw(state: dict) -> int:
#     """
#     Old classifier: averaged ratio of AG share and V-vs-T balance.
#     Non-monotone at R2→R3 (Resource Competition is civilian-driven, not AG-driven).
#     """
#     total_pop = sum(state.get(sp, 0.0)
#                     for sp in ('SR_RL', 'SR_SL', 'WR_RL', 'WR_SL', 'AG_RL', 'AG_SL'))
#     ag  = state.get('AG_RL', 0.0) + state.get('AG_SL', 0.0)
#     v   = state.get('V',     0.0)
#     t   = state.get('T',     0.0)
#     ratio_ag = ag / (total_pop + 1e-9)
#     ratio_vt = v  / (v + t + 1e-9)
#     avg = (ratio_ag + ratio_vt) / 2.0
#     if avg < _RATIO_IDEAL:       return 1
#     if avg < _RATIO_WORRYING:    return 2
#     if avg < _RATIO_PROBLEMATIC: return 3
#     mx = max(ratio_ag, ratio_vt)
#     if mx < _CHAOS_MILD:  return 4
#     if mx < _CHAOS_HEAVY: return 5
#     return 6
# ─────────────────────────────────────────────────────────────────────────────


# ── NEW classifier: Conflict Index (CI) — indicator-based ────────────────────
#
# Three monotone indicators are combined into a single Conflict Index (CI):
#
#   war_economy     = E_conflict / (E_peace + E_conflict)
#                     Structural economic lock-in.  Cleanly separates all 6 regimes.
#
#   poverty         = WR_SL / (SR_RL + WR_SL)
#                     Livelihood collapse.  Captures civilian deterioration in R3
#                     even when AG is still small (the raw classifier missed this).
#
#   social_cohesion = T / (T + V)
#                     Social fabric erosion.  Trust–violence balance.
#
#   CI = war_economy + poverty + (1 - social_cohesion)   ∈ [0, 3]
#
# CI values at t=0 computed from regime IC means:
#   R1 ≈ 0.26  |  R2 ≈ 0.64  |  R3 ≈ 1.25
#   R4 ≈ 1.53  |  R5 ≈ 2.15  |  R6 ≈ 2.47
#
# Thresholds placed at midpoints between successive regime CI values:
_CI_THRESHOLDS = (0.45, 0.94, 1.39, 1.84, 2.31)
#  CI < 0.45  → R1   (Functional Stability)
#  CI < 0.94  → R2   (Fragile Governance)
#  CI < 1.39  → R3   (Resource Competition)
#  CI < 1.84  → R4   (Insurgency Periphery)
#  CI < 2.31  → R5   (Insurgency Stronghold)
#  CI ≥ 2.31  → R6   (State Fragmentation)


def conflict_index(state: dict) -> float:
    """
    Compute the Conflict Index (CI) for a state dict.  CI ∈ [0, 3].

    CI = war_economy + poverty + (1 - social_cohesion)

    Higher CI = more severe conflict.  The three components are all
    structural ratios that are monotone across the six regime ICs.
    """
    eps = 1e-9
    Ep  = state.get('E_peace',    0.0)
    Ec  = state.get('E_conflict', 0.0)
    SR  = state.get('SR_RL',      0.0)
    WS  = state.get('WR_SL',      0.0)
    T   = state.get('T',          0.0)
    V   = state.get('V',          0.0)
    war_economy     = Ec  / (Ep + Ec + eps)
    poverty         = WS  / (SR + WS  + eps)
    social_cohesion = T   / (T  + V   + eps)
    return war_economy + poverty + (1.0 - social_cohesion)


def classify_regime(state: dict) -> int:
    """
    Classify a state dict into regime 1–6 using the Conflict Index (CI).

    CI = war_economy + poverty + (1 - social_cohesion)

    Thresholds: 0.45 / 0.94 / 1.39 / 1.84 / 2.31
    These are mid-points between regime CI values computed from IC means.
    """
    ci = conflict_index(state)
    for regime, threshold in enumerate(_CI_THRESHOLDS, start=1):
        if ci < threshold:
            return regime
    return 6


def make_simulator(kinetics: Optional[dict] = None,
                   dt: float = 0.05,
                   climate: float = 0.3,
                   rng_seed: Optional[int] = None) -> BudgetAllocationSimulator:
    """
    Return a configured BudgetAllocationSimulator for the conflict model.

    Parameters
    ----------
    kinetics  : dict, optional   Override DEFAULT_KINETICS (partial or full).
    dt        : float            Time step in years.
    climate   : float            [0=dry … 1=wet].
    rng_seed  : int, optional
    """
    kin = dict(DEFAULT_KINETICS)
    if kinetics:
        for rxn_id, params in kinetics.items():
            kin[rxn_id] = dict(kin.get(rxn_id, {}))
            kin[rxn_id].update(params)
    return BudgetAllocationSimulator(
        reaction_defs=REACTION_DEFS,
        kinetics=kin,
        dt=dt,
        climate=climate,
        rng_seed=rng_seed,
    )
