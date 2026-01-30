"""
pyCOT Kinetics Module

Comprehensive collection of kinetic rate laws for reaction network simulations.

Organization:
-------------
- deterministic_basic: MAK, MMK, Hill (canonical forms)
- deterministic_advanced: Threshold-memory, saturation, oscillatory, inhibition
- stochastic_events: Poisson events, Gillespie SSA, tau-leaping, Langevin

Usage:
------
# Import specific kinetics
from pyCOT.kinetics.deterministic_basic import rate_mak, rate_mmk
from pyCOT.kinetics.deterministic_advanced import rate_threshold_memory
from pyCOT.kinetics.stochastic_events import PoissonEventGenerator

# Or import full registry
from pyCOT.kinetics import KINETIC_REGISTRY

# Use in simulation
from pyCOT.simulations import simulation
time_series, flux = simulation(rn, rate='threshold_memory', ...)
"""

from .deterministic_basic import (
    rate_mak,
    rate_mmk,
    rate_hill,
    BASIC_KINETICS
)

from .deterministic_advanced import (
    rate_cosine,
    rate_saturated,
    rate_threshold_memory,
    rate_inhibition,
    ADVANCED_KINETICS
)

from .stochastic_events import (
    PoissonEventGenerator,
    concentration_shock,
    rate_perturbation_event,
    GillespieSSA,
    TauLeaping,
    LangevinDynamics,
    STOCHASTIC_METHODS
)


# Unified kinetic registry for simulation dispatch
KINETIC_REGISTRY = {
    **BASIC_KINETICS,
    **ADVANCED_KINETICS,
    # Stochastic methods are class-based, handled separately
}


def get_kinetic_function(kinetic_name: str):
    """
    Retrieve kinetic function by name.
    
    Parameters:
    -----------
    kinetic_name : str
        Name of kinetic law ('mak', 'mmk', 'threshold_memory', etc.)
        
    Returns:
    --------
    callable
        Kinetic rate function
        
    Raises:
    -------
    KeyError
        If kinetic_name not in registry
    """
    if kinetic_name not in KINETIC_REGISTRY:
        available = ', '.join(KINETIC_REGISTRY.keys())
        raise KeyError(
            f"Kinetic law '{kinetic_name}' not found. "
            f"Available: {available}"
        )
    return KINETIC_REGISTRY[kinetic_name]


def list_kinetics(category=None):
    """
    List available kinetic laws, optionally filtered by category.
    
    Parameters:
    -----------
    category : str, optional
        'basic', 'advanced', or 'stochastic'
        
    Returns:
    --------
    list
        Names of available kinetic laws
    """
    if category == 'basic':
        return list(BASIC_KINETICS.keys())
    elif category == 'advanced':
        return list(ADVANCED_KINETICS.keys())
    elif category == 'stochastic':
        return list(STOCHASTIC_METHODS.keys())
    elif category is None:
        return list(KINETIC_REGISTRY.keys())
    else:
        raise ValueError(f"Unknown category '{category}'")


__all__ = [
    # Basic kinetics
    'rate_mak',
    'rate_mmk',
    'rate_hill',
    
    # Advanced kinetics
    'rate_cosine',
    'rate_saturated',
    'rate_threshold_memory',
    'rate_inhibition',
    
    # Stochastic events
    'PoissonEventGenerator',
    'concentration_shock',
    'rate_perturbation_event',
    'GillespieSSA',
    'TauLeaping',
    'LangevinDynamics',
    
    # Registries and utilities
    'KINETIC_REGISTRY',
    'BASIC_KINETICS',
    'ADVANCED_KINETICS',
    'STOCHASTIC_METHODS',
    'get_kinetic_function',
    'list_kinetics'
]
