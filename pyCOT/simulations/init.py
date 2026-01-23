"""
pyCOT Simulations Module

Unified interface for reaction network simulations across different spatial structures.

Simulation Types:
-----------------
1. ODE (Well-Mixed): Homogeneous systems, no spatial structure
2. Spatial (PDE): Reaction-diffusion on continuous domains
3. Metapopulation: Discrete patches with connectivity

All simulation types share:
- Common parameter interface
- Kinetic rule dispatch from pyCOT.kinetics
- Consistent output format
- Stoichiometric matrix integration

Usage:
------
# Standard ODE simulation
from pyCOT.simulations import simulation
ts, flux = simulation(rn, rate='mak', t_span=(0, 100))

# Spatial simulation (future)
from pyCOT.simulations import simulate_spatial_dynamics
ts, spatial_fields = simulate_spatial_dynamics(rn, grid_shape=(50, 50), ...)

# Metapopulation (future)
from pyCOT.simulations import simulate_metapopulation_dynamics
ts, patch_data = simulate_metapopulation_dynamics(rn, num_patches=10, ...)
"""

# Primary ODE simulation (most common use case)
from .ode import simulation

# Import core utilities for advanced users
from .core import (
    generate_random_vector,
    validate_rate_list,
    build_reaction_dict,
    parse_parameters,
    generate_default_parameters,
    print_differential_equations,
    print_velocity_expressions,
    time_series_dataframe,
    flux_vector_dataframe
)

# Spatial and metapopulation simulations (to be fully implemented)
# For now, keep placeholders for backward compatibility
try:
    from .spatial import simulate_spatial_dynamics
except ImportError:
    def simulate_spatial_dynamics(*args, **kwargs):
        raise NotImplementedError(
            "Spatial (PDE) simulations are being refactored. "
            "Use legacy simulations.py temporarily or implement spatial.py"
        )

try:
    from .metapopulation import simulate_metapopulation_dynamics
except ImportError:
    def simulate_metapopulation_dynamics(*args, **kwargs):
        raise NotImplementedError(
            "Metapopulation simulations are being refactored. "
            "Use legacy simulations.py temporarily or implement metapopulation.py"
        )


__all__ = [
    # Primary interface
    'simulation',
    
    # Spatial variants (future)
    'simulate_spatial_dynamics',
    'simulate_metapopulation_dynamics',
    
    # Core utilities
    'generate_random_vector',
    'validate_rate_list',
    'build_reaction_dict',
    'parse_parameters',
    'generate_default_parameters',
    'print_differential_equations',
    'print_velocity_expressions',
    'time_series_dataframe',
    'flux_vector_dataframe'
]


# Module version
__version__ = '2.0.0'  # Major refactor with kinetics separation
