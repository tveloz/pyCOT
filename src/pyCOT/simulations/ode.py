"""
ODE Simulation Module

Implements well-mixed (homogeneous) reaction network simulations using
ordinary differential equations (ODEs).

Assumptions:
- Spatial homogeneity (no diffusion/transport)
- Continuous concentrations (deterministic dynamics)
- Smooth time evolution

Integration methods:
- solve_ivp (default: LSODA, adaptive)
- odeint (fallback)

Supports all kinetic laws from pyCOT.kinetics module.
"""

import numpy as np
from scipy.integrate import solve_ivp, odeint
from typing import List, Dict, Optional, Union, Tuple
import pandas as pd

# Import kinetics registry
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pyCOT.kinetics import KINETIC_REGISTRY

# Import core utilities
from pyCOT.simulations.core import (
    validate_rate_list,
    build_reaction_dict,
    parse_parameters,
    generate_default_parameters,
    print_differential_equations,
    print_velocity_expressions,
    time_series_dataframe,
    flux_vector_dataframe,
    generate_random_vector
)


def simulation(rn,
               rate: Union[str, List[str]] = 'mak',
               spec_vector: Optional[List] = None,
               x0: Optional[List] = None,
               t_span: Tuple[float, float] = (0, 100),
               n_steps: int = 1000,
               additional_laws: Optional[Dict] = None,
               method: str = 'LSODA',
               rtol: float = 1e-8,
               atol: float = 1e-10,
               verbose: bool = True) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Simulate well-mixed reaction network using ODEs.
    
    Parameters:
    -----------
    rn : ReactionNetwork
        pyCOT reaction network object
    rate : str or list
        Kinetic law name(s) - single string or list per reaction
        Available: 'mak', 'mmk', 'hill', 'cosine', 'saturated', 'threshold_memory'
    spec_vector : list, optional
        Parameters for each reaction. If None, generates random defaults.
        Format: [[params_r1], [params_r2], ...]
    x0 : list, optional
        Initial concentrations. If None, generates random values [50-100].
    t_span : tuple
        Time interval (t_start, t_end)
    n_steps : int
        Number of time points for output
    additional_laws : dict, optional
        Custom kinetic functions {name: function}
    method : str
        Integration method ('LSODA', 'RK45', 'BDF', etc.)
    rtol : float
        Relative tolerance for integrator
    atol : float
        Absolute tolerance for integrator
    verbose : bool
        Print equations and parameters
        
    Returns:
    --------
    time_series_df : DataFrame
        Columns: Time, species_1, species_2, ...
    flux_vector_df : DataFrame
        Columns: Time, reaction_1, reaction_2, ...
        
    Examples:
    ---------
    # Basic simulation with mass action kinetics
    >>> ts, flux = simulation(rn, rate='mak', t_span=(0, 100))
    
    # Mixed kinetics
    >>> ts, flux = simulation(rn, rate=['mak', 'mak', 'mmk', 'hill'])
    
    # Threshold-memory kinetics
    >>> params = [[0.5, 1.0, 0.6, 0.1, 5.0]]  # threshold, base_rate, mem_str, decay, steep
    >>> ts, flux = simulation(rn, rate='threshold_memory', spec_vector=params)
    """
    # Extract species and reactions
    species = [specie.name for specie in rn.species()]
    reactions = [reaction.name() for reaction in rn.reactions()]
    rn_dict = build_reaction_dict(rn)
    
    # Validate and normalize rate specification
    rate = validate_rate_list(rate, len(reactions))
    
    # Generate initial conditions if not provided
    if x0 is None:
        x0 = generate_random_vector(len(species), min_value=50, max_value=100.0).astype(int).tolist()
        if verbose:
            print(f"\nx0 = {x0}")
    
    # Generate parameters if not provided
    if spec_vector is None:
        spec_vector = generate_default_parameters(rate, len(reactions), additional_laws)
        parameters = parse_parameters(rn, rn_dict, rate, spec_vector)
        if verbose:
            print("spec_vector =", spec_vector)
    else:
        parameters = parse_parameters(rn, rn_dict, rate, spec_vector)
        if verbose:
            print("spec_vector =", spec_vector)
    
    # Build kinetic law registry
    rate_laws = dict(KINETIC_REGISTRY)
    if additional_laws:
        for name in rate:
            if name not in rate_laws:
                if name in additional_laws:
                    rate_laws[name] = additional_laws[name]
                else:
                    raise NotImplementedError(
                        f"Kinetic law '{name}' not defined or registered"
                    )
    
    # Create species index mapping
    species_idx = {s: i for i, s in enumerate(species)}
    
    # Define ODE system with non-negativity constraint
    def rates_fn_constrained(t, x):
        # Enforce non-negative concentrations
        x_constrained = np.maximum(x, 0)
        dxdt = np.zeros_like(x_constrained)
        concentrations = x_constrained
        rates = np.zeros(len(reactions))
        
        # Compute reaction rates
        for i, reaction in enumerate(reactions):
            kinetic = rate[i]
            reactants = rn_dict[reaction][0]
            
            # Skip MMK with no reactants
            if kinetic == 'mmk' and (not reactants or len(reactants) == 0):
                continue
            
            param_vector = parameters[reaction]
            rate_fn = rate_laws[kinetic]
            
            # Handle time-dependent kinetics
            if kinetic in ['cosine', 'threshold_memory']:
                rates[i] = rate_fn(reactants, concentrations, species_idx, param_vector, t=t)
            else:
                rates[i] = rate_fn(reactants, concentrations, species_idx, param_vector)
        
        # Apply stoichiometry
        for i, reaction in enumerate(reactions):
            # Consumption (reactants)
            for sp, coef in rn_dict[reaction][0]:
                dxdt[species_idx[sp]] -= coef * rates[i]
            # Production (products)
            for sp, coef in rn_dict[reaction][1]:
                dxdt[species_idx[sp]] += coef * rates[i]
        
        return dxdt
    
    # Print system information if verbose
    if verbose:
        print_differential_equations(rn, rn_dict)
        print_velocity_expressions(rn, rn_dict, rate, rate_laws, additional_laws)
    
    # Time points for output
    time_points = np.linspace(t_span[0], t_span[1], n_steps)
    
    # Integrate ODEs using solve_ivp
    try:
        sol = solve_ivp(
            rates_fn_constrained,
            t_span,
            x0,
            t_eval=time_points,
            method=method,
            rtol=rtol,
            atol=atol
        )
        
        if not sol.success:
            print(f"Warning: Integration failed: {sol.message}")
            print("Trying with LSODA method...")
            sol = solve_ivp(
                rates_fn_constrained,
                t_span,
                x0,
                t_eval=time_points,
                method='LSODA',
                rtol=1e-6,
                atol=1e-8
            )
        
        # Extract and clean results
        result = sol.y.T
        result = np.maximum(result, 0)  # Remove small negative artifacts
        time_points = sol.t  # Use actual time points from solver
        
        # Report if significant negative values were corrected
        negative_mask = sol.y.T < -1e-10
        if np.any(negative_mask):
            neg_count = np.sum(negative_mask)
            print(f"Warning: {neg_count} negative values detected and corrected to 0")
    
    except Exception as e:
        print(f"Error during integration: {e}")
        print("Using fallback method (odeint)...")
        
        # Fallback to odeint
        def rates_fn_odeint(x, t):
            return rates_fn_constrained(t, x)
        
        result = odeint(rates_fn_odeint, x0, time_points)
        result = np.maximum(result, 0)
    
    # Create output DataFrames
    time_series_df = time_series_dataframe(rn, result, time_points)
    flux_vector_df = flux_vector_dataframe(
        rn, rn_dict, rate, result, time_points, spec_vector, rate_laws
    )
    
    return time_series_df, flux_vector_df


__all__ = ['simulation']
