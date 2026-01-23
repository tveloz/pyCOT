"""
Simulation Core Infrastructure

Shared utilities for all simulation types (ODE, PDE, metapopulation):
- Parameter parsing and validation
- Kinetic rule dispatch
- Reaction network processing
- Output formatting
- Common helper functions

This module provides the foundation that all simulation variants build upon.
"""

import numpy as np
import pandas as pd
from typing import List, Dict, Union, Optional, Callable


def generate_random_vector(n: int, seed=None, min_value=0.1, max_value=1.0):
    """
    Generate random parameter vector.
    
    Parameters:
    -----------
    n : int
        Vector size
    seed : int, optional
        Random seed for reproducibility
    min_value : float
        Minimum value
    max_value : float
        Maximum value
        
    Returns:
    --------
    array
        Random values rounded to 2 decimal places
    """
    rng = np.random.default_rng(seed)
    values = rng.uniform(min_value, max_value, n)
    return np.round(values, 2)


def validate_rate_list(rate: Union[str, List[str]], n_reactions: int) -> List[str]:
    """
    Validate and normalize rate specification.
    
    Parameters:
    -----------
    rate : str or list
        Single kinetic law name or list of names (one per reaction)
    n_reactions : int
        Expected number of reactions
        
    Returns:
    --------
    list
        Validated list of kinetic law names
        
    Raises:
    -------
    ValueError
        If list length doesn't match n_reactions
    """
    if isinstance(rate, str):
        return [rate] * n_reactions
    elif len(rate) != n_reactions:
        raise ValueError(
            f"Rate list length ({len(rate)}) must equal "
            f"number of reactions ({n_reactions})"
        )
    return rate


def build_reaction_dict(rn):
    """
    Build dictionary mapping reactions to their reactants and products.
    
    Parameters:
    -----------
    rn : ReactionNetwork
        pyCOT reaction network object
        
    Returns:
    --------
    dict
        {reaction_name: (reactants_list, products_list)}
        where each list contains (species, coefficient) tuples
    """
    rn_dict = {}
    species = rn.stoichiometry_matrix().species
    reactions = rn.stoichiometry_matrix().reactions
    reactants_vectors = rn.reactants_matrix().T
    products_vectors = rn.products_matrix().T
    
    for i in range(len(reactions)):
        reactants = [
            (species[j], reactants_vectors[i][j])
            for j in range(len(species))
            if reactants_vectors[i][j] > 0
        ]
        products = [
            (species[j], products_vectors[i][j])
            for j in range(len(species))
            if products_vectors[i][j] > 0
        ]
        rn_dict[reactions[i]] = (reactants, products)
    
    return rn_dict


def parse_parameters(rn, rn_dict, rate, spec_vector):
    """
    Parse parameter specifications into reaction-indexed dictionary.
    
    Parameters:
    -----------
    rn : ReactionNetwork
        Network structure
    rn_dict : dict
        Reaction dictionary from build_reaction_dict
    rate : list
        Kinetic law names per reaction
    spec_vector : list
        Parameters for each reaction
        
    Returns:
    --------
    dict
        {reaction_name: parameter_list}
    """
    reactions = rn.stoichiometry_matrix().reactions
    parameters = {}
    for i, reaction in enumerate(reactions):
        parameters[reaction] = spec_vector[i]
    return parameters


def generate_default_parameters(rate_list: List[str], 
                               n_reactions: int,
                               additional_laws: Optional[Dict] = None) -> List[List[float]]:
    """
    Generate default parameters for kinetic laws.
    
    Parameters:
    -----------
    rate_list : list
        Kinetic law names
    n_reactions : int
        Number of reactions
    additional_laws : dict, optional
        Custom kinetic law registry
        
    Returns:
    --------
    list
        Parameter specifications for each reaction
    """
    spec_vector = []
    
    for kinetic in rate_list:
        if kinetic == 'mak':
            params = generate_random_vector(1, min_value=0.01, max_value=1.0)
        elif kinetic == 'mmk':
            Vmax = generate_random_vector(1, min_value=1, max_value=1.5)
            Km = generate_random_vector(1, min_value=5, max_value=10)
            params = np.round([float(Vmax), float(Km)], 2)
        elif kinetic == 'hill':
            Vmax = generate_random_vector(1, min_value=1, max_value=1.5)
            Kd = generate_random_vector(1, min_value=5, max_value=10)
            n = generate_random_vector(1, min_value=1, max_value=4)
            params = np.round([float(Vmax), float(Kd), float(n)], 2)
        elif kinetic == 'cosine':
            A = generate_random_vector(1, min_value=0.5, max_value=2.0)
            w = generate_random_vector(1, min_value=0.1, max_value=1.0)
            params = np.round([float(A), float(w)], 2)
        elif kinetic == 'saturated':
            Vmax = generate_random_vector(1, min_value=1, max_value=1.5)
            Km = generate_random_vector(1, min_value=5, max_value=10)
            params = np.round([float(Vmax), float(Km)], 2)
        elif kinetic == 'threshold_memory':
            threshold = generate_random_vector(1, min_value=0.5, max_value=2.0)
            base_rate = generate_random_vector(1, min_value=0.5, max_value=1.5)
            memory_strength = generate_random_vector(1, min_value=0.1, max_value=0.8)
            decay_rate = generate_random_vector(1, min_value=0.1, max_value=0.5)
            steepness = 5.0
            params = np.round([float(threshold), float(base_rate), 
                             float(memory_strength), float(decay_rate), steepness], 2)
        elif kinetic in (additional_laws or {}):
            # Default 3 parameters for custom laws
            params = np.round(np.random.uniform(0.1, 1.0, 3), 3)
        else:
            raise ValueError(f"Unknown kinetic law: {kinetic}")
        
        spec_vector.append(params.tolist())
    
    return spec_vector


def print_differential_equations(rn, rn_dict):
    """
    Print symbolic differential equations for the system.
    
    Shows dx[i]/dt = ... for each species.
    """
    print("\nDifferential equations of the system:")
    species = rn.stoichiometry_matrix().species
    reactions = rn.stoichiometry_matrix().reactions
    
    for sp in species:
        term_dict = {}
        for i, reaction in enumerate(reactions):
            rate_expr = f"v{i+1}"
            coef = 0
            
            # Count consumption
            for reactant, coeff in rn_dict[reaction][0]:
                if reactant == sp:
                    coef -= coeff
            
            # Count production
            for product, coeff in rn_dict[reaction][1]:
                if product == sp:
                    coef += coeff
            
            if coef != 0:
                term_dict[rate_expr] = coef
        
        # Format equation
        eq_terms = []
        for v, coef in term_dict.items():
            if coef == 1:
                eq_terms.append(f"+{v}")
            elif coef == -1:
                eq_terms.append(f"-{v}")
            elif coef > 1:
                eq_terms.append(f"+{coef}*{v}")
            else:
                eq_terms.append(f"{coef}*{v}")
        
        rhs = " ".join(eq_terms) if eq_terms else "0"
        print(f"d[{sp}]/dt = {rhs}")


def print_velocity_expressions(rn, rn_dict, rate, rate_laws, additional_laws=None):
    """
    Print symbolic rate expressions for each reaction.
    
    Shows v_i = f([species]) based on kinetic law.
    """
    reactions = rn.stoichiometry_matrix().reactions
    print("\nAnalytical velocity expressions v(x, param):")
    
    for i, reaction in enumerate(reactions):
        kinetic = rate[i]
        vlabel = f"v{i+1}"
        reactants = rn_dict[reaction][0]
        
        rate_func = rate_laws.get(kinetic)
        if rate_func and hasattr(rate_func, 'expression'):
            expr = rate_func.expression(reactants, reaction)
        else:
            expr = f"undefined kinetics '{kinetic}' or missing expression"
        
        print(f"{vlabel} = {expr}  ({kinetic})")


def time_series_dataframe(rn, result, time_points):
    """
    Create time series DataFrame from simulation results.
    
    Parameters:
    -----------
    rn : ReactionNetwork
        Network structure
    result : array
        Simulation results (n_steps × n_species)
    time_points : array
        Time values
        
    Returns:
    --------
    DataFrame
        Columns: Time, species_1, species_2, ...
    """
    species = [specie.name for specie in rn.species()]
    data = {'Time': time_points}
    
    for i, sp in enumerate(species):
        data[sp] = result[:, i]
    
    return pd.DataFrame(data)


def flux_vector_dataframe(rn, rn_dict, rate, result, time_points, 
                         spec_vector, rate_laws):
    """
    Compute flux vector (reaction rates) over time.
    
    Parameters:
    -----------
    rn : ReactionNetwork
    rn_dict : dict
        Reaction dictionary
    rate : list
        Kinetic law names
    result : array
        Concentration time series
    time_points : array
        Time values
    spec_vector : list
        Parameters
    rate_laws : dict
        Kinetic function registry
        
    Returns:
    --------
    DataFrame
        Columns: Time, reaction_1, reaction_2, ...
    """
    species = [specie.name for specie in rn.species()]
    reactions = [reaction.name() for reaction in rn.reactions()]
    species_idx = {s: i for i, s in enumerate(species)}
    parameters = parse_parameters(rn, rn_dict, rate, spec_vector)
    
    flux_data = {'Time': time_points}
    
    for i, reaction in enumerate(reactions):
        kinetic = rate[i]
        reactants = rn_dict[reaction][0]
        param_vector = parameters[reaction]
        rate_fn = rate_laws[kinetic]
        
        rates = []
        for t_idx, t in enumerate(time_points):
            concentrations = result[t_idx, :]
            
            # Handle time-dependent kinetics
            if kinetic in ['cosine']:
                v = rate_fn(reactants, concentrations, species_idx, param_vector, t=t)
            else:
                v = rate_fn(reactants, concentrations, species_idx, param_vector)
            
            rates.append(v)
        
        flux_data[reaction] = rates
    
    return pd.DataFrame(flux_data)


__all__ = [
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
