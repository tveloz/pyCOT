"""
Basic Deterministic Kinetic Laws

This module contains fundamental, well-established kinetic rate laws:
- Mass Action Kinetics (MAK)
- Michaelis-Menten Kinetics (MMK)
- Hill Equation

These are canonical forms used for teaching and as the foundation for more
complex kinetic behaviors.
"""

import numpy as np


def rate_mak(reactants, concentrations, species_idx, spec_vector):
    """
    Mass Action Kinetics (MAK)
    
    Rate is proportional to the product of reactant concentrations raised
    to their stoichiometric coefficients.
    
    Parameters:
    -----------
    reactants : list of tuple
        List of (species_name, stoichiometric_coefficient) tuples
    concentrations : array-like
        Current concentrations of all species
    species_idx : dict
        Mapping from species names to concentration array indices
    spec_vector : list
        Parameters [k] where k is the rate constant
        
    Returns:
    --------
    float
        Reaction rate
        
    Example:
    --------
    For reaction A + 2B -> C with k=0.5:
    rate = 0.5 * [A] * [B]^2
    """
    k = spec_vector[0]
    rate = k
    for sp, coef in reactants:
        rate *= concentrations[species_idx[sp]] ** coef
    return rate


rate_mak.expression = lambda reactants, reaction: (
    f"k_{reaction} * " + " * ".join(
        f"[{reactant}]^{coef}" if coef != 1 else f"[{reactant}]"
        for reactant, coef in reactants
    )
)


def rate_mmk(reactants, concentrations, species_idx, spec_vector):
    """
    Michaelis-Menten Kinetics (MMK)
    
    Describes enzymatic reactions with substrate saturation.
    Rate approaches Vmax as substrate concentration increases.
    
    Parameters:
    -----------
    reactants : list of tuple
        List of (species_name, stoichiometric_coefficient) tuples
        First reactant is treated as the substrate
    concentrations : array-like
        Current concentrations of all species
    species_idx : dict
        Mapping from species names to concentration array indices
    spec_vector : list
        Parameters [Vmax, Km] where:
        - Vmax: maximum reaction rate
        - Km: Michaelis constant (substrate conc at half-maximal rate)
        
    Returns:
    --------
    float
        Reaction rate
        
    Example:
    --------
    For enzyme-substrate reaction:
    rate = Vmax * [S] / (Km + [S])
    """
    if not reactants:
        return 0.0
        
    Vmax, Km = spec_vector
    substrate = reactants[0][0]
    S = concentrations[species_idx[substrate]]
    
    return Vmax * S / (Km + S)


rate_mmk.expression = lambda reactants, reaction: (
    f"(Vmax_{reactants[0][0]} * [{reactants[0][0]}]) / "
    f"(Km_{reactants[0][0]} + [{reactants[0][0]}])"
    if reactants else "0  (mmk without defined substrate)"
)


def rate_hill(substrate, concentrations, species_idx, spec_vector):
    """
    Hill Equation
    
    Generalization of Michaelis-Menten for cooperative binding.
    Exhibits sigmoidal kinetics with Hill coefficient n controlling steepness.
    
    Parameters:
    -----------
    substrate : list of tuple or list
        Substrate species (can be single or multiple)
    concentrations : array-like
        Current concentrations of all species
    species_idx : dict
        Mapping from species names to concentration array indices
    spec_vector : list
        Parameters [Vmax, K, n] where:
        - Vmax: maximum reaction rate
        - K: dissociation constant
        - n: Hill coefficient (cooperativity)
        
    Returns:
    --------
    float
        Reaction rate
        
    Example:
    --------
    For cooperative binding:
    rate = Vmax * [S]^n / (K^n + [S]^n)
    
    n > 1: positive cooperativity (sigmoidal)
    n = 1: reduces to Michaelis-Menten
    n < 1: negative cooperativity
    """
    if not substrate:
        return 0.0
        
    Vmax, K, n = spec_vector
    
    # Handle different substrate input formats
    if isinstance(substrate, list) and isinstance(substrate[0], tuple):
        substrate_names = [s[0] for s in substrate]
    elif isinstance(substrate, list):
        substrate_names = substrate
    else:
        substrate_names = [substrate]
    
    # Sum concentrations of all substrate species
    S = sum(concentrations[species_idx[s]] for s in substrate_names)
    
    return Vmax * (S**n) / (K**n + S**n)


rate_hill.expression = lambda reactants, reaction: (
    f"(Vmax_{reactants[0][0]} * [{reactants[0][0]}]^n) / "
    f"(Km_{reactants[0][0]}^n + [{reactants[0][0]}]^n)"
    if reactants else "0  (hill without defined substrate)"
)


# Registry of basic kinetic laws
BASIC_KINETICS = {
    'mak': rate_mak,
    'mmk': rate_mmk,
    'hill': rate_hill
}


__all__ = [
    'rate_mak',
    'rate_mmk', 
    'rate_hill',
    'BASIC_KINETICS'
]
