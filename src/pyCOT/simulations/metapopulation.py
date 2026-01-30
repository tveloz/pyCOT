"""
Metapopulation Simulation Module

Implements discrete patch dynamics with local reactions and inter-patch dispersal.

Features:
- Multiple discrete patches/populations
- Local reaction dynamics in each patch
- Dispersal between patches via connectivity matrix
- Flexible connectivity (migration networks, metapopulations)

Mathematical form:
dx_i^(p)/dt = f_i(x^(p)) + Σ_q [D_i * C_qp * x_i^(q) - D_i * C_pq * x_i^(p)]

where:
- p, q are patch indices
- f_i(x^(p)) are local reactions in patch p
- D_i is dispersal rate for species i
- C_pq is connectivity from patch p to q
"""

import numpy as np
from scipy.integrate import solve_ivp
from typing import List, Dict, Optional, Tuple
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from kinetics import KINETIC_REGISTRY

from .core import (
    validate_rate_list,
    generate_default_parameters,
    generate_random_vector
)


def get_reaction_components(reaction, species):
    """
    Extract reactants, products, and stoichiometry from pyCOT Reaction object.
    
    Parameters:
    -----------
    reaction : Reaction
        pyCOT reaction object
    species : list
        List of species names
        
    Returns:
    --------
    reactants : list of tuple
        [(species_idx, coefficient), ...]
    products : list of tuple
        [(species_idx, coefficient), ...]
    stoichiometry : list
        Net stoichiometric coefficients for each species
    """
    try:
        reactants = [(species.index(edge.species_name), edge.coefficient) 
                     for edge in reaction.support_edges()]
        products = [(species.index(edge.species_name), edge.coefficient) 
                    for edge in reaction.products_edges()]
        
        stoichiometry = [0] * len(species)
        for sp_idx, coeff in reactants:
            stoichiometry[sp_idx] -= coeff
        for sp_idx, coeff in products:
            stoichiometry[sp_idx] += coeff
        
        return reactants, products, stoichiometry
    except Exception as e:
        print(f"Error extracting reaction components for {reaction.name()}: {e}")
        raise


def simulate_metapopulation_dynamics(rn,
                                     rate: str = 'mak',
                                     num_patches: Optional[int] = None,
                                     D_dict: Optional[Dict] = None,
                                     x0_dict: Optional[Dict] = None,
                                     spec_vector: Optional[List] = None,
                                     t_span: Tuple[float, float] = (0, 20),
                                     n_steps: int = 500,
                                     connectivity_matrix: Optional[np.ndarray] = None,
                                     additional_laws: Optional[Dict] = None,
                                     method: str = 'LSODA',
                                     rtol: float = 1e-8,
                                     atol: float = 1e-10):
    """
    Simulate metapopulation dynamics with local reactions and dispersal.
    
    Parameters:
    -----------
    rn : ReactionNetwork
        pyCOT reaction network object
    rate : str or list
        Kinetic law name(s)
    num_patches : int, optional
        Number of discrete patches/populations. Default: 3
    D_dict : dict, optional
        Dispersal rates {species_name: D_value}
        Default: random uniform [0.01, 0.2]
    x0_dict : dict, optional
        Initial patch concentrations {species_name: 1D_array(num_patches)}
        Default: random uniform [0, 2.0]
    spec_vector : list, optional
        Kinetic parameters per reaction
    t_span : tuple
        Time interval (t_start, t_end)
    n_steps : int
        Number of time points
    connectivity_matrix : array, optional
        Connectivity matrix (num_patches × num_patches)
        C[i,j] = probability/rate of dispersal from patch i to j
        Rows must sum to 1 (conservation of mass)
        Default: random normalized matrix
    additional_laws : dict, optional
        Custom kinetic functions
    method : str
        Integration method ('LSODA', 'RK45', etc.)
    rtol, atol : float
        Relative and absolute tolerances
        
    Returns:
    --------
    t : array
        Time points (n_steps,)
    X_out : dict
        Patch concentration time series {species: array(n_steps, num_patches)}
    flux_out : dict
        Patch flux time series {reaction: array(n_steps, num_patches)}
        
    Notes:
    ------
    - Each patch evolves with local reaction dynamics
    - Dispersal couples patches: species move according to connectivity
    - Connectivity matrix normalized so rows sum to 1
    - Dispersal flux: inflow from patch j to i minus outflow from i to j
    
    Example:
    --------
    # 3 patches with asymmetric connectivity
    connectivity = np.array([
        [0.7, 0.2, 0.1],  # Patch 0: 70% stay, 20% to 1, 10% to 2
        [0.1, 0.8, 0.1],  # Patch 1: 10% to 0, 80% stay, 10% to 2
        [0.2, 0.2, 0.6]   # Patch 2: 20% to 0, 20% to 1, 60% stay
    ])
    
    ts, X, flux = simulate_metapopulation_dynamics(
        rn, num_patches=3, connectivity_matrix=connectivity
    )
    """
    np.random.seed(seed=42)
    species = [specie.name for specie in rn.species()]
    reactions = [reaction.name() for reaction in rn.reactions()]
    
    rate = validate_rate_list(rate, len(reactions))
    
    # Default number of patches
    if num_patches is None:
        num_patches = 3
    
    # Default dispersal rates
    if D_dict is None:
        D_dict = {sp: np.round(np.random.uniform(0.01, 0.2), 3) for sp in species}
    
    # Default initial conditions
    if x0_dict is None:
        x0_dict = {sp: np.round(np.random.uniform(0, 2.0, size=num_patches), 2) 
                   for sp in species}
    else:
        # Validate x0_dict
        for sp in species:
            if sp not in x0_dict:
                raise ValueError(f"Missing initial condition for species '{sp}'")
            if len(x0_dict[sp]) != num_patches:
                raise ValueError(
                    f"Initial condition for '{sp}' must have length {num_patches}"
                )
    
    # Default connectivity matrix (random, row-normalized)
    if connectivity_matrix is None:
        connectivity_matrix = np.random.uniform(0, 1, size=(num_patches, num_patches))
        connectivity_matrix = np.round(connectivity_matrix, 2)
        
        # Normalize rows to sum to 1
        row_sums = np.sum(connectivity_matrix, axis=1, keepdims=True)
        connectivity_matrix = connectivity_matrix / row_sums
        connectivity_matrix = np.round(connectivity_matrix, 2)
        
        # Fix rounding errors
        for i in range(num_patches):
            current_sum = np.sum(connectivity_matrix[i])
            if current_sum != 1.0:
                diff = 1.0 - current_sum
                max_idx = np.argmax(connectivity_matrix[i])
                connectivity_matrix[i, max_idx] += diff
                connectivity_matrix[i, max_idx] = np.round(
                    connectivity_matrix[i, max_idx], 2
                )
    else:
        # Validate connectivity matrix
        if connectivity_matrix.shape != (num_patches, num_patches):
            raise ValueError(
                f"Connectivity matrix must have shape ({num_patches}, {num_patches})"
            )
        if not np.allclose(np.sum(connectivity_matrix, axis=1), 1.0, atol=1e-6):
            print("Warning: Normalizing connectivity matrix so rows sum to 1")
            connectivity_matrix = connectivity_matrix / np.sum(
                connectivity_matrix, axis=1, keepdims=True
            )
            connectivity_matrix = np.round(connectivity_matrix, 2)
    
    print("\nConnectivity matrix:\n", connectivity_matrix)
    print(f"Row sums: {np.sum(connectivity_matrix, axis=1)}")
    
    # Generate kinetic parameters
    if spec_vector is None:
        spec_vector = generate_default_parameters(rate, len(reactions), additional_laws)
    
    # State flattening/reshaping utilities
    def flatten_state(x_dict):
        return np.concatenate([x_dict[sp] for sp in species])
    
    def reshape_state(x):
        return {sp: x[i*num_patches:(i+1)*num_patches] 
                for i, sp in enumerate(species)}
    
    x0 = flatten_state(x0_dict)
    
    # Build kinetic law registry
    rate_laws = dict(KINETIC_REGISTRY)
    if additional_laws:
        rate_laws.update(additional_laws)
    
    # Local reaction dynamics in each patch
    def reaction_dynamics(Xdict, rate, spec_vector):
        dxdt_dict = {sp: np.zeros(num_patches) for sp in species}
        flux_dict = {r: np.zeros(num_patches) for r in reactions}
        
        # Iterate over patches
        for patch_idx in range(num_patches):
            local_x = [Xdict[sp][patch_idx] for sp in species]
            
            # Compute reaction rates
            for r_idx, reaction in enumerate(rn.reactions()):
                reactants, products, stoichiometry = get_reaction_components(
                    reaction, species
                )
                
                # Compute rate based on kinetic law
                kinetic = rate[r_idx]
                params = spec_vector[r_idx]
                
                if kinetic == 'mak':
                    k = params[0]
                    v_r = k
                    if reactants:
                        for sp_idx, stoich in reactants:
                            v_r *= local_x[sp_idx] ** abs(stoich)
                elif kinetic == 'mmk':
                    Vmax, Km = params
                    sp_idx = reactants[0][0] if reactants else 0
                    v_r = (Vmax * local_x[sp_idx] / (Km + local_x[sp_idx]) 
                           if reactants else 0)
                elif kinetic == 'hill':
                    Vmax, Kd, n = params
                    sp_idx = reactants[0][0] if reactants else 0
                    v_r = (Vmax * (local_x[sp_idx] ** n) / 
                           (Kd ** n + local_x[sp_idx] ** n) if reactants else 0)
                elif kinetic in rate_laws:
                    # Use registry for advanced kinetics
                    species_idx = {sp: idx for idx, sp in enumerate(species)}
                    v_r = rate_laws[kinetic](
                        [(species[sp_idx], coef) for sp_idx, coef in reactants],
                        local_x,
                        species_idx,
                        params
                    )
                else:
                    v_r = 0
                
                flux_dict[reaction.name()][patch_idx] = v_r
                
                # Apply stoichiometry
                for sp_idx, stoich in enumerate(stoichiometry):
                    if stoich != 0:
                        dxdt_dict[species[sp_idx]][patch_idx] += stoich * v_r
        
        return dxdt_dict, flux_dict
    
    # Dispersal term (connectivity-based)
    def dispersal_term(Xdict):
        dxdt_dict = {sp: np.zeros(num_patches) for sp in species}
        
        for sp in species:
            D = D_dict.get(sp, 0.0)  # Dispersal rate
            X = Xdict[sp]  # Concentration vector
            
            for i in range(num_patches):
                dxdt = 0
                for j in range(num_patches):
                    # Inflow from patch j to patch i
                    dxdt += D * connectivity_matrix[j, i] * X[j]
                    # Outflow from patch i to patch j
                    dxdt -= D * connectivity_matrix[i, j] * X[i]
                dxdt_dict[sp][i] = dxdt
        
        return dxdt_dict
    
    # Combined ODE system
    def combined_ode(t, x):
        x = np.maximum(x, 0)  # Non-negativity
        Xdict = reshape_state(x)
        dxdt_reac, _ = reaction_dynamics(Xdict, rate, spec_vector)
        dxdt_disperse = dispersal_term(Xdict)
        dxdt_total = {sp: dxdt_reac[sp] + dxdt_disperse[sp] for sp in species}
        return flatten_state(dxdt_total)
    
    # Integrate
    t_eval = np.linspace(t_span[0], t_span[1], n_steps)
    sol = solve_ivp(combined_ode, t_span, x0, t_eval=t_eval,
                    method=method, rtol=rtol, atol=atol)
    
    # Format output: patch time series
    X_out = {sp: np.zeros((n_steps, num_patches)) for sp in species}
    for i, xt in enumerate(sol.y.T):
        xt_dict = reshape_state(xt)
        for sp in species:
            X_out[sp][i] = xt_dict[sp]
    
    # Compute flux time series
    flux_out = {r: np.zeros((n_steps, num_patches)) for r in reactions}
    for time_idx, t_val in enumerate(sol.t):
        x_val = sol.y[:, time_idx]
        Xdict = reshape_state(x_val)
        _, flux_dict = reaction_dynamics(Xdict, rate, spec_vector)
        
        for r in reactions:
            flux_out[r][time_idx] = flux_dict[r]
    
    return sol.t, X_out, flux_out


__all__ = ['simulate_metapopulation_dynamics']
