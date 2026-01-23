"""
Spatial (PDE) Simulation Module

Implements reaction-diffusion dynamics on 2D grids using finite difference methods.

Features:
- Reaction networks on spatial grids
- Diffusion via Laplacian operators
- Periodic boundary conditions
- All kinetic laws supported

Mathematical form:
∂x_i/∂t = f_i(x) + D_i ∇²x_i

where f_i(x) are reaction terms and D_i are diffusion coefficients.
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
        # Reactants from support_edges()
        reactants = [(species.index(edge.species_name), edge.coefficient) 
                     for edge in reaction.support_edges()]
        # Products from products_edges()
        products = [(species.index(edge.species_name), edge.coefficient) 
                    for edge in reaction.products_edges()]
        
        # Build stoichiometry manually
        stoichiometry = [0] * len(species)
        for sp_idx, coeff in reactants:
            stoichiometry[sp_idx] -= coeff  # Reactants negative
        for sp_idx, coeff in products:
            stoichiometry[sp_idx] += coeff  # Products positive
        
        return reactants, products, stoichiometry
    except Exception as e:
        print(f"Error extracting reaction components for {reaction.name()}: {e}")
        raise


def simulate_spatial_dynamics(rn,
                              rate: str = 'mak',
                              grid_shape: Optional[Tuple[int, int]] = None,
                              D_dict: Optional[Dict] = None,
                              x0_dict: Optional[Dict] = None,
                              spec_vector: Optional[List] = None,
                              t_span: Tuple[float, float] = (0, 20),
                              n_steps: int = 500,
                              additional_laws: Optional[Dict] = None,
                              method: str = 'RK45',
                              rtol: float = 1e-6,
                              atol: float = 1e-8):
    """
    Simulate reaction-diffusion dynamics on a 2D grid.
    
    Parameters:
    -----------
    rn : ReactionNetwork
        pyCOT reaction network object
    rate : str or list
        Kinetic law name(s)
    grid_shape : tuple, optional
        Grid dimensions (rows, cols). Default: (2, 2)
    D_dict : dict, optional
        Diffusion coefficients {species_name: D_value}
        Default: random uniform [0.01, 0.2]
    x0_dict : dict, optional
        Initial spatial concentrations {species_name: 2D_array}
        Default: random uniform [0, 2.0]
    spec_vector : list, optional
        Kinetic parameters per reaction
    t_span : tuple
        Time interval (t_start, t_end)
    n_steps : int
        Number of time points
    additional_laws : dict, optional
        Custom kinetic functions
    method : str
        Integration method ('RK45', 'LSODA', etc.)
    rtol, atol : float
        Relative and absolute tolerances
        
    Returns:
    --------
    t : array
        Time points (n_steps,)
    X_out : dict
        Spatial concentration fields {species: array(n_steps, rows, cols)}
    flux_out : dict
        Spatial flux fields {reaction: array(n_steps, rows, cols)}
        
    Notes:
    ------
    - Uses periodic boundary conditions (wrap-around)
    - Laplacian approximated with 5-point stencil
    - Diffusion term: D * ∇²x ≈ D * (4 neighbors - 4*center)
    """
    np.random.seed(seed=42)
    species = [specie.name for specie in rn.species()]
    reactions = [reaction.name() for reaction in rn.reactions()]
    
    rate = validate_rate_list(rate, len(reactions))
    
    # Default grid
    if grid_shape is None:
        grid_shape = (2, 2)
    
    rows, cols = grid_shape
    num_cells = rows * cols
    
    # Default diffusion coefficients
    if D_dict is None:
        D_dict = {sp: np.round(np.random.uniform(0.01, 0.2), 3) for sp in species}
    
    # Default initial conditions
    if x0_dict is None:
        x0_dict = {sp: np.round(np.random.uniform(0, 2.0, size=(rows, cols)), 2) 
                   for sp in species}
    
    # State flattening/reshaping utilities
    def flatten_state(x_dict):
        return np.concatenate([x_dict[sp].flatten() for sp in species])
    
    def reshape_state(x):
        return {sp: x[i*num_cells:(i+1)*num_cells].reshape((rows, cols)) 
                for i, sp in enumerate(species)}
    
    x0 = flatten_state(x0_dict)
    
    # Generate kinetic parameters
    if spec_vector is None:
        spec_vector = generate_default_parameters(rate, len(reactions), additional_laws)
    
    # Build kinetic law registry
    rate_laws = dict(KINETIC_REGISTRY)
    if additional_laws:
        rate_laws.update(additional_laws)
    
    # Local reaction dynamics at each grid point
    def reaction_dynamics(Xdict, rate, spec_vector):
        dxdt_dict = {sp: np.zeros((rows, cols)) for sp in species}
        flux_dict = {r: np.zeros((rows, cols)) for r in reactions}
        
        # Iterate over grid points
        for i in range(rows):
            for j in range(cols):
                local_x = [Xdict[sp][i, j] for sp in species]
                
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
                    
                    flux_dict[reaction.name()][i, j] = v_r
                    
                    # Apply stoichiometry
                    for sp_idx, stoich in enumerate(stoichiometry):
                        if stoich != 0:
                            dxdt_dict[species[sp_idx]][i, j] += stoich * v_r
        
        return dxdt_dict, flux_dict
    
    # Diffusion term with Laplacian (periodic boundary conditions)
    def diffusion_term(Xdict):
        def laplacian(D):
            # 5-point stencil with periodic boundaries
            return (-4*D +
                    np.roll(D, 1, axis=0) +   # Up
                    np.roll(D, -1, axis=0) +  # Down
                    np.roll(D, 1, axis=1) +   # Right
                    np.roll(D, -1, axis=1))   # Left
        
        return {sp: D_dict.get(sp, 0.0) * laplacian(Xdict[sp]) for sp in species}
    
    # Combined ODE system
    def combined_ode(t, x):
        Xdict = reshape_state(x)
        dxdt_reac, _ = reaction_dynamics(Xdict, rate, spec_vector)
        dxdt_diff = diffusion_term(Xdict)
        dxdt_total = {sp: dxdt_reac[sp] + dxdt_diff[sp] for sp in species}
        return flatten_state(dxdt_total)
    
    # Integrate
    t_eval = np.linspace(t_span[0], t_span[1], n_steps)
    sol = solve_ivp(combined_ode, t_span, x0, t_eval=t_eval, 
                    method=method, rtol=rtol, atol=atol)
    
    # Format output: spatial fields
    X_out = {sp: np.zeros((n_steps, rows, cols)) for sp in species}
    for i, xt in enumerate(sol.y.T):
        xt_dict = reshape_state(xt)
        for sp in species:
            X_out[sp][i] = xt_dict[sp]
    
    # Compute flux time series
    flux_out = {r: np.zeros((n_steps, rows, cols)) for r in reactions}
    for time_idx, t_val in enumerate(sol.t):
        x_val = sol.y[:, time_idx]
        Xdict = reshape_state(x_val)
        _, flux_dict = reaction_dynamics(Xdict, rate, spec_vector)
        
        for r in reactions:
            flux_out[r][time_idx] = flux_dict[r]
    
    return sol.t, X_out, flux_out


__all__ = ['simulate_spatial_dynamics']
