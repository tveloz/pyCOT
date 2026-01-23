"""
Advanced Deterministic Kinetic Laws

This module contains sophisticated kinetic rate laws for specialized scenarios:
- Time-dependent kinetics (cosine/oscillatory)
- Saturation with multiple substrates
- Threshold-activated kinetics with memory
- Inhibition kinetics
- State-dependent switching

These extend beyond canonical forms for modeling complex real-world phenomena.
"""

import numpy as np


def rate_cosine(reactants, concentrations, species_idx, spec_vector, t=0):
    """
    Cosine Oscillatory Kinetics
    
    Models time-varying inflows or seasonal processes with smooth oscillation.
    Always positive, varies between 0 and amplitude A.
    
    Parameters:
    -----------
    reactants : list of tuple
        (not used, kept for interface consistency)
    concentrations : array-like
        (not used)
    species_idx : dict
        (not used)
    spec_vector : list
        Parameters [A, w] where:
        - A: amplitude (maximum rate)
        - w: angular frequency (2π/period)
    t : float, optional
        Current simulation time
        
    Returns:
    --------
    float
        Instantaneous rate at time t
        
    Example:
    --------
    For seasonal water inflow with 12-month period:
    A = 10 (max inflow)
    w = 2π/12 = 0.524 (monthly frequency)
    rate = 10 * (1 + cos(0.524*t)) / 2
    
    Oscillates between 0 (dry season) and 10 (wet season)
    """
    A, w = spec_vector
    rate = A * (1 + np.cos(w * t / (2 * np.pi))) / 2
    return rate


rate_cosine.expression = lambda reactants, reaction: (
    f"A_{reaction} * (1 + cos(w_{reaction} * t)) / 2"
)


def rate_saturated(substrates, concentrations, species_idx, spec_vector):
    """
    Saturated Multi-Substrate Kinetics
    
    Combines Michaelis-Menten saturation with multi-substrate dependence.
    One substrate exhibits saturation, others enter multiplicatively.
    
    Parameters:
    -----------
    substrates : list of tuple
        List of (species_name, stoichiometric_coefficient) tuples
        First substrate: resource with saturation
        Additional substrates: catalysts/populations (linear dependence)
    concentrations : array-like
        Current concentrations of all species
    species_idx : dict
        Mapping from species names to concentration array indices
    spec_vector : list
        Parameters [Vmax, Km] where:
        - Vmax: maximum rate
        - Km: half-saturation constant for first substrate
        
    Returns:
    --------
    float
        Reaction rate
        
    Examples:
    ---------
    Single substrate (pure saturation):
    rate = Vmax * [S] / (Km + [S])
    
    Two substrates (resource processing by population):
    rate = Vmax * [Resource] * [Population] / (Km + [Resource])
    
    Use case: Conflict model where violence production depends on
    grievances (saturating) and population (linear).
    """
    Vmax, Km = spec_vector
    
    if len(substrates) == 1:
        # Pure saturation
        substrate_name = substrates[0][0]
        substrate_conc = concentrations[species_idx[substrate_name]]
        return (Vmax * substrate_conc) / (Km + substrate_conc)
    
    elif len(substrates) >= 2:
        # Saturation on first, multiply by others
        substrate1_name = substrates[0][0]
        substrate1_conc = concentrations[species_idx[substrate1_name]]
        
        # Multiply by all other substrate concentrations
        other_product = 1.0
        for i in range(1, len(substrates)):
            substrate_i_name = substrates[i][0]
            other_product *= concentrations[species_idx[substrate_i_name]]
        
        return (Vmax * substrate1_conc * other_product) / (Km + substrate1_conc)
    
    else:
        return 0.0


rate_saturated.expression = lambda substrates, reaction: (
    f"(Vmax_{reaction} * [{substrates[0][0]}] * [{substrates[1][0] if len(substrates) > 1 else '1'}]) / "
    f"(Km_{reaction} + [{substrates[0][0]}])"
    if len(substrates) >= 1 else "0"
)


def rate_threshold_memory(reactants, concentrations, species_idx, spec_vector, 
                          memory_state=None, t=0, dt=1.0):
    """
    Threshold-Activated Kinetics with Memory
    
    IMPORTANT FOR CONSERVATION:
    - Threshold detection uses PRODUCT of reactant concentrations
    - Rate magnitude uses appropriate concentration terms for stoichiometry
    
    Reactions activate only when reactant concentrations exceed thresholds.
    Memory accumulates with activation history, influencing future rates.
    Enables modeling of rare events, critical transitions, and hysteresis.
    
    Parameters:
    -----------
    reactants : list of tuple
        List of (species_name, stoichiometric_coefficient) tuples
    concentrations : array-like
        Current concentrations of all species
    species_idx : dict
        Mapping from species names to concentration array indices
    spec_vector : list
        Parameters [threshold, base_rate, memory_strength, decay_rate, steepness]:
        - threshold: concentration threshold for activation
        - base_rate: baseline rate when activated
        - memory_strength: influence of past activations (0-1)
        - decay_rate: memory decay rate (>0)
        - steepness: sigmoid steepness (default 5)
    memory_state : dict, optional
        Current memory value for this reaction (external tracking required)
    t : float, optional
        Current time (for memory updates)
    dt : float, optional
        Time step (for memory integration)
        
    Returns:
    --------
    float
        Reaction rate including memory effects
        
    Mathematical Form:
    ------------------
    rate = base_rate * magnitude * sigmoid(product - threshold, steepness) 
           * (1 + memory_strength * memory)
    
    where:
        product = concentration product for threshold detection
        magnitude = concentration term for rate (respects stoichiometry)
    
    """
    # Robust unpacking with sensible defaults
    if len(spec_vector) >= 5:
        threshold, base_rate, memory_strength, decay_rate, steepness = spec_vector[:5]
    elif len(spec_vector) == 4:
        threshold, base_rate, memory_strength, decay_rate = spec_vector[:4]
        steepness = 5.0
    elif len(spec_vector) == 3:
        threshold, base_rate, memory_strength = spec_vector[:3]
        decay_rate = 1.0
        steepness = 5.0
    elif len(spec_vector) == 2:
        threshold, base_rate = spec_vector[:2]
        memory_strength = 0.0
        decay_rate = 1.0
        steepness = 5.0
    elif len(spec_vector) == 1:
        threshold = 0.0
        base_rate = float(spec_vector[0])
        memory_strength = 0.0
        decay_rate = 1.0
        steepness = 5.0
    else:
        # No parameters provided: use safe defaults
        threshold = 0.0
        base_rate = 1.0
        memory_strength = 0.0
        decay_rate = 1.0
        steepness = 5.0
    
    # Handle case with no reactants
    if not reactants:
        return 0.0
    
    # CRITICAL FIX: Separate product (for threshold) and magnitude (for rate)
    
    if len(reactants) == 1:
        # Single reactant: e.g., 2IP => 2IP + V
        species_name, coeff = reactants[0]
        conc = concentrations[species_idx[species_name]]
        
        # Product for threshold detection
        product = conc ** coeff
        
        # Magnitude for rate (mass action: conc^coeff)
        magnitude = conc ** coeff
        
    else:
        # Multiple reactants: e.g., IP + V => DP + V
        # Product for threshold = product of all reactants
        product = 1.0
        for species_name, coeff in reactants:
            product *= concentrations[species_idx[species_name]] ** coeff
        
        # Magnitude for rate = ONLY FIRST REACTANT (for conservation)
        # This ensures population conservation in compartmental models
        # First reactant is the one being transformed (IP in IP+V=>DP+V)
        # Other reactants are catalysts for threshold detection
        first_species, first_coeff = reactants[0]
        magnitude = concentrations[species_idx[first_species]] ** first_coeff
    
    # Sigmoid activation based on product
    activation = 1.0 / (1.0 + np.exp(-steepness * (product - threshold)))
    
    # Memory contribution (default to 0 if not tracked)
    memory_factor = 1.0
    if memory_state is not None:
        memory_value = memory_state.get('value', 0.0)
        memory_factor = 1.0 + memory_strength * memory_value
    
    # Combined rate: INCLUDES MAGNITUDE TERM!
    rate = base_rate * magnitude * activation * memory_factor
    if reactants:
        min_available = float('inf')
        for species_name, coeff in reactants:
            conc = concentrations[species_idx[species_name]]
            if coeff > 0 and conc > 0:
                # Maximum rate this species can support (assuming it's consumed)
                available_rate = conc / coeff
                min_available = min(min_available, available_rate)
        
        # Cap the rate at what's physically available
        if min_available < float('inf'):
            rate = min(rate, min_available)


    return rate


rate_threshold_memory.expression = lambda reactants, reaction: (
    f"base_rate_{reaction} * [{reactants[0][0]}]^{reactants[0][1]} * "
    f"sigmoid(product - threshold_{reaction}) * (1 + memory_strength_{reaction} * memory_{reaction})"
    if len(reactants) == 1 else
    f"base_rate_{reaction} * [{reactants[0][0]}]^{reactants[0][1]} * "
    f"sigmoid({' * '.join(f'[{r[0]}]^{r[1]}' if r[1] != 1 else f'[{r[0]}]' for r in reactants)} - threshold_{reaction}) * "
    f"(1 + memory_strength_{reaction} * memory_{reaction})"
    if reactants else "0"
)


def rate_inhibition(substrates, concentrations, species_idx, spec_vector, inhibitor_type='competitive'):
    """
    Inhibition Kinetics (Competitive and Non-Competitive)
    
    Models reactions with inhibitory species that reduce reaction rates.
    
    Parameters:
    -----------
    substrates : list of tuple
        [(substrate, coef), (inhibitor, coef)]
        First is substrate, second is inhibitor
    concentrations : array-like
        Current concentrations
    species_idx : dict
        Species index mapping
    spec_vector : list
        [Vmax, Km, Ki] where Ki is inhibition constant
    inhibitor_type : str
        'competitive' or 'noncompetitive'
        
    Returns:
    --------
    float
        Inhibited reaction rate
        
    Forms:
    ------
    Competitive: rate = Vmax*[S] / (Km*(1 + [I]/Ki) + [S])
    Non-competitive: rate = Vmax*[S] / ((Km + [S])*(1 + [I]/Ki))
    """
    if len(substrates) < 2 or len(spec_vector) < 3:
        return 0.0
        
    Vmax, Km, Ki = spec_vector[:3]
    
    substrate_name = substrates[0][0]
    inhibitor_name = substrates[1][0]
    
    S = concentrations[species_idx[substrate_name]]
    I = concentrations[species_idx[inhibitor_name]]
    
    if inhibitor_type == 'competitive':
        return Vmax * S / (Km * (1 + I/Ki) + S)
    else:  # non-competitive
        return Vmax * S / ((Km + S) * (1 + I/Ki))


rate_inhibition.expression = lambda substrates, reaction: (
    f"Vmax_{reaction} * [{substrates[0][0]}] / "
    f"(Km_{reaction}*(1+[{substrates[1][0]}]/Ki_{reaction}) + [{substrates[0][0]}])"
    if len(substrates) >= 2 else "0"
)


# Registry of advanced kinetic laws
ADVANCED_KINETICS = {
    'cosine': rate_cosine,
    'saturated': rate_saturated,
    'threshold_memory': rate_threshold_memory,
    'inhibition': rate_inhibition
}


__all__ = [
    'rate_cosine',
    'rate_saturated',
    'rate_threshold_memory',
    'rate_inhibition',
    'ADVANCED_KINETICS'
]
