import numpy as np 
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

###################################################################################
# Function to calculate the reaction rate
def reaction_rate(reaction_name, species_concentrations, RN_dict, k):
    """
    Calculate the reaction rate for a given reaction.

    Parameters:
    reaction_name (str): The name of the reaction.
    species_concentrations (dict): Current concentrations of species.
    RN_dict (dict): Reaction network dictionary mapping reactions to stoichiometry.
    k (dict): Dictionary of reaction rate constants.

    Returns:
    float: Reaction rate for the specified reaction. 
    """
    reactants_coeff, _ = RN_dict[reaction_name]
    rate = k[reaction_name]
    
    # Multiply the concentrations of reactants raised to their stoichiometric coefficients
    for reactant, coef in reactants_coeff:
        rate *= species_concentrations.get(reactant, 0) ** coef
    
    return rate

###################################################################################
# Generate ODEs for the reaction network
def generate_odes(species_concentrations, RN_dict, k):
    """
    Generate the ODE system for the reaction network.

    Parameters:
    species_concentrations (dict): Mapping of species to indices.
    RN_dict (dict): Reaction network dictionary mapping reactions to stoichiometry.
    k (dict): Dictionary of reaction rate constants.

    Returns:
    function: A function representing the ODE system.
    """
    def odes(t, y):
        # Initialize derivatives for all species
        dSdt = {species: 0 for species in species_concentrations.keys()}
        
        # Map current species concentrations
        current_concentrations = dict(zip(species_concentrations.keys(), y))
        
        # Compute derivatives for each species
        for species in current_concentrations:
            for reaction_name, (reactants, products) in RN_dict.items():
                # Reactants contribution
                for reactant, coef in reactants:
                    if reactant == species:
                        dSdt[species] -= coef * reaction_rate(reaction_name, current_concentrations, RN_dict, k)
                # Products contribution
                for product, coef in products:
                    if product == species:
                        dSdt[species] += coef * reaction_rate(reaction_name, current_concentrations, RN_dict, k)
        
        # Ensure all species are accounted for in dSdt
        for key in species_concentrations.keys():
            if key not in dSdt:
                raise KeyError(f"Unexpected key '{key}' not found in dSdt.")
        
        return np.array([dSdt[species] for species in species_concentrations.keys()])
    
    return odes

###################################################################################
# Solve the ODE system and plot results
def solve_ode(RN, k, t_span, y0,n_steps=100):
    """
    Solve the ODE system for the reaction network and plot the results.

    Parameters:
    RN: Reaction network object containing species and reaction details.
    k (dict): Dictionary of reaction rate constants.
    t_span (tuple): Time span for the simulation (start, end).
    y0 (list): Initial concentrations of species.
    """
    # Build the reaction network dictionary
    RN_dict = {}
    for i in range(len(RN.RnStr)):
        reactants_coeff = []
        products_coeff = []
        for j in range(len(RN.SpStr)):
            coef_r = RN.RnMsupp[i][j]
            if coef_r > 0:
                reactants_coeff.append((RN.SpStr[j], coef_r))
            coef_p = RN.RnMprod[i][j]
            if coef_p > 0:
                products_coeff.append((RN.SpStr[j], coef_p))
        RN_dict[RN.RnStr[i]] = (reactants_coeff, products_coeff)
    
    # Map species to indices
    species_concentrations = {species: idx for idx, species in enumerate(RN.SpStr)}
    
    # Generate ODE system
    odes = generate_odes(species_concentrations, RN_dict, k)
    
    # Solve the ODE system
    sol = solve_ivp(odes, t_span, y0, t_eval=np.linspace(t_span[0], t_span[1], n_steps))
    return sol