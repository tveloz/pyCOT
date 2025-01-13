import numpy as np 
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt 
import pandas as pd

###################################################################################
# Function to calculate the reaction rate 
def reaction_rate_mak(reaction_name, species_concentrations, RN_dict, k):
    """
    Calculate the reaction rate for a given reaction.
    """
    if k is None:
        raise ValueError("The speed constants ‘k’ are not correctly defined.")
    
    # print("Constantes de velocidad:", k)  # Verifica el valor de k
    reactants_coeff, _ = RN_dict[reaction_name]
    
    rate = k[reaction_name]
    
    # Multiply the concentrations of reactants raised to their stoichiometric coefficients
    for reactant, coef in reactants_coeff:
        rate *= species_concentrations.get(reactant, 0) ** coef
    
    return rate


###################################################################################
# Generate ODEs for the reaction network
def generate_odes_mak(species_concentrations, RN_dict, k):
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
                        dSdt[species] -= coef * reaction_rate_mak(reaction_name, current_concentrations, RN_dict, k)
                # Products contribution
                for product, coef in products:
                    if product == species:
                        dSdt[species] += coef * reaction_rate_mak(reaction_name, current_concentrations, RN_dict, k)
        
        # Ensure all species are accounted for in dSdt
        for key in species_concentrations.keys():
            if key not in dSdt:
                raise KeyError(f"Unexpected key '{key}' not found in dSdt.")
        
        return np.array([dSdt[species] for species in species_concentrations.keys()])
    
    return odes

###################################################################################
# Solve the ODE system
def simulate_ode_mak(RN, x0 = None, t_span = None, n_steps = None, k0 = None):
    """
    Solve the ODE system for the given reaction network and return the results as a time series.

    This function constructs the reaction network, generates the ODE system based on the stoichiometry 
    of the reactions, and solves the system using an ODE solver. It returns the time series of species 
    concentrations over the specified time span.

    Parameters:
    RN (object): Reaction network object containing species and reaction details.
    x0 (list, optional): Initial concentrations of species. If not provided, they will be generated automatically.
    k0 (dict, optional): Dictionary of reaction rate constants. If not provided, they will be generated automatically.
    t_span (tuple, optional): Time span for the simulation (start, end). Defaults to (0, 100).
    n_steps (int, optional): Number of time steps for the simulation. Defaults to 100.

    Returns:
    pd.DataFrame: A pandas DataFrame containing the time series of species concentrations over time.
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

    # Define the initial state vector
    if x0 is None:
        x0 = generate_random_vector(len(RN.SpStr), seed=3).round(1)

    # Number of steps
    if n_steps is None:
        n_steps = 100
        
    if t_span is None:    
        t_span = [0, 100]

    # Automatically generate rate constants if not provided
    if k0 is None: 
        k0 = generate_random_vector(len(RN.RnStr), seed=3).round(2) 
        k = {RN.RnStr[i]: float(k0[i]) for i in range(len(RN.RnStr))}        
    # Asegúrate de que `k` siempre esté definido
    k = {RN.RnStr[i]: float(k0[i]) for i in range(len(RN.RnStr))}

    # Generate ODE system
    odes = generate_odes_mak(species_concentrations, RN_dict, k)

    # Solve the ODE system
    sol = solve_ivp(odes, t_span, x0, t_eval=np.linspace(t_span[0],t_span[1], n_steps))

    # Extract the vector of times and the solutions
    data_concentrations = {"Time": sol.t}
    for i, species in enumerate(RN.SpStr):
        data_concentrations[species] = sol.y[i]

    time_series = pd.DataFrame(data_concentrations)
    # return time_series

    # Calculate flux vector at each time step
    flux_vector = []
    for t_idx in range(len(sol.t)):
        concentrations = sol.y[:, t_idx]
        flux = []
        for rxn, (reactants, products) in RN_dict.items():
            rate = k[rxn]
            for species, coeff in reactants:
                rate *= concentrations[species_concentrations[species]] ** coeff
            flux.append(rate)
        flux_vector.append(flux)

    # Calculate flux vector at each time step
    flux_data = {"Time": sol.t}
    for i, rxn in enumerate(RN.RnStr):
        flux_data[f"Flux_{i}"] = []

    for t_idx in range(len(sol.t)):
        concentrations = sol.y[:, t_idx]
        for i, rxn in enumerate(RN.RnStr):
            rate = k[rxn]
            reactants = RN_dict[rxn][0]
            for species, coeff in reactants:
                rate *= concentrations[species_concentrations[species]] ** coeff
            flux_data[f"Flux_{i}"].append(rate)

    flux_vector = pd.DataFrame(flux_data)
    
    return time_series, flux_vector    
###################################################################################








###################################################################################
# Function to generate the state vector with random positive values
# def generate_random_vector(n_species):
#     """
#     Generates a state vector with random values between 0 and 1 for each species.
    
#     Parameters:
#     n_species (int): Number of species in the system.

#     Returns:
#     np.ndarray: A vector of size n_species with values between 0 and 1.
#     """
#     return np.random.rand(n_species)  # Generates a vector of size n_species with values between 0 and 1

def generate_random_vector(n_species, seed=None):
    """
    Generates a state vector with random values between 0 and 1 for each species.
    The same seed will always produce the same vector.

    Parameters:
    n_species (int): Number of species in the system.
    seed (int or None): Seed for random number generator. If None, the seed is not set.

    Returns:
    np.ndarray: A vector of size n_species with values between 0 and 1.
    """
    if seed is not None:
        np.random.seed(seed)  # Set the seed for reproducibility
    return np.random.rand(n_species)  # Generates a vector of size n_species with values between 0 and 1

###################################################################################
# Function to create the stoichiometric matrix 
def universal_stoichiometric_matrix(RN):
    """
    Creates the stoichiometric matrix based on the reaction network (RN).
    
    Parameters:
    RN (object): The reaction network object that contains species and reaction details.
    
    Returns:
    np.ndarray: The stoichiometric matrix representing the reactions, with integer values.
    """
    # Initialize the stoichiometric matrix with zeros
    matrix = np.zeros((len(RN.SpStr), len(RN.RnStr)))  # #species x #reactions

    # Iterate over reactions
    for i in range(len(RN.RnStr)):  # For each reaction
        # Process reactants
        for j in range(len(RN.SpStr)):  # For each species
            coef_r = RN.RnMsupp[i][j]  # Coefficient of the reactant in reaction i for species j
            if coef_r > 0:  # If greater than 0, it's a reactant
                matrix[j, i] = -coef_r  # Reactant coefficients are negative (transposition here)

        # Process products
        for j in range(len(RN.SpStr)):  # Same for products
            coef_r = RN.RnMprod[i][j]  # Coefficient of the product in reaction i for species j
            if coef_r > 0:  # If greater than 0, it's a product
                matrix[j, i] = coef_r  # Product coefficients are positive (transposition here)

    # Convert the matrix to integer type
    matrix = matrix.astype(int)  # Convert all elements to integers

    return matrix  # Returns the stoichiometric matrix with integer values


###################################################################################
# Function to update the state vector at each time step 
def simulate_discrete_random(RN, S, x, n_iter=10):
    """
    Iterates the state vector over a number of iterations, updating the state 
    based on the stoichiometric matrix and random values, and stores the flux history.

    Parameters:
    RN (object): The reaction network object.
    S (np.ndarray): The stoichiometric matrix.
    x (np.ndarray): The initial state vector. 
    n_iter (int): Number of iterations to run.
    
    Returns:
    df_state: DataFrame with the time and updated species concentrations.
    df_flux: DataFrame with the time and flux vectors v(x,k) for each iteration.
    """
    history = [x.copy()]  # Store the initial state in the history
    flux_history = []  # Store the flux vectors (v)
    
    for i in range(n_iter):  # Iterate for n_iter steps
        v = generate_random_vector(S.shape[1])  # Generate a new random state vector for the reaction network
        x = x + S @ v  # Update the state vector using the stoichiometric matrix and random values
        x = np.maximum(x, 0)  # Ensure all values are positive (no negative concentrations)
        
        history.append(x.copy())  # Store the updated state in the history
        flux_history.append(v.copy())  # Store the flux vector
    
    # Convert the history into a numpy array
    history = np.array(history)
    
    # Create a time vector t
    t = np.arange(n_iter)  # Time steps based on the history shape
    
    # Concatenate the time vector with the history of states
    history = np.column_stack((t, history[:n_iter]))  # Combine time with the state history
    
    # Create DataFrame for species concentrations
    columns = ['Time'] + RN.SpStr  # The first column is 'Time', followed by species names
    df_state = pd.DataFrame(history, columns=columns)  # Create the DataFrame with time and species concentrations
    
    # Create DataFrame for flux history
    flux_columns = ['Time'] + [f'Flux_{i}' for i in range(S.shape[1])]  # Column names for fluxes
    flux_history = np.column_stack((np.arange(0, n_iter), flux_history))  # Combine time with flux history
    df_flux = pd.DataFrame(flux_history, columns=flux_columns)  # Create DataFrame for fluxes

    return df_state, df_flux # Return the state and flux DataFrames
