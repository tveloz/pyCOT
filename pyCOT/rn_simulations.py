import numpy as np
import pandas as pd

###################################################################################
# Function to generate the state vector with random positive values
def generate_state_vector(n_species):
    """
    Generates a state vector with random values between 0 and 1 for each species.
    
    Parameters:
    n_species (int): Number of species in the system.

    Returns:
    np.ndarray: A vector of size n_species with values between 0 and 1.
    """
    return np.random.rand(n_species)  # Generates a vector of size n_species with values between 0 and 1

###################################################################################
# Function to create the stoichiometric matrix
def stoichiometric_matrix(RN):
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
def iterate_state_vector(S, x, species, n_iter=10):
    """
    Iterates the state vector over a number of iterations, updating the state 
    based on the stoichiometric matrix and random values, and stores the flux history.

    Parameters:
    S (np.ndarray): The stoichiometric matrix.
    x (np.ndarray): The initial state vector.
    species (list): List of species names.
    n_iter (int): Number of iterations to run.

    Returns:
    df_state: DataFrame with the time and updated species concentrations.
    df_flux: DataFrame with the time and flux vectors v(x,k) for each iteration.
    """
    history = [x.copy()]  # Store the initial state in the history
    flux_history = []  # Store the flux vectors (v)
    
    for i in range(n_iter):  # Iterate for n_iter steps
        v = generate_state_vector(S.shape[1])  # Generate a new random state vector for the reaction network
        x = x + S @ v  # Update the state vector using the stoichiometric matrix and random values
        x = np.maximum(x, 0)  # Ensure all values are positive (no negative concentrations)
        
        history.append(x.copy())  # Store the updated state in the history
        flux_history.append(v.copy())  # Store the flux vector
    
    # Convert the history into a numpy array
    history = np.array(history)
    
    # Create a time vector t
    t = np.arange(history.shape[0])  # Time steps based on the history shape
    
    # Concatenate the time vector with the history of states
    history = np.column_stack((t, history))  # Combine time with the state history
    
    # Create DataFrame for species concentrations
    columns = ['Time'] + species  # The first column is 'Time', followed by species names
    df_state = pd.DataFrame(history, columns=columns)  # Create the DataFrame with time and species concentrations
    
    # Create DataFrame for flux history
    flux_columns = ['Time'] + [f'Flux_{i}' for i in range(S.shape[1])]  # Column names for fluxes
    flux_history = np.column_stack((np.arange(n_iter), flux_history))  # Combine time with flux history
    df_flux = pd.DataFrame(flux_history, columns=flux_columns)  # Create DataFrame for fluxes

    return df_state, df_flux  # Return the DataFrame with the time series of species states
