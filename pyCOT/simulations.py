# Import libraries
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt 
import re
from scipy.integrate import solve_ivp
from scipy.integrate import odeint
from matplotlib.widgets import Slider
from scipy.optimize import linprog  
import matplotlib.animation as animation  

# -----------------------------
# USER-DEFINED KINETIC FUNCTIONS
# -----------------------------
# Function to calculate the rate of reaction of the type 'mmk'
def rate_mak(reactants, concentrations, species_idx, spec_vector):
    """
    Calculate the rate of reaction using the mass action kinetics.
    The mass action kinetics describes the rate of a reaction as a function of the concentrations of the reactants.
    Parameters:
    - reactants: List of tuples (species, stoichiometric coefficient) for the reactants.
    - concentrations: List of current concentrations of the species.
    - species_idx: Dictionary mapping species names to their indices in the concentrations list.
    - spec_vector: List of parameters for the reaction, including the rate constants k.
    Returns:
    - rate: The rate of the reaction.
    """
    # Extract the rate constant and the stoichiometric coefficients
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

# Function to calculate the rate of reaction of the type 'mmk'
def rate_mmk(reactants, concentrations, species_idx, spec_vector):
    """
    Calculate the rate of reaction using the Michaelis-Menten equation.
    The Michaelis-Menten equation describes the rate of enzymatic reactions as a function of substrate concentration.
    Parameters:
    - reactants: List of tuples (species, stoichiometric coefficient) for the reactants.
    - concentrations: List of current concentrations of the species.
    - species_idx: Dictionary mapping species names to their indices in the concentrations list.
    - spec_vector: List of parameters for the reaction, including Vmax and Km.
    Returns:
    - rate: The rate of the reaction.
    """
    Vmax, Km = spec_vector
    s = reactants[0][0]
    S = concentrations[species_idx[s]]
    return Vmax * S / (Km + S)
rate_mmk.expression = lambda reactants, reaction: (
    f"(Vmax_{reactants[0][0]} * [{reactants[0][0]}]) / "
    f"(Km_{reactants[0][0]} + [{reactants[0][0]}])"
    if reactants else "0  (mmk without defined substrate)"
)

# Function to calculate the rate of reaction of the type 'hill'
def rate_hill(substrate, concentrations, species_idx, spec_vector):
    """
    Calculate the rate of reaction using the Hill equation.
    The Hill equation is a generalization of the Michaelis-Menten equation.
    Parameters:
    - substrate: List of tuples (species, stoichiometric coefficient) for the substrate.
    - concentrations: List of current concentrations of the species.
    - species_idx: Dictionary mapping species names to their indices in the concentrations list.
    - spec_vector: List of parameters for the reaction, including Vmax, K, and n.
    Returns:
    - rate: The rate of the reaction.
    """
    Vmax, K, n = spec_vector

    if isinstance(substrate, list) and isinstance(substrate[0], tuple):
        substrate_names = [s[0] for s in substrate]
    elif isinstance(substrate, list):
        substrate_names = substrate
    else:
        substrate_names = [substrate]

    S = sum(concentrations[species_idx[s]] for s in substrate_names)
    return Vmax * (S**n) / (K**n + S**n)
rate_hill.expression = lambda reactants, reaction: (
    f"(Vmax_{reactants[0][0]} * [{reactants[0][0]}]^n) / "
    f"(Km_{reactants[0][0]}^n + [{reactants[0][0]}]^n)"
    if reactants else "0  (hill without defined substrate)"
)

# Function to calculate the rate of reaction of the type 'ping_pong' 
def rate_ping_pong(substrates, concentrations, species_idx, spec_vector):
    """
    Calculate the rate of reaction using the ping-pong mechanism.
    The ping-pong mechanism is a type of enzyme kinetics where the enzyme alternates between two states.
    It assumes two substrates and a single product.
    Parameters:
    - substrates: List of tuples (species, stoichiometric coefficient) for the substrates.
    - concentrations: List of current concentrations of the species.
    - species_idx: Dictionary mapping species names to their indices in the concentrations list.
    - spec_vector: List of parameters for the reaction, including Vmax, KmA, and KmB.
    Returns:
    - rate: The rate of the reaction.
    """
    # Extract the parameters from the spec_vector
    Vmax, KmA, KmB = spec_vector
    if len(substrates) < 2:
        return 0  # No suficientes sustratos definidos
    
    substrateA = substrates[0][0]
    substrateB = substrates[1][0]
    A = concentrations[species_idx[substrateA]]
    B = concentrations[species_idx[substrateB]]
    return Vmax * A * B / (KmA * B + KmB * A + A * B)
rate_ping_pong.expression = lambda substrates, reaction: (
    f"(Vmax_{reaction} * [{substrates[0][0]}] * [{substrates[1][0]}]) / "
    f"(Km_{substrates[0][0]} * [{substrates[1][0]}] + Km_{substrates[1][0]} * [{substrates[0][0]}] + [{substrates[0][0]}] * [{substrates[1][0]}])"
    if len(substrates) >= 2 else "0 (ping-pong requires two substrates)"
)

def update_rate_laws(rate_list, additional_laws=None):
    # Diccionario base
    rate_laws = {
        'mak': rate_mak,
        'mmk': rate_mmk,
        'hill': rate_hill
    }

    # Añadir cinéticas adicionales si se entregan
    if additional_laws:
        for name in rate_list:
            if name not in rate_laws:
                if name in additional_laws:
                    rate_laws[name] = additional_laws[name]
                else:
                    raise NotImplementedError(f"La cinética '{name}' no está definida ni registrada.")

    return rate_laws 

# Function to generate a random vector of size n
def generate_random_vector(n, seed=None):
    rng = np.random.default_rng(seed)
    return rng.random(n)

# Function to validate the rate list
def validate_rate_list(rate, n_reactions):
    if isinstance(rate, str):
        return [rate] * n_reactions
    elif len(rate) != n_reactions:
        raise ValueError(f"The size of the list rate={len(rate)} must be equal to the number of reactions={n_reactions}.")
    return rate

def build_reaction_dict(rn):
    rn_dict = {} 
    species =   rn.stoichiometry_matrix().species
    reactions = rn.stoichiometry_matrix().reactions
    reactants_vectors = rn.reactants_matrix().T
    products_vectors = rn.products_matrix().T
    for i in range(len(reactions)):
        reactants = [(species[j], reactants_vectors[i][j]) for j in range(len(species)) if reactants_vectors[i][j] > 0]
        products  = [(species[j], products_vectors[i][j]) for j in range(len(species)) if products_vectors[i][j] > 0]
        rn_dict[reactions[i]] = (reactants, products) 
    return rn_dict



# Function to parse the parameters from the reaction network
def parse_parameters(rn, rn_dict, rate, spec_vector):
    reactions = rn.stoichiometry_matrix().reactions
    parameters = {}
    for i, reaction in enumerate(reactions):
        parameters[reaction] = spec_vector[i]
    return parameters

# Function to compute the reaction rates 
def compute_reaction_rates(rn, rn_dict, rate, parameters, rate_laws):
    species = rn.stoichiometry_matrix().species
    reactions = rn.stoichiometry_matrix().reactions
    species_idx = {s: i for i, s in enumerate(species)}

    def rates_fn(x, t):
        dxdt = np.zeros_like(x)
        concentrations = x
        rates = np.zeros(len(reactions))

        for i, reaction in enumerate(reactions):
            kinetic = rate[i]
            reactants = rn_dict[reaction][0]
            param_vector = parameters[reaction]
            rate_fn = rate_laws[kinetic]
            rates[i] = rate_fn(reactants, concentrations, species_idx, param_vector)

        for i, reaction in enumerate(reactions):
            for sp, coef in rn_dict[reaction][0]:
                dxdt[species_idx[sp]] -= coef * rates[i]
            for sp, coef in rn_dict[reaction][1]:
                dxdt[species_idx[sp]] += coef * rates[i]

        return dxdt

    return rates_fn

def print_differential_equations(rn, rn_dict):
    print("\nDifferential equations of the system:")
    species = rn.stoichiometry_matrix().species
    reactions = rn.stoichiometry_matrix().reactions 
    for sp in species:
        term_dict = {}
        for i, reaction in enumerate(reactions):
            rate_expr = f"v{i+1}"
            coef = 0
            for reactant, coeff in rn_dict[reaction][0]:
                if reactant == sp:
                    coef -= coeff
            for product, coeff in rn_dict[reaction][1]:
                if product == sp:
                    coef += coeff
            if coef != 0:
                term_dict[rate_expr] = coef

        eq_terms = []
        for v, coef in term_dict.items():
            if coef == 1:
                eq_terms.append(f"+{v}")
            elif coef == -1:
                eq_terms.append(f"-{v}")
            elif coef > 1:
                eq_terms.append(f"+{coef}*{v}")
            else:
                eq_terms.append(f"{coef}*{v}")  # coef negativo

        rhs = " ".join(eq_terms)
        print(f"d[{sp}]/dt = {rhs}")

def print_velocity_expressions(rn, rn_dict, rate, additional_laws=None):
    reactions = rn.stoichiometry_matrix().reactions
    print("\nAnalytical expressions of velocities v(x,param):")

    # Dictionary of rate laws
    rate_laws = {
        'mak': rate_mak,
        'mmk': rate_mmk,
        'hill': rate_hill
    }
    # Register the rate laws new
    rate_laws = update_rate_laws(rate, rate_laws if additional_laws is None else {**rate_laws, **additional_laws})

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

# Function to create the time series DataFrame
def time_series(rn, result, time_points):
    species = rn.stoichiometry_matrix().species
    data_concentrations = {"Time": time_points}
    for i, species in enumerate(species):
        data_concentrations[species] = result[:, i]
    df = pd.DataFrame(data_concentrations)
    time_series = df[["Time"] + sorted(df.columns.difference(["Time"]))]
    return time_series 

# Function to create the flux vector DataFrame
def flux_vector(rn, rn_dict, rate, result, time_points, spec_vector, additional_laws=None):
    species = rn.stoichiometry_matrix().species
    reactions = rn.stoichiometry_matrix().reactions
    species_idx = {s: i for i, s in enumerate(species)}
    
    # Validar y expandir lista de tipos de cinéticas
    rate = validate_rate_list(rate, len(reactions))
    
    # Registrar las leyes de velocidad
    rate_laws = update_rate_laws(rate, additional_laws)

    # Parsear parámetros para cada reacción
    parameters = parse_parameters(rn, rn_dict, rate, spec_vector)

    # Crear estructura para guardar flujos
    flux_data = {"Time": time_points}
    for i in range(len(reactions)):
        flux_data[f"Flux_r{i+1}"] = []

    # Calcular flujos en cada punto de tiempo
    for concentrations in result:
        flux_row = []
        for i, reaction in enumerate(reactions):
            kinetic = rate[i]
            reactants = rn_dict[reaction][0]
            param_vector = parameters[reaction]
            rate_fn = rate_laws[kinetic]
            rate_value = rate_fn(reactants, concentrations, species_idx, param_vector)
            flux_row.append(rate_value)

        # Guardar los valores de flujo
        for i, value in enumerate(flux_row):
            flux_data[f"Flux_r{i+1}"].append(value)

    # Convertir a DataFrame
    flux_vector = pd.DataFrame(flux_data)
    return flux_vector

# Function to simulate the reaction network
def simulation(rn, rate='mak', spec_vector=None, x0=None, t_span=(0, 100), n_steps=300, additional_laws=None):
    species = rn.stoichiometry_matrix().species
    reactions = rn.stoichiometry_matrix().reactions
    rn_dict = build_reaction_dict(rn)

    rate = validate_rate_list(rate, len(reactions))

    # Initial conditions vector by default random
    if x0 is None:
        x0 = generate_random_vector(len(species), seed=3).round(1) * 20
        print(f"\nx0 = {x0.tolist()}")
    else:
        print(f"\nx0 = {x0}") 

    if spec_vector is None:
        spec_vector = []
        for kinetic in rate:
            if kinetic == 'mak':
                params = np.round(np.random.uniform(0.1, 1.0, 1), 3)
            elif kinetic == 'mmk':
                params = np.round(np.random.uniform(0.1, 2.0, 2), 3)
            elif kinetic == 'hill':
                Vmax = np.random.uniform(0.5, 2.0)
                Kd = np.random.uniform(0.1, 1.5)
                n = np.random.randint(1, 4)
                params = np.round([Vmax, Kd, n], 3)
            elif kinetic in (additional_laws or {}):  # soporte para cinéticas extra
                # Generar 3 parámetros aleatorios genéricos por defecto
                params = np.round(np.random.uniform(0.1, 1.0, 3), 3)
            else:
                raise ValueError(f"Unknown kinetic law: {kinetic}")
            spec_vector.append(params.tolist())


        parameters = parse_parameters(rn, rn_dict, rate, spec_vector) 
        print("spec_vector =", spec_vector)
    else:
        parameters = parse_parameters(rn, rn_dict, rate, spec_vector)
        print("spec_vector =", spec_vector)

    rate_laws = {
        'mak': rate_mak,
        'mmk': rate_mmk,
        'hill': rate_hill
    }
    # Register the rate laws
    rate_laws = update_rate_laws(rate, rate_laws if additional_laws is None else {**rate_laws, **additional_laws})

    rates_fn = compute_reaction_rates(rn, rn_dict, rate, parameters, rate_laws)

    print_differential_equations(rn, rn_dict)
    print_velocity_expressions(rn, rn_dict, rate, additional_laws)

    time_points = np.linspace(t_span[0], t_span[1], n_steps)
    result = odeint(rates_fn, x0, time_points)

    time_series_df = time_series(rn, result, time_points) 
    flux_vector_df = flux_vector(rn, rn_dict, rate, result, time_points, spec_vector, additional_laws)
    return time_series_df, flux_vector_df

###############################################################################
#################################################################################
#########################################################################
# Function to create the stoichiometric matrix 
def universal_stoichiometric_matrix(rn):
    """
    Creates the stoichiometric matrix based on the reaction network (RN).
    
    Parameters:
    RN (object): The reaction network object that contains species and reaction details.
    
    Returns:
    np.ndarray: The stoichiometric matrix representing the reactions, with integer values.
    """
    species =   rn.stoichiometry_matrix().species
    reactions = rn.stoichiometry_matrix().reactions
    reactants_vectors, products_vectors=build_stoichiometric_vectors(file_path, species)
    # Initialize the stoichiometric matrix with zeros
    matrix = np.zeros((len(species), len(reactions)))  # #species x #reactions

    # Iterate over reactions
    for i in range(len(reactions)):  # For each reaction 
        for j in range(len(species)):  # Same for products
            coef_p = products_vectors[i][j]  # Coefficient of the product in reaction i for species j
            coef_r = reactants_vectors[i][j]  # Coefficient of the reactant in reaction i for species j
            matrix[j, i] = coef_p - coef_r  # Product coefficients are positive (transposition here)  

    return matrix  # Returns the stoichiometric matrix with integer values

# Function to update the state vector at each time step 
def simulate_discrete_random(rn, S, x, n_iter=10):
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
    species =   rn.stoichiometry_matrix().species  # List of species names 
    
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
    columns = ['Time'] + species  # The first column is 'Time', followed by species names
    df_state = pd.DataFrame(history, columns=columns)  # Create the DataFrame with time and species concentrations
    
    # Create DataFrame for flux history
    flux_columns = ['Time'] + [f'Flux_{i}' for i in range(S.shape[1])]  # Column names for fluxes
    flux_history = np.column_stack((np.arange(0, n_iter), flux_history))  # Combine time with flux history
    df_flux = pd.DataFrame(flux_history, columns=flux_columns)  # Create DataFrame for fluxes

    return df_state, df_flux # Return the state and flux DataFrames


###################################################################################
###################################################################################
###################################################################################
###########################################################################################
# Function that computes the time series of a metapopulation reaction-diffusion model
def simulate_odes_metapop_mak(RN, x0=None, t_span=None, n_steps=None, k=None, exchange_rates=None, D=None, grid_shape=None):
    """
    Simulates the dynamics of a spatiotemporal metapopulation with reaction-diffusion terms.
    Specifically, it models the evolution of species in a reaction network distributed over a 
    discrete space (a grid of patches), considering both chemical interactions within each patch 
    and exchange between patches through diffusion and exchange rates.

    Parameters:
        RN (dict): Object describing the reaction network.
        x0 (array-like, optional): Initial conditions of the species. Default is None.
        t_span (tuple, optional): Time range for the simulation. Default is None.
        n_steps (int, optional): Number of steps for the simulation. Default is None.
        k (array-like, optional): Reaction rate constants. Default is None.
        exchange_rates (array-like, optional): Exchange rates between modules (diffusion between nodes). Default is None.
        D (array-like, optional): Spatial diffusion coefficients for each species. Default is None.
        grid_shape (tuple, optional): Shape of the spatial grid (rows, columns). Default is None.

    Returns:
        dict: Simulation results for each patch (key: patch index).
    """
    # Validate that RN and its attributes are not None
    if RN is None or not hasattr(RN, 'SpStr') or not hasattr(RN, 'RnStr'):
        raise ValueError("The reaction network (RN) is invalid or incomplete.")
    
    # Default parameters if not specified
    if grid_shape is None:
        grid_shape = (2, 2)                   # Default configuration if not provided
    if x0 is None:
        x0 = np.random.rand(grid_shape[0], grid_shape[1], len(RN.SpStr)) * 2 # Random initial conditions
    if n_steps is None:
        n_steps = 500                         # Default number of steps    
    if t_span is None:
        t_span = (0, 100)                      # Default simulation time
    if k is None:
        k = np.random.rand(len(RN.RnStr)) * 1 # Random rate constants
    if exchange_rates is None:
        exchange_rates = np.ones((grid_shape[0], grid_shape[1], len(RN.SpStr)))  # No exchange by default 
    if D is None:
        D = np.random.rand(len(RN.SpStr))     # Default uniform diffusion coefficients
    
    # Function describing the system of differential equations with diffusion
    def ode_system(x_flat, t, k, RN, exchange_rates, D, grid_shape):
        # Reconstruct the spatial grid
        x = x_flat.reshape(grid_shape[0], grid_shape[1], -1)
        dxdt = np.zeros_like(x)
        
        # Stoichiometric matrices
        S = np.array(RN.RnMsupp)
        P = np.array(RN.RnMprod)
        
        # Reaction terms
        for i in range(len(RN.SpStr)):      # For each species
            for r in range(len(RN.RnStr)):  # For each reaction
                rate = k[r] * np.prod(x[:, :, :] ** S[r, :], axis=-1)
                dxdt[:, :, i] += rate * (P[r, i] - S[r, i])
        
        # Diffusion terms (diffusion in the grid)
        for i in range(len(RN.SpStr)):      # For each species
            dxdt[:, :, i] += D[i] * (
                np.roll(x[:, :, i],  1, axis=0) +  # Upper neighbor
                np.roll(x[:, :, i], -1, axis=0) +  # Lower neighbor
                np.roll(x[:, :, i],  1, axis=1) +  # Left neighbor
                np.roll(x[:, :, i], -1, axis=1) -  # Right neighbor
                4 * x[:, :, i]  # Center
            )
        
        # Exchange rates between nodes (patches)
        for i in range(len(RN.SpStr)):
            dxdt[:, :, i] += exchange_rates[:, :, i]
        
        return dxdt.flatten()
    # Verify that t_span is a tuple with two values
    if not (isinstance(t_span, tuple) and len(t_span) == 2):
        raise TypeError(f"t_span must be a tuple (t0, tf), but received {t_span} of type {type(t_span)}.")
    
    # Define the time range
    t = np.linspace(t_span[0], t_span[1], n_steps)
    
    # Solve the system of ODEs
    result = odeint(ode_system, x0.flatten(), t, args=(k, RN, exchange_rates, D, grid_shape))
    
    # Reconstruct results as a dictionary 
    modules = {
        f'Patch ({i+1}, {j+1})': result[:, (i * grid_shape[1] + j    ) * len(RN.SpStr):
                                           (i * grid_shape[1] + j + 1) * len(RN.SpStr)]
        for i in range(grid_shape[0]) for j in range(grid_shape[1])
    }
    # Number of patches
    num_patches = grid_shape[0] * grid_shape[1]
    print(f"Number of patches: {num_patches}")
    
    return modules

###########################################################################################
# Function to simulate the dynamics of a metapopulation with Michaelis-Menten kinetics and diffusion
def simulate_odes_metapop_mmk(RN, x0=None, t_span=None, n_steps=None, Vmax=None, Km=None, D=None, grid_shape=None):
    """
    Simulates the dynamics of a metapopulation with Michaelis-Menten kinetics and diffusion on a spatial grid.

    Parameters:
        RN (dict): Object describing the reaction network.
        x0 (array-like, optional): Initial conditions of the species.
        t_span (tuple, optional): Simulation time range.
        n_steps (int, optional): Number of simulation steps.
        Vmax (array-like, optional): Maximum reaction rate for each reaction.
        Km (array-like, optional): Michaelis-Menten constant for each reaction.
        D (array-like, optional): Diffusion coefficients for each species.
        grid_shape (tuple, optional): Shape of the spatial grid (rows, columns).

    Returns:
        dict: Simulation results for each patch (key: patch index).
    """
    if RN is None or not hasattr(RN, 'SpStr') or not hasattr(RN, 'RnStr'):
        raise ValueError("The reaction network (RN) is invalid or incomplete.")
    
    if grid_shape is None:
        grid_shape = (2, 2)
    if x0 is None:
        x0 = np.random.rand(grid_shape[0], grid_shape[1], len(RN.SpStr)) * 2
    if n_steps is None:
        n_steps = 500
    if t_span is None:
        t_span = (0, 100)
    if Vmax is None:
        Vmax = np.random.rand(len(RN.RnStr))
    if Km is None:
        Km = np.random.rand(len(RN.RnStr))
    if D is None:
        D = np.random.rand(len(RN.SpStr))
    
    def ode_system(x_flat, t, Vmax, Km, RN, D, grid_shape):
        x = x_flat.reshape(grid_shape[0], grid_shape[1], -1)
        dxdt = np.zeros_like(x)
        
        S = np.array(RN.RnMsupp)
        P = np.array(RN.RnMprod)
        
        for i in range(len(RN.SpStr)):
            for r in range(len(RN.RnStr)):
                rate = Vmax[r] * (x[:, :, :] ** S[r, :]) / (Km[r] + x[:, :, :])
                dxdt[:, :, i] += np.sum(rate * (P[r, i] - S[r, i]), axis=-1)
        
        for i in range(len(RN.SpStr)):
            dxdt[:, :, i] += D[i] * (
                np.roll(x[:, :, i],  1, axis=0) + 
                np.roll(x[:, :, i], -1, axis=0) + 
                np.roll(x[:, :, i],  1, axis=1) + 
                np.roll(x[:, :, i], -1, axis=1) - 
                4 * x[:, :, i]
            )
        
        return dxdt.flatten()
    
    t = np.linspace(t_span[0], t_span[1], n_steps)
    result = odeint(ode_system, x0.flatten(), t, args=(Vmax, Km, RN, D, grid_shape))
    
    modules = {
        f'Patch ({i+1}, {j+1})': result[:, (i * grid_shape[1] + j) * len(RN.SpStr):
                                           (i * grid_shape[1] + j + 1) * len(RN.SpStr)]
        for i in range(grid_shape[0]) for j in range(grid_shape[1])
    }
    
    print(f"Number of patches: {grid_shape[0] * grid_shape[1]}")
    return modules

###########################################################################################
# Function to create the simulations of the metapopulation dynamics of a reaction-diffusion model PDE
def simulate_pde_rd(RN, x0=None, t_span=None, n_steps=None, k=None, exchange_rates=None, D=None, grid_shape=None):
    if RN is None or not hasattr(RN, 'SpStr') or not hasattr(RN, 'RnStr'):
        raise ValueError("The reaction network (RN) is invalid or incomplete.")
    
    if grid_shape is None:
        grid_shape = (2, 2)
    if x0 is None:
        x0 = np.random.rand(grid_shape[0], grid_shape[1], len(RN.SpStr)) * 2
        print(f"Values of x0: {x0}")
    if n_steps is None:
        n_steps = 500
    if t_span is None:
        t_span = (0, 100)
    if k is None:
        k = np.random.rand(len(RN.RnStr)) * 1 # Random rate constants into the interval [0, 1]
    if exchange_rates is None:
        exchange_rates = np.ones((grid_shape[0], grid_shape[1], len(RN.SpStr)))  # No exchange by default        
    if D is None:
        D = np.random.rand(len(RN.SpStr)) # Default diffusion coefficients 
    
    def ode_system(x_flat, t, k, RN, D, grid_shape):
        x = x_flat.reshape(grid_shape[0], grid_shape[1], -1)
        dxdt = np.zeros_like(x)
        
        S = np.array(RN.RnMsupp)
        P = np.array(RN.RnMprod)
        
        for i in range(len(RN.SpStr)):
            for r in range(len(RN.RnStr)):
                rate = k[r] * np.prod(x[:, :, :] ** S[r, :], axis=-1) # Calculate the reaction rate 
                dxdt[:, :, i] += rate * (P[r, i] - S[r, i])           # Mass action kinetics
            
            # Difusión con diferencias finitas centradas en una malla 2D con coeficientes D_i específicos para cada especie.
            dxdt[:, :, i] += D[i] * (
                np.roll(x[:, :, i],  1, axis=0) +  # Vecino superior
                np.roll(x[:, :, i], -1, axis=0) +  # Vecino inferior
                np.roll(x[:, :, i],  1, axis=1) +  # Vecino derecho
                np.roll(x[:, :, i], -1, axis=1) -  # Vecino izquierdo
                4 * x[:, :, i]  # Centro (multiplicado por 4)
            )

        # Exchange rates between nodes (patches)
        for i in range(len(RN.SpStr)):
            dxdt[:, :, i] += exchange_rates[:, :, i]            
        
        return dxdt.flatten()
    
    if not (isinstance(t_span, tuple) and len(t_span) == 2):
        raise TypeError(f"t_span must be a tuple (t0, tf), but received {t_span} of type {type(t_span)}.")
    
    t = np.linspace(t_span[0], t_span[1], n_steps)
    
    result = odeint(ode_system, x0.flatten(), t, args=(k, RN, D, grid_shape))
    
    return result.reshape(n_steps, grid_shape[0], grid_shape[1], len(RN.SpStr))

###########################################################################################
# Function for solve the Linear Programming Problem: S.v>=0, v>0
def minimize_sv(S, epsilon=1,method='highs'):
    """
    Minimizes c @ v subject to S @ v >= 0 and v > 0.
    
    Parameters:
    ----------
    S : numpy.ndarray
        Stoichiometric matrix of the reaction network (RN).
    epsilon : float, optional
        Small positive value to ensure all coordinates of v remain strictly positive (default is 1).

    Returns:
    -------
    numpy.ndarray
        Optimal process vector v with all strictly positive coordinates.

    Raises:
    ------
    ValueError
        If no feasible solution is found.
    """
    n_species, n_reactions = S.shape  # Dimensions of S
    c = np.ones(n_reactions)          # Objective function: minimize the sum of v
    A_ub = -S                         # Reformulate S @ v >= 0 as -S @ v <= 0
    b_ub = np.zeros(n_species)        # Inequality constraints
    bounds = [(epsilon, None)] * n_reactions  # v > 0 (avoiding exact zero values) 
    
    # Solve the linear programming problem: minimize c @ v subject to (A_ub @ v <= b_ub)
    result = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=bounds, 
                     method=method) # 'highs' uses the Dual Simplex method (good for many constraints)
                                     # 'highs-ipm' uses the Interior Point method (faster for large sparse problems)
                    
    
    if result.success:
        return result.x
    else:
        raise ValueError("No feasible solution was found.")
    
# Function for solve the Linear Programming Problem: S.v>=0, v>0
def maximize_sv(S, epsilon=1, method='highs'):
    """
    Maximizes c @ v subject to S @ v >= 0 and v > 0.
    
    Parameters:
    ----------
    S : numpy.ndarray
        Stoichiometric matrix of the reaction network (RN).
    epsilon : float, optional
        Small positive value to ensure all coordinates of v remain strictly positive (default is 1).

    Returns:
    -------
    numpy.ndarray
        Optimal process vector v with all strictly positive coordinates.

    Raises:
    ------
    ValueError
        If no feasible solution is found.
    """
    n_species, n_reactions = S.shape  # Dimensions of S
    c = -np.ones(n_reactions)         # Objective function: maximize the sum of -v (equivalent to minimize c @ v)
    A_ub = -S                         # Reformulate S @ v >= 0 as -S @ v <= 0
    b_ub = np.zeros(n_species)        # Inequality constraints
    bounds = [( None,epsilon)] * n_reactions  # v > 0 (avoiding exact zero values) 
    
    # Solve the linear programming problem: maximize c @ v subject to (A_ub @ v <= b_ub)
    result = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method=method)
    
    if result.success:
        return result.x
    else:
        raise ValueError("No feasible solution was found.")

    
########################################################################################### 
# Function to calculate the echenlon form of a matrix
def row_echelon_form(A):
    A = A.astype(float)  # Asegurar que los cálculos sean en punto flotante
    rows, cols = A.shape
    for i in range(min(rows, cols)):
        # Encontrar el pivote y hacer intercambio de filas si es necesario
        max_row = np.argmax(abs(A[i:, i])) + i
        if A[max_row, i] != 0:
            A[[i, max_row]] = A[[max_row, i]]  # Intercambio de filas
        
        # Hacer ceros debajo del pivote
        for j in range(i + 1, rows):
            if A[i, i] != 0:
                factor = A[j, i] / A[i, i]
                A[j, i:] -= factor * A[i, i:]
    
    return A