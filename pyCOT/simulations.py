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
    if not reactants:  # Handle case where no reactants are defined
            return 0.0    
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
    if not substrate:  # Handle case where no substrate is defined
            return 0.0    
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

# def rate_cosine(substrates, concentrations, species_idx, spec_vector, t=None):
#     """
#     Cosine inflow kinetics: A*(1 + cos(w*t))/2 (always positive, oscillates between 0 and A)
#     Used for: R1 (seasonal water inflow)
#     Parameters: [A, w] where A is amplitude, w is frequency
#     """
#     A, w = spec_vector
    
#     # Handle time parameter - it should be passed from the modified simulation
#     if t is not None:
#         current_t = float(t)
#     else:
#         print("WARNING: No time parameter received, using t=0")
#         current_t = 0.0
    
#     # Use (1 + cos(wt))/2 to ensure always positive and smooth oscillation
#     rate = A * (1 + np.cos(w * current_t/(2*np.pi))) / 2
    
#     # Debug print
#     if int(current_t * 100) % 100 == 0:  # Print every 1 time unit
#         print(f"Time: {current_t:.2f}, Rate: {rate:.4f}")
    
#     return rate
# rate_cosine.expression = lambda substrates, reaction: (
#     f"A_{reaction} * (1 + cos(w_{reaction} * t)) / 2"
# )
# Function to calculate the rate of reaction of the type 'cosine'
def rate_cosine(reactants, concentrations, species_idx, spec_vector, t=0):
    """
    Calculate the rate of an inflow using a cosine-based oscillatory function.
    This kinetic law models an oscillating inflow (e.g., seasonal water inflow)
    that smoothly varies between 0 and amplitude A.
    
    Parameters:
    - reactants: List of tuples (species, stoichiometric coefficient) for the reactants
      (not used here but kept for consistency with the interface).
    - concentrations: List of current concentrations of all species (not used here).
    - species_idx: Dictionary mapping species names to their indices (not used here).
    - spec_vector: List of parameters [A, w], where
      A = amplitude (maximum inflow),
      w = frequency of oscillation.
    - t: Current simulation time (float). Default is 0.
    
    Returns:
    - rate: The instantaneous inflow rate at time t.
    """
    # Extract amplitude (A) and frequency (w)
    A, w = spec_vector
    # Cosine inflow formula (always positive)
    # Oscillates smoothly between 0 and A
    rate = A * (1 + np.cos(w * t / (2 * np.pi))) / 2
    return rate
rate_cosine.expression = lambda reactants, reaction: (
    f"A_{reaction} * (1 + cos(w_{reaction} * t)) / 2"
)

def rate_saturated(substrates, concentrations, species_idx, spec_vector):
    """
    Saturated kinetics: (Vmax * [substrate1] * [substrate2]) / (Km + [substrate1])
    Used for: R2, R3, R5, R6, R7, R9 (resource processing with saturation)
    Parameters: [Vmax, Km] where Vmax is max rate, Km is half-saturation
    """
    Vmax, Km = spec_vector
    
    if len(substrates) == 1:
        # Single substrate reaction
        substrate_name = substrates[0][0]
        substrate_conc = concentrations[species_idx[substrate_name]]
        return (Vmax * substrate_conc) / (Km + substrate_conc)
    
    elif len(substrates) >= 2:
        # Two substrate reaction - use first for saturation, multiply by both
        substrate1_name = substrates[0][0]  # Resource being processed
        substrate2_name = substrates[1][0]  # Population/catalyst
        
        substrate1_conc = concentrations[species_idx[substrate1_name]]
        substrate2_conc = concentrations[species_idx[substrate2_name]]
        
        return (Vmax * substrate1_conc * substrate2_conc) / (Km + substrate1_conc)
    
    else:
        return 0
rate_saturated.expression = lambda substrates, reaction: (
    f"(Vmax_{reaction} * [{substrates[0][0]}] * [{substrates[1][0] if len(substrates) > 1 else '1'}]) / (Km_{reaction} + [{substrates[0][0]}])"
    if len(substrates) >= 1 else "0"
)

def update_rate_laws(rate_list, additional_laws=None):
    # Diccionario base
    rate_laws = {
        'mak': rate_mak,
        'mmk': rate_mmk,
        'hill': rate_hill,
        'cosine': rate_cosine,
        'saturated': rate_saturated
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
def generate_random_vector(n, seed=None, min_value=0.1, max_value=1.0):
    rng = np.random.default_rng(seed)             # Generator with optional seed
    values = rng.uniform(min_value, max_value, n) # Random vector between [0.1, 1.0)
    return np.round(values, 2)                    # Round to 2 decimal places

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
def simulation(rn, rate='mak', spec_vector=None, x0=None, t_span=(0, 200),
               n_steps=200, additional_laws=None, method='LSODA', rtol=1e-8, atol=1e-10):
    
    species = [specie.name for specie in rn.species()]
    reactions = [reaction.name() for reaction in rn.reactions()]
    rn_dict = build_reaction_dict(rn)
    rate = validate_rate_list(rate, len(reactions))
    
    # Initial conditions vector by default random
    if x0 is None:
        x0 = generate_random_vector(len(species), min_value=50, max_value=100.0).astype(int).tolist()
        print(f"\nx0 = {x0}")
    
    # If spec_vector is None, generate random parameters for each reaction
    if spec_vector is None:
        spec_vector = []
        for kinetic in rate:
            if kinetic == 'mak':
                params = generate_random_vector(1, min_value=0.01, max_value=1.0)
            elif kinetic == 'mmk':
                Vmax = generate_random_vector(1, min_value=1, max_value=1.5)
                Km = generate_random_vector(1, min_value=5, max_value=10)
                params = np.round([float(Vmax), float(Km)], 2)
            elif kinetic == 'hill':
                Vmax = generate_random_vector(1, min_value=1, max_value=1.5)
                Kd = generate_random_vector(1, min_value=5, max_value=10)
                n = generate_random_vector(1, min_value=1, max_value=4) #np.random.randint(1, 4) #
                params = np.round([float(Vmax), float(Kd), float(n)], 2)
            elif kinetic in (additional_laws or {}):
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
    
    # Crear índice de especies
    species_idx = {s: i for i, s in enumerate(species)}
    
    # Función de velocidades CON RESTRICCIÓN DE NO-NEGATIVIDAD
    def rates_fn_constrained(t, x):
        # Asegurar que las concentraciones no sean negativas
        x_constrained = np.maximum(x, 0)
        dxdt = np.zeros_like(x_constrained)
        concentrations = x_constrained
        rates = np.zeros(len(reactions))
        
        for i, reaction in enumerate(reactions):
            kinetic = rate[i]
            reactants = rn_dict[reaction][0]
            
            if kinetic == 'mmk' and (not reactants or len(reactants) == 0):
                continue
            
            param_vector = parameters[reaction]
            rate_fn = rate_laws[kinetic]
            
            # *** KEY FIX: Pass time parameter for time-dependent kinetics ***
            if kinetic in ['cosine']:  # Add other time-dependent kinetics here
                rates[i] = rate_fn(reactants, concentrations, species_idx, param_vector, t=t)
            else:
                rates[i] = rate_fn(reactants, concentrations, species_idx, param_vector)
        
        for i, reaction in enumerate(reactions):
            for sp, coef in rn_dict[reaction][0]:
                dxdt[species_idx[sp]] -= coef * rates[i]
            for sp, coef in rn_dict[reaction][1]:
                dxdt[species_idx[sp]] += coef * rates[i]
        
        return dxdt
    
    # Imprimir ecuaciones y expresiones
    print_differential_equations(rn, rn_dict)
    print_velocity_expressions(rn, rn_dict, rate, additional_laws)
    
    # Puntos de tiempo
    time_points = np.linspace(t_span[0], t_span[1], n_steps)
    
    # INTEGRACIÓN CON solve_ivp
    try:
        sol = solve_ivp(rates_fn_constrained, t_span, x0, t_eval=time_points,
                       method=method, rtol=rtol, atol=atol)
        
        if not sol.success:
            print(f"Warning: Integration failed: {sol.message}")
            # Fallback a método más robusto
            print("Trying with LSODA method...")
            sol = solve_ivp(rates_fn_constrained, t_span, x0, t_eval=time_points,
                           method='LSODA', rtol=1e-6, atol=1e-8)
        
        # Extraer resultado y transponer para mantener compatibilidad
        result = sol.y.T
        
        # Asegurar que no hay valores negativos pequeños
        result = np.maximum(result, 0)
        
        # Reportar valores negativos corregidos
        negative_mask = sol.y.T < -1e-10  # Solo reportar negativos significativos
        if np.any(negative_mask):
            neg_count = np.sum(negative_mask)
            print(f"Warning: {neg_count} valores negativos detectados y corregidos a 0")
            
    except Exception as e:
        print(f"Error durante la integración: {e}")
        print("Usando método de respaldo (odeint)...")
        from scipy.integrate import odeint
        
        # Función original para odeint (t, x en orden diferente)
        def rates_fn_odeint(x, t):
            return rates_fn_constrained(t, x)
        
        result = odeint(rates_fn_odeint, x0, time_points)
        result = np.maximum(result, 0)  # Limpiar negativos
    
    # Crear DataFrames de salida
    time_series_df = time_series(rn, result, time_points)
    flux_vector_df = flux_vector(rn, rn_dict, rate, result, time_points, spec_vector, additional_laws)
    
    return time_series_df, flux_vector_df

###############################################################################
#################################################################################
#########################################################################
# # Function to create the stoichiometric matrix 
# def universal_stoichiometric_matrix(rn):
#     """
#     Creates the stoichiometric matrix based on the reaction network (RN).
    
#     Parameters:
#     RN (object): The reaction network object that contains species and reaction details.
    
#     Returns:
#     np.ndarray: The stoichiometric matrix representing the reactions, with integer values.
#     """
#     species =   rn.stoichiometry_matrix().species
#     reactions = rn.stoichiometry_matrix().reactions
#     reactants_vectors, products_vectors=build_stoichiometric_vectors(file_path, species)
#     # Initialize the stoichiometric matrix with zeros
#     matrix = np.zeros((len(species), len(reactions)))  # #species x #reactions

#     # Iterate over reactions
#     for i in range(len(reactions)):  # For each reaction 
#         for j in range(len(species)):  # Same for products
#             coef_p = products_vectors[i][j]  # Coefficient of the product in reaction i for species j
#             coef_r = reactants_vectors[i][j]  # Coefficient of the reactant in reaction i for species j
#             matrix[j, i] = coef_p - coef_r  # Product coefficients are positive (transposition here)  

#     return matrix  # Returns the stoichiometric matrix with integer values

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

####################################################################################
# Function to simulate diffusion dynamics 2D 
def simulate_diffusion_dynamics_2D(rn, rate='mak', grid_shape=None, D_dict=None, 
                                     x0_dict=None, spec_vector=None, t_span=(0, 20), n_steps=500, additional_laws=None):
    """
    Simulates the reaction-diffusion dynamics of a reaction network in a 2D grid.
    Parameters:
    - rn: Reaction network object containing species and reactions.
    - rate: List of kinetic laws for each reaction (default is 'mak').
    - grid_shape: Tuple defining the shape of the grid (rows, cols).
    - D_dict: Dictionary with diffusion coefficients for each species (default is random).
    - x0_dict: Dictionary with initial conditions for each species (default is random).
    - spec_vector: List of parameters for each reaction (default is random).
    - t_span: Tuple defining the time span for the simulation (default is (0, 20)).
    - n_steps: Number of time steps for the simulation (default is 500).
    - additional_laws: Dictionary with additional kinetic laws (default is None).
    Returns:
    - t: array of time points
    - X_out: dictionary with concentration time series for each species
    - flux_out: dictionary with flux time series for each reaction
    """
    np.random.seed(seed=42)  # Para reproducibilidad 
    species = [specie.name for specie in rn.species()]
    reactions = [reaction.name() for reaction in rn.reactions()] 

    rate = validate_rate_list(rate, len(reactions)) 

    if grid_shape is None:
        grid_shape = (2, 2)  # Por defecto, una cuadrícula de 2x2

    rows, cols = grid_shape
    num_cells = rows * cols

    if D_dict is None:
        D_dict = {sp: np.round(np.random.uniform(0.01, 0.2), 3) for sp in species}

    if x0_dict is None:
        x0_dict = {sp: np.round(np.random.uniform(0, 2.0, size=(rows, cols)), 2) for sp in species}

    # Flatten the initial state dictionary to a 1D array
    def flatten_state(x0_dict):
        return np.concatenate([x0_dict[sp].flatten() for sp in species])

    # Reshape the initial state to match the grid shape
    def reshape_state(x):
        return {sp: x[i*num_cells:(i+1)*num_cells].reshape((rows, cols)) for i, sp in enumerate(species)}

    x0 = flatten_state(x0_dict)

    # If spec_vector is None, generate random parameters for each reaction
    if spec_vector is None:
        spec_vector = []
        for kinetic in rate:
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
            elif kinetic in (additional_laws or {}):
                params = np.round(np.random.uniform(0.1, 1.0, 3), 3)
            else:
                raise ValueError(f"Unknown kinetic law: {kinetic}")
            spec_vector.append(params.tolist())

    ## Define the ODE system for the reaction-diffusion dynamics
    # reaction_dynamics computes the local reaction dynamics for each species
    def reaction_dynamics(Xdict, rate=rate, spec_vector=spec_vector):
        dxdt_dict = {sp: np.zeros((rows, cols)) for sp in species}
        flux_dict = {r: np.zeros((rows, cols)) for r in reactions}
        
        # Each point (i,j) of the grid is traversed.
        for i in range(rows):
            for j in range(cols):
                local_x = [Xdict[sp][i, j] for sp in species] # Current local concentrations
                
                # Calculate reaction rates and fluxes directly
                for r_idx, reaction in enumerate(rn.reactions()):
                    # Obtener reactantes, productos y estequiometría
                    reactants, products, stoichiometry = get_reaction_components(reaction, species)
                    
                    # Calcular tasa de reacción según la ley cinética
                    v_r = 0
                    if rate[r_idx] == 'mak':
                        k = spec_vector[r_idx][0]
                        v_r = k
                        if reactants:  # Reacciones con reactantes
                            for sp_idx, stoich in reactants:
                                v_r *= local_x[sp_idx] ** abs(stoich)
                    elif rate[r_idx] == 'mmk':
                        Vmax, Km = spec_vector[r_idx]
                        sp_idx = reactants[0][0] if reactants else 0
                        v_r = Vmax * local_x[sp_idx] / (Km + local_x[sp_idx]) if reactants else 0
                    elif rate[r_idx] == 'hill':
                        Vmax, Kd, n = spec_vector[r_idx]
                        sp_idx = reactants[0][0] if reactants else 0
                        v_r = Vmax * (local_x[sp_idx] ** n) / (Kd ** n + local_x[sp_idx] ** n) if reactants else 0
                    elif rate[r_idx] in (additional_laws or {}):
                        v_r = additional_laws[rate[r_idx]](local_x, spec_vector[r_idx])
                    
                    flux_dict[reaction.name()][i, j] = v_r
                    for sp_idx, stoich in enumerate(stoichiometry):
                        if stoich != 0:
                            dxdt_dict[species[sp_idx]][i, j] += stoich * v_r
        
        return dxdt_dict, flux_dict

    # The Laplacian is computed as the difference between the central value and the average of the neighbors
    def diffusion_term(Xdict):
        def laplacian(D):
            # Return a simple finite difference approximation of the diffusion term
            return (-4*D +
                    np.roll(D, 1, axis=0) + # Arriba
                    np.roll(D, -1, axis=0)+ # Abajo
                    np.roll(D, 1, axis=1) + # Derecha
                    np.roll(D, -1, axis=1)) # Izquierda
        return {sp: D_dict.get(sp, 0.0) * laplacian(Xdict[sp]) for sp in species}

    # Combined ODE function for reaction and diffusion dynamics
    def combined_ode(t, x):
        Xdict = reshape_state(x)
        dxdt_reac, _ = reaction_dynamics(Xdict, rate=rate, spec_vector=spec_vector)
        dxdt_diff = diffusion_term(Xdict)
        dxdt_total = {sp: dxdt_reac[sp] + dxdt_diff[sp] for sp in species}
        return flatten_state(dxdt_total)

    # ODE integration
    t_eval = np.linspace(*t_span, n_steps)
    sol = solve_ivp(combined_ode, t_span, x0, t_eval=t_eval, method='RK45', rtol=1e-6)

    # Formatear salida de concentraciones
    X_out = {sp: np.zeros((n_steps, rows, cols)) for sp in species}
    for i, xt in enumerate(sol.y.T):
        xt_dict = reshape_state(xt)
        for sp in species:
            X_out[sp][i] = xt_dict[sp]
    
    # Calcular series temporales de flujos
    flux_out = {r: np.zeros((n_steps, rows, cols)) for r in reactions}
    
    for time_idx, t_val in enumerate(sol.t):
        x_val = sol.y[:, time_idx]
        Xdict = reshape_state(x_val)
        _, flux_dict = reaction_dynamics(Xdict, rate, spec_vector)
        
        for r in reactions:
            flux_out[r][time_idx] = flux_dict[r]
    
    return sol.t, X_out, flux_out  

###################################################################################
# METAPOPULATION SIMULATION
###################################################################################

def get_reaction_components(reaction, species):
    """
    Obtiene reactantes, productos y estequiometría de un objeto Reaction de pyCOT.
    Construye la estequiometría manualmente usando support_edges() y products_edges().
    """
    try:
        # Reactantes desde support_edges()
        reactants = [(species.index(edge.species_name), edge.coefficient) for edge in reaction.support_edges()]
        # Productos desde products_edges()
        products = [(species.index(edge.species_name), edge.coefficient) for edge in reaction.products_edges()]
        # Construir estequiometría manualmente
        stoichiometry = [0] * len(species)
        for sp_idx, coeff in reactants:
            stoichiometry[sp_idx] -= coeff  # Reactantes tienen coeficientes negativos
        for sp_idx, coeff in products:
            stoichiometry[sp_idx] += coeff  # Productos tienen coeficientes positivos
        return reactants, products, stoichiometry
    except Exception as e:
        print(f"Error al obtener componentes de la reacción {reaction.name()}: {e}")
        print(f"Atributos disponibles: {dir(reaction)}")
        raise
def simulate_metapopulation_dynamics(rn, rate='mak', num_patches=None, D_dict=None, 
                                    x0_dict=None, spec_vector=None, t_span=(0, 20), 
                                    n_steps=500, connectivity_matrix=None, additional_laws=None):
    """
    Simula la dinámica de metapoblaciones con reacciones locales y dispersión entre parches.
    Implementación vectorial: las poblaciones se tratan como elementos de un vector.
    
    Parámetros:
    - rn: Objeto ReactionNetwork con especies y reacciones.
    - rate: Lista de leyes cinéticas para cada reacción ('mak', 'mmk', 'hill' o adicionales).
    - num_patches: Número de parches/poblaciones (e.g., 3 para tres poblaciones).
    - D_dict: Diccionario con tasas de dispersión global por especie (default: aleatorio).
    - x0_dict: Diccionario con condiciones iniciales como vectores 1D de longitud num_patches
               para cada especie (e.g., {"species1": np.array([0.5, 1.0, 0.3])}).
    - spec_vector: Lista de parámetros para cada reacción (default: aleatorio).
    - t_span: Tupla con el intervalo de tiempo de simulación (default: (0, 20)).
    - n_steps: Número de pasos temporales (default: 500).
    - connectivity_matrix: Matriz de conectividad (num_patches × num_patches).
                           Las filas deben sumar 1 (default: aleatorio).
    - additional_laws: Diccionario con leyes cinéticas adicionales (default: None).
    
    Retorna:
    - t: Array de tiempos de forma (n_steps,)
    - X_out: Diccionario con arrays de concentraciones de forma (n_steps, num_patches) para cada especie
    - flux_out: Diccionario con arrays de flujos de forma (n_steps, num_patches) para cada reacción
    """
    np.random.seed(seed=42)  # Para reproducibilidad
    species = [specie.name for specie in rn.species()]
    reactions = [reaction.name() for reaction in rn.reactions()]
    
    rate = validate_rate_list(rate, len(reactions))
    
    # Configurar número de parches
    if num_patches is None:
        num_patches = 3  # Default: 3 poblaciones
    
    # Tasas de dispersión por especie
    if D_dict is None:
        D_dict = {sp: np.round(np.random.uniform(0.01, 0.2), 3) for sp in species}
    
    # Condiciones iniciales (vectores 1D)
    if x0_dict is None:
        x0_dict = {sp: np.round(np.random.uniform(0, 2.0, size=num_patches), 2) for sp in species}
    else:
        # Validar que x0_dict contiene vectores de la longitud correcta
        for sp in species:
            if sp not in x0_dict:
                raise ValueError(f"Falta condición inicial para la especie '{sp}'")
            if len(x0_dict[sp]) != num_patches:
                raise ValueError(f"La condición inicial para '{sp}' debe tener longitud {num_patches}")
    
    # Matriz de conectividad o probabilidad de dispersión
    if connectivity_matrix is None:
        # Generar matriz aleatoria con dos decimales
        connectivity_matrix = np.random.uniform(0, 1, size=(num_patches, num_patches))
        connectivity_matrix = np.round(connectivity_matrix, 2)
        
        # Asegurar que las filas sumen 1
        row_sums = np.sum(connectivity_matrix, axis=1, keepdims=True)
        connectivity_matrix = connectivity_matrix / row_sums
        connectivity_matrix = np.round(connectivity_matrix, 2)
        
        # Ajustar para compensar errores de redondeo
        for i in range(num_patches):
            current_sum = np.sum(connectivity_matrix[i])
            if current_sum != 1.0:
                diff = 1.0 - current_sum
                max_idx = np.argmax(connectivity_matrix[i])
                connectivity_matrix[i, max_idx] += diff
                connectivity_matrix[i, max_idx] = np.round(connectivity_matrix[i, max_idx], 2)
        
        
    else:
        # Validar connectivity_matrix
        if connectivity_matrix.shape != (num_patches, num_patches):
            raise ValueError(f"La matriz de conectividad debe tener forma ({num_patches}, {num_patches})")
        if not np.allclose(np.sum(connectivity_matrix, axis=1), 1.0, atol=1e-6):
            print("Advertencia: Normalizando la matriz de conectividad para que las filas sumen 1")
            connectivity_matrix = connectivity_matrix / np.sum(connectivity_matrix, axis=1, keepdims=True)
            connectivity_matrix = np.round(connectivity_matrix, 2)
    
    # Imprimir matriz final
    print("\nMatriz de conectividad o probabilidad:\n", connectivity_matrix)
    print(f"Sumas por fila: {np.sum(connectivity_matrix, axis=1)}")
    
    # Parámetros de reacción
    if spec_vector is None:
        spec_vector = []
        for kinetic in rate:
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
            elif kinetic in (additional_laws or {}):
                params = np.round(np.random.uniform(0.1, 1.0, 3), 3)
            else:
                raise ValueError(f"Ley cinética desconocida: {kinetic}")
            spec_vector.append(params.tolist())
    
    # Aplanar el estado inicial (concatenar vectores de todas las especies)
    def flatten_state(x_dict):
        return np.concatenate([x_dict[sp] for sp in species])
    
    # Reconstruir el estado (deshacer el flatten)
    def reshape_state(x):
        return {sp: x[i*num_patches:(i+1)*num_patches] for i, sp in enumerate(species)}
    
    x0 = flatten_state(x0_dict)
    
    # Término de reacciones locales
    def reaction_dynamics(Xdict, rate, spec_vector):
        dxdt_dict = {sp: np.zeros(num_patches) for sp in species}
        flux_dict = {r: np.zeros(num_patches) for r in reactions}
        
        # Para cada parche
        for patch_idx in range(num_patches):
            local_x = [Xdict[sp][patch_idx] for sp in species]
            
            # Para cada reacción
            for r_idx, reaction in enumerate(rn.reactions()):
                # Obtener reactantes, productos y estequiometría
                reactants, products, stoichiometry = get_reaction_components(reaction, species)
                
                # Calcular tasa de reacción según la ley cinética
                v_r = 0
                if rate[r_idx] == 'mak':
                    k = spec_vector[r_idx][0]
                    v_r = k
                    if reactants:  # Reacciones con reactantes
                        for sp_idx, stoich in reactants:
                            v_r *= local_x[sp_idx] ** abs(stoich)
                elif rate[r_idx] == 'mmk':
                    Vmax, Km = spec_vector[r_idx]
                    sp_idx = reactants[0][0] if reactants else 0
                    v_r = Vmax * local_x[sp_idx] / (Km + local_x[sp_idx]) if reactants else 0
                elif rate[r_idx] == 'hill':
                    Vmax, Kd, n = spec_vector[r_idx]
                    sp_idx = reactants[0][0] if reactants else 0
                    v_r = Vmax * (local_x[sp_idx] ** n) / (Kd ** n + local_x[sp_idx] ** n) if reactants else 0
                elif rate[r_idx] in (additional_laws or {}):
                    v_r = additional_laws[rate[r_idx]](local_x, spec_vector[r_idx])
                
                flux_dict[reaction.name()][patch_idx] = v_r
                for sp_idx, stoich in enumerate(stoichiometry):
                    if stoich != 0:
                        dxdt_dict[species[sp_idx]][patch_idx] += stoich * v_r
        
        return dxdt_dict, flux_dict
    
    # Término de dispersión
    def dispersal_term(Xdict):
        dxdt_dict = {sp: np.zeros(num_patches) for sp in species}
        
        for sp in species:
            D = D_dict.get(sp, 0.0)  # Tasa de dispersión para la especie
            X = Xdict[sp]  # Vector de concentraciones de la especie
            
            for i in range(num_patches):
                dxdt = 0
                for j in range(num_patches):
                    # Flujo entrante desde parche j hacia parche i
                    dxdt += D * connectivity_matrix[j, i] * X[j]
                    # Flujo saliente desde parche i hacia parche j
                    dxdt -= D * connectivity_matrix[i, j] * X[i]
                dxdt_dict[sp][i] = dxdt
        
        return dxdt_dict
    
    # Sistema ODE combinado
    def combined_ode(t, x):
        x = np.maximum(x, 0)
        Xdict = reshape_state(x)
        dxdt_reac, _ = reaction_dynamics(Xdict, rate, spec_vector)
        dxdt_disperse = dispersal_term(Xdict)
        dxdt_total = {sp: dxdt_reac[sp] + dxdt_disperse[sp] for sp in species}
        return flatten_state(dxdt_total)
    
    # Integración de ODEs
    t_eval = np.linspace(t_span[0], t_span[1], n_steps)
    sol = solve_ivp(combined_ode, t_span, x0, t_eval=t_eval, method='LSODA', rtol=1e-8, atol=1e-10)
    
    # Formatear salida
    X_out = {sp: np.zeros((n_steps, num_patches)) for sp in species}
    for i, xt in enumerate(sol.y.T):
        xt_dict = reshape_state(xt)
        for sp in species:
            X_out[sp][i] = xt_dict[sp]
    
    # Calcular series temporales de flujos
    flux_out = {r: np.zeros((n_steps, num_patches)) for r in reactions}
    
    for time_idx, t_val in enumerate(sol.t):
        x_val = sol.y[:, time_idx]
        Xdict = reshape_state(x_val)
        _, flux_dict = reaction_dynamics(Xdict, rate, spec_vector)
        
        for r in reactions:
            flux_out[r][time_idx] = flux_dict[r]
    
    return sol.t, X_out, flux_out

###################################################################################

###################################################################################

########################################################################################### 
# LINEAR PROGRAMMING 
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