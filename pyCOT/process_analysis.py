"""
Process Analysis Library for Reaction Networks

This module provides analysis functions for classifying processes, analyzing cones,
and processing time series of reaction network dynamics.

Main components:
1. Process Classification - classify process vectors based on stoichiometric effects
2. Time Series Processing - rescaling and batch classification
3. Cone Analysis - nullspace, feasible regions, and comprehensive cone analysis
4. Utilities - helper functions for interval extraction and data export
"""

import numpy as np
import pandas as pd
import sympy as sp
from typing import List, Dict, Tuple, Optional, Union
from collections import Counter
import os


######################################################################################
# 1. PROCESS CLASSIFICATION
######################################################################################

def classify_process_mode(v, S, tol=1e-3):
    """
    Clasifica un vector de proceso v con respecto a la matriz estequiométrica S.

    Categorías base:
    - Stationary Mode: El proceso no produce cambios netos en las concentraciones (Sv=0).
    - Problem: El proceso consume especies pero no produce ninguna (Sv <= 0 y al menos un Sv < 0).
    - Challenge: El proceso consume al menos una especie (al menos un Sv < 0).
    - Cognitive Domain: El proceso mantiene o aumenta todas las especies (Sv >= 0) y utiliza todas las reacciones de la red (v > 0).
    
    Categorías extendidas (requieren un proceso previo 'v_pert'):
    - Counteraction: Un proceso 'v' que, combinado con un 'Challenge' previo 'v_pert', resulta en un no-consumo neto (S(v + v_pert) >= 0).
    - Solution: Un proceso 'v' del espacio 'Challenge' que sirve como solución para un 'Problem' previo 'v_pert' (S(v+v_pert) >= 0).
    - Cognitive Control: Un proceso 'v' del 'Cognitive Domain' que sirve como solución para un 'Problem' previo 'v_pert' (S(v+v_pert) >= 0, v > 0).

    Args:
        v (np.ndarray): El vector de proceso a clasificar. Debe ser un array de NumPy que representa los flujos.
        S (np.ndarray): La matriz estequiométrica, donde las filas son especies y las columnas son reacciones.
        tol (float): Tolerancia para comparaciones numéricas con cero para manejar la imprecisión de punto flotante. Por defecto es 1e-3.
        
    Returns:
        list[str]: Una lista ordenada de las categorías a las que pertenece el proceso v.
                   Si no coincide con ninguna categoría específica, se clasifica como "Other".
                   Un proceso puede pertenecer a múltiples categorías.
    
    Raises:
        TypeError: Si 'v' o 'S' no son arrays de NumPy.
    """
    
    # Validaciones de tipo para los inputs
    if not isinstance(v, np.ndarray) or not isinstance(S, np.ndarray):
        raise TypeError("Los inputs 'v' y 'S' deben ser arrays de NumPy.")

    # Calcular el cambio neto en las concentraciones de las especies (Sv)
    Sv = S @ v
    classifications = []

    # --- Propiedades fundamentales de Sv y v (calculadas una vez para eficiencia) ---
    is_stationary = np.all((-tol <= Sv) & (Sv <= tol))
    is_overproduced = np.all(Sv >= -tol) and np.any(Sv > tol)
    is_challenge = np.any(Sv < -tol) and np.any(Sv > tol)
    is_problem = np.all(Sv <= tol) and np.any(Sv <= -tol)
    is_complete = np.all(v > tol) 
    
    # --- Clasificaciones de Modo y Completitud (basadas solo en 'v' y 'S') --- 
    if is_stationary:
        classifications = ["Stationary Mode"]
        
    elif is_overproduced:
        classifications = ["Overproduction Mode"]
    elif is_challenge:
        classifications = ["Challenge"]
    elif is_problem:
        classifications = ["Problem"]
    else:
        classifications = ["Other"]
    
    if is_complete:
        classifications.append("Complete Process")
    else:
        classifications.append("Incomplete Process")
    
    # Special case: Cognitive Control = Overproduction + Complete
    if is_overproduced and is_complete:
        # Replace "Overproduction Mode" with "Cognitive Control"
        classifications[0] = "Cognitive Control"
    
    return classifications


def is_cognitive_domain(v, S, tol=1e-3):
    """
    Verifica si un proceso v está en el dominio cognitivo.
    
    Un proceso está en el dominio cognitivo si:
    - Es Stationary Mode o Overproduction Mode
    - Es un proceso completo (v > 0 para todas las reacciones)
    
    Args:
        v (np.ndarray): El vector de proceso a verificar.
        S (np.ndarray): La matriz estequiométrica.
        tol (float): Tolerancia para comparaciones numéricas.
        
    Returns:
        bool: True si el proceso está en el dominio cognitivo, False en caso contrario.
    """
    v_class = classify_process_mode(v, S, tol=tol) 
    if (v_class[0] == "Stationary Mode" or v_class[0] == "Overproduction Mode" or v_class[0] == "Cognitive Control") and v_class[1] == "Complete Process":
        return True 
    else:
        return False


def classify_response_to_disturbance(v, S, v_pert=None, x_pert=None, tol=1e-3, verbose=False):
    """
    Clasifica la respuesta de un proceso a una perturbación.
    
    Analiza cómo un proceso v responde a una perturbación, que puede ser:
    - v_pert: Una perturbación en el espacio de procesos
    - x_pert: Una perturbación en el espacio de estados
    
    Args:
        v (np.ndarray): El vector de proceso de respuesta.
        S (np.ndarray): La matriz estequiométrica.
        v_pert (np.ndarray, optional): Perturbación en el espacio de procesos.
        x_pert (np.ndarray, optional): Perturbación en el espacio de estados.
        tol (float): Tolerancia para comparaciones numéricas.
        verbose (bool): Si True, imprime información de debug.
        
    Returns:
        list[str]: Clasificación de la respuesta a la perturbación.
        
    Raises:
        TypeError: Si los tipos de entrada no son correctos.
        ValueError: Si las dimensiones no coinciden.
    """
    if v_pert is not None and not isinstance(v_pert, np.ndarray):
        raise TypeError("El input 'v_pert', si se proporciona, debe ser un array de NumPy.")
    
    if v_pert is not None and v_pert.shape[0] != S.shape[1]:
        raise ValueError("El input 'v_pert' debe tener la misma dimensión que el número de reacciones en 'S'.")
    
    if x_pert is not None and not isinstance(x_pert, np.ndarray):
        raise TypeError("El input 'x_pert', si se proporciona, debe ser un array de NumPy.")
    
    if x_pert is not None and x_pert.shape[0] != S.shape[0]:
        raise ValueError("El input 'x_pert' debe tener la misma dimensión que el número de especies en 'S'.")

    v_mode = classify_process_mode(v, S, tol=tol)
    
    if v_pert is None and x_pert is None:
        if verbose:
            print("Null disturbance")
        return v_mode

    elif v_pert is not None and x_pert is None:
        if verbose:
            print("Process disturbance")
        v_pert_class = classify_process_mode(v_pert, S, tol=tol)
        v_cognitive_domain = is_cognitive_domain(v, S, tol=tol)
        v_combined = (v + v_pert)
        v_combined_class = classify_process_mode(v_combined, S, tol=tol)
        v_combined_cognitive_domain = is_cognitive_domain(v_combined, S, tol=tol)
        
        # Cognitive Control Situations
        if v_cognitive_domain:
            if verbose:
                print("In Cognitive Domain")
            # Cognitive Domain controla Challenge
            if v_pert_class[0] == "Challenge":
                if v_combined_cognitive_domain:
                    return ["Cognitive Controls Challenge"]
                else:
                    return ["Cognitive Breakdown by Challenge"]
            # Cognitive Domain controla Problem
            if v_pert_class[0] == "Problem": 
                if v_combined_cognitive_domain:
                    return ["Cognitive Controls Problem"]
                else:
                    return ["Cognitive Breakdown by Problem"]
            # Cognitive Domain se sostiene por proceso automantenido o sobreproducido
            if (v_pert_class[0] == "Stationary Mode" or v_pert_class[0] == "Overproduction Mode"): 
                if v_combined_cognitive_domain:
                    return ["Cognitive Domain Sustained"]
            else:
                if verbose:
                    print("Check Cognitive Breakdown by Glitch!!" + str(v_mode) + " + " + str(v_pert_class) + " =>" + str(v_combined_class))
                    return ["Cognitive Breakdown by Glitch"]
        else:
            if verbose:
                print("Out of Cognitive Domain")
                print(f"v_mode: {v_mode}, v_pert_class: {v_pert_class}, v_combined_class: {v_combined_class}")
                print(f"v_combined_cognitive_domain: {v_combined_cognitive_domain}")
            if v_combined_cognitive_domain:
                return ["Cognitive Domain Recovered"]
            else:
                return ["Cognitive Breakdown Sustained"]
                
    elif v_pert is None and x_pert is not None:
        # State Disturbance and response result
        v_cognitive_domain = is_cognitive_domain(v, S, tol=tol)
        x_next = x_pert + S @ v
        
        if v_cognitive_domain:
            if np.any(x_pert < -tol) and np.any(x_pert > tol):
                if verbose:
                    print("State disturbance is a challenge")
                if np.all(x_next >= -tol):
                    return ["Cognitive Controls Challenge"]
                else:
                    if np.any(x_next < -tol):
                        return ["Cognitive Breakdown by Challenge"]
            elif np.all(x_pert <= -tol) and np.any(x_pert < tol):
                if verbose:
                    print("State disturbance is a problem")
                if np.all(x_next >= -tol):
                    return ["Cognitive Controls Problem"]
                else:
                    if np.any(x_next < -tol):
                        return ["Cognitive Breakdown by Problem"]
            else:
                if verbose:
                    print("State disturbance is resources incoming")
                if np.all(x_next >= -tol):
                    return ["Cognitive Domain Sustained"]
                else:
                    return ["Cognitive Breakdown by Glitch"]
        else:
            if verbose:
                print("Out of Cognitive Domain")
                print(f"v_mode: {v_cognitive_domain} x_next: {x_next}")
            if np.all(x_next >= -tol):
                return ["Cognitive Domain Recovered"]
            else:
                return ["Cognitive Breakdown Sustained"]


######################################################################################
# 2. TIME SERIES PROCESSING
######################################################################################

def rescale_process_time_series(process_series, window_size=1):
    """
    Rescales a process time series by aggregating consecutive vectors using a rolling sum.
    
    This function performs a sliding window sum over process vectors, creating a new 
    time series where each point represents the sum of 'window_size' consecutive processes.
    
    Parameters:
    -----------
    process_series : pd.DataFrame
        DataFrame containing:
        - 'Time' column: timestamps for each process
        - Process columns: vector components (e.g., 'Flux_r1', 'Flux_r2', ...)
    
    window_size : int, default=1
        Number of consecutive processes to sum together.
        - window_size=1: Returns copy of original (no rescaling)
        - window_size>1: Applies rolling sum
    
    Returns:
    --------
    pd.DataFrame
        Rescaled time series with:
        - Length: len(process_series) - window_size + 1
        - 'Time' column: adjusted to last time in each window
        - Process columns: rolling sum of original vectors
    
    Examples:
    ---------
    >>> # Original time series (5 steps)
    >>> process_series = pd.DataFrame({
    ...     'Time': [0, 1, 2, 3, 4],
    ...     'Flux_r1': [1, 2, 3, 4, 5],
    ...     'Flux_r2': [0.5, 1.0, 1.5, 2.0, 2.5]
    ... })
    >>> 
    >>> # Rescale with window_size=3
    >>> rescaled = rescale_process_time_series(process_series, window_size=3)
    >>> # Result has 3 rows: times [2, 3, 4]
    >>> # Flux_r1: [6, 9, 12] = [1+2+3, 2+3+4, 3+4+5]
    """
    if not isinstance(process_series, pd.DataFrame):
        raise TypeError("process_series debe ser un DataFrame de pandas.")
    
    if 'Time' not in process_series.columns:
        raise ValueError("El DataFrame debe contener una columna 'Time'.")
    
    if window_size == 1:
        # No rescaling needed
        return process_series.copy()
    
    # Separate time and process data
    time_col = process_series['Time']
    process_numeric = process_series.drop(columns=['Time']).copy()
    
    # Apply rolling sum
    process_rolling = (
        process_numeric
        .rolling(window=window_size, min_periods=window_size)
        .sum()
        .dropna()
        .reset_index(drop=True)
    )
    
    # Adjust time to the last time in each window
    time_adjusted = time_col.iloc[window_size - 1:].reset_index(drop=True)
    
    # Combine time and rescaled processes
    rescaled_series = pd.concat([time_adjusted, process_rolling], axis=1)
    
    return rescaled_series


def classify_process_series(process_series, S, tol=1e-3):
    """
    Clasifica todos los procesos en una serie temporal.
    
    Esta función aplica classify_process_mode a cada vector de proceso en la serie,
    proporcionando una clasificación batch eficiente.
    
    Parameters:
    -----------
    process_series : pd.DataFrame
        DataFrame con columna 'Time' y columnas de proceso.
    S : np.ndarray
        Matriz estequiométrica.
    tol : float, default=1e-3
        Tolerancia para clasificación.
        
    Returns:
    --------
    list[list[str]]
        Lista de clasificaciones, una por cada paso temporal.
    """
    if not isinstance(process_series, pd.DataFrame):
        raise TypeError("process_series debe ser un DataFrame de pandas.")
    
    if not isinstance(S, np.ndarray):
        raise TypeError("S debe ser un array de NumPy.")
    
    # Extract process vectors (exclude Time column)
    if 'Time' in process_series.columns:
        process_values = process_series.drop(columns=['Time'])
    else:
        process_values = process_series
    
    # Classify each process
    classifications = []
    for _, row in process_values.iterrows():
        v = row.to_numpy()
        cat = classify_process_mode(v, S, tol=tol)
        classifications.append(cat)
    
    return classifications


######################################################################################
# 3. CONE ANALYSIS
######################################################################################

def compute_nullspace_vectors(S):
    """
    Calcula los vectores del espacio nulo de la matriz estequiométrica S.
    
    El espacio nulo consiste en todos los vectores v tales que Sv = 0,
    representando modos estacionarios del sistema.
    
    Parameters:
    -----------
    S : np.ndarray
        Matriz estequiométrica (especies × reacciones).
        
    Returns:
    --------
    list[np.ndarray]
        Lista de vectores base del espacio nulo.
    """
    if not isinstance(S, np.ndarray):
        raise TypeError("S debe ser un array de NumPy.")
    
    S_sym = sp.Matrix(S)
    null_basis = S_sym.nullspace()
    null_vectors = [np.array(v, dtype=float).flatten() for v in null_basis]
    
    return null_vectors


def compute_feasible_region(S, grid_max=None, grid_res=5, auto_scale=True):
    """
    Calcula la región factible donde Sv >= 0 (región de no-consumo).
    
    Genera una cuadrícula de puntos en el espacio de procesos y filtra
    aquellos que satisfacen Sv >= 0 para todas las especies.
    
    Parameters:
    -----------
    S : np.ndarray
        Matriz estequiométrica (especies × reacciones).
    grid_max : float, optional
        Valor máximo para la cuadrícula. Si None, se calcula automáticamente.
    grid_res : int, default=5
        Resolución de la cuadrícula (puntos por dimensión).
    auto_scale : bool, default=True
        Si True, ajusta grid_max basándose en los vectores del espacio nulo.
        
    Returns:
    --------
    np.ndarray
        Array de puntos factibles (shape: n_points × n_reactions).
    """
    if not isinstance(S, np.ndarray):
        raise TypeError("S debe ser un array de NumPy.")
    
    n = S.shape[1]  # Number of reactions
    
    # Auto-scale grid_max if needed
    if grid_max is None or auto_scale:
        null_vectors = compute_nullspace_vectors(S)
        if null_vectors:
            base_max = np.max([np.max(np.abs(v)) for v in null_vectors])
        else:
            base_max = 1.0
        
        if grid_max is None:
            grid_max = base_max
        else:
            grid_max = max(grid_max, base_max)
    
    # Generate grid
    grid = np.linspace(0, grid_max, grid_res)
    V = np.array(np.meshgrid(*([grid] * n))).T.reshape(-1, n)
    
    # Filter feasible points (Sv >= 0)
    mask = np.all(S @ V.T >= -1e-10, axis=0)
    feasible_points = V[mask]
    
    return feasible_points


def classify_feasible_points(points, S, tol=1e-3):
    """
    Clasifica cada punto en la región factible.
    
    Parameters:
    -----------
    points : np.ndarray
        Array de vectores de proceso (shape: n_points × n_reactions).
    S : np.ndarray
        Matriz estequiométrica.
    tol : float, default=1e-3
        Tolerancia para clasificación.
        
    Returns:
    --------
    list[str]
        Lista de clasificaciones (una por punto).
    """
    if not isinstance(points, np.ndarray) or not isinstance(S, np.ndarray):
        raise TypeError("points y S deben ser arrays de NumPy.")
    
    classifications = []
    for v in points:
        cat = classify_process_mode(v, S, tol=tol)
        cat_str = ",".join(cat) if cat else "None"
        classifications.append(cat_str)
    
    return classifications


def analyze_cone(S, grid_max=None, grid_res=5, classify=True, tol=1e-3):
    """
    Análisis completo del cono definido por la matriz estequiométrica S.
    
    Combina el cálculo del espacio nulo, la región factible, y opcionalmente
    la clasificación de procesos en una sola función conveniente.
    
    Parameters:
    -----------
    S : np.ndarray
        Matriz estequiométrica (especies × reacciones).
    grid_max : float, optional
        Valor máximo para la cuadrícula.
    grid_res : int, default=5
        Resolución de la cuadrícula.
    classify : bool, default=True
        Si True, clasifica los puntos factibles.
    tol : float, default=1e-3
        Tolerancia para clasificación.
        
    Returns:
    --------
    dict
        Diccionario con claves:
        - 'nullspace_vectors': List[np.ndarray] - vectores del espacio nulo
        - 'feasible_points': np.ndarray - puntos factibles
        - 'Sv_values': np.ndarray - valores de S @ v para cada punto
        - 'grid_max': float - valor máximo usado
        - 'classifications': List[str] - clasificaciones (si classify=True)
        - 'classification_counts': dict - conteo por categoría (si classify=True)
    """
    if not isinstance(S, np.ndarray):
        raise TypeError("S debe ser un array de NumPy.")
    
    # Compute nullspace
    nullspace_vectors = compute_nullspace_vectors(S)
    
    # Compute feasible region
    feasible_points = compute_feasible_region(S, grid_max=grid_max, grid_res=grid_res, auto_scale=True)
    
    # Compute Sv for each point
    if feasible_points.shape[0] > 0:
        Sv_values = (S @ feasible_points.T).T
    else:
        Sv_values = np.array([])
    
    # Determine actual grid_max used
    if grid_max is None:
        if nullspace_vectors:
            grid_max = np.max([np.max(np.abs(v)) for v in nullspace_vectors])
        else:
            grid_max = 1.0
    
    # Build result dictionary
    result = {
        'nullspace_vectors': nullspace_vectors,
        'feasible_points': feasible_points,
        'Sv_values': Sv_values,
        'grid_max': grid_max
    }
    
    # Classify if requested
    if classify and feasible_points.shape[0] > 0:
        classifications = classify_feasible_points(feasible_points, S, tol=tol)
        result['classifications'] = classifications
        
        # Count classifications
        classification_counts = Counter(classifications)
        result['classification_counts'] = dict(classification_counts)
    
    return result


######################################################################################
# 4. UTILITIES
######################################################################################

def get_intervals_by_category(times, classifications, category):
    """
    Extrae intervalos temporales donde los procesos pertenecen a una categoría específica.
    
    Esta función identifica períodos continuos donde los procesos están clasificados
    en la categoría dada, útil para visualización y análisis de dinámicas.
    
    Parameters:
    -----------
    times : np.ndarray or pd.Series
        Timestamps para cada proceso.
    classifications : list[list[str]] or list[str]
        Clasificaciones para cada proceso. Puede ser:
        - List[List[str]]: clasificaciones múltiples por proceso
        - List[str]: clasificación única (string con categorías separadas por comas)
    category : str
        Categoría a buscar (e.g., "Cognitive Control", "Problem").
        
    Returns:
    --------
    list[tuple]
        Lista de intervalos (start_time, end_time) donde se encuentra la categoría.
        
    Examples:
    ---------
    >>> times = np.array([0, 1, 2, 3, 4, 5])
    >>> classifications = [["Problem"], ["Problem"], ["Cognitive Control"], 
    ...                    ["Cognitive Control"], ["Problem"], ["Problem"]]
    >>> intervals = get_intervals_by_category(times, classifications, "Cognitive Control")
    >>> # Returns: [(2, 4)]
    """
    if isinstance(times, pd.Series):
        times = times.to_numpy()
    
    # Create mask for category presence
    mask = np.array([category in cat if isinstance(cat, list) else category in cat 
                     for cat in classifications])
    
    # Extract intervals
    intervals = []
    in_interval = False
    start = None
    
    for t, flag in zip(times, mask):
        if flag and not in_interval:
            in_interval = True
            start = t
        elif not flag and in_interval:
            in_interval = False
            intervals.append((start, t))
    
    # Close last interval if still open
    if in_interval:
        intervals.append((start, times[-1]))
    
    return intervals


def export_classified_processes(process_data, S, filepath, 
                                format='excel', 
                                include_Sv=True,
                                separate_sheets=True):
    """
    Exporta procesos clasificados a Excel o CSV.
    
    Parameters:
    -----------
    process_data : pd.DataFrame
        DataFrame con vectores de proceso (puede incluir columna 'Time').
    S : np.ndarray
        Matriz estequiométrica para calcular Sv.
    filepath : str
        Ruta del archivo de salida.
    format : str, default='excel'
        Formato de salida: 'excel' o 'csv'.
    include_Sv : bool, default=True
        Si True, incluye columnas S*v.
    separate_sheets : bool, default=True
        Si True (solo Excel), crea hojas separadas por categoría.
        
    Returns:
    --------
    str
        Ruta del archivo guardado.
    """
    if not isinstance(process_data, pd.DataFrame):
        raise TypeError("process_data debe ser un DataFrame.")
    
    if not isinstance(S, np.ndarray):
        raise TypeError("S debe ser un array de NumPy.")
    
    # Extract process vectors
    if 'Time' in process_data.columns:
        process_values = process_data.drop(columns=['Time'])
    else:
        process_values = process_data
    
    # Classify processes
    process_types = []
    for _, row in process_values.iterrows():
        v = row.to_numpy()
        cat = classify_process_mode(v, S)
        cat_str = ",".join(cat) if cat else "None"
        process_types.append(cat_str)
    
    # Create output DataFrame
    output_df = process_data.copy()
    
    # Add Sv columns if requested
    if include_Sv:
        Sv_matrix = process_values.apply(lambda v: S @ v.to_numpy(), axis=1)
        Sv_expanded = pd.DataFrame(Sv_matrix.tolist(),
                                   columns=[f"S*v_{i+1}" for i in range(S.shape[0])],
                                   index=process_data.index)
        output_df = pd.concat([output_df, Sv_expanded], axis=1)
    
    # Add classification column
    output_df["Process_Type"] = process_types
    
    # Create output directory if needed
    os.makedirs(os.path.dirname(filepath) if os.path.dirname(filepath) else '.', exist_ok=True)
    
    # Export based on format
    if format.lower() == 'excel':
        with pd.ExcelWriter(filepath, engine="openpyxl") as writer:
            # Main sheet with all data
            output_df.to_excel(writer, sheet_name="All_Processes", index=False)
            
            # Separate sheets by category if requested
            if separate_sheets:
                for category, group in output_df.groupby("Process_Type", sort=False):
                    # Excel sheet names limited to 31 characters
                    sheet_name = category[:31]
                    group.to_excel(writer, sheet_name=sheet_name, index=False)
        
        print(f"Procesos clasificados guardados en: {filepath}")
    
    elif format.lower() == 'csv':
        output_df.to_csv(filepath, index=False)
        print(f"Procesos clasificados guardados en: {filepath}")
    
    else:
        raise ValueError(f"Formato no soportado: {format}. Use 'excel' o 'csv'.")
    
    return filepath


######################################################################################
# SUMMARY
######################################################################################

__all__ = [
    # Process Classification
    'classify_process_mode',
    'is_cognitive_domain',
    'classify_response_to_disturbance',
    
    # Time Series Processing
    'rescale_process_time_series',
    'classify_process_series',
    
    # Cone Analysis
    'compute_nullspace_vectors',
    'compute_feasible_region',
    'classify_feasible_points',
    'analyze_cone',
    
    # Utilities
    'get_intervals_by_category',
    'export_classified_processes'
]