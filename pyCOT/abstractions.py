import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt 

###################################################################################
# Function to generate the abstraction of a time series based on a threshold
def abstraction_ordinary(time_series, threshold=None, title="Abstraction"):
    """
    Builds an abstraction of a time series based on a threshold.

    Parameters:
    time_series (pd.DataFrame): DataFrame where each column represents a species and rows correspond to time points.
                                Must include a 'Time' column with time values.
    threshold (float or dict): Threshold to determine the existence of a species.
                               If a dictionary, defines a different threshold for each species.
    title (string): Name of the abstraction column.

    Returns:
    abstract_time_series: Abstract time series where each row specifies the species meeting the existence condition.
    """

    # Ensure the 'Time' column is present
    if 'Time' not in time_series.columns:
        raise ValueError("The DataFrame must include a 'Time' column for time values.")

    # Assign a default threshold if none is provided
    if threshold is None:
        threshold = 0.5

    # Check if the threshold is a dictionary or a single value
    if isinstance(threshold, dict):
        # Validate that all species have a defined threshold
        for species in time_series.columns:
            if species != 'Time' and species not in threshold:
                raise ValueError(f"The threshold for species '{species}' is not defined.")
    else:
        # Convert the threshold into a uniform dictionary
        threshold = {species: threshold for species in time_series.columns if species != 'Time'}
    
    # Build the abstraction
    abstract_series = []
    for idx, row in time_series.iterrows():
        # List of species meeting the existence condition at this time point
        active_species = [species for species in time_series.columns if species != 'Time' and row[species] >= threshold[species]]
        # Save the time and active species for this time point
        abstract_series.append({"Time": row['Time'], title: active_species})
    
    # Convert to DataFrame
    abstract_time_series = pd.DataFrame(abstract_series)
    return abstract_time_series