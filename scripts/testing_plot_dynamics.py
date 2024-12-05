#######################################################################
import os
import sys
# Add the root directory to the PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pyCOT.abstractions import abstraction_ordinary
from pyCOT.plot_dynamics import *

#########################################################################
# Time series of species concentrations
data_concentrations = {
    "Time": [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0],
    "A": [0.8953792537708223, 0.26486054747788695, 0.7814640741363406, 0.4804573804243072, 0.4485337249500629, 0.4629068467286066, 0.2454758964704813, 0.5405723189720857, 0.6591293320543123, 0.49285372534630484],
    "B": [0.2649035608580179, 0.5035839068288304, 0.44496204400906814, 0.1131711028253838, 0.5859495343929363, 0.5620931096408655, 0.5474624930955375, 0.48221947601400783, 0.40781057137176213, 0.2017722955105233],
    "C": [0.543215778204073, 0.9256885528146014, 0.5166438102659423, 0.37965876777582463, 0.2846592197931128, 0.19985473169947682, 0.3628683818281516, 0.9723371590709325, 0.7496519452985205, 0.5110804891140438]
}

time_series = pd.DataFrame(data_concentrations)
print(time_series)

# Define a fixed threshold or per-species thresholds
threshold = 0.5  # Fixed threshold
# threshold = {"A": 0.4, "B": 0.5, "C": 0.4}  # Per-species thresholds

# Calculate abstractions
# abstract_time_series = abstraction_ordinary(time_series)
abstract_time_series = abstraction_ordinary(time_series, threshold)
print(abstract_time_series)
#########################################################################

#########################################################################
# Examples of Plots
#########################################################################
# # Plot of the time series
# Plot the time series of species concentrations
plot_series_ode(time_series)

# # Plot of the number of abstractions
plot_abstraction_size(abstract_time_series)

# # Plot the abstractions sorted by size
plot_abstraction_sets(abstract_time_series)

###########################################################################
# # Graphs
# # Static graph
# plot_static_abstraction_graph(abstract_time_series)

# # Movie
# plot_abstraction_graph_movie(abstract_time_series, interval=1000)

###########################################################################
# # Histograms
# plot_join_concentration_histogram(time_series, bins=5, alpha=0.6)

# # Call the function to plot the histograms
# plot_species_histograms(time_series, species_names=["A", "B", "C"], bins=5, alpha=0.6)

# # Call the function to plot the combined histogram
# plot_combined_species_histogram(time_series, species_names=["A", "B", "C"], bins=5, alpha=0.6)
