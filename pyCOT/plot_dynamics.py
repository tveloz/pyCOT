from pyvis.network import Network        # Enables creating and visualizing interactive networks in HTML.
import networkx as nx                    # Library for creating, manipulating, and analyzing graphs.
import pandas as pd                      # Efficient handling and analysis of tabular data.
import matplotlib.pyplot as plt          # For generating static plots and visualizations.
import matplotlib.animation as animation # For creating animations using Matplotlib.
from matplotlib.colors import to_hex     # Converts colors to hexadecimal format (#RRGGBB).
import mplcursors                        # Adds interactive cursors to Matplotlib plots.
import webbrowser                        # Opens URLs or local files in the system's default web browser.
import os                                # Handles file and directory operations in the operating system.
from collections import defaultdict      # Dictionary that provides default values for missing keys.
import plotly.graph_objects as go        # Generates interactive plots and advanced visualizations using Plotly.
from itertools import combinations       # Generates all possible combinations of elements in an iterable.

import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout

import sys                               # Provides access to system-specific parameters and functions.
sys.stdout.reconfigure(encoding='utf-8') # Reconfigures the standard output to use UTF-8 encoding, ensuring proper handling of special characters.
import tempfile                          # Provides utilities for creating temporary files and directories.

 


######################################################################################
# Plots the time series of ODE concentrations and abstractions
######################################################################################

def plot_series_ode(time_series, xlabel="Time", ylabel="Concentration", title="Time Series of Concentrations"):
    """
    Plots the time series of ODE concentrations.

    Parameters:
    time_series (pd.DataFrame): Time series with a 'Time' column and species concentrations as columns.
    xlabel (str): Label for the x-axis. Default is "Time".
    ylabel (str): Label for the y-axis. Default is "Concentration".
    title (str): Title of the plot. Default is "Time Series of Concentrations".

    Raises:
    ValueError: If the DataFrame does not contain a 'Time' column.

    Returns:
    None: Displays a line plot for each species in the time series.
    """
    if 'Time' not in time_series.columns:
        # Check if the 'Time' column exists, raise error if not found
        raise ValueError("The DataFrame must include a 'Time' column for time values.")
    
    # Create a new figure and axis object for the plot
    fig, ax = plt.subplots(figsize=(10, 6))  # Plot size 10x6 inches
    
    # Iterate over the columns to plot each species
    for species in time_series.columns:
        if species != 'Time':  # Skip the 'Time' column
            ax.plot(time_series['Time'], time_series[species], label=species)  # Plot species concentration
    
    # Set axis labels and plot title
    ax.set_xlabel(xlabel)  # Set x-axis label
    ax.set_ylabel(ylabel)  # Set y-axis label
    ax.set_title(title)    # Set the title of the plot
    ax.grid()              # Enable grid on the plot
    ax.legend()            # Add a legend for species
    
    plt.tight_layout()  # Adjust layout to fit all elements
    plt.show()          # Display the plot

######################################################################################## 

def plot_abstraction_size(abstract_time_series, xlabel="Time", ylabel="Number of Species", title="Number of species per abstraction over time", marker='o', label="Abstraction Size"):
    """
    Plots the number of abstractions over the time series.

    Parameters:
    abstract_time_series (pd.DataFrame): Time series with a 'Time' column and a column for species abstractions.
    xlabel (str): Label for the x-axis. Default is "Time".
    ylabel (str): Label for the y-axis. Default is "Number of Species".
    title (str): Title of the plot. Default is "Number of species per abstraction over time".
    marker (str): Marker style for the plot. Default is 'o'.
    label (str): Legend label for the plot. Default is "Abstraction Size".

    Raises:
    ValueError: If the DataFrame does not contain a 'Time' column.

    Returns:
    None: Displays a line plot of abstraction sizes over time.
    """
    if 'Time' not in abstract_time_series.columns:
        # Check if the 'Time' column exists, raise error if not found
        raise ValueError("The DataFrame must include a 'Time' column for time values.")
    
    # Create a new figure and axis object for the plot
    fig, ax = plt.subplots(figsize=(10, 6))  # Plot size 10x6 inches
    
    # Create a mapping of time and number of active species (abstraction size)
    abstraction_mapping = {
        idx: (row["Time"], len(row["Abstraction"]))  # Map index to (time, abstraction size)
        for idx, row in abstract_time_series.iterrows()
    }
    
    # Extract times and abstraction sizes
    times = [abstraction_mapping[idx][0] for idx in abstraction_mapping]  # Extract time values
    sizes = [abstraction_mapping[idx][1] for idx in abstraction_mapping]  # Extract abstraction sizes

    # Plot abstractions as ordered points
    ax.plot(times, sizes, marker=marker, linestyle='-', label=label)  # Plot with line and marker
    ax.set_xlabel(xlabel)  # Set x-axis label
    ax.set_ylabel(ylabel)  # Set y-axis label
    ax.set_title(title)  # Set plot title
    ax.grid()  # Enable grid on the plot
    ax.legend()  # Add a legend for the abstraction size
    
    plt.tight_layout()  # Adjust layout to fit all elements
    plt.show()  # Display the plot

######################################################################################## 

def plot_abstraction_sets(abstract_time_series, xlabel="Time",ylabel="Species",title="Abstraction of Time Series"):
    """
    Plots the time series abstraction generated by the abstraction_ordinary function.

    Parameters:
    abstract_time_series (pd.DataFrame): DataFrame containing 'Time' and 'Abstraction' columns.
                                         'Abstraction' should be a list of species present at each time point.
    xlabel (str): Label for the x-axis. Default is "Time".
    ylabel (str): Label for the y-axis. Default is "Presence of Species".
    title (str): Title of the plot. Default is "Abstraction of ODE Time Series".

    Returns:
    None: Displays a stacked bar chart showing the presence of species over time.
    """
    # Extract all unique species from the 'Abstraction' column
    all_species = sorted({species for row in abstract_time_series["Abstraction"] for species in row})  # Unique sorted species
    
    # Build a binary matrix indicating the presence of species over time
    binary_matrix = pd.DataFrame(
        [
            {species: (species in row) for species in all_species}  # Check presence of each species in a row
            for row in abstract_time_series["Abstraction"]
        ],
        index=abstract_time_series["Time"]  # Use time as the index
    ).astype(int)  # Convert to integer (0 or 1)
    
    # Create a new figure and axis object for the plot
    fig, ax = plt.subplots(figsize=(10, 6))  # Plot size 10x6 inches
    
    # Plot each species as a stacked bar chart
    for i, species in enumerate(all_species):
        ax.bar(
            binary_matrix.index,  # Time points
            binary_matrix[species],  # Presence of species
            label=species,  # Label for legend
            bottom=binary_matrix.iloc[:, :i].sum(axis=1)  # Stacked position
        )
    
    # Configure plot labels and title
    ax.set_xlabel(xlabel)  # Set x-axis label
    ax.set_ylabel(ylabel)  # Set y-axis label
    ax.set_title(title)  # Set plot title
    ax.legend(title="Species")  # Add legend with title
    
    plt.show()  # Display the plot

######################################################################################
# Abstraction graph static and movie 
######################################################################################
# # # Original: plot_static_abstraction_graph
def plot_static_abstraction_graph(abstract_time_series, title="Static Abstraction Graph"):
    """
    Plots a static abstraction graph with nodes, edges, and associated attributes.

    The graph represents:
    - Node size: Depends on the frequency of occurrences of each abstraction.
    - Node color: Represents the weighted average time of occurrence (normalized).
    - Edge thickness: Indicates the frequency of transitions between abstractions.

    Parameters:
    abstract_time_series (pd.DataFrame): A DataFrame with the following columns:
        - 'Time': Time of occurrence for each abstraction (numeric).
        - 'Abstraction': Abstraction sets represented as iterables.
    title (str): Title of the graph. Defaults to "Static Abstraction Graph".

    Returns:
    None. Displays a static plot of the abstraction graph.
    """
    abstractions = abstract_time_series["Abstraction"]  # Extracts the abstraction column
    times = abstract_time_series["Time"]  # Extracts the time column
    nodes = abstractions.apply(tuple).value_counts()  # Counts the occurrences of each abstraction

    # Compute transition frequencies
    transitions = [(tuple(abstractions[i]), tuple(abstractions[i + 1])) for i in range(len(abstractions) - 1)]  # Transitions between consecutive abstractions
    transitions_freq = pd.Series(transitions).value_counts()  # Counts the frequency of each transition

    # Create the graph
    G = nx.DiGraph()  # Creates a directed graph

    # Add nodes with attributes
    for node, freq in nodes.items():
        node_times = times[abstractions.apply(tuple) == node]  # Filters the times associated with the current node
        weighted_time = (node_times * freq).mean()  # Calculates the weighted average time
        normalized_time = (weighted_time - times.min()) / (times.max() - times.min())  # Normalizes the time
        # color = plt.cm.coolwarm(1 - normalized_time)  # Assigns a color based on the normalized time
        color = to_hex(plt.cm.coolwarm(1 - normalized_time))  # Convierte el color a formato hexadecimal
        G.add_node(node, size=freq, color=color)  # Adds the node to the graph with attributes

    # Add edges with weights
    for (source, target), weight in transitions_freq.items():
        G.add_edge(source, target, weight=weight)  # Adds edges to the graph with the corresponding weight

    # Hierarchical layout
    pos = nx.drawing.nx_agraph.graphviz_layout(G, prog="dot")  # Calculates the node positions in a hierarchical layout
    pos = {node: (x, -y) for node, (x, y) in pos.items()}  # Inverts the y-axis so the graph is drawn from bottom to top

    # Create the plot
    plt.figure(figsize=(12, 8))  # Creates a figure with a specific size

    # Draw nodes
    node_sizes = [G.nodes[node]["size"] * 100 for node in G.nodes]  # Scales the node sizes
    node_colors = [G.nodes[node]["color"] for node in G.nodes]  # Gets the node colors
    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors)  # Draws the nodes

    # Draw edges
    edge_weights = [G.edges[edge]["weight"] for edge in G.edges]  # Gets the edge weights
    nx.draw_networkx_edges(G, pos, width=[w / 2 for w in edge_weights], alpha=0.7)  # Draws the edges with adjusted opacity

    # Draw node labels
    # nx.draw_networkx_labels(G, pos, font_size=8, font_color="black")  # Draws node labels 
    node_labels = {node: f"({', '.join(map(str, node))})" for node in G.nodes}  # Formatea las etiquetas
    nx.draw_networkx_labels(G, pos, labels=node_labels, font_size=8, font_color="black")  # Usa las etiquetas formateadas

    # Title and display
    plt.title(title)  # Adds a title to the graph
    plt.axis("off")  # Hides the axes
    plt.show()  # Displays the graph

########################################################################################

def get_plot_static_abstraction_graph_html(
    abstract_time_series, 
    node_size=20,  
    edge_color="gray", 
    shape_node='dot', 
    node_font_size=14, 
    edge_width=2,  
    filename="Static_Abstraction_Graph.html"
):
    """
    Generates an interactive HTML visualization of a directed graph.

    The graph is constructed from a time series of abstractions, with nodes representing unique abstractions
    and edges representing transitions between consecutive abstractions. Node attributes are customized
    based on their frequency, and edge attributes are customized based on transition frequencies.

    Parameters:
    ----------
    abstract_time_series : pd.DataFrame
        A DataFrame containing two columns:
        - "Abstraction": A column with abstraction data (e.g., lists or tuples).
        - "Time": A column with corresponding time values for each abstraction.
    node_size : int, optional
        The size of the nodes in the visualization (default is 20).
    edge_color : str, optional
        The color of the edges in the graph (default is "gray").
    shape_node : str, optional
        The shape of the nodes (e.g., 'dot', 'square'; default is 'dot').
    node_font_size : int, optional
        The font size for node labels (default is 14).
    edge_width : int, optional
        The base width of the edges (default is 2).
    filename : str, optional
        The name of the output HTML file (default is "Static_Abstraction_Graph.html").

    Returns:
    -------
    pyvis.network.Network
        A PyVis Network object containing the interactive visualization.

    Notes:
    -----
    - The node color encodes its position in the sequence:
      * Blue: Initial node
      * Red: Final node
      * Gradient colors (coolwarm colormap): Intermediate nodes based on normalized time.
    - Nodes and edges are added to the visualization with specific attributes such as size, color, 
      label, and weight.
    - The output file is an interactive HTML visualization saved at the specified filename.

    Example:
    -------
    >>> df = pd.DataFrame({
    >>>     "Time": [0, 1, 2, 3, 4],
    >>>     "Abstraction": [[A], [A,B], [C], [A,C], [A]]    
    >>> })
    >>> get_plot_static_abstraction_graph_html(df) 
    """
    # Extract abstraction and time columns
    abstractions = abstract_time_series["Abstraction"]  
    times = abstract_time_series["Time"]  
    
    # Count node frequencies
    nodes = abstractions.apply(tuple).value_counts(sort=False)  
    # print("Nodes de", nodes)

    # Compute transition frequencies
    transitions = [(tuple(abstractions[i]), tuple(abstractions[i + 1])) for i in range(len(abstractions) - 1)]  
    transitions_freq = pd.Series(transitions).value_counts()  
    # print("transitions", transitions)
    # print("transitions_freq", transitions_freq)
     
    # Create the graph
    G = nx.DiGraph()  

    # Add nodes with attributes
    for node, freq in nodes.items():
        node_labels = {node: f"S{idx + 1}" for idx, node in enumerate(nodes.keys())}
        node_times = times[abstractions.apply(tuple) == node]  
        weighted_time = (node_times * freq).mean()  
        normalized_time = (weighted_time - times.min()) / (times.max() - times.min())  

        if node == tuple(abstractions.iloc[0]):  
            color = "#0000FF"  
        elif node == tuple(abstractions.iloc[-1]):  
            color = "#FF0000"  
        else:
            color = to_hex(plt.cm.coolwarm(1 - normalized_time))  

        label = f"{node_labels[node]} ({freq})"  
        G.add_node(node, size=freq, color=color, label=label)  

    # Add edges with attributes
    for (source, target), freq in transitions_freq.items():
        G.add_edge(
            source, 
            target, 
            weight=freq,  
            color=edge_color,  
            width=edge_width + freq / max(transitions_freq) * edge_width  
        )

    # Convert to PyVis visualization 
    net = Network(height="750px", width="100%", directed=True, notebook=False)
    net.set_options(f"""
    {{
      "nodes": {{
        "shape": "{shape_node}",  
        "font": {{
        "size": {node_font_size},  
        "align": "center"  
        }},
        "borderWidth": 2,  
        "borderColor": "black"  
      }},
      "edges": {{
        "smooth": false,
        "color": "{edge_color}",
        "width": {edge_width}  
      }},
      "physics": {{
        "enabled": false,
        "stabilization": {{
          "enabled": false
        }},
        "hierarchicalRepulsion": {{
          "nodeDistance": 150
        }}
      }},
      "layout": {{
        "hierarchical": {{
          "enabled": false,
          "direction": "DU",  
          "sortMethod": "directed"
        }}
      }}
    }}
    """) 

    # Add nodes to the PyVis network
    for node, attrs in G.nodes(data=True):
        net.add_node(
            node_labels[node],  
            label=attrs['label'],  
            title=f"{node}",
            color=attrs['color'],
            size=node_size
        )

    # Add edges to the PyVis network
    for source, target, attrs in G.edges(data=True):
        net.add_edge(
            node_labels[source],  
            node_labels[target],  
            width=edge_width
        )
        
    # Save the visualization to an HTML file 
    net.html = net.generate_html() 
    with open(filename, "w", encoding="utf-8") as f:
        f.write(net.html)
    return net

########################################################################################

def plot_static_abstraction_graph_html(
    abstract_time_series, 
    node_size=20,  
    edge_color="gray", 
    shape_node='dot', 
    node_font_size=14, 
    edge_width=2,  
    filename="Static_Abstraction_Graph.html"
):
    """
    Generate and visualize the static abstraction graph as an interactive HTML file.

    This function is a wrapper to visualize the containment hierarchy among sets
    using a static abstraction graph. The graph is saved as an HTML file and can be 
    opened in a web browser.

    Args:
        abstract_time_series (pd.DataFrame): DataFrame containing abstraction and time series information.
            It must include at least the following columns:
            - "Abstraction": Represents the sets or elements to visualize.
            - "Time": Represents the corresponding time for the abstractions.
        node_size (int, optional): Size of the nodes in the graph. Default is 20.
        edge_color (str, optional): Color of the edges connecting nodes. Default is "gray".
        shape_node (str, optional): Shape of the nodes. Options include 'dot', 'square', 
            'triangle', 'star', or 'diamond'. Default is 'dot'.
        node_font_size (int, optional): Font size for node labels. Default is 14.
        edge_width (int, optional): Thickness of the edges. Default is 2.
        filename (str, optional): Name of the output HTML file. Default is "Static_Abstraction_Graph.html".

    Raises:
        FileNotFoundError: If the generated HTML file is not found in the specified location.

    Side Effects:
        - Saves the generated graph as an HTML file at the specified location.
        - Opens the HTML file in the default web browser.

    Example:
        ```python
        plot_static_abstraction_graph_html(
            abstract_time_series=my_dataframe,
            node_size=30,
            edge_color="blue",
            shape_node="triangle",
            filename="My_Abstraction_Graph.html"
        )
        ```
    """
    # Call the `get_plot_static_abstraction_graph_html` function to generate the graph
    get_plot_static_abstraction_graph_html(
        abstract_time_series,  
        node_size=node_size,  
        edge_color=edge_color, 
        shape_node=shape_node, 
        node_font_size=node_font_size, 
        edge_width=edge_width,  
        filename=filename
    )
    
    # Convert to an absolute path
    abs_path = os.path.abspath(filename) 
    
    # Check if the file was created correctly
    if not os.path.isfile(abs_path):
        print(f"File not found at {abs_path}")  # Additional message for debugging
        raise FileNotFoundError(f"The file {abs_path} was not found. Check if `get_plot_static_abstraction_graph_html` generated the file correctly.")
    
    # Inform the user about the file's location
    print(f"The hierarchy visualization was saved to: {abs_path}")
    
    # Open the HTML file in the default browser
    webbrowser.open(f"file://{abs_path}")

########################################################################################

def get_plot_static_abstraction_graph_hierarchy_html(
    abstract_time_series, 
    node_size=20,  
    edge_color="gray", 
    shape_node='dot', 
    node_font_size=14, 
    edge_width=2,  
    filename="Static_Abstraction_Graph_Hierarchy.html"
):
    """
    Generates a hierarchical abstraction graph from time series data and saves it as an interactive HTML file.

    Args:
        abstract_time_series (pd.DataFrame): A DataFrame with two columns:
            - "Abstraction": Contains sets or lists representing the abstract states.
            - "Time": Contains the corresponding timestamps for each state.
        node_size (int): Size of the nodes in the graph visualization.
        edge_color (str): Color of the edges in the graph.
        shape_node (str): Shape of the nodes; options include 'dot', 'square', 'triangle', etc.
        node_font_size (int): Font size for the node labels.
        edge_width (int): Base width of the edges in the graph.
        filename (str): Name of the output HTML file where the visualization will be saved.

    Returns:
        net (pyvis.network.Network): A PyVis Network object representing the generated graph.

    Description:
        - Extracts nodes and their frequencies from the "Abstraction" column of the time series data.
        - Computes containment relationships among nodes to build a directed graph.
        - Nodes are colored based on their temporal positions (initial, final, or normalized weighted time).
        - Saves the graph as an interactive HTML file with hierarchical layout and opens it in the browser.

    Example:
        ```python
        abstract_time_series = pd.DataFrame({
            "Abstraction": [{"A", "B"}, {"A"}, {"A", "C"}, {"C"}],
            "Time": [1, 2, 3, 4]
        })
        get_plot_static_abstraction_graph_hierarchy_html(
            abstract_time_series,
            node_size=30,
            edge_color="blue",
            shape_node="square",
            filename="Hierarchy_Graph.html"
        )
        ```
    """
    abstractions = abstract_time_series["Abstraction"]  # Extracts the abstraction column
    times = abstract_time_series["Time"]  # Extracts the time column
    nodes = abstractions.apply(tuple).value_counts(sort=False)  # Counts the occurrences of each abstraction    
    # print("nodes", nodes)

    # Compute transition frequencies
    transitions = [(tuple(abstractions[i]), tuple(abstractions[i + 1])) for i in range(len(abstractions) - 1)]  # Transitions between consecutive abstractions
    transitions_freq = pd.Series(transitions).value_counts()  # Counts the frequency of each transition
    # print("transitions", transitions)
    # print("transitions_freq", transitions_freq)
    
    # Build the containment graph
    containment_graph = {node: [] for node in nodes.index}  # Initialize containment graph
    for i, child_set in enumerate(nodes.index):  # Iterate over sets as children
        for j, parent_set in enumerate(nodes.index):  # Iterate over sets as parents
            if i != j and set(child_set).issubset(parent_set):  # If they are not the same and containment exists
                if parent_set not in containment_graph[child_set]:  # Check subsets of child_set
                    containment_graph[child_set].append(parent_set)  # Add the relationship
    # print("containment_graph", containment_graph)

    # Create the graph
    G = nx.DiGraph()  # Directed graph

    # Create unique labels for nodes (e.g., S1, S2, S3...)
    node_labels = {node: f"S{idx + 1}" for idx, node in enumerate(nodes.keys())}

    # Add nodes with attributes 
    for node, freq in nodes.items():  
        node_times = times[abstractions.apply(tuple) == node]  # Filters the times associated with the current node
        weighted_time = (node_times * freq).mean()  # Calculates the weighted average time
        normalized_time = (weighted_time - times.min()) / (times.max() - times.min())  # Normalizes the time

        if node == tuple(abstractions.iloc[0]):  # Initial node
            color = "#0000FF"  # Blue
        elif node == tuple(abstractions.iloc[-1]):  # Final node
            color = "#FF0000"  # Red
        else:
            color = to_hex(plt.cm.coolwarm(1 - normalized_time))  # Color based on normalized time

        label = f"{node_labels[node]} ({freq})"  # Label includes the unique label, frequency, and elements
        G.add_node(node, size=freq, color=color, label=label)  # Add node with attributes

    # print("nodes items:", nodes)

    # Add edges based on containment graph (parents -> children) 
    for child, parents in containment_graph.items(): 
        for parent in parents:
            G.add_edge(child, parent, weight=1)  # Add directed edges from parent to child

    # Initialize the PyVis network
    net = Network(height="750px", width="100%", directed=True, notebook=False)
    net.set_options(f"""
    {{
      "nodes": {{
        "shape": "{shape_node}",  
        "font": {{
          "size": {node_font_size},  
          "align": "center"  
        }},
        "borderWidth": 2,  
        "borderColor": "black"  
      }},
      "edges": {{
        "smooth": false,
        "color": "{edge_color}",
        "width": {edge_width}  
      }},
      "physics": {{
        "enabled": false,
        "stabilization": {{
          "enabled": false
        }},
        "hierarchicalRepulsion": {{
          "nodeDistance": 150
        }}
      }},
      "layout": {{
        "hierarchical": {{
          "enabled": true,
          "direction": "DU",  
          "sortMethod": "directed"
        }}
      }}
    }}
    """)

    # Add nodes to the PyVis network
    for node, attrs in G.nodes(data=True):
        net.add_node(
            node_labels[node],  # Use the unique label for the node
            label=attrs['label'],  # Set the label with frequency and elements
            title=f"{node}",
            color=attrs['color'],
            size=node_size
        )

    # Add edges to the PyVis network
    for source, target, attrs in G.edges(data=True):
        net.add_edge(
            node_labels[source],  # Use unique labels for source
            node_labels[target],  # Use unique labels for target
            width=edge_width
        )

    # Save the visualization to an HTML file 
    net.html = net.generate_html() 
    with open(filename, "w", encoding="utf-8") as f:
        f.write(net.html)
    return net

########################################################################################

def plot_static_abstraction_graph_hierarchy_html(
    abstract_time_series, 
    node_size=20,  
    edge_color="gray", 
    shape_node='dot', 
    node_font_size=14, 
    edge_width=2,  
    filename="Static_Abstraction_Graph_Hierarchy.html"
):        
    """
    Generates and visualizes the containment hierarchy of abstractions as an interactive HTML file.
    
    This function serves as a wrapper to call the underlying function 
    `get_plot_static_abstraction_graph_hierarchy_html`, which creates the graph visualization. 
    Once the graph is generated, it saves the visualization to an HTML file and optionally opens it 
    in the user's default web browser.
    
    Args:
        abstract_time_series (pd.DataFrame): A DataFrame containing columns "Abstraction" (list or set of 
            elements) and "Time" (time associated with the abstractions). 
        node_size (int, optional): Size of the nodes in the visualization. Default is 20.
        edge_color (str, optional): Color of the edges in the graph. Default is "gray".
        shape_node (str, optional): Shape of the nodes in the graph. Options include 'dot', 'square', 
            'triangle', 'star', and 'diamond'. Default is 'dot'.
        node_font_size (int, optional): Font size of the labels displayed on nodes. Default is 14.
        edge_width (int, optional): Width of the edges connecting nodes. Default is 2.
        filename (str, optional): Name of the output HTML file to save the graph visualization. Default is 
            "Static_Abstraction_Graph_Hierarchy.html".
    
    Raises:
        FileNotFoundError: If the HTML file is not found after attempting to generate it.
    
    Outputs:
        The function generates an interactive HTML file that represents the abstraction graph hierarchy 
        and opens it in the default web browser.
    
    Example:
        >>> plot_static_abstraction_graph_hierarchy_html(
        ...     abstract_time_series=data,
        ...     node_size=30,
        ...     edge_color="blue",
        ...     shape_node="triangle",
        ...     filename="AbstractionGraph.html"
        ... )
    """
    # Call the function to generate the visualization
    get_plot_static_abstraction_graph_hierarchy_html(
        abstract_time_series,  
        node_size=node_size,  
        edge_color=edge_color, 
        shape_node=shape_node, 
        node_font_size=node_font_size, 
        edge_width=edge_width,  
        filename=filename
    )
    
    # Convert to an absolute path
    abs_path = os.path.abspath(filename) 
    
    # Check if the file was created correctly
    if not os.path.isfile(abs_path):
        print(f"File not found at {abs_path}")  # Additional message for debugging
        raise FileNotFoundError(f"The file {abs_path} was not found. Check if `get_plot_static_abstraction_graph_hierarchy_html` generated the file correctly.")
    
    # Inform the user about the file's location
    print(f"The hierarchy visualization was saved to: {abs_path}")
    
    # Open the HTML file in the default browser
    webbrowser.open(f"file://{abs_path}")

########################################################################################

def plot_static_abstraction_graph_hierarchy(abstract_time_series, title="Hierarchical Abstraction Graph"):
    abstractions = abstract_time_series["Abstraction"]
    times = abstract_time_series["Time"]
    nodes = abstractions.apply(tuple).value_counts()

    # Compute transition frequencies
    transitions = [(tuple(abstractions[i]), tuple(abstractions[i + 1])) for i in range(len(abstractions) - 1)]
    transitions_freq = pd.Series(transitions).value_counts()

    # Create a directed graph
    G = nx.DiGraph()

    # Add nodes with attributes
    for node, freq in nodes.items():
        node_times = times[abstractions.apply(tuple) == node]
        weighted_time = (node_times * freq).mean()
        normalized_time = (weighted_time - times.min()) / (times.max() - times.min())
        color = to_hex(plt.cm.coolwarm(1 - normalized_time))
        G.add_node(node, size=freq, color=color)

    # Add edges for transitions
    for (source, target), weight in transitions_freq.items():
        G.add_edge(source, target, weight=weight)

    # Add edges based on containment relationships
    node_list = list(G.nodes)
    for i, smaller in enumerate(node_list):
        for j, larger in enumerate(node_list):
            if i != j and set(smaller).issubset(set(larger)):  # Containment relationship
                if not G.has_edge(smaller, larger):  # Avoid duplicate edges
                    G.add_edge(smaller, larger, weight=1)  # Default weight

    # Compute hierarchical positions based on containment
    levels = {}
    for node in G.nodes:
        level = sum(1 for other in G.nodes if set(node).issubset(set(other)) and node != other)
        levels[node] = level

    pos = {node: (0, -level) for node, level in levels.items()}

    # Create the plot
    plt.figure(figsize=(12, 8))

    # Draw nodes
    node_sizes = [G.nodes[node]["size"] * 100 for node in G.nodes]
    node_colors = [G.nodes[node]["color"] for node in G.nodes]
    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors)

    # Draw edges
    edge_weights = [G.edges[edge].get("weight", 1) for edge in G.edges]  # Use default weight
    nx.draw_networkx_edges(G, pos, width=[w / 2 for w in edge_weights], alpha=0.7)

    # Draw node labels
    node_labels = {node: f"({', '.join(map(str, node))})" for node in G.nodes}
    nx.draw_networkx_labels(G, pos, labels=node_labels, font_size=8, font_color="black")

    # Title and display
    plt.title(title)
    plt.axis("off")
    plt.show()

########################################################################################

def plot_static_abstraction_graph_html_with_hierarchy(
    abstract_time_series,
    filename="static_abstraction_graph_hierarchy.html",
    node_size_scale=5,
    edge_width_scale=0.5,
    default_color="cyan",
    title="plot_static_abstraction_graph_html_with_hierarchy"
):
    """
    Generates an HTML visualization of a static abstraction graph considering hierarchy based on subset relationships.

    Parameters:
    abstract_time_series (pd.DataFrame): DataFrame with columns:
        - 'Time': Time of occurrence for each abstraction (numeric).
        - 'Abstraction': Abstraction sets represented as iterables.
    filename (str): Name of the output HTML file. Defaults to "static_abstraction_graph_hierarchy.html".
    node_size_scale (float): Scaling factor for node sizes.
    edge_width_scale (float): Scaling factor for edge widths.
    default_color (str): Default color for nodes.
    title (str): Title of the graph. Defaults to "Static Abstraction Graph with Hierarchy".

    Returns:
    str: Filename of the generated HTML file.
    """
    # Extract data from DataFrame
    abstractions = abstract_time_series["Abstraction"]
    times = abstract_time_series["Time"]
    nodes = abstractions.apply(tuple).value_counts()

    # Compute transition frequencies
    transitions = [(tuple(abstractions[i]), tuple(abstractions[i + 1])) for i in range(len(abstractions) - 1)]
    transitions_freq = pd.Series(transitions).value_counts()

    # Create hierarchical levels based on subset relationships
    hierarchy = {}
    sorted_nodes = sorted(nodes.keys(), key=lambda x: (len(x), x))  # Sort by size of the set and lexicographically
    for node in sorted_nodes:
        hierarchy[node] = []
        for potential_parent in sorted_nodes:
            if set(node).issubset(set(potential_parent)) and node != potential_parent:
                hierarchy[node].append(potential_parent)

    # Initialize the PyVis network
    net = Network(height="750px", width="100%", directed=True, notebook=False)
    net.set_options(f"""
    {{
      "edges": {{
        "smooth": false,
        "color": "gray"
      }},
      "physics": {{
        "enabled": false,
        "stabilization": {{
          "enabled": true
        }}
      }},
      "layout": {{
        "hierarchical": {{
          "enabled": false,
          "direction": "DU",
          "sortMethod": "directed"
        }}
      }}
    }}
    """)

    # Add nodes
    for node, freq in nodes.items():
        node_times = times[abstractions.apply(tuple) == node]
        weighted_time = (node_times * freq).mean()
        normalized_time = (weighted_time - times.min()) / (times.max() - times.min())
        color = to_hex(plt.cm.coolwarm(1 - normalized_time)) if not pd.isna(normalized_time) else default_color
        net.add_node(
            str(node),
            label=f"{node}",
            title=f"Freq: {freq}, Avg Time: {weighted_time:.2f}" if not pd.isna(weighted_time) else "N/A",
            color=color,
            size=freq * node_size_scale
        )

    # # Add edges based on hierarchy
    # for node, parents in hierarchy.items():
    #     for parent in parents:
    #         net.add_edge(str(parent), str(node), width=1, color="black", title="Subset Relationship")

    # Add transitions as additional edges
    for (source, target), weight in transitions_freq.items():
        net.add_edge(str(source), str(target), weight=weight * edge_width_scale, title=f"Freq: {weight}", color="gray")

    # Save the visualization to an HTML file 
    net.html = net.generate_html()
    with open(filename, "w", encoding="utf-8") as f:
        f.write(net.html)
    return filename

########################################################################################
########################################################################################
########################################################################################
# # # Original: plot_abstraction_graph_movie
def plot_abstraction_graph_movie(abstract_time_series, interval=400, title="Abstraction Graph - Time"):
    """
    Creates an animated abstraction graph that evolves over time.

    The animation highlights:
    - The current node in red.
    - Past nodes gradually fading (transparent) as time progresses.
    - The last three unique abstraction sets in varying sizes and shades of blue.
    - Node size: Depends on the frequency of occurrences of each abstraction.
    - Node color: Represents the weighted average time of occurrence (normalized).
    - Edge thickness: Indicates the frequency of transitions between abstractions.

    Parameters:
    abstract_time_series (pd.DataFrame): A DataFrame with the following columns:
        - 'Time': Time of occurrence for each abstraction (numeric).
        - 'Abstraction': Abstraction sets represented as iterables.
    interval (int): Interval in milliseconds between frames of the animation. Defaults to 400.
    title (str): Title of the graph animation. Defaults to "Abstraction Graph - Time".

    Returns:
    None. Displays an animation of the abstraction graph.
    """
    abstractions = abstract_time_series["Abstraction"]
    times = abstract_time_series["Time"]
    unique_abstractions = abstractions.apply(tuple).drop_duplicates(keep='last')  # Keep last occurrences of unique abstractions
    last_three_unique = unique_abstractions.tail(3)  # Last three unique sets

    nodes = abstractions.apply(tuple).value_counts(sort=False)
    # Ordenar el Series por la longitud de la primera coordenada de las claves
    # nodes= nodes.sort_index(key=lambda x: [len(k) if isinstance(k, (list, tuple)) else 0 for k in x])

    print(nodes)

    # Compute transition frequencies
    transitions = [(tuple(abstractions[i]), tuple(abstractions[i + 1])) for i in range(len(abstractions) - 1)]
    transitions_freq = pd.Series(transitions).value_counts()

    # Create the graph
    G = nx.DiGraph()

    # Add nodes with attributes
    for node, freq in nodes.items():
        node_times = times[abstractions.apply(tuple) == node]
        weighted_time = (node_times * freq).mean()
        normalized_time = (weighted_time - times.min()) / (times.max() - times.min())
        color = to_hex(plt.cm.Blues(normalized_time))  # Convert the color to hexadecimal
        G.add_node(node, size=freq, color=color)

    # Add edges with weights
    for (source, target), weight in transitions_freq.items():
        G.add_edge(source, target, weight=weight)

    # Hierarchical layout
    pos = nx.drawing.nx_agraph.graphviz_layout(G, prog="dot")
    pos = {node: (x, -y) for node, (x, y) in pos.items()}

    # Create the figure and axis for animation
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.axis("off")

    # Sizes and colors for the last three nodes
    highlight_sizes = [1 / 4, 1 / 2, 1]
    highlight_colors = ['lightblue', 'dodgerblue', 'blue']  # Different shades of blue

    # Create a function for drawing each frame
    def update_frame(i):
        ax.clear()  # Clear the current frame

        # Draw the nodes
        current_abstraction = abstractions.iloc[i]
        visited_nodes = abstractions.iloc[:i + 1].apply(tuple).value_counts().index
        node_sizes = [G.nodes[node]["size"] * 100 for node in G.nodes]

        # Adjust node colors based on visited state, with transparency for past nodes
        node_colors = [
            G.nodes[node]["color"] if node not in visited_nodes else (0.8, 0.8, 1, 0.3)  # Fade visited nodes
            for node in G.nodes
        ]

        # Highlight the last three unique nodes
        for idx, last_node in enumerate(last_three_unique):
            if last_node in G.nodes:
                node_index = list(G.nodes).index(last_node)
                node_sizes[node_index] = highlight_sizes[idx] * 500  # Update size
                node_colors[node_index] = highlight_colors[idx]  # Update color

        # Draw the nodes
        nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors, ax=ax)

        # Draw the edges with width based on transition frequency
        edge_weights = [G.edges[edge]["weight"] for edge in G.edges]
        nx.draw_networkx_edges(G, pos, width=[w / 2 for w in edge_weights], alpha=0.7, ax=ax)

        # Draw the current node as a larger, fully opaque red node
        current_node = tuple(current_abstraction)
        nx.draw_networkx_nodes(G, pos, nodelist=[current_node], node_size=500, node_color='red', ax=ax)

        # Add labels for nodes
        node_labels = {node: f"({', '.join(map(str, node))})" for node in G.nodes}
        nx.draw_networkx_labels(G, pos, labels=node_labels, font_size=8, font_color="black", ax=ax)

        # Title
        ax.set_title(f"{title} {i}", fontsize=15)

        return ax,

    # Create the animation by updating the frame for each timestep
    ani = animation.FuncAnimation(fig, update_frame, frames=len(abstractions), interval=interval, blit=False, repeat=False)

    # Show the animation
    plt.show()

##############################################################################################

def plot_abstraction_graph_movie_3last_nodes(abstract_time_series, 
                                             first_node_size=125, opacity_first_node=0.5,
                                             last_node_size=500, last_adjust_sizes = [1/4, 1/2, 1], 
                                             colour_last_nodes = ['blue', 'blue', 'blue'], 
                                             interval=400, 
                                             title="Abstraction Graph - Time"):
    """
    Creates an animated abstraction graph showing the evolution of abstraction sets over time.

    Parameters:
    ----------
    abstract_time_series : pd.DataFrame
        A DataFrame with the following columns:
        - 'Time': Numeric, indicating the occurrence time of each abstraction.
        - 'Abstraction': Iterable, representing abstraction sets.
    
    first_node_size : int, optional
        Size of the nodes that are not part of the last three unique abstractions. Defaults to 125.
    
    opacity_first_node : float, optional
        Opacity of nodes representing visited states in the animation. Range [0, 1]. Defaults to 0.5.
    
    last_node_size : int, optional
        Base size for the last three unique nodes. Their sizes are scaled by `last_adjust_sizes`. Defaults to 500.
    
    last_adjust_sizes : list of float, optional
        Scaling factors for the sizes of the last three unique nodes. Defaults to [1/4, 1/2, 1].
    
    colour_last_nodes : list of str, optional
        Colors for the last three unique nodes. Should contain three color specifications. Defaults to ['blue', 'blue', 'blue'].
    
    interval : int, optional
        Interval in milliseconds between frames of the animation. Defaults to 400.
    
    title : str, optional
        Title of the graph animation. Defaults to "Abstraction Graph - Time".

    Returns:
    -------
    None
        Displays an animation of the abstraction graph.
    """

    abstractions = abstract_time_series["Abstraction"]
    abstractions= pd.concat([abstractions, pd.Series([abstractions.iloc[-1]])], ignore_index=True)
    times = abstract_time_series["Time"]
    unique_abstractions = abstractions.apply(tuple).drop_duplicates(keep='last')  # Keep last occurrences of unique abstractions
    last_three_unique = unique_abstractions.tail(3)  # Last three unique sets
    nodes = abstractions.apply(tuple).value_counts(sort=False) # Count node frequencies and without sort (sort=False).                 # nodes = nodes.sort_index(key=lambda x: [len(k) if isinstance(k, (list, tuple)) else 0 for k in x])
    
    # Create unique labels for nodes (S1, S2, S3, ...)
    node_labels = {node: f"S{idx + 1}" for idx, node in enumerate(nodes.keys())}
    # print(node_labels)

    # Compute transition frequencies
    transitions = [(tuple(abstractions[i]), tuple(abstractions[i + 1])) for i in range(len(abstractions) - 1)]
    transitions_freq = pd.Series(transitions).value_counts()

    # Create the graph
    G = nx.DiGraph()

    # Add nodes with attributes
    for node, freq in nodes.items():
        nodee = ', '.join(node)                    # Delete inverted commas from the elements of the vector nodes
        node_times = times[abstractions.apply(tuple) == node]
        weighted_time = (node_times * freq).mean()
        normalized_time = (weighted_time - times.min()) / (times.max() - times.min())
        color = to_hex(plt.cm.Blues(normalized_time))       # Convert the color to hexadecimal
        label = f"{node_labels[node]} ({freq}) \n[{nodee}]" # Node labels: S1 (2) [a, b,c]
        G.add_node(node, size=freq, color=color, label=label)

    # Add edges with weights
    for (source, target), weight in transitions_freq.items():
        G.add_edge(source, target, weight=weight)

    # Hierarchical layout
    pos = nx.drawing.nx_agraph.graphviz_layout(G, prog="dot")
    pos = {node: (x, -y) for node, (x, y) in pos.items()}

    # Create the figure and axis for animation
    fig, ax = plt.subplots(figsize=(12, 9))
    ax.axis("off")

    # Create a function for drawing each frame
    def update_frame(i):
        ax.clear()     # Limpia el eje para la nueva iteración
        ax.axis("off") # Desactiva el borde y las marcas del eje

        # Draw the nodes
        current_abstraction = abstractions.iloc[i]
        visited_nodes = abstractions.iloc[:i + 1].apply(tuple).value_counts().index
        node_sizes = [G.nodes[node]["size"] * first_node_size for node in G.nodes]

        # Adjust node colors based on visited state, with transparency for past nodes
        node_colors = [
            G.nodes[node]["color"] if node not in visited_nodes else (0.8, 0.8, 1, opacity_first_node)  # Fade visited nodes. Azul claro:0.8, 0.8, 1. Opacidad: 0.5
            for node in G.nodes
        ]

        for idx, last_node in enumerate(last_three_unique):
            if last_node in G.nodes:
                node_index = list(G.nodes).index(last_node)
                node_sizes[node_index] = last_adjust_sizes[idx] * last_node_size  # Adjust size
                node_colors[node_index] = colour_last_nodes[idx]      # Uniform blue color

        # Draw the nodes
        nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors, ax=ax)

        # Draw the edges with width based on transition frequency
        edge_weights = [G.edges[edge]["weight"] for edge in G.edges]
        nx.draw_networkx_edges(G, pos, width=[w / 2 for w in edge_weights], alpha=0.9, ax=ax)

        # Draw the current node as a larger, fully opaque red node, except for the last frame
        if i < len(abstractions)-1:  # Evita el último frame
            current_node = tuple(current_abstraction)
            nx.draw_networkx_nodes(G, pos, nodelist=[current_node], node_size=500, node_color='red', ax=ax)
        
        # Add labels for nodes 
        label_pos = {node: (x, y + 10) for node, (x, y) in pos.items()}  # Ajusta '10' según el desplazamiento deseado        
        nx.draw_networkx_labels(
            G, label_pos, labels=nx.get_node_attributes(G, 'label'), 
            font_size=8, font_color="black", ax=ax, 
            bbox=dict(facecolor="white", edgecolor="white", boxstyle="round,pad=0.3", alpha=0)
        )
        
        # Title
        if i < len(abstractions) - 1:  # Para todos los frames excepto el último
            ax.set_title(f"{title} {i}", fontsize=15)
        else:  # Para el último frame
            j=len(abstractions)-2
            ax.set_title(f"{title} {j}", fontsize=15)

        return ax,

    # Create the animation by updating the frame for each timestep
    ani = animation.FuncAnimation(fig, update_frame, frames=len(abstractions), interval=interval, blit=False, repeat=False)

    # Show the animation 
    plt.show()

##############################################################################################    

def get_plot_abstraction_graph_movie_html(abstract_time_series, 
                                          first_node_size=25, colour_first_nodes="lightcyan",
                                          last_node_size=35, last_adjust_sizes=[1/2, 3/4, 1], 
                                          colour_last_nodes="red", 
                                          interval=1000, showlegend=True,
                                          title="Abstraction Graph Movie",
                                          filename="abstraction_graph_movie.html"):
    """
    Generates an animated HTML file visualizing the abstraction graph over time.

    Parameters:
    ----------
    abstract_time_series : pd.DataFrame
        A DataFrame containing two columns:
        - "Time": Corresponding timestamps for each abstraction state.        
        - "Abstraction": A list of abstraction states over time.
    first_node_size : int, optional. The base size for nodes representing the abstraction states, by default 25.
    colour_first_nodes : str, optional. Color for the nodes, by default "lightcyan".
    last_node_size : int, optional. The size of the most recent nodes in the animation, by default 35.
    last_adjust_sizes : list of float, optional. Adjustment multipliers for the size of the last three unique nodes, by default [1/4, 1/2, 1].
    colour_last_nodes : list of str, optional. Colors for the last three unique nodes, by default ['blue', 'blue', 'blue'].
    interval : int, optional. Time interval between frames in milliseconds, by default 1000.
    showlegend : bool, optional. Whether to display the legend in the animation, by default True.
    title : str, optional. Title of the animation, by default "Abstraction Graph Movie".
    filename : str, optional. Name of the output HTML file, by default "abstraction_graph_movie.html".

    Returns:
    -------
    str
        The filename of the generated HTML file.
    """
    abstractions = abstract_time_series["Abstraction"]                           # Obtiene la columna "Abstraction" del DataFrame 'abstract_time_series'.
    abstractions = pd.concat([abstractions, pd.Series([abstractions.iloc[-1]])], ignore_index=True)  # Añade el último valor de la columna al final, para garantizar que se considere la transición final.
    times = abstract_time_series["Time"]                                         # Obtiene la columna "Time" del DataFrame 'abstract_time_series'.
    unique_abstractions = abstractions.apply(tuple).drop_duplicates(keep='last') # Convierte los elementos de 'abstractions' a tuplas y elimina duplicados, manteniendo el último.
    last_three_unique = unique_abstractions.tail(3)                              # Selecciona las tres últimas abstracciones únicas.
    nodes = abstractions.apply(tuple).value_counts(sort=False)                   # Cuenta las ocurrencias de cada tupla única en 'abstractions'.

    # Create unique labels for nodes
    node_labels = {node: f"S{idx + 1}" for idx, node in enumerate(nodes.keys())}  # Crea etiquetas únicas para los nodos, como "S1", "S2", etc.

    # Compute transition frequencies
    transitions = [(tuple(abstractions[i]), tuple(abstractions[i + 1])) for i in range(len(abstractions) - 1)]  # Crea las transiciones entre los elementos consecutivos de 'abstractions'.
    transitions_freq = pd.Series(transitions).value_counts()  # Cuenta la frecuencia de cada transición.

    # Create the graph
    G = nx.DiGraph()  # Crea un grafo dirigido vacío.

    for node, freq in nodes.items():                          # Itera sobre los nodos y sus frecuencias.
        nodee = ', '.join(node)                               # Convierte la tupla del nodo en una cadena de texto.
        color = colour_first_nodes                            # Define el color de los nodos iniciales.
        label = f"{node_labels[node]} ({freq})"               # Etiqueta del nodo con su frecuencia.
        hover_info = f"[{nodee}]"                             # Información detallada para el hover.
        G.add_node(node, size=1, color=color, label=label, hovertext=hover_info)  # Añade el nodo al grafo con atributos adicionales.

    for (source, target), weight in transitions_freq.items(): # Itera sobre las transiciones y sus frecuencias.
        G.add_edge(source, target, weight=weight)             # Añade las transiciones como aristas con su peso (frecuencia).

    pos = nx.drawing.nx_agraph.graphviz_layout(G, prog="dot") # Calcula las posiciones de los nodos utilizando Graphviz.
    pos = {node: (x, -y) for node, (x, y) in pos.items()}     # Invierte la coordenada y para ajustar la visualización.

    frames = []                                               # Lista para almacenar los frames de la animación.
    sliders_steps = []                                        # Lista para los pasos del control deslizante de la animación.
    
    for i in range(len(abstractions)-1):                      # Itera sobre las abstracciones (menos el último valor).
        node_sizes = [G.nodes[node]["size"] * first_node_size for node in G.nodes]  # Define el tamaño de los nodos.
        node_colors = [G.nodes[node]["color"] for node in G.nodes]                  # Define el color de los nodos.

        # # Destacar los nodos actuales y ajustar tamaños progresivamente
        # for offset, size_factor in enumerate(last_adjust_sizes[::-1]):  # Ajustes de tamaño
        #     # Calcular el índice máximo de los nodos a graficar
        #     max_nodes_to_plot = i + 1  # Número de nodos a graficar en función del tiempo `i`
            
        #     # Verificar que el nodo actual esté dentro del rango permitido
        #     if offset < max_nodes_to_plot:
        #         current_node = tuple(abstractions.iloc[i - offset])  # Obtener el nodo actual
                
        #         # Verificar si el nodo existe en el grafo
        #         if current_node in G.nodes:
        #             index = list(G.nodes).index(current_node)        # Índice del nodo
        #             node_sizes[index] = last_node_size * size_factor # Ajustar tamaño
        #             node_colors[index] = colour_last_nodes           # Cambiar color del nodo
        #         else:
        #             continue  # Saltar si el nodo no existe en el grafo
        # # Destacar los nodos actuales y ajustar tamaños progresivamente
        # for offset, size_factor in enumerate(last_adjust_sizes[::-1]):  # Ajustes de tamaño
        #     # Calcular el índice máximo de los nodos a graficar
        #     max_nodes_to_plot = i + 1  # Número de nodos a graficar en función del tiempo `i`
            
        #     # Verificar que el nodo actual esté dentro del rango permitido
        #     if offset < max_nodes_to_plot:
        #         current_node = tuple(abstractions.iloc[i - offset])  # Obtener el nodo actual
                
        #         # Verificar si el nodo existe en el grafo
        #         if current_node in G.nodes:
        #             index = list(G.nodes).index(current_node)  # Índice del nodo
                    
        #             # Encontrar el tamaño máximo de nodos repetidos
        #             max_size = max(node_sizes[index], last_node_size * size_factor)
        #             node_sizes[index] = max_size  # Ajustar tamaño al máximo
        #             node_colors[index] = colour_last_nodes  # Cambiar color del nodo
        #         else:
        #             continue  # Saltar si el nodo no existe en el grafo
        # Destacar los nodos actuales y ajustar tamaños progresivamente
        for offset, size_factor in enumerate(last_adjust_sizes[::-1]):  # Ajustes de tamaño
            # Calcular el índice máximo de los nodos a graficar
            max_nodes_to_plot = i + 1  # Número de nodos a graficar en función del tiempo `i`
            
            # Verificar que el nodo actual esté dentro del rango permitido
            if offset < max_nodes_to_plot:
                current_node = tuple(abstractions.iloc[i - offset])  # Obtener el nodo actual
                
                # Verificar si el nodo existe en el grafo
                if current_node in G.nodes:
                    index = list(G.nodes).index(current_node)  # Índice del nodo
                    
                    # Condición para tres nodos consecutivos diferentes
                    if offset >= 2:  # Verificar si hay al menos tres nodos consecutivos
                        prev_node_1 = tuple(abstractions.iloc[i - offset + 1])
                        prev_node_2 = tuple(abstractions.iloc[i - offset + 2])
                        if (
                            prev_node_1 in G.nodes
                            and prev_node_2 in G.nodes
                            and current_node != prev_node_1
                            and current_node != prev_node_2
                            and prev_node_1 != prev_node_2
                        ):
                            node_sizes[index] = last_node_size * size_factor
                        else:
                            # Ajustar tamaño al máximo existente
                            max_size = max(node_sizes[index], last_node_size * size_factor)
                            node_sizes[index] = max_size
                    else:
                        # Ajustar tamaño al máximo existente
                        max_size = max(node_sizes[index], last_node_size * size_factor)
                        node_sizes[index] = max_size
                    
                    node_colors[index] = colour_last_nodes  # Cambiar color del nodo
                else:
                    continue  # Saltar si el nodo no existe en el grafo

 
        edge_x = []                    # Lista para las coordenadas x de las aristas.
        edge_y = []                    # Lista para las coordenadas y de las aristas.
        annotations = []               # Lista para las anotaciones de las aristas.

        for source, target in G.edges: # Itera sobre las aristas.
            x0, y0 = pos[source]       # Obtiene las coordenadas del nodo fuente.
            x1, y1 = pos[target]       # Obtiene las coordenadas del nodo objetivo.
            edge_x += [x0, x1, None]   # Añade las coordenadas de la arista.
            edge_y += [y0, y1, None]   # Añade las coordenadas de la arista.

            annotations.append(dict(   # Añade anotaciones para cada arista.
                ax=x0, ay=y0,
                x=x1, y=y1,
                xref="x", yref="y",
                axref="x", ayref="y",
                showarrow=True,
                arrowhead=3,
                arrowsize=2,
                arrowwidth=1,
                arrowcolor="gray"
            ))

        frames.append(go.Frame(                                  # Añade un nuevo frame para la animación.
            data=[
                go.Scatter(x=[pos[node][0] for node in G.nodes], # Posiciones x de los nodos.
                        y=[pos[node][1] for node in G.nodes],    # Posiciones y de los nodos.
                        mode='markers+text',                     # Modo de visualización de los nodos.
                        marker=dict(size=node_sizes, color=node_colors, opacity=1), # Tamaño y color de los nodos.
                        text=[G.nodes[node]["label"] for node in G.nodes],          # Etiquetas de los nodos.
                        textposition="bottom center"),           # Posición del texto.
                go.Scatter(x=edge_x, y=edge_y,                   # Coordenadas de las aristas.
                        mode='lines',                            # Modo de visualización de las aristas.
                        line=dict(width=1, color='gray'))        # Estilo de las aristas.
            ],
            layout=dict(annotations=annotations),                # Añade las anotaciones a la visualización.
            name=f"Frame {i}"                                    # Nombre del frame para la animación.
        ))
        sliders_steps.append(dict(                               # Añade un paso al control deslizante de la animación.
            args=[[f"Frame {i}"], dict(frame=dict(duration=0, redraw=True), mode="immediate")],  # Define el paso.
            label=f"{i}",                                        # Etiqueta del paso.
            method="animate"                                     # Define el método de animación.
        ))

    legend_text = "<br>".join([f"{node_labels[node]}= [{', '.join(node)}]" for node in G.nodes])  # Genera el texto de la leyenda con etiquetas de nodos (S1, S2, ...) y sus elementos.

    fig = go.Figure(
        data=[
            go.Scatter(
                x=[pos[node][0] for node in G.nodes],  # Coordenadas X de los nodos.
                y=[pos[node][1] for node in G.nodes],  # Coordenadas Y de los nodos.
                mode='markers+text',                   # Modo de visualización: muestra puntos (markers) y texto (text).
                marker=dict(size=node_sizes, color=node_colors, opacity=1),  # Configuración de los nodos: tamaño, color, y opacidad.
                text=[G.nodes[node]["label"] for node in G.nodes],  # Etiquetas que se muestran en cada nodo.
                textposition="bottom center",          # Posición del texto respecto al nodo (debajo y centrado).
                hovertext=[G.nodes[node]["hovertext"] for node in G.nodes],  # Texto que aparece al pasar el mouse sobre un nodo.
                hoverinfo="text",                      # Especifica que se muestra el texto de hover al interactuar.
                name=legend_text                       # Asocia el texto de la leyenda con el gráfico.
            )
        ],
        layout=go.Layout(
            title=title,                      # Título del gráfico.
            updatemenus=[                     # Configuración de los botones de control (play y pause).
                dict(
                    type="buttons",           # Define un grupo de botones.
                    showactive=False,         # Desactiva el resaltado del botón seleccionado.
                    buttons=[
                        dict(
                            label="Play",     # Etiqueta del botón "Play".
                            method="animate", # Método que activa la animación.
                            args=[None, dict(frame=dict(duration=interval, redraw=True), fromcurrent=True)]  # Configuración de la animación.
                        ),
                        dict(
                            label="Pause",    # Etiqueta del botón "Pause".
                            method="animate", # Método que pausa la animación.
                            args=[[None], dict(frame=dict(duration=0, redraw=False), mode="immediate")]  # Configuración de pausa.
                        )
                    ]
                )
            ],
            sliders=[{  # Configuración del control deslizante para avanzar manualmente entre los cuadros de la animación.
                'steps': sliders_steps,   # Pasos del deslizador.
                'currentvalue': {
                    'prefix': "Time: ",   # Prefijo que aparece junto al valor actual del deslizador.
                    'font': {'size': 16}, # Tamaño de fuente del prefijo.
                    'visible': True,      # Muestra el valor actual.
                },
                'x': 0.1,   # Posición horizontal del deslizador.
                'len': 0.9, # Longitud del deslizador en proporción al gráfico.
            }],
            xaxis=dict(showgrid=False, zeroline=False, visible=False),  # Oculta la cuadrícula y el eje X.
            yaxis=dict(showgrid=False, zeroline=False, visible=False),  # Oculta la cuadrícula y el eje Y.
            plot_bgcolor='rgba(0,0,0,0)', # Fondo transparente del gráfico.
            width=1400,                   # Ancho del gráfico.
            height=800,                   # Altura del gráfico.
            annotations=annotations,      # Anotaciones adicionales para destacar elementos.
            showlegend=showlegend,        # Habilita la visualización de la leyenda.
            legend=dict( 
                title=dict(text="Node Legend", font=dict(family="Arial", size=14, color="darkblue")),  # Título de la leyenda.
                bgcolor="white",  # Fondo blanco.
                bordercolor="white",                               # Borde blanco.
                borderwidth=2,                                     # Grosor del borde.
                font=dict(family="Arial", size=12, color="black"), # Fuente de los elementos.
                x=1,                                               # Posición horizontal a la derecha.
                y=1.0                                              # Posición vertical en la parte superior.
            )
        ),
        frames=frames  # Cuadros de la animación.
    )

    fig.write_html(filename)
    return filename


############################################################################################## 
 
def plot_abstraction_graph_movie_html(
                abstract_time_series, 
                first_node_size=25, colour_first_nodes="lightcyan",
                last_node_size=35, last_adjust_sizes=[1/2, 3/4, 1], 
                colour_last_nodes="red", 
                interval=1000, showlegend=True,
                title="Abstraction Graph Movie",
                filename="abstraction_graph_movie.html"):
    """
    Generates an interactive animation of an abstraction graph from an abstract time series, 
    saving the animation as an HTML file. The animation shows how the graph evolves over time, 
    highlighting the most important nodes.

    Parameters:
    ----------
    abstract_time_series : pd.DataFrame
        A DataFrame containing two columns:
        - "Abstraction": A list of abstraction states over time.
        - "Time": Corresponding timestamps for each abstraction state.
    first_node_size : int, optional. The base size for nodes representing the abstraction states, by default 25.
    colour_first_nodes : str, optional. Color for the nodes, by default "lightcyan".
    last_node_size : int, optional. The size of the most recent nodes in the animation, by default 35.
    last_adjust_sizes : list of float, optional. Adjustment multipliers for the size of the last three unique nodes, by default [1/4, 1/2, 1].
    colour_last_nodes : list of str, optional. Colors for the last three unique nodes, by default ['blue', 'blue', 'blue'].
    interval : int, optional. Time interval between frames in milliseconds, by default 1000.
    showlegend : bool, optional. Whether to display the legend in the animation, by default True.
    title : str, optional. Title of the animation, by default "Abstraction Graph Movie".
    filename : str, optional. Name of the output HTML file, by default "abstraction_graph_movie.html".

    Returns:
    -------
    str
        The filename of the generated HTML file.
    Example:
    -------
    # Load an abstraction time series
    abstract_time_series = pd.DataFrame({
        'Time': [0, 1, 2],
        'Abstraction': [['A', 'B'], ['B'], ['C', 'D']]        
    })

    # Call the function to generate the animation
    plot_abstraction_graph_movie_html(abstract_time_series, filename="graph_movie.html")
    """    
    # Call the `rn_get_visualization` function to generate the HTML file
    get_plot_abstraction_graph_movie_html(
        abstract_time_series,
        first_node_size=first_node_size, colour_first_nodes=colour_first_nodes,
        last_node_size=last_node_size, last_adjust_sizes=last_adjust_sizes, 
        colour_last_nodes=colour_last_nodes, 
        interval=interval, showlegend=showlegend,
        title=title,
        filename=filename
    )
    
    # Convert to an absolute path
    abs_path = os.path.abspath(filename) 
    # print("absolute path")
    # print(abs_path)
    # Check if the file was created correctly
    if not os.path.isfile(abs_path):
        print(f"File not found at {abs_path}")  # Additional message for debugging
        raise FileNotFoundError(f"The file {abs_path} was not found. Check if `rn_get_visualization` generated the file correctly.")
    
    # Inform the user about the file's location
    print(f"The animation was saved to: {abs_path}")
    
    # Open the HTML file in the default browser
    webbrowser.open(f"file://{abs_path}")

##############################################################################################    













##############################################################################################

def plot_abstraction_graph_movie_html0(abstract_time_series, interval=400, title="Abstraction Graph - Time", 
                                      filename="abstraction_graph_movie0.html"):
    """
    Creates an interactive abstraction graph that evolves over time and saves it as an HTML file.

    Parameters:
    abstract_time_series (pd.DataFrame): A DataFrame with the following columns:
        - 'Time': Time of occurrence for each abstraction (numeric).
        - 'Abstraction': Abstraction sets represented as iterables.
    interval (int): Interval in milliseconds between frames of the animation. Defaults to 400.
    title (str): Title of the graph animation. Defaults to "Abstraction Graph - Time".
    filename (str): The name of the HTML file to save the animation.

    Returns:
    None. Saves the interactive abstraction graph animation to an HTML file.
    """
    abstractions = abstract_time_series["Abstraction"]
    times = abstract_time_series["Time"]
    unique_abstractions = abstractions.apply(tuple).drop_duplicates(keep='last')  # Keep last occurrences of unique abstractions
    last_three_unique = unique_abstractions.tail(3)  # Last three unique sets

    nodes = abstractions.apply(tuple).value_counts()

    # Compute transition frequencies
    transitions = [(tuple(abstractions[i]), tuple(abstractions[i + 1])) for i in range(len(abstractions) - 1)]
    transitions_freq = pd.Series(transitions).value_counts()

    # Create the graph
    G = nx.DiGraph()

    # Add nodes with attributes
    for node, freq in nodes.items():
        node_times = times[abstractions.apply(tuple) == node]
        weighted_time = (node_times * freq).mean()
        normalized_time = (weighted_time - times.min()) / (times.max() - times.min())
        color = to_hex(plt.cm.Blues(normalized_time))  # Convert the color to hexadecimal
        G.add_node(node, size=freq, color=color)

    # Add edges with weights
    for (source, target), weight in transitions_freq.items():
        G.add_edge(source, target, weight=weight)

    # Hierarchical layout
    pos = nx.drawing.nx_agraph.graphviz_layout(G, prog="dot")
    pos = {node: (x, -y) for node, (x, y) in pos.items()}

    # Prepare the PyVis network
    net = Network(height="750px", width="100%", directed=True)

    # Add nodes to the PyVis network
    for node in G.nodes:
        node_size = G.nodes[node]["size"] * 100
        color = G.nodes[node]["color"]
        
        # Add node to the network (convert to string for PyVis compatibility)
        net.add_node(str(node), label=f"({', '.join(map(str, node))})", title=f"({', '.join(map(str, node))})", color=color, size=node_size, x=pos[node][0], y=pos[node][1])

    # Add edges to the PyVis network
    for source, target in G.edges:
        weight = G.edges[source, target].get("weight", 1)
        # Add edges to network (convert nodes to string)
        net.add_edge(str(source), str(target), width=weight)

    # Set up the layout options for the network
    net.set_options("""
    {
        "nodes": {
            "shape": "dot",
            "font": {
                "size": 14,
                "align": "center"
            },
            "borderWidth": 2,
            "borderColor": "black"
        },
        "edges": {
            "smooth": false,
            "color": "#888888",
            "width": 1
        },
        "physics": {
            "enabled": false
        },
        "layout": {
            "hierarchical": {
                "enabled": true,
                "direction": "DU",
                "sortMethod": "directed"
            }
        }
    }
    """) 
        # Save the visualization to an HTML file 
    net.html = net.generate_html()
    with open(filename, "w", encoding="utf-8") as f:
        f.write(net.html)
    return filename

##############################################################################################
# # Películas animadas de gráficos de abstracción y semiorganizaciones
##############################################################################################

# # Function to generate a Hasse diagram of the set of species in a reaction network.
def plot_hasse_diagram(species_set, node_size=25, color_node='cyan', arrowcolor='black', bgcolor_legend="white", bordercolor_legend="white", name_legend="Legend:", title="Hasse diagram"):
    """
    Generates and plots a Hasse diagram for a given set of species. The diagram includes node labels 
    and a legend, showing the relationships between different subsets of species. 

    Parameters:
    ----------
    species_set : set
        A set of species that will be used to generate the Hasse diagram.
    node_size : int, optional
        The size of the nodes in the diagram, by default 25.
    color_node : str, optional
        The color for the nodes in the diagram, by default 'cyan'.
    arrowcolor : str, optional
        The color of the arrows connecting the nodes, by default 'black'.
    bgcolor_legend : str, optional
        The background color of the legend, by default 'white'.
    bordercolor_legend : str, optional
        The border color of the legend, by default 'white'.
    name_legend : str, optional
        The title for the legend, by default 'Legend:'.
    title : str, optional
        The title of the diagram, by default 'Hasse diagram'.

    Returns:
    -------
    None
        The function displays the Hasse diagram plot using Plotly.
    
    Example:
    -------
    species_set = {"A", "B", "C", "D"}
    plot_hasse_diagram(species_set)
    """
    # Generar todos los subconjuntos del conjunto dado
    elementos = list(species_set)
    subconjuntos = [set(comb) for i in range(len(elementos) + 1) for comb in combinations(elementos, i)]
    
    # Crear nodos y relaciones
    nodos = []
    edges = []
    niveles = {}

    for subconjunto in subconjuntos:
        nodo = frozenset(subconjunto)
        nodos.append(nodo)
        niveles[nodo] = len(subconjunto)
    
    for subconjunto1 in subconjuntos:
        for subconjunto2 in subconjuntos:
            if subconjunto1 < subconjunto2 and len(subconjunto2 - subconjunto1) == 1:
                edges.append((frozenset(subconjunto1), frozenset(subconjunto2)))

    # Etiquetar nodos
    etiquetas = {nodo: f"S{idx + 1}" for idx, nodo in enumerate(nodos)}
    
    # Posiciones de los nodos (usando el nivel como coordenada y distribuyendo en X)
    posiciones = {}
    nivel_max = max(niveles.values())
    for nivel in range(nivel_max + 1):
        nodos_nivel = [nodo for nodo in nodos if niveles[nodo] == nivel]
        x_coords = range(-len(nodos_nivel) // 2, len(nodos_nivel) // 2 + 1)
        for x, nodo in zip(x_coords, nodos_nivel):
            posiciones[nodo] = (x, nivel)

    # Datos para Plotly
    x_coords = []
    y_coords = []
    labels = []
    hovertexts = []
    for nodo, (x, y) in posiciones.items():
        x_coords.append(x)
        y_coords.append(y)  # Sin invertir Y, el vacío queda abajo
        labels.append(etiquetas[nodo])
        hovertexts.append(f"{[', '.join(nodo) if nodo else '∅']}")

    # Crear las aristas con flechas
    annotations = []
    for nodo1, nodo2 in edges:
        x0, y0 = posiciones[nodo1]
        x1, y1 = posiciones[nodo2]
        annotations.append(dict(
            ax=x0, ay=y0,
            x=x1, y=y1,
            xref='x', yref='y',
            axref='x', ayref='y',
            showarrow=True,
            arrowhead=2,
            arrowsize=1,
            arrowwidth=2,
            arrowcolor=arrowcolor
        ))

    # Trazar los nodos
    node_trace = go.Scatter(
        x=x_coords, y=y_coords,
        mode='markers+text',
        text=labels,
        hovertext=hovertexts,
        textposition="top center",
        marker=dict(
            size=node_size,
            color=color_node,
            line=dict(width=2, color=color_node)
        )
    )

    # Crear la leyenda
    leyenda_texto = "<br>".join([f"{etiquetas[nodo]}= [{', '.join(nodo) if nodo else '∅'}]" for nodo in nodos])

    # Crear la figura
    fig = go.Figure(data=[node_trace])
    fig.update_layout(
        title=title,
        showlegend=False,
        xaxis=dict(
            showgrid=False, zeroline=False, showticklabels=False
        ),
        yaxis=dict(
            showgrid=False, zeroline=False, showticklabels=False
        ),
        plot_bgcolor='white',
        annotations=annotations
    )
    
    # Agregar la leyenda como anotación
    fig.add_annotation(
        x=1.05, y=1, xref="paper", yref="paper",
        text=f"<b>{name_legend}</b><br>{leyenda_texto}",
        showarrow=False,
        align="left",
        bordercolor=bordercolor_legend,
        borderwidth=1,
        bgcolor=bgcolor_legend, 
        font=dict(size=12)
    )

    fig.show()

##############################################################################################
# # Function to generate an interactive graph from a Hasse diagram of Semi-Organisations and Abstractions.
def get_film_semiorganizations_abstractions_html(abstract_time_series, input_data,
                                          first_node_size=25, 
                                          colour_abs_nodes="lightcyan", colour_semiorg_nodes="lightgreen", colour_cap_nodes="orange",
                                          last_node_size=35, last_adjust_sizes=[1/2, 3/4, 1], colour_last_nodes="red", 
                                          interval=100, showlegend=True,
                                          title="Semi-organisations and Abstractions Film",
                                          filename="semiorg_abstraction_film.html"):
    """
    Generates an animated HTML file visualizing the abstraction graph over time.

    Parameters:
    ----------
    abstract_time_series : pd.DataFrame
        A DataFrame containing two columns:
        - "Time": Corresponding timestamps for each abstraction state.        
        - "Abstraction": A list of abstraction states over time.
    input_data : list of lists
        A list of subsets to visualize as additional nodes and intersections.
    first_node_size : int, optional. The base size for nodes representing the abstraction states, by default 25.
    colour_abs_nodes : str, optional. Color for the nodes, by default "lightcyan".
    colour_semiorg_nodes : str, optional. Color for semi-organization nodes, by default "lightgreen".
    colour_cap_nodes : str, optional. Color for the intersection nodes, by default "orange".
    last_node_size : int, optional. The size of the most recent nodes in the animation, by default 35.
    last_adjust_sizes : list of float, optional. Adjustment multipliers for the size of the last three unique nodes, by default [1/4, 1/2, 1].
    colour_last_nodes : str, optional. Colors for the last three unique nodes, by default ['red', 'red', 'red'].
    interval : int, optional. Time interval between frames in milliseconds, by default 100.
    showlegend : bool, optional. Whether to display the legend in the animation, by default True.
    title : str, optional. Title of the animation, by default "Semi-organisations and Abstractions Film".
    filename : str, optional. Name of the output HTML file, by default "semiorg_abstraction_film.html".

    Returns:
    -------
    str
        The filename of the generated HTML file.
    """
    # Convert input_data to unique sets
    unique_subsets = []
    for sublist in input_data:
        if set(sublist) not in [set(x) for x in unique_subsets]:
            unique_subsets.append(sublist) 

    # Create a list of unique sets for reference
    set_of_sets = [tuple(subset) for subset in unique_subsets] 
    # Ordenar set_of_sets por la cantidad de elementos en cada subconjunto
    set_of_sets = sorted([tuple(sorted(subset)) for subset in set_of_sets], key=len) 

    set_of_sets2 = pd.Series(set_of_sets) # Convertir a Series para aplicar .apply() y .value_counts()
    set_of_sets2 = set_of_sets2.apply(tuple).value_counts(sort=False) # Luego puedes aplicar .apply() y .value_counts()
    # print("set_of_sets2")
    # print(set_of_sets2)
    # Create set names based on their level
    Set_names = {node: f"S{idx + 1}" for idx, node in enumerate(set_of_sets2.keys())}
    # Set_names = [f"S{i+1}" for i in range(len(set_of_sets2))]

    # Create set labels to display when hovering over nodes
    abstractions_0 = abstract_time_series["Abstraction"] # Abstraction states over time # abstractions = pd.concat([abstractions, pd.Series([abstractions.iloc[-1]])], ignore_index=True) # Add the last state to the end
    abstractions_0 = abstractions_0.apply(sorted)
    # print("abstractions_0")
    # print(abstractions_0)
    
    abstractions = pd.concat([abstractions_0, pd.Series(set_of_sets)], ignore_index=True)
    # unique_abstractions = abstractions.apply(tuple).drop_duplicates(keep='last') # Unique abstractions times = abstract_time_series["Time"] # Time values for each abstraction state # last_three_unique = unique_abstractions.tail(3) # Get the last three unique abstractions

    # Ordenar los elementos de cada subconjunto en unique_subsets
    unique_subsets = [sorted(subset) for subset in unique_subsets]
    # print()
    # print("unique_abstractions", unique_abstractions)

    nodes = abstractions.apply(tuple).value_counts(sort=False) 
    nodes = nodes[~nodes.index.duplicated()]

    nodes_0 = abstractions_0.apply(tuple).value_counts(sort=False)
    nodes_1 = nodes.loc[~nodes.index.isin(nodes_0.index)] 

    # print("nodes")
    # print(nodes)  
    # print("nodes_0")  
    # print(nodes_0)
    # print("nodes_1")    
    # print(nodes_1)     

    # Create unique labels for nodes
    # node_labels = {node: f"S{idx + 1}" for idx, node in enumerate(nodes.keys())}
    
    # Crear node_labels para nodes_0 y nodes_1
    node_labels_0 = {node: f"A{idx + 1}" for idx, node in enumerate(nodes_0.keys())}
    node_labels_1 = {node: f"S{idx + 1}" for idx, node in enumerate(nodes_1.keys())}

    # Fusionar los diccionarios
    node_labels = {**node_labels_0, **node_labels_1}
    # print()
    # print(node_labels)

    # Compute transition frequencies
    transitions = [(tuple(abstractions[i]), tuple(abstractions[i + 1])) for i in range(len(abstractions) - 1)]
    transitions_freq = pd.Series(transitions).value_counts()

    # Create the graph
    G = nx.DiGraph()

    for node, freq in nodes_0.items():
        nodee = ', '.join(node)
        color = colour_abs_nodes
        label = f"{node_labels_0[node]} ({freq})"
        hover_info = f"[{nodee}]"
        G.add_node(node, size=1, color=color, label=label, hovertext=hover_info)

    for node, freq in set_of_sets2.items():
        nodee = ', '.join(node)
        color = colour_semiorg_nodes
        label = f"{Set_names[node]} ({freq})"
        hover_info = f"[{nodee}]"
        G.add_node(node, size=1, color=color, label=label, hovertext=hover_info)

    # Graficar nodos de intersección entre nodes_0 y set_of_sets2
    intersection_nodes = set(nodes_0.keys()).intersection(set_of_sets2.keys())
    for node in intersection_nodes:
        freq_nodes_0 = nodes_0[node]
        # freq_set_of_sets2 = set_of_sets2[node]
        nodee = ', '.join(node)
        color =  colour_cap_nodes
        label = f"{node_labels_0[node]} ∩ {Set_names[node]} ({freq_nodes_0})" #, {freq_set_of_sets2})"
        hover_info = f"[{nodee}]"
        G.add_node(node, size=1, color=color, label=label, hovertext=hover_info)

    for (source, target), weight in transitions_freq.items():
        G.add_edge(source, target, weight=weight)

    # pos = nx.drawing.nx_agraph.graphviz_layout(G, prog="dot")
    # pos = {node: (x, -y) for node, (x, y) in pos.items()}

    # Usar disposición con Graphviz y ajustar con un desplazamiento
    # pos = nx.shell_layout(G)
    pos = graphviz_layout(G, prog="dot")
    # # pos = nx.drawing.nx_agraph.graphviz_layout(G, prog="dot")
    offset_x, offset_y = 10, 10
    pos = {node: (x + offset_x, y + offset_y) for node, (x, y) in pos.items()}


    frames = []
    sliders_steps = []

    for i in range(len(abstract_time_series["Abstraction"])):
        node_sizes = [G.nodes[node]["size"] * first_node_size for node in G.nodes]
        node_colors = [G.nodes[node]["color"] for node in G.nodes]

        # Destacar los nodos actuales y ajustar tamaños progresivamente
        for offset, size_factor in enumerate(last_adjust_sizes[::-1]):
            max_nodes_to_plot = i + 1

            if offset < max_nodes_to_plot:
                current_node = tuple(abstractions.iloc[i - offset])

                if current_node in G.nodes:
                    index = list(G.nodes).index(current_node)

                    if offset >= 2:
                        prev_node_1 = tuple(abstractions.iloc[i - offset + 1])
                        prev_node_2 = tuple(abstractions.iloc[i - offset + 2])
                        if (
                            prev_node_1 in G.nodes
                            and prev_node_2 in G.nodes
                            and current_node != prev_node_1
                            and current_node != prev_node_2
                            and prev_node_1 != prev_node_2
                        ):
                            node_sizes[index] = last_node_size * size_factor
                        else:
                            max_size = max(node_sizes[index], last_node_size * size_factor)
                            node_sizes[index] = max_size
                    else:
                        max_size = max(node_sizes[index], last_node_size * size_factor)
                        node_sizes[index] = max_size

                    node_colors[index] = colour_last_nodes
                else:
                    continue
 
        edge_x = []                    # Lista para las coordenadas x de las aristas.
        edge_y = []                    # Lista para las coordenadas y de las aristas.
        annotations = []               # Lista para las anotaciones de las aristas.

        for source, target in G.edges: # Itera sobre las aristas.
            x0, y0 = pos[source]       # Obtiene las coordenadas del nodo fuente.
            x1, y1 = pos[target]       # Obtiene las coordenadas del nodo objetivo.
            edge_x += [x0, x1, None]   # Añade las coordenadas de la arista.
            edge_y += [y0, y1, None]   # Añade las coordenadas de la arista.

            annotations.append(dict(   # Añade anotaciones para cada arista.
                ax=x0, ay=y0,
                x=x1, y=y1,
                xref="x", yref="y",
                axref="x", ayref="y",
                showarrow=True,
                arrowhead=3,
                arrowsize=2,
                arrowwidth=1,
                arrowcolor="gray"
            ))

        frames.append(go.Frame(                                  # Añade un nuevo frame para la animación.
            data=[
                go.Scatter(x=[pos[node][0] for node in G.nodes], # Posiciones x de los nodos.
                        y=[pos[node][1] for node in G.nodes],    # Posiciones y de los nodos.
                        mode='markers+text',                     # Modo de visualización de los nodos.
                        marker=dict(size=node_sizes, color=node_colors, opacity=1), # Tamaño y color de los nodos.
                        text=[G.nodes[node]["label"] for node in G.nodes],          # Etiquetas de los nodos.
                        textposition="bottom center"),           # Posición del texto.
                go.Scatter(x=edge_x, y=edge_y,                   # Coordenadas de las aristas.
                        mode='lines',                            # Modo de visualización de las aristas.
                        line=dict(width=1, color='gray'))        # Estilo de las aristas.
            ],
            layout=dict(annotations=annotations),                # Añade las anotaciones a la visualización.
            name=f"Frame {i}"                                    # Nombre del frame para la animación.
        ))
        sliders_steps.append(dict(                               # Añade un paso al control deslizante de la animación.
            args=[[f"Frame {i}"], dict(frame=dict(duration=0, redraw=True), mode="immediate")],  # Define el paso.
            label=f"{i}",                                        # Etiqueta del paso.
            method="animate"                                     # Define el método de animación.
        ))

    legend_text = "<br>".join([f"{node_labels[node]}= [{', '.join(node)}]" for node in G.nodes])  # Genera el texto de la leyenda con etiquetas de nodos (S1, S2, ...) y sus elementos.
    # legend_text = "<br>".join([f"[{node}]" for node in G.nodes])  # Genera el texto de la leyenda con etiquetas de nodos (S1, S2, ...) y sus elementos.

    fig = go.Figure(
        data=[
            go.Scatter(
                x=[pos[node][0] for node in G.nodes],  # Coordenadas X de los nodos.
                y=[pos[node][1] for node in G.nodes],  # Coordenadas Y de los nodos.
                mode='markers+text',                   # Modo de visualización: muestra puntos (markers) y texto (text).
                marker=dict(size=node_sizes, color=node_colors, opacity=1),  # Configuración de los nodos: tamaño, color, y opacidad.
                text=[G.nodes[node]["label"] for node in G.nodes],  # Etiquetas que se muestran en cada nodo.
                textposition="bottom center",          # Posición del texto respecto al nodo (debajo y centrado).
                hovertext=[G.nodes[node]["hovertext"] for node in G.nodes],  # Texto que aparece al pasar el mouse sobre un nodo.
                hoverinfo="text",                      # Especifica que se muestra el texto de hover al interactuar.
                name=legend_text                       # Asocia el texto de la leyenda con el gráfico.
            )
        ],
        layout=go.Layout(
            title=title,                      # Título del gráfico.
            updatemenus=[                     # Configuración de los botones de control (play y pause).
                dict(
                    type="buttons",           # Define un grupo de botones.
                    showactive=False,         # Desactiva el resaltado del botón seleccionado.
                    buttons=[
                        dict(
                            label="Play",     # Etiqueta del botón "Play".
                            method="animate", # Método que activa la animación.
                            args=[None, dict(frame=dict(duration=interval, redraw=True), fromcurrent=True)]  # Configuración de la animación.
                        ),
                        dict(
                            label="Pause",    # Etiqueta del botón "Pause".
                            method="animate", # Método que pausa la animación.
                            args=[[None], dict(frame=dict(duration=0, redraw=False), mode="immediate")]  # Configuración de pausa.
                        )
                    ]
                )
            ],
            sliders=[{  # Configuración del control deslizante para avanzar manualmente entre los cuadros de la animación.
                'steps': sliders_steps,   # Pasos del deslizador.
                'currentvalue': {
                    'prefix': "Time: ",   # Prefijo que aparece junto al valor actual del deslizador.
                    'font': {'size': 16}, # Tamaño de fuente del prefijo.
                    'visible': True,      # Muestra el valor actual.
                },
                'x': 0.1,   # Posición horizontal del deslizador.
                'len': 0.9, # Longitud del deslizador en proporción al gráfico.
            }],
            xaxis=dict(showgrid=False, zeroline=False, visible=False),  # Oculta la cuadrícula y el eje X.
            yaxis=dict(showgrid=False, zeroline=False, visible=False),  # Oculta la cuadrícula y el eje Y.
            plot_bgcolor='rgba(0,0,0,0)', # Fondo transparente del gráfico.
            width=1400,                   # Ancho del gráfico.
            height=800,                   # Altura del gráfico.
            annotations=annotations,      # Anotaciones adicionales para destacar elementos.
            showlegend=showlegend,        # Habilita la visualización de la leyenda.
            legend=dict( 
                title=dict(text="Node Legend", font=dict(family="Arial", size=14, color="darkblue")),  # Título de la leyenda.
                bgcolor="white",  # Fondo blanco.
                bordercolor="white",                               # Borde blanco.
                borderwidth=2,                                     # Grosor del borde.
                font=dict(family="Arial", size=12, color="black"), # Fuente de los elementos.
                x=1,                                               # Posición horizontal a la derecha.
                y=1.0                                              # Posición vertical en la parte superior.
            )
        ),
        frames=frames  # Cuadros de la animación.
    )

    fig.write_html(filename)
    return filename

############################################################################################## 
 
def film_semiorganizations_abstractions_html(
                abstract_time_series, input_data,
                first_node_size=25, 
                colour_abs_nodes="lightcyan",  colour_semiorg_nodes="lightgreen", colour_cap_nodes="orange",
                last_node_size=35, last_adjust_sizes=[1/2, 3/4, 1], 
                colour_last_nodes="red", 
                interval=100, showlegend=True,
                title="Semi-organisations and Abstractions Film",
                filename="semiorg_abstraction_film.html"):
    """
    Generates an interactive animation of an abstraction graph from an abstract time series, 
    saving the animation as an HTML file. The animation shows how the graph evolves over time, 
    highlighting the most important nodes.

    Parameters:
    ----------
    abstract_time_series : pd.DataFrame
        A DataFrame containing two columns:
        - "Abstraction": A list of abstraction states over time.
        - "Time": Corresponding timestamps for each abstraction state.
    first_node_size : int, optional. The base size for nodes representing the abstraction states, by default 25.
    colour_abs_nodes : str, optional. Color for the nodes, by default "lightcyan".
    colour_semiorg_nodes : str, optional. Color for semi-organization nodes, by default "lightgreen".
    colour_cap_nodes : str, optional. Color for the intersection nodes, by default "orange".
    last_node_size : int, optional. The size of the most recent nodes in the animation, by default 35.
    last_adjust_sizes : list of float, optional. Adjustment multipliers for the size of the last three unique nodes, by default [1/4, 1/2, 1].
    colour_last_nodes : str, optional. Colors for the last three unique nodes, by default 'red'.
    interval : int, optional. Time interval between frames in milliseconds, by default 100.
    showlegend : bool, optional. Whether to display the legend in the animation, by default True.
    title : str, optional. Title of the animation, by default "Semi-organisations and Abstractions Film".
    filename : str, optional. Name of the output HTML file, by default "semiorg_abstraction_film.html".

    Returns:
    -------
    str
        The filename of the generated HTML file.
        
    Example:
    -------
    # Load an abstraction time series
    abstract_time_series = pd.DataFrame({
        'Time': [0, 1, 2],
        'Abstraction': [['A', 'B'], ['B'], ['C', 'D']]        
    })

    # Call the function to generate the animation
    plot_abstraction_graph_movie_html(abstract_time_series, filename="semiorg_abstraction_film.html")
    """
    # Call the `get_film_semiorganizations_abstractions_html` function to generate the HTML file
    get_film_semiorganizations_abstractions_html(
        abstract_time_series, input_data,
        first_node_size=first_node_size, 
        colour_abs_nodes=colour_abs_nodes, colour_semiorg_nodes=colour_semiorg_nodes, colour_cap_nodes=colour_cap_nodes,
        last_node_size=last_node_size, last_adjust_sizes=last_adjust_sizes, 
        colour_last_nodes=colour_last_nodes, 
        interval=interval, showlegend=showlegend,
        title=title,
        filename=filename
    )
    
    # Convert to an absolute path
    abs_path = os.path.abspath(filename) 
    # print("absolute path")
    # print(abs_path)
    # Check if the file was created correctly
    if not os.path.isfile(abs_path):
        print(f"File not found at {abs_path}")  # Additional message for debugging
        raise FileNotFoundError(f"The file {abs_path} was not found. Check if `rn_get_visualization` generated the file correctly.")
    
    # Inform the user about the file's location
    print(f"The animation was saved to: {abs_path}")
    
    # Open the HTML file in the default browser
    webbrowser.open(f"file://{abs_path}")

######################################################################################
# Histograms
######################################################################################

def plot_join_concentration_histogram(time_series, bins=30, alpha=0.7, color='skyblue', edgecolor='black', xlabel="Concentration", ylabel="Frequency", title="Histogram of Species Concentrations"):
    """
    Plots a histogram of combined species concentrations from the time series.

    Parameters:
    ----------
    time_series : pd.DataFrame
        Time series with a 'Time' column and species concentrations as additional columns.
    bins : int, optional
        Number of bins for the histogram (default is 30).
    alpha : float, optional
        Transparency level for the histogram bars (default is 0.7).
    color : str, optional
        Color of the histogram bars (default is 'skyblue').
    edgecolor : str, optional
        Color of the edges of the histogram bars (default is 'black').
    xlabel : str, optional
        Label for the x-axis (default is "Concentration").
    ylabel : str, optional
        Label for the y-axis (default is "Frequency").
    title : str, optional
        Title of the plot (default is "Histogram of Species Concentrations").
    """
    if 'Time' not in time_series.columns:
        raise ValueError("The DataFrame must include a 'Time' column for time values.")
    
    concentrations = time_series.drop(columns='Time').values.flatten()
    
    plt.figure(figsize=(10, 6))
    plt.hist(concentrations, bins=bins, alpha=alpha, color=color, edgecolor=edgecolor)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.show()


def plot_reaction_rate_mak_histogram(reaction_rate_maks, bins=10, color='skyblue', edgecolor='black', alpha=0.7, xlabel="Reaction Rate", ylabel="Frequency", title="Histogram of Reaction Rates"):
    """
    Generates a histogram for the reaction rates of species.

    Parameters:
    ----------
    reaction_rate_maks : pd.DataFrame
        DataFrame with reaction rates for species. Each column represents a species, and rows represent different time points.
    bins : int, optional
        Number of bins in the histogram (default is 10).
    color : str, optional
        Color of the histogram bars (default is 'skyblue').
    edgecolor : str, optional
        Color of the edges of the histogram bars (default is 'black').
    alpha : float, optional
        Transparency level for the histogram bars (default is 0.7).
    xlabel : str, optional
        Label for the x-axis (default is "Reaction Rate").
    ylabel : str, optional
        Label for the y-axis (default is "Frequency").
    title : str, optional
        Title of the histogram (default is "Histogram of Reaction Rates").
    """
    all_rates = reaction_rate_maks.values.flatten()
    
    plt.figure(figsize=(10, 6))
    plt.hist(all_rates, bins=bins, color=color, edgecolor=edgecolor, alpha=alpha)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.show()


def plot_species_histograms(reaction_rate_maks, species_names, bins=10, alpha=0.7, color="skyblue", edgecolor="black", xlabel="Concentration", ylabel="Frequency", title_plot="Histogram of"):
    """
    Generates individual histograms for each species in the reaction rates DataFrame.

    Parameters:
    ----------
    reaction_rate_maks : pd.DataFrame
        DataFrame with reaction rates for species.
    species_names : list
        List of species names (column names in the DataFrame).
    bins : int, optional
        Number of bins for the histograms (default is 10).
    alpha : float, optional
        Transparency level for the histogram bars (default is 0.7).
    color : str, optional
        Color of the histogram bars (default is "skyblue").
    edgecolor : str, optional
        Color of the edges of the histogram bars (default is "black").
    xlabel : str, optional
        Label for the x-axis (default is "Concentration").
    ylabel : str, optional
        Label for the y-axis (default is "Frequency").
    title_plot : str, optional
        Prefix for the titles of each subplot (default is "Histogram of").
    """
    num_species = len(species_names)
    fig, axes = plt.subplots(1, num_species, figsize=(5 * num_species, 4), sharey=True)

    if num_species == 1:
        axes = [axes]

    for ax, species in zip(axes, species_names):
        ax.hist(reaction_rate_maks[species], bins=bins, alpha=alpha, color=color, edgecolor=edgecolor, label=species)
        ax.set_title(f"{title_plot} {species}")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.legend()  # Adds the legend to the subplot.

    plt.tight_layout()
    plt.show()


def plot_combined_species_histogram(reaction_rate_maks, species_names, bins=10, alpha=0.7, edgecolor="black", xlabel="Concentration", ylabel="Frequency", title_plot="Combined Histogram"):
    """
    Generates a combined histogram of reaction rates for a set of species.

    Parameters:
    ----------
    reaction_rate_maks : pd.DataFrame
        DataFrame with reaction rates for species.
    species_names : list
        List of species names (columns in the DataFrame).
    bins : int, optional
        Number of bins for the histograms (default is 10).
    alpha : float, optional
        Transparency level for the histogram bars (default is 0.7).
    edgecolor : str, optional
        Color of the edges of the histogram bars (default is "black").
    xlabel : str, optional
        Label for the x-axis (default is "Concentration").
    ylabel : str, optional
        Label for the y-axis (default is "Frequency").
    title_plot : str, optional
        Title of the plot (default is "Combined Histogram").
    """
    plt.figure(figsize=(8, 6))

    for species in species_names:
        plt.hist(reaction_rate_maks[species], bins=bins, alpha=alpha, label=species, edgecolor=edgecolor)
    
    plt.title(title_plot)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(title="Species")
    plt.grid(axis="y", linestyle="--", alpha=alpha)
    plt.tight_layout()
    plt.show()

    #############################################################################################
def plot_series_MP(modules, RN, t_span=None, n_steps=None, species_labels=None):
    """
    Grafica los resultados de la simulación espaciotemporal de la metapoblación en subfiguras por especie.

    Parámetros:
        modules (dict): Resultados de la simulación para cada parche.
        grid_shape (tuple): Forma de la malla espacial (filas, columnas).
        t_span (tuple, opcional): Rango de tiempo para la simulación. Por defecto es None.
        n_steps (int, opcional): Número de pasos para la simulación. Por defecto es None.
        species_labels (list, opcional): Nombres de las especies. Por defecto es None.
        RN (dict, opcional): Red de reacciones, para acceder a `SpStr`.
    """
    if modules is None or len(modules) == 0:
        raise ValueError("Debe proporcionar los resultados de la simulación.")

    if n_steps is None:
        n_steps = 500  # Pasos por defecto  

    if t_span is None:
        t_span = (0, 20)  # Tiempo de simulación por defecto

    # Usar los nombres de las especies desde RN.SpStr si no se proporcionan labels
    if species_labels is None and RN is not None:
        species_labels = RN.SpStr
    elif species_labels is None:
        species_labels = [f"Especie {i + 1}" for i in range(next(iter(modules.values())).shape[1])]

    t = np.linspace(t_span[0], t_span[1], n_steps)
    n_species = len(species_labels)

    # Crear una figura con subplots para cada especie
    fig, axes = plt.subplots(n_species, 1, figsize=(10, 6 * n_species), sharex=True)

    if n_species == 1:
        axes = [axes]  # Asegurarse de que 'axes' sea iterable si hay una sola especie

    for idx, ax in enumerate(axes):
        for patch, values in modules.items():
            ax.plot(t, values[:, idx], label=f"{patch}") 
        ax.set_xlabel("Time")
        ax.set_ylabel(f"{species_labels[idx]}")
        ax.legend(loc="upper right", fontsize="small", ncol=1)
        ax.grid(True)

    plt.tight_layout()
    plt.show()