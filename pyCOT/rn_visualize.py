from collections import Counter
from collections import defaultdict
from graphviz import Digraph
from IPython.display import Image
# from PIL import Image
from pyvis.network import Network
import networkx as nx
import rustworkx as rx
import matplotlib.pyplot as plt
import mplcursors 
import os
import webbrowser
import sys
sys.stdout.reconfigure(encoding='utf-8')
# import codecs  # Agregado para manejar codificación UTF-8
    
import tempfile

from pyCOT.simulations import build_reaction_dict

##################################################################
# # Plot a reaction network in HTML:
##################################################################
def rn_get_visualization(rn, lst_color_spcs=None, lst_color_reacs=None, 
                         global_species_color=None, global_reaction_color=None,
                         global_input_edge_color=None, global_output_edge_color=None, 
                         node_size=20, shape_species_node='dot', shape_reactions_node='box', 
                         curvature=None, physics_enabled=False, 
                         filename="reaction_network.html"):
    """
    Visualizes a reaction network as an interactive HTML file.

    This function uses a given reaction network (RN) to generate a visualization in which 
    species and reactions are represented as nodes, and interactions are represented as edges. 
    The resulting visualization is saved to an HTML file that can be opened in a browser.

    Parameters:
    ----------
    rn : ReactionNetwork of rn_rustwork.py
        An object representing the reaction network, which includes species, 
        reactions, and their stoichiometric relationships.
    
    lst_color_spcs : list of tuples, optional
        A list of tuples specifying colors for specific species. Each tuple should have the 
        format (color, species_list), where species_list is a list of species names 
        to be assigned the given color. Example: [('cyan', ['A', 'B'])].
    
    lst_color_reacs : list of tuples, optional
        A list of tuples specifying colors for specific reactions. Each tuple should have 
        the format (color, reaction_list), where reaction_list is a list of reaction 
        names to be assigned the given color. Example: [('gray', ['R1', 'R2'])].
    
    global_species_color : str, optional
        Default color ('cyan') for all species nodes if no specific color is assigned in lst_color_spcs.
    
    global_reaction_color : str, optional
        Default color ('lightgray') for all reaction nodes if no specific color is assigned in lst_color_reacs.
    
    global_input_edge_color : str, optional
        Color for edges representing inputs (species consumed in reactions). Default is 'red'.
    
    global_output_edge_color : str, optional
        Color for edges representing outputs (species produced in reactions). Default is 'green'.
    
    node_size : int, optional
        The size of the nodes in the visualization. Default is 20.
    
    shape_species_node : str, optional
        Shape of the nodes representing species. Common options include 'dot', 'circle', 'box', 'ellipse', etc. Default is 'dot'.
    
    shape_reactions_node : str, optional
        Shape of the nodes representing reactions. Common options include 'box', 'dot',  'ellipse', etc. Default is 'box'.
    
    curvature : str or None, optional
        Determines whether edges are curved or straight. If 'curvedCCW' (curved counter-clockwise), 
        edges will be curved. Default is None (straight edges).
    
    physics_enabled : bool, optional
        If True, enables physics-based positioning for the nodes in the visualization. Default is False.
    
    filename : str, optional
        Name of the output HTML file where the visualization will be saved. Default is "reaction_network.html".

    Returns:
    -------
    str
        The filename of the saved HTML file.

    Example:
    -------
    # Load a reaction network object
    file_path = 'Txt/Farm.txt'
    rn = read_txt(file_path)

    # Visualize the reaction network, assigning yellow to species 's1'
    rn_get_visualization(rn, lst_color_spcs=[('yellow', ['s1'])], filename="reaction_network.html")
    """
    # Initialize the network visualization with specific size, settings, and directionality
    net = Network(height='750px', width='100%', notebook=True, directed=True, cdn_resources='in_line')
 
    ######################################     
    # Set physics options depending on ‘physics_enabled’.
    options = f"""
    var options = {{
        "physics": {{
            "enabled": {str(physics_enabled).lower()}
        }}
    }}
    """
    net.set_options(options)
    ######################################
    # Default colors: Define fallback colors for species, reactions, and edges if not provided
    default_species_color = global_species_color or 'cyan'
    default_reaction_color = global_reaction_color or 'lightgray'
    input_edge_color = global_input_edge_color or 'red'     # Color of input arrows to the reaction (Consume)
    output_edge_color = global_output_edge_color or 'green'  # Color of output arrows from the reaction (Produce)
    
    # Map specific colors to certain species or reactions if lists provided, otherwise use defaults
    species_colors = {species: color for color, species_list in (lst_color_spcs or []) for species in species_list}
    reaction_colors = {reaction: color for color, reaction_list in (lst_color_reacs or []) for reaction in reaction_list}
    
    #Construct RN dictionary from reactionnetwork object 
    RN_dict = build_reaction_dict(rn)

    ###################################### 
    # Identify all species in the network to check their presence in lst_color_spcs
    species_vector = sorted(set([spcs for reactants, products in RN_dict.values() for spcs, _ in reactants + products]))
    species_set = set(species_vector)  # Convert to a set for easier comparison
    
    # Validate species in lst_color_spcs to ensure they belong to the network's species
    if lst_color_spcs:
        for color, species_list in lst_color_spcs:
            for species in species_list:
                if species not in species_set:
                    print(f"Warning: The species '{species}' specified in lst_color_spcs does not belong to the species of the network.")

    # Identify all reactions in the network to check their presence in lst_color_reacs
    reaction_vector = list(RN_dict.keys())
    reaction_set = set(reaction_vector)  # Convert to a set for easier comparison

    # Validate reactions in lst_color_reacs to ensure they belong to the network's reactions
    if lst_color_reacs:
        for color, reaction_list in lst_color_reacs:
            for reaction in reaction_list:
                if reaction not in reaction_set:
                    print(f"Warning: The reaction '{reaction}' specified in lst_color_reacs does not belong to the network reactions.")
    
    # Calculate node size based on character length of species name for better visualization 
    for species in species_vector:
        color = species_colors.get(species, default_species_color)
        net.add_node(species, shape=shape_species_node, label=species, color=color, 
                     size=node_size, font={'size': 14, 'color': 'black'})
 
    ######################################

    # Add reaction nodes with color defined in lst_color_reacs or global/default color
    for reaction in reaction_vector:
        color = reaction_colors.get(reaction, default_reaction_color)
        net.add_node(reaction, shape=shape_reactions_node, label=reaction, color=color, size=10, font={'size': 14, 'color': 'black'})

    ######################################
    # Add edges and labels
    connections = set() # Set of connections to track opposite edges and avoid redundant edges
    for reaction, (inputs, outputs) in RN_dict.items():
        for species, coef in inputs:
            # coef = int(coef)  # Ensure coefficient is an integer
            if coef.is_integer():
                coef = int(coef)
            else:
                coef = float(coef)
                        
            edge_id = (species, reaction)
            default_curve='curvedCCW'
            smooth_type = {'type': default_curve} if curvature else None
            # Check for opposite edges and set edge direction and label accordingly
            if (reaction, species) in connections:
                if coef != 1:
                    net.add_edge(species, reaction, arrows='to', label=str(coef), color=input_edge_color, 
                                 font_color=input_edge_color, smooth=smooth_type)
                else:
                    net.add_edge(species, reaction, arrows='to', color=input_edge_color, 
                                 font_color=input_edge_color, smooth=smooth_type)
            else:
                if coef != 1:
                    net.add_edge(species, reaction, arrows='to', label=str(coef), color=input_edge_color, 
                                 font_color=input_edge_color)
                else:
                    net.add_edge(species, reaction, arrows='to', color=input_edge_color, 
                                 font_color=input_edge_color)

            # Add edge to connections set to avoid redundant edges
            connections.add(edge_id)

        for species, coef in outputs:
            # coef = int(coef)  # Ensure coefficient is an integer
            if coef.is_integer():
                coef = int(coef)
            else:
                coef = float(coef)

            edge_id = (reaction, species)
            default_curve='curvedCCW'
            smooth_type = {'type': default_curve} if curvature else None
            # Check for opposite edges and set edge direction and label accordingly
            if (species, reaction) in connections:
                if coef != 1:
                    net.add_edge(reaction, species, arrows='to', label=str(coef), color=output_edge_color, 
                                 font_color=output_edge_color, smooth=smooth_type)
                else:
                    net.add_edge(reaction, species, arrows='to', color=output_edge_color, 
                                 font_color=output_edge_color, smooth=smooth_type)
            else:
                if coef != 1:
                    net.add_edge(reaction, species, arrows='to', label=str(coef), color=output_edge_color, 
                                 font_color=output_edge_color)
                else:
                    net.add_edge(reaction, species, arrows='to', color=output_edge_color, 
                                 font_color=output_edge_color)

            # Add edge to connections set to avoid redundant edges
            connections.add(edge_id)
    
    ######################################            
    # Save the visualization to an HTML file 
    net.html = net.generate_html()
    with open(filename, "w", encoding="utf-8") as f:
        f.write(net.html)
    return filename

###########################################################################################
# Function to visualize the reaction network and open the HTML file in a browser
def rn_visualize_html(rn, lst_color_spcs=None, lst_color_reacs=None, 
                 global_species_color=None, global_reaction_color=None,
                 global_input_edge_color=None, global_output_edge_color=None, 
                 node_size=20, shape_species_node='dot', shape_reactions_node='box', 
                 curvature=None, physics_enabled=False, 
                 filename="reaction_network.html"):
    """
    Example:
    -------
    # Load a reaction network object
    file_path = 'Txt/Farm.txt'
    rn = read_txt(file_path)

    # Visualize the reaction network and open the HTML file in a browser
    rn_visualize_html(rn, lst_color_spcs=[('yellow', ['s1'])], filename="reaction_network.html")
    """    
    # Call the rn_get_visualization function to generate the HTML file
    rn_get_visualization(
        rn, lst_color_spcs, lst_color_reacs, 
        global_species_color, global_reaction_color,
        global_input_edge_color, global_output_edge_color, 
        node_size=node_size, shape_species_node=shape_species_node, shape_reactions_node=shape_reactions_node,
        curvature=curvature, physics_enabled=physics_enabled, filename=filename
    )
    
    # Convert to an absolute path
    abs_path = os.path.abspath(filename)  
    # Check if the file was created correctly
    if not os.path.isfile(abs_path):
        print(f"\nFile not found at {abs_path}")  # Additional message for debugging
        raise FileNotFoundError(f"\nThe file {abs_path} was not found. Check if rn_get_visualization generated the file correctly.")
    
    # Inform the user about the file's location
    print(f"\nThe visualization in HTML format of the reaction network was saved in:\n{abs_path}\n")
    
    # Open the HTML file in the default browser
    webbrowser.open(f"file://{abs_path}")

###########################################################################################
# # Plot the Hierarchy
###########################################################################################

# from pyCOT.rn_hierarchy import *

# def hierarchy_visualize(input_data):
#     """
#     Visualizes the containment hierarchy among sets.
    
#     Args:
#         input_data (list): List of lists or sets, where each element represents a set.
        
#     Displays a graph with nodes and containment relationships between sets. Smaller sets 
#     will be at lower levels, and interactive cursors are used to view set elements 
#     when hovering over nodes.
#     """
#     # Build the graph with species sets stored in nodes
#     G = hierarchy_build(input_data)
    
#     # Extract labels for visualization
#     labels = nx.get_node_attributes(G, 'label')
    
#     # Automatically assign positions by levels in the graph
#     levels = {}
#     for node in nx.topological_sort(G):  # Topologically sort the nodes
#         depth = 0
#         for parent in G.predecessors(node):  # Find predecessors of the current node
#             depth = max(depth, levels.get(parent, 0) + 1)  # Define the depth of the current node
#         levels[node] = depth
    
#     # Arrange nodes by levels and center them on the X-axis
#     pos = {}  # Dictionary to store node positions
#     nodes_by_level = defaultdict(list)  # Group nodes by containment levels
    
#     for node, level in levels.items():
#         nodes_by_level[level].append(node)  # Add nodes to the corresponding level
    
#     max_nodes_in_level = max(len(nodes) for nodes in nodes_by_level.values()) if nodes_by_level else 1
    
#     # Do not invert levels on the Y-axis; smaller sets will be at the bottom
#     for level, nodes in nodes_by_level.items():
#         num_nodes = len(nodes)  # Number of nodes at the current level
#         horizontal_spacing = max_nodes_in_level / (num_nodes + 1)  # Horizontal spacing between nodes
#         start_x = -((num_nodes - 1) * horizontal_spacing) / 2  # Initial X-coordinate to center nodes
        
#         # Calculate vertical position without inverting the Y-axis
#         for i, node in enumerate(nodes):
#             pos[node] = (start_x + i * horizontal_spacing, level)  # Assign position (X, Y) to each node
    
#     try:
#         import matplotlib.pyplot as plt
#         import mplcursors

#         # Draw the graph with labels and arrows, adjusting visual properties
#         plt.figure(figsize=(10, 8))
#         nx.draw(
#             G, pos, labels=labels, arrows=True, arrowstyle='-|>', arrowsize=20, 
#             node_size=1000, node_color="skyblue", font_size=12, font_weight="bold", 
#             edge_color="gray"
#         )

#         # Create cursor labels using the species sets stored in nodes
#         cursor_texts = {}
#         for node in G.nodes():
#             species_set = G.nodes[node]['species_set']
#             cursor_texts[node] = ', '.join(sorted(list(species_set))) if species_set else '∅'
        
#         # Configure interactive cursor to display set labels
#         cursor = mplcursors.cursor(plt.gca(), hover=True)
        
#         @cursor.connect("add")
#         def on_add(sel):
#             """Displays the content of the set when hovering over a node."""
#             node = list(G.nodes())[sel.index]
#             sel.annotation.set_text(cursor_texts.get(node, ''))
#             sel.annotation.set_visible(True)

#         @cursor.connect("remove")
#         def on_remove(sel):
#             """Hides the annotation when the cursor leaves the node."""
#             sel.annotation.set_visible(False)
        
#         plt.show()  # Display the graph
#     except ImportError:
#         print("Matplotlib or mplcursors not available, skipping visualization")
    
#     # Extract positions and store in a vector
#     node_positions = {node: pos[node] for node in G.nodes()}
    
#     return G, node_positions  # Return both the graph and the positions

# Function to access species sets from nodes
def get_species_from_node(G, node_name):
    """
    Get the species set stored in a node.
    
    Args:
        G (nx.DiGraph): The hierarchy graph
        node_name (str): The name of the node (e.g., "S1")
        
    Returns:
        set: The set of species stored in the node
    """
    if node_name in G:
        return G.nodes[node_name]['species_set']
    else:
        raise ValueError(f"Node {node_name} not found in the graph")

##################################################################
# Function to visualize the hierarchy of sets and save it as an HTML file
def hierarchy_get_visualization_html(
    input_data, 
    node_size=20, 
    node_color="cyan", 
    edge_color="gray", 
    shape_node='dot',
    lst_color_subsets=None,
    node_font_size=14, 
    edge_width=2,  
    filename="hierarchy_visualization.html"
):
    """
    Visualizes the containment hierarchy among sets with automatic positions (inverted hierarchy).
    
    Args:
        input_data (list): List of lists or sets, where each element represents a set.
        node_size (int): Size of the nodes in the visualization.
        node_color (str): Default color of the nodes in the visualization.
        edge_color (str): Color of the edges in the visualization.
        shape_node (str): Specifies the shape of the nodes in the visualization, options include 'dot', 'square', 'triangle', 'star', 'diamond' (default is 'dot').
        lst_color_subsets (list): List of tuples with color and subsets to highlight specific nodes. Example: lst_color_subsets = [("red", [{"A", "B"}, {"A", "C"}])]
        node_font_size (int): Font size of the node labels.
        edge_width (int): Width of the edges. 
        filename (str): Name of the output HTML file.
    
    Saves an interactive HTML graph with automatic positioning (smallest sets at the bottom).
    """
    # Convert input data to unique sets
    unique_subsets = []
    for sublist in input_data:
        if set(sublist) not in [set(x) for x in unique_subsets]:
            unique_subsets.append(sublist)

    # Create a list of unique sets for reference
    Set_of_sets = [set(s) for s in unique_subsets]

    # Sort the sets by size (from smallest to largest)
    Set_of_sets.sort(key=lambda x: len(x))
    
    # Create set names based on their level
    Set_names = [f"S{i+1}" for i in range(len(Set_of_sets))]
    
    # Create a dictionary of labels for the nodes
    labels = {f"S{i+1}": f"{', '.join(sorted(s))}" for i, s in enumerate(Set_of_sets)}

    # Create set labels to display when hovering over nodes
    cursor_labels = [
        ', '.join(sorted(list(s))) if s else '∅'  # Use ∅ for empty sets
        for s in Set_of_sets
    ]

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

    # Assign colors to nodes based on lst_color_subsets
    color_map = {}
    if lst_color_subsets:
        for color, subsets in lst_color_subsets:
            for subset in subsets:
                for i, s in enumerate(Set_of_sets):
                    if s == set(subset):
                        color_map[Set_names[i]] = color

    # Add nodes to the graph
    for name, hover_text in zip(Set_names, cursor_labels):
        color = color_map.get(name, node_color)  # Use the specific color if defined, otherwise use the default
        net.add_node(
            name,  # Node identifier
            label=name,  # Node label to be displayed inside the node
            title=hover_text,  # Text shown on hover
            color=color,  # Node color
            size=node_size,  # Node size
            font={"size": node_font_size},  # Font size of the label
            shape=shape_node  # Node shape (e.g., circle, to ensure the label fits inside)
        )
 

    # Add edges based on containment relationships
    for i, child_set in enumerate(Set_of_sets):
        for j, parent_set in enumerate(Set_of_sets):
            if i < j and child_set.issubset(parent_set):  # Add edge if child is a subset of parent
                # Ensure no intermediate set exists between child and parent
                is_direct = True
                for k, intermediate_set in enumerate(Set_of_sets):
                    if i < k < j and child_set.issubset(intermediate_set) and intermediate_set.issubset(parent_set):
                        is_direct = False
                        break
                if is_direct:
                    net.add_edge(Set_names[i], Set_names[j])

    # Print tuple with labels and subsets
    label_set_pairs = [(label, s) for label, s in zip(Set_names, Set_of_sets)]
    print(f"\nTuples of Labels and Subsets:\n{label_set_pairs}")

    # Save the visualization to an HTML file 
    net.html = net.generate_html() 
    with open(filename, "w", encoding="utf-8") as f:
        f.write(net.html)
    return net, label_set_pairs

##################################################################
# Function to visualize the hierarchy of sets and open the HTML file in a browser   
def hierarchy_visualize_html(
    input_data, 
    node_size=20, 
    node_color="cyan", 
    edge_color="gray", 
    shape_node='dot',    
    lst_color_subsets=None, 
    node_font_size=14, 
    edge_width=2,  
    filename="hierarchy_visualization.html"
):        
    """
    Wrapper function to generate and visualize the containment hierarchy among sets as an HTML file.
    
    Args:
        input_data (list): List of lists or sets, where each element represents a set.
        node_size (int): Size of the nodes in the visualization.
        node_color (str): Default color of the nodes in the visualization.
        edge_color (str): Color of the edges in the visualization.
        shape_node (str): Specifies the shape of the nodes in the visualization, options include 'dot', 'square', 'triangle', 'star', 'diamond' (default is 'dot').
        lst_color_subsets (list): List of tuples with color and subsets to highlight specific nodes. Example: lst_color_subsets = [("red", [{"A", "B"}, {"A", "C"}])]
        node_font_size (int): Font size of the node labels.
        edge_width (int): Width of the edges. 
        filename (str): Name of the output HTML file.
    
    Generates the visualization and opens it in the default web browser.
    """
    # Call the hierarchy_get_visualization_html function to generate the HTML file with the additional parameters
    hierarchy_get_visualization_html(
        input_data,  
        node_size=node_size, 
        node_color=node_color, 
        edge_color=edge_color, 
        shape_node=shape_node,
        lst_color_subsets=lst_color_subsets,  
        node_font_size=node_font_size, 
        edge_width=edge_width,  
        filename=filename
    )
    
    # Convert to an absolute path
    abs_path = os.path.abspath(filename) 
    
    # Check if the file was created correctly
    if not os.path.isfile(abs_path):
        print(f"File not found at {abs_path}")  # Additional message for debugging
        raise FileNotFoundError(f"The file {abs_path} was not found. Check if hierarchy_get_visualization_html generated the file correctly.")
    
    # Inform the user about the file's location
    print(f"\nThe hierarchy visualization was saved to:\n{abs_path}\n")
    
    # Open the HTML file in the default browser
    webbrowser.open(f"file://{abs_path}")

def rn_get_string(rn):
    # Print basic information about the loaded network
    print(f"Loaded reaction network with:")
    print(f"  - {len(rn.species())} species")
    print(f"  - {len(rn.reactions())} reactions")

    # Print the species
    print("\nSpecies:")
    for species in rn.species():
        print(f"  - {species.name}")

    # Print the reactions
    print("\nReactions:")
    for reaction in rn.reactions():
        support = " + ".join([f"{edge.coefficient}*{edge.species_name}" if edge.coefficient != 1 else edge.species_name 
                                for edge in reaction.support_edges()])
        products = " + ".join([f"{edge.coefficient}*{edge.species_name}" if edge.coefficient != 1 else edge.species_name 
                                for edge in reaction.products_edges()])
        
        if not support:
            support = "∅"  # Empty set symbol for inflow reactions
        if not products:
            products = "∅"  # Empty set symbol for outflow reactions
            
        print(f"  - {reaction.name()}: {support} => {products}")

##################################################################
# # Plot a bipartite metabolic network graph from a ReactionNetwork object
##################################################################
# Function to create a bipartite graph from a ReactionNetwork object
def create_bipartite_graph_from_rn(rn):
    graph = rx.PyDiGraph()
    metabolite_nodes = {}
    reaction_nodes = {} 

    for specie in rn.species():
        idx = graph.add_node(('specie', specie.name))
        metabolite_nodes[specie.name] = idx

    for reaction in rn.reactions():
        idx = graph.add_node(('reaction', reaction.name()))
        reaction_nodes[reaction.name()] = idx 

    # Agregar aristas según edges de cada reacción
    for reaction in rn.reactions():
        rxn_name = reaction.name()
        rxn_node = reaction_nodes[rxn_name]

        for edge in reaction.edges:
            met_name = edge.species_name
            coeff = edge.coefficient
            if edge.type == 'reactant':
                graph.add_edge(metabolite_nodes[met_name], rxn_node, coeff)
            elif edge.type == 'product':
                graph.add_edge(rxn_node, metabolite_nodes[met_name], coeff)

    return graph, metabolite_nodes, reaction_nodes

# Function to plot the bipartite graph using Graphviz library
def plot_graph_with_graphviz(
    graph,
    lst_color_spcs=None,
    lst_color_reacs=None,
    global_species_color='cyan',
    global_reaction_color='lightgray',
    global_input_edge_color='red',
    global_output_edge_color='green',
    node_size=20,
    shape_species_node='circle',
    shape_reactions_node='box',
    filename="metabolic_network"
):
    dot = Digraph(comment="Bipartite Metabolic Network")

    # Crear diccionarios para color específico por nombre
    species_colors = {species: color for color, species_list in (lst_color_spcs or []) for species in species_list}
    reaction_colors = {reaction: color for color, reaction_list in (lst_color_reacs or []) for reaction in reaction_list}

    for idx, (tipo, nombre) in enumerate(graph.nodes()):
        if tipo == 'specie':
            color = species_colors.get(nombre, global_species_color)
            shape = shape_species_node
        else:
            color = reaction_colors.get(nombre, global_reaction_color)
            shape = shape_reactions_node

        dot.node(
            str(idx),
            nombre,
            shape=shape,
            style='filled',
            fillcolor=color,
            width=str(node_size / 72),  # Aproximación para escalar (Graphviz usa pulgadas)
            fontsize='10'
        )

    for src, dst in graph.edge_list():
        data = graph.get_edge_data(src, dst)
        label = str(data)

        # Determinar tipo de arista
        src_tipo = graph.nodes()[src][0]
        color = global_input_edge_color if src_tipo == 'specie' else global_output_edge_color

        # No mostrar etiqueta si peso = 1
        if label == '1':
            dot.edge(str(src), str(dst), color=color)
        else:
            dot.edge(str(src), str(dst), label=label, color=color)

    dot.render(filename, format='png', cleanup=True)
    full_path = os.path.abspath(f"{filename}.png")
    print(f"Reaction network saved as: {full_path}")

    return Image(filename + '.png')

# def plot_graph_with_graphviz(graph, filename="metabolic_network"):
#     dot = Digraph(comment="Bipartite Metabolic Network")

#     for idx, (tipo, nombre) in enumerate(graph.nodes()):
#         if tipo == 'specie':
#             dot.node(str(idx), nombre, shape='circle', style='filled', fillcolor='cyan')
#         else:
#             dot.node(str(idx), nombre, shape='box', style='filled', fillcolor='lightgray')

#     for src, dst in graph.edge_list():
#         data = graph.get_edge_data(src, dst)
#         color = 'red' if graph.nodes()[src][0] == 'specie' else 'green'
#         # dot.edge(str(src), str(dst), label=str(data), color=color)
#         label = str(data)
#         if label == '1':
#             dot.edge(str(src), str(dst), color=color)
#         else:
#             dot.edge(str(src), str(dst), label=label, color=color)


#     dot.render(filename, format='png', cleanup=True)
#     full_path = os.path.abspath(f"{filename}.png")
#     print(f"Reaction network saved as: {full_path}")

#     return Image(filename + '.png')

# Function to get the visualization of the bipartite metabolic network graph as an interactive HTML file using the pyvis library
def get_visualize_metabolic_network(graph, lst_color_spcs=None, lst_color_reacs=None, 
                         global_species_color=None, global_reaction_color=None,
                         global_input_edge_color=None, global_output_edge_color=None, 
                         node_size=20, shape_species_node='dot', shape_reactions_node='box', 
                         curvature=None, physics_enabled=False, filename="metabolic_network.html"):
    """
    Visualizes a metabolic network graph as an interactive HTML file.

    This function uses a given graph to generate a visualization in which 
    species and reactions are represented as nodes, and interactions are represented as edges. 
    The resulting visualization is saved to an HTML file that can be opened in a browser.

    Parameters:
    ----------
    graph : rustworkx.PyDiGraph
        A directed graph representing the metabolic network, where nodes contain 
        (tipo, nombre) tuples indicating node type and name.
    
    lst_color_spcs : list of tuples, optional
        A list of tuples specifying colors for specific species. Each tuple should have the 
        format (color, species_list), where species_list is a list of species names 
        to be assigned the given color. Example: [('cyan', ['A', 'B'])].
    
    lst_color_reacs : list of tuples, optional
        A list of tuples specifying colors for specific reactions. Each tuple should have 
        the format (color, reaction_list), where reaction_list is a list of reaction 
        names to be assigned the given color. Example: [('gray', ['R1', 'R2'])].
    
    global_species_color : str, optional
        Default color for all species nodes if no specific color is assigned in lst_color_spcs.
        Default is 'cyan'.
    
    global_reaction_color : str, optional
        Default color for all reaction nodes if no specific color is assigned in lst_color_reacs.
        Default is 'lightgray'.
    
    global_input_edge_color : str, optional
        Color for edges representing inputs (species consumed in reactions). Default is 'red'.
    
    global_output_edge_color : str, optional
        Color for edges representing outputs (species produced in reactions). Default is 'green'.
    
    node_size : int, optional
        The size of the nodes in the visualization. Default is 20.
    
    shape_species_node : str, optional
        Shape of the nodes representing species. Common options include 'dot', 'circle', 'box', 'ellipse', etc. Default is 'dot'.
    
    shape_reactions_node : str, optional
        Shape of the nodes representing reactions. Common options include 'box', 'dot', 'ellipse', etc. Default is 'box'.
    
    curvature : str or None, optional
        Determines whether edges are curved or straight. If 'curvedCCW' (curved counter-clockwise), 
        edges will be curved. Default is None (straight edges).
    
    physics_enabled : bool, optional
        If True, enables physics-based positioning for the nodes in the visualization. Default is False.
    
    filename : str, optional
        Name of the output HTML file where the visualization will be saved. Default is "metabolic_network.html".

    Returns:
    -------
    str
        The filename of the saved HTML file.
    """
    # Initialize the network visualization
    
    net = Network(height='750px', width='100%', notebook=True, directed=True, cdn_resources='in_line')
    # net = Network(height='100%', width='100%', notebook=True, directed=True, cdn_resources='in_line')
    # net = Network(height='750px', width='100%', notebook=True, directed=True, cdn_resources='in_line')
    # net = Network(height="100%", width="100%", directed=True)
    
    # Configure physics
    if physics_enabled:
        net.barnes_hut()
    else:
        net.toggle_physics(False)
    
    # Default colors
    default_species_color = global_species_color or 'cyan'
    default_reaction_color = global_reaction_color or 'lightgray'
    input_edge_color = global_input_edge_color or 'red'
    output_edge_color = global_output_edge_color or 'green'
    
    # Map specific colors to certain species or reactions
    species_colors = {species: color for color, species_list in (lst_color_spcs or []) for species in species_list}
    reaction_colors = {reaction: color for color, reaction_list in (lst_color_reacs or []) for reaction in reaction_list}
    
    # Convert rustworkx graph to networkx for easier manipulation
    g_nx = nx.DiGraph()
    node_mapping = {}  # Map rustworkx indices to node names
    species_set = set()
    reaction_set = set()
    
    for idx, (tipo, nombre) in enumerate(graph.nodes()):
        g_nx.add_node(idx, tipo=tipo, nombre=nombre)
        node_mapping[idx] = nombre
        if tipo == 'specie':
            species_set.add(nombre)
        else:
            reaction_set.add(nombre)

    # Validate species in lst_color_spcs
    if lst_color_spcs:
        for color, species_list in lst_color_spcs:
            for species in species_list:
                if species not in species_set:
                    print(f"Warning: The species '{species}' specified in lst_color_spcs does not belong to the species of the network.")

    # Validate reactions in lst_color_reacs
    if lst_color_reacs:
        for color, reaction_list in lst_color_reacs:
            for reaction in reaction_list:
                if reaction not in reaction_set:
                    print(f"Warning: The reaction '{reaction}' specified in lst_color_reacs does not belong to the network reactions.")

    # Count edge occurrences
    edge_counts = Counter((src, dst) for src, dst in graph.edge_list())
    
    for src, dst in graph.edge_list():
        g_nx.add_edge(src, dst, weight=graph.get_edge_data(src, dst))

    # Get positions using Graphviz layout if physics is disabled
    if not physics_enabled:
        try:
            pos = nx.nx_pydot.graphviz_layout(g_nx, prog="dot")
        except:
            # Fallback to spring layout if graphviz is not available
            pos = nx.spring_layout(g_nx)
    else:
        pos = {}

    # Add nodes to the visualization
    for idx in g_nx.nodes():
        tipo = g_nx.nodes[idx]['tipo']
        nombre = g_nx.nodes[idx]['nombre']
        
        # Determine node color
        if tipo == 'specie':
            color = species_colors.get(nombre, default_species_color)
            shape = shape_species_node
        else:
            color = reaction_colors.get(nombre, default_reaction_color)
            shape = shape_reactions_node
        
        # Add node with or without fixed position
        if not physics_enabled and idx in pos:
            x, y = pos[idx]
            net.add_node(n_id=idx, label=nombre, shape=shape, color=color, 
                        size=node_size, x=x, y=-y, physics=False,
                        font={'size': 14, 'color': 'black'})
        else:
            net.add_node(n_id=idx, label=nombre, shape=shape, color=color, 
                        size=node_size, font={'size': 14, 'color': 'black'})

    # Add edges with proper styling
    edge_usage = Counter()
    
    connections = set()
    
    for src, dst in g_nx.edges():
        weight = g_nx.edges[src, dst]['weight']
        count = edge_counts[(src, dst)]
        edge_usage[(src, dst)] += 1
        
        # Determine edge color based on node types
        src_tipo = g_nx.nodes[src]['tipo']
        dst_tipo = g_nx.nodes[dst]['tipo']
        
        # Species -> Reaction (input/consumption)
        if src_tipo == 'specie' and dst_tipo == 'reaction':
            edge_color = input_edge_color
        # Reaction -> Species (output/production)
        elif src_tipo == 'reaction' and dst_tipo == 'specie':
            edge_color = output_edge_color
        else:
            edge_color = 'gray'  # Default for other cases
        
        # Configure edge smoothing based on curvature and edge usage
        smooth_config = {}

        if curvature:
            if count != 1 or (dst, src) in connections:
                # Caso con conexión opuesta y múltiples aristas
                curve_type = "cubicBezier" if edge_usage[(src, dst)] % 2 == 0 else "curvedCCW"
                smooth_config = {
                    "smooth": {
                        "type": curve_type,
                        "roundness": 0.3
                    }
                }
            else:
                # Caso único o sin conexión opuesta
                smooth_config = {
                    "smooth": {
                        "type": "cubicBezier",
                        "forceDirection": "vertical",
                        "roundness": 0.4
                    }
                }
        elif count != 1 or (dst, src) in g_nx.edges():
            # Caso sin curvature, pero con múltiples aristas o arista opuesta
            curve_type = "curvedCW" if edge_usage[(src, dst)] % 2 == 0 else "curvedCCW"
            smooth_config = {
                "smooth": {
                    "type": curve_type,
                    "forceDirection": "vertical",
                    "roundness": 0.4
                }
            }
        
        # Add edge with weight label if greater than 1
        if weight != 1: 
            net.add_edge(src, dst, label=str(weight), 
                        arrows='to', 
                        color=edge_color, font_color=edge_color, **smooth_config)
        else: 
            net.add_edge(src, dst, arrows='to', 
                        color=edge_color, font_color=edge_color, **smooth_config)

        connections.add((src, dst)) 

    ######################################            
    # Save the visualization to an HTML file 
    net.html = net.generate_html()
    with open(filename, "w", encoding="utf-8") as f:
        f.write(net.html)
    return filename 

    # # Ensure filename has .html extension
    # if not filename.endswith('.html'):
    #     filename += '.html'
    
    # # SOLUCIÓN AL PROBLEMA DE CODIFICACIÓN: Reemplazar net.write_html(filename) con:
    # try:
    #     # Generar contenido HTML
    #     html_content = net.generate_html()
        
    #     # Escribir con codificación UTF-8
    #     with codecs.open(filename, 'w', encoding='utf-8') as f:
    #         f.write(html_content)
        
    #     print(f"Visualización guardada exitosamente en: {filename}")
        
    # except Exception as e:
    #     print(f"Error al guardar el archivo: {e}")
    #     # Respaldo: intentar con el método original
    #     try:
    #         net.write_html(filename)
    #         print(f"Guardado con método alternativo en: {filename}")
    #     except:
    #         print("No se pudo guardar el archivo HTML")
    #         raise
    
    # return filename   
 
    # # Ensure filename has .html extension
    # if not filename.endswith('.html'):
    #     filename += '.html'
    
    # # Save the visualization
    # net.write_html(filename)
    
    # return filename

# Function to visualize the metabolic network graph as an interactive HTML file using the pyvis library
def visualize_metabolic_network(graph, lst_color_spcs=None, lst_color_reacs=None, 
                         global_species_color=None, global_reaction_color=None,
                         global_input_edge_color=None, global_output_edge_color=None, 
                         node_size=20, shape_species_node='dot', shape_reactions_node='box', 
                         curvature=None, physics_enabled=False, filename="metabolic_network.html"):
    """
    Example:
    -------
    # Visualize an interactive graph object 
    graph = create_bipartite_graph_from_rn()
    get_plot_graph_interactive(graph, lst_color_spcs=[('blue', ['node1'])], filename="interactive_graph.html")
    """    
    # Call a function to generate the interactive graph HTML file
    get_visualize_metabolic_network(
        graph, lst_color_spcs, lst_color_reacs, 
        global_species_color, global_reaction_color,
        global_input_edge_color, global_output_edge_color, 
        node_size=node_size, shape_species_node=shape_species_node, shape_reactions_node=shape_reactions_node,
        curvature=curvature, physics_enabled=physics_enabled, filename=filename
    )
    
    # Convert the filename to an absolute path
    abs_path = os.path.abspath(filename)  
    
    # Check if the file was created successfully
    if not os.path.isfile(abs_path):
        print(f"\nFile not found at {abs_path}")  # Debugging message
        raise FileNotFoundError(f"\nThe file {abs_path} was not found. Check if get_visualize_metabolic_network created the file correctly.")
    
    # Inform the user about the file's location
    print(f"\nThe visualization of the bipartite graph of the metabolic network was saved in HTML format:\n{abs_path}\n")
    
    # Open the HTML file in the default browser
    webbrowser.open(f"file://{abs_path}")