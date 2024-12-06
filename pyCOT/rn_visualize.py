from pyvis.network import Network 
import networkx as nx
import matplotlib.pyplot as plt
import mplcursors
import webbrowser  # Allows opening URLs or local files in the system's default browser
import os  # For handling paths and checking file existence
from collections import defaultdict
import sys
sys.stdout.reconfigure(encoding='utf-8')

##################################################################
# # Plot a reaction network in HTML:
##################################################################
def rn_get_visualization(RN, lst_color_spcs=None, lst_color_reacs=None, 
                 global_species_color=None, global_reaction_color=None,
                 global_input_edge_color=None, global_output_edge_color=None, 
                 node_size=1000, curvature=None, physics_enabled=False, 
                 filename="reaction_network.html"):
    """
    Visualizes a reaction network in an interactive HTML file. 

    Parameters:
    - RN: ReactionNetwork object 
    - lst_color_spcs (list, optional): List of tuples that associates specific colors with species (e.g., [('cyan', ['A', 'B'])]).
    - lst_color_reacs (list, optional): List of tuples that associates specific colors with reactions (e.g., [('gray', ['R1', 'R2'])]).
    - global_species_color (str, optional): Default color for species, if no other color is specified.
    - global_reaction_color (str, optional): Default color for reactions, if no other color is specified.
    - global_input_edge_color (str, optional): Color of input arrows (consumption) to reactions.
    - global_output_edge_color (str, optional): Color of output arrows (production) from reactions.
    - node_size (int, optional): Base size for network nodes.
    - curvature (str, optional): None=False (Without curvature), True=curvedCCW (Counter Clockwise)
    - filename (str, optional): Name of the output HTML file.    

    Returns:
    - None. The function saves the visualization to an HTML file specified by `filename`.

    Example of usage: 
    # Define the file path where the reaction network model is located.
    file_path = 'Txt/Farm.txt'  

    # Load the reaction network from the specified file in `file_path`.
    testRN = load_pyCOT_from_file(file_path)  

    # Visualize the reaction network `testRN`, assigning the color yellow to species 's1' and saving the result in an HTML file.
    rn_visualize(testRN.RN, lst_color_spcs=[('yellow', ['s1'])], filename="Reaction_network.html")  
   
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

    # Default colors: Define fallback colors for species, reactions, and edges if not provided
    default_species_color = global_species_color or 'cyan'
    default_reaction_color = global_reaction_color or 'lightgray'
    input_edge_color = global_input_edge_color or 'red'     # Color of input arrows to the reaction (Consume)
    output_edge_color = global_output_edge_color or 'green'  # Color of output arrows from the reaction (Produce)
    
    # Map specific colors to certain species or reactions if lists provided, otherwise use defaults
    species_colors = {species: color for color, species_list in (lst_color_spcs or []) for species in species_list}
    reaction_colors = {reaction: color for color, reaction_list in (lst_color_reacs or []) for reaction in reaction_list}
    
    #Construct RN dictionary from reactionnetwork object 
    RN_dict = {}  
    
    # Lista para almacenar los nombres de las reacciones
    reaction_names = []  

    # Iterar sobre las reacciones
    for i in range(len(RN.RnStr)):  # Usamos RnStr directamente, ya que RN aún no está completamente construido
        reactants_coeff = []
        products_coeff = []

        # Procesar reactivos
        for j in range(len(RN.SpStr)):  # Usamos SpStr, ya que SpStr contiene las especies
            coef_r = RN.RnMsupp[i][j]  # Coeficiente del reactivo en la reacción i para la especie j
            if coef_r > 0:
                reactants_coeff.append((RN.SpStr[j], coef_r))
        
        # Procesar productos
        for j in range(len(RN.SpStr)):  # Lo mismo para los productos
            coef_r = RN.RnMprod[i][j]  # Coeficiente del producto en la reacción i para la especie j
            if coef_r > 0:
                products_coeff.append((RN.SpStr[j], coef_r))

        # Crear el nombre de la reacción
        reaction_name = RN.RnStr[i] #f'R{i+1}'

        # Guardar la reacción en el diccionario
        RN_dict[reaction_name] = (reactants_coeff, products_coeff)

        # Agregar el nombre de la reacción a la lista de nombres de reacciones
        reaction_names.append(reaction_name)

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
    longest_species = max(RN.SpStr, key=len)
    node_size = len(longest_species) * 10   
    # Add species nodes with color defined in `lst_color_spcs` or global/default color
    for species in species_vector:
        color = species_colors.get(species, default_species_color)
        net.add_node(species, shape='circle', label=species, color=color, 
                     size=node_size, font={'size': 14, 'color': 'black'})
 
################################################################################################ 

    # Add reaction nodes with color defined in `lst_color_reacs` or global/default color
    for reaction in reaction_vector:
        color = reaction_colors.get(reaction, default_reaction_color)
        net.add_node(reaction, shape='box', label=reaction, color=color, size=10, font={'size': 14, 'color': 'black'})

    ######################################
    # Add edges and labels
    connections = set() # Set of connections to track opposite edges and avoid redundant edges
    for reaction, (inputs, outputs) in RN_dict.items():
        for species, coef in inputs:
            edge_id = (species, reaction)
            default_curve='curvedCCW'
            smooth_type = {'type': default_curve} if curvature else None
            # Check for opposite edges and set edge direction and label accordingly
            if (reaction, species) in connections:
                if coef > 1:
                    net.add_edge(species, reaction, arrows='to', label=str(coef), color=input_edge_color, 
                                 font_color=input_edge_color, smooth=smooth_type)
                else:
                    net.add_edge(species, reaction, arrows='to', color=input_edge_color, 
                                 font_color=input_edge_color, smooth=smooth_type)
            else:
                if coef > 1:
                    net.add_edge(species, reaction, arrows='to', label=str(coef), color=input_edge_color, 
                                 font_color=input_edge_color)
                else:
                    net.add_edge(species, reaction, arrows='to', color=input_edge_color, 
                                 font_color=input_edge_color)

            # Add edge to connections set to avoid redundant edges
            connections.add(edge_id)

        for species, coef in outputs:
            edge_id = (reaction, species)
            default_curve='curvedCCW'
            smooth_type = {'type': default_curve} if curvature else None
            # Check for opposite edges and set edge direction and label accordingly
            if (species, reaction) in connections:
                if coef > 1:
                    net.add_edge(reaction, species, arrows='to', label=str(coef), color=output_edge_color, 
                                 font_color=output_edge_color, smooth=smooth_type)
                else:
                    net.add_edge(reaction, species, arrows='to', color=output_edge_color, 
                                 font_color=output_edge_color, smooth=smooth_type)
            else:
                if coef > 1:
                    net.add_edge(reaction, species, arrows='to', label=str(coef), color=output_edge_color, 
                                 font_color=output_edge_color)
                else:
                    net.add_edge(reaction, species, arrows='to', color=output_edge_color, 
                                 font_color=output_edge_color)

            # Add edge to connections set to avoid redundant edges
            connections.add(edge_id)

    # Save the visualization to an HTML file 
    with open(filename, "w", encoding="utf-8") as f:
        f.write(net.html)
    return f
 
##################################################################

def rn_visualize(RN, lst_color_spcs=None, lst_color_reacs=None, 
                 global_species_color=None, global_reaction_color=None,
                 global_input_edge_color=None, global_output_edge_color=None, 
                 node_size=1000, curvature=None, physics_enabled=False, 
                 filename="reaction_network.html"):
    
    # Call the `rn_get_visualization` function to generate the HTML file
    rn_get_visualization(
        RN, lst_color_spcs, lst_color_reacs, 
        global_species_color, global_reaction_color,
        global_input_edge_color, global_output_edge_color, 
        node_size=node_size, curvature=curvature, 
        physics_enabled=physics_enabled, filename=filename
    )
    
    # Convert to an absolute path
    abs_path = os.path.abspath(filename) 
    
    # Check if the file was created correctly
    if not os.path.isfile(abs_path):
        print(f"File not found at {abs_path}")  # Additional message for debugging
        raise FileNotFoundError(f"The file {abs_path} was not found. Check if `rn_get_visualization` generated the file correctly.")
    
    # Inform the user about the file's location
    print(f"The reaction network visualization was saved to: {abs_path}")
    
    # Open the HTML file in the default browser
    webbrowser.open(f"file://{abs_path}") 

##################################################################
# # Plot the Hierarchy
##################################################################

def hierarchy_visualize(input_data):
    """
    Visualizes the containment hierarchy among sets.
    
    Args:
        input_data (list): List of lists or sets, where each element represents a set.
        
    Displays a graph with nodes and containment relationships between sets. Smaller sets 
    will be at lower levels, and interactive cursors are used to view set elements 
    when hovering over nodes.
    """
    
    # Create a list of tuples (set name, set) for reference
    Set_of_sets = [(set(s)) for s in input_data]
    
    # Dictionary to store the containment graph
    containment_graph = defaultdict(list)
    
    #######################################################################
    # Helper function to check if one set is contained within another
    def is_subset(subset, superset):
        """Checks if `subset` is contained within `superset`."""
        return subset.issubset(superset) 
    # issubset() is a method that returns:
    # True if all elements of the subset are present in the superset. 
    # Otherwise, it returns False.
    #######################################################################

    # Sort the sets by size (from smallest to largest)
    Set_of_sets.sort(key=lambda x: len(x))  # Sort `Set_of_sets` by the size of each set

    # Create set names based on their level
    Set_names = [f"S{i+1}" for i in range(len(Set_of_sets))]
    
    # Create a dictionary of labels for the nodes
    labels = {f"S{i+1}": f"S{i+1}" for i in range(len(Set_of_sets))}
    
    # Build the containment graph
    for i, child_set in enumerate(Set_of_sets):  # Iterate over sets as children
        for j, parent_set in enumerate(Set_of_sets):  # Iterate over sets as parents
            if i != j and is_subset(child_set, parent_set):  # If they are not the same and containment exists
                # Add only if there is no reverse edge (avoid indirect edges)
                if Set_names[j] not in containment_graph[Set_names[i]]:  # Check subsets of Set_names[i]
                    containment_graph[Set_names[i]].append(Set_names[j]) # Add the relationship
    
    # Print tuple with labels and subsets
    label_set_pairs = [(label, s) for label, s in zip(Set_names, Set_of_sets)]
    print("Tuples of Labels and Subsets:", label_set_pairs)

    # Create a directed graph
    G = nx.DiGraph()
    
    # Add nodes to the graph
    for name in Set_names:
        G.add_node(name)
    #######################################################################
    # Add nodes and direct containment relationships
    for i, child_set in enumerate(Set_of_sets):
        for j in range(i + 1, len(Set_of_sets)):
            parent_set = Set_of_sets[j]
            if child_set.issubset(parent_set):
                # Check if an intermediate set exists
                is_direct = True
                for k in range(i + 1, j):
                    intermediate_set = Set_of_sets[k]
                    if child_set.issubset(intermediate_set) and intermediate_set.issubset(parent_set):
                        is_direct = False
                        break
                if is_direct:
                    G.add_edge(Set_names[i], Set_names[j])  # Add direct edges only
    #######################################################################        
    # Automatically assign positions by levels in the graph
    levels = {}
    for node in nx.topological_sort(G):  # Topologically sort the nodes
        depth = 0
        for parent in G.predecessors(node):  # Find predecessors of the current node
            depth = max(depth, levels[parent] + 1)  # Define the depth of the current node
        levels[node] = depth
    
    # Arrange nodes by levels and center them on the X-axis
    pos = {}  # Dictionary to store node positions
    nodes_by_level = defaultdict(list)  # Group nodes by containment levels
    
    for node, level in levels.items():
        nodes_by_level[level].append(node)  # Add nodes to the corresponding level
    
    max_nodes_in_level = max(len(nodes) for nodes in nodes_by_level.values())  # Maximum number of nodes at any level
    
    # Do not invert levels on the Y-axis; smaller sets will be at the bottom
    for level, nodes in nodes_by_level.items():
        num_nodes = len(nodes)  # Number of nodes at the current level
        horizontal_spacing = max_nodes_in_level / (num_nodes + 1)  # Horizontal spacing between nodes
        start_x = -((num_nodes - 1) * horizontal_spacing) / 2  # Initial X-coordinate to center nodes
        
        # Calculate vertical position without inverting the Y-axis
        for i, node in enumerate(nodes):
            pos[node] = (start_x + i * horizontal_spacing, level)  # Assign position (X, Y) to each node

    # Draw the graph with labels and arrows, adjusting visual properties
    plt.figure(figsize=(10, 8))
    nx.draw(
        G, pos, labels=labels, arrows=True, arrowstyle='-|>', arrowsize=20, 
        node_size=1000, node_color="skyblue", font_size=12, font_weight="bold", 
        edge_color="gray"
    )  

    #######################################################################
    # Create set labels to display when hovering over nodes
    cursor_labels = [
        ', '.join(sorted(list(s))) if s else '∅'
        for s in Set_of_sets
    ]
    
    # Configure interactive cursor to display set labels
    cursor = mplcursors.cursor(plt.gca(), hover=True)
    
    @cursor.connect("add")
    def on_add(sel):
        """Displays the content of the set when hovering over a node."""
        index = sel.index
        if index < len(cursor_labels):
            sel.annotation.set_text(cursor_labels[index])
            sel.annotation.set_visible(True)

    @cursor.connect("remove")
    def on_remove(sel):
        """Hides the annotation when the cursor leaves the node."""
        sel.annotation.set_visible(False)
    
    plt.show()  # Display the graph

##################################################################

