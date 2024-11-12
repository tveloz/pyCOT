import networkx as nx
import matplotlib.pyplot as plt
from mplcursors import cursor
import numpy as np
import re 

##################################################################################### 
def Hierarchy_visualize(input_data):
    """
    Creates and visualizes a hierarchical graph based on input sets and levels.

    Parameters:
    - input_data (list of tuples): List of tuples containing the name, elements, and level of each node.
    
    Returns:
    - G (networkx.DiGraph): Directed graph with the hierarchical structure based on the input levels.

    Description:
    This function creates a hierarchical directed graph representing inclusion relationships across different levels.
    It assigns nodes and establishes edges based on the provided sets and levels. Nodes at lower levels
    are connected to upper levels if the set of elements is contained within the set of a higher level.
    """

    # Sort sets by the third element
    input_data.sort(key=lambda x: x[2])  # Sort by the third element

    # Create a directed graph
    G = nx.DiGraph()

    # Add nodes and their elements, assigning their level
    for name, elements, level in input_data:
        G.add_node(name, elements=elements, level=level)

    # Set to keep track of nodes with incoming edges
    nodes_with_incoming_edges = set()

    # Automatically define inclusion relationships between levels
    for i in range(len(input_data)):
        for j in range(len(input_data)):
            if i != j:  # Ensure nodes are not compared to themselves
                set_i = set(input_data[i][1])  # Convert the second coordinate of entry 'i' to a set
                set_j = set(input_data[j][1])  # Convert the second coordinate of entry 'j' to a set

                # Calculate level difference
                k = input_data[j][2] - input_data[i][2]
                
                # Add edges from level 0 to level 1 only if the node has no incoming edges
                for node in input_data:
                    if node[2] == k and k > 3:  # Level 1
                        if len(G.predecessors(node[0])) == 0:  # Check if there are no incoming edges
                            G.add_edge(input_data[0][0], node[0])  # Add edge from "S0" to each node at level j
                                
                # Check if 'set_i' is contained within 'set_j' and if k >= 1
                if set_i.issubset(set_j) and k >= 1:
                    # Check if an indirect edge has already been added
                    indirect_edge_exists = False
                    for neighbor in G.successors(input_data[i][0]):
                        if set(input_data[G.nodes[neighbor]['level']][1]).issubset(set_j):
                            indirect_edge_exists = True  # Corrected to reflect an indirect edge exists
                            break 
                    if not indirect_edge_exists:
                        G.add_edge(input_data[i][0], input_data[j][0])  # Add directed edge between 'i' and 'j'
                        nodes_with_incoming_edges.add(input_data[j][0])  # Mark that node j has an incoming edge
                        break 

    # Generate positions
    pos_hierarchy = generate_positions(G)

    # Draw the graph
    plt.figure(figsize=(8, 6))

    # Draw nodes
    nx.draw_networkx_nodes(G, pos=pos_hierarchy, node_shape='s', node_size=350, node_color='lightblue')

    # Draw edges and labels with straight arrows
    nx.draw_networkx_edges(G, pos=pos_hierarchy, arrows=True, arrowstyle='-|>', arrowsize=20, edge_color='gray', width=2,
                         connectionstyle='arc3,rad=0.2') #rad=0.2: changes arrow curvature

    nx.draw_networkx_labels(G, pos=pos_hierarchy, font_size=10, font_weight='bold')

    #####################################################################################
    # Create a cursor to show only set elements when hovering over nodes
    cursor_labels = [
        ', '.join(elements) if elements else '∅'
        for name, elements, _ in input_data
    ]

    cursor_obj = cursor(plt.gca(), hover=True)

    @cursor_obj.connect("add")
    def on_add(sel):
        index = sel.index
        if index < len(cursor_labels):
            sel.annotation.set_text(cursor_labels[index])
            sel.annotation.set_visible(True)

    @cursor_obj.connect("remove")
    def on_remove(sel):
        sel.annotation.set_visible(False)

    plt.title("Hierarchy")
    plt.gca().set_frame_on(False)

    plt.show()
    return G

#######################################################################################
def generate_positions(graph):
    """
    Automatically generates positions for a hierarchical graph.

    Parameters:
    - graph (networkx.DiGraph): Directed graph for which positions will be generated.

    Returns:
    - dict: Position dictionary in the format {node: (x, y)}, where 'y' represents the hierarchical level.
    """
    pos = {}
    levels = {}

    for node in graph.nodes:
        level = graph.nodes[node]['level']
        if level not in levels:
            levels[level] = []
        levels[level].append(node)

    for level, nodes in levels.items():
        y_offset = level

        if len(nodes) == 1:
            x_offset = [0]
        else:
            x_offset = np.linspace(-1, 1, len(nodes))

        for i, node in enumerate(nodes):
            pos[node] = (x_offset[i], y_offset)

    return pos


#######################################################################################
def convert_to_input_data(reactive_sorgs):
    """
    Converts a list of subsets into a format suitable for hierarchical graph visualization.

    Parameters:
    - reactive_sorgs (list of lists): List of subsets, where each subset represents a set of reactive elements.

    Returns:
    - list of lists: List of nodes in the format [name, set of elements, level], where the level is the set cardinality.
    """
    input_data = []

    # Remove duplicates by converting each subset to a tuple, then back to a list
    reactive_sorgs = [list(subset) for subset in set(tuple(sorted(subset)) for subset in reactive_sorgs)]

    # Sort each subset within reactive_sorgs
    reactive_sorgs = [sorted(subset) for subset in reactive_sorgs]

    # Check if the empty set is already in reactive_sorgs
    if [] in reactive_sorgs:
        input_data.append(["S1", [], 0])
        start_index = 1  # Start index for remaining elements from 1

    # Create nodes for each subset
    for i, subset in enumerate(reactive_sorgs, start=1):  # start=1 to start indices from 1
        if subset:  # Ignore empty set if already added
            level = len(subset)  # Calculate level as the subset size
            input_data.append(["", subset, level])  # Assign name 'S1', 'S2', etc. at creation time

    input_data.sort(key=lambda x: (x[2], x[1]))

    # Assign sequential node names (S1, S2, ...)
    for i, node in enumerate(input_data):
        node[0] = f"S{i + 1}"

    return input_data
