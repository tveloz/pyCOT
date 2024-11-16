from collections import defaultdict
import networkx as nx
import matplotlib.pyplot as plt
import mplcursors

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