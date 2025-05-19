from pyvis.network import Network 
import networkx as nx
import matplotlib.pyplot as plt
import mplcursors
import webbrowser  # Allows opening URLs or local files in the system's default browser
import os  # For handling paths and checking file existence
from collections import defaultdict
import sys
sys.stdout.reconfigure(encoding='utf-8')
import tempfile
from pyCOT.closure_structure_b import *
import networkx as nx
from collections import defaultdict

def hierarchy_build(input_data):
    """
    Builds a hierarchy graph where each node contains both its name and the set of species.
    
    Args:
        input_data (list): List of lists, where each element represents a set of species.
        
    Returns:
        nx.DiGraph: A directed graph with node attributes including the set of species.
    """
    # Convert each subset to a set and then to a list to remove duplicates 
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
    labels = {f"S{i+1}": f"S{i+1}" for i in range(len(Set_of_sets))}
    
    # Create a directed graph
    G = nx.DiGraph()
    
    # Add nodes to the graph with attributes
    for i, name in enumerate(Set_names):
        # Store both the name and the set of species in the node attributes
        G.add_node(name, species_set=Set_of_sets[i], label=name)
    
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
    
    # Print information about the node attributes
    print("Nodes with their species sets:")
    for node in G.nodes():
        species_set = G.nodes[node]['species_set']
        print(f"{node}: {sorted(list(species_set)) if species_set else '∅'}")
    
    return G
def get_species_from_node(G, node_name):
    if node_name in G:
        return G.nodes[node_name]['species_set']
    else:
        raise ValueError(f"Node {node_name} not found in the graph")
    
def join(RN, G, node1, node2):
    # Get species sets from the nodes
    if nx.has_path(G,node1,node2):
        return get_species_from_node(G,node2), True, node2
    elif nx.has_path(G,node2,node1):
        return get_species_from_node(G,node1), True, node1
    else:
        species_set1 = G.nodes[node1]['species_set']
        species_set2 = G.nodes[node2]['species_set']

        # Compute the union and its closure
        union_set = list(species_set1.union(species_set2))
        closure_set = closure(RN, union_set)
        
        # Check if semi-self-maintaining
        is_ssm = RN.is_semi_self_maintaining(closure_set)
        
        if is_ssm:
            for node in G.nodes():
                X=get_species_from_node(G,node)
                if set(X)==set(closure_set):
                    return closure_set, is_ssm, node
        else:
            return closure_set, is_ssm, "Not Found"
def meet(RN, G, node1, node2):
    # Get species sets from the nodes
    if nx.has_path(G,node1,node2):
        return get_species_from_node(G,node1), True, node1
    elif nx.has_path(G,node2,node1):
        return get_species_from_node(G,node2), True, node2
    else:
        species_set1 = G.nodes[node1]['species_set']
        species_set2 = G.nodes[node2]['species_set']

        # Compute the union and its closure
        intersection_set = list(species_set1.intersection(species_set2))
        
        # print("intersection_set meet")
        # print(intersection_set)
        reacs=RN.get_reactions_from_species(intersection_set)
        intersection_set=RN.get_supp_from_reactions(reacs)
        closure_set = closure(RN, intersection_set)
        # print("closure = "+str(closure_set))
        # # Check if semi-self-maintaining
        ##get the reactive part

        is_ssm = RN.is_semi_self_maintaining(closure_set)
    
        if is_ssm:
            for node in G.nodes():
                X=get_species_from_node(G,node)
                if set(X)==set(closure_set):
                    return closure_set, is_ssm, node
            return closure_set, is_ssm, "Not Found" 
        else:
            return closure_set, is_ssm, "Not Found"

def find_node_for_species(G, species_set):
    """
    Find the node in G that corresponds to the given species set.
    
    Parameters:
    -----------
    G : nx.DiGraph
        The hierarchy graph with nodes containing species sets
    species_set : list or set
        The species set to find in the graph
        
    Returns:
    --------
    str or None
        The node label if found, None otherwise
    """
    species_set = set(species_set)  # Convert to set for comparison
    
    for node in G.nodes():
        node_species = set(G.nodes[node]['species_set'])
        if node_species == species_set:
            return node
    
    return None
def test_completeness(RN, G):
    """
    Tests whether the organizational structure forms a complete lattice.
    
    A lattice is complete if every subset of organizations has both a join (least upper bound)
    and a meet (greatest lower bound).
    
    Parameters:
    -----------
    RN : ReactionNetwork
        The reaction network object.
    G : networkx.DiGraph
        The directed graph representing the hierarchy of organizations.
        
    Returns:
    --------
    dict
        A dictionary containing:
        - 'is_complete': Boolean indicating if the structure is a complete lattice
        - 'completeness_ratio': Ratio of successful join/meet operations to total possible
        - 'missing_joins': List of pairs of nodes that don't have a valid join
        - 'missing_meets': List of pairs of nodes that don't have a valid meet
        - 'problematic_nodes': Dictionary of nodes ranked by how often they appear in missing operations
        - 'stats': Dictionary with more detailed statistics
    """
    nodes = list(G.nodes())
    n = len(nodes)
    
    # Initialize counters
    total_pairs = n * (n - 1) // 2  # Number of unique pairs
    successful_joins = 0
    successful_meets = 0
    missing_joins = []
    missing_meets = []
    node_missing_join_count = {node: 0 for node in nodes}
    node_missing_meet_count = {node: 0 for node in nodes}
    
    for i in range(n):
        for j in range(i+1, n):
            node1 = nodes[i]
            node2 = nodes[j]
            
            # Test join
            join_result = test_join_exists(RN, G, node1, node2)
            if join_result['exists']:
                successful_joins += 1
            else:
                missing_joins.append((node1, node2, join_result['computed_species']))
                node_missing_join_count[node1] += 1
                node_missing_join_count[node2] += 1
            
            # Test meet
            meet_result = test_meet_exists(RN, G, node1, node2)
            if meet_result['exists']:
                successful_meets += 1
            else:
                missing_meets.append((node1, node2, meet_result['computed_species']))
                node_missing_meet_count[node1] += 1
                node_missing_meet_count[node2] += 1
    
    # Calculate completeness ratio (average of join and meet success rates)
    join_ratio = successful_joins / total_pairs if total_pairs > 0 else 1
    meet_ratio = successful_meets / total_pairs if total_pairs > 0 else 1
    completeness_ratio = (join_ratio + meet_ratio) / 2
    
    # Identify problematic nodes (those that appear most frequently in missing operations)
    total_node_problems = {}
    for node in nodes:
        total_node_problems[node] = node_missing_join_count[node] + node_missing_meet_count[node]
    
    # Sort nodes by how problematic they are
    problematic_nodes = {
        'overall': sorted([(node, count) for node, count in total_node_problems.items()], 
                         key=lambda x: x[1], reverse=True),
        'joins': sorted([(node, count) for node, count in node_missing_join_count.items()], 
                        key=lambda x: x[1], reverse=True),
        'meets': sorted([(node, count) for node, count in node_missing_meet_count.items()], 
                        key=lambda x: x[1], reverse=True)
    }
    
    # Add species information for the most problematic nodes
    top_problematic = problematic_nodes['overall'][:min(5, len(nodes))]
    problematic_with_species = []
    for node, count in top_problematic:
        if count > 0:  # Only include actually problematic nodes
            species = get_species_from_node(G, node)
            problematic_with_species.append((node, count, species))
    
    return {
        'is_complete': len(missing_joins) == 0 and len(missing_meets) == 0,
        'completeness_ratio': completeness_ratio,
        'missing_joins': missing_joins,
        'missing_meets': missing_meets,
        'problematic_nodes': problematic_nodes,
        'top_problematic_nodes': problematic_with_species,
        'stats': {
            'total_pairs': total_pairs,
            'successful_joins': successful_joins,
            'successful_meets': successful_meets,
            'join_ratio': join_ratio,
            'meet_ratio': meet_ratio
        }
    }

def test_join_exists(RN, G, node1, node2):
    """
    Tests whether a join exists for two nodes in the organizational hierarchy.
    
    Parameters:
    -----------
    RN : ReactionNetwork
        The reaction network object.
    G : networkx.DiGraph
        The directed graph representing the hierarchy of organizations.
    node1, node2 : str
        Names of the nodes to test.
        
    Returns:
    --------
    dict
        A dictionary containing:
        - 'exists': Boolean indicating if a valid join exists
        - 'join_node': The node representing the join if it exists, None otherwise
        - 'computed_species': The set of species in the computed join
    """
    try:
        # Get species sets for the nodes
        species1 = get_species_from_node(G, node1)
        species2 = get_species_from_node(G, node2)
        
        # Compute the union of species
        union_species = set(species1).union(set(species2))
        
        # Calculate closure
        closure_species = compute_closure(RN, list(union_species))
        
        # Check if this set is self-maintaining
        is_self_maintaining = RN.is_semi_self_maintaining(list(closure_species))
        
        if not is_self_maintaining:
            return {'exists': False, 'join_node': None, 'computed_species': closure_species}
        
        # Find the node that corresponds to this closure
        join_node = find_node_for_species(G, closure_species)
        
        # Check if a node with exactly these species exists
        if join_node is None:
            return {'exists': False, 'join_node': None, 'computed_species': closure_species}
        
        return {'exists': True, 'join_node': join_node, 'computed_species': closure_species}
    
    except Exception as e:
        print(f"Error testing join for {node1} and {node2}: {e}")
        return {'exists': False, 'join_node': None, 'computed_species': set()}

def test_meet_exists(RN, G, node1, node2):
    """
    Tests whether a meet exists for two nodes in the organizational hierarchy.
    
    Parameters:
    -----------
    RN : ReactionNetwork
        The reaction network object.
    G : networkx.DiGraph
        The directed graph representing the hierarchy of organizations.
    node1, node2 : str
        Names of the nodes to test.
        
    Returns:
    --------
    dict
        A dictionary containing:
        - 'exists': Boolean indicating if a valid meet exists
        - 'meet_node': The node representing the meet if it exists, None otherwise
        - 'computed_species': The set of species in the computed meet
    """
    try:
        # Get species sets for the nodes
        species1 = get_species_from_node(G, node1)
        species2 = get_species_from_node(G, node2)
        
        # Compute the intersection of species
        intersection_species = set(species1).intersection(set(species2))
        
        # For meet, we need to find the largest organization contained in the intersection
        # This is more complex than join and may not exist for all pairs
        
        # First, check if the intersection itself is an organization
        is_closed = RN.is_closed(list(intersection_species))
        is_self_maintaining = RN.is_semi_self_maintaining(list(intersection_species))
        
        if is_closed and is_self_maintaining:
            # The intersection itself is an organization
            meet_node = find_node_for_species(G, intersection_species)
            if meet_node is not None:
                return {'exists': True, 'meet_node': meet_node, 'computed_species': intersection_species}
        
        # If the intersection is not an organization, find the largest organization contained in it
        candidates = []
        for node in G.nodes():
            node_species = get_species_from_node(G, node)
            if set(node_species).issubset(intersection_species):
                candidates.append((node, node_species))
        
        if not candidates:
            return {'exists': False, 'meet_node': None, 'computed_species': intersection_species}
        
        # Find the candidate with the most species (the largest organization in the intersection)
        largest_candidate = max(candidates, key=lambda x: len(x[1]))
        
        return {'exists': True, 'meet_node': largest_candidate[0], 'computed_species': largest_candidate[1]}
    
    except Exception as e:
        print(f"Error testing meet for {node1} and {node2}: {e}")
        return {'exists': False, 'meet_node': None, 'computed_species': set()}

def compute_closure(RN, species_list):
    """
    Computes the closure of a set of species in the reaction network.
    
    Parameters:
    -----------
    RN : ReactionNetwork
        The reaction network object.
    species_list : list
        List of species names.
        
    Returns:
    --------
    set
        The closure of the input species set.
    """
    current_set = set(species_list)
    
    while True:
        # Get reactions that can be triggered by the current set
        reactions = RN.get_reactions_from_species(list(current_set))
        
        # Get products of these reactions
        products = RN.get_prod_from_reactions(reactions)
        
        # Union with the current set
        new_set = current_set.union(products)
        
        # If no new species were added, we've reached closure
        if new_set == current_set:
            break
        
        current_set = new_set
    
    return current_set

def safe_join(RN, G, node1, node2):
    """
    A safe version of the join operation that checks if the join exists first.
    
    Parameters:
    -----------
    RN : ReactionNetwork
        The reaction network object.
    G : networkx.DiGraph
        The directed graph representing the hierarchy of organizations.
    node1, node2 : str
        Names of the nodes to join.
        
    Returns:
    --------
    tuple
        (species_set, success, node)
        - species_set: The set of species in the join
        - success: Boolean indicating if the join was successful
        - node: The node representing the join if it exists, None otherwise
    """
    join_result = test_join_exists(RN, G, node1, node2)
    
    if join_result['exists']:
        return join_result['computed_species'], True, join_result['join_node']
    else:
        # Return the computed species even if no node exists for it
        return join_result['computed_species'], False, None

def safe_meet(RN, G, node1, node2):
    """
    A safe version of the meet operation that checks if the meet exists first.
    
    Parameters:
    -----------
    RN : ReactionNetwork
        The reaction network object.
    G : networkx.DiGraph
        The directed graph representing the hierarchy of organizations.
    node1, node2 : str
        Names of the nodes to meet.
        
    Returns:
    --------
    tuple
        (species_set, success, node)
        - species_set: The set of species in the meet
        - success: Boolean indicating if the meet was successful
        - node: The node representing the meet if it exists, None otherwise
    """
    meet_result = test_meet_exists(RN, G, node1, node2)
    
    if meet_result['exists']:
        return meet_result['computed_species'], True, meet_result['meet_node']
    else:
        # Return the computed species even if no node exists for it
        return meet_result['computed_species'], False, None

def test_distributive_law1_detailed(RN, G, x, y, z):
    """
    Tests the first distributive law: x ∧ (y ∨ z) = (x ∧ y) ∨ (x ∧ z)
    with detailed information about where the law fails.
    
    Parameters:
    -----------
    RN : ReactionNetwork
        The reaction network object.
    G : networkx.DiGraph
        The directed graph representing the hierarchy of organizations.
    x, y, z : str
        Names of the nodes to test.
        
    Returns:
    --------
    dict
        A dictionary containing:
        - 'success': Boolean indicating if the distributive law holds
        - 'operation_failed': Boolean indicating if any operation failed
        - 'failure_step': String describing where the law failed (if it failed)
        - 'left_side': The left side of the equation (x ∧ (y ∨ z))
        - 'right_side': The right side of the equation ((x ∧ y) ∨ (x ∧ z))
        - 'steps': Dictionary with intermediate steps and results
    """
    steps = {}
    operation_failed = False
    failure_step = None
    
    # Compute y ∨ z
    y_join_z = test_join_exists(RN, G, y, z)
    steps['y_join_z'] = y_join_z
    
    if not y_join_z['exists']:
        operation_failed = True
        failure_step = "y ∨ z doesn't exist"
        return {
            'success': False,
            'operation_failed': operation_failed,
            'failure_step': failure_step,
            'left_side': None,
            'right_side': None,
            'steps': steps
        }
    
    # Compute x ∧ (y ∨ z)
    left_side = test_meet_exists(RN, G, x, y_join_z['join_node'])
    steps['x_meet_y_join_z'] = left_side
    
    if not left_side['exists']:
        operation_failed = True
        failure_step = "x ∧ (y ∨ z) doesn't exist"
        return {
            'success': False,
            'operation_failed': operation_failed,
            'failure_step': failure_step,
            'left_side': left_side['computed_species'],
            'right_side': None,
            'steps': steps
        }
    
    # Compute x ∧ y
    x_meet_y = test_meet_exists(RN, G, x, y)
    steps['x_meet_y'] = x_meet_y
    
    if not x_meet_y['exists']:
        operation_failed = True
        failure_step = "x ∧ y doesn't exist"
        return {
            'success': False,
            'operation_failed': operation_failed,
            'failure_step': failure_step,
            'left_side': left_side['computed_species'],
            'right_side': None,
            'steps': steps
        }
    
    # Compute x ∧ z
    x_meet_z = test_meet_exists(RN, G, x, z)
    steps['x_meet_z'] = x_meet_z
    
    if not x_meet_z['exists']:
        operation_failed = True
        failure_step = "x ∧ z doesn't exist"
        return {
            'success': False,
            'operation_failed': operation_failed,
            'failure_step': failure_step,
            'left_side': left_side['computed_species'],
            'right_side': None,
            'steps': steps
        }
    
    # Compute (x ∧ y) ∨ (x ∧ z)
    right_side = test_join_exists(RN, G, x_meet_y['meet_node'], x_meet_z['meet_node'])
    steps['x_meet_y_join_x_meet_z'] = right_side
    
    if not right_side['exists']:
        operation_failed = True
        failure_step = "(x ∧ y) ∨ (x ∧ z) doesn't exist"
        return {
            'success': False,
            'operation_failed': operation_failed,
            'failure_step': failure_step,
            'left_side': left_side['computed_species'],
            'right_side': right_side['computed_species'],
            'steps': steps
        }
    
    # Compare the species sets of left_side and right_side
    left_species = left_side['computed_species']
    right_species = right_side['computed_species']
    
    if left_species == right_species:
        return {
            'success': True,
            'operation_failed': False,
            'failure_step': None,
            'left_side': left_species,
            'right_side': right_species,
            'steps': steps
        }
    else:
        failure_step = "Left side != Right side"
        return {
            'success': False,
            'operation_failed': False,
            'failure_step': failure_step,
            'left_side': left_species,
            'right_side': right_species,
            'steps': steps
        }

def test_distributive_law1(RN, G, x, y, z):
    """
    Tests the first distributive law: x ∧ (y ∨ z) = (x ∧ y) ∨ (x ∧ z)
    
    Parameters:
    -----------
    RN : ReactionNetwork
        The reaction network object.
    G : networkx.DiGraph
        The directed graph representing the hierarchy of organizations.
    x, y, z : str
        Names of the nodes to test.
        
    Returns:
    --------
    bool
        True if the distributive law holds, False otherwise.
    """
    result = test_distributive_law1_detailed(RN, G, x, y, z)
    return result['success']

def test_distributive_law2_detailed(RN, G, x, y, z):
    """
    Tests the second distributive law: x ∨ (y ∧ z) = (x ∨ y) ∧ (x ∨ z)
    with detailed information about where the law fails.
    
    Parameters:
    -----------
    RN : ReactionNetwork
        The reaction network object.
    G : networkx.DiGraph
        The directed graph representing the hierarchy of organizations.
    x, y, z : str
        Names of the nodes to test.
        
    Returns:
    --------
    dict
        A dictionary containing:
        - 'success': Boolean indicating if the distributive law holds
        - 'operation_failed': Boolean indicating if any operation failed
        - 'failure_step': String describing where the law failed (if it failed)
        - 'left_side': The left side of the equation (x ∨ (y ∧ z))
        - 'right_side': The right side of the equation ((x ∨ y) ∧ (x ∨ z))
        - 'steps': Dictionary with intermediate steps and results
    """
    steps = {}
    operation_failed = False
    failure_step = None
    
    # Compute y ∧ z
    y_meet_z = test_meet_exists(RN, G, y, z)
    steps['y_meet_z'] = y_meet_z
    
    if not y_meet_z['exists']:
        operation_failed = True
        failure_step = "y ∧ z doesn't exist"
        return {
            'success': False,
            'operation_failed': operation_failed,
            'failure_step': failure_step,
            'left_side': None,
            'right_side': None,
            'steps': steps
        }
    
    # Compute x ∨ (y ∧ z)
    left_side = test_join_exists(RN, G, x, y_meet_z['meet_node'])
    steps['x_join_y_meet_z'] = left_side
    
    if not left_side['exists']:
        operation_failed = True
        failure_step = "x ∨ (y ∧ z) doesn't exist"
        return {
            'success': False,
            'operation_failed': operation_failed,
            'failure_step': failure_step,
            'left_side': left_side['computed_species'],
            'right_side': None,
            'steps': steps
        }
    
    # Compute x ∨ y
    x_join_y = test_join_exists(RN, G, x, y)
    steps['x_join_y'] = x_join_y
    
    if not x_join_y['exists']:
        operation_failed = True
        failure_step = "x ∨ y doesn't exist"
        return {
            'success': False,
            'operation_failed': operation_failed,
            'failure_step': failure_step,
            'left_side': left_side['computed_species'],
            'right_side': None,
            'steps': steps
        }
    
    # Compute x ∨ z
    x_join_z = test_join_exists(RN, G, x, z)
    steps['x_join_z'] = x_join_z
    
    if not x_join_z['exists']:
        operation_failed = True
        failure_step = "x ∨ z doesn't exist"
        return {
            'success': False,
            'operation_failed': operation_failed,
            'failure_step': failure_step,
            'left_side': left_side['computed_species'],
            'right_side': None,
            'steps': steps
        }
    
    # Compute (x ∨ y) ∧ (x ∨ z)
    right_side = test_meet_exists(RN, G, x_join_y['join_node'], x_join_z['join_node'])
    steps['x_join_y_meet_x_join_z'] = right_side
    
    if not right_side['exists']:
        operation_failed = True
        failure_step = "(x ∨ y) ∧ (x ∨ z) doesn't exist"
        return {
            'success': False,
            'operation_failed': operation_failed,
            'failure_step': failure_step,
            'left_side': left_side['computed_species'],
            'right_side': right_side['computed_species'],
            'steps': steps
        }
    
    # Compare the species sets of left_side and right_side
    left_species = left_side['computed_species']
    right_species = right_side['computed_species']
    
    if left_species == right_species:
        return {
            'success': True,
            'operation_failed': False,
            'failure_step': None,
            'left_side': left_species,
            'right_side': right_species,
            'steps': steps
        }
    else:
        failure_step = "Left side != Right side"
        return {
            'success': False,
            'operation_failed': False,
            'failure_step': failure_step,
            'left_side': left_species,
            'right_side': right_species,
            'steps': steps
        }

def test_distributive_law2(RN, G, x, y, z):
    """
    Tests the second distributive law: x ∨ (y ∧ z) = (x ∨ y) ∧ (x ∨ z)
    
    Parameters:
    -----------
    RN : ReactionNetwork
        The reaction network object.
    G : networkx.DiGraph
        The directed graph representing the hierarchy of organizations.
    x, y, z : str
        Names of the nodes to test.
        
    Returns:
    --------
    bool
        True if the distributive law holds, False otherwise.
    """
    result = test_distributive_law2_detailed(RN, G, x, y, z)
    return result['success']

def test_distributed_joins_and_meets(RN, G, node_triples=None):
    """
    Tests the distributive property for joins and meets in the organizational hierarchy.
    
    For a lattice to be distributive, for all elements x, y, z:
    x ∧ (y ∨ z) = (x ∧ y) ∨ (x ∧ z)
    x ∨ (y ∧ z) = (x ∨ y) ∧ (x ∨ z)
    
    Parameters:
    -----------
    RN : ReactionNetwork
        The reaction network object.
    G : networkx.DiGraph
        The directed graph representing the hierarchy of organizations.
    node_triples : list, optional
        Specific triples of nodes to test. If None, all possible triples will be tested.
        
    Returns:
    --------
    dict
        A dictionary containing:
        - 'is_distributive': Boolean indicating if the structure is distributive
        - 'distributivity_ratio': Ratio of distributive triples to total triples
        - 'non_distributive_triples': List of triples that violate distributivity with details
        - 'problematic_nodes': Dictionary of nodes ranked by how often they appear in non-distributive triples
        - 'failure_types': Dictionary counting different types of failures
        - 'stats': Dictionary with more detailed statistics
    """
    nodes = list(G.nodes())
    
    if node_triples is None:
        # Generate all possible triples
        node_triples = []
        for i in range(len(nodes)):
            for j in range(len(nodes)):
                for k in range(len(nodes)):
                    if i != j and i != k and j != k:
                        node_triples.append((nodes[i], nodes[j], nodes[k]))
    
    total_triples = len(node_triples)
    distributive_count = 0
    dist1_success = 0
    dist2_success = 0
    non_distributive_triples = []
    node_failure_count = {node: 0 for node in nodes}
    failure_types = {
        "First law failed": 0,
        "Second law failed": 0,
        "Both laws failed": 0,
        "Operation failed": 0,
        "Error": 0
    }
    
    for x, y, z in node_triples:
        try:
            # Detailed results for each law
            dist1_result = test_distributive_law1_detailed(RN, G, x, y, z)
            dist2_result = test_distributive_law2_detailed(RN, G, x, y, z)
            
            # Track failures
            if not dist1_result['success'] and not dist2_result['success']:
                if dist1_result['operation_failed'] or dist2_result['operation_failed']:
                    failure_types["Operation failed"] += 1
                else:
                    failure_types["Both laws failed"] += 1
                non_distributive_triples.append((x, y, z, "Both laws failed", 
                                               {"law1": dist1_result, "law2": dist2_result}))
                node_failure_count[x] += 1
                node_failure_count[y] += 1
                node_failure_count[z] += 1
            elif not dist1_result['success']:
                if dist1_result['operation_failed']:
                    failure_types["Operation failed"] += 1
                else:
                    failure_types["First law failed"] += 1
                non_distributive_triples.append((x, y, z, "First law failed", dist1_result))
                node_failure_count[x] += 1
                node_failure_count[y] += 1
                node_failure_count[z] += 1
                dist2_success += 1
            elif not dist2_result['success']:
                if dist2_result['operation_failed']:
                    failure_types["Operation failed"] += 1
                else:
                    failure_types["Second law failed"] += 1
                non_distributive_triples.append((x, y, z, "Second law failed", dist2_result))
                node_failure_count[x] += 1
                node_failure_count[y] += 1
                node_failure_count[z] += 1
                dist1_success += 1
            else:
                distributive_count += 1
                dist1_success += 1
                dist2_success += 1
        
        except Exception as e:
            failure_types["Error"] += 1
            print(f"Error testing distributivity for {x}, {y}, {z}: {e}")
            non_distributive_triples.append((x, y, z, f"Error: {e}", None))
            node_failure_count[x] += 1
            node_failure_count[y] += 1
            node_failure_count[z] += 1
    
    distributivity_ratio = distributive_count / total_triples if total_triples > 0 else 1
    
    # Sort nodes by how problematic they are
    problematic_nodes = sorted([(node, count) for node, count in node_failure_count.items()], 
                              key=lambda x: x[1], reverse=True)
    
    # Add species information for the most problematic nodes
    top_problematic = problematic_nodes[:min(5, len(nodes))]
    problematic_with_species = []
    for node, count in top_problematic:
        if count > 0:  # Only include actually problematic nodes
            species = get_species_from_node(G, node)
            problematic_with_species.append((node, count, species))
    
    return {
        'is_distributive': len(non_distributive_triples) == 0,
        'distributivity_ratio': distributivity_ratio,
        'non_distributive_triples': non_distributive_triples,
        'problematic_nodes': problematic_nodes,
        'top_problematic_nodes': problematic_with_species,
        'failure_types': failure_types,
        'stats': {
            'total_triples': total_triples,
            'distributive_count': distributive_count,
            'dist1_success': dist1_success,
            'dist2_success': dist2_success,
            'dist1_ratio': dist1_success / total_triples if total_triples > 0 else 1,
            'dist2_ratio': dist2_success / total_triples if total_triples > 0 else 1
        }
    }

def safe_test_modularity(RN, G):
    """
    Tests whether the organizational structure satisfies the modular law.
    
    For a lattice to be modular, for all elements x, y, z:
    If x ≤ z, then x ∨ (y ∧ z) = (x ∨ y) ∧ z
    
    Parameters:
    -----------
    RN : ReactionNetwork
        The reaction network object.
    G : networkx.DiGraph
        The directed graph representing the hierarchy of organizations.
        
    Returns:
    --------
    dict
        A dictionary containing:
        - 'is_modular': Boolean indicating if the structure is modular
        - 'modularity_ratio': Ratio of modular triples to total applicable triples
        - 'non_modular_triples': List of triples that violate modularity
        - 'problematic_nodes': Dictionary of nodes ranked by how often they appear in non-modular triples
        - 'failure_types': Dictionary counting different types of failures
    """
    nodes = list(G.nodes())
    total_applicable = 0
    modular_count = 0
    non_modular_triples = []
    node_failure_count = {node: 0 for node in nodes}
    failure_types = {
        "y∧z doesn't exist": 0,
        "x∨(y∧z) doesn't exist": 0,
        "x∨y doesn't exist": 0,
        "(x∨y)∧z doesn't exist": 0,
        "Equation doesn't hold": 0,
        "Error": 0
    }
    
    for x in nodes:
        for y in nodes:
            for z in nodes:
                if x == y or x == z or y == z:
                    continue
                    
                # Check if x ≤ z (x is contained in z)
                x_species = get_species_from_node(G, x)
                z_species = get_species_from_node(G, z)
                
                if not set(x_species).issubset(set(z_species)):
                    continue
                
                total_applicable += 1
                
                try:
                    # Compute y ∧ z
                    y_meet_z_species, y_meet_z_success, y_meet_z_node = safe_meet(RN, G, y, z)
                    
                    if not y_meet_z_success:
                        failure_types["y∧z doesn't exist"] += 1
                        non_modular_triples.append((x, y, z, "y∧z doesn't exist", y_meet_z_species))
                        node_failure_count[y] += 1
                        node_failure_count[z] += 1
                        continue
                    
                    # Compute x ∨ (y ∧ z)
                    left_species, left_success, left_node = safe_join(RN, G, x, y_meet_z_node)
                    
                    if not left_success:
                        failure_types["x∨(y∧z) doesn't exist"] += 1
                        non_modular_triples.append((x, y, z, "x∨(y∧z) doesn't exist", left_species))
                        node_failure_count[x] += 1
                        node_failure_count[y_meet_z_node] += 1
                        continue
                    
                    # Compute x ∨ y
                    x_join_y_species, x_join_y_success, x_join_y_node = safe_join(RN, G, x, y)
                    
                    if not x_join_y_success:
                        failure_types["x∨y doesn't exist"] += 1
                        non_modular_triples.append((x, y, z, "x∨y doesn't exist", x_join_y_species))
                        node_failure_count[x] += 1
                        node_failure_count[y] += 1
                        continue
                    
                    # Compute (x ∨ y) ∧ z
                    right_species, right_success, right_node = safe_meet(RN, G, x_join_y_node, z)
                    
                    if not right_success:
                        failure_types["(x∨y)∧z doesn't exist"] += 1
                        non_modular_triples.append((x, y, z, "(x∨y)∧z doesn't exist", right_species))
                        node_failure_count[x_join_y_node] += 1
                        node_failure_count[z] += 1
                        continue
                    
                    # Check if x ∨ (y ∧ z) = (x ∨ y) ∧ z
                    if set(left_species) == set(right_species):
                        modular_count += 1
                    else:
                        failure_types["Equation doesn't hold"] += 1
                        non_modular_triples.append((x, y, z, "Equation doesn't hold", 
                                                  {"left": left_species, "right": right_species}))
                        node_failure_count[x] += 1
                        node_failure_count[y] += 1
                        node_failure_count[z] += 1
                
                except Exception as e:
                    failure_types["Error"] += 1
                    print(f"Error in modularity test for {x}, {y}, {z}: {e}")
                    non_modular_triples.append((x, y, z, f"Error: {e}", None))
                    node_failure_count[x] += 1
                    node_failure_count[y] += 1
                    node_failure_count[z] += 1
    
    modularity_ratio = modular_count / total_applicable if total_applicable > 0 else 1
    
    # Sort nodes by how problematic they are
    problematic_nodes = sorted([(node, count) for node, count in node_failure_count.items()], 
                              key=lambda x: x[1], reverse=True)
    
    # Add species information for the most problematic nodes
    top_problematic = problematic_nodes[:min(5, len(nodes))]
    problematic_with_species = []
    for node, count in top_problematic:
        if count > 0:  # Only include actually problematic nodes
            species = get_species_from_node(G, node)
            problematic_with_species.append((node, count, species))
    
    return {
        'is_modular': len(non_modular_triples) == 0,
        'modularity_ratio': modularity_ratio,
        'non_modular_triples': non_modular_triples,
        'problematic_nodes': problematic_nodes,
        'top_problematic_nodes': problematic_with_species,
        'failure_types': failure_types,
        'stats': {
            'total_applicable': total_applicable,
            'modular_count': modular_count
        }
    }

def safe_test_distributivity(RN, G):
    """
    Tests whether the organizational structure satisfies the distributive laws.
    
    For a lattice to be distributive, for all elements x, y, z:
    x ∧ (y ∨ z) = (x ∧ y) ∨ (x ∧ z)
    
    Parameters:
    -----------
    RN : ReactionNetwork
        The reaction network object.
    G : networkx.DiGraph
        The directed graph representing the hierarchy of organizations.
        
    Returns:
    --------
    dict
        A dictionary containing:
        - 'is_distributive': Boolean indicating if the structure is distributive
        - 'distributivity_ratio': Ratio of distributive triples to total triples
        - 'non_distributive_triples': List of triples that violate distributivity with details
        - 'problematic_nodes': List of nodes ranked by how often they appear in non-distributive triples
        - 'top_problematic_nodes': List of most problematic nodes with their species
        - 'failure_types': Dictionary counting different types of failures
    """
    return test_distributed_joins_and_meets(RN, G)

def find_node_for_species(G, species_set):
    """
    Find a node in the graph that corresponds to a given set of species.
    
    Parameters:
    -----------
    G : networkx.DiGraph
        The directed graph representing the hierarchy of organizations.
    species_set : set or list
        The set of species to find a node for.
    
    Returns:
    --------
    str or None
        The name of the node if found, None otherwise.
    """
    species_set = set(species_set)
    for node in G.nodes():
        node_species = set(get_species_from_node(G, node))
        if node_species == species_set:
            return node
    return None

# Add user-friendly reporting functions
def print_completeness_report(completeness_result):
    """
    Prints a user-friendly report of the completeness analysis.
    
    Parameters:
    -----------
    completeness_result : dict
        The result from test_completeness function.
    """
    print("\n===== COMPLETENESS ANALYSIS =====")
    print(f"Is the structure complete? {'Yes' if completeness_result['is_complete'] else 'No'}")
    print(f"Completeness ratio: {completeness_result['completeness_ratio']:.2f}")
    print(f"Join success ratio: {completeness_result['stats']['join_ratio']:.2f}")
    print(f"Meet success ratio: {completeness_result['stats']['meet_ratio']:.2f}")
    
    if not completeness_result['is_complete']:
        print("\nMost problematic nodes:")
        for node, count, species in completeness_result['top_problematic_nodes']:
            if count > 0:
                print(f"- Node {node}: involved in {count} failures")
                print(f"  Species: {species}")
        
        print(f"\nMissing joins: {len(completeness_result['missing_joins'])}")
        if len(completeness_result['missing_joins']) > 0:
            print("Sample of missing joins:")
            for i, (node1, node2, species) in enumerate(completeness_result['missing_joins'][:3]):
                print(f"- {node1} ∨ {node2} -> computed species: {species}")
        
        print(f"\nMissing meets: {len(completeness_result['missing_meets'])}")
        if len(completeness_result['missing_meets']) > 0:
            print("Sample of missing meets:")
            for i, (node1, node2, species) in enumerate(completeness_result['missing_meets'][:3]):
                print(f"- {node1} ∧ {node2} -> computed species: {species}")

def print_modularity_report(modularity_result):
    """
    Prints a user-friendly report of the modularity analysis.
    
    Parameters:
    -----------
    modularity_result : dict
        The result from safe_test_modularity function.
    """
    print("\n===== MODULARITY ANALYSIS =====")
    print(f"Is the structure modular? {'Yes' if modularity_result['is_modular'] else 'No'}")
    print(f"Modularity ratio: {modularity_result['modularity_ratio']:.2f}")
    print(f"Total applicable triples: {modularity_result['stats']['total_applicable']}")
    print(f"Modular triples: {modularity_result['stats']['modular_count']}")
    
    if not modularity_result['is_modular']:
        print("\nFailure type distribution:")
        for failure_type, count in modularity_result['failure_types'].items():
            if count > 0:
                print(f"- {failure_type}: {count} instances")
        
        print("\nMost problematic nodes:")
        for node, count, species in modularity_result['top_problematic_nodes']:
            if count > 0:
                print(f"- Node {node}: involved in {count} failures")
                print(f"  Species: {species}")
        
        print(f"\nNon-modular triples: {len(modularity_result['non_modular_triples'])}")
        if len(modularity_result['non_modular_triples']) > 0:
            print("Sample of non-modular triples:")
            for i, triple in enumerate(modularity_result['non_modular_triples'][:3]):
                x, y, z, failure_type, details = triple
                print(f"- ({x}, {y}, {z}): {failure_type}")

def print_distributivity_report(distributivity_result):
    """
    Prints a user-friendly report of the distributivity analysis.
    
    Parameters:
    -----------
    distributivity_result : dict
        The result from safe_test_distributivity function.
    """
    print("\n===== DISTRIBUTIVITY ANALYSIS =====")
    print(f"Is the structure distributive? {'Yes' if distributivity_result['is_distributive'] else 'No'}")
    print(f"Distributivity ratio: {distributivity_result['distributivity_ratio']:.2f}")
    print(f"First law success ratio: {distributivity_result['stats']['dist1_ratio']:.2f}")
    print(f"Second law success ratio: {distributivity_result['stats']['dist2_ratio']:.2f}")
    
    if not distributivity_result['is_distributive']:
        print("\nFailure type distribution:")
        for failure_type, count in distributivity_result['failure_types'].items():
            if count > 0:
                print(f"- {failure_type}: {count} instances")
        
        print("\nMost problematic nodes:")
        for node, count, species in distributivity_result['top_problematic_nodes']:
            if count > 0:
                print(f"- Node {node}: involved in {count} failures")
                print(f"  Species: {species}")
        
        print(f"\nNon-distributive triples: {len(distributivity_result['non_distributive_triples'])}")
        if len(distributivity_result['non_distributive_triples']) > 0:
            print("Sample of non-distributive triples:")
            for i, triple in enumerate(distributivity_result['non_distributive_triples'][:3]):
                x, y, z, failure_type, details = triple
                print(f"- ({x}, {y}, {z}): {failure_type}")
                
                # Print more detailed information for the first example
                if i == 0 and details is not None:
                    if failure_type == "First law failed":
                        left = details.get('left_side')
                        right = details.get('right_side')
                        print(f"  Left side (x ∧ (y ∨ z)): {left}")
                        print(f"  Right side ((x ∧ y) ∨ (x ∧ z)): {right}")
                    elif failure_type == "Second law failed":
                        left = details.get('left_side')
                        right = details.get('right_side')
                        print(f"  Left side (x ∨ (y ∧ z)): {left}")
                        print(f"  Right side ((x ∨ y) ∧ (x ∨ z)): {right}")
                    elif failure_type == "Both laws failed":
                        if 'law1' in details and 'law2' in details:
                            law1 = details['law1']
                            print(f"  First law: {law1.get('failure_step')}")
                            law2 = details['law2']
                            print(f"  Second law: {law2.get('failure_step')}")

def analyze_lattice_properties(RN, G):
    """
    Comprehensive analysis of lattice properties for the organizational structure.
    
    Parameters:
    -----------
    RN : ReactionNetwork
        The reaction network object.
    G : networkx.DiGraph
        The directed graph representing the hierarchy of organizations.
        
    Returns:
    --------
    dict
        A dictionary containing all analysis results.
    """
    print("\n===== LATTICE PROPERTIES ANALYSIS =====")
    
    # Check completeness first
    completeness_result = test_completeness(RN, G)
    print_completeness_report(completeness_result)
    
    # Check modularity
    modularity_result = safe_test_modularity(RN, G)
    print_modularity_report(modularity_result)
    
    # Check distributivity
    distributivity_result = safe_test_distributivity(RN, G)
    print_distributivity_report(distributivity_result)
    
    return {
        'completeness': completeness_result,
        'modularity': modularity_result,
        'distributivity': distributivity_result
    }
def identify_top_bottom_elements(G):
    """
    Identifies the top and bottom elements of the lattice, if they exist.
    
    A top element (⊤) is an element that contains all other elements (is a superset of all other elements).
    A bottom element (⊥) is an element that is contained in all other elements (is a subset of all other elements).
    
    Parameters:
    -----------
    G : networkx.DiGraph
        The directed graph representing the hierarchy of organizations.
        
    Returns:
    --------
    dict
        A dictionary containing:
        - 'has_top': Boolean indicating if a top element exists
        - 'has_bottom': Boolean indicating if a bottom element exists
        - 'top_element': The node name of the top element if it exists, None otherwise
        - 'bottom_element': The node name of the bottom element if it exists, None otherwise
        - 'missing_relations': List of missing containment relations for incomplete lattices
    """
    nodes = list(G.nodes())
    
    # Find candidates for top and bottom elements
    top_candidates = []
    bottom_candidates = []
    
    # Check for nodes with no outgoing edges (potential top elements)
    for node in nodes:
        if G.out_degree(node) == 0:
            top_candidates.append(node)
    
    # Check for nodes with no incoming edges (potential bottom elements)
    for node in nodes:
        if G.in_degree(node) == 0:
            bottom_candidates.append(node)
    
    # For a true top element, it must be reachable from all other nodes
    true_top = None
    for candidate in top_candidates:
        is_true_top = True
        missing = []
        for node in nodes:
            if node != candidate and not nx.has_path(G, node, candidate):
                is_true_top = False
                missing.append((node, candidate))
        if is_true_top:
            true_top = candidate
            break
    
    # For a true bottom element, all other nodes must be reachable from it
    true_bottom = None
    for candidate in bottom_candidates:
        is_true_bottom = True
        missing = []
        for node in nodes:
            if node != candidate and not nx.has_path(G, candidate, node):
                is_true_bottom = False
                missing.append((candidate, node))
        if is_true_bottom:
            true_bottom = candidate
            break
    
    # Collect missing relations for incomplete lattices
    missing_relations = []
    if true_top is None and len(top_candidates) > 0:
        # Find out which nodes aren't connected to our top candidates
        for candidate in top_candidates:
            for node in nodes:
                if node != candidate and not nx.has_path(G, node, candidate):
                    missing_relations.append((node, candidate))
    
    if true_bottom is None and len(bottom_candidates) > 0:
        # Find out which nodes aren't connected from our bottom candidates
        for candidate in bottom_candidates:
            for node in nodes:
                if node != candidate and not nx.has_path(G, candidate, node):
                    missing_relations.append((candidate, node))
    
    return {
        'has_top': true_top is not None,
        'has_bottom': true_bottom is not None,
        'top_element': true_top,
        'bottom_element': true_bottom,
        'top_candidates': top_candidates,
        'bottom_candidates': bottom_candidates,
        'missing_relations': missing_relations
    }

def print_top_bottom_report(top_bottom_result, G):
    """
    Prints a user-friendly report of the top and bottom element analysis.
    
    Parameters:
    -----------
    top_bottom_result : dict
        The result from identify_top_bottom_elements function.
    G : networkx.DiGraph
        The directed graph representing the hierarchy of organizations.
    """
    print("\n===== TOP AND BOTTOM ELEMENT ANALYSIS =====")
    
    # Report on top element
    if top_bottom_result['has_top']:
        top_element = top_bottom_result['top_element']
        top_species = get_species_from_node(G, top_element)
        print(f"Top element (⊤): {top_element}")
        print(f"Top element species: {sorted(list(top_species))}")
    else:
        print("No top element (⊤) found")
        if top_bottom_result['top_candidates']:
            print("Potential top candidates:")
            for candidate in top_bottom_result['top_candidates']:
                species = get_species_from_node(G, candidate)
                print(f"- {candidate}: {sorted(list(species))}")
    
    # Report on bottom element
    if top_bottom_result['has_bottom']:
        bottom_element = top_bottom_result['bottom_element']
        bottom_species = get_species_from_node(G, bottom_element)
        print(f"Bottom element (⊥): {bottom_element}")
        print(f"Bottom element species: {sorted(list(bottom_species))}")
    else:
        print("No bottom element (⊥) found")
        if top_bottom_result['bottom_candidates']:
            print("Potential bottom candidates:")
            for candidate in top_bottom_result['bottom_candidates']:
                species = get_species_from_node(G, candidate)
                print(f"- {candidate}: {sorted(list(species))}")
    
    # Report on overall completeness
    is_bounded_lattice = top_bottom_result['has_top'] and top_bottom_result['has_bottom']
    print(f"\nIs a bounded lattice: {is_bounded_lattice}")
    
    # If not a bounded lattice, report on what's missing
    if not is_bounded_lattice:
        print("\nMissing relations for completeness:")
        if top_bottom_result['missing_relations']:
            print(f"Number of missing relations: {len(top_bottom_result['missing_relations'])}")
            if len(top_bottom_result['missing_relations']) <= 10:
                for source, target in top_bottom_result['missing_relations']:
                    source_species = get_species_from_node(G, source)
                    target_species = get_species_from_node(G, target)
                    print(f"- Missing path from {source} to {target}")
                    print(f"  Source species: {sorted(list(source_species))}")
                    print(f"  Target species: {sorted(list(target_species))}")
            else:
                print("Too many missing relations to display")

def propose_negation_candidates(RN, G, node, top_element=None, bottom_element=None):
    """
    Proposes candidate nodes that could serve as negations for a given node.
    
    Parameters:
    -----------
    RN : ReactionNetwork
        The reaction network object.
    G : networkx.DiGraph
        The directed graph representing the hierarchy of organizations.
    node : str
        The name of the node to find negation candidates for.
    top_element : str, optional
        The name of the top element of the lattice, if known.
    bottom_element : str, optional
        The name of the bottom element of the lattice, if known.
        
    Returns:
    --------
    dict
        A dictionary containing:
        - 'orthogonal_complement': Set of possible orthogonal complements (no common species)
        - 'relative_complement': Set of possible relative complements (join = top, meet = bottom)
        - 'pseudo_complement': Set of possible pseudo-complements (meet = bottom)
        - 'score': Dictionary scoring how well each candidate satisfies negation axioms
    """
    # Get the species for the node
    node_species = get_species_from_node(G, node)
    nodes = list(G.nodes())
    
    # Find top and bottom if not provided
    if top_element is None or bottom_element is None:
        top_bottom = identify_top_bottom_elements(G)
        top_element = top_bottom['top_element']
        bottom_element = top_bottom['bottom_element']
    
    # Get species for top and bottom if they exist
    top_species = set()
    if top_element:
        top_species = get_species_from_node(G, top_element)
    
    bottom_species = set()
    if bottom_element:
        bottom_species = get_species_from_node(G, bottom_element)
    
    # Initialize result containers
    orthogonal_complements = []
    relative_complements = []
    pseudo_complements = []
    
    candidate_scores = {}
    
    for candidate in nodes:
        if candidate == node:
            continue
        
        candidate_species = get_species_from_node(G, candidate)
        
        # Score tracking
        score = {'orthogonality': 0, 'double_negation': 0, 'complement': 0, 'demorgan': 0}
        
        # Check orthogonality: no common species
        if not node_species.intersection(candidate_species):
            orthogonal_complements.append(candidate)
            score['orthogonality'] = 1
        
        # Check relative complement property
        if top_element and bottom_element:
            # Join (node ∨ candidate) should be top
            join_species, join_success, join_node = safe_join(RN, G, node, candidate)
            if join_success and join_node == top_element:
                # Meet (node ∧ candidate) should be bottom
                meet_species, meet_success, meet_node = safe_meet(RN, G, node, candidate)
                if meet_success and meet_node == bottom_element:
                    relative_complements.append(candidate)
                    score['complement'] = 1
        
        # Check pseudo-complement property
        if bottom_element:
            # Meet (node ∧ candidate) should be bottom
            meet_species, meet_success, meet_node = safe_meet(RN, G, node, candidate)
            if meet_success and meet_node == bottom_element:
                pseudo_complements.append(candidate)
                score['complement'] += 0.5
        
        # Check double negation property (need to see if candidate's negation is node)
        # This is harder to check precisely, but we can heuristically score it
        # based on similarity of species sets
        
        # De Morgan's laws would require checking many combinations
        # For now, we just compute an overall score
        candidate_scores[candidate] = score
    
    # Rank candidates by combined score
    ranked_candidates = []
    for candidate in nodes:
        if candidate == node:
            continue
        
        if candidate in candidate_scores:
            score = candidate_scores[candidate]
            total_score = sum(score.values())
            ranked_candidates.append((candidate, total_score, score))
    
    ranked_candidates.sort(key=lambda x: x[1], reverse=True)
    
    return {
        'orthogonal_complements': orthogonal_complements,
        'relative_complements': relative_complements,
        'pseudo_complements': pseudo_complements,
        'ranked_candidates': ranked_candidates,
        'scores': candidate_scores
    }

def print_negation_report(negation_results, G, node):
    """
    Prints a user-friendly report of the negation analysis.
    
    Parameters:
    -----------
    negation_results : dict
        The result from propose_negation_candidates function.
    G : networkx.DiGraph
        The directed graph representing the hierarchy of organizations.
    node : str
        The node for which negation candidates were found.
    """
    node_species = get_species_from_node(G, node)
    
    print(f"\n===== NEGATION ANALYSIS FOR {node} =====")
    print(f"Species: {sorted(list(node_species))}")
    
    print("\nTop ranked negation candidates:")
    
    for i, (candidate, total_score, score_breakdown) in enumerate(negation_results['ranked_candidates'][:5]):
        if i >= 5:  # Limit to top 5
            break
        
        candidate_species = get_species_from_node(G, candidate)
        print(f"\n{i+1}. {candidate} (Score: {total_score:.2f})")
        print(f"   Species: {sorted(list(candidate_species))}")
        print(f"   Score breakdown: {score_breakdown}")
        
        # Check if it's in special categories
        categories = []
        if candidate in negation_results['orthogonal_complements']:
            categories.append("Orthogonal complement")
        if candidate in negation_results['relative_complements']:
            categories.append("Relative complement")
        if candidate in negation_results['pseudo_complements']:
            categories.append("Pseudo-complement")
        
        if categories:
            print(f"   Categories: {', '.join(categories)}")
    
    print(f"\nTotal orthogonal complements: {len(negation_results['orthogonal_complements'])}")
    print(f"Total relative complements: {len(negation_results['relative_complements'])}")
    print(f"Total pseudo-complements: {len(negation_results['pseudo_complements'])}")

def test_demorgan_laws(RN, G, node1, node2, negation_map):
    """
    Tests De Morgan's laws using the provided negation mapping.
    
    De Morgan's laws:
    1. ¬(A ∧ B) = ¬A ∨ ¬B
    2. ¬(A ∨ B) = ¬A ∧ ¬B
    
    Parameters:
    -----------
    RN : ReactionNetwork
        The reaction network object.
    G : networkx.DiGraph
        The directed graph representing the hierarchy of organizations.
    node1, node2 : str
        Names of the nodes to test.
    negation_map : dict
        A dictionary mapping each node to its negation.
        
    Returns:
    --------
    dict
        A dictionary containing:
        - 'law1_success': Boolean indicating if the first law holds
        - 'law2_success': Boolean indicating if the second law holds
        - 'details': Dictionary with detailed results of each test
    """
    results = {
        'law1_success': False,
        'law2_success': False,
        'details': {}
    }
    
    # Check if negations exist for both nodes
    if node1 not in negation_map or node2 not in negation_map:
        results['details']['error'] = "Negation not defined for one or both nodes"
        return results
    
    neg_node1 = negation_map[node1]
    neg_node2 = negation_map[node2]
    
    # Test Law 1: ¬(A ∧ B) = ¬A ∨ ¬B
    
    # Compute A ∧ B
    meet_result = test_meet_exists(RN, G, node1, node2)
    results['details']['meet'] = meet_result
    
    if not meet_result['exists']:
        results['details']['law1_error'] = "Meet doesn't exist"
    else:
        meet_node = meet_result['meet_node']
        
        # Compute ¬(A ∧ B)
        if meet_node not in negation_map:
            results['details']['law1_error'] = "Negation not defined for meet"
        else:
            neg_meet = negation_map[meet_node]
            
            # Compute ¬A ∨ ¬B
            join_result = test_join_exists(RN, G, neg_node1, neg_node2)
            results['details']['neg_join'] = join_result
            
            if not join_result['exists']:
                results['details']['law1_error'] = "Join of negations doesn't exist"
            else:
                neg_join_node = join_result['join_node']
                
                # Check if ¬(A ∧ B) = ¬A ∨ ¬B
                left_species = get_species_from_node(G, neg_meet)
                right_species = get_species_from_node(G, neg_join_node)
                
                if set(left_species) == set(right_species):
                    results['law1_success'] = True
                else:
                    results['details']['law1_error'] = "Left side != Right side"
                    results['details']['law1_left'] = left_species
                    results['details']['law1_right'] = right_species
    
    # Test Law 2: ¬(A ∨ B) = ¬A ∧ ¬B
    
    # Compute A ∨ B
    join_result = test_join_exists(RN, G, node1, node2)
    results['details']['join'] = join_result
    
    if not join_result['exists']:
        results['details']['law2_error'] = "Join doesn't exist"
    else:
        join_node = join_result['join_node']
        
        # Compute ¬(A ∨ B)
        if join_node not in negation_map:
            results['details']['law2_error'] = "Negation not defined for join"
        else:
            neg_join = negation_map[join_node]
            
            # Compute ¬A ∧ ¬B
            meet_result = test_meet_exists(RN, G, neg_node1, neg_node2)
            results['details']['neg_meet'] = meet_result
            
            if not meet_result['exists']:
                results['details']['law2_error'] = "Meet of negations doesn't exist"
            else:
                neg_meet_node = meet_result['meet_node']
                
                # Check if ¬(A ∨ B) = ¬A ∧ ¬B
                left_species = get_species_from_node(G, neg_join)
                right_species = get_species_from_node(G, neg_meet_node)
                
                if set(left_species) == set(right_species):
                    results['law2_success'] = True
                else:
                    results['details']['law2_error'] = "Left side != Right side"
                    results['details']['law2_left'] = left_species
                    results['details']['law2_right'] = right_species
    
    return results

def generate_best_negation_map(RN, G):
    """
    Generates the best possible negation map for the given lattice.
    
    Parameters:
    -----------
    RN : ReactionNetwork
        The reaction network object.
    G : networkx.DiGraph
        The directed graph representing the hierarchy of organizations.
        
    Returns:
    --------
    tuple
        (negation_map, score_info)
        - negation_map: Dictionary mapping each node to its best negation candidate
        - score_info: Information about the quality of the negation map
    """
    # First, identify top and bottom elements
    top_bottom = identify_top_bottom_elements(G)
    top_element = top_bottom['top_element']
    bottom_element = top_bottom['bottom_element']
    
    # If we don't have both top and bottom elements, we can't create a proper negation map
    if not (top_bottom['has_top'] and top_bottom['has_bottom']):
        return {}, {'error': 'Lattice must have both top and bottom elements for negation'}
    
    nodes = list(G.nodes())
    negation_map = {}
    score_info = {
        'complement_score': 0,
        'involution_score': 0,
        'contrapositive_score': 0,
        'demorgan_score': 0,
        'node_details': {}
    }
    
    # For each node, find the best negation candidate
    for node in nodes:
        # Skip top and bottom elements (they map to each other)
        if node == top_element:
            negation_map[node] = bottom_element
            continue
        if node == bottom_element:
            negation_map[node] = top_element
            continue
        
        # Find negation candidates
        candidates = propose_negation_candidates(RN, G, node, top_element, bottom_element)
        
        # Select the best candidate
        if candidates['ranked_candidates']:
            best_candidate, best_score, score_breakdown = candidates['ranked_candidates'][0]
            negation_map[node] = best_candidate
            score_info['node_details'][node] = {
                'negation': best_candidate,
                'score': best_score,
                'breakdown': score_breakdown
            }
            
            # Update overall scores
            score_info['complement_score'] += score_breakdown.get('complement', 0)
            score_info['involution_score'] += score_breakdown.get('double_negation', 0)
            # Other scores will be computed after the full map is created
    
    # Compute involution score - check if double negation returns to original
    successful_involutions = 0
    for node in nodes:
        if node in negation_map:
            neg = negation_map[node]
            if neg in negation_map and negation_map[neg] == node:
                successful_involutions += 1
                if node in score_info['node_details']:
                    score_info['node_details'][node]['breakdown']['double_negation'] = 1
    
    score_info['involution_score'] = successful_involutions / len(nodes) if nodes else 0
    
    # Check De Morgan's laws for random pairs
    demorgan_tests = 0
    successful_demorgan = 0
    
    # Limit to 20 random tests to avoid excessive computation
    import random
    test_pairs = []
    if len(nodes) >= 2:
        for _ in range(min(20, len(nodes))):
            node1, node2 = random.sample(nodes, 2)
            test_pairs.append((node1, node2))
        
        for node1, node2 in test_pairs:
            demorgan_result = test_demorgan_laws(RN, G, node1, node2, negation_map)
            demorgan_tests += 1
            if demorgan_result['law1_success']:
                successful_demorgan += 0.5
            if demorgan_result['law2_success']:
                successful_demorgan += 0.5
    
    score_info['demorgan_score'] = successful_demorgan / demorgan_tests if demorgan_tests > 0 else 0
    
    # Compute overall negation quality score
    overall_score = (
        score_info['complement_score'] / len(nodes) +
        score_info['involution_score'] +
        score_info['demorgan_score']
    ) / 3 if nodes else 0
    
    score_info['overall_score'] = overall_score
    
    return negation_map, score_info

def print_negation_map_report(negation_map, score_info, G):
    """
    Prints a report on the quality of the negation map.
    
    Parameters:
    -----------
    negation_map : dict
        Dictionary mapping each node to its negation.
    score_info : dict
        Information about the quality of the negation map.
    G : networkx.DiGraph
        The directed graph representing the hierarchy of organizations.
    """
    print("\n===== NEGATION MAP ANALYSIS =====")
    
    if 'error' in score_info:
        print(f"Error: {score_info['error']}")
        return
    
    print(f"Overall negation quality score: {score_info['overall_score']:.2f}")
    print(f"Complement score: {score_info['complement_score']:.2f}")
    print(f"Involution score: {score_info['involution_score']:.2f}")
    print(f"De Morgan score: {score_info['demorgan_score']:.2f}")
    
    print("\nNegation mappings:")
    for node, neg_node in negation_map.items():
        node_species = get_species_from_node(G, node)
        neg_species = get_species_from_node(G, neg_node)
        print(f"¬{node} = {neg_node}")
        print(f"  {node} species: {sorted(list(node_species))}")
        print(f"  {neg_node} species: {sorted(list(neg_species))}")
        
        # Check if double negation works
        double_neg = negation_map.get(neg_node, None)
        if double_neg == node:
            print(f"  ✓ Double negation works: ¬¬{node} = {node}")
        else:
            print(f"  ✗ Double negation fails: ¬¬{node} = {double_neg} ≠ {node}")

def analyze_negation_properties(RN, G):
    """
    Comprehensive analysis of negation properties for the organizational structure.
    
    Parameters:
    -----------
    RN : ReactionNetwork
        The reaction network object.
    G : networkx.DiGraph
        The directed graph representing the hierarchy of organizations.
        
    Returns:
    --------
    dict
        A dictionary containing all analysis results.
    """
    print("\n===== NEGATION PROPERTIES ANALYSIS =====")
    
    # Check for top and bottom elements first
    top_bottom = identify_top_bottom_elements(G)
    print_top_bottom_report(top_bottom, G)
    
    # If we don't have both top and bottom, we can't define proper negation
    if not (top_bottom['has_top'] and top_bottom['has_bottom']):
        print("\nCannot perform full negation analysis without both top and bottom elements.")
        print("Will analyze partial negation properties instead.")
    
    # Generate best possible negation map
    negation_map, score_info = generate_best_negation_map(RN, G)
    
    # Print report on the negation map
    print_negation_map_report(negation_map, score_info, G)
    
    return {
        'top_bottom': top_bottom,
        'negation_map': negation_map,
        'score_info': score_info
    }