# Script: Organization-Based Abstraction Analysis (IMPROVED VERSION)
# Key improvements:
# 1. Stability-based simulation stopping (distribution convergence)
# 2. Analysis of both abstractions AND their closures
# 3. Proper subnetwork simulation (restricts dynamics to organization)

# ========================================
# 1. LIBRARY LOADING AND CONFIGURATION
# ========================================
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
from collections import defaultdict, Counter
from tqdm import tqdm
from scipy.stats import entropy

# Add the project root directory to the PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import pyCOT modules
from pyCOT.io.functions import read_txt
from pyCOT.visualization.rn_visualize import rn_visualize_html
from pyCOT.analysis.Persistent_Modules_Generator import compute_all_organizations
from pyCOT.simulations import simulation
from pyCOT.simulations.abstractions import abstraction_ordinary
from pyCOT.analysis.ERC_Hierarchy import ERC_Hierarchy

# ========================================
# 2. CONFIGURATION PARAMETERS
# ========================================
# Network file path
#file_path = 'networks/testing/Farm.txt'  # Specify the network file
#file_path = 'networks/Marine_Ecosystem/Las_Cruces_251021.txt'  # Specify the network file
#file_path = 'networks/testing/autopoietic.txt'  # Specify the network file
file_path = 'networks/Conflict_Theory/Resource_Community_Insurgency_Loops_model1.txt'  # Conflict model


# Simulation parameters
N_sim = 100          # Number of simulations per organization
t_max = 200          # Maximum simulation time
t_transient_tmax = 180  # If t_max reached, use last t_transient_tmax for frequencies
n_steps_per_unit = 5   # Time steps per time unit

# Stability detection parameters
N_buffer = 40        # Number of recent abstractions to keep in buffer
N_abs_threshold = 2  # Maximum number of unique abstractions for stability
stability_checks = 10 # Number of consecutive checks needed to confirm stability
check_interval = 20  # Check stability every N time steps

# Parameter ranges for random conditions
sp_min = 0.5         # Minimum value for species parameters
sp_max = 1.5         # Maximum value for species parameters
k_min = 0.5          # Minimum value for rate constants
k_max = 1.5          # Maximum value for rate constants

# Abstraction threshold
threshold = 0.0001     # Threshold for abstraction_ordinary

# ========================================
# 3. HELPER FUNCTIONS
# ========================================

def compute_closure(rn, species_set):
    """
    Compute the closure of a species set using the reaction network.
    Returns the closed set of species names.
    """
    closure_species = rn.generated_closure(list(species_set))
    return set(sp.name for sp in closure_species)

def check_stability(abstraction_buffer, n_abs_threshold):
    """
    Check if the distribution of abstractions has stabilized.
    Returns True if the buffer contains N_abs_threshold or fewer unique abstractions.
    """
    if len(abstraction_buffer) < N_buffer:
        return False
    
    unique_abstractions = len(set(tuple(sorted(abs)) for abs in abstraction_buffer))
    return unique_abstractions <= n_abs_threshold

def compute_distribution(abstraction_list):
    """
    Compute the frequency distribution of abstractions.
    Returns a Counter object.
    """
    return Counter(tuple(sorted(abs)) for abs in abstraction_list)

# ========================================
# 4. LOAD AND VISUALIZE NETWORK
# ========================================
print("="*80)
print("LOADING REACTION NETWORK")
print("="*80)
rn = read_txt(file_path)
print(f"Network loaded from: {file_path}")
print(f"Number of species: {len(rn.species())}")
print(f"Number of reactions: {len(rn.reactions())}")

# Visualize the network
print("\nGenerating network visualization...")
rn_visualize_html(rn, filename="network_visualization.html")
print("Network visualization saved as: network_visualization.html")

# ========================================
# 5. CALCULATE AND PRINT ERC HIERARCHY
# ========================================
print("\n" + "="*80)
print("CALCULATING ERC HIERARCHY")
print("="*80)

try:
    hierarchy = ERC_Hierarchy(rn, verbose=True)
    print("\nERC Hierarchy Structure:")
    print(f"Number of ERCs: {len(hierarchy.ercs)}")
    print(f"Number of trees in forest: {len(hierarchy.trees)}")
    
    if hierarchy.hierarchy.graph:
        print(f"Hierarchy graph has {hierarchy.hierarchy.graph.number_of_nodes()} nodes")
        print(f"Hierarchy graph has {hierarchy.hierarchy.graph.number_of_edges()} edges")
    else:
        print("No hierarchical structure (single level)")
        
except Exception as e:
    print(f"Warning: Could not compute ERC hierarchy: {e}")
    hierarchy = None

# ========================================
# 6. COMPUTE ORGANIZATIONS
# ========================================
print("\n" + "="*80)
print("COMPUTING ORGANIZATIONS")
print("="*80)

results = compute_all_organizations(
    rn, 
    max_generator_size=8,
    max_organization_size=10,
    verbose=True
)

# Extract organizations
all_organizations = results['elementary_organizations']
all_organizations_sets = []
for org in all_organizations:
    if hasattr(org, "combined_closure"):
        sset = set(sp.name for sp in org.combined_closure)
    elif hasattr(org, "closure_species"):
        sset = set(sp.name for sp in org.closure_species)
    else:
        sset = set()
    all_organizations_sets.append(sset)

orgs_sets = sorted(all_organizations_sets, key=len)
print(f"\nFound {len(orgs_sets)} organizations:")
for i, org in enumerate(orgs_sets):
    print(f"  Organization {i}: {org} (size: {len(org)})")

# ========================================
# 7. SIMULATION WITH STABILITY DETECTION
# ========================================

def run_simulation_with_stability(rn_sim, org_species, n_sim, t_max, t_transient_tmax, 
                                   n_steps_per_unit, sp_min, sp_max, k_min, k_max, 
                                   threshold, N_buffer, N_abs_threshold, stability_checks, 
                                   check_interval):
    """
    Run simulations with stability-based stopping criterion.
    Simulates on the SUBNETWORK defined by org_species, not the full network.
    
    Returns:
        abstraction_data: List of dicts with abstraction frequencies and metadata
        closure_data: List of dicts with closure frequencies and metadata
    """
    abstraction_data = []
    closure_data = []
    
    print(f"\n  Running {n_sim} simulations with stability detection...")
    
    for sim_idx in tqdm(range(n_sim), desc="  Simulations"):
        try:
            # Create subnetwork for this organization
            rn_sub = rn_sim.sub_reaction_network(list(org_species))
            
            num_species_sub = len(rn_sub.species())
            num_reactions_sub = len(rn_sub.reactions())
            
            if num_reactions_sub == 0:
                print(f"\n  Warning: Organization has no reactions, skipping simulation {sim_idx}")
                continue
            
            # Generate random initial conditions for subnetwork species
            x0 = np.random.uniform(sp_min, sp_max, num_species_sub)
            
            # Generate random rate constants for subnetwork reactions
            spec_vector = [[np.random.uniform(k_min, k_max)] for _ in range(num_reactions_sub)]
            
            # Initialize tracking
            abstraction_buffer = []
            stability_counter = 0
            current_time = 0
            all_abstractions = []
            all_closures = []
            
            # Simulation loop with adaptive stopping
            reached_stability = False
            
            while current_time < t_max:
                # Determine time span for this chunk
                chunk_duration = check_interval
                t_end = min(current_time + chunk_duration, t_max)
                n_steps = int((t_end - current_time) * n_steps_per_unit)
                
                # Run simulation chunk
                time_series, _ = simulation(
                    rn_sub,
                    rate='mak',
                    spec_vector=spec_vector,
                    x0=x0,
                    t_span=(current_time, t_end),
                    n_steps=max(n_steps, 10)  # At least 10 steps
                )
                
                # Update x0 for next chunk (continuity)
                x0 = time_series.iloc[-1, 1:].values
                
                # Compute abstractions for this chunk
                abstract_ts = abstraction_ordinary(time_series, threshold=threshold)
                
                # Track abstractions and closures
                for abstraction in abstract_ts['Abstraction']:
                    all_abstractions.append(abstraction)
                    abstraction_buffer.append(abstraction)
                    
                    # Compute closure for this abstraction
                    closure_set = compute_closure(rn_sim, abstraction)
                    all_closures.append(closure_set)
                
                # Keep buffer size manageable
                if len(abstraction_buffer) > N_buffer:
                    abstraction_buffer = abstraction_buffer[-N_buffer:]
                
                # Check stability
                if check_stability(abstraction_buffer, N_abs_threshold):
                    stability_counter += 1
                    if stability_counter >= stability_checks:
                        reached_stability = True
                        break
                else:
                    stability_counter = 0
                
                current_time = t_end
            
            # Determine which data to use for frequency calculation
            if reached_stability:
                # Use all data (it's all stationary after stabilization)
                abstractions_to_count = all_abstractions
                closures_to_count = all_closures
                status = 'stabilized'
            else:
                # Use last t_transient_tmax portion
                time_points = len(all_abstractions)
                cutoff_idx = max(0, time_points - int(t_transient_tmax * n_steps_per_unit * check_interval / 10))
                abstractions_to_count = all_abstractions[cutoff_idx:]
                closures_to_count = all_closures[cutoff_idx:]
                status = 'max_time_reached'
            
            # Compute abstraction frequencies
            abs_freq = compute_distribution(abstractions_to_count)
            
            # Compute closure frequencies
            closure_freq = compute_distribution(closures_to_count)
            
            # Store results
            abstraction_data.append({
                'sim_idx': sim_idx,
                'frequencies': abs_freq,
                'total_time': current_time,
                'status': status,
                'num_unique': len(abs_freq)
            })
            
            closure_data.append({
                'sim_idx': sim_idx,
                'frequencies': closure_freq,
                'total_time': current_time,
                'status': status,
                'num_unique': len(closure_freq)
            })
            
        except Exception as e:
            print(f"\n  Warning: Simulation {sim_idx} failed: {e}")
            continue
    
    return abstraction_data, closure_data

# ========================================
# 8. ANALYZE EACH ORGANIZATION
# ========================================
print("\n" + "="*80)
print("ANALYZING ORGANIZATIONS")
print("="*80)

# Store results for all organizations
all_org_results = []

for org_idx, org_set in enumerate(orgs_sets):
    print(f"\n{'='*80}")
    print(f"ORGANIZATION {org_idx}: {org_set}")
    print(f"Size: {len(org_set)} species")
    print(f"{'='*80}")
    
    # Check if organization has closure equal to itself
    org_closure = compute_closure(rn, org_set)
    print(f"Closure size: {len(org_closure)} species")
    print(f"Is closed: {org_closure == org_set}")
    
    # Run simulations with stability detection
    abs_data, closure_data = run_simulation_with_stability(
        rn, org_set, N_sim, t_max, t_transient_tmax, n_steps_per_unit,
        sp_min, sp_max, k_min, k_max, threshold,
        N_buffer, N_abs_threshold, stability_checks, check_interval
    )
    
    # Aggregate abstraction frequencies across all simulations
    aggregated_abs_freq = Counter()
    for sim_data in abs_data:
        aggregated_abs_freq.update(sim_data['frequencies'])
    
    # Aggregate closure frequencies across all simulations
    aggregated_closure_freq = Counter()
    for sim_data in closure_data:
        aggregated_closure_freq.update(sim_data['frequencies'])
    
    # Compute transition frequencies for abstractions
    abs_transitions = Counter()
    for sim_data in abs_data:
        freq_items = list(sim_data['frequencies'].items())
        for i in range(len(freq_items)):
            for j in range(len(freq_items)):
                if i != j:
                    source, source_freq = freq_items[i]
                    target, target_freq = freq_items[j]
                    # Estimate transitions (proportional to product of frequencies)
                    abs_transitions[(source, target)] += source_freq * target_freq
    
    # Compute transition frequencies for closures
    closure_transitions = Counter()
    for sim_data in closure_data:
        freq_items = list(sim_data['frequencies'].items())
        for i in range(len(freq_items)):
            for j in range(len(freq_items)):
                if i != j:
                    source, source_freq = freq_items[i]
                    target, target_freq = freq_items[j]
                    closure_transitions[(source, target)] += source_freq * target_freq
    
    # Store results
    all_org_results.append({
        'org_idx': org_idx,
        'org_set': org_set,
        'org_closure': org_closure,
        'abstraction_frequencies': aggregated_abs_freq,
        'closure_frequencies': aggregated_closure_freq,
        'abstraction_transitions': abs_transitions,
        'closure_transitions': closure_transitions,
        'simulation_data': {
            'abstractions': abs_data,
            'closures': closure_data
        }
    })
    
    # Print statistics
    print(f"\n  Simulation Statistics:")
    stabilized = sum(1 for d in abs_data if d['status'] == 'stabilized')
    print(f"    Stabilized: {stabilized}/{len(abs_data)}")
    print(f"    Reached max time: {len(abs_data) - stabilized}/{len(abs_data)}")
    
    print(f"\n  Abstraction Analysis:")
    print(f"    Total unique abstractions: {len(aggregated_abs_freq)}")
    print(f"    Top 5 abstractions:")
    for abstraction, freq in aggregated_abs_freq.most_common(5):
        is_org = set(abstraction) in orgs_sets
        org_marker = " [ORG]" if is_org else ""
        print(f"      {set(abstraction)}: {freq}{org_marker}")
    
    print(f"\n  Closure Analysis:")
    print(f"    Total unique closures: {len(aggregated_closure_freq)}")
    print(f"    Top 5 closures:")
    for closure_tuple, freq in aggregated_closure_freq.most_common(5):
        closure_set = set(closure_tuple)
        is_org = closure_set in orgs_sets
        org_marker = " [ORG]" if is_org else ""
        print(f"      {closure_set}: {freq}{org_marker}")

# ========================================
# 9. VISUALIZATION FUNCTIONS
# ========================================

def create_frequency_histogram(frequencies, orgs_sets, title, filename, data_type="Abstraction"):
    """Create histogram with organization highlighting"""
    if len(frequencies) == 0:
        print(f"  Warning: No {data_type.lower()}s to plot")
        return
    
    # Prepare data
    items = []
    counts = []
    colors = []
    
    for item_tuple, count in frequencies.most_common(30):
        items.append(str(set(item_tuple)))
        counts.append(count)
        is_org = set(item_tuple) in orgs_sets
        colors.append('green' if is_org else 'skyblue')
    
    # Create plot
    fig, ax = plt.subplots(figsize=(16, 6))
    bars = ax.bar(range(len(counts)), counts, color=colors)
    
    ax.set_xlabel(data_type, fontsize=12)
    ax.set_ylabel('Frequency', fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.set_xticks(range(len(items)))
    ax.set_xticklabels(items, rotation=45, ha='right', fontsize=7)
    ax.grid(axis='y', alpha=0.3)
    
    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='green', label='Organization'),
        Patch(facecolor='skyblue', label='Non-organization')
    ]
    ax.legend(handles=legend_elements, loc='upper right')
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

def create_transition_graph(transitions, frequencies, orgs_sets, title, filename):
    """Create transition graph visualization"""
    G = nx.DiGraph()
    
    # Add organization nodes
    org_nodes = []
    for i, org in enumerate(orgs_sets):
        node_label = f"Org_{i}"
        G.add_node(node_label, label=str(org), is_org=True, org_set=org)
        org_nodes.append((node_label, org))
    
    # Add non-organization node
    G.add_node("Non-Org", label="Non-Organization", is_org=False)
    
    # Map items to nodes
    def find_node(item_set):
        for node_label, org in org_nodes:
            if item_set == org:
                return node_label
        return "Non-Org"
    
    # Process transitions
    transition_counts = defaultdict(int)
    for (source_tuple, target_tuple), weight in transitions.items():
        source_set = set(source_tuple)
        target_set = set(target_tuple)
        
        source_node = find_node(source_set)
        target_node = find_node(target_set)
        
        if source_node != target_node:
            transition_counts[(source_node, target_node)] += weight
    
    # Add edges
    for (source, target), weight in transition_counts.items():
        G.add_edge(source, target, weight=weight)
    
    if G.number_of_edges() == 0:
        print(f"  No transitions to plot")
        return
    
    # Visualization
    fig, ax = plt.subplots(figsize=(14, 10))
    pos = nx.spring_layout(G, k=3, iterations=50, seed=42)
    
    # Draw nodes
    org_node_list = [n for n in G.nodes() if n != "Non-Org"]
    non_org_node_list = ["Non-Org"] if "Non-Org" in G.nodes() else []
    
    nx.draw_networkx_nodes(G, pos, nodelist=org_node_list,
                          node_color='lightgreen', node_size=2500,
                          alpha=0.8, ax=ax)
    nx.draw_networkx_nodes(G, pos, nodelist=non_org_node_list,
                          node_color='lightcoral', node_size=2500,
                          alpha=0.8, ax=ax)
    
    # Draw edges
    edges = G.edges()
    weights = [G[u][v]['weight'] for u, v in edges]
    max_weight = max(weights) if weights else 1
    
    nx.draw_networkx_edges(G, pos, edgelist=edges,
                          width=[5 * w / max_weight for w in weights],
                          alpha=0.6, edge_color='gray',
                          arrows=True, arrowsize=25, ax=ax,
                          connectionstyle='arc3,rad=0.1')
    
    # Labels
    nx.draw_networkx_labels(G, pos, font_size=9, font_weight='bold', ax=ax)
    
    # Edge labels (show only significant ones)
    significant_edges = [(u, v) for u, v in edges if G[u][v]['weight'] > max_weight * 0.1]
    edge_labels = {(u, v): f"{int(G[u][v]['weight'])}" for u, v in significant_edges}
    nx.draw_networkx_edge_labels(G, pos, edge_labels, font_size=8, ax=ax)
    
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.axis('off')
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

# ========================================
# 10. GENERATE ALL VISUALIZATIONS
# ========================================
print("\n" + "="*80)
print("GENERATING VISUALIZATIONS")
print("="*80)

os.makedirs("visualizations/organization_analysis", exist_ok=True)

for result in all_org_results:
    org_idx = result['org_idx']
    org_set = result['org_set']
    
    print(f"\nGenerating plots for Organization {org_idx}...")
    
    # Abstraction frequency histogram
    create_frequency_histogram(
        result['abstraction_frequencies'],
        orgs_sets,
        f'Abstraction Frequencies - Organization {org_idx}\n{org_set}',
        f"visualizations/organization_analysis/org_{org_idx}_abstraction_freq.png",
        data_type="Abstraction"
    )
    
    # Closure frequency histogram
    create_frequency_histogram(
        result['closure_frequencies'],
        orgs_sets,
        f'Closure Frequencies - Organization {org_idx}\n{org_set}',
        f"visualizations/organization_analysis/org_{org_idx}_closure_freq.png",
        data_type="Closure"
    )
    
    # Abstraction transition graph
    create_transition_graph(
        result['abstraction_transitions'],
        result['abstraction_frequencies'],
        orgs_sets,
        f'Abstraction Transitions - Organization {org_idx}',
        f"visualizations/organization_analysis/org_{org_idx}_abstraction_transitions.png"
    )
    
    # Closure transition graph
    create_transition_graph(
        result['closure_transitions'],
        result['closure_frequencies'],
        orgs_sets,
        f'Closure Transitions - Organization {org_idx}',
        f"visualizations/organization_analysis/org_{org_idx}_closure_transitions.png"
    )

# ========================================
# 11. FINAL SUMMARY
# ========================================
print("\n" + "="*80)
print("FINAL SUMMARY")
print("="*80)

for result in all_org_results:
    org_idx = result['org_idx']
    print(f"\nOrganization {org_idx}:")
    print(f"  Unique abstractions: {len(result['abstraction_frequencies'])}")
    print(f"  Unique closures: {len(result['closure_frequencies'])}")
    print(f"  Abstraction transitions: {len(result['abstraction_transitions'])}")
    print(f"  Closure transitions: {len(result['closure_transitions'])}")

print("\n" + "="*80)
print("ANALYSIS COMPLETE")
print("="*80)
print(f"\nAll visualizations saved in: visualizations/organization_analysis/")
print("\nGenerated files:")
print("  - Abstraction frequency histograms")
print("  - Closure frequency histograms")
print("  - Abstraction transition graphs")
print("  - Closure transition graphs")
