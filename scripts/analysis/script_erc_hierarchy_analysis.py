#!/usr/bin/env python3
"""
Enhanced ERC Hierarchy Analysis Script

This script analyzes ERC hierarchies as partial ordered sets (posets),
extracting structural properties relevant for understanding their 
organization and computational efficiency.

Usage:
    python enhanced_hierarchy_benchmark.py <network_file>
    python enhanced_hierarchy_benchmark.py <network_folder>
"""

import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
from collections import defaultdict, Counter, deque
import seaborn as sns
from pathlib import Path
import os
import sys
import glob

# Import required modules
from pyCOT.io.functions import read_txt
from pyCOT.analysis.ERC_Hierarchy import ERC, species_list_to_names

def count_cliques_in_hierarchy(graph):
    """Count the number of connected components (cliques) in the hierarchy graph"""
    if len(graph.nodes()) == 0:
        return 0
        
    if graph.is_directed():
        # For directed graphs, use weakly connected components
        return len(list(nx.weakly_connected_components(graph)))
    else:
        # For undirected graphs, use connected components
        return len(list(nx.connected_components(graph)))

def analyze_chains_in_hierarchy(graph):
    """Analyze chain structure in the hierarchy DAG"""
    if len(graph.nodes()) == 0:
        return []
        
    if not graph.is_directed():
        # Convert to directed if needed
        dag = nx.DiGraph()
        dag.add_edges_from(graph.edges())
    else:
        dag = graph
    
    # Find all maximal chains (paths)
    chains = []
    
    try:
        # Find all simple paths in the DAG
        # A chain is a maximal path where each internal node has exactly one predecessor and one successor
        in_degrees = dict(dag.in_degree())
        out_degrees = dict(dag.out_degree())
        
        # Find all potential chain start nodes (in_degree <= 1)
        potential_starts = [node for node in dag.nodes() if in_degrees[node] == 0]
        
        visited_in_chains = set()
        
        for start_node in potential_starts:
            if start_node in visited_in_chains:
                continue
                
            # Build chain starting from this node
            chain = [start_node]
            current = start_node
            
            # Follow the unique path as long as current node has exactly one successor
            # and that successor has exactly one predecessor
            while (out_degrees[current] == 1):
                successors = list(dag.successors(current))
                next_node = successors[0]
                
                if in_degrees[next_node] == 1:  # Next node has exactly one predecessor
                    chain.append(next_node)
                    current = next_node
                else:
                    break
            
            # Only consider it a chain if it has more than 1 node
            if len(chain) > 1:
                chains.append(chain)
                visited_in_chains.update(chain)
        
    except Exception as e:
        print(f"Warning: Could not analyze chains: {e}")
        return []
    
    return chains

def analyze_cliques_depth_and_branching(graph):
    """Analyze depth and branching factor of each clique (connected component)"""
    if len(graph.nodes()) == 0:
        return []
        
    # Get connected components (cliques)
    if graph.is_directed():
        cliques = list(nx.weakly_connected_components(graph))
    else:
        cliques = list(nx.connected_components(graph))
    
    clique_stats = []
    
    for clique_nodes in cliques:
        if len(clique_nodes) == 1:
            # Single node clique
            clique_stats.append({
                'nodes': list(clique_nodes),
                'depth': 1,
                'avg_branching': 0
            })
            continue
            
        # Extract subgraph for this clique
        subgraph = graph.subgraph(clique_nodes)
        
        # Convert to directed for depth analysis
        if not subgraph.is_directed():
            dag_subgraph = nx.DiGraph()
            dag_subgraph.add_edges_from(subgraph.edges())
        else:
            dag_subgraph = subgraph
        
        # Calculate depth (longest path in this clique)
        try:
            if nx.is_directed_acyclic_graph(dag_subgraph):
                depth = nx.dag_longest_path_length(dag_subgraph) + 1  # +1 to convert path length to depth
            else:
                # If it has cycles, use number of nodes as approximation
                depth = len(clique_nodes)
        except:
            depth = len(clique_nodes)
        
        # Calculate average branching factor within this clique
        if dag_subgraph.number_of_nodes() > 0:
            out_degrees = [deg for node, deg in dag_subgraph.out_degree()]
            avg_branching = np.mean(out_degrees) if out_degrees else 0
        else:
            avg_branching = 0
        
        clique_stats.append({
            'nodes': list(clique_nodes),
            'depth': depth,
            'avg_branching': avg_branching
        })
    
    return clique_stats

def calculate_chain_statistics(chains, total_ercs):
    """Calculate statistics about chains"""
    if not chains or total_ercs == 0:
        return {}, {}
        
    chain_lengths = [len(chain) for chain in chains]
    chain_length_counts = Counter(chain_lengths)
    
    # Calculate proportion of ERCs in chains of each length
    chain_proportions = {}
    for length, count in chain_length_counts.items():
        ercs_in_chains_of_this_length = count * length
        proportion = ercs_in_chains_of_this_length / total_ercs if total_ercs > 0 else 0
        chain_proportions[length] = proportion
    
    return chain_length_counts, chain_proportions

def analyze_hierarchy_structure(graph, ercs, rn):
    """Analyze the structural properties of an ERC hierarchy graph"""
    
    if len(graph.nodes()) == 0:
        return {}
    
    # Basic graph properties
    n_nodes = len(graph.nodes())
    n_edges = len(graph.edges())
    
    # Convert to directed graph for proper analysis
    if not graph.is_directed():
        # Assume it's a DAG representing hierarchy relationships
        dag = nx.DiGraph()
        dag.add_edges_from(graph.edges())
    else:
        dag = graph
    
    # Count cliques (connected components)
    n_cliques = count_cliques_in_hierarchy(graph)
    
    # Analyze chains
    chains = analyze_chains_in_hierarchy(graph)
    chain_length_counts, chain_proportions = calculate_chain_statistics(chains, n_nodes)
    
    # Analyze cliques (connected components) depth and branching
    clique_stats = analyze_cliques_depth_and_branching(graph)
    
    # Hierarchy levels using topological sorting
    try:
        levels = list(nx.topological_generations(dag))
        n_levels = max(1, len(levels))  # Ensure minimum of 1 level
        level_sizes = [len(level) for level in levels]
        avg_level_size = np.mean(level_sizes) if level_sizes else 1  # Changed from 0 to 1
        max_level_size = max(level_sizes) if level_sizes else 1  # Changed from 0 to 1
        level_variance = np.var(level_sizes) if len(level_sizes) > 1 else 0
    except nx.NetworkXError:
        # Graph has cycles or other issues
        n_levels = 1  # Changed from 0 to 1
        level_sizes = []
        avg_level_size = 1  # Changed from 0 to 1
        max_level_size = 1  # Changed from 0 to 1
        level_variance = 0
    
    # Connectivity properties
    if n_edges > 0 and n_nodes > 1:
        max_possible_edges = n_nodes * (n_nodes - 1) / 2  # For undirected
        hierarchy_density = n_edges / max_possible_edges
        avg_degree = 2 * n_edges / n_nodes  # For undirected graph
    else:
        hierarchy_density = 0
        avg_degree = 0
    
    # Leaf and root nodes
    if dag.number_of_nodes() > 0:
        in_degrees = dict(dag.in_degree())
        out_degrees = dict(dag.out_degree())
        
        root_nodes = sum(1 for node, degree in in_degrees.items() if degree == 0)
        leaf_nodes = sum(1 for node, degree in out_degrees.items() if degree == 0)
        
        # Ratios
        root_ratio = root_nodes / n_nodes if n_nodes > 0 else 0
        leaf_ratio = leaf_nodes / n_nodes if n_nodes > 0 else 0
    else:
        root_nodes = 0
        leaf_nodes = 0
        root_ratio = 0
        leaf_ratio = 0
    
    # Chain analysis (longest paths)
    try:
        if nx.is_directed_acyclic_graph(dag):
            # Find longest path
            longest_path_length = nx.dag_longest_path_length(dag)
            max_chain_length = longest_path_length
        else:
            max_chain_length = 0
    except:
        max_chain_length = 0
    
    # Branching factor analysis
    if dag.number_of_nodes() > 0:
        out_degrees_values = list(out_degrees.values())
        avg_branching = np.mean(out_degrees_values) if out_degrees_values else 0
        max_branching = max(out_degrees_values) if out_degrees_values else 0
    else:
        avg_branching = 0
        max_branching = 0
    
    # ERC-specific properties
    avg_generators = 0  # Placeholder - would need specific ERC analysis
    max_generators = 0  # Placeholder
    
    # Calculate average closure size with error handling
    try:
        closure_sizes = []
        for erc in ercs:
            if hasattr(erc, 'get_closure') and rn:
                closure_sizes.append(len(erc.get_closure(rn)))
        avg_closure_size = np.mean(closure_sizes) if closure_sizes else 0
    except Exception as e:
        print(f"Warning: Could not calculate closure sizes: {e}")
        avg_closure_size = 0
    
    ercs_per_species = len(ercs) / len(rn.species()) if ercs and rn and len(rn.species()) > 0 else 0
    
    # Hierarchy shape metrics
    if n_levels > 0:
        hierarchy_depth_ratio = max_chain_length / n_levels if n_levels > 0 else 0
    else:
        hierarchy_depth_ratio = 0
    
    return {
        # Basic properties
        'n_edges': n_edges,
        'hierarchy_density': hierarchy_density,
        'avg_degree': avg_degree,
        'n_cliques': n_cliques,
        
        # Level structure
        'n_levels': n_levels,
        'avg_level_size': avg_level_size,
        'max_level_size': max_level_size,
        'level_variance': level_variance,
        
        # Node types
        'leaf_nodes': leaf_nodes,
        'root_nodes': root_nodes,
        'leaf_ratio': leaf_ratio,
        'root_ratio': root_ratio,
        
        # Chain properties
        'max_chain_length': max_chain_length,
        'hierarchy_depth_ratio': hierarchy_depth_ratio,
        'chains': chains,
        'chain_length_counts': chain_length_counts,
        'chain_proportions': chain_proportions,
        'clique_stats': clique_stats,
        
        # Branching properties
        'avg_branching': avg_branching,
        'max_branching': max_branching,
        
        # ERC properties
        'avg_generators': avg_generators,
        'max_generators': max_generators,
        'avg_closure_size': avg_closure_size,
        'ercs_per_species': ercs_per_species
    }

def benchmark_single_network(network_file):
    """Benchmark hierarchy building on a single network file"""
    print("="*80)
    print(f"BENCHMARKING: {os.path.basename(network_file)}")
    print("="*80)
    
    try:
        # Load the reaction network
        print("Loading reaction network...")
        rn = read_txt(network_file)
        print(f"Network loaded: {len(rn.species())} species, {len(rn.reactions())} reactions")
        
        # Generate ERCs
        print("Generating ERCs...")
        start_erc_time = time.time()
        ercs = ERC.ERCs(rn)
        erc_time = time.time() - start_erc_time
        print(f"Generated {len(ercs)} ERCs in {erc_time:.4f}s")
        
        if len(ercs) == 0:
            print("⚠️ No ERCs found, skipping hierarchy tests")
            return None
            
    except Exception as e:
        print(f"❌ Error loading network: {e}")
        return None
    
    # Initialize success flags and results
    success= False
    graph= None
    
    # Test methods based on input
    print("\n" + "-"*50)
    print("TESTING ORIGINAL METHOD...")
    print("-"*50)
    start_time = time.time()
    try:
        graph = ERC.build_hierarchy_graph(ercs, rn)
        hierarchy_time = time.time() - start_time
        print(f"✓ Hierarchy calculation: {hierarchy_time:.4f}s, {len(graph.edges())} edges")
        success= True
    except Exception as e:
        print(f"✗ ERC hierarchy calculation failed: {e}")
        hierarchy_time = float('inf')
        graph = None
        success = False

    successful_graphs = []
    if success and graph is not None:
        successful_graphs.append(("ERC_hierarchy", set(graph.edges())))
        results_match = True
        print(f"✓ {successful_graphs[0][0]} succeeded")
    else:
        print("✗ The hierarchy computation did not succeed")

    # Choose best available graph for structural analysis

    # Analyze hierarchy structure
    hierarchy_props = {}
    if graph:
        print(f"\n" + "-"*50)
        print(f"ANALYZING HIERARCHY STRUCTURE)")
        print("-"*50)
        hierarchy_props = analyze_hierarchy_structure(graph, ercs, rn)
        
        print(f"Hierarchy Analysis:")
        print(f"  Levels: {hierarchy_props.get('n_levels', 0)}")
        print(f"  Max chain length: {hierarchy_props.get('max_chain_length', 0)}")
        print(f"  Average branching: {hierarchy_props.get('avg_branching', 0):.2f}")
        print(f"  Leaf ratio: {hierarchy_props.get('leaf_ratio', 0):.3f}")
        print(f"  Root ratio: {hierarchy_props.get('root_ratio', 0):.3f}")
        print(f"  Hierarchy density: {hierarchy_props.get('hierarchy_density', 0):.4f}")
        print(f"  Number of cliques: {hierarchy_props.get('n_cliques', 0)}")
        print(f"  Number of chains: {len(hierarchy_props.get('chains', []))}")
        print(f"  Number of clique stats: {len(hierarchy_props.get('clique_stats', []))}")

    # Performance summary
    print("\n" + "="*80)
    print("PERFORMANCE SUMMARY")
    print("="*80)
    
    methods_data = []
    if success:
        methods_data.append(("ERC hierarchy", hierarchy_time))
    
    if methods_data:
        fastest_time = methods_data[0][1]
        
        print(f"{'Hierarchy':<20} {'Time (s)':<12} {'Speedup':<10}")
        print("-" * 45)
        
        for method_name, exec_time in methods_data:
            speedup = fastest_time / exec_time if exec_time > 0 else 1.0
            speedup_str = f"{speedup:.1f}x" if speedup != 1.0 else "baseline"
            print(f"{method_name:<20} {exec_time:<12.4f} {speedup_str:<10}")
        
        print(f"\n🏆 Computation: {methods_data[0][0]} ({methods_data[0][1]:.4f}s)")
            
    # Build comprehensive result dictionary
    result = {
        'file': os.path.basename(network_file),
        'n_species': len(rn.species()),
        'n_reactions': len(rn.reactions()),
        'n_ercs': len(ercs),
        'erc_generation_time': erc_time,
        'results_match': results_match,
        #Structure
        'ERC_hierarchy': graph,
        'RN':rn,
        # Algorithm timing and success
        'time': hierarchy_time if success else None,
        'edges': len(graph.edges()) if success and graph else None,
        'success': success,
        }
    
    # Add hierarchy structural properties
    result.update(hierarchy_props)
    
    
    return result

def benchmark_multiple_networks(network_files, output_dir="benchmark_results", method="all", 
                                output_file="erc_analysis.csv"):
    """Benchmark hierarchy building methods on multiple network files"""
    print("="*80)
    print(f"MULTI-NETWORK ERC HIERARCHY BENCHMARK")
    print("="*80)
    print(f"Testing {len(network_files)} networks with method(s): {method}")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    results = []
    
    for i, network_file in enumerate(network_files, 1):
        print(f"\n{'='*20} NETWORK {i}/{len(network_files)} {'='*20}")
        result = benchmark_single_network(network_file, method)
        if result:
            results.append(result)
        
        # Small delay between tests
        time.sleep(0.1)
    
    if not results:
        print("❌ No successful benchmarks to analyze")
        return None
    
    # Create summary DataFrame
    df = pd.DataFrame(results)
    
    # Save detailed results
    csv_file = os.path.join(output_dir, output_file)
    df.to_csv(csv_file, index=False)
    print(f"\n💾 Detailed results saved to: {csv_file}")
    
    # Create visualizations
    create_hierarchy_visualizations(df, output_dir)
    create_custom_analysis_plots(df, output_dir)
    
    # Print summary statistics
    print_hierarchy_statistics(df)
    
    return df

def create_custom_analysis_plots(df, output_dir):
    """Create the three specific custom analysis plots requested"""
    
    plt.style.use('seaborn-v0_8' if 'seaborn-v0_8' in plt.style.available else 'default')
    
    fig, axes = plt.subplots(1, 3, figsize=(21, 6))
    fig.suptitle('Custom ERC Hierarchy Analysis', fontsize=16, fontweight='bold')
    
    # Plot 1: Number of reactions vs ratios
    ax1 = axes[0]
    
    # Calculate ratios
    erc_reaction_ratio = df['n_ercs'] / df['n_reactions']
    clique_reaction_ratio = df['n_cliques'] / df['n_reactions'] if 'n_cliques' in df.columns else None
    
    # Plot ERCs/Reactions ratio
    ax1.scatter(df['n_reactions'], erc_reaction_ratio, alpha=0.6, label='ERCs/Reactions', color='blue', s=50)
    
    # Plot Cliques/Reactions ratio
    if clique_reaction_ratio is not None:
        ax1.scatter(df['n_reactions'], clique_reaction_ratio, alpha=0.6, label='Cliques/Reactions', color='green', s=50)
    
    ax1.set_xlabel('Number of Reactions')
    ax1.set_ylabel('Ratio')
    ax1.set_title('Network Component Ratios vs Reactions')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_xscale('log')
    
    # Plot 2: Clique size analysis with proportional coloring
    ax2 = axes[1]
    
    # Prepare data for clique size vs ERCs plot
    plot_data = []
    
    for idx, row in df.iterrows():
        n_ercs = row['n_ercs']
        clique_stats = row.get('clique_stats', [])
        
        if isinstance(clique_stats, list) and clique_stats:
            # Calculate clique size distribution
            clique_sizes = []
            for clique_info in clique_stats:
                if isinstance(clique_info, dict) and 'nodes' in clique_info:
                    clique_size = len(clique_info['nodes'])
                    clique_sizes.append(clique_size)
            
            if clique_sizes:
                # Calculate proportion of ERCs in cliques of each size
                from collections import Counter
                size_counts = Counter(clique_sizes)
                
                for clique_size, count in size_counts.items():
                    ercs_in_cliques_of_this_size = count * clique_size
                    proportion = ercs_in_cliques_of_this_size / n_ercs if n_ercs > 0 else 0
                    
                    plot_data.append({
                        'n_ercs': n_ercs,
                        'clique_size': clique_size,
                        'proportion': proportion
                    })
    
    if plot_data:
        plot_df = pd.DataFrame(plot_data)
        
        # Create scatter plot with color based on proportion (red=0, blue=1)
        scatter = ax2.scatter(plot_df['n_ercs'], plot_df['clique_size'], 
                             c=plot_df['proportion'], cmap='coolwarm', 
                             alpha=0.7, s=60, vmin=0, vmax=1)
        
        # Add colorbar
        cbar = plt.colorbar(scatter, ax=ax2)
        cbar.set_label('Proportion of ERCs in cliques of this size')
        
        ax2.set_xlabel('Number of ERCs in Network')
        ax2.set_ylabel('Clique Size')
        ax2.set_title('Clique Size Distribution (Color = Proportion)')
        ax2.grid(True, alpha=0.3)
        ax2.set_xscale('log')
    else:
        ax2.text(0.5, 0.5, 'No clique data available', 
                ha='center', va='center', transform=ax2.transAxes)
        ax2.set_xlabel('Number of ERCs in Network')
        ax2.set_ylabel('Clique Size')
        ax2.set_title('Clique Size Distribution (No data)')
    
    # Plot 3: Clique depth vs average branching factor within cliques
    ax3 = axes[2]
    
    # Prepare data for clique depth vs branching
    clique_data = []
    
    for idx, row in df.iterrows():
        clique_stats = row.get('clique_stats', [])
        
        if isinstance(clique_stats, list) and clique_stats:
            for clique_info in clique_stats:
                if isinstance(clique_info, dict):
                    depth = clique_info.get('depth', 0)
                    avg_branching = clique_info.get('avg_branching', 0)
                    
                    if depth > 0:  # Only include valid cliques
                        clique_data.append({
                            'depth': depth,
                            'avg_branching': avg_branching
                        })
    
    if clique_data:
        clique_df = pd.DataFrame(clique_data)
        
        # Group by depth and calculate statistics
        grouped = clique_df.groupby('depth')['avg_branching'].agg(['mean', 'std', 'count']).reset_index()
        
        # Plot with error bars
        ax3.errorbar(grouped['depth'], grouped['mean'], 
                    yerr=grouped['std'], fmt='o-', capsize=5, alpha=0.7, label='Average ± Std')
        
        # Also show individual points
        ax3.scatter(clique_df['depth'], clique_df['avg_branching'], 
                   alpha=0.3, s=20, color='red', label='Individual cliques')
        
        ax3.set_xlabel('Clique Depth (Longest chain within clique)')
        ax3.set_ylabel('Average Branching Factor within Clique')
        ax3.set_title('Clique Depth vs Internal Branching Factor')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
    else:
        ax3.text(0.5, 0.5, 'No clique data available', 
                ha='center', va='center', transform=ax3.transAxes)
        ax3.set_xlabel('Clique Depth')
        ax3.set_ylabel('Average Branching Factor')
        ax3.set_title('Clique Depth vs Internal Branching Factor (No data)')
    
    plt.tight_layout()
    plot_file = os.path.join(output_dir, "custom_erc_analysis.png")
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"📊 Custom analysis plots saved to: {plot_file}")
    plt.show()

def create_hierarchy_visualizations(df, output_dir):
    """Create hierarchy structure visualization plots"""
    plt.style.use('seaborn-v0_8' if 'seaborn-v0_8' in plt.style.available else 'default')
    
    fig, axes = plt.subplots(3, 3, figsize=(18, 15))
    fig.suptitle('ERC Hierarchy Structure Analysis', fontsize=16, fontweight='bold')
    
    # 1. Hierarchy depth distribution
    if 'n_levels' in df.columns:
        axes[0, 0].hist(df['n_levels'].dropna(), bins=20, alpha=0.7, edgecolor='black')
        axes[0, 0].set_xlabel('Number of Hierarchy Levels')
        axes[0, 0].set_ylabel('Frequency')
        axes[0, 0].set_title('Hierarchy Depth Distribution')
    
    # 2. Chain length distribution  
    if 'max_chain_length' in df.columns:
        axes[0, 1].hist(df['max_chain_length'].dropna(), bins=20, alpha=0.7, edgecolor='black')
        axes[0, 1].set_xlabel('Maximum Chain Length')
        axes[0, 1].set_ylabel('Frequency')
        axes[0, 1].set_title('Chain Length Distribution')
    
    # 3. Branching factor distribution
    if 'avg_branching' in df.columns:
        axes[0, 2].hist(df['avg_branching'].dropna(), bins=20, alpha=0.7, edgecolor='black')
        axes[0, 2].set_xlabel('Average Branching Factor')
        axes[0, 2].set_ylabel('Frequency')
        axes[0, 2].set_title('Branching Distribution')
    
    # 4. Network size vs hierarchy depth
    if 'n_ercs' in df.columns and 'n_levels' in df.columns:
        axes[1, 0].scatter(df['n_ercs'], df['n_levels'], alpha=0.6)
        axes[1, 0].set_xlabel('Number of ERCs')
        axes[1, 0].set_ylabel('Hierarchy Levels')
        axes[1, 0].set_title('Network Size vs Hierarchy Depth')
        axes[1, 0].set_xscale('log')
    
    # 5. Branching vs depth relationship
    if 'avg_branching' in df.columns and 'n_levels' in df.columns:
        axes[1, 1].scatter(df['avg_branching'], df['n_levels'], alpha=0.6)
        axes[1, 1].set_xlabel('Average Branching Factor')
        axes[1, 1].set_ylabel('Hierarchy Levels')
        axes[1, 1].set_title('Branching vs Depth')
    
    # 6. Hierarchy density distribution
    if 'hierarchy_density' in df.columns:
        axes[1, 2].hist(df['hierarchy_density'].dropna(), bins=20, alpha=0.7, edgecolor='black')
        axes[1, 2].set_xlabel('Hierarchy Density')
        axes[1, 2].set_ylabel('Frequency') 
        axes[1, 2].set_title('Connectivity Density Distribution')
    
    # 7. Leaf vs root ratio
    if 'leaf_ratio' in df.columns and 'root_ratio' in df.columns:
        axes[2, 0].scatter(df['leaf_ratio'], df['root_ratio'], alpha=0.6)
        axes[2, 0].set_xlabel('Leaf Node Ratio')
        axes[2, 0].set_ylabel('Root Node Ratio')
        axes[2, 0].set_title('Structural Balance (Leaf vs Root)')
    
    # 8. Performance vs structure
    if 'ultra_optimized_time' in df.columns and 'n_levels' in df.columns:
        valid_data = df[df['ultra_optimized_time'].notna()]
        if len(valid_data) > 0:
            axes[2, 1].scatter(valid_data['n_levels'], valid_data['ultra_optimized_time'], alpha=0.6)
            axes[2, 1].set_xlabel('Hierarchy Levels')
            axes[2, 1].set_ylabel('Computation Time (s)')
            axes[2, 1].set_title('Complexity vs Performance')
            axes[2, 1].set_yscale('log')
    
    # 9. ERC size vs closure size
    if 'n_ercs' in df.columns and 'avg_closure_size' in df.columns:
        axes[2, 2].scatter(df['n_ercs'], df['avg_closure_size'], alpha=0.6)
        axes[2, 2].set_xlabel('Number of ERCs')
        axes[2, 2].set_ylabel('Average Closure Size')
        axes[2, 2].set_title('ERC Count vs Closure Size')
        axes[2, 2].set_xscale('log')
    
    plt.tight_layout()
    plot_file = os.path.join(output_dir, "erc_hierarchy_analysis.png")
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"📊 Hierarchy analysis plots saved to: {plot_file}")
    plt.show()

def print_hierarchy_statistics(df):
    """Print comprehensive hierarchy statistics"""
    print("\n" + "="*80)
    print("ERC HIERARCHY STRUCTURE STATISTICS")
    print("="*80)
    
    print(f"Dataset Summary:")
    print(f"  Total networks analyzed: {len(df)}")
    print(f"  Network size range: {df['n_ercs'].min()}-{df['n_ercs'].max()} ERCs")
    print(f"  Species range: {df['n_species'].min()}-{df['n_species'].max()}")
    print(f"  Reactions range: {df['n_reactions'].min()}-{df['n_reactions'].max()}")
    
    # Hierarchy structure statistics
    if 'n_levels' in df.columns:
        levels_data = df['n_levels'].dropna()
        if len(levels_data) > 0:
            print(f"\nHierarchy Depth Statistics:")
            print(f"  Average levels: {levels_data.mean():.2f}")
            print(f"  Median levels: {levels_data.median():.1f}")
            print(f"  Max levels: {levels_data.max()}")
            print(f"  Min levels: {levels_data.min()}")
    
    # Chain analysis
    if 'max_chain_length' in df.columns:
        chain_data = df['max_chain_length'].dropna()
        if len(chain_data) > 0:
            print(f"\nChain Length Statistics:")
            print(f"  Average max chain: {chain_data.mean():.2f}")
            print(f"  Median max chain: {chain_data.median():.1f}")
            print(f"  Longest chain: {chain_data.max()}")
    
    # Branching analysis
    if 'avg_branching' in df.columns:
        branching_data = df['avg_branching'].dropna()
        if len(branching_data) > 0:
            print(f"\nBranching Statistics:")
            print(f"  Average branching factor: {branching_data.mean():.3f}")
            print(f"  Median branching factor: {branching_data.median():.3f}")
            print(f"  Max branching factor: {branching_data.max():.3f}")
    
    # Clique analysis
    if 'n_cliques' in df.columns:
        clique_data = df['n_cliques'].dropna()
        if len(clique_data) > 0:
            print(f"\nClique Statistics:")
            print(f"  Average cliques: {clique_data.mean():.2f}")
            print(f"  Median cliques: {clique_data.median():.1f}")
            print(f"  Max cliques: {clique_data.max()}")
    
    # Clique depth and branching analysis
    all_clique_depths = []
    all_clique_branchings = []
    for idx, row in df.iterrows():
        clique_stats = row.get('clique_stats', [])
        if isinstance(clique_stats, list):
            for clique_info in clique_stats:
                if isinstance(clique_info, dict):
                    depth = clique_info.get('depth', 0)
                    branching = clique_info.get('avg_branching', 0)
                    if depth > 0:
                        all_clique_depths.append(depth)
                        all_clique_branchings.append(branching)
    
    if all_clique_depths:
        print(f"\nClique Depth Statistics:")
        print(f"  Average clique depth: {np.mean(all_clique_depths):.2f}")
        print(f"  Median clique depth: {np.median(all_clique_depths):.1f}")
        print(f"  Max clique depth: {max(all_clique_depths)}")
        
    if all_clique_branchings:
        print(f"\nClique Internal Branching Statistics:")
        print(f"  Average branching within cliques: {np.mean(all_clique_branchings):.3f}")
        print(f"  Median branching within cliques: {np.median(all_clique_branchings):.3f}")
        print(f"  Max branching within cliques: {max(all_clique_branchings):.3f}")
    
    # Structural balance
    if 'leaf_ratio' in df.columns and 'root_ratio' in df.columns:
        leaf_data = df['leaf_ratio'].dropna()
        root_data = df['root_ratio'].dropna()
        if len(leaf_data) > 0 and len(root_data) > 0:
            print(f"\nStructural Balance:")
            print(f"  Average leaf ratio: {leaf_data.mean():.3f}")
            print(f"  Average root ratio: {root_data.mean():.3f}")
            print(f"  Leaf/Root balance: {leaf_data.mean()/root_data.mean():.2f}")
    
    # Connectivity
    if 'hierarchy_density' in df.columns:
        density_data = df['hierarchy_density'].dropna()
        if len(density_data) > 0:
            print(f"\nConnectivity Statistics:")
            print(f"  Average hierarchy density: {density_data.mean():.4f}")
            print(f"  Median hierarchy density: {density_data.median():.4f}")
    
    # Performance statistics (if available)
    perf_methods = [
        ('ultra_optimized_time', 'Ultra-optimized'),
        ('optimized_time', 'Optimized'),
        ('original_time', 'Original')
    ]
    
    print(f"\nPerformance Statistics:")
    for col, name in perf_methods:
        if col in df.columns:
            data = df[col].dropna()
            if len(data) > 0:
                print(f"  {name}: {data.mean():.4f}s avg, {data.median():.4f}s median")

def find_network_files(path):
    """Find all network files in a given path"""
    if os.path.isfile(path):
        return [path]
    elif os.path.isdir(path):
        patterns = ['*.txt', '*.net', '*.rn']
        files = []
        for pattern in patterns:
            files.extend(glob.glob(os.path.join(path, pattern)))
        return sorted(files)
    else:
        print(f"❌ Path not found: {path}")
        return []

def main():
    """Main function to run the enhanced benchmark"""
    # Configuration - modify these paths as needed
    input_path = "networks/RandomAlife"  # Your networks directory
    output_file ="erc_analysis_random.csv"
    input_path = "networks/RandomAlife/RN_Ns_40_Norg_14_id_356.txt"
    # Find network files
    network_files = find_network_files(input_path)
    
    
    if not network_files:
        print(f"❌ No network files found in: {input_path}")
        return
    
    print(f"Found {len(network_files)} network file(s)")
    
    if len(network_files) == 1:
        # Single network benchmark
        result = benchmark_single_network(network_files[0])
        # Visualize network
        hierarchy=result["ERC_hierarchy"]
        rn=result["RN"]
        ERC.plot_hierarchy(RN=rn, graph=hierarchy)

        if result:
            print(f"\n✅ Benchmark completed successfully")
        else:
            print(f"\n❌ Benchmark failed")
    else:
        # Multiple network benchmark
        df = benchmark_multiple_networks(network_files)
        if df is not None:
            print(f"\n✅ Multi-network benchmark completed successfully")
            print(f"📁 Results saved in 'benchmark_results/' directory")
        else:
            print(f"\n❌ Multi-network benchmark failed")

if __name__ == "__main__":
    main()