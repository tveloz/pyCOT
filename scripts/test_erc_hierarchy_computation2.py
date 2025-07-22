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
from pyCOT.ERC_Hierarchy import ERC, species_list_to_names

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
        
        # Branching properties
        'avg_branching': avg_branching,
        'max_branching': max_branching,
        
        # ERC properties
        'avg_generators': avg_generators,
        'max_generators': max_generators,
        'avg_closure_size': avg_closure_size,
        'ercs_per_species': ercs_per_species
    }

def benchmark_single_network(network_file, method="all"):
    """Benchmark hierarchy building methods on a single network file"""
    print("="*80)
    print(f"BENCHMARKING: {os.path.basename(network_file)}")
    print("Methods:", method)
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
    success1, success2, success3 = False, False, False
    graph1, graph2, graph3 = None, None, None
    time1, time2, time3 = float('inf'), float('inf'), float('inf')
    
    # Test methods based on input
    methods_to_test = method
    if "all" in method or "all"==method:
        methods_to_test.append("Original")
        methods_to_test.append("Optimized")  
        methods_to_test.append("Ultra-optimized")
    
    # Test Original method
    if "Original" in methods_to_test:
        print("\n" + "-"*50)
        print("TESTING ORIGINAL METHOD...")
        print("-"*50)
        start_time = time.time()
        try:
            graph1 = ERC.build_hierarchy_graph(ercs, rn)
            time1 = time.time() - start_time
            print(f"✓ Original: {time1:.4f}s, {len(graph1.edges())} edges")
            success1 = True
        except Exception as e:
            print(f"✗ Original failed: {e}")
            time1 = float('inf')
            graph1 = None
            success1 = False

    # Test Optimized method
    if "Optimized" in methods_to_test:
        print("\n" + "-"*50)
        print("TESTING OPTIMIZED METHOD")
        print("-"*50)
        start_time = time.time()
        try:
            graph2 = ERC.build_hierarchy_graph2(ercs, rn)
            time2 = time.time() - start_time
            print(f"✓ Optimized: {time2:.4f}s, {len(graph2.edges())} edges")
            success2 = True
        except Exception as e:
            print(f"✗ Optimized failed: {e}")
            time2 = float('inf')
            graph2 = None
            success2 = False

    # Test Ultra-optimized method  
    if "Ultra-optimized" in methods_to_test:
        print("\n" + "-"*50)
        print("TESTING ULTRA-OPTIMIZED METHOD")
        print("-"*50)
        start_time = time.time()
        try:
            graph3 = ERC.build_hierarchy_graph3(ercs, rn)
            time3 = time.time() - start_time
            print(f"✓ Ultra-optimized: {time3:.4f}s, {len(graph3.edges())} edges")
            success3 = True
        except Exception as e:
            print(f"✗ Ultra-optimized failed: {e}")
            time3 = float('inf')
            graph3 = None
            success3 = False

    # Verify results match
    print("\n" + "-"*50)
    print("VERIFYING RESULTS")
    print("-"*50)
    
    successful_graphs = []
    if success1 and graph1 is not None:
        successful_graphs.append(("Original", set(graph1.edges())))
    if success2 and graph2 is not None:
        successful_graphs.append(("Optimized", set(graph2.edges())))
    if success3 and graph3 is not None:
        successful_graphs.append(("Ultra-optimized", set(graph3.edges())))
    
    results_match = True
    if len(successful_graphs) > 1:
        reference_edges = successful_graphs[0][1]
        for name, edges in successful_graphs[1:]:
            if edges != reference_edges:
                results_match = False
                print(f"✗ WARNING: {name} produces different results!")
                print(f"  Reference: {len(reference_edges)} edges")
                print(f"  {name}: {len(edges)} edges")
                break
        
        if results_match:
            print("✓ All successful methods produce identical results")
    elif len(successful_graphs) == 1:
        print(f"✓ Only {successful_graphs[0][0]} method succeeded")
    else:
        print("✗ No methods succeeded")

    # Choose best available graph for structural analysis
    best_graph = None
    if success3 and graph3:
        best_graph = graph3
        best_method = "Ultra-optimized"
    elif success2 and graph2:
        best_graph = graph2  
        best_method = "Optimized"
    elif success1 and graph1:
        best_graph = graph1
        best_method = "Original"
    
    # Analyze hierarchy structure
    hierarchy_props = {}
    if best_graph:
        print(f"\n" + "-"*50)
        print(f"ANALYZING HIERARCHY STRUCTURE ({best_method})")
        print("-"*50)
        hierarchy_props = analyze_hierarchy_structure(best_graph, ercs, rn)
        
        print(f"Hierarchy Analysis:")
        print(f"  Levels: {hierarchy_props.get('n_levels', 0)}")
        print(f"  Max chain length: {hierarchy_props.get('max_chain_length', 0)}")
        print(f"  Average branching: {hierarchy_props.get('avg_branching', 0):.2f}")
        print(f"  Leaf ratio: {hierarchy_props.get('leaf_ratio', 0):.3f}")
        print(f"  Root ratio: {hierarchy_props.get('root_ratio', 0):.3f}")
        print(f"  Hierarchy density: {hierarchy_props.get('hierarchy_density', 0):.4f}")

    # Performance summary
    print("\n" + "="*80)
    print("PERFORMANCE SUMMARY")
    print("="*80)
    
    methods_data = []
    if success1:
        methods_data.append(("Original", time1))
    if success2:
        methods_data.append(("Optimized", time2))
    if success3:
        methods_data.append(("Ultra-optimized", time3))
    
    if methods_data:
        methods_data.sort(key=lambda x: x[1])
        fastest_time = methods_data[0][1]
        
        print(f"{'Method':<20} {'Time (s)':<12} {'Speedup':<10}")
        print("-" * 45)
        
        for method_name, exec_time in methods_data:
            speedup = fastest_time / exec_time if exec_time > 0 else 1.0
            speedup_str = f"{speedup:.1f}x" if speedup != 1.0 else "baseline"
            print(f"{method_name:<20} {exec_time:<12.4f} {speedup_str:<10}")
        
        print(f"\n🏆 Fastest method: {methods_data[0][0]} ({methods_data[0][1]:.4f}s)")
        
        if len(methods_data) > 1 and methods_data[0][1] > 0:
            overall_speedup = methods_data[-1][1] / methods_data[0][1]
        else:
            overall_speedup = 100
        print(f"📈 Maximum speedup: {overall_speedup:.1f}x improvement")
    
    # Build comprehensive result dictionary
    result = {
        'file': os.path.basename(network_file),
        'n_species': len(rn.species()),
        'n_reactions': len(rn.reactions()),
        'n_ercs': len(ercs),
        'erc_generation_time': erc_time,
        'results_match': results_match,
        
        # Algorithm timing and success
        'original_time': time1 if success1 else None,
        'optimized_time': time2 if success2 else None,
        'ultra_optimized_time': time3 if success3 else None,
        'original_edges': len(graph1.edges()) if success1 and graph1 else None,
        'optimized_edges': len(graph2.edges()) if success2 and graph2 else None,
        'ultra_optimized_edges': len(graph3.edges()) if success3 and graph3 else None,
        'success1': success1,
        'success2': success2, 
        'success3': success3
    }
    
    # Add hierarchy structural properties
    result.update(hierarchy_props)
    
    return result

def benchmark_multiple_networks(network_files, output_dir="benchmark_results", method="all"):
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
    csv_file = os.path.join(output_dir, "erc_hierarchy_benchmark.csv")
    df.to_csv(csv_file, index=False)
    print(f"\n💾 Detailed results saved to: {csv_file}")
    
    # Create visualizations
    create_hierarchy_visualizations(df, output_dir)
    
    # Print summary statistics
    print_hierarchy_statistics(df)
    
    return df

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
    input_path = "networks/biomodels_all_txt"  # Your networks directory
    methods_to_test = ["Original", "Ultra-optimized"]  # Focus on best performing method

    # Find network files
    network_files = find_network_files(input_path)
    
    if not network_files:
        print(f"❌ No network files found in: {input_path}")
        return
    
    print(f"Found {len(network_files)} network file(s)")
    
    if len(network_files) == 1:
        # Single network benchmark
        result = benchmark_single_network(network_files[0], method=methods_to_test)
        if result:
            print(f"\n✅ Benchmark completed successfully")
        else:
            print(f"\n❌ Benchmark failed")
    else:
        # Multiple network benchmark
        df = benchmark_multiple_networks(network_files, method=methods_to_test)
        if df is not None:
            print(f"\n✅ Multi-network benchmark completed successfully")
            print(f"📁 Results saved in 'benchmark_results/' directory")
        else:
            print(f"\n❌ Multi-network benchmark failed")

if __name__ == "__main__":
    main()