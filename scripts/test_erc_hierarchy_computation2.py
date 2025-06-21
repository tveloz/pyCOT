#!/usr/bin/env python3
"""
Hierarchy Graph Building Methods Performance Comparison Script

This script compares the performance of three different methods for building 
ERC hierarchy graphs: original, optimized, and ultra-optimized versions.

Usage:
    python hierarchy_benchmark.py <network_file>
    python hierarchy_benchmark.py <network_folder>
"""

import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
from collections import defaultdict, Counter
import seaborn as sns
from pathlib import Path
import os
import sys
import glob

# Import required modules
from pyCOT.io.functions import read_txt
from pyCOT.ERC_Hierarchy import ERC, species_list_to_names

def benchmark_single_network(network_file):
    """Benchmark hierarchy building methods on a single network file"""
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
    
    # Test original method
    print("\n" + "-"*50)
    print("TESTING ORIGINAL METHOD...skipped")
    print("-"*50)
    # start_time = time.time()
    # try:
    #     graph1 = ERC.build_hierarchy_graph(ercs, rn)
    #     time1 = time.time() - start_time
    #     print(f"✓ Original: {time1:.4f}s, {len(graph1.edges())} edges")
    #     success1 = True
    # except Exception as e:
    #     print(f"✗ Original failed: {e}")
    #     time1 = float('inf')
    #     graph1 = None
    #     success1 = False

    # Test optimized method
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

    # Test ultra-optimized method
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

    # Verify results match (only for successful methods)
    print("\n" + "-"*50)
    print("VERIFYING RESULTS")
    print("-"*50)
    
    successful_graphs = []
    # if success1 and graph1 is not None:
    #     successful_graphs.append(("Original", set(graph1.edges())))
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

    # Performance summary
    print("\n" + "="*80)
    print("PERFORMANCE SUMMARY")
    print("="*80)
    
    methods_data = []
    # if success1:
    #     methods_data.append(("Original", time1))
    if success2:
        methods_data.append(("Optimized", time2))
    if success3:
        methods_data.append(("Ultra-optimized", time3))
    
    if methods_data:
        # Sort by time
        methods_data.sort(key=lambda x: x[1])
        fastest_time = methods_data[0][1]
        
        print(f"{'Method':<20} {'Time (s)':<12} {'Speedup':<10} {'Status'}")
        print("-" * 55)
        
        for method, exec_time in methods_data:
            speedup = fastest_time / exec_time if exec_time > 0 else 1.0
            speedup_str = f"{speedup:.1f}x" if speedup != 1.0 else "baseline"
            print(f"{method:<20} {exec_time:<12.4f} {speedup_str:<10}")
        
        print(f"\n🏆 Fastest method: {methods_data[0][0]} ({methods_data[0][1]:.4f}s)")
        
        if len(methods_data) > 1:
            overall_speedup = methods_data[-1][1] / methods_data[0][1]
            print(f"📈 Maximum speedup: {overall_speedup:.1f}x improvement")
    
    # Return benchmark results
    return {
        'network_file': os.path.basename(network_file),
        'n_species': len(rn.species()),
        'n_reactions': len(rn.reactions()),
        'n_ercs': len(ercs),
        'erc_generation_time': erc_time,
        # 'original_time': time1 if success1 else None,
        'optimized_time': time2 if success2 else None,
        'ultra_optimized_time': time3 if success3 else None,
        # 'original_edges': len(graph1.edges()) if success1 and graph1 else None,
        'optimized_edges': len(graph2.edges()) if success2 and graph2 else None,
        'ultra_optimized_edges': len(graph3.edges()) if success3 and graph3 else None,
        'results_match': results_match,
        # 'success1': success1,
        'success2': success2,
        'success3': success3
    }

def benchmark_multiple_networks(network_files, output_dir="benchmark_results"):
    """Benchmark hierarchy building methods on multiple network files"""
    print("="*80)
    print(f"MULTI-NETWORK HIERARCHY BENCHMARK")
    print("="*80)
    print(f"Testing {len(network_files)} networks")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    results = []
    
    for i, network_file in enumerate(network_files, 1):
        print(f"\n{'='*20} NETWORK {i}/{len(network_files)} {'='*20}")
        result = benchmark_single_network(network_file)
        if result:
            results.append(result)
        
        # Small delay between tests
        time.sleep(0.1)
    
    if not results:
        print("❌ No successful benchmarks to analyze")
        return
    
    # Create summary DataFrame
    df = pd.DataFrame(results)
    
    # Save detailed results
    csv_file = os.path.join(output_dir, "hierarchy_benchmark_results.csv")
    df.to_csv(csv_file, index=False)
    print(f"\n💾 Detailed results saved to: {csv_file}")
    
    # Create visualizations
    create_benchmark_visualizations(df, output_dir)
    
    # Print summary statistics
    print_summary_statistics(df)
    
    return df

def create_benchmark_visualizations(df, output_dir):
    """Create performance visualization plots"""
    plt.style.use('seaborn-v0_8' if 'seaborn-v0_8' in plt.style.available else 'default')
    
    # Filter successful results for plotting
    df_clean = df.copy()
    
    # 1. Performance comparison plot
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('Hierarchy Building Methods Performance Comparison', fontsize=16, fontweight='bold')
    
    # Time comparison
    ax1 = axes[0, 0]
    #methods = ['original_time', 'optimized_time', 'ultra_optimized_time']
    methods = [ 'optimized_time', 'ultra_optimized_time']
    #method_labels = ['Original', 'Optimized', 'Ultra-optimized']
    method_labels = [ 'Optimized', 'Ultra-optimized']
    
    valid_data = []
    valid_labels = []
    for method, label in zip(methods, method_labels):
        if method in df_clean.columns:
            data = df_clean[method].dropna()
            if len(data) > 0:
                valid_data.append(data)
                valid_labels.append(label)
    
    if valid_data:
        ax1.boxplot(valid_data, labels=valid_labels)
        ax1.set_ylabel('Time (seconds)')
        ax1.set_title('Execution Time Distribution')
        ax1.set_yscale('log')
    
    # Speedup vs network size
    ax2 = axes[0, 1]
    if 'original_time' in df_clean.columns and 'ultra_optimized_time' in df_clean.columns:
        mask = (df_clean['original_time'].notna()) & (df_clean['ultra_optimized_time'].notna())
        if mask.any():
            speedup = df_clean.loc[mask, 'original_time'] / df_clean.loc[mask, 'ultra_optimized_time']
            ercs = df_clean.loc[mask, 'n_ercs']
            ax2.scatter(ercs, speedup, alpha=0.7)
            ax2.set_xlabel('Number of ERCs')
            ax2.set_ylabel('Speedup Factor')
            ax2.set_title('Speedup vs Network Size')
            ax2.set_yscale('log')
    
    # Success rate
    ax3 = axes[1, 0]
    success_rates = []
    method_names = []
    if 'success1' in df_clean.columns:
        success_rates.append(df_clean['success1'].mean() * 100)
        method_names.append('Original')
    if 'success2' in df_clean.columns:
        success_rates.append(df_clean['success2'].mean() * 100)
        method_names.append('Optimized')
    if 'success3' in df_clean.columns:
        success_rates.append(df_clean['success3'].mean() * 100)
        method_names.append('Ultra-optimized')
    
    if success_rates:
        bars = ax3.bar(method_names, success_rates, color=['red', 'orange', 'green'][:len(success_rates)])
        ax3.set_ylabel('Success Rate (%)')
        ax3.set_title('Method Success Rate')
        ax3.set_ylim(0, 105)
        
        # Add value labels on bars
        for bar, rate in zip(bars, success_rates):
            height = bar.get_height()
            ax3.text(bar.get_x() + bar.get_width()/2., height + 1,
                    f'{rate:.1f}%', ha='center', va='bottom')
    
    # Network size distribution
    ax4 = axes[1, 1]
    ax4.hist(df_clean['n_ercs'], bins=min(20, len(df_clean)), alpha=0.7, edgecolor='black')
    ax4.set_xlabel('Number of ERCs')
    ax4.set_ylabel('Count')
    ax4.set_title('Network Size Distribution')
    
    plt.tight_layout()
    plot_file = os.path.join(output_dir, "hierarchy_benchmark_plots.png")
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"📊 Visualization plots saved to: {plot_file}")
    plt.show()

def print_summary_statistics(df):
    """Print summary statistics from benchmark results"""
    print("\n" + "="*80)
    print("SUMMARY STATISTICS")
    print("="*80)
    
    print(f"Total networks tested: {len(df)}")
    print(f"Network size range: {df['n_ercs'].min()}-{df['n_ercs'].max()} ERCs")
    
    # Success rates
    if 'success1' in df.columns:
        print(f"Original method success rate: {df['success1'].mean()*100:.1f}%")
    if 'success2' in df.columns:
        print(f"Optimized method success rate: {df['success2'].mean()*100:.1f}%")
    if 'success3' in df.columns:
        print(f"Ultra-optimized success rate: {df['success3'].mean()*100:.1f}%")
    
    # Performance statistics
    methods = [
        ('original_time', 'Original'),
        ('optimized_time', 'Optimized'), 
        ('ultra_optimized_time', 'Ultra-optimized')
    ]
    
    print(f"\nPerformance Statistics:")
    for col, name in methods:
        if col in df.columns:
            data = df[col].dropna()
            if len(data) > 0:
                print(f"  {name}: {data.mean():.4f}s avg, {data.median():.4f}s median, {data.std():.4f}s std")
    
    # Overall speedup
    if 'original_time' in df.columns and 'ultra_optimized_time' in df.columns:
        mask = (df['original_time'].notna()) & (df['ultra_optimized_time'].notna())
        if mask.any():
            speedups = df.loc[mask, 'original_time'] / df.loc[mask, 'ultra_optimized_time']
            print(f"\nSpeedup Statistics:")
            print(f"  Average speedup: {speedups.mean():.1f}x")
            print(f"  Median speedup: {speedups.median():.1f}x")
            print(f"  Maximum speedup: {speedups.max():.1f}x")

def find_network_files(path):
    """Find all network files in a given path"""
    if os.path.isfile(path):
        return [path]
    elif os.path.isdir(path):
        # Look for common network file extensions
        patterns = ['*.txt', '*.net', '*.rn']
        files = []
        for pattern in patterns:
            files.extend(glob.glob(os.path.join(path, pattern)))
        return sorted(files)
    else:
        print(f"❌ Path not found: {path}")
        return []

def main():
    """Main function to run the benchmark"""
    # if len(sys.argv) < 2:
    #     print("Usage:")
    #     print("  python hierarchy_benchmark.py <network_file>")
    #     print("  python hierarchy_benchmark.py <network_folder>")
    #     print("\nExample:")
    #     print("  python hierarchy_benchmark.py ERCs_test2.txt")
    #     print("  python hierarchy_benchmark.py networks/")
    #     return
    
    #input_path = sys.argv[1]
    input_path= 'networks/testing/ERCs_test2.txt'
    input_path= 'networks/Navarino/RN_IN_05.txt'
    input_path= 'networks/biomodels_interesting/bigg_iAF692.txt'  # Adjust path as needed
    network_files = find_network_files(input_path)
    
    if not network_files:
        print(f"❌ No network files found in: {input_path}")
        return
    
    print(f"Found {len(network_files)} network file(s)")
    
    if len(network_files) == 1:
        # Single network benchmark
        result = benchmark_single_network(network_files[0])
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