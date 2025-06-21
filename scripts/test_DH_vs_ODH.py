#!/usr/bin/env python3
"""
Comparison script for original vs optimized synergy algorithms.
Validates correctness and measures performance improvements.
Can handle single files or folders of files.
"""

import time
import os
import glob
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from pyCOT.io.functions import read_txt
from pyCOT.ERC_Hierarchy import ERC
from pyCOT.Dynamical_Hierarchy_Synergy_Algorithm import EfficientFundamentalSynergyFinder
from pyCOT.Opt_Dynamical_Hierarchy_Synergy import OptimizedFundamentalSynergyFinder


def compare_algorithms_on_file(rn_file: str, verbose: bool = True):
    """Compare original and optimized algorithms on a single network file."""
    if verbose:
        print(f"\n{'='*70}")
        print(f"Processing: {os.path.basename(rn_file)}")
        print(f"{'='*70}")
    
    try:
        # Load reaction network
        rn = read_txt(rn_file)
        
        # Compute ERCs and hierarchy once
        if verbose:
            print("Computing ERCs and hierarchy...")
        ercs = ERC.ERCs(rn)
        hierarchy = ERC.build_hierarchy_graph(ercs, rn)
        
        if verbose:
            print(f"Network has {len(ercs)} ERCs, {len(rn.species())} species, {len(rn.reactions())} reactions")
        
        # Run original algorithm
        # if verbose:
        #     print("\nRunning ORIGINAL algorithm...")
        
        # original_finder = EfficientFundamentalSynergyFinder(rn, ercs, hierarchy)
        # start_time = time.time()
        # original_results = original_finder.find_fundamental_synergies(verbose=False)
        # original_time = time.time() - start_time
        
        # if verbose:
        #     print(f"Original algorithm completed in {original_time:.3f} seconds")
        #     print(f"Found {len(original_results)} fundamental synergies")
        
        # Run optimized algorithm
        if verbose:
            print("\nRunning OPTIMIZED algorithm...")
        
        optimized_finder = OptimizedFundamentalSynergyFinder(rn, ercs, hierarchy)
        start_time = time.time()
        optimized_results = optimized_finder.find_fundamental_synergies(verbose=True)
        optimized_time = time.time() - start_time
        
        if verbose:
            print(f"Optimized algorithm completed in {optimized_time:.3f} seconds")
            print(f"Found {len(optimized_results)} fundamental synergies")
        
        # # Validate results
        # def synergy_to_tuple(syn):
        #     return tuple(sorted([syn['base1'].label, syn['base2'].label]) + [syn['target'].label])
        
        # original_set = set(synergy_to_tuple(s) for s in original_results)
        # optimized_set = set(synergy_to_tuple(s) for s in optimized_results)
        
        # results_match = original_set == optimized_set
        
        # if verbose:
        #     print("\nVALIDATION:")
        #     if results_match:
        #         print("✓ Results match! Both algorithms found the same fundamental synergies.")
        #     else:
        #         print("✗ Results differ!")
        #         if original_set - optimized_set:
        #             print(f"  Missing in optimized: {original_set - optimized_set}")
        #         if optimized_set - original_set:
        #             print(f"  Extra in optimized: {optimized_set - original_set}")
        
        # # Performance comparison
        # speedup = original_time / optimized_time if optimized_time > 0 else float('inf')
        
        # if verbose:
        #     print(f"\nPERFORMANCE:")
        #     print(f"  Original time: {original_time:.3f} seconds")
        #     print(f"  Optimized time: {optimized_time:.3f} seconds")
        #     print(f"  Speedup: {speedup:.2f}x")
        
        # return {
        #     'file': os.path.basename(rn_file),
        #     'path': rn_file,
        #     'n_ercs': len(ercs),
        #     'n_species': len(rn.species()),
        #     'n_reactions': len(rn.reactions()),
        #     'n_synergies': len(original_results),
        #     'original_time': original_time,
        #     'optimized_time': optimized_time,
        #     'speedup': speedup,
        #     'results_match': results_match,
        #     'original_stats': dict(original_finder.stats),
        #     'optimized_stats': dict(optimized_finder.stats),
        #     'error': None
        # }
        
    except Exception as e:
        print(f"Error processing {rn_file}: {str(e)}")
        return {
            'file': os.path.basename(rn_file),
            'path': rn_file,
            'error': str(e),
            'results_match': False
        }


def compare_algorithms(path: str, verbose: bool = True, file_pattern: str = "*.txt"):
    """
    Compare algorithms on a single file or all files in a folder.
    
    Parameters:
    - path: File path or folder path
    - verbose: Print detailed output
    - file_pattern: Pattern for files to process (if path is a folder)
    """
    results = []
    
    # Determine if path is file or folder
    if os.path.isfile(path):
        # Single file
        result = compare_algorithms_on_file(path, verbose=verbose)
        results.append(result)
        
    elif os.path.isdir(path):
        # Folder of files
        pattern = os.path.join(path, file_pattern)
        files = sorted(glob.glob(pattern))
        
        if not files:
            print(f"No files matching pattern '{file_pattern}' found in {path}")
            return results
            
        print(f"\nFound {len(files)} files to process in {path}")
        
        for i, file_path in enumerate(files, 1):
            print(f"\n[{i}/{len(files)}] ", end="")
            result = compare_algorithms_on_file(file_path, verbose=verbose)
            results.append(result)
            
    else:
        print(f"Path '{path}' is neither a file nor a directory")
        return results
    
    # Generate summary
    if len(results) > 1:
        generate_summary(results, verbose)
    
    return results


def generate_summary(results: list, verbose: bool = True):
    """Generate summary statistics and plots for multiple files."""
    # Filter out errors
    valid_results = [r for r in results if r.get('error') is None]
    error_results = [r for r in results if r.get('error') is not None]
    
    if not valid_results:
        print("\nNo valid results to summarize")
        return
    
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    
    print(f"\nProcessed {len(results)} files:")
    print(f"  Successful: {len(valid_results)}")
    print(f"  Errors: {len(error_results)}")
    
    if error_results:
        print("\nFiles with errors:")
        for r in error_results:
            print(f"  - {r['file']}: {r['error']}")
    
    # Check correctness
    matching = [r for r in valid_results if r['results_match']]
    print(f"\nCorrectness:")
    print(f"  Matching results: {len(matching)}/{len(valid_results)}")
    
    non_matching = [r for r in valid_results if not r['results_match']]
    if non_matching:
        print(f"  Non-matching files:")
        for r in non_matching:
            print(f"    - {r['file']}")
    
    # Performance statistics
    speedups = [r['speedup'] for r in valid_results]
    print(f"\nPerformance:")
    print(f"  Average speedup: {np.mean(speedups):.2f}x")
    print(f"  Median speedup: {np.median(speedups):.2f}x")
    print(f"  Min speedup: {np.min(speedups):.2f}x")
    print(f"  Max speedup: {np.max(speedups):.2f}x")
    
    # Network size statistics
    print(f"\nNetwork sizes:")
    print(f"  ERCs: {np.min([r['n_ercs'] for r in valid_results])}-{np.max([r['n_ercs'] for r in valid_results])}")
    print(f"  Synergies found: {np.min([r['n_synergies'] for r in valid_results])}-{np.max([r['n_synergies'] for r in valid_results])}")
    
    # Create summary plots if requested
    if verbose and len(valid_results) > 2:
        plot_comparison_results(valid_results)


def plot_comparison_results(results: list):
    """Create visualization of comparison results."""
    import numpy as np
    
    # Prepare data
    df = pd.DataFrame(results)
    df = df[df['error'].isna()].sort_values('n_ercs')
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
    
    # Plot 1: Execution times vs ERCs
    ax1.scatter(df['n_ercs'], df['original_time'], 
                label='Original', alpha=0.7, s=100, c='blue')
    ax1.scatter(df['n_ercs'], df['optimized_time'], 
                label='Optimized', alpha=0.7, s=100, c='red')
    ax1.set_xlabel('Number of ERCs')
    ax1.set_ylabel('Time (seconds)')
    ax1.set_title('Execution Time vs Network Size')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    if len(df) > 10:
        ax1.set_yscale('log')
    
    # Plot 2: Speedup vs ERCs
    ax2.scatter(df['n_ercs'], df['speedup'], alpha=0.7, s=100, c='green')
    ax2.axhline(y=1, color='black', linestyle='--', alpha=0.5)
    ax2.set_xlabel('Number of ERCs')
    ax2.set_ylabel('Speedup Factor')
    ax2.set_title('Speedup vs Network Size')
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Speedup distribution
    ax3.hist(df['speedup'], bins=min(20, len(df)//2), alpha=0.7, color='purple', edgecolor='black')
    ax3.axvline(x=df['speedup'].mean(), color='red', linestyle='--', 
                label=f'Mean: {df["speedup"].mean():.2f}x')
    ax3.axvline(x=df['speedup'].median(), color='green', linestyle='--',
                label=f'Median: {df["speedup"].median():.2f}x')
    ax3.set_xlabel('Speedup Factor')
    ax3.set_ylabel('Frequency')
    ax3.set_title('Distribution of Speedups')
    ax3.legend()
    ax3.grid(True, alpha=0.3, axis='y')
    
    # Plot 4: Network characteristics
    sizes = df['n_ercs'].values
    synergies = df['n_synergies'].values
    speedups = df['speedup'].values
    
    # Create scatter plot with color based on speedup
    scatter = ax4.scatter(sizes, synergies, c=speedups, s=100, 
                         cmap='viridis', alpha=0.7, edgecolors='black')
    ax4.set_xlabel('Number of ERCs')
    ax4.set_ylabel('Number of Synergies')
    ax4.set_title('Network Characteristics (color = speedup)')
    ax4.grid(True, alpha=0.3)
    
    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax4)
    cbar.set_label('Speedup Factor')
    
    # Add file labels for outliers
    for idx, row in df.iterrows():
        if row['speedup'] > df['speedup'].mean() + 2 * df['speedup'].std():
            ax4.annotate(row['file'], (row['n_ercs'], row['n_synergies']),
                        fontsize=8, alpha=0.7)
    
    plt.suptitle('Algorithm Comparison Results', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.show()


def save_results_to_csv(results: list, output_file: str = "comparison_results.csv"):
    """Save comparison results to CSV file."""
    df = pd.DataFrame(results)
    
    # Flatten the stats dictionaries
    for prefix in ['original', 'optimized']:
        stats_col = f'{prefix}_stats'
        if stats_col in df.columns:
            stats_df = pd.json_normalize(df[stats_col])
            stats_df.columns = [f'{prefix}_{col}' for col in stats_df.columns]
            df = pd.concat([df.drop(columns=[stats_col]), stats_df], axis=1)
    
    df.to_csv(output_file, index=False)
    print(f"\nResults saved to {output_file}")


if __name__ == "__main__":
    import argparse
    import numpy as np
    
    # parser = argparse.ArgumentParser(description='Compare synergy algorithms')
    # parser.add_argument('path', help='networks/testing/Farm.txt')
    # #parser.add_argument('--pattern', default='*.txt', help='File pattern for folders (default: *.txt)')
    # #parser.add_argument('--quiet', action='store_true', help='Reduce output verbosity')
    # #parser.add_argument('--save-csv', help='Save results to CSV file')
    
    # args = parser.parse_args()
    # print(args)
    # # Run comparison
    #results = compare_algorithms(args.path, verbose=not args.quiet, file_pattern=args.pattern)
    #results = compare_algorithms('networks/Navarino/RN_IN_05.txt', verbose=True)
    #results = compare_algorithms('networks/biomodels_interesting/central_ecoli.txt')
    results = compare_algorithms('networks/biomodels_interesting/bigg_iAF692.txt')
    
    
    # Save results if requested
    #if args.save_csv and results:
    #    save_results_to_csv(results, args.save_csv)