#!/usr/bin/env python3
"""
Multi-Network Synergy Algorithm Benchmarking Script

This script compares the brute force and dynamical hierarchy algorithms across
multiple reaction networks to understand performance scaling and identify factors
that contribute to efficiency improvements.

@author: Algorithm Benchmarking Analysis
"""

import os
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from collections import defaultdict
import networkx as nx
from typing import List, Dict, Tuple, Optional
import warnings
warnings.filterwarnings('ignore')

# Import required modules
from pyCOT.io.functions import read_txt
from pyCOT.ERC_Hierarchy import ERC, species_list_to_names
from pyCOT.Synergy_Brute_Force_Algorithm import BruteForceSynergyCalculator
from pyCOT.Dynamical_Hierarchy_Synergy_Algorithm import EfficientFundamentalSynergyFinder


class NetworkBenchmarker:
    """Benchmarks synergy algorithms across multiple reaction networks."""
    
    def __init__(self, network_folder: str, output_folder: str = "benchmark_results"):
        """
        Initialize benchmarker.
        
        Args:
            network_folder: Path to folder containing network files
            output_folder: Path to save results and figures
        """
        self.network_folder = network_folder
        self.output_folder = output_folder
        self.results = []
        self.network_data = {}
        
        # Create output folder if it doesn't exist
        os.makedirs(output_folder, exist_ok=True)
    
    def analyze_network_properties(self, rn, ercs, hierarchy_graph) -> Dict:
        """Extract detailed network properties for analysis."""
        # Basic network properties
        n_species = len(rn.species())
        n_reactions = len(rn.reactions())
        n_ercs = len(ercs)
        
        # Hierarchy properties
        levels = ERC.get_node_levels(hierarchy_graph)
        level_distribution = defaultdict(int)
        for erc, level in levels.items():
            level_distribution[level] += 1
        
        n_levels = len(level_distribution)
        avg_level_size = n_ercs / n_levels if n_levels > 0 else 0
        max_level_size = max(level_distribution.values()) if level_distribution else 0
        
        # Hierarchy structure metrics
        n_edges = hierarchy_graph.number_of_edges()
        avg_degree = 2 * n_edges / n_ercs if n_ercs > 0 else 0
        
        # Count leaf and root nodes
        leaf_nodes = sum(1 for node in hierarchy_graph.nodes() 
                        if hierarchy_graph.out_degree(node) == 0)
        root_nodes = sum(1 for node in hierarchy_graph.nodes() 
                        if hierarchy_graph.in_degree(node) == 0)
        
        # Maximum chain length
        max_chain_length = 0
        for node in hierarchy_graph.nodes():
            ancestors = nx.ancestors(hierarchy_graph, node)
            if len(ancestors) > max_chain_length:
                max_chain_length = len(ancestors)
        
        # Average generators per ERC
        avg_generators = np.mean([len(erc.min_generators) for erc in ercs])
        max_generators = max([len(erc.min_generators) for erc in ercs])
        
        # Closure size statistics
        closure_sizes = [len(species_list_to_names(erc.get_closure(rn))) for erc in ercs]
        avg_closure_size = np.mean(closure_sizes)
        
        # Hierarchy density (edges / possible edges)
        possible_edges = n_ercs * (n_ercs - 1) / 2
        hierarchy_density = n_edges / possible_edges if possible_edges > 0 else 0
        
        # Level variance (measure of hierarchy balance)
        level_counts = list(level_distribution.values())
        level_variance = np.var(level_counts) if level_counts else 0
        
        # Branching factor
        branching_factors = []
        for node in hierarchy_graph.nodes():
            out_degree = hierarchy_graph.out_degree(node)
            if out_degree > 0:
                branching_factors.append(out_degree)
        avg_branching = np.mean(branching_factors) if branching_factors else 0
        
        properties = {
            # Basic properties
            'n_species': n_species,
            'n_reactions': n_reactions,
            'n_ercs': n_ercs,
            
            # Hierarchy structure
            'n_levels': n_levels,
            'avg_level_size': avg_level_size,
            'max_level_size': max_level_size,
            'level_variance': level_variance,
            'n_edges': n_edges,
            'avg_degree': avg_degree,
            'hierarchy_density': hierarchy_density,
            
            # Node distribution
            'leaf_nodes': leaf_nodes,
            'root_nodes': root_nodes,
            'max_chain_length': max_chain_length,
            'avg_branching': avg_branching,
            
            # Generator and closure properties
            'avg_generators': avg_generators,
            'max_generators': max_generators,
            'avg_closure_size': avg_closure_size,
            
            # Derived metrics
            'ercs_per_species': n_ercs / n_species if n_species > 0 else 0,
            'hierarchy_depth_ratio': max_chain_length / n_ercs if n_ercs > 0 else 0,
            'leaf_ratio': leaf_nodes / n_ercs if n_ercs > 0 else 0,
            'root_ratio': root_nodes / n_ercs if n_ercs > 0 else 0,
        }
        
        return properties
    
    def benchmark_network(self, network_file: str, timeout: Optional[float] = 300) -> Optional[Dict]:
        """
        Benchmark both algorithms on a single network.
        
        Args:
            network_file: Path to network file
            timeout: Maximum time allowed per algorithm (seconds)
            
        Returns:
            Dictionary with benchmark results or None if error/timeout
        """
        print(f"\nBenchmarking: {os.path.basename(network_file)}")
        
        try:
            # Load network
            rn = read_txt(network_file)
            ercs = ERC.ERCs(rn)
            hierarchy_graph = ERC.build_hierarchy_graph(ercs, rn)
            
            # Get network properties
            properties = self.analyze_network_properties(rn, ercs, hierarchy_graph)
            print(f"  Network: {properties['n_species']} species, {properties['n_reactions']} reactions, {properties['n_ercs']} ERCs")
            print(f"  Hierarchy: {properties['n_levels']} levels, max chain {properties['max_chain_length']}")
            
            # Run brute force algorithm
            print("  Running brute force...", end="", flush=True)
            brute_calculator = BruteForceSynergyCalculator(rn, ercs, hierarchy_graph)
            
            brute_start = time.time()
            all_synergies, brute_fundamental = brute_calculator.brute_force(fundamental=True, verbose=False)
            brute_time = time.time() - brute_start
            
            if brute_time > timeout:
                print(f" TIMEOUT ({brute_time:.1f}s)")
                return None
            
            print(f" {brute_time:.3f}s")
            
            # Run efficient algorithm
            print("  Running efficient algorithm...", end="", flush=True)
            efficient_finder = EfficientFundamentalSynergyFinder(rn, ercs, hierarchy_graph)
            
            efficient_start = time.time()
            efficient_fundamental = efficient_finder.find_fundamental_synergies(verbose=False)
            efficient_time = time.time() - efficient_start
            
            print(f" {efficient_time:.3f}s")
            
            # Verify correctness
            brute_set = {(s['base1'].label, s['base2'].label, s['target'].label) 
                        for s in brute_fundamental}
            efficient_set = {(s['base1'].label, s['base2'].label, s['target'].label) 
                            for s in efficient_fundamental}
            
            results_match = (brute_set == efficient_set)
            if not results_match:
                print("  WARNING: Results don't match!")
            
            # Calculate performance metrics
            speedup = brute_time / efficient_time if efficient_time > 0 else float('inf')
            
            # Algorithm-specific metrics
            brute_checks = brute_calculator.stats['total_checks']
            efficient_ops = efficient_finder.stats['base_pairs_generated']
            
            # Theoretical maximum checks
            n = properties['n_ercs']
            max_possible_checks = n * (n - 1) * (n - 2) // 2  # All possible 3-ERC combinations
            
            # Efficiency metrics
            brute_efficiency = brute_checks / max_possible_checks if max_possible_checks > 0 else 0
            efficient_reduction = 1 - (efficient_ops / brute_checks) if brute_checks > 0 else 0
            
            # Steps per ERC metrics
            brute_steps_per_erc = brute_checks / n if n > 0 else 0
            efficient_steps_per_erc = efficient_ops / n if n > 0 else 0
            
            # Pruning effectiveness
            pruning_total = (efficient_finder.stats['within_target_pruning'] + 
                           efficient_finder.stats['cross_target_pruning'])
            pruning_per_synergy = (pruning_total / efficient_finder.stats['synergies_found'] 
                                 if efficient_finder.stats['synergies_found'] > 0 else 0)
            
            result = {
                'file': os.path.basename(network_file),
                'results_match': results_match,
                
                # Network properties
                **properties,
                
                # Performance metrics
                'brute_time': brute_time,
                'efficient_time': efficient_time,
                'speedup': speedup,
                
                # Algorithm metrics
                'brute_checks': brute_checks,
                'efficient_ops': efficient_ops,
                'max_possible_checks': max_possible_checks,
                'brute_efficiency': brute_efficiency,
                'efficient_reduction': efficient_reduction,
                
                # Per-ERC metrics
                'brute_steps_per_erc': brute_steps_per_erc,
                'efficient_steps_per_erc': efficient_steps_per_erc,
                
                # Synergy counts
                'total_synergies': len(all_synergies),
                'fundamental_synergies': len(brute_fundamental),
                'fundamentality_rate': len(brute_fundamental) / len(all_synergies) if all_synergies else 0,
                
                # Efficient algorithm specific
                'valid_base_target_pairs': efficient_finder.stats['valid_base_target_pairs'],
                'partial_overlap_checks': efficient_finder.stats['partial_overlap_checks'],
                'within_target_pruning': efficient_finder.stats['within_target_pruning'],
                'cross_target_pruning': efficient_finder.stats['cross_target_pruning'],
                'pruning_per_synergy': pruning_per_synergy,
            }
            
            # Store network data for later analysis
            self.network_data[network_file] = {
                'rn': rn,
                'ercs': ercs,
                'hierarchy_graph': hierarchy_graph
            }
            
            return result
            
        except Exception as e:
            print(f"  ERROR: {str(e)}")
            return None
    
    def run_benchmarks(self, file_pattern: str = "*.txt", max_files: Optional[int] = None,
                      timeout: float = 300) -> pd.DataFrame:
        """
        Run benchmarks on all matching files in the network folder.
        
        Args:
            file_pattern: Pattern to match network files
            max_files: Maximum number of files to process
            timeout: Maximum time per algorithm
            
        Returns:
            DataFrame with all results
        """
        print("="*80)
        print("MULTI-NETWORK SYNERGY ALGORITHM BENCHMARKING")
        print("="*80)
        
        # Find all matching files
        import glob
        pattern = os.path.join(self.network_folder, file_pattern)
        files = sorted(glob.glob(pattern))
        
        if max_files:
            files = files[:max_files]
        
        print(f"Found {len(files)} network files to benchmark")
        
        # Run benchmarks
        for i, file in enumerate(files, 1):
            print(f"\n[{i}/{len(files)}]", end="")
            result = self.benchmark_network(file, timeout)
            if result:
                self.results.append(result)
        
        # Convert to DataFrame
        df = pd.DataFrame(self.results)
        
        # Save results
        csv_path = os.path.join(self.output_folder, "benchmark_results.csv")
        df.to_csv(csv_path, index=False)
        print(f"\n\nResults saved to: {csv_path}")
        
        return df
    
    def analyze_scaling(self, df: pd.DataFrame) -> Dict:
        """Analyze how performance scales with network size."""
        print("\n" + "="*50)
        print("SCALING ANALYSIS")
        print("="*50)
        
        # Create figure for scaling plots
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        fig.suptitle('Algorithm Scaling Analysis', fontsize=16, fontweight='bold')
        
        # 1. Speedup vs network size (ERCs)
        ax = axes[0, 0]
        ax.scatter(df['n_ercs'], df['speedup'], alpha=0.6)
        ax.set_xlabel('Number of ERCs')
        ax.set_ylabel('Speedup Factor')
        ax.set_title('Speedup vs Network Size')
        ax.set_yscale('log')
        
        # Fit power law: speedup = a * n^b
        if len(df) > 3:
            x = np.log(df['n_ercs'])
            y = np.log(df['speedup'])
            slope, intercept, r_value, _, _ = stats.linregress(x, y)
            
            x_fit = np.linspace(df['n_ercs'].min(), df['n_ercs'].max(), 100)
            y_fit = np.exp(intercept) * x_fit**slope
            ax.plot(x_fit, y_fit, 'r-', label=f'y = {np.exp(intercept):.2f} * x^{slope:.2f}')
            ax.legend()
            
            print(f"Speedup scaling: speedup ∝ n_ercs^{slope:.2f} (R² = {r_value**2:.3f})")
        
        # 2. Time complexity analysis
        ax = axes[0, 1]
        ax.scatter(df['n_ercs'], df['brute_time'], label='Brute Force', alpha=0.6)
        ax.scatter(df['n_ercs'], df['efficient_time'], label='Efficient', alpha=0.6)
        ax.set_xlabel('Number of ERCs')
        ax.set_ylabel('Time (seconds)')
        ax.set_title('Time Complexity')
        ax.set_yscale('log')
        ax.legend()
        
        # 3. Operations per ERC
        ax = axes[0, 2]
        ax.scatter(df['n_ercs'], df['brute_steps_per_erc'], label='Brute Force', alpha=0.6)
        ax.scatter(df['n_ercs'], df['efficient_steps_per_erc'], label='Efficient', alpha=0.6)
        ax.set_xlabel('Number of ERCs')
        ax.set_ylabel('Operations per ERC')
        ax.set_title('Algorithm Efficiency')
        ax.set_yscale('log')
        ax.legend()
        
        # 4. Reduction ratio vs size
        ax = axes[1, 0]
        ax.scatter(df['n_ercs'], df['efficient_reduction'] * 100, alpha=0.6)
        ax.set_xlabel('Number of ERCs')
        ax.set_ylabel('Check Reduction (%)')
        ax.set_title('Efficiency Improvement with Size')
        
        # 5. Pruning effectiveness
        ax = axes[1, 1]
        total_pruning = df['within_target_pruning'] + df['cross_target_pruning']
        ax.scatter(df['n_ercs'], total_pruning, alpha=0.6)
        ax.set_xlabel('Number of ERCs')
        ax.set_ylabel('Total Pruning Operations')
        ax.set_title('Pruning Scalability')
        
        # 6. Theoretical vs actual operations
        ax = axes[1, 2]
        ax.scatter(df['n_ercs'], df['max_possible_checks'], label='Theoretical Max', alpha=0.6)
        ax.scatter(df['n_ercs'], df['brute_checks'], label='Brute Force', alpha=0.6)
        ax.scatter(df['n_ercs'], df['efficient_ops'], label='Efficient', alpha=0.6)
        ax.set_xlabel('Number of ERCs')
        ax.set_ylabel('Number of Operations')
        ax.set_title('Operations Scaling')
        ax.set_yscale('log')
        ax.legend()
        
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_folder, 'scaling_analysis.png'), dpi=300, bbox_inches='tight')
        
        # Calculate scaling exponents
        scaling_results = {}
        if len(df) > 3:
            # Brute force time scaling
            x = np.log(df['n_ercs'])
            y = np.log(df['brute_time'])
            slope, _, r_value, _, _ = stats.linregress(x, y)
            scaling_results['brute_time_exponent'] = slope
            scaling_results['brute_time_r2'] = r_value**2
            print(f"Brute force time: O(n^{slope:.2f})")
            
            # Efficient time scaling
            y = np.log(df['efficient_time'])
            slope, _, r_value, _, _ = stats.linregress(x, y)
            scaling_results['efficient_time_exponent'] = slope
            scaling_results['efficient_time_r2'] = r_value**2
            print(f"Efficient time: O(n^{slope:.2f})")
        
        return scaling_results
    
    def analyze_hierarchy_factors(self, df: pd.DataFrame):
        """Analyze which hierarchy properties affect performance."""
        print("\n" + "="*50)
        print("HIERARCHY FACTOR ANALYSIS")
        print("="*50)
        
        # Create correlation matrix
        factors = ['n_levels', 'avg_level_size', 'max_level_size', 'level_variance',
                  'hierarchy_density', 'max_chain_length', 'avg_branching', 
                  'leaf_ratio', 'root_ratio', 'hierarchy_depth_ratio']
        
        performance_metrics = ['speedup', 'efficient_reduction', 'pruning_per_synergy']
        
        # Calculate correlations
        correlations = {}
        for factor in factors:
            correlations[factor] = {}
            for metric in performance_metrics:
                if factor in df.columns and metric in df.columns:
                    corr, p_value = stats.pearsonr(df[factor], df[metric])
                    correlations[factor][metric] = (corr, p_value)
        
        # Create visualization
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle('Hierarchy Factors vs Performance', fontsize=16, fontweight='bold')
        
        # 1. Speedup vs hierarchy depth
        ax = axes[0, 0]
        scatter = ax.scatter(df['max_chain_length'], df['speedup'], 
                           c=df['n_ercs'], cmap='viridis', alpha=0.6)
        ax.set_xlabel('Maximum Chain Length')
        ax.set_ylabel('Speedup Factor')
        ax.set_title('Speedup vs Hierarchy Depth')
        plt.colorbar(scatter, ax=ax, label='Number of ERCs')
        
        # 2. Speedup vs number of levels
        ax = axes[0, 1]
        scatter = ax.scatter(df['n_levels'], df['speedup'],
                           c=df['n_ercs'], cmap='viridis', alpha=0.6)
        ax.set_xlabel('Number of Levels')
        ax.set_ylabel('Speedup Factor')
        ax.set_title('Speedup vs Number of Levels')
        plt.colorbar(scatter, ax=ax, label='Number of ERCs')
        
        # 3. Speedup vs hierarchy density
        ax = axes[1, 0]
        scatter = ax.scatter(df['hierarchy_density'], df['speedup'],
                           c=df['n_ercs'], cmap='viridis', alpha=0.6)
        ax.set_xlabel('Hierarchy Density')
        ax.set_ylabel('Speedup Factor')
        ax.set_title('Speedup vs Hierarchy Density')
        plt.colorbar(scatter, ax=ax, label='Number of ERCs')
        
        # 4. Correlation heatmap
        ax = axes[1, 1]
        corr_matrix = []
        factor_labels = []
        metric_labels = []
        
        for factor in factors[:6]:  # Limit to 6 for readability
            if factor in df.columns:
                row = []
                factor_labels.append(factor)
                for metric in performance_metrics:
                    if metric in df.columns:
                        corr, _ = correlations[factor][metric]
                        row.append(corr)
                        if factor == factors[0]:
                            metric_labels.append(metric)
                corr_matrix.append(row)
        
        if corr_matrix:
            im = ax.imshow(corr_matrix, cmap='RdBu', vmin=-1, vmax=1, aspect='auto')
            ax.set_xticks(range(len(metric_labels)))
            ax.set_yticks(range(len(factor_labels)))
            ax.set_xticklabels(metric_labels, rotation=45, ha='right')
            ax.set_yticklabels(factor_labels)
            ax.set_title('Factor-Performance Correlations')
            plt.colorbar(im, ax=ax)
            
            # Add correlation values
            for i in range(len(factor_labels)):
                for j in range(len(metric_labels)):
                    text = ax.text(j, i, f'{corr_matrix[i][j]:.2f}',
                                 ha="center", va="center", color="black" if abs(corr_matrix[i][j]) < 0.5 else "white")
        
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_folder, 'hierarchy_factors.png'), dpi=300, bbox_inches='tight')
        
        # Print significant correlations
        print("\nSignificant correlations (p < 0.05):")
        for factor in factors:
            for metric in performance_metrics:
                if factor in correlations and metric in correlations[factor]:
                    corr, p_value = correlations[factor][metric]
                    if p_value < 0.05:
                        print(f"  {factor} vs {metric}: r = {corr:.3f} (p = {p_value:.3f})")
        
        return correlations
    
    def create_summary_report(self, df: pd.DataFrame):
        """Create a comprehensive summary report."""
        print("\n" + "="*50)
        print("SUMMARY REPORT")
        print("="*50)
        
        # Overall statistics
        print(f"\nDataset Overview:")
        print(f"  Total networks analyzed: {len(df)}")
        print(f"  Network size range: {df['n_ercs'].min()}-{df['n_ercs'].max()} ERCs")
        print(f"  All results matched: {df['results_match'].all()}")
        
        print(f"\nPerformance Summary:")
        print(f"  Average speedup: {df['speedup'].mean():.2f}x (±{df['speedup'].std():.2f})")
        print(f"  Median speedup: {df['speedup'].median():.2f}x")
        print(f"  Maximum speedup: {df['speedup'].max():.2f}x")
        print(f"  Average check reduction: {(df['efficient_reduction'].mean()*100):.1f}%")
        
        print(f"\nScaling Behavior:")
        print(f"  Speedup increases with network size: ", end="")
        corr, p = stats.pearsonr(df['n_ercs'], df['speedup'])
        print(f"r = {corr:.3f} (p = {p:.3f})")
        
        # Create summary visualization
        fig = plt.figure(figsize=(16, 10))
        
        # 1. Speedup distribution
        ax1 = plt.subplot(2, 3, 1)
        ax1.hist(df['speedup'], bins=20, edgecolor='black', alpha=0.7)
        ax1.axvline(df['speedup'].mean(), color='red', linestyle='--', label=f'Mean: {df["speedup"].mean():.2f}x')
        ax1.set_xlabel('Speedup Factor')
        ax1.set_ylabel('Count')
        ax1.set_title('Speedup Distribution')
        ax1.legend()
        
        # 2. Network size distribution
        ax2 = plt.subplot(2, 3, 2)
        ax2.hist(df['n_ercs'], bins=20, edgecolor='black', alpha=0.7)
        ax2.set_xlabel('Number of ERCs')
        ax2.set_ylabel('Count')
        ax2.set_title('Network Size Distribution')
        
        # 3. Speedup vs size with trend
        ax3 = plt.subplot(2, 3, 3)
        ax3.scatter(df['n_ercs'], df['speedup'], alpha=0.6)
        z = np.polyfit(df['n_ercs'], df['speedup'], 2)
        p = np.poly1d(z)
        x_trend = np.linspace(df['n_ercs'].min(), df['n_ercs'].max(), 100)
        ax3.plot(x_trend, p(x_trend), "r--", alpha=0.8, label='Trend')
        ax3.set_xlabel('Number of ERCs')
        ax3.set_ylabel('Speedup Factor')
        ax3.set_title('Speedup Scaling Trend')
        ax3.legend()
        
        # 4. Time comparison
        ax4 = plt.subplot(2, 3, 4)
        ax4.boxplot([df['brute_time'], df['efficient_time']], 
                   labels=['Brute Force', 'Efficient'])
        ax4.set_ylabel('Time (seconds)')
        ax4.set_title('Algorithm Time Comparison')
        ax4.set_yscale('log')
        
        # 5. Efficiency metrics
        ax5 = plt.subplot(2, 3, 5)
        metrics = ['Speedup', 'Check\nReduction', 'Pruning\nEfficiency']
        values = [
            df['speedup'].mean(),
            df['efficient_reduction'].mean() * 100,
            (df['pruning_per_synergy'].mean() / df['n_ercs'].mean()) * 100
        ]
        bars = ax5.bar(metrics, values, color=['green', 'blue', 'orange'])
        ax5.set_ylabel('Value')
        ax5.set_title('Average Efficiency Metrics')
        
        # Add value labels
        for bar, value in zip(bars, values):
            height = bar.get_height()
            ax5.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                    f'{value:.1f}', ha='center', va='bottom')
        
        # 6. Network properties impact
        ax6 = plt.subplot(2, 3, 6)
        # Show top factors affecting speedup
        factors = ['n_levels', 'max_chain_length', 'hierarchy_density', 'avg_branching']
        correlations = []
        for factor in factors:
            if factor in df.columns:
                corr, _ = stats.pearsonr(df[factor], df['speedup'])
                correlations.append(abs(corr))
            else:
                correlations.append(0)
        
        bars = ax6.bar(range(len(factors)), correlations, color='purple', alpha=0.7)
        ax6.set_xticks(range(len(factors)))
        ax6.set_xticklabels(factors, rotation=45, ha='right')
        ax6.set_ylabel('|Correlation| with Speedup')
        ax6.set_title('Hierarchy Factors Impact')
        
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_folder, 'summary_report.png'), dpi=300, bbox_inches='tight')
        
        # Save detailed statistics
        stats_file = os.path.join(self.output_folder, 'summary_statistics.txt')
        with open(stats_file, 'w') as f:
            f.write("SYNERGY ALGORITHM BENCHMARKING SUMMARY\n")
            f.write("="*50 + "\n\n")
            f.write(f"Total networks analyzed: {len(df)}\n")
            f.write(f"Network size range: {df['n_ercs'].min()}-{df['n_ercs'].max()} ERCs\n")
            f.write(f"Average speedup: {df['speedup'].mean():.2f}x (±{df['speedup'].std():.2f})\n")
            f.write(f"Maximum speedup: {df['speedup'].max():.2f}x\n")
            f.write(f"Average check reduction: {(df['efficient_reduction'].mean()*100):.1f}%\n")
            f.write(f"\nScaling behavior: speedup ∝ n_ercs^{slope:.2f}\n" if 'slope' in locals() else "")
        
        print(f"\nResults saved to: {self.output_folder}")
    
    def run_complete_benchmark(self, file_pattern: str = "*.txt", 
                              max_files: Optional[int] = None,
                              timeout: float = 300):
        """Run complete benchmarking analysis."""
        # Run benchmarks
        df = self.run_benchmarks(file_pattern, max_files, timeout)
        
        if len(df) == 0:
            print("No successful benchmarks to analyze!")
            return None
        
        # Perform analyses
        scaling_results = self.analyze_scaling(df)
        hierarchy_correlations = self.analyze_hierarchy_factors(df)
        self.create_summary_report(df)
        
        # Show all plots
        plt.show()
        
        return {
            'dataframe': df,
            'scaling_results': scaling_results,
            'hierarchy_correlations': hierarchy_correlations,
            'network_data': self.network_data
        }


def benchmark_synergy_algorithms(network_folder: str, 
                               output_folder: str = "benchmark_results",
                               file_pattern: str = "*.txt",
                               max_files: Optional[int] = None,
                               timeout: float = 300):
    """
    Main function to benchmark synergy algorithms across multiple networks.
    
    Args:
        network_folder: Path to folder containing network files
        output_folder: Path to save results
        file_pattern: Pattern to match network files
        max_files: Maximum number of files to process
        timeout: Maximum time per algorithm (seconds)
    
    Returns:
        Dictionary with all results and analyses
    """
    benchmarker = NetworkBenchmarker(network_folder, output_folder)
    results = benchmarker.run_complete_benchmark(file_pattern, max_files, timeout)
    return benchmarker, results


if __name__ == "__main__":
    # Example usage
    import sys
    
    if len(sys.argv) > 1:
        network_folder = sys.argv[1]
    else:
        # Default folder
        network_folder = "networks/RandomAlife/"
    
    print(f"Benchmarking networks in: {network_folder}")
    
    # Run benchmarks
    benchmarker, results = benchmark_synergy_algorithms(
        network_folder,
        output_folder="benchmark_results",
        file_pattern="*.txt",
        max_files=None,  # Process all files
        timeout=600  # 5 minute timeout per algorithm
    )
    
    # Keep plots open9
    plt.show(block=True)