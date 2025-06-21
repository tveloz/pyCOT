#!/usr/bin/env python3
"""
Dynamical Hierarchy Algorithm Analysis Tool

This script analyzes the performance and efficiency of the Dynamical Hierarchy 
synergy algorithm across multiple reaction networks, providing detailed insights
into network structure, hierarchy properties, and algorithm behavior.

@author: DH Algorithm Analysis
"""

import os
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from collections import defaultdict, Counter
import networkx as nx
from typing import List, Dict, Tuple, Optional
import warnings
warnings.filterwarnings('ignore')

# Import required modules
from pyCOT.io.functions import read_txt
from pyCOT.ERC_Hierarchy import ERC, species_list_to_names
from pyCOT.Dynamical_Hierarchy_Synergy_Algorithm import EfficientFundamentalSynergyFinder


class DHNetworkAnalyzer:
    """Analyzes Dynamical Hierarchy algorithm performance across multiple networks."""
    
    def __init__(self, network_folder: str, output_folder: str = "dh_analysis_results"):
        """
        Initialize analyzer.
        
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
    
    def analyze_network_structure(self, rn, ercs, hierarchy_graph) -> Dict:
        """Extract comprehensive network structure properties."""
        # Basic network properties
        n_species = len(rn.species())
        n_reactions = len(rn.reactions())
        n_ercs = len(ercs)
        
        # Reaction complexity metrics
        reaction_sizes = [len(r.support_names()) + len(r.products_names()) for r in rn.reactions()]
        avg_reaction_size = np.mean(reaction_sizes) if reaction_sizes else 0
        max_reaction_size = max(reaction_sizes) if reaction_sizes else 0
        
        # Species connectivity
        species_reaction_count = defaultdict(int)
        for reaction in rn.reactions():
            for species in reaction.species_names():
                species_reaction_count[species] += 1
        
        connectivity_values = list(species_reaction_count.values())
        avg_species_connectivity = np.mean(connectivity_values) if connectivity_values else 0
        max_species_connectivity = max(connectivity_values) if connectivity_values else 0
        
        # Network topology metrics
        reaction_graph = nx.Graph()
        for reaction in rn.reactions():
            species = reaction.species_names()
            for i, s1 in enumerate(species):
                for s2 in species[i+1:]:
                    reaction_graph.add_edge(s1, s2)
        
        # Graph connectivity metrics
        is_connected = nx.is_connected(reaction_graph) if reaction_graph.nodes() else False
        n_components = nx.number_connected_components(reaction_graph) if reaction_graph.nodes() else 0
        
        # Clustering coefficient
        clustering = nx.average_clustering(reaction_graph) if reaction_graph.nodes() else 0
        
        # Density
        density = nx.density(reaction_graph) if reaction_graph.nodes() else 0
        
        return {
            'n_species': n_species,
            'n_reactions': n_reactions,
            'n_ercs': n_ercs,
            'avg_reaction_size': avg_reaction_size,
            'max_reaction_size': max_reaction_size,
            'avg_species_connectivity': avg_species_connectivity,
            'max_species_connectivity': max_species_connectivity,
            'is_connected': is_connected,
            'n_components': n_components,
            'clustering_coefficient': clustering,
            'network_density': density,
            'species_per_reaction': n_species / n_reactions if n_reactions > 0 else 0,
            'reactions_per_species': n_reactions / n_species if n_species > 0 else 0,
            'ercs_per_species': n_ercs / n_species if n_species > 0 else 0,
            'ercs_per_reaction': n_ercs / n_reactions if n_reactions > 0 else 0,
        }
    
    def analyze_hierarchy_structure(self, ercs, hierarchy_graph) -> Dict:
        """Extract detailed hierarchy structure properties."""
        n_ercs = len(ercs)
        
        # Level distribution
        levels = ERC.get_node_levels(hierarchy_graph)
        level_distribution = defaultdict(int)
        for erc_label, level in levels.items():
            level_distribution[level] += 1
        
        n_levels = len(level_distribution)
        level_counts = list(level_distribution.values())
        
        # Hierarchy shape metrics
        avg_level_size = np.mean(level_counts) if level_counts else 0
        max_level_size = max(level_counts) if level_counts else 0
        min_level_size = min(level_counts) if level_counts else 0
        level_variance = np.var(level_counts) if level_counts else 0
        level_skewness = stats.skew(level_counts) if len(level_counts) > 2 else 0
        
        # Hierarchy connectivity
        n_edges = hierarchy_graph.number_of_edges()
        avg_degree = 2 * n_edges / n_ercs if n_ercs > 0 else 0
        
        # Node role distribution
        leaf_nodes = sum(1 for node in hierarchy_graph.nodes() 
                        if hierarchy_graph.out_degree(node) == 0)
        root_nodes = sum(1 for node in hierarchy_graph.nodes() 
                        if hierarchy_graph.in_degree(node) == 0)
        intermediate_nodes = n_ercs - leaf_nodes - root_nodes
        
        # Branching analysis
        out_degrees = [hierarchy_graph.out_degree(node) for node in hierarchy_graph.nodes()]
        in_degrees = [hierarchy_graph.in_degree(node) for node in hierarchy_graph.nodes()]
        
        avg_out_degree = np.mean(out_degrees) if out_degrees else 0
        max_out_degree = max(out_degrees) if out_degrees else 0
        avg_in_degree = np.mean(in_degrees) if in_degrees else 0
        max_in_degree = max(in_degrees) if in_degrees else 0
        
        # Hierarchy depth metrics
        max_chain_length = 0
        total_path_lengths = []
        for node in hierarchy_graph.nodes():
            ancestors = nx.ancestors(hierarchy_graph, node)
            path_length = len(ancestors)
            total_path_lengths.append(path_length)
            if path_length > max_chain_length:
                max_chain_length = path_length
        
        avg_depth = np.mean(total_path_lengths) if total_path_lengths else 0
        
        # Hierarchy balance metrics
        hierarchy_density = n_edges / (n_ercs * (n_ercs - 1) / 2) if n_ercs > 1 else 0
        
        # Width-to-height ratio
        width_height_ratio = max_level_size / n_levels if n_levels > 0 else 0
        
        return {
            'n_levels': n_levels,
            'avg_level_size': avg_level_size,
            'max_level_size': max_level_size,
            'min_level_size': min_level_size,
            'level_variance': level_variance,
            'level_skewness': level_skewness,
            'n_edges': n_edges,
            'avg_degree': avg_degree,
            'hierarchy_density': hierarchy_density,
            'leaf_nodes': leaf_nodes,
            'root_nodes': root_nodes,
            'intermediate_nodes': intermediate_nodes,
            'max_chain_length': max_chain_length,
            'avg_depth': avg_depth,
            'avg_out_degree': avg_out_degree,
            'max_out_degree': max_out_degree,
            'avg_in_degree': avg_in_degree,
            'max_in_degree': max_in_degree,
            'width_height_ratio': width_height_ratio,
            'leaf_ratio': leaf_nodes / n_ercs if n_ercs > 0 else 0,
            'root_ratio': root_nodes / n_ercs if n_ercs > 0 else 0,
            'intermediate_ratio': intermediate_nodes / n_ercs if n_ercs > 0 else 0,
            'hierarchy_depth_ratio': max_chain_length / n_ercs if n_ercs > 0 else 0,
        }
    
    def analyze_erc_properties(self, ercs, rn) -> Dict:
        """Analyze ERC-specific properties."""
        # Generator analysis
        generator_counts = [len(erc.min_generators) for erc in ercs]
        avg_generators = np.mean(generator_counts) if generator_counts else 0
        max_generators = max(generator_counts) if generator_counts else 0
        min_generators = min(generator_counts) if generator_counts else 0
        generator_variance = np.var(generator_counts) if generator_counts else 0
        
        # Closure size analysis
        closure_sizes = []
        for erc in ercs:
            closure_species = erc.get_closure(rn)
            closure_sizes.append(len(closure_species))
        
        avg_closure_size = np.mean(closure_sizes) if closure_sizes else 0
        max_closure_size = max(closure_sizes) if closure_sizes else 0
        min_closure_size = min(closure_sizes) if closure_sizes else 0
        closure_variance = np.var(closure_sizes) if closure_sizes else 0
        
        # Generator size distribution
        generator_sizes = []
        for erc in ercs:
            for gen in erc.min_generators:
                generator_sizes.append(len(gen))
        
        avg_generator_size = np.mean(generator_sizes) if generator_sizes else 0
        max_generator_size = max(generator_sizes) if generator_sizes else 0
        
        # ERC complexity metrics
        complexity_ratios = []
        for i, erc in enumerate(ercs):
            closure_size = closure_sizes[i]
            generator_count = generator_counts[i]
            if generator_count > 0:
                complexity_ratios.append(closure_size / generator_count)
        
        avg_complexity_ratio = np.mean(complexity_ratios) if complexity_ratios else 0
        
        return {
            'avg_generators': avg_generators,
            'max_generators': max_generators,
            'min_generators': min_generators,
            'generator_variance': generator_variance,
            'avg_closure_size': avg_closure_size,
            'max_closure_size': max_closure_size,
            'min_closure_size': min_closure_size,
            'closure_variance': closure_variance,
            'avg_generator_size': avg_generator_size,
            'max_generator_size': max_generator_size,
            'avg_complexity_ratio': avg_complexity_ratio,
        }
    
    def analyze_network(self, network_file: str, timeout: Optional[float] = 300) -> Optional[Dict]:
        """
        Analyze a single network with the DH algorithm.
        
        Args:
            network_file: Path to network file
            timeout: Maximum time allowed
            
        Returns:
            Dictionary with analysis results or None if error/timeout
        """
        print(f"\nAnalyzing: {os.path.basename(network_file)}")
        
        try:
            # Load network
            rn = read_txt(network_file)
            ercs = ERC.ERCs(rn)
            hierarchy_graph = ERC.build_hierarchy_graph(ercs, rn)
            
            # Get comprehensive network properties
            network_props = self.analyze_network_structure(rn, ercs, hierarchy_graph)
            hierarchy_props = self.analyze_hierarchy_structure(ercs, hierarchy_graph)
            erc_props = self.analyze_erc_properties(ercs, rn)
            
            print(f"  Network: {network_props['n_species']} species, {network_props['n_reactions']} reactions, {network_props['n_ercs']} ERCs")
            print(f"  Hierarchy: {hierarchy_props['n_levels']} levels, max chain {hierarchy_props['max_chain_length']}")
            
            # Run DH algorithm
            print("  Running DH algorithm...", end="", flush=True)
            dh_finder = EfficientFundamentalSynergyFinder(rn, ercs, hierarchy_graph)
            
            start_time = time.time()
            fundamental_synergies = dh_finder.find_fundamental_synergies(verbose=False)
            execution_time = time.time() - start_time
            
            if execution_time > timeout:
                print(f" TIMEOUT ({execution_time:.1f}s)")
                return None
            
            print(f" {execution_time:.3f}s")
            
            # Analyze DH algorithm performance
            stats = dh_finder.stats
            
            # Calculate efficiency metrics
            n_ercs = len(ercs)
            theoretical_max = n_ercs * (n_ercs - 1) * (n_ercs - 2) // 2
            
            # Operations efficiency
            total_operations = (stats['partial_overlap_checks'] + 
                              stats['base_pairs_generated'] + 
                              stats['synergies_found'])
            
            efficiency_ratio = total_operations / theoretical_max if theoretical_max > 0 else 0
            operations_per_erc = total_operations / n_ercs if n_ercs > 0 else 0
            
            # Pruning effectiveness
            total_pruning = stats['within_target_pruning'] + stats['cross_target_pruning']
            pruning_ratio = total_pruning / stats['base_pairs_generated'] if stats['base_pairs_generated'] > 0 else 0
            
            # Synergy density
            synergy_density = len(fundamental_synergies) / theoretical_max if theoretical_max > 0 else 0
            synergies_per_erc = len(fundamental_synergies) / n_ercs if n_ercs > 0 else 0
            
            # Base-target relationship efficiency
            base_target_efficiency = (stats['valid_base_target_pairs'] / 
                                    stats['partial_overlap_checks'] if stats['partial_overlap_checks'] > 0 else 0)
            
            # Success rates
            synergy_success_rate = (stats['synergies_found'] / 
                                  stats['base_pairs_generated'] if stats['base_pairs_generated'] > 0 else 0)
            
            result = {
                'file': os.path.basename(network_file),
                
                # Network structure properties
                **network_props,
                
                # Hierarchy properties
                **hierarchy_props,
                
                # ERC properties
                **erc_props,
                
                # DH Algorithm performance
                'execution_time': execution_time,
                'fundamental_synergies': len(fundamental_synergies),
                
                # Raw statistics
                'partial_overlap_checks': stats['partial_overlap_checks'],
                'valid_base_target_pairs': stats['valid_base_target_pairs'],
                'base_pairs_generated': stats['base_pairs_generated'],
                'synergies_found': stats['synergies_found'],
                'within_target_pruning': stats['within_target_pruning'],
                'cross_target_pruning': stats['cross_target_pruning'],
                'final_fundamental': stats['final_fundamental'],
                
                # Efficiency metrics
                'theoretical_max_checks': theoretical_max,
                'total_operations': total_operations,
                'efficiency_ratio': efficiency_ratio,
                'operations_per_erc': operations_per_erc,
                'total_pruning': total_pruning,
                'pruning_ratio': pruning_ratio,
                'synergy_density': synergy_density,
                'synergies_per_erc': synergies_per_erc,
                'base_target_efficiency': base_target_efficiency,
                'synergy_success_rate': synergy_success_rate,
                
                # Performance ratios
                'time_per_erc': execution_time / n_ercs if n_ercs > 0 else 0,
                'time_per_operation': execution_time / total_operations if total_operations > 0 else 0,
                'reduction_factor': theoretical_max / total_operations if total_operations > 0 else 0,
            }
            
            # Store network data for later analysis
            self.network_data[network_file] = {
                'rn': rn,
                'ercs': ercs,
                'hierarchy_graph': hierarchy_graph,
                'fundamental_synergies': fundamental_synergies,
                'dh_finder': dh_finder
            }
            
            return result
            
        except Exception as e:
            print(f"  ERROR: {str(e)}")
            import traceback
            traceback.print_exc()
            return None
    
    def run_analysis(self, file_pattern: str = "*.txt", max_files: Optional[int] = None,
                    timeout: float = 300) -> pd.DataFrame:
        """
        Run analysis on all matching files in the network folder.
        
        Args:
            file_pattern: Pattern to match network files
            max_files: Maximum number of files to process
            timeout: Maximum time per analysis
            
        Returns:
            DataFrame with all results
        """
        print("="*80)
        print("DYNAMICAL HIERARCHY ALGORITHM NETWORK ANALYSIS")
        print("="*80)
        
        # Find all matching files
        import glob
        pattern = os.path.join(self.network_folder, file_pattern)
        files = sorted(glob.glob(pattern))
        
        if max_files:
            files = files[:max_files]
        
        print(f"Found {len(files)} network files to analyze")
        
        # Run analyses
        for i, file in enumerate(files, 1):
            print(f"\n[{i}/{len(files)}]", end="")
            result = self.analyze_network(file, timeout)
            if result:
                self.results.append(result)
        
        # Convert to DataFrame
        df = pd.DataFrame(self.results)
        
        if len(df) == 0:
            print("\nNo successful analyses to process!")
            return df
        
        # Save results
        csv_path = os.path.join(self.output_folder, "dh_analysis_results.csv")
        df.to_csv(csv_path, index=False)
        print(f"\n\nResults saved to: {csv_path}")
        
        return df
    
    def analyze_performance_scaling(self, df: pd.DataFrame) -> Dict:
        """Analyze how DH performance scales with network properties."""
        print("\n" + "="*60)
        print("PERFORMANCE SCALING ANALYSIS")
        print("="*60)
        
        if len(df) == 0:
            print("No data to analyze!")
            return {}
        
        # Create comprehensive scaling analysis
        fig, axes = plt.subplots(3, 3, figsize=(20, 15))
        fig.suptitle('DH Algorithm Performance Scaling Analysis', fontsize=16, fontweight='bold')
        
        # 1. Execution time vs network size
        ax = axes[0, 0]
        ax.scatter(df['n_ercs'], df['execution_time'], alpha=0.6, color='blue')
        ax.set_xlabel('Number of ERCs')
        ax.set_ylabel('Execution Time (seconds)')
        ax.set_title('Execution Time vs Network Size')
        ax.set_yscale('log')
        
        # 2. Operations vs network size
        ax = axes[0, 1]
        ax.scatter(df['n_ercs'], df['total_operations'], alpha=0.6, color='green')
        ax.set_xlabel('Number of ERCs')
        ax.set_ylabel('Total Operations')
        ax.set_title('Operations vs Network Size')
        ax.set_yscale('log')
        
        # 3. Efficiency ratio vs network size
        ax = axes[0, 2]
        ax.scatter(df['n_ercs'], df['efficiency_ratio'] * 100, alpha=0.6, color='red')
        ax.set_xlabel('Number of ERCs')
        ax.set_ylabel('Efficiency Ratio (%)')
        ax.set_title('Algorithm Efficiency vs Size')
        
        # 4. Synergy discovery vs hierarchy depth
        ax = axes[1, 0]
        scatter = ax.scatter(df['max_chain_length'], df['fundamental_synergies'], 
                           c=df['n_ercs'], cmap='viridis', alpha=0.6)
        ax.set_xlabel('Maximum Chain Length')
        ax.set_ylabel('Fundamental Synergies Found')
        ax.set_title('Synergy Discovery vs Hierarchy Depth')
        plt.colorbar(scatter, ax=ax, label='Number of ERCs')
        
        # 5. Pruning effectiveness vs hierarchy structure
        ax = axes[1, 1]
        ax.scatter(df['hierarchy_density'], df['pruning_ratio'] * 100, alpha=0.6, color='orange')
        ax.set_xlabel('Hierarchy Density')
        ax.set_ylabel('Pruning Ratio (%)')
        ax.set_title('Pruning Effectiveness vs Hierarchy Density')
        
        # 6. Base-target efficiency vs network connectivity
        ax = axes[1, 2]
        ax.scatter(df['clustering_coefficient'], df['base_target_efficiency'] * 100, alpha=0.6, color='purple')
        ax.set_xlabel('Network Clustering Coefficient')
        ax.set_ylabel('Base-Target Efficiency (%)')
        ax.set_title('Base-Target Efficiency vs Network Structure')
        
        # 7. Time per ERC vs complexity
        ax = axes[2, 0]
        ax.scatter(df['avg_complexity_ratio'], df['time_per_erc'], alpha=0.6, color='brown')
        ax.set_xlabel('Average ERC Complexity Ratio')
        ax.set_ylabel('Time per ERC (seconds)')
        ax.set_title('Time Efficiency vs ERC Complexity')
        ax.set_yscale('log')
        
        # 8. Synergy density vs network properties
        ax = axes[2, 1]
        ax.scatter(df['network_density'], df['synergy_density'] * 100, alpha=0.6, color='pink')
        ax.set_xlabel('Network Density')
        ax.set_ylabel('Synergy Density (%)')
        ax.set_title('Synergy Density vs Network Density')
        
        # 9. Reduction factor vs hierarchy levels
        ax = axes[2, 2]
        ax.scatter(df['n_levels'], df['reduction_factor'], alpha=0.6, color='cyan')
        ax.set_xlabel('Number of Hierarchy Levels')
        ax.set_ylabel('Computational Reduction Factor')
        ax.set_title('Reduction Factor vs Hierarchy Complexity')
        ax.set_yscale('log')
        
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_folder, 'performance_scaling.png'), 
                   dpi=300, bbox_inches='tight')
        
        # Calculate scaling relationships
        scaling_results = {}
        if len(df) > 3:
            # Time complexity
            x = np.log(df['n_ercs'])
            y = np.log(df['execution_time'])
            slope, _, r_value, _, _ = stats.linregress(x, y)
            scaling_results['time_scaling_exponent'] = slope
            scaling_results['time_scaling_r2'] = r_value**2
            print(f"Time scaling: O(n^{slope:.3f}) with R² = {r_value**2:.3f}")
            
            # Operations scaling
            y = np.log(df['total_operations'])
            slope, _, r_value, _, _ = stats.linregress(x, y)
            scaling_results['operations_scaling_exponent'] = slope
            scaling_results['operations_scaling_r2'] = r_value**2
            print(f"Operations scaling: O(n^{slope:.3f}) with R² = {r_value**2:.3f}")
        
        return scaling_results
    
    def analyze_hierarchy_impact(self, df: pd.DataFrame):
        """Analyze how hierarchy properties impact algorithm performance."""
        print("\n" + "="*60)
        print("HIERARCHY IMPACT ANALYSIS")
        print("="*60)
        
        if len(df) == 0:
            print("No data to analyze!")
            return {}
        
        # Define hierarchy factors and performance metrics
        hierarchy_factors = [
            'n_levels', 'max_chain_length', 'hierarchy_density', 'avg_depth',
            'width_height_ratio', 'level_variance', 'avg_out_degree', 'leaf_ratio'
        ]
        
        performance_metrics = [
            'execution_time', 'efficiency_ratio', 'pruning_ratio', 
            'synergy_success_rate', 'reduction_factor'
        ]
        
        # Calculate correlations
        correlations = {}
        for factor in hierarchy_factors:
            correlations[factor] = {}
            for metric in performance_metrics:
                if factor in df.columns and metric in df.columns:
                    corr, p_value = stats.pearsonr(df[factor], df[metric])
                    correlations[factor][metric] = (corr, p_value)
        
        # Create visualization
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle('Hierarchy Impact on DH Algorithm Performance', fontsize=16, fontweight='bold')
        
        # 1. Correlation heatmap
        ax = axes[0, 0]
        corr_matrix = []
        factor_labels = []
        metric_labels = []
        
        for factor in hierarchy_factors:
            if factor in df.columns:
                row = []
                factor_labels.append(factor.replace('_', '\n'))
                for metric in performance_metrics:
                    if metric in df.columns:
                        corr, _ = correlations[factor][metric]
                        row.append(corr)
                        if factor == hierarchy_factors[0]:
                            metric_labels.append(metric.replace('_', '\n'))
                corr_matrix.append(row)
        
        if corr_matrix:
            im = ax.imshow(corr_matrix, cmap='RdBu', vmin=-1, vmax=1, aspect='auto')
            ax.set_xticks(range(len(metric_labels)))
            ax.set_yticks(range(len(factor_labels)))
            ax.set_xticklabels(metric_labels, rotation=45, ha='right')
            ax.set_yticklabels(factor_labels)
            ax.set_title('Hierarchy-Performance Correlations')
            plt.colorbar(im, ax=ax)
            
            # Add correlation values
            for i in range(len(factor_labels)):
                for j in range(len(metric_labels)):
                    value = corr_matrix[i][j]
                    color = "white" if abs(value) > 0.5 else "black"
                    ax.text(j, i, f'{value:.2f}', ha="center", va="center", color=color)
        
        # 2. Key relationships
        ax = axes[0, 1]
        ax.scatter(df['max_chain_length'], df['execution_time'], alpha=0.6)
        ax.set_xlabel('Maximum Chain Length')
        ax.set_ylabel('Execution Time (seconds)')
        ax.set_title('Performance vs Hierarchy Depth')
        ax.set_yscale('log')
        
        # 3. Efficiency vs structure
        ax = axes[1, 0]
        ax.scatter(df['hierarchy_density'], df['efficiency_ratio'], alpha=0.6)
        ax.set_xlabel('Hierarchy Density')
        ax.set_ylabel('Efficiency Ratio')
        ax.set_title('Efficiency vs Hierarchy Density')
        
        # 4. Pruning vs hierarchy balance
        ax = axes[1, 1]
        ax.scatter(df['level_variance'], df['pruning_ratio'], alpha=0.6)
        ax.set_xlabel('Level Variance (Hierarchy Balance)')
        ax.set_ylabel('Pruning Ratio')
        ax.set_title('Pruning Effectiveness vs Hierarchy Balance')
        
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_folder, 'hierarchy_impact.png'), 
                   dpi=300, bbox_inches='tight')
        
        # Print significant correlations
        print("\nSignificant correlations (p < 0.05):")
        for factor in hierarchy_factors:
            for metric in performance_metrics:
                if factor in correlations and metric in correlations[factor]:
                    corr, p_value = correlations[factor][metric]
                    if p_value < 0.05:
                        strength = "strong" if abs(corr) > 0.7 else "moderate" if abs(corr) > 0.5 else "weak"
                        direction = "positive" if corr > 0 else "negative"
                        print(f"  {factor} vs {metric}: {direction} {strength} correlation "
                              f"(r = {corr:.3f}, p = {p_value:.3f})")
        
        return correlations
    
    def analyze_algorithm_efficiency(self, df: pd.DataFrame):
        """Analyze detailed algorithm efficiency patterns."""
        print("\n" + "="*60)
        print("ALGORITHM EFFICIENCY ANALYSIS")
        print("="*60)
        
        if len(df) == 0:
            print("No data to analyze!")
            return
        
        # Create efficiency analysis plots
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        fig.suptitle('DH Algorithm Efficiency Analysis', fontsize=16, fontweight='bold')
        
        # 1. Algorithm stage breakdown
        ax = axes[0, 0]
        stages = ['Partial Overlap\nChecks', 'Base Pairs\nGenerated', 'Synergies\nFound', 
                 'Pruning\nOperations']
        avg_operations = [
            df['partial_overlap_checks'].mean(),
            df['base_pairs_generated'].mean(),
            df['synergies_found'].mean(),
            df['total_pruning'].mean()
        ]
        bars = ax.bar(stages, avg_operations, color=['blue', 'green', 'red', 'orange'])
        ax.set_ylabel('Average Operations')
        ax.set_title('Algorithm Stage Breakdown')
        ax.set_yscale('log')
        
        # Add value labels
        for bar, value in zip(bars, avg_operations):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height * 1.1,
                   f'{value:.1f}', ha='center', va='bottom')
        
        # 2. Efficiency distribution
        ax = axes[0, 1]
        ax.hist(df['efficiency_ratio'] * 100, bins=20, edgecolor='black', alpha=0.7)
        ax.axvline(df['efficiency_ratio'].mean() * 100, color='red', linestyle='--', 
                  label=f'Mean: {df["efficiency_ratio"].mean()*100:.2f}%')
        ax.set_xlabel('Efficiency Ratio (%)')
        ax.set_ylabel('Count')
        ax.set_title('Algorithm Efficiency Distribution')
        ax.legend()
        
        # 3. Success rate analysis
        ax = axes[0, 2]
        ax.scatter(df['base_pairs_generated'], df['synergy_success_rate'] * 100, alpha=0.6)
        ax.set_xlabel('Base Pairs Generated')
        ax.set_ylabel('Synergy Success Rate (%)')
        ax.set_title('Synergy Discovery Success Rate')
        ax.set_xscale('log')
        
        # 4. Pruning effectiveness
        ax = axes[1, 0]
        pruning_data = [df['within_target_pruning'], df['cross_target_pruning']]
        labels = ['Within Target', 'Cross Target']
        bp = ax.boxplot(pruning_data, labels=labels, patch_artist=True)
        colors = ['lightblue', 'lightgreen']
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
        ax.set_ylabel('Pruning Operations')
        ax.set_title('Pruning Strategy Effectiveness')
        ax.set_yscale('log')
        
        # 5. Time vs operations relationship
        ax = axes[1, 1]
        ax.scatter(df['total_operations'], df['execution_time'], alpha=0.6)
        ax.set_xlabel('Total Operations')
        ax.set_ylabel('Execution Time (seconds)')
        ax.set_title('Time vs Operations Relationship')
        ax.set_xscale('log')
        ax.set_yscale('log')
        
        # Add trend line
        if len(df) > 2:
            z = np.polyfit(np.log(df['total_operations']), np.log(df['execution_time']), 1)
            p = np.poly1d(z)
            x_trend = np.logspace(np.log10(df['total_operations'].min()), 
                                 np.log10(df['total_operations'].max()), 100)
            y_trend = np.exp(p(np.log(x_trend)))
            ax.plot(x_trend, y_trend, "r--", alpha=0.8, label=f'Slope: {z[0]:.2f}')
            ax.legend()
        
        # 6. Reduction factor distribution
        ax = axes[1, 2]
        ax.hist(df['reduction_factor'], bins=20, edgecolor='black', alpha=0.7, color='green')
        ax.axvline(df['reduction_factor'].mean(), color='red', linestyle='--',
                  label=f'Mean: {df["reduction_factor"].mean():.1f}x')
        ax.set_xlabel('Computational Reduction Factor')
        ax.set_ylabel('Count')
        ax.set_title('Computational Reduction Distribution')
        ax.set_xscale('log')
        ax.legend()
        
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_folder, 'algorithm_efficiency.png'), 
                   dpi=300, bbox_inches='tight')
        
        # Print efficiency statistics
        print(f"\nAlgorithm Efficiency Statistics:")
        print(f"  Average efficiency ratio: {df['efficiency_ratio'].mean()*100:.2f}%")
        print(f"  Average reduction factor: {df['reduction_factor'].mean():.1f}x")
        print(f"  Average synergy success rate: {df['synergy_success_rate'].mean()*100:.2f}%")
        print(f"  Average pruning ratio: {df['pruning_ratio'].mean()*100:.2f}%")
        print(f"  Average time per ERC: {df['time_per_erc'].mean()*1000:.2f} ms")
    
    def create_comprehensive_report(self, df: pd.DataFrame):
        """Create a comprehensive analysis report."""
        print("\n" + "="*60)
        print("COMPREHENSIVE ANALYSIS REPORT")
        print("="*60)
        
        if len(df) == 0:
            print("No data to analyze!")
            return
        
        # Overall statistics
        print(f"\nDataset Overview:")
        print(f"  Networks analyzed: {len(df)}")
        print(f"  Network size range: {df['n_ercs'].min()}-{df['n_ercs'].max()} ERCs")
        print(f"  Species range: {df['n_species'].min()}-{df['n_species'].max()}")
        print(f"  Reactions range: {df['n_reactions'].min()}-{df['n_reactions'].max()}")
        
        print(f"\nPerformance Summary:")
        print(f"  Average execution time: {df['execution_time'].mean():.3f}s (±{df['execution_time'].std():.3f})")
        print(f"  Median execution time: {df['execution_time'].median():.3f}s")
        print(f"  Total synergies found: {df['fundamental_synergies'].sum()}")
        print(f"  Average synergies per network: {df['fundamental_synergies'].mean():.1f}")
        
        print(f"\nEfficiency Metrics:")
        print(f"  Average efficiency ratio: {df['efficiency_ratio'].mean()*100:.2f}%")
        print(f"  Average reduction factor: {df['reduction_factor'].mean():.1f}x")
        print(f"  Best reduction factor: {df['reduction_factor'].max():.1f}x")
        print(f"  Average operations per ERC: {df['operations_per_erc'].mean():.1f}")
        
        # Network structure insights
        print(f"\nNetwork Structure Insights:")
        connected_networks = df['is_connected'].sum()
        print(f"  Connected networks: {connected_networks}/{len(df)} ({connected_networks/len(df)*100:.1f}%)")
        print(f"  Average clustering coefficient: {df['clustering_coefficient'].mean():.3f}")
        print(f"  Average network density: {df['network_density'].mean():.3f}")
        
        # Hierarchy insights
        print(f"\nHierarchy Structure Insights:")
        print(f"  Average hierarchy levels: {df['n_levels'].mean():.1f}")
        print(f"  Average max chain length: {df['max_chain_length'].mean():.1f}")
        print(f"  Average hierarchy density: {df['hierarchy_density'].mean():.3f}")
        print(f"  Average leaf ratio: {df['leaf_ratio'].mean()*100:.1f}%")
        
        # Create summary visualization
        fig = plt.figure(figsize=(20, 12))
        
        # Network size distribution
        ax1 = plt.subplot(2, 4, 1)
        ax1.hist(df['n_ercs'], bins=15, edgecolor='black', alpha=0.7)
        ax1.set_xlabel('Number of ERCs')
        ax1.set_ylabel('Count')
        ax1.set_title('Network Size Distribution')
        
        # Performance distribution
        ax2 = plt.subplot(2, 4, 2)
        ax2.hist(df['execution_time'], bins=15, edgecolor='black', alpha=0.7, color='green')
        ax2.set_xlabel('Execution Time (seconds)')
        ax2.set_ylabel('Count')
        ax2.set_title('Performance Distribution')
        
        # Synergies found
        ax3 = plt.subplot(2, 4, 3)
        ax3.hist(df['fundamental_synergies'], bins=15, edgecolor='black', alpha=0.7, color='red')
        ax3.set_xlabel('Fundamental Synergies Found')
        ax3.set_ylabel('Count')
        ax3.set_title('Synergy Discovery Distribution')
        
        # Efficiency metrics
        ax4 = plt.subplot(2, 4, 4)
        metrics = ['Efficiency\nRatio', 'Pruning\nRatio', 'Success\nRate']
        values = [
            df['efficiency_ratio'].mean() * 100,
            df['pruning_ratio'].mean() * 100,
            df['synergy_success_rate'].mean() * 100
        ]
        bars = ax4.bar(metrics, values, color=['blue', 'orange', 'green'])
        ax4.set_ylabel('Percentage')
        ax4.set_title('Average Efficiency Metrics')
        
        for bar, value in zip(bars, values):
            height = bar.get_height()
            ax4.text(bar.get_x() + bar.get_width()/2., height + 1,
                    f'{value:.1f}%', ha='center', va='bottom')
        
        # Performance vs size
        ax5 = plt.subplot(2, 4, 5)
        ax5.scatter(df['n_ercs'], df['execution_time'], alpha=0.6)
        ax5.set_xlabel('Number of ERCs')
        ax5.set_ylabel('Execution Time (seconds)')
        ax5.set_title('Performance Scaling')
        ax5.set_yscale('log')
        
        # Hierarchy depth impact
        ax6 = plt.subplot(2, 4, 6)
        ax6.scatter(df['max_chain_length'], df['reduction_factor'], alpha=0.6, color='purple')
        ax6.set_xlabel('Maximum Chain Length')
        ax6.set_ylabel('Reduction Factor')
        ax6.set_title('Hierarchy Impact on Efficiency')
        ax6.set_yscale('log')
        
        # Network complexity
        ax7 = plt.subplot(2, 4, 7)
        ax7.scatter(df['network_density'], df['efficiency_ratio'], alpha=0.6, color='brown')
        ax7.set_xlabel('Network Density')
        ax7.set_ylabel('Efficiency Ratio')
        ax7.set_title('Network Complexity vs Efficiency')
        
        # Algorithm stages
        ax8 = plt.subplot(2, 4, 8)
        stage_data = [
            df['partial_overlap_checks'].sum(),
            df['base_pairs_generated'].sum(),
            df['synergies_found'].sum(),
            df['total_pruning'].sum()
        ]
        stage_labels = ['Overlap\nChecks', 'Base\nPairs', 'Synergies\nFound', 'Pruning\nOps']
        colors = ['blue', 'green', 'red', 'orange']
        
        wedges, texts, autotexts = ax8.pie(stage_data, labels=stage_labels, colors=colors,
                                          autopct='%1.1f%%', startangle=90)
        ax8.set_title('Algorithm Operation Distribution')
        
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_folder, 'comprehensive_report.png'), 
                   dpi=300, bbox_inches='tight')
        
        # Save detailed report
        report_file = os.path.join(self.output_folder, 'analysis_report.txt')
        with open(report_file, 'w') as f:
            f.write("DYNAMICAL HIERARCHY ALGORITHM ANALYSIS REPORT\n")
            f.write("="*60 + "\n\n")
            
            f.write(f"Dataset Overview:\n")
            f.write(f"  Networks analyzed: {len(df)}\n")
            f.write(f"  Network size range: {df['n_ercs'].min()}-{df['n_ercs'].max()} ERCs\n")
            f.write(f"  Average network size: {df['n_ercs'].mean():.1f} ERCs\n\n")
            
            f.write(f"Performance Summary:\n")
            f.write(f"  Average execution time: {df['execution_time'].mean():.3f}s\n")
            f.write(f"  Average efficiency ratio: {df['efficiency_ratio'].mean()*100:.2f}%\n")
            f.write(f"  Average reduction factor: {df['reduction_factor'].mean():.1f}x\n")
            f.write(f"  Total synergies found: {df['fundamental_synergies'].sum()}\n\n")
            
            # Best and worst performers
            best_idx = df['reduction_factor'].idxmax()
            worst_idx = df['reduction_factor'].idxmin()
            
            f.write(f"Best Performer:\n")
            f.write(f"  File: {df.loc[best_idx, 'file']}\n")
            f.write(f"  ERCs: {df.loc[best_idx, 'n_ercs']}\n")
            f.write(f"  Reduction factor: {df.loc[best_idx, 'reduction_factor']:.1f}x\n")
            f.write(f"  Execution time: {df.loc[best_idx, 'execution_time']:.3f}s\n\n")
            
            f.write(f"Worst Performer:\n")
            f.write(f"  File: {df.loc[worst_idx, 'file']}\n")
            f.write(f"  ERCs: {df.loc[worst_idx, 'n_ercs']}\n")
            f.write(f"  Reduction factor: {df.loc[worst_idx, 'reduction_factor']:.1f}x\n")
            f.write(f"  Execution time: {df.loc[worst_idx, 'execution_time']:.3f}s\n\n")
        
        print(f"\nDetailed report saved to: {report_file}")
    
    def run_complete_analysis(self, file_pattern: str = "*.txt", 
                             max_files: Optional[int] = None,
                             timeout: float = 300):
        """Run complete DH algorithm analysis."""
        # Run network analyses
        df = self.run_analysis(file_pattern, max_files, timeout)
        
        if len(df) == 0:
            print("No successful analyses to process!")
            return None
        
        # Perform detailed analyses
        scaling_results = self.analyze_performance_scaling(df)
        hierarchy_correlations = self.analyze_hierarchy_impact(df)
        self.analyze_algorithm_efficiency(df)
        self.create_comprehensive_report(df)
        
        # Show all plots
        plt.show()
        
        return {
            'dataframe': df,
            'scaling_results': scaling_results,
            'hierarchy_correlations': hierarchy_correlations,
            'network_data': self.network_data
        }


def analyze_dh_algorithm(network_folder: str, 
                        output_folder: str = "dh_analysis_results",
                        file_pattern: str = "*.txt",
                        max_files: Optional[int] = None,
                        timeout: float = 300):
    """
    Main function to analyze DH algorithm performance across multiple networks.
    
    Args:
        network_folder: Path to folder containing network files
        output_folder: Path to save results
        file_pattern: Pattern to match network files
        max_files: Maximum number of files to process
        timeout: Maximum time per analysis (seconds)
    
    Returns:
        Tuple of (analyzer, results)
    """
    analyzer = DHNetworkAnalyzer(network_folder, output_folder)
    results = analyzer.run_complete_analysis(file_pattern, max_files, timeout)
    return analyzer, results


if __name__ == "__main__":
    # Example usage
    import sys
    
    if len(sys.argv) > 1:
        network_folder = sys.argv[1]
    else:
        # Default folder - change this to your network folder
        network_folder = "networks/biomodels_all_txt"
    
    print(f"Analyzing DH algorithm performance on networks in: {network_folder}")
    
    # Run analysis
    analyzer, results = analyze_dh_algorithm(
        network_folder,
        output_folder="networks/biomodels_all_txt/dh_analysis_results",
        file_pattern="*.txt",
        max_files=None,  # Process all files
        timeout=600  # 10 minute timeout per analysis
    )
    
    if results:
        print(f"\nAnalysis complete! Results saved to: {analyzer.output_folder}")
        print(f"Total networks analyzed: {len(results['dataframe'])}")
        print(f"Check the output folder for detailed results and visualizations.")
    
    # Keep plots open
    plt.show(block=True)