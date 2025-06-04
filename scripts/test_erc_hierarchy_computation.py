#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ERC Hierarchy Analysis Script

This script builds ERC hierarchies from reaction networks and performs
comprehensive structural analysis including levels, chains, degrees,
and other topological measures.

@author: Analysis Script
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
from collections import defaultdict, Counter
import seaborn as sns
from pathlib import Path

# Import required modules (assuming they're in the path)
from pyCOT.io.functions import read_txt
from pyCOT.ERC_Hierarchy import ERC, species_list_to_names

class ERCHierarchyAnalyzer:
    """Class for analyzing ERC hierarchy structures."""
    
    def __init__(self, reaction_network_file):
        """Initialize analyzer with reaction network file."""
        self.rn_file = reaction_network_file
        self.rn = None
        self.ercs = None
        self.hierarchy_graph = None
        self.analysis_results = {}
        
    def load_reaction_network(self):
        """Load reaction network from file."""
        print(f"Loading reaction network from: {self.rn_file}")
        self.rn = read_txt(self.rn_file)
        print(f"Loaded {len(self.rn.species())} species and {len(self.rn.reactions())} reactions")
        
    def generate_ercs(self):
        """Generate ERCs from the reaction network."""
        print("Generating ERCs...")
        self.ercs = ERC.ERCs(self.rn)
        print(f"Generated {len(self.ercs)} ERCs")
        
        # Print ERC details
        for i, erc in enumerate(self.ercs):
            closure = species_list_to_names(erc.get_closure(self.rn))
            min_gens = [species_list_to_names(gen) for gen in erc.min_generators]
            print(f"  {erc.label}: Closure={closure}, Min_generators={min_gens}")
    
    def build_hierarchy(self):
        """Build the ERC hierarchy graph."""
        print("Building hierarchy graph...")
        self.hierarchy_graph = ERC.build_hierarchy_graph(self.ercs, self.rn)
        print(f"Hierarchy has {self.hierarchy_graph.number_of_nodes()} nodes and {self.hierarchy_graph.number_of_edges()} edges")
        
    def plot_hierarchy(self, figsize=(12, 8), save_path=None):
        """
        Plot the ERC hierarchy.
        
        Parameters
        ----------
        figsize : tuple
            Figure size (width, height) in inches
        save_path : str, optional
            Path where to save the plot
            
        Returns
        -------
        matplotlib.figure.Figure
            The generated figure
        """
        print("Plotting hierarchy...")
        fig = ERC.plot_hierarchy(self.ercs, self.rn, self.hierarchy_graph, figsize=figsize)
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Hierarchy plot saved to: {save_path}")
            
        return fig
    
    def analyze_hierarchy_structure(self):
        """Perform comprehensive structural analysis of the hierarchy."""
        print("Analyzing hierarchy structure...")
        
        G = self.hierarchy_graph
        levels = ERC.get_node_levels(G)
        
        # Basic metrics
        self.analysis_results['basic_metrics'] = {
            'total_nodes': G.number_of_nodes(),
            'total_edges': G.number_of_edges(),
            'max_level': max(levels.values()) if levels else 0,
            'min_level': min(levels.values()) if levels else 0,
            'hierarchy_depth': max(levels.values()) - min(levels.values()) + 1 if levels else 0
        }
        
        # Level distribution
        level_distribution = Counter(levels.values())
        self.analysis_results['level_distribution'] = dict(level_distribution)
        
        # Nodes per level
        nodes_per_level = defaultdict(list)
        for node, level in levels.items():
            nodes_per_level[level].append(node)
        self.analysis_results['nodes_per_level'] = dict(nodes_per_level)
        
        # Degree analysis
        in_degrees = dict(G.in_degree())
        out_degrees = dict(G.out_degree())
        
        # Degree statistics per level
        degree_stats_per_level = {}
        for level in sorted(set(levels.values())):
            level_nodes = nodes_per_level[level]
            level_in_degrees = [in_degrees[node] for node in level_nodes]
            level_out_degrees = [out_degrees[node] for node in level_nodes]
            
            degree_stats_per_level[level] = {
                'avg_in_degree': np.mean(level_in_degrees) if level_in_degrees else 0,
                'avg_out_degree': np.mean(level_out_degrees) if level_out_degrees else 0,
                'max_in_degree': max(level_in_degrees) if level_in_degrees else 0,
                'max_out_degree': max(level_out_degrees) if level_out_degrees else 0,
                'std_in_degree': np.std(level_in_degrees) if level_in_degrees else 0,
                'std_out_degree': np.std(level_out_degrees) if level_out_degrees else 0,
                'node_count': len(level_nodes)
            }
        
        self.analysis_results['degree_stats_per_level'] = degree_stats_per_level
        
        # Root and leaf analysis
        root_nodes = [node for node in G.nodes() if G.in_degree(node) == 0]
        leaf_nodes = [node for node in G.nodes() if G.out_degree(node) == 0]
        
        self.analysis_results['root_leaf_analysis'] = {
            'root_nodes': root_nodes,
            'leaf_nodes': leaf_nodes,
            'num_roots': len(root_nodes),
            'num_leaves': len(leaf_nodes),
            'root_levels': [levels[node] for node in root_nodes],
            'leaf_levels': [levels[node] for node in leaf_nodes]
        }
        
        # Chain analysis
        all_paths = []
        for root in root_nodes:
            for leaf in leaf_nodes:
                paths = list(nx.all_simple_paths(G, root, leaf))
                all_paths.extend(paths)
        
        chain_lengths = [len(path) - 1 for path in all_paths]  # -1 because path length is edges, not nodes
        
        self.analysis_results['chain_analysis'] = {
            'total_chains': len(all_paths),
            'chain_lengths': chain_lengths,
            'avg_chain_length': np.mean(chain_lengths) if chain_lengths else 0,
            'max_chain_length': max(chain_lengths) if chain_lengths else 0,
            'min_chain_length': min(chain_lengths) if chain_lengths else 0,
            'std_chain_length': np.std(chain_lengths) if chain_lengths else 0
        }
        
        # Branching analysis
        branching_factors = []
        internal_nodes = [node for node in G.nodes() if G.out_degree(node) > 0]
        
        for node in internal_nodes:
            branching_factors.append(G.out_degree(node))
        
        self.analysis_results['branching_analysis'] = {
            'internal_nodes': len(internal_nodes),
            'avg_branching_factor': np.mean(branching_factors) if branching_factors else 0,
            'max_branching_factor': max(branching_factors) if branching_factors else 0,
            'std_branching_factor': np.std(branching_factors) if branching_factors else 0
        }
        
        # Width analysis (max nodes at any level)
        level_widths = [len(nodes_per_level[level]) for level in sorted(set(levels.values()))]
        
        self.analysis_results['width_analysis'] = {
            'level_widths': level_widths,
            'max_width': max(level_widths) if level_widths else 0,
            'avg_width': np.mean(level_widths) if level_widths else 0,
            'width_at_each_level': {level: len(nodes_per_level[level]) for level in sorted(set(levels.values()))}
        }
        
        # Connectivity metrics
        if G.number_of_nodes() > 0:
            density = nx.density(G)
            # For directed graphs, we can check weak connectivity
            is_weakly_connected = nx.is_weakly_connected(G)
            num_weakly_connected_components = nx.number_weakly_connected_components(G)
        else:
            density = 0
            is_weakly_connected = False
            num_weakly_connected_components = 0
            
        self.analysis_results['connectivity_metrics'] = {
            'density': density,
            'is_weakly_connected': is_weakly_connected,
            'num_weakly_connected_components': num_weakly_connected_components
        }
        
        # Topological properties
        # Transitivity (clustering coefficient for directed graphs)
        try:
            transitivity = nx.transitivity(G)
        except:
            transitivity = 0
            
        self.analysis_results['topological_properties'] = {
            'transitivity': transitivity,
            'reciprocity': nx.reciprocity(G) if G.number_of_edges() > 0 else 0
        }
    
    def print_analysis_summary(self):
        """Print a comprehensive summary of the hierarchy analysis."""
        print("\n" + "="*70)
        print("ERC HIERARCHY ANALYSIS SUMMARY")
        print("="*70)
        
        # Basic metrics
        basic = self.analysis_results['basic_metrics']
        print(f"\nBASIC METRICS:")
        print(f"  Total ERCs: {basic['total_nodes']}")
        print(f"  Total containment relations: {basic['total_edges']}")
        print(f"  Hierarchy depth: {basic['hierarchy_depth']} levels")
        print(f"  Level range: {basic['min_level']} to {basic['max_level']}")
        
        # Level distribution
        print(f"\nLEVEL DISTRIBUTION:")
        for level in sorted(self.analysis_results['level_distribution'].keys()):
            count = self.analysis_results['level_distribution'][level]
            print(f"  Level {level}: {count} ERCs")
        
        # Root and leaf analysis
        root_leaf = self.analysis_results['root_leaf_analysis']
        print(f"\nROOT AND LEAF ANALYSIS:")
        print(f"  Root ERCs (most general): {root_leaf['num_roots']} ({root_leaf['root_nodes']})")
        print(f"  Leaf ERCs (most specific): {root_leaf['num_leaves']} ({root_leaf['leaf_nodes']})")
        
        # Chain analysis
        chain = self.analysis_results['chain_analysis']
        print(f"\nCHAIN ANALYSIS:")
        print(f"  Total containment chains: {chain['total_chains']}")
        print(f"  Average chain length: {chain['avg_chain_length']:.2f}")
        print(f"  Chain length range: {chain['min_chain_length']} to {chain['max_chain_length']}")
        
        # Degree analysis summary
        print(f"\nDEGREE ANALYSIS BY LEVEL:")
        degree_stats = self.analysis_results['degree_stats_per_level']
        for level in sorted(degree_stats.keys()):
            stats = degree_stats[level]
            print(f"  Level {level}: {stats['node_count']} ERCs, "
                  f"avg in-degree: {stats['avg_in_degree']:.2f}, "
                  f"avg out-degree: {stats['avg_out_degree']:.2f}")
        
        # Branching analysis
        branch = self.analysis_results['branching_analysis']
        print(f"\nBRANCHING ANALYSIS:")
        print(f"  Internal nodes: {branch['internal_nodes']}")
        print(f"  Average branching factor: {branch['avg_branching_factor']:.2f}")
        print(f"  Maximum branching factor: {branch['max_branching_factor']}")
        
        # Width analysis
        width = self.analysis_results['width_analysis']
        print(f"\nWIDTH ANALYSIS:")
        print(f"  Maximum level width: {width['max_width']} ERCs")
        print(f"  Average level width: {width['avg_width']:.2f}")
        
        # Connectivity
        connectivity = self.analysis_results['connectivity_metrics']
        print(f"\nCONNECTIVITY METRICS:")
        print(f"  Graph density: {connectivity['density']:.3f}")
        print(f"  Weakly connected: {connectivity['is_weakly_connected']}")
        print(f"  Connected components: {connectivity['num_weakly_connected_components']}")
        
        print("\n" + "="*70)
    
    def create_analysis_plots(self):
        """Create comprehensive analysis plots."""
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        fig.suptitle('ERC Hierarchy Analysis', fontsize=16, fontweight='bold')
        
        # 1. Level distribution
        levels = list(self.analysis_results['level_distribution'].keys())
        counts = list(self.analysis_results['level_distribution'].values())
        axes[0, 0].bar(levels, counts, alpha=0.7, color='skyblue')
        axes[0, 0].set_xlabel('Hierarchy Level')
        axes[0, 0].set_ylabel('Number of ERCs')
        axes[0, 0].set_title('ERCs per Level')
        axes[0, 0].grid(True, alpha=0.3)
        
        # 2. Chain length distribution
        chain_lengths = self.analysis_results['chain_analysis']['chain_lengths']
        if chain_lengths:
            axes[0, 1].hist(chain_lengths, bins=max(1, len(set(chain_lengths))), 
                           alpha=0.7, color='lightgreen', edgecolor='black')
        axes[0, 1].set_xlabel('Chain Length')
        axes[0, 1].set_ylabel('Frequency')
        axes[0, 1].set_title('Chain Length Distribution')
        axes[0, 1].grid(True, alpha=0.3)
        
        # 3. In-degree vs Out-degree scatter
        in_degrees = [self.hierarchy_graph.in_degree(node) for node in self.hierarchy_graph.nodes()]
        out_degrees = [self.hierarchy_graph.out_degree(node) for node in self.hierarchy_graph.nodes()]
        axes[0, 2].scatter(in_degrees, out_degrees, alpha=0.7, color='orange')
        axes[0, 2].set_xlabel('In-degree')
        axes[0, 2].set_ylabel('Out-degree')
        axes[0, 2].set_title('In-degree vs Out-degree')
        axes[0, 2].grid(True, alpha=0.3)
        
        # 4. Average degrees per level
        degree_stats = self.analysis_results['degree_stats_per_level']
        levels = sorted(degree_stats.keys())
        avg_in_degrees = [degree_stats[level]['avg_in_degree'] for level in levels]
        avg_out_degrees = [degree_stats[level]['avg_out_degree'] for level in levels]
        
        axes[1, 0].plot(levels, avg_in_degrees, 'o-', label='In-degree', color='red')
        axes[1, 0].plot(levels, avg_out_degrees, 's-', label='Out-degree', color='blue')
        axes[1, 0].set_xlabel('Hierarchy Level')
        axes[1, 0].set_ylabel('Average Degree')
        axes[1, 0].set_title('Average Degrees by Level')
        axes[1, 0].legend()
        axes[1, 0].grid(True, alpha=0.3)
        
        # 5. Level width
        width_analysis = self.analysis_results['width_analysis']
        level_widths = [width_analysis['width_at_each_level'][level] for level in sorted(width_analysis['width_at_each_level'].keys())]
        levels = sorted(width_analysis['width_at_each_level'].keys())
        
        axes[1, 1].plot(levels, level_widths, 'o-', color='purple', linewidth=2)
        axes[1, 1].set_xlabel('Hierarchy Level')
        axes[1, 1].set_ylabel('Number of ERCs')
        axes[1, 1].set_title('Level Width')
        axes[1, 1].grid(True, alpha=0.3)
        
        # 6. Degree distribution
        all_degrees = in_degrees + out_degrees
        if all_degrees:
            axes[1, 2].hist(all_degrees, bins=max(1, len(set(all_degrees))), 
                           alpha=0.7, color='pink', edgecolor='black')
        axes[1, 2].set_xlabel('Degree')
        axes[1, 2].set_ylabel('Frequency')
        axes[1, 2].set_title('Overall Degree Distribution')
        axes[1, 2].grid(True, alpha=0.3)
        
        plt.tight_layout()
        return fig
    
    def export_analysis_to_csv(self, save_dir=None):
        """Export detailed analysis results to CSV files."""
        if save_dir:
            save_dir = Path(save_dir)
            save_dir.mkdir(exist_ok=True)
        else:
            save_dir = Path(".")
        
        # ERC details
        erc_data = []
        levels = ERC.get_node_levels(self.hierarchy_graph)
        
        for erc in self.ercs:
            closure = species_list_to_names(erc.get_closure(self.rn))
            min_gens = [species_list_to_names(gen) for gen in erc.min_generators]
            
            erc_data.append({
                'ERC_Label': erc.label,
                'Level': levels.get(erc.label, 0),
                'Closure_Size': len(closure),
                'Closure_Species': str(closure),
                'Min_Generators_Count': len(min_gens),
                'Min_Generators': str(min_gens),
                'In_Degree': self.hierarchy_graph.in_degree(erc.label),
                'Out_Degree': self.hierarchy_graph.out_degree(erc.label)
            })
        
        erc_df = pd.DataFrame(erc_data)
        erc_csv_path = save_dir / "erc_details.csv"
        erc_df.to_csv(erc_csv_path, index=False)
        print(f"ERC details exported to: {erc_csv_path}")
        
        # Level statistics
        level_stats_data = []
        for level, stats in self.analysis_results['degree_stats_per_level'].items():
            level_stats_data.append({
                'Level': level,
                'Node_Count': stats['node_count'],
                'Avg_In_Degree': stats['avg_in_degree'],
                'Avg_Out_Degree': stats['avg_out_degree'],
                'Max_In_Degree': stats['max_in_degree'],
                'Max_Out_Degree': stats['max_out_degree'],
                'Std_In_Degree': stats['std_in_degree'],
                'Std_Out_Degree': stats['std_out_degree']
            })
        
        level_stats_df = pd.DataFrame(level_stats_data)
        level_csv_path = save_dir / "level_statistics.csv"
        level_stats_df.to_csv(level_csv_path, index=False)
        print(f"Level statistics exported to: {level_csv_path}")
        
        return erc_df, level_stats_df
    
    def run_full_analysis(self, plot=True, save_dir=None):
        """Run the complete ERC hierarchy analysis pipeline."""
        print("Starting ERC Hierarchy Analysis...")
        print("="*50)
        
        # Load and process data
        self.load_reaction_network()
        self.generate_ercs()
        self.build_hierarchy()
        
        # Perform analysis
        self.analyze_hierarchy_structure()
        
        # Generate outputs
        if plot:
            if save_dir:
                Path(save_dir).mkdir(exist_ok=True)
                hierarchy_plot_path = str(Path(save_dir) / "hierarchy_plot.png")
            else:
                hierarchy_plot_path = None
            self.print_analysis_summary()    
            self.plot_hierarchy(save_path=hierarchy_plot_path)
            analysis_fig = self.create_analysis_plots()
            
            if save_dir:
                analysis_plot_path = str(Path(save_dir) / "analysis_plots.png")
                analysis_fig.savefig(analysis_plot_path, dpi=300, bbox_inches='tight')
                print(f"Analysis plots saved to: {analysis_plot_path}")
        
        
        self.export_analysis_to_csv(save_dir=save_dir)
        
        print(f"\nAnalysis complete!")
        return self.analysis_results


def main():
    """Main function to run the ERC hierarchy analysis."""
    
    # Configuration
    file_path = 'networks/testing/Farm.txt'
    file_path = 'networks/RandomAlife/RN_Ns_20_Norg_6_id_137.txt'
    #file_path = 'networks/Navarino/RN_IN_05.txt'
    #file_path = 'networks/biomodels_interesting/central_ecoli.txt'
    #file_path = 'networks/biomodels_interesting/BIOMD0000000237_manyOrgs.txt'    

    SAVE_DIR = "hierarchy_analysis_results"  # Directory to save results
    
    # Create analyzer and run analysis
    analyzer = ERCHierarchyAnalyzer(file_path)
    results = analyzer.run_full_analysis(plot=True, save_dir=SAVE_DIR)
    
    return analyzer, results


if __name__ == "__main__":
    analyzer, results = main()