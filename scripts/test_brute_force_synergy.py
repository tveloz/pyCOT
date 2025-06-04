#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Comprehensive Synergy Test Script

This script provides detailed testing and analysis of the brute force synergy
calculator, showing exactly why each synergy is kept or discarded.

@author: Synergy Test Analysis
"""

import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from collections import defaultdict
from typing import List, Dict, Set, Tuple, Optional
import os

# Import your modules
from pyCOT.io.functions import read_txt
from pyCOT.ERC_Hierarchy import ERC, species_list_to_names
from pyCOT.Synergy_Brute_Force_Algorithm import BruteForceSynergyCalculator


class SynergyTestAnalyzer:
    """Comprehensive analyzer for synergy test results."""
    
    def __init__(self, reaction_network_file: str):
        """Initialize with reaction network file."""
        self.rn_file = reaction_network_file
        self.rn = None
        self.ercs = None
        self.hierarchy_graph = None
        self.calculator = None
        self.test_results = {}
        
    def setup(self):
        """Load reaction network and setup calculator."""
        print("="*80)
        print("COMPREHENSIVE SYNERGY TEST ANALYSIS")
        print("="*80)
        print(f"Loading reaction network from: {self.rn_file}")
        
        # Load reaction network
        self.rn = read_txt(self.rn_file)
        
        # Build ERCs and hierarchy
        print("Building ERC hierarchy...")
        self.ercs = ERC.ERCs(self.rn)
        self.hierarchy_graph = ERC.build_hierarchy_graph(self.ercs, self.rn)
        
        # Initialize calculator
        self.calculator = BruteForceSynergyCalculator(self.rn, self.ercs, self.hierarchy_graph)
        
        print(f"Setup complete:")
        print(f"  - {len(self.ercs)} ERCs found")
        print(f"  - {self.hierarchy_graph.number_of_edges()} containment relations")
        print(f"  - {len(self.rn.species())} species")
        print(f"  - {len(self.rn.reactions())} reactions")
        
        # Print ERC details
        print(f"\nERC Details:")
        for erc in self.ercs:
            closure_species = species_list_to_names(erc.get_closure(self.rn))
            print(f"  {erc.label}: closure = {closure_species}")
        
    def analyze_hierarchy_structure(self):
        """Analyze and display hierarchy structure."""
        print(f"\n" + "="*80)
        print("HIERARCHY STRUCTURE ANALYSIS")
        print("="*80)
        
        # Get hierarchy levels
        levels = ERC.get_node_levels(self.hierarchy_graph)
        
        print(f"Hierarchy levels:")
        level_groups = defaultdict(list)
        for erc_label, level in levels.items():
            level_groups[level].append(erc_label)
        
        for level in sorted(level_groups.keys()):
            ercs_at_level = sorted(level_groups[level])
            print(f"  Level {level}: {ercs_at_level}")
        
        # Show containment relationships
        print(f"\nContainment relationships:")
        for parent, child in self.hierarchy_graph.edges():
            print(f"  {parent} contains {child}")
            
        return levels, level_groups
    
    def test_individual_synergies(self, verbose: bool = True):
        """Test each potential synergy combination individually with detailed output."""
        print(f"\n" + "="*80)
        print("INDIVIDUAL SYNERGY TESTING")
        print("="*80)
        
        all_combinations = []
        
        # Test all possible combinations
        for target in self.ercs:
            other_ercs = [erc for erc in self.ercs if erc != target]
            
            for i, base1 in enumerate(other_ercs):
                for j, base2 in enumerate(other_ercs[i+1:], i+1):
                    # Order bases consistently
                    if base1.label < base2.label:
                        ordered_base1, ordered_base2 = base1, base2
                    else:
                        ordered_base1, ordered_base2 = base2, base1
                    
                    combination = (ordered_base1, ordered_base2, target)
                    if combination not in all_combinations:
                        all_combinations.append(combination)
        
        print(f"Testing {len(all_combinations)} unique ERC combinations...")
        
        test_results = []
        synergies_found = []
        
        for i, (base1, base2, target) in enumerate(all_combinations):
            print(f"\n--- Test {i+1}/{len(all_combinations)} ---")
            print(f"Testing: {base1.label} + {base2.label} → {target.label}")
            
            # Get closures for analysis
            base1_closure = set(species_list_to_names(base1.get_closure(self.rn)))
            base2_closure = set(species_list_to_names(base2.get_closure(self.rn)))
            target_closure = set(species_list_to_names(target.get_closure(self.rn)))
            combined_closure = self.calculator.compute_combined_closure(base1, base2)
            
            print(f"  Base1 ({base1.label}) closure: {sorted(base1_closure)}")
            print(f"  Base2 ({base2.label}) closure: {sorted(base2_closure)}")
            print(f"  Target ({target.label}) closure: {sorted(target_closure)}")
            print(f"  Combined closure: {sorted(combined_closure)}")
            
            # Test synergy
            is_synergy, details = self.calculator.validates_synergy(base1, base2, target, verbose=False)
            
            if is_synergy:
                print(f"  ✓ SYNERGY DETECTED!")
                print(f"    Coverage ratio: {details['coverage_ratio']:.1%}")
                print(f"    Synergy ratio: {details['synergy_ratio']:.1%}")
                print(f"    Synergistic generators: {len(details['synergistic_generators'])}")
                
                synergy_dict = {
                    'base1': base1,
                    'base2': base2,
                    'target': target,
                    'details': details
                }
                synergies_found.append(synergy_dict)
            else:
                print(f"  ✗ No synergy")
                
                # Analyze why no synergy
                reasons = []
                if not combined_closure.issuperset(target_closure):
                    reasons.append("combined closure doesn't contain target closure")
                if base1_closure.issuperset(target_closure):
                    reasons.append("base1 alone contains target closure")
                if base2_closure.issuperset(target_closure):
                    reasons.append("base2 alone contains target closure")
                    
                # Check for synergistic generators
                has_synergistic = False
                for gen in target.min_generators:
                    gen_species = set(species_list_to_names(gen))
                    if (gen_species.issubset(combined_closure) and
                        not gen_species.issubset(base1_closure) and
                        not gen_species.issubset(base2_closure) and
                        len(gen_species & base1_closure - base2_closure) > 0 and
                        len(gen_species & base2_closure - base1_closure) > 0):
                        has_synergistic = True
                        break
                
                if not has_synergistic:
                    reasons.append("no generators require true synergy")
                
                print(f"    Reasons: {', '.join(reasons)}")
            
            test_result = {
                'base1': base1,
                'base2': base2,
                'target': target,
                'is_synergy': is_synergy,
                'details': details if is_synergy else None,
                'closures': {
                    'base1': base1_closure,
                    'base2': base2_closure,
                    'target': target_closure,
                    'combined': combined_closure
                }
            }
            test_results.append(test_result)
        
        self.test_results['individual_tests'] = test_results
        self.test_results['synergies_found'] = synergies_found
        
        print(f"\n" + "-"*50)
        print(f"INDIVIDUAL TESTING SUMMARY")
        print(f"-"*50)
        print(f"Total combinations tested: {len(all_combinations)}")
        print(f"Synergies found: {len(synergies_found)}")
        print(f"Success rate: {len(synergies_found)/len(all_combinations)*100:.1f}%")
        
        return synergies_found
    
    def test_minimality_analysis(self, synergies: List[Dict], verbose: bool = True):
        """Perform detailed minimality analysis for each synergy."""
        print(f"\n" + "="*80)
        print("MINIMALITY ANALYSIS")
        print("="*80)
        
        if not synergies:
            print("No synergies to analyze for minimality.")
            return []
        
        minimal_synergies = []
        minimality_details = []
        
        for i, synergy in enumerate(synergies):
            print(f"\n--- Minimality Test {i+1}/{len(synergies)} ---")
            base1 = synergy['base1']
            base2 = synergy['base2']
            target = synergy['target']
            
            print(f"Analyzing: {base1.label} + {base2.label} → {target.label}")
            
            # Test all three minimality criteria
            criterion1, base1_violations = self.calculator.check_minimal_base1_criterion(synergy, verbose=True)
            criterion2, base2_violations = self.calculator.check_minimal_base2_criterion(synergy, verbose=True)
            criterion3, target_violations = self.calculator.check_maximal_target_criterion(synergy, verbose=True)
            
            is_minimal = criterion1 and criterion2 and criterion3
            
            minimality_info = {
                'base1_minimal': criterion1,
                'base2_minimal': criterion2,
                'target_maximal': criterion3,
                'is_minimal': is_minimal,
                'base1_violations': base1_violations,
                'base2_violations': base2_violations,
                'target_violations': target_violations
            }
            
            # Add minimality info to synergy
            synergy_with_minimality = synergy.copy()
            synergy_with_minimality['minimality'] = minimality_info
            
            print(f"\n  MINIMALITY SUMMARY:")
            print(f"    Base1 ({base1.label}) minimal: {'✓' if criterion1 else '✗'}")
            if base1_violations:
                print(f"      Violations: {base1_violations}")
            print(f"    Base2 ({base2.label}) minimal: {'✓' if criterion2 else '✗'}")
            if base2_violations:
                print(f"      Violations: {base2_violations}")
            print(f"    Target ({target.label}) maximal: {'✓' if criterion3 else '✗'}")
            if target_violations:
                print(f"      Violations: {target_violations}")
            print(f"    OVERALL: {'✓ MINIMAL' if is_minimal else '✗ NOT MINIMAL'}")
            
            minimality_details.append(synergy_with_minimality)
            
            if is_minimal:
                minimal_synergies.append(synergy_with_minimality)
        
        self.test_results['minimality_analysis'] = minimality_details
        self.test_results['minimal_synergies'] = minimal_synergies
        
        print(f"\n" + "-"*50)
        print(f"MINIMALITY ANALYSIS SUMMARY")
        print(f"-"*50)
        print(f"Total synergies analyzed: {len(synergies)}")
        print(f"Minimal synergies: {len(minimal_synergies)}")
        print(f"Minimality rate: {len(minimal_synergies)/len(synergies)*100:.1f}%")
        
        # Summary of why synergies were filtered
        base1_failures = sum(1 for s in minimality_details if not s['minimality']['base1_minimal'])
        base2_failures = sum(1 for s in minimality_details if not s['minimality']['base2_minimal'])
        target_failures = sum(1 for s in minimality_details if not s['minimality']['target_maximal'])
        
        print(f"\nFailure reasons:")
        print(f"  Base1 not minimal: {base1_failures}")
        print(f"  Base2 not minimal: {base2_failures}")
        print(f"  Target not maximal: {target_failures}")
        
        return minimal_synergies
    
    def compare_algorithms(self):
        """Compare manual analysis with calculator results."""
        print(f"\n" + "="*80)
        print("ALGORITHM COMPARISON")
        print("="*80)
        
        # Run calculator
        print("Running brute force calculator...")
        all_synergies_calc, minimal_synergies_calc = self.calculator.brute_force(minimal=True, verbose=False)
        
        # Compare with manual analysis
        manual_synergies = self.test_results.get('synergies_found', [])
        manual_minimal = self.test_results.get('minimal_synergies', [])
        
        print(f"\nComparison Results:")
        print(f"  Manual analysis - All synergies: {len(manual_synergies)}")
        print(f"  Calculator - All synergies: {len(all_synergies_calc)}")
        print(f"  Match: {'✓' if len(manual_synergies) == len(all_synergies_calc) else '✗'}")
        
        print(f"  Manual analysis - Minimal synergies: {len(manual_minimal)}")
        print(f"  Calculator - Minimal synergies: {len(minimal_synergies_calc)}")
        print(f"  Match: {'✓' if len(manual_minimal) == len(minimal_synergies_calc) else '✗'}")
        
        # Detailed comparison
        if len(manual_synergies) != len(all_synergies_calc):
            print(f"\n  DISCREPANCY DETECTED in all synergies!")
            manual_set = {(s['base1'].label, s['base2'].label, s['target'].label) for s in manual_synergies}
            calc_set = {(s['base1'].label, s['base2'].label, s['target'].label) for s in all_synergies_calc}
            
            only_manual = manual_set - calc_set
            only_calc = calc_set - manual_set
            
            if only_manual:
                print(f"    Only in manual: {only_manual}")
            if only_calc:
                print(f"    Only in calculator: {only_calc}")
        
        if len(manual_minimal) != len(minimal_synergies_calc):
            print(f"\n  DISCREPANCY DETECTED in minimal synergies!")
            manual_min_set = {(s['base1'].label, s['base2'].label, s['target'].label) for s in manual_minimal}
            calc_min_set = {(s['base1'].label, s['base2'].label, s['target'].label) for s in minimal_synergies_calc}
            
            only_manual_min = manual_min_set - calc_min_set
            only_calc_min = calc_min_set - manual_min_set
            
            if only_manual_min:
                print(f"    Only in manual minimal: {only_manual_min}")
            if only_calc_min:
                print(f"    Only in calculator minimal: {only_calc_min}")
        
        return all_synergies_calc, minimal_synergies_calc
    
    def create_comprehensive_visualization(self, figsize=(16, 12)):
        """Create comprehensive visualization of results."""
        if not self.test_results:
            print("No test results to visualize.")
            return None
        
        fig, axes = plt.subplots(2, 2, figsize=figsize)
        fig.suptitle('Comprehensive Synergy Analysis Results', fontsize=16, fontweight='bold')
        
        # 1. Hierarchy with synergies
        ax1 = axes[0, 0]
        self._plot_hierarchy_with_synergies(ax1)
        
        # 2. Synergy statistics
        ax2 = axes[0, 1]
        self._plot_synergy_statistics(ax2)
        
        # 3. Minimality analysis
        ax3 = axes[1, 0]
        self._plot_minimality_analysis(ax3)
        
        # 4. Detailed breakdown
        ax4 = axes[1, 1]
        self._plot_detailed_breakdown(ax4)
        
        plt.tight_layout()
        return fig
    
    def _plot_hierarchy_with_synergies(self, ax):
        """Plot hierarchy structure with synergies highlighted."""
        # Use existing hierarchy layout
        levels = ERC.get_node_levels(self.hierarchy_graph)
        pos = {}
        level_nodes = defaultdict(list)
        
        for node, level in levels.items():
            level_nodes[level].append(node)
        
        for level, nodes in level_nodes.items():
            n_nodes = len(nodes)
            nodes.sort()
            for i, node in enumerate(nodes):
                x = (i - (n_nodes-1)/2) * 2.0
                y = level * 1.5
                pos[node] = (x, y)
        
        # Draw hierarchy edges
        nx.draw_networkx_edges(self.hierarchy_graph, pos, ax=ax,
                             edge_color='lightgray', arrows=True, arrowsize=10)
        
        # Color nodes based on involvement in synergies
        synergies = self.test_results.get('synergies_found', [])
        involved_nodes = set()
        for s in synergies:
            involved_nodes.update([s['base1'].label, s['base2'].label, s['target'].label])
        
        node_colors = ['lightcoral' if node in involved_nodes else 'lightblue' 
                      for node in self.hierarchy_graph.nodes()]
        
        nx.draw_networkx_nodes(self.hierarchy_graph, pos, ax=ax,
                             node_color=node_colors, node_size=500)
        nx.draw_networkx_labels(self.hierarchy_graph, pos, ax=ax, font_size=8)
        
        ax.set_title('ERC Hierarchy\n(Red: involved in synergies)')
        ax.axis('off')
    
    def _plot_synergy_statistics(self, ax):
        """Plot synergy statistics."""
        if 'synergies_found' not in self.test_results:
            ax.text(0.5, 0.5, 'No synergy data', ha='center', va='center', transform=ax.transAxes)
            ax.set_title('Synergy Statistics')
            return
        
        total_combinations = len(self.test_results.get('individual_tests', []))
        synergies_found = len(self.test_results.get('synergies_found', []))
        minimal_synergies = len(self.test_results.get('minimal_synergies', []))
        
        categories = ['Total\nCombinations', 'Synergies\nFound', 'Minimal\nSynergies']
        values = [total_combinations, synergies_found, minimal_synergies]
        colors = ['lightgray', 'lightcoral', 'lightgreen']
        
        bars = ax.bar(categories, values, color=colors)
        ax.set_title('Synergy Statistics')
        ax.set_ylabel('Count')
        
        # Add value labels on bars
        for bar, value in zip(bars, values):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                   f'{value}', ha='center', va='bottom')
    
    def _plot_minimality_analysis(self, ax):
        """Plot minimality failure analysis."""
        if 'minimality_analysis' not in self.test_results:
            ax.text(0.5, 0.5, 'No minimality data', ha='center', va='center', transform=ax.transAxes)
            ax.set_title('Minimality Analysis')
            return
        
        minimality_data = self.test_results['minimality_analysis']
        
        base1_failures = sum(1 for s in minimality_data if not s['minimality']['base1_minimal'])
        base2_failures = sum(1 for s in minimality_data if not s['minimality']['base2_minimal'])
        target_failures = sum(1 for s in minimality_data if not s['minimality']['target_maximal'])
        
        categories = ['Base1\nNot Minimal', 'Base2\nNot Minimal', 'Target\nNot Maximal']
        values = [base1_failures, base2_failures, target_failures]
        
        bars = ax.bar(categories, values, color='orange', alpha=0.7)
        ax.set_title('Minimality Failure Reasons')
        ax.set_ylabel('Number of Failures')
        
        # Add value labels
        for bar, value in zip(bars, values):
            height = bar.get_height()
            if height > 0:
                ax.text(bar.get_x() + bar.get_width()/2., height + 0.05,
                       f'{value}', ha='center', va='bottom')
    
    def _plot_detailed_breakdown(self, ax):
        """Plot detailed breakdown of results."""
        if 'synergies_found' not in self.test_results:
            ax.text(0.5, 0.5, 'No data available', ha='center', va='center', transform=ax.transAxes)
            ax.set_title('Detailed Breakdown')
            return
        
        total_synergies = len(self.test_results.get('synergies_found', []))
        minimal_synergies = len(self.test_results.get('minimal_synergies', []))
        filtered_synergies = total_synergies - minimal_synergies
        
        if total_synergies == 0:
            ax.text(0.5, 0.5, 'No synergies found', ha='center', va='center', transform=ax.transAxes)
            ax.set_title('Synergy Breakdown')
            return
        
        labels = ['Minimal', 'Filtered']
        sizes = [minimal_synergies, filtered_synergies]
        colors = ['lightgreen', 'lightcoral']
        
        # Only plot if we have data
        if sum(sizes) > 0:
            wedges, texts, autotexts = ax.pie(sizes, labels=labels, colors=colors, 
                                            autopct='%1.1f%%', startangle=90)
        
        ax.set_title('Synergy Breakdown')
    
    def run_comprehensive_test(self, test_file_path: str = None):
        """Run the complete comprehensive test suite."""
        if test_file_path:
            self.rn_file = test_file_path
        
        print("Starting comprehensive synergy test...")
        
        # Setup
        self.setup()
        
        # Analyze hierarchy
        self.analyze_hierarchy_structure()
        
        # Test individual synergies
        synergies_found = self.test_individual_synergies(verbose=True)
        
        # Test minimality
        minimal_synergies = self.test_minimality_analysis(synergies_found, verbose=True)
        
        # Compare with calculator
        calc_all, calc_minimal = self.compare_algorithms()
        
        # Create visualization
        fig = self.create_comprehensive_visualization()
        
        # Final summary
        print(f"\n" + "="*80)
        print("COMPREHENSIVE TEST COMPLETE")
        print("="*80)
        
        print(f"Final Results:")
        print(f"  Total ERC combinations tested: {len(self.test_results.get('individual_tests', []))}")
        print(f"  Synergies found: {len(synergies_found)}")
        print(f"  Minimal synergies: {len(minimal_synergies)}")
        print(f"  Calculator matches manual analysis: {'✓' if len(calc_all) == len(synergies_found) and len(calc_minimal) == len(minimal_synergies) else '✗'}")
        
        if minimal_synergies:
            print(f"\nFinal Minimal Synergies:")
            for i, synergy in enumerate(minimal_synergies, 1):
                print(f"  {i}. {synergy['base1'].label} + {synergy['base2'].label} → {synergy['target'].label}")
        else:
            print(f"\nNo minimal synergies found.")
        
        return {
            'manual_synergies': synergies_found,
            'manual_minimal': minimal_synergies,
            'calculator_synergies': calc_all,
            'calculator_minimal': calc_minimal,
            'test_results': self.test_results,
            'visualization': fig
        }


def test_with_specific_network(file_path: str):
    """Test the comprehensive analyzer with a specific network file."""
    print(f"Testing comprehensive synergy analysis with: {file_path}")
    
    if not os.path.exists(file_path):
        print(f"Error: File {file_path} not found.")
        return None
    
    analyzer = SynergyTestAnalyzer(file_path)
    results = analyzer.run_comprehensive_test()
    
    # Show visualization
    if results['visualization']:
        plt.show()
    
    return analyzer, results


if __name__ == "__main__":
    # Test with your reaction network
    # Replace with the path to your test file
    test_file = "networks/RandomAlife/RN_Ns_25_Norg_14_id_206.txt"  # This should match the file in your documents
    
    try:
        analyzer, results = test_with_specific_network(test_file)
        
        # Keep plots open
        plt.show(block=True)
        
    except Exception as e:
        print(f"Error during testing: {str(e)}")
        import traceback
        traceback.print_exc()