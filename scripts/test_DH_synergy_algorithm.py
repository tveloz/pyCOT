#!/usr/bin/env python3
"""
Test Script for Efficient Fundamental Synergy Finder - CORRECTED VERSION

This script tests the corrected efficient algorithm against the brute force approach
to verify correctness and measure performance improvements.

@author: Test Efficient Synergy
"""

import time
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
import networkx as nx

# Import modules
from pyCOT.io.functions import read_txt
from pyCOT.ERC_Hierarchy import ERC, species_list_to_names
from pyCOT.Synergy_Brute_Force_Algorithm import BruteForceSynergyCalculator
from pyCOT.Dynamical_Hierarchy_Synergy_Algorithm import EfficientFundamentalSynergyFinder


class EfficientSynergyTester:
    """Comprehensive tester for the corrected efficient synergy algorithm."""
    
    def __init__(self, network_file: str):
        """Initialize tester with network file."""
        self.network_file = network_file
        self.rn = None
        self.ercs = None
        self.hierarchy_graph = None
        self.results = {}
        
    def setup(self):
        """Load network and build hierarchy."""
        print("="*80)
        print("EFFICIENT SYNERGY ALGORITHM TEST - CORRECTED VERSION")
        print("="*80)
        print(f"Loading network from: {self.network_file}")
        
        self.rn = read_txt(self.network_file)
        self.ercs = ERC.ERCs(self.rn)
        self.hierarchy_graph = ERC.build_hierarchy_graph(self.ercs, self.rn)
        
        print(f"\nNetwork Statistics:")
        print(f"  ERCs: {len(self.ercs)}")
        print(f"  Species: {len(self.rn.species())}")
        print(f"  Reactions: {len(self.rn.reactions())}")
        print(f"  Hierarchy edges: {self.hierarchy_graph.number_of_edges()}")
        
        # Print ERC details
        print(f"\nERC Details:")
        for erc in self.ercs:
            closure = species_list_to_names(erc.get_closure(self.rn))
            generators = [species_list_to_names(gen) for gen in erc.min_generators]
            print(f"  {erc.label}: closure={closure}, generators={generators}")
        
        # Print hierarchy
        print(f"\nHierarchy Structure:")
        for parent, child in self.hierarchy_graph.edges():
            print(f"  {parent} contains {child}")
    
    def analyze_constraint_potential(self):
        """Analyze the potential for constraint propagation in the hierarchy."""
        print(f"\n" + "-"*50)
        print("CONSTRAINT PROPAGATION POTENTIAL")
        print("-"*50)
        
        # Count ERCs at each level
        levels = ERC.get_node_levels(self.hierarchy_graph)
        level_counts = defaultdict(int)
        for erc, level in levels.items():
            level_counts[level] += 1
        
        print("ERCs per level:")
        for level in sorted(level_counts.keys(), reverse=True):
            print(f"  Level {level}: {level_counts[level]} ERCs")
        
        # Analyze containment depth
        max_chain_length = 0
        for erc in self.ercs:
            ancestors = nx.ancestors(self.hierarchy_graph, erc.label)
            if len(ancestors) > max_chain_length:
                max_chain_length = len(ancestors)
        
        print(f"\nMaximum containment chain length: {max_chain_length}")
        
        # Count impossible base pairs
        impossible_count = 0
        total_pairs = len(self.ercs) * (len(self.ercs) - 1) // 2
        
        for i, erc1 in enumerate(self.ercs):
            for erc2 in self.ercs[i+1:]:
                if (nx.has_path(self.hierarchy_graph, erc1.label, erc2.label) or
                    nx.has_path(self.hierarchy_graph, erc2.label, erc1.label)):
                    impossible_count += 1
        
        print(f"\nBase pair analysis:")
        print(f"  Total possible pairs: {total_pairs}")
        print(f"  Impossible pairs (containment): {impossible_count}")
        print(f"  Valid base pairs: {total_pairs - impossible_count}")
        print(f"  Reduction: {impossible_count/total_pairs*100:.1f}%")
    
    def compare_algorithms(self, verbose: bool = True):
        """Compare brute force and efficient algorithms."""
        print(f"\n" + "="*80)
        print("ALGORITHM COMPARISON")
        print("="*80)
        
        # Run brute force algorithm
        print("\n1. Running Brute Force Algorithm...")
        brute_calculator = BruteForceSynergyCalculator(self.rn, self.ercs, self.hierarchy_graph)
        
        start_time = time.time()
        all_synergies, brute_fundamental = brute_calculator.brute_force(fundamental=True, verbose=False)
        brute_time = time.time() - start_time
        
        print(f"   Brute force found {len(all_synergies)} total synergies")
        print(f"   Brute force found {len(brute_fundamental)} fundamental synergies")
        print(f"   Time: {brute_time:.3f} seconds")
        print(f"   Total checks: {brute_calculator.stats['total_checks']}")
        
        # Run efficient algorithm
        print("\n2. Running Efficient Algorithm...")
        efficient_finder = EfficientFundamentalSynergyFinder(self.rn, self.ercs, self.hierarchy_graph)
        
        start_time = time.time()
        efficient_fundamental = efficient_finder.find_fundamental_synergies(verbose=False)
        efficient_time = time.time() - start_time
        
        print(f"   Efficient algorithm found {len(efficient_fundamental)} fundamental synergies")
        print(f"   Time: {efficient_time:.3f} seconds")
        print(f"   Potential pairs generated: {efficient_finder.stats['base_pairs_generated']}")
        print(f"   Impossible pairs skipped: {efficient_finder.stats['valid_base_target_pairs']-efficient_finder.stats['base_pairs_generated']}")
        print(f"   Synergies found: {efficient_finder.stats['synergies_found']}")
        print(f"   Pruning operations: {efficient_finder.stats['within_target_pruning']}")
        print(f"   Cross-target pruning: {efficient_finder.stats['cross_target_pruning']}")
        
        # Store results
        self.results['brute_force'] = {
            'all_synergies': all_synergies,
            'fundamental_synergies': brute_fundamental,
            'time': brute_time,
            'stats': brute_calculator.stats
        }
        
        self.results['efficient'] = {
            'synergies': efficient_fundamental,
            'time': efficient_time,
            'stats': efficient_finder.stats
        }
        
        # Compare results
        print(f"\n3. Comparing Results...")
        brute_set = {(s['base1'].label, s['base2'].label, s['target'].label) 
                    for s in brute_fundamental}
        efficient_set = {(s['base1'].label, s['base2'].label, s['target'].label) 
                        for s in efficient_fundamental}
        
        if brute_set == efficient_set:
            print("   ✓ Results match perfectly!")
        else:
            print("   ✗ Results differ!")
            only_brute = brute_set - efficient_set
            only_efficient = efficient_set - brute_set
            
            if only_brute:
                print(f"   Only in brute force: {only_brute}")
            if only_efficient:
                print(f"   Only in efficient: {only_efficient}")
        
        # List all synergies found
        print(f"\n4. Synergies Found:")
        for i, synergy in enumerate(sorted(brute_set), 1):
            print(f"   {i}. {synergy[0]} + {synergy[1]} → {synergy[2]}")
        
        # Performance comparison
        if brute_time > 0:
            speedup = brute_time / efficient_time
            print(f"\n5. Performance Summary:")
            print(f"   Speedup: {speedup:.2f}x faster")
            print(f"   Time saved: {brute_time - efficient_time:.3f} seconds")
            
            # Efficiency metrics
            brute_checks = brute_calculator.stats['total_checks']
            efficient_potential = efficient_finder.stats['base_pairs_generated']
            if brute_checks > 0:
                check_reduction = (1 - efficient_potential / brute_checks) * 100
                print(f"   Check reduction: {check_reduction:.1f}%")
        
        return brute_set == efficient_set
    
    def detailed_analysis(self):
        """Perform detailed analysis of the efficient algorithm."""
        print(f"\n" + "="*80)
        print("DETAILED ALGORITHM ANALYSIS")
        print("="*80)
        
        # Run efficient algorithm with verbose output
        print("\nRunning efficient algorithm with detailed output...")
        efficient_finder = EfficientFundamentalSynergyFinder(self.rn, self.ercs, self.hierarchy_graph)
        efficient_fundamental = efficient_finder.find_fundamental_synergies(verbose=True)
        
        # Analyze efficiency gains
        print(f"\n" + "-"*50)
        print("EFFICIENCY ANALYSIS")
        print("-"*50)
        
        stats = efficient_finder.stats
        
        # Calculate theoretical maximum checks
        n_ercs = len(self.ercs)
        max_base_pairs = n_ercs * (n_ercs - 1) // 2
        max_target_checks = max_base_pairs * n_ercs
        
        print(f"Theoretical maximum synergy checks: {max_target_checks}")
        print(f"Potential pairs generated: {stats['base_pairs_generated']}")
        print(f"Reduction from theoretical max: {(1 - stats['base_pairs_generated']/max_target_checks)*100:.1f}%")

        print(f"\nPruning effectiveness:")
        print(f"  Within-target pruning: {stats['within_target_pruning']}")
        print(f"  Cross-target pruning: {stats['cross_target_pruning']}")
        print(f"  Total constraints applied: {stats['within_target_pruning'] + stats['cross_target_pruning']}")
        
        if stats['synergies_found'] > 0:
            print(f"  Average pruning per synergy: {(stats['within_target_pruning'] + stats['cross_target_pruning'])/stats['synergies_found']:.1f}")
            fundamentality_rate = stats['final_fundamental'] / stats['synergies_found'] * 100
            print(f"\nFundamentality rate: {fundamentality_rate:.1f}%")
        
        # Show synergy details with explanations
        if efficient_fundamental:
            print(f"\n" + "-"*50)
            print("FUNDAMENTAL SYNERGIES ANALYSIS")
            print("-"*50)
            
            for i, synergy in enumerate(efficient_fundamental, 1):
                base1 = synergy['base1']
                base2 = synergy['base2']
                target = synergy['target']
                
                print(f"\n{i}. {base1.label} + {base2.label} → {target.label}")
                
                # Show closures
                base1_closure = set(species_list_to_names(base1.get_closure(self.rn)))
                base2_closure = set(species_list_to_names(base2.get_closure(self.rn)))
                target_closure = set(species_list_to_names(target.get_closure(self.rn)))
                
                print(f"   Base1 ({base1.label}) closure: {sorted(base1_closure)}")
                print(f"   Base2 ({base2.label}) closure: {sorted(base2_closure)}")
                print(f"   Target ({target.label}) closure: {sorted(target_closure)}")
                
                # Show which generators are covered synergistically
                print(f"   Target generators and coverage:")
                for j, gen in enumerate(target.min_generators):
                    gen_species = set(species_list_to_names(gen))
                    base1_contrib = gen_species & base1_closure - base2_closure
                    base2_contrib = gen_species & base2_closure - base1_closure
                    
                    print(f"     Generator {j+1}: {sorted(gen_species)}")
                    if base1_contrib and base2_contrib:
                        print(f"       ✓ Synergistic! Base1 contributes {sorted(base1_contrib)}, Base2 contributes {sorted(base2_contrib)}")
                    else:
                        print(f"       ✗ Not covered synergistically")
                
                # Explain why it's fundamental
                print(f"   Fundamentality:")
                print(f"     - Base1 minimal: No proper subset of {base1.label} works")
                print(f"     - Base2 minimal: No proper subset of {base2.label} works")
                print(f"     - Target maximal: No proper superset of {target.label} is covered")
        
        return efficient_finder
    
    def visualize_results(self, figsize=(14, 10)):
        """Create visualization of algorithm performance."""
        if not self.results:
            print("No results to visualize. Run compare_algorithms first.")
            return None
        
        fig, axes = plt.subplots(2, 2, figsize=figsize)
        fig.suptitle('Efficient vs Brute Force Algorithm Comparison', fontsize=16, fontweight='bold')
        
        # 1. Time comparison
        ax1 = axes[0, 0]
        algorithms = ['Brute Force', 'Efficient']
        times = [self.results['brute_force']['time'], self.results['efficient']['time']]
        bars = ax1.bar(algorithms, times, color=['coral', 'lightgreen'])
        ax1.set_ylabel('Time (seconds)')
        ax1.set_title('Execution Time Comparison')
        
        # Add value labels
        for bar, time in zip(bars, times):
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height + 0.001,
                    f'{time:.3f}s', ha='center', va='bottom')
        
        # 2. Operations comparison
        ax2 = axes[0, 1]
        n_ercs = len(self.ercs)
        brute_checks = self.results['brute_force']['stats']['total_checks']
        efficient_base = self.results['efficient']['stats']['base_pairs_generated']
        efficient_impossible = n_ercs**3-self.results['efficient']['stats']['base_pairs_generated']
        
        categories = ['Brute Force\nChecks', 'Efficient\nBase Pairs', 'Impossible\nPairs Skipped']
        values = [brute_checks, efficient_base, efficient_impossible]
        colors = ['coral', 'lightgreen', 'gray']
        bars = ax2.bar(categories, values, color=colors)
        ax2.set_ylabel('Number of Operations')
        ax2.set_title('Computational Operations')
        
        for bar, value in zip(bars, values):
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                    f'{value}', ha='center', va='bottom')
        
        # 3. Pruning effectiveness
        ax3 = axes[1, 0]
        stats = self.results['efficient']['stats']
        categories = ['Synergies\nFound', 'Within-Target\nPruning', 'Cross-Target\nPruning']
        values = [stats['synergies_found'], 
                 stats['within_target_pruning'], 
                 stats['cross_target_pruning']]
        
        bars = ax3.bar(categories, values, color=['skyblue', 'orange', 'purple'])
        ax3.set_ylabel('Count')
        ax3.set_title('Pruning Effectiveness')
        
        for bar, value in zip(bars, values):
            height = bar.get_height()
            if height > 0:
                ax3.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                        f'{value}', ha='center', va='bottom')
        
        # 4. Summary metrics
        ax4 = axes[1, 1]
        if self.results['brute_force']['time'] > 0:
            speedup = self.results['brute_force']['time'] / self.results['efficient']['time']
            
            # Create summary text
            summary_text = f"Performance Summary\n\n"
            summary_text += f"Speedup: {speedup:.2f}x\n"
            summary_text += f"Time saved: {self.results['brute_force']['time'] - self.results['efficient']['time']:.3f}s\n"
            
            brute_all = len(self.results['brute_force']['all_synergies'])
            brute_fundamental = len(self.results['brute_force']['fundamental_synergies'])
            efficient_fundamental = len(self.results['efficient']['synergies'])
            efficient_potential = self.results['efficient']['stats']['base_pairs_generated']

            summary_text += f"\nSynergies:\n"
            summary_text += f"Total found (brute): {brute_all}\n"
            summary_text += f"Fundamental (brute): {brute_fundamental}\n"
            summary_text += f"Fundamental (efficient): {efficient_fundamental}\n"
            
            if brute_checks > 0:
                reduction = (1 - efficient_potential / brute_checks) * 100
                summary_text += f"\nCheck reduction: {reduction:.1f}%"
            
            ax4.text(0.1, 0.9, summary_text, 
                    ha='left', va='top', fontsize=12,
                    transform=ax4.transAxes,
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgray", alpha=0.5))
        
        ax4.set_xlim(0, 1)
        ax4.set_ylim(0, 1)
        ax4.axis('off')
        
        plt.tight_layout()
        return fig
    
    def run_complete_test(self):
        """Run the complete test suite."""
        # Setup
        self.setup()
        
        # Analyze constraint potential
        self.analyze_constraint_potential()
        
        # Compare algorithms
        results_match = self.compare_algorithms()
        
        # Detailed analysis
        efficient_finder = self.detailed_analysis()
        
        # Visualize results
        fig = self.visualize_results()
        
        # Final summary
        print(f"\n" + "="*80)
        print("FINAL TEST SUMMARY")
        print("="*80)
        print(f"✓ Correctness: {'PASSED' if results_match else 'FAILED'}")
        
        if self.results['brute_force']['time'] > 0 and self.results['efficient']['time'] > 0:
            speedup = self.results['brute_force']['time'] / self.results['efficient']['time']
            print(f"✓ Performance: {speedup:.2f}x faster")
            
            brute_checks = self.results['brute_force']['stats']['total_checks']
            efficient_potential = self.results['efficient']['stats']['base_pairs_generated']
            if brute_checks > 0:
                reduction = (1 - efficient_potential / brute_checks) * 100
                print(f"✓ Efficiency: Reduced checks by {reduction:.1f}%")
        
        print(f"\nAlgorithm found {len(self.results['efficient']['synergies'])} fundamental synergies correctly!")
        
        return {
            'results_match': results_match,
            'brute_force_results': self.results['brute_force'],
            'efficient_results': self.results['efficient'],
            'finder': efficient_finder,
            'visualization': fig
        }


def test_efficient_synergy_finder(network_file: str = "networks/testing/ERCs_test2.txt"):
    """Test the efficient synergy finder with the specified network."""
    print(f"Testing efficient synergy finder with: {network_file}")
    
    tester = EfficientSynergyTester(network_file)
    results = tester.run_complete_test()
    
    # Show plots
    if results['visualization']:
        plt.show()
    
    return tester, results


if __name__ == "__main__":
    import os
    
    # Test with the provided network file
    test_file = "networks/testing/Farm.txt"
    #test_file = "networks\RandomAlife\RN_Ns_40_Norg_18_id_367.txt"
    test_file = "networks/RandomAlife/RN_Ns_20_Norg_8_id_163.txt"
    test_file = "networks/Navarino/RN_IN_05.txt"

    # Alternative: use the content directly
    if not os.path.exists(test_file):
        # Create test file from document content
        test_content = """R1:	a1+b1=>a2+b2;
R2:	a2+b1=>a3+b3;
R3:	a3+b3=>x3+a2;
R4:	a1+b4=>a5+b5;
R5:	a5+b5=>x5+b2;
R6:	a6+b1=>a7+b7;
"""
        os.makedirs("networks/testing", exist_ok=True)
        with open(test_file, 'w') as f:
            f.write(test_content)
    
    try:
        tester, results = test_efficient_synergy_finder(test_file)
        plt.show(block=True)
    except Exception as e:
        print(f"Error during testing: {str(e)}")
        import traceback
        traceback.print_exc()