#!/usr/bin/env python3
"""
Efficient Fundamental Synergy Finder Analysis Script

This script provides detailed analysis of the efficient synergy algorithm,
showing how it works internally, its optimization strategies, and performance characteristics.

@author: Efficient Synergy Analysis
"""

import time
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
import networkx as nx

# Import modules
from pyCOT.io.functions import read_txt
from pyCOT.ERC_Hierarchy import ERC, species_list_to_names
from pyCOT.Dynamical_Hierarchy_Synergy_Algorithm import EfficientFundamentalSynergyFinder


class EfficientSynergyAnalyzer:
    """Comprehensive analyzer for the efficient synergy algorithm."""
    
    def __init__(self, network_file: str):
        """Initialize analyzer with network file."""
        self.network_file = network_file
        self.rn = None
        self.ercs = None
        self.hierarchy_graph = None
        self.finder = None
        self.results = {}
        
    def setup(self):
        """Load network and build hierarchy."""
        print("="*80)
        print("EFFICIENT SYNERGY ALGORITHM ANALYSIS")
        print("="*80)
        print(f"Loading network from: {self.network_file}")
        
        self.rn = read_txt(self.network_file)
        "File loaded"
        self.ercs = ERC.ERCs(self.rn)
        self.hierarchy_graph = ERC.build_hierarchy_graph2(self.ercs, self.rn)
        
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
            print(f"  {erc.label}: closure={closure}")
            print(f"    generators={generators}")
        
        # Print hierarchy
        print(f"\nHierarchy Structure:")
        if self.hierarchy_graph.number_of_edges() > 0:
            for parent, child in self.hierarchy_graph.edges():
                print(f"  {parent} contains {child}")
        else:
            print("  No containment relationships (flat hierarchy)")
    
    def analyze_hierarchy_structure(self):
        """Analyze hierarchy structure and its impact on algorithm efficiency."""
        print(f"\n" + "="*60)
        print("HIERARCHY STRUCTURE ANALYSIS")
        print("="*60)
        
        # Count ERCs at each level
        levels = ERC.get_node_levels(self.hierarchy_graph)
        level_counts = defaultdict(int)
        level_ercs = defaultdict(list)
        
        for erc_label, level in levels.items():
            level_counts[level] += 1
            level_ercs[level].append(erc_label)
        
        print("ERCs per hierarchy level:")
        for level in sorted(level_counts.keys(), reverse=True):
            ercs_at_level = sorted(level_ercs[level])
            print(f"  Level {level}: {level_counts[level]} ERCs - {ercs_at_level}")
        
        # Analyze containment structure
        max_chain_length = 0
        longest_chains = []
        for erc_label in self.hierarchy_graph.nodes():
            ancestors = nx.ancestors(self.hierarchy_graph, erc_label)
            chain_length = len(ancestors)
            if chain_length > max_chain_length:
                max_chain_length = chain_length
                longest_chains = [erc_label]
            elif chain_length == max_chain_length:
                longest_chains.append(erc_label)
        
        print(f"\nContainment analysis:")
        print(f"  Maximum chain length: {max_chain_length}")
        print(f"  ERCs at maximum depth: {longest_chains}")
        
        # Count branches and roots
        roots = [n for n in self.hierarchy_graph.nodes() if self.hierarchy_graph.in_degree(n) == 0]
        leaves = [n for n in self.hierarchy_graph.nodes() if self.hierarchy_graph.out_degree(n) == 0]
        
        print(f"  Root ERCs (no parents): {len(roots)} - {roots}")
        print(f"  Leaf ERCs (no children): {len(leaves)} - {leaves}")
        
        # Calculate potential base pair reductions
        impossible_pairs = 0
        total_pairs = len(self.ercs) * (len(self.ercs) - 1) // 2
        
        for i, erc1 in enumerate(self.ercs):
            for erc2 in self.ercs[i+1:]:
                if (nx.has_path(self.hierarchy_graph, erc1.label, erc2.label) or
                    nx.has_path(self.hierarchy_graph, erc2.label, erc1.label)):
                    impossible_pairs += 1
        
        print(f"\nBase pair feasibility:")
        print(f"  Total possible pairs: {total_pairs}")
        print(f"  Impossible pairs (containment): {impossible_pairs}")
        print(f"  Valid base pairs: {total_pairs - impossible_pairs}")
        print(f"  Hierarchy reduction: {impossible_pairs/total_pairs*100:.1f}%")
        
        return {
            'levels': levels,
            'max_chain_length': max_chain_length,
            'roots': roots,
            'leaves': leaves,
            'impossible_pairs': impossible_pairs,
            'total_pairs': total_pairs
        }
    
    def analyze_partial_overlap_strategy(self):
        """Analyze the partial overlap filtering strategy."""
        print(f"\n" + "="*60)
        print("PARTIAL OVERLAP FILTERING ANALYSIS")
        print("="*60)
        
        print("Analyzing partial overlap between ERCs and target generators...")
        
        overlap_analysis = {}
        total_base_target_checks = 0
        valid_base_target_pairs = 0
        
        for target in self.ercs:
            overlap_analysis[target.label] = {
                'generators': [species_list_to_names(gen) for gen in target.min_generators],
                'valid_bases': [],
                'rejected_bases': []
            }
            
            for base in self.ercs:
                if base.label == target.label:
                    continue
                    
                total_base_target_checks += 1
                
                # Check if base contains target (impossible synergy)
                ancestors = nx.ancestors(self.hierarchy_graph, target.label)
                if base.label in ancestors:
                    overlap_analysis[target.label]['rejected_bases'].append({
                        'base': base.label,
                        'reason': 'containment'
                    })
                    continue
                
                # Check partial overlap with generators
                base_closure = set(species_list_to_names(base.get_closure(self.rn)))
                has_partial_overlap = False
                
                for gen in target.min_generators:
                    gen_species = set(species_list_to_names(gen))
                    intersection = base_closure & gen_species
                    
                    if intersection and not gen_species.issubset(base_closure):
                        has_partial_overlap = True
                        break
                
                if has_partial_overlap:
                    valid_base_target_pairs += 1
                    overlap_analysis[target.label]['valid_bases'].append(base.label)
                else:
                    overlap_analysis[target.label]['rejected_bases'].append({
                        'base': base.label,
                        'reason': 'no_partial_overlap'
                    })
        
        # Print results
        print(f"Partial overlap filtering results:")
        print(f"  Total base-target checks: {total_base_target_checks}")
        print(f"  Valid base-target pairs: {valid_base_target_pairs}")
        print(f"  Rejection rate: {(1 - valid_base_target_pairs/total_base_target_checks)*100:.1f}%")
        
        print(f"\nDetailed analysis per target:")
        for target_label, analysis in overlap_analysis.items():
            valid_count = len(analysis['valid_bases'])
            rejected_count = len(analysis['rejected_bases'])
            total_count = valid_count + rejected_count
            
            print(f"\n  Target {target_label}:")
            print(f"    Generators: {analysis['generators']}")
            print(f"    Valid bases: {valid_count}/{total_count} - {analysis['valid_bases']}")
            
            if rejected_count > 0:
                containment_rejects = sum(1 for r in analysis['rejected_bases'] if r['reason'] == 'containment')
                overlap_rejects = sum(1 for r in analysis['rejected_bases'] if r['reason'] == 'no_partial_overlap')
                print(f"    Rejected: {rejected_count} (containment: {containment_rejects}, no overlap: {overlap_rejects})")
        
        return overlap_analysis
    
    def run_algorithm_with_analysis(self):
        """Run the efficient algorithm with detailed step-by-step analysis."""
        print(f"\n" + "="*60)
        print("RUNNING EFFICIENT ALGORITHM WITH DETAILED ANALYSIS")
        print("="*60)
        
        self.finder = EfficientFundamentalSynergyFinder(self.rn, self.ercs, self.hierarchy_graph)
        
        print("Running algorithm with verbose output...")
        start_time = time.time()
        fundamental_synergies = self.finder.find_fundamental_synergies(verbose=True)
        end_time = time.time()
        
        execution_time = end_time - start_time
        stats = self.finder.stats
        
        print(f"\n" + "-"*50)
        print("ALGORITHM EXECUTION SUMMARY")
        print("-"*50)
        
        print(f"Execution time: {execution_time:.3f} seconds")
        print(f"Fundamental synergies found: {len(fundamental_synergies)}")
        
        print(f"\nDetailed statistics:")
        for key, value in stats.items():
            print(f"  {key}: {value}")
        
        # Calculate efficiency metrics
        n_ercs = len(self.ercs)
        theoretical_max_checks = n_ercs * (n_ercs - 1) * (n_ercs - 2) // 2
        
        print(f"\nEfficiency metrics:")
        print(f"  Theoretical maximum checks: {theoretical_max_checks}")
        print(f"  Base pairs generated: {stats['base_pairs_generated']}")
        print(f"  Reduction from maximum: {(1 - stats['base_pairs_generated']/theoretical_max_checks)*100:.1f}%")
        
        if stats['synergies_found'] > 0:
            pruning_per_synergy = (stats['within_target_pruning'] + stats['cross_target_pruning']) / stats['synergies_found']
            print(f"  Average pruning operations per synergy: {pruning_per_synergy:.1f}")
            
            fundamentality_rate = len(fundamental_synergies) / stats['synergies_found'] * 100
            print(f"  Fundamentality rate: {fundamentality_rate:.1f}%")
        
        self.results = {
            'fundamental_synergies': fundamental_synergies,
            'execution_time': execution_time,
            'stats': stats,
            'finder': self.finder
        }
        
        return fundamental_synergies
    
    def analyze_synergy_details(self, fundamental_synergies):
        """Provide detailed analysis of each discovered synergy."""
        if not fundamental_synergies:
            print(f"\nNo fundamental synergies found to analyze.")
            return
        
        print(f"\n" + "="*60)
        print("DETAILED SYNERGY ANALYSIS")
        print("="*60)
        
        for i, synergy in enumerate(fundamental_synergies, 1):
            base1 = synergy['base1']
            base2 = synergy['base2']
            target = synergy['target']
            
            print(f"\n{i}. SYNERGY: {base1.label} + {base2.label} → {target.label}")
            print(f"   " + "-"*50)
            
            # Show closures
            base1_closure = set(species_list_to_names(base1.get_closure(self.rn)))
            base2_closure = set(species_list_to_names(base2.get_closure(self.rn)))
            target_closure = set(species_list_to_names(target.get_closure(self.rn)))
            
            print(f"   Base1 ({base1.label}) closure: {sorted(base1_closure)}")
            print(f"   Base2 ({base2.label}) closure: {sorted(base2_closure)}")
            print(f"   Target ({target.label}) closure: {sorted(target_closure)}")
            
            # Calculate combined closure
            combined_species = list(set(base1.get_closure(self.rn)) | set(base2.get_closure(self.rn)))
            from pyCOT.ERC_Hierarchy import closure
            combined_closure = set(species_list_to_names(closure(self.rn, combined_species)))
            print(f"   Combined closure: {sorted(combined_closure)}")
            
            # Verify basic synergy conditions
            print(f"\n   Synergy validation:")
            covers_target = combined_closure.issuperset(target_closure)
            base1_insufficient = not base1_closure.issuperset(target_closure)
            base2_insufficient = not base2_closure.issuperset(target_closure)
            
            print(f"     ✓ Combined covers target: {covers_target}")
            print(f"     ✓ Base1 alone insufficient: {base1_insufficient}")
            print(f"     ✓ Base2 alone insufficient: {base2_insufficient}")
            
            # Analyze generator coverage
            print(f"\n   Generator analysis:")
            synergistic_count = 0
            for j, gen in enumerate(target.min_generators):
                gen_species = set(species_list_to_names(gen))
                
                covered = gen_species.issubset(combined_closure)
                base1_contrib = gen_species & base1_closure - base2_closure
                base2_contrib = gen_species & base2_closure - base1_closure
                
                is_synergistic = (covered and 
                                not gen_species.issubset(base1_closure) and
                                not gen_species.issubset(base2_closure) and
                                len(base1_contrib) > 0 and len(base2_contrib) > 0)
                
                if is_synergistic:
                    synergistic_count += 1
                
                print(f"     Generator {j+1}: {sorted(gen_species)}")
                print(f"       Covered by combination: {covered}")
                if covered:
                    print(f"       Base1 unique contribution: {sorted(base1_contrib) if base1_contrib else 'none'}")
                    print(f"       Base2 unique contribution: {sorted(base2_contrib) if base2_contrib else 'none'}")
                    print(f"       Truly synergistic: {is_synergistic}")
            
            print(f"\n   Summary: {synergistic_count}/{len(target.min_generators)} generators require true synergy")
            
            # Explain fundamentality
            print(f"\n   Fundamentality explanation:")
            print(f"     - Base1 minimal: No ERC contained in {base1.label} can replace it")
            print(f"     - Base2 minimal: No ERC contained in {base2.label} can replace it")
            print(f"     - Target maximal: No ERC containing {target.label} can be reached instead")
    
    def analyze_pruning_effectiveness(self):
        """Analyze how effective the pruning strategies are."""
        if not self.finder:
            print("Run algorithm first to analyze pruning effectiveness.")
            return
        
        print(f"\n" + "="*60)
        print("PRUNING STRATEGY ANALYSIS")
        print("="*60)
        
        stats = self.finder.stats
        
        print(f"Pruning statistics:")
        print(f"  Within-target pruning operations: {stats['within_target_pruning']}")
        print(f"  Cross-target pruning operations: {stats['cross_target_pruning']}")
        print(f"  Total constraint applications: {stats['within_target_pruning'] + stats['cross_target_pruning']}")
        
        if stats['synergies_found'] > 0:
            avg_pruning = (stats['within_target_pruning'] + stats['cross_target_pruning']) / stats['synergies_found']
            print(f"  Average pruning per synergy discovered: {avg_pruning:.1f}")
        
        print(f"\nPruning strategy effectiveness:")
        print(f"  Within-target pruning: Eliminates non-minimal base combinations for same target")
        print(f"  Cross-target pruning: Eliminates redundant checks for contained targets")
        print(f"  Dynamic constraint propagation: Guides search away from dominated solutions")
        
        # Estimate computational savings
        n_ercs = len(self.ercs)
        total_pruning = stats['within_target_pruning'] + stats['cross_target_pruning']
        base_pairs_generated = stats['base_pairs_generated']
        
        print(f"\nComputational impact:")
        print(f"  Base pairs that would be generated without pruning: {base_pairs_generated + total_pruning}")
        print(f"  Base pairs actually generated: {base_pairs_generated}")
        print(f"  Pruning saved: {total_pruning} operations ({total_pruning/(base_pairs_generated + total_pruning)*100:.1f}%)")
    
    def visualize_algorithm_performance(self, figsize=(16, 12)):
        """Create comprehensive visualization of algorithm performance."""
        if not self.results:
            print("No results to visualize. Run algorithm first.")
            return None
        
        fig, axes = plt.subplots(2, 3, figsize=figsize)
        fig.suptitle('Efficient Synergy Algorithm Performance Analysis', fontsize=16, fontweight='bold')
        
        stats = self.results['stats']
        
        # 1. Algorithm steps breakdown
        ax1 = axes[0, 0]
        steps = ['Partial Overlap\nChecks', 'Valid Base-Target\nPairs', 'Base Pairs\nGenerated', 'Synergies\nFound']
        values = [stats['partial_overlap_checks'], stats['valid_base_target_pairs'], 
                 stats['base_pairs_generated'], stats['synergies_found']]
        
        bars = ax1.bar(steps, values, color=['lightcoral', 'lightblue', 'lightgreen', 'gold'])
        ax1.set_ylabel('Count')
        ax1.set_title('Algorithm Steps')
        plt.setp(ax1.get_xticklabels(), rotation=45, ha='right')
        
        # Add value labels
        for bar, value in zip(bars, values):
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height + max(values)*0.01,
                    f'{value}', ha='center', va='bottom')
        
        # 2. Pruning effectiveness
        ax2 = axes[0, 1]
        pruning_categories = ['Within-Target\nPruning', 'Cross-Target\nPruning', 'Final\nFundamental']
        pruning_values = [stats['within_target_pruning'], stats['cross_target_pruning'], 
                         stats['final_fundamental']]
        
        colors = ['orange', 'purple', 'green']
        bars = ax2.bar(pruning_categories, pruning_values, color=colors)
        ax2.set_ylabel('Count')
        ax2.set_title('Pruning and Results')
        plt.setp(ax2.get_xticklabels(), rotation=45, ha='right')
        
        for bar, value in zip(bars, pruning_values):
            height = bar.get_height()
            if height > 0:
                ax2.text(bar.get_x() + bar.get_width()/2., height + max(pruning_values)*0.01,
                        f'{value}', ha='center', va='bottom')
        
        # 3. Efficiency funnel
        ax3 = axes[0, 2]
        n_ercs = len(self.ercs)
        theoretical_max = n_ercs * (n_ercs - 1) * (n_ercs - 2) // 2
        
        funnel_steps = ['Theoretical\nMaximum', 'After Hierarchy\nFiltering', 'Base Pairs\nGenerated', 'Synergies\nFound']
        funnel_values = [theoretical_max, stats['valid_base_target_pairs'], 
                        stats['base_pairs_generated'], stats['synergies_found']]
        
        bars = ax3.bar(funnel_steps, funnel_values, color=['red', 'orange', 'yellow', 'green'])
        ax3.set_ylabel('Count (log scale)')
        ax3.set_yscale('log')
        ax3.set_title('Efficiency Funnel')
        plt.setp(ax3.get_xticklabels(), rotation=45, ha='right')
        
        for bar, value in zip(bars, funnel_values):
            if value > 0:
                height = bar.get_height()
                ax3.text(bar.get_x() + bar.get_width()/2., height * 1.1,
                        f'{value}', ha='center', va='bottom')
        
        # 4. Hierarchy structure impact
        ax4 = axes[1, 0]
        levels = ERC.get_node_levels(self.hierarchy_graph)
        level_counts = defaultdict(int)
        for level in levels.values():
            level_counts[level] += 1
        
        if level_counts:
            level_nums = sorted(level_counts.keys(), reverse=True)
            level_ercs = [level_counts[level] for level in level_nums]
            
            bars = ax4.bar([f'Level {l}' for l in level_nums], level_ercs, color='skyblue')
            ax4.set_ylabel('Number of ERCs')
            ax4.set_title('Hierarchy Distribution')
            plt.setp(ax4.get_xticklabels(), rotation=45, ha='right')
            
            for bar, value in zip(bars, level_ercs):
                height = bar.get_height()
                ax4.text(bar.get_x() + bar.get_width()/2., height + 0.05,
                        f'{value}', ha='center', va='bottom')
        
        # 5. Time complexity analysis
        ax5 = axes[1, 1]
        complexity_components = ['Hierarchy\nTraversal', 'Partial Overlap\nFiltering', 'Dynamic\nPruning', 'Verification']
        # Estimate relative time contributions (this is illustrative)
        time_percentages = [15, 40, 30, 15]
        
        wedges, texts, autotexts = ax5.pie(time_percentages, labels=complexity_components, 
                                          autopct='%1.1f%%', startangle=90,
                                          colors=['lightcoral', 'lightblue', 'lightgreen', 'gold'])
        ax5.set_title('Estimated Time Distribution')
        
        # 6. Performance summary
        ax6 = axes[1, 2]
        execution_time = self.results['execution_time']
        fundamental_count = len(self.results['fundamental_synergies'])
        
        summary_text = f"Performance Summary\n\n"
        summary_text += f"Execution time: {execution_time:.3f}s\n"
        summary_text += f"Fundamental synergies: {fundamental_count}\n"
        summary_text += f"ERCs processed: {n_ercs}\n"
        summary_text += f"Base pairs generated: {stats['base_pairs_generated']}\n"
        
        if theoretical_max > 0:
            reduction = (1 - stats['base_pairs_generated'] / theoretical_max) * 100
            summary_text += f"\nEfficiency gain: {reduction:.1f}%\n"
        
        summary_text += f"\nKey optimizations:\n"
        summary_text += f"• Hierarchy-based filtering\n"
        summary_text += f"• Partial overlap detection\n"
        summary_text += f"• Dynamic constraint propagation\n"
        summary_text += f"• Bottom-up exploration"
        
        ax6.text(0.05, 0.95, summary_text, 
                ha='left', va='top', fontsize=10,
                transform=ax6.transAxes,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgray", alpha=0.8))
        
        ax6.set_xlim(0, 1)
        ax6.set_ylim(0, 1)
        ax6.axis('off')
        
        plt.tight_layout()
        return fig
    
    def run_complete_analysis(self):
        """Run the complete analysis of the efficient algorithm."""
        # Setup
        self.setup()
        
        # Analyze hierarchy structure
        hierarchy_analysis = self.analyze_hierarchy_structure()
        
        # Analyze partial overlap strategy
        overlap_analysis = self.analyze_partial_overlap_strategy()
        
        # Run algorithm with detailed analysis
        fundamental_synergies = self.run_algorithm_with_analysis()
        
        # Analyze discovered synergies
        self.analyze_synergy_details(fundamental_synergies)
        
        # Analyze pruning effectiveness
        self.analyze_pruning_effectiveness()
        
        # Create visualization
        fig = self.visualize_algorithm_performance()
        
        # Final summary
        print(f"\n" + "="*80)
        print("COMPLETE ANALYSIS SUMMARY")
        print("="*80)
        
        stats = self.results['stats']
        n_ercs = len(self.ercs)
        theoretical_max = n_ercs * (n_ercs - 1) * (n_ercs - 2) // 2
        
        print(f"Algorithm Performance:")
        print(f"  Execution time: {self.results['execution_time']:.3f} seconds")
        print(f"  Fundamental synergies found: {len(fundamental_synergies)}")
        
        print(f"\nEfficiency Achievements:")
        if theoretical_max > 0:
            reduction = (1 - stats['base_pairs_generated'] / theoretical_max) * 100
            print(f"  Reduced search space by: {reduction:.1f}%")
        print(f"  Hierarchy filtering eliminated: {hierarchy_analysis['impossible_pairs']} impossible pairs")
        print(f"  Partial overlap filtering: {(1 - stats['valid_base_target_pairs']/stats['partial_overlap_checks'])*100:.1f}% rejection rate")
        print(f"  Dynamic pruning operations: {stats['within_target_pruning'] + stats['cross_target_pruning']}")
        
        print(f"\nKey Algorithm Insights:")
        print(f"  • Hierarchy structure provides {hierarchy_analysis['impossible_pairs']/(hierarchy_analysis['total_pairs'])*100:.1f}% base pair reduction")
        print(f"  • Partial overlap filtering is the primary efficiency driver")
        print(f"  • Dynamic pruning prevents redundant computation")
        print(f"  • Bottom-up exploration ensures minimality")
        
        return {
            'fundamental_synergies': fundamental_synergies,
            'hierarchy_analysis': hierarchy_analysis,
            'overlap_analysis': overlap_analysis,
            'results': self.results,
            'visualization': fig
        }


def analyze_efficient_synergy_algorithm(network_file: str):
    """Analyze the efficient synergy algorithm with the specified network."""
    print(f"Analyzing efficient synergy algorithm with: {network_file}")
    
    analyzer = EfficientSynergyAnalyzer(network_file)
    results = analyzer.run_complete_analysis()
    
    # Show plots
    if results['visualization']:
        plt.show()
    
    return analyzer, results


if __name__ == "__main__":
    import os
    
    # Test with various network files
    test_files = [
        "networks/biomodels_interesting/bigg_iAF692.txt"
        #'networks/Navarino/RN_IN_05.txt'
    ]
    
    # Use the first available file or create a default one
    test_file = None
    for file in test_files:
        if os.path.exists(file):
            test_file = file
            break
    
    if not test_file:
        # Create test file from document content
        
        test_content = """R1:	a1+b1=>a2+b2;
R2:	a2+b1=>a3+b3;
R3:	a3+b3=>x3+a2;
R4:	a1+b4=>a5+b5;
R5:	a5+b5=>x5+b2;
R6:	a6+b1=>a7+b7;
"""
        with open(test_file, 'w') as f:
            f.write(test_content)
    
    try:
        analyzer, results = analyze_efficient_synergy_algorithm(test_file)
        plt.show(block=True)
    except Exception as e:
        print(f"Error during analysis: {str(e)}")
        import traceback
        traceback.print_exc()