#!/usr/bin/env python3
"""
Synergy Algorithm Comparison Script

This script compares the results from brute force and hierarchy-aware synergy algorithms,
measuring performance and identifying any differences in the synergies found.

@author: ERC Synergy Comparison Framework
"""

import time
import numpy as np
from collections import defaultdict
from typing import Dict, List, Set, Tuple
import pandas as pd
from pathlib import Path

# Import required modules
from pyCOT.io.functions import read_txt
from pyCOT.ERC_Hierarchy import ERC, species_list_to_names
from pyCOT.Synergy_Brute_Force_Algorithm import BruteForceSynergyCalculator
from pyCOT.ERC_Exploration import ERCExplorer
from pyCOT.Synergy_Hierarchy_Aware_Algorithm import SynergyHarvester


class SynergyAlgorithmComparator:
    """Compare synergy results from different algorithms."""
    
    def __init__(self, reaction_network_file: str):
        """Initialize comparator with a reaction network file."""
        print(f"Loading reaction network from: {reaction_network_file}")
        self.rn = read_txt(reaction_network_file)
        
        # Generate ERCs and hierarchy
        print("Generating ERCs...")
        self.ercs = ERC.ERCs(self.rn)
        print(f"  Found {len(self.ercs)} ERCs")
        
        print("Building hierarchy...")
        self.hierarchy = ERC.build_hierarchy_graph(self.ercs, self.rn)
        
        # Results storage
        self.brute_force_results = None
        self.hierarchy_aware_results = None
        self.comparison_results = None
        
    def normalize_synergy_format(self, synergy_data):
        """
        Normalize synergy data to a common format for comparison.
        Returns a set of tuples: (base1_label, base2_label, target_label)
        """
        normalized_synergies = set()
        
        # Handle brute force format
        if isinstance(synergy_data, list) and synergy_data and 'base1' in synergy_data[0]:
            for synergy in synergy_data:
                base1 = synergy['base1'].label
                base2 = synergy['base2'].label
                target = synergy['target'].label
                
                # Ensure consistent ordering of base1 and base2
                if base1 > base2:
                    base1, base2 = base2, base1
                
                normalized_synergies.add((base1, base2, target))
        
        # Handle hierarchy-aware format (target-centric)
        elif isinstance(synergy_data, dict):
            for target_erc, synergy_list in synergy_data.items():
                target = target_erc.label
                for synergy in synergy_list:
                    base1 = synergy.base1_erc.label
                    base2 = synergy.base2_erc.label
                    
                    # Ensure consistent ordering
                    if base1 > base2:
                        base1, base2 = base2, base1
                    
                    normalized_synergies.add((base1, base2, target))
        
        return normalized_synergies
    
    def run_brute_force_algorithm(self, minimal=False, verbose=False):
        """Run the brute force synergy algorithm and measure time."""
        print("\n" + "="*70)
        print("RUNNING BRUTE FORCE ALGORITHM")
        print("="*70)
        
        calculator = BruteForceSynergyCalculator(self.rn, self.ercs, self.hierarchy)
        
        start_time = time.time()
        
        if minimal:
            all_synergies, minimal_synergies = calculator.brute_force(minimal=True, verbose=verbose)
            self.brute_force_results = {
                'all': all_synergies,
                'minimal': minimal_synergies
            }
        else:
            all_synergies = calculator.brute_force(minimal=False, verbose=verbose)
            self.brute_force_results = {
                'all': all_synergies,
                'minimal': []
            }
        
        end_time = time.time()
        
        execution_time = end_time - start_time
        
        print(f"\nBrute Force Results:")
        print(f"  Execution time: {execution_time:.3f} seconds")
        print(f"  Total synergies found: {len(all_synergies)}")
        if minimal:
            print(f"  Minimal synergies: {len(minimal_synergies)}")
        print(f"  Statistics: {calculator.stats}")
        
        return {
            'synergies': self.brute_force_results,
            'execution_time': execution_time,
            'stats': calculator.stats
        }
    
    def run_hierarchy_aware_algorithm(self, verbose=False):
        """Run the hierarchy-aware synergy algorithm and measure time."""
        print("\n" + "="*70)
        print("RUNNING HIERARCHY-AWARE ALGORITHM")
        print("="*70)
        
        # Phase 1: ERC Exploration
        print("\nPhase 1: ERC Exploration")
        explorer = ERCExplorer(self.rn, self.ercs, self.hierarchy)
        
        start_time = time.time()
        
        # Identify potential synergy targets
        potential_synergies = explorer.identify_potential_synergy_targets()
        
        # Reorganize to target-centric
        target_centric_synergies = explorer.reorganize_to_target_centric(potential_synergies)
        
        exploration_time = time.time() - start_time
        
        print(f"  Exploration time: {exploration_time:.3f} seconds")
        print(f"  Targets identified: {len(target_centric_synergies)}")
        
        # Phase 2: Synergy Harvesting
        print("\nPhase 2: Synergy Harvesting")
        harvester = SynergyHarvester(self.rn, self.hierarchy)
        
        harvest_start = time.time()
        
        effective_synergies = harvester.harvest_effective_synergies(target_centric_synergies, verbose=verbose)
        
        harvest_time = time.time() - harvest_start
        total_time = time.time() - start_time
        
        self.hierarchy_aware_results = effective_synergies
        
        # Count total synergies
        total_synergies = sum(len(synergies) for synergies in effective_synergies.values())
        
        print(f"\nHierarchy-Aware Results:")
        print(f"  Total execution time: {total_time:.3f} seconds")
        print(f"    - Exploration: {exploration_time:.3f} seconds")
        print(f"    - Harvesting: {harvest_time:.3f} seconds")
        print(f"  Total effective synergies found: {total_synergies}")
        print(f"  Harvesting statistics: {harvester.harvesting_statistics}")
        
        return {
            'synergies': effective_synergies,
            'execution_time': total_time,
            'exploration_time': exploration_time,
            'harvest_time': harvest_time,
            'stats': harvester.harvesting_statistics
        }
    
    def compare_results(self, brute_force_data, hierarchy_aware_data):
        """Compare the results from both algorithms."""
        print("\n" + "="*70)
        print("COMPARING ALGORITHM RESULTS")
        print("="*70)
        
        # Normalize synergy formats
        bf_all = self.normalize_synergy_format(brute_force_data['synergies']['all'])
        bf_minimal = self.normalize_synergy_format(brute_force_data['synergies']['minimal'])
        ha_synergies = self.normalize_synergy_format(hierarchy_aware_data['synergies'])
        
        # Basic comparison
        print(f"\n📊 Basic Comparison:")
        print(f"  Brute Force (all): {len(bf_all)} synergies")
        print(f"  Brute Force (minimal): {len(bf_minimal)} synergies")
        print(f"  Hierarchy-Aware: {len(ha_synergies)} synergies")
        
        # Time comparison
        print(f"\n⏱️  Time Comparison:")
        print(f"  Brute Force: {brute_force_data['execution_time']:.3f} seconds")
        print(f"  Hierarchy-Aware: {hierarchy_aware_data['execution_time']:.3f} seconds")
        speedup = brute_force_data['execution_time'] / hierarchy_aware_data['execution_time']
        print(f"  Speedup: {speedup:.2f}x")
        
        # Set operations for detailed comparison
        print(f"\n🔍 Detailed Comparison:")
        
        # Compare with all brute force synergies
        common_with_all = bf_all & ha_synergies
        bf_only_all = bf_all - ha_synergies
        ha_only_vs_all = ha_synergies - bf_all
        
        print(f"\n  Compared to Brute Force (all):")
        print(f"    Common synergies: {len(common_with_all)}")
        print(f"    Only in Brute Force: {len(bf_only_all)}")
        print(f"    Only in Hierarchy-Aware: {len(ha_only_vs_all)}")
        
        # Compare with minimal brute force synergies
        if bf_minimal:
            common_with_minimal = bf_minimal & ha_synergies
            bf_only_minimal = bf_minimal - ha_synergies
            ha_only_vs_minimal = ha_synergies - bf_minimal
            
            print(f"\n  Compared to Brute Force (minimal):")
            print(f"    Common synergies: {len(common_with_minimal)}")
            print(f"    Only in Brute Force (minimal): {len(bf_only_minimal)}")
            print(f"    Only in Hierarchy-Aware: {len(ha_only_vs_minimal)}")
        
        # Store comparison results
        self.comparison_results = {
            'bf_all': bf_all,
            'bf_minimal': bf_minimal,
            'ha_synergies': ha_synergies,
            'common_with_all': common_with_all,
            'bf_only_all': bf_only_all,
            'ha_only_vs_all': ha_only_vs_all,
            'time_speedup': speedup
        }
        
        # Investigate differences if any
        if bf_only_all or ha_only_vs_all:
            self.investigate_differences(bf_only_all, ha_only_vs_all, brute_force_data, hierarchy_aware_data)
        else:
            print("\n✅ Both algorithms found exactly the same synergies!")
        
        return self.comparison_results
    
    def investigate_differences(self, bf_only, ha_only, bf_data, ha_data):
        """Investigate why certain synergies appear in only one algorithm."""
        print("\n" + "="*70)
        print("INVESTIGATING DIFFERENCES")
        print("="*70)
        
        # Helper to get ERC by label
        def get_erc(label):
            return next(erc for erc in self.ercs if erc.label == label)
        
        # Investigate synergies only in brute force
        if bf_only:
            print(f"\n🔍 Synergies only in Brute Force ({len(bf_only)}):")
            for i, (base1_label, base2_label, target_label) in enumerate(sorted(bf_only)[:5]):  # Show first 5
                print(f"\n  [{i+1}] {base1_label} + {base2_label} → {target_label}")
                
                # Get ERCs
                base1 = get_erc(base1_label)
                base2 = get_erc(base2_label)
                target = get_erc(target_label)
                
                # Check why it might be missing from hierarchy-aware
                base1_closure = set(species_list_to_names(base1.get_closure(self.rn)))
                base2_closure = set(species_list_to_names(base2.get_closure(self.rn)))
                target_closure = set(species_list_to_names(target.get_closure(self.rn)))
                combined = base1_closure | base2_closure
                
                print(f"    Base1 closure: {sorted(base1_closure)}")
                print(f"    Base2 closure: {sorted(base2_closure)}")
                print(f"    Target closure: {sorted(target_closure)}")
                print(f"    Combined covers target: {target_closure.issubset(combined)}")
                
                # Check generator coverage
                print(f"    Target generators: {len(target.min_generators)}")
                for j, gen in enumerate(target.min_generators):
                    gen_species = set(species_list_to_names(gen))
                    in_combined = gen_species.issubset(combined)
                    in_base1 = gen_species.issubset(base1_closure)
                    in_base2 = gen_species.issubset(base2_closure)
                    print(f"      Gen{j+1}: {sorted(gen_species)} - Combined:{in_combined}, Base1:{in_base1}, Base2:{in_base2}")
            
            if len(bf_only) > 5:
                print(f"\n  ... and {len(bf_only) - 5} more")
        
        # Investigate synergies only in hierarchy-aware
        if ha_only:
            print(f"\n🔍 Synergies only in Hierarchy-Aware ({len(ha_only)}):")
            for i, (base1_label, base2_label, target_label) in enumerate(sorted(ha_only)[:5]):
                print(f"\n  [{i+1}] {base1_label} + {base2_label} → {target_label}")
                
                # This case is less likely but could happen due to algorithm differences
                base1 = get_erc(base1_label)
                base2 = get_erc(base2_label)
                target = get_erc(target_label)
                
                print(f"    Base1: {base1_label} - {species_list_to_names(base1.get_closure(self.rn))}")
                print(f"    Base2: {base2_label} - {species_list_to_names(base2.get_closure(self.rn))}")
                print(f"    Target: {target_label} - {species_list_to_names(target.get_closure(self.rn))}")
            
            if len(ha_only) > 5:
                print(f"\n  ... and {len(ha_only) - 5} more")
    
    def generate_detailed_report(self, output_file=None):
        """Generate a detailed comparison report."""
        if not self.comparison_results:
            print("No comparison results available. Run comparison first.")
            return
        
        report_lines = []
        report_lines.append("="*80)
        report_lines.append("SYNERGY ALGORITHM COMPARISON REPORT")
        report_lines.append("="*80)
        
        # Summary statistics
        report_lines.append("\nSUMMARY STATISTICS")
        report_lines.append("-"*40)
        
        results = self.comparison_results
        report_lines.append(f"Brute Force (all): {len(results['bf_all'])} synergies")
        report_lines.append(f"Brute Force (minimal): {len(results['bf_minimal'])} synergies")
        report_lines.append(f"Hierarchy-Aware: {len(results['ha_synergies'])} synergies")
        report_lines.append(f"Time speedup: {results['time_speedup']:.2f}x")
        
        # Agreement metrics
        report_lines.append("\nAGREEMENT METRICS")
        report_lines.append("-"*40)
        
        if len(results['bf_all']) > 0:
            agreement_all = len(results['common_with_all']) / len(results['bf_all']) * 100
            report_lines.append(f"Agreement with BF(all): {agreement_all:.1f}%")
        
        if len(results['bf_minimal']) > 0:
            common_with_minimal = results['bf_minimal'] & results['ha_synergies']
            agreement_minimal = len(common_with_minimal) / len(results['bf_minimal']) * 100
            report_lines.append(f"Agreement with BF(minimal): {agreement_minimal:.1f}%")
        
        # Differences
        report_lines.append("\nDIFFERENCES")
        report_lines.append("-"*40)
        report_lines.append(f"Only in Brute Force: {len(results['bf_only_all'])}")
        report_lines.append(f"Only in Hierarchy-Aware: {len(results['ha_only_vs_all'])}")
        
        # List all differences
        if results['bf_only_all']:
            report_lines.append("\nSynergies only in Brute Force:")
            for base1, base2, target in sorted(results['bf_only_all']):
                report_lines.append(f"  {base1} + {base2} → {target}")
        
        if results['ha_only_vs_all']:
            report_lines.append("\nSynergies only in Hierarchy-Aware:")
            for base1, base2, target in sorted(results['ha_only_vs_all']):
                report_lines.append(f"  {base1} + {base2} → {target}")
        
        # Print or save report
        report_text = "\n".join(report_lines)
        
        if output_file:
            with open(output_file, 'w') as f:
                f.write(report_text)
            print(f"\nDetailed report saved to: {output_file}")
        else:
            print("\n" + report_text)
        
        return report_text
    
    def export_comparison_data(self, output_dir="comparison_results"):
        """Export comparison data to CSV files for further analysis."""
        Path(output_dir).mkdir(exist_ok=True)
        
        # Export all synergies found
        all_synergies = []
        
        # Add brute force synergies
        for synergy in self.brute_force_results['all']:
            base1 = synergy['base1'].label
            base2 = synergy['base2'].label
            target = synergy['target'].label
            all_synergies.append({
                'base1': base1,
                'base2': base2,
                'target': target,
                'algorithm': 'brute_force',
                'is_minimal': synergy in self.brute_force_results['minimal']
            })
        
        # Add hierarchy-aware synergies
        for target_erc, synergy_list in self.hierarchy_aware_results.items():
            for synergy in synergy_list:
                all_synergies.append({
                    'base1': synergy.base1_erc.label,
                    'base2': synergy.base2_erc.label,
                    'target': target_erc.label,
                    'algorithm': 'hierarchy_aware',
                    'is_minimal': True  # Hierarchy-aware only finds minimal
                })
        
        # Create DataFrame and save
        df = pd.DataFrame(all_synergies)
        df.to_csv(f"{output_dir}/all_synergies.csv", index=False)
        
        # Export comparison summary
        summary = {
            'metric': ['BF_all', 'BF_minimal', 'HA_synergies', 'Common', 'BF_only', 'HA_only', 'Time_speedup'],
            'value': [
                len(self.comparison_results['bf_all']),
                len(self.comparison_results['bf_minimal']),
                len(self.comparison_results['ha_synergies']),
                len(self.comparison_results['common_with_all']),
                len(self.comparison_results['bf_only_all']),
                len(self.comparison_results['ha_only_vs_all']),
                self.comparison_results['time_speedup']
            ]
        }
        
        summary_df = pd.DataFrame(summary)
        summary_df.to_csv(f"{output_dir}/comparison_summary.csv", index=False)
        
        print(f"\nComparison data exported to: {output_dir}/")


def run_comprehensive_comparison(reaction_network_file: str, 
                                verbose: bool = False,
                                check_minimal: bool = True,
                                export_results: bool = True):
    """
    Run a comprehensive comparison of synergy algorithms.
    
    Parameters
    ----------
    reaction_network_file : str
        Path to the reaction network file
    verbose : bool
        Print detailed output during algorithm execution
    check_minimal : bool
        Also compute minimal synergies for brute force
    export_results : bool
        Export results to CSV files
    """
    print("🚀 Starting Synergy Algorithm Comparison")
    print("="*80)
    
    # Initialize comparator
    comparator = SynergyAlgorithmComparator(reaction_network_file)
    
    # Run brute force algorithm
    bf_results = comparator.run_brute_force_algorithm(minimal=check_minimal, verbose=verbose)
    
    # Run hierarchy-aware algorithm
    ha_results = comparator.run_hierarchy_aware_algorithm(verbose=verbose)
    
    # Compare results
    comparison = comparator.compare_results(bf_results, ha_results)
    
    # Generate detailed report
    comparator.generate_detailed_report()
    
    # Export results if requested
    if export_results:
        comparator.export_comparison_data()
    
    print("\n✅ Comparison Complete!")
    
    return comparator


if __name__ == "__main__":
    # Example usage
    reaction_file = "ERCs_test2.txt"  # Update with your file path
    
    # Run comparison
    comparator = run_comprehensive_comparison(
        reaction_file,
        verbose=False,  # Set to True for detailed output
        check_minimal=True,
        export_results=True
    )
    
    # Access results programmatically if needed
    # print(comparator.comparison_results)