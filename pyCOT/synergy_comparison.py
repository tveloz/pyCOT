#!/usr/bin/env python3
"""
Synergy Algorithm Comparison Script

Compares the Brute Force and Hierarchy Aware synergy algorithms on:
1. Execution time
2. Number of synergies found
3. Exact synergy matches
4. Detailed analysis of differences

@author: ERC Analysis Framework
"""

import time
import sys
import os
from pathlib import Path
from collections import defaultdict, Counter
from dataclasses import dataclass
from typing import List, Set, Tuple, Dict, Optional

# Import the required modules
from pyCOT.io.functions import read_txt, from_string
from pyCOT.ERC_Hierarchy import ERC
from pyCOT.Synergy_Brute_Force_Algorithm import BruteForceSynergyCalculator
from pyCOT.Synergy_Hierarchy_Aware_Algorithm import DynamicSynergyExplorer


@dataclass
class SynergyTriple:
    """Normalized synergy representation for comparison."""
    base1: str
    base2: str
    target: str
    
    def __post_init__(self):
        # Ensure consistent ordering for comparison
        if self.base1 > self.base2:
            self.base1, self.base2 = self.base2, self.base1
    
    def __hash__(self):
        return hash((self.base1, self.base2, self.target))
    
    def __eq__(self, other):
        return (self.base1, self.base2, self.target) == (other.base1, other.base2, other.target)
    
    def __str__(self):
        return f"{self.base1} + {self.base2} → {self.target}"


@dataclass
class AlgorithmResult:
    """Results from running a synergy algorithm."""
    algorithm_name: str
    execution_time: float
    synergies: List[SynergyTriple]
    raw_results: List  # Original algorithm output
    error: Optional[str] = None
    
    @property
    def synergy_count(self):
        return len(self.synergies)
    
    @property
    def synergy_set(self):
        return set(self.synergies)


@dataclass
class ComparisonResults:
    """Results of comparing two algorithms."""
    brute_force: AlgorithmResult
    hierarchy_aware: AlgorithmResult
    common_synergies: Set[SynergyTriple]
    brute_force_only: Set[SynergyTriple]
    hierarchy_aware_only: Set[SynergyTriple]
    
    @property
    def perfect_match(self):
        return len(self.brute_force_only) == 0 and len(self.hierarchy_aware_only) == 0
    
    @property
    def speedup_factor(self):
        if self.hierarchy_aware.execution_time > 0:
            return self.brute_force.execution_time / self.hierarchy_aware.execution_time
        return float('inf')


class SynergyAlgorithmComparator:
    """Compares synergy algorithms on timing and accuracy."""
    
    def __init__(self, verbose=True):
        self.verbose = verbose
        self.comparison_history = []
    
    def load_reaction_network(self, source):
        """
        Load reaction network from various sources.
        
        Parameters
        ----------
        source : str or ReactionNetwork
            Path to file, reaction string, or ReactionNetwork object
        """
        if hasattr(source, 'reactions'):  # Already a ReactionNetwork
            return source
        elif isinstance(source, str):
            if os.path.exists(source):
                # File path
                if source.endswith('.txt'):
                    return read_txt(source)
                else:
                    raise ValueError(f"Unsupported file format: {source}")
            else:
                # Treat as reaction string
                return from_string(source)
        else:
            raise ValueError(f"Unsupported source type: {type(source)}")
    
    def normalize_brute_force_results(self, raw_results):
        """Convert brute force results to normalized synergy triples."""
        synergies = []
        for result in raw_results:
            triple = SynergyTriple(
                base1=result['base1'].label,
                base2=result['base2'].label,
                target=result['target'].label
            )
            synergies.append(triple)
        return synergies
    
    def normalize_hierarchy_aware_results(self, raw_results):
        """Convert hierarchy aware results to normalized synergy triples."""
        synergies = []
        for result in raw_results:
            triple = SynergyTriple(
                base1=result.base1_erc.label,
                base2=result.base2_erc.label,
                target=result.target_erc.label
            )
            synergies.append(triple)
        return synergies
    
    def run_brute_force_algorithm(self, rn, ercs, hierarchy_graph):
        """Run the brute force algorithm and measure performance."""
        if self.verbose:
            print("🔄 Running Brute Force Algorithm...")
        
        try:
            start_time = time.perf_counter()
            
            calculator = BruteForceSynergyCalculator(rn, ercs, hierarchy_graph)
            raw_results = calculator.find_all_synergies(verbose=False)
            
            end_time = time.perf_counter()
            execution_time = end_time - start_time
            
            synergies = self.normalize_brute_force_results(raw_results)
            
            return AlgorithmResult(
                algorithm_name="Brute Force",
                execution_time=execution_time,
                synergies=synergies,
                raw_results=raw_results
            )
            
        except Exception as e:
            return AlgorithmResult(
                algorithm_name="Brute Force",
                execution_time=0.0,
                synergies=[],
                raw_results=[],
                error=str(e)
            )
    
    def run_hierarchy_aware_algorithm(self, rn, ercs, hierarchy_graph):
        """Run the hierarchy aware algorithm and measure performance."""
        if self.verbose:
            print("🔄 Running Hierarchy Aware Algorithm...")
        
        try:
            start_time = time.perf_counter()
            
            explorer = DynamicSynergyExplorer(rn, ercs, hierarchy_graph)
            raw_results = explorer.run_complete_analysis(verbose=False)
            
            end_time = time.perf_counter()
            execution_time = end_time - start_time
            
            synergies = self.normalize_hierarchy_aware_results(raw_results)
            
            return AlgorithmResult(
                algorithm_name="Hierarchy Aware",
                execution_time=execution_time,
                synergies=synergies,
                raw_results=raw_results
            )
            
        except Exception as e:
            return AlgorithmResult(
                algorithm_name="Hierarchy Aware",
                execution_time=0.0,
                synergies=[],
                raw_results=[],
                error=str(e)
            )
    
    def compare_algorithms(self, reaction_network_source, test_name="Test"):
        """
        Compare both algorithms on a given reaction network.
        
        Parameters
        ----------
        reaction_network_source : str or ReactionNetwork
            Source for the reaction network
        test_name : str
            Name for this test case
            
        Returns
        -------
        ComparisonResults
            Detailed comparison results
        """
        if self.verbose:
            print(f"\n{'='*80}")
            print(f"🧪 SYNERGY ALGORITHM COMPARISON: {test_name}")
            print(f"{'='*80}")
        
        # Load reaction network
        try:
            rn = self.load_reaction_network(reaction_network_source)
            if self.verbose:
                print(f"✅ Loaded reaction network: {len(rn.species())} species, {len(rn.reactions())} reactions")
        except Exception as e:
            if self.verbose:
                print(f"❌ Failed to load reaction network: {e}")
            return None
        
        # Generate ERCs and hierarchy
        try:
            if self.verbose:
                print("🔧 Computing ERCs and hierarchy...")
            start_time = time.perf_counter()
            ercs = ERC.ERCs(rn)
            hierarchy_graph = ERC.build_hierarchy_graph(ercs, rn)
            setup_time = time.perf_counter() - start_time
            
            if self.verbose:
                print(f"✅ Generated {len(ercs)} ERCs in {setup_time:.3f}s")
                print(f"   Hierarchy: {hierarchy_graph.number_of_nodes()} nodes, {hierarchy_graph.number_of_edges()} edges")
        except Exception as e:
            if self.verbose:
                print(f"❌ Failed to compute ERCs: {e}")
            return None
        
        # Run both algorithms
        brute_force_result = self.run_brute_force_algorithm(rn, ercs, hierarchy_graph)
        hierarchy_aware_result = self.run_hierarchy_aware_algorithm(rn, ercs, hierarchy_graph)
        
        # Analyze results
        if brute_force_result.error or hierarchy_aware_result.error:
            if self.verbose:
                if brute_force_result.error:
                    print(f"❌ Brute Force Error: {brute_force_result.error}")
                if hierarchy_aware_result.error:
                    print(f"❌ Hierarchy Aware Error: {hierarchy_aware_result.error}")
            return None
        
        # Compare synergies
        brute_force_set = brute_force_result.synergy_set
        hierarchy_aware_set = hierarchy_aware_result.synergy_set
        
        common_synergies = brute_force_set.intersection(hierarchy_aware_set)
        brute_force_only = brute_force_set - hierarchy_aware_set
        hierarchy_aware_only = hierarchy_aware_set - brute_force_set
        
        comparison = ComparisonResults(
            brute_force=brute_force_result,
            hierarchy_aware=hierarchy_aware_result,
            common_synergies=common_synergies,
            brute_force_only=brute_force_only,
            hierarchy_aware_only=hierarchy_aware_only
        )
        
        self.comparison_history.append((test_name, comparison))
        
        if self.verbose:
            self.print_comparison_results(comparison)
        
        return comparison
    
    def print_comparison_results(self, comparison):
        """Print detailed comparison results."""
        bf = comparison.brute_force
        ha = comparison.hierarchy_aware
        
        print(f"\n📊 TIMING COMPARISON")
        print(f"{'='*50}")
        print(f"Brute Force:     {bf.execution_time:.4f}s")
        print(f"Hierarchy Aware: {ha.execution_time:.4f}s")
        if ha.execution_time > 0:
            print(f"Speedup Factor:  {comparison.speedup_factor:.2f}x")
        
        print(f"\n🎯 SYNERGY COUNT COMPARISON")
        print(f"{'='*50}")
        print(f"Brute Force:     {bf.synergy_count} synergies")
        print(f"Hierarchy Aware: {ha.synergy_count} synergies")
        print(f"Common:          {len(comparison.common_synergies)} synergies")
        
        print(f"\n🔍 RESULT ACCURACY")
        print(f"{'='*50}")
        if comparison.perfect_match:
            print("✅ PERFECT MATCH - Both algorithms found identical synergies!")
        else:
            print("❌ MISMATCH DETECTED")
            
            if comparison.brute_force_only:
                print(f"\n🔸 Synergies found ONLY by Brute Force ({len(comparison.brute_force_only)}):")
                for synergy in sorted(comparison.brute_force_only, key=str):
                    print(f"   {synergy}")
            
            if comparison.hierarchy_aware_only:
                print(f"\n🔹 Synergies found ONLY by Hierarchy Aware ({len(comparison.hierarchy_aware_only)}):")
                for synergy in sorted(comparison.hierarchy_aware_only, key=str):
                    print(f"   {synergy}")
    
    def print_detailed_investigation(self, comparison, rn, ercs):
        """Print detailed information for investigating differences."""
        if comparison.perfect_match:
            print("\n✅ No differences to investigate - perfect match!")
            return
        
        print(f"\n🔬 DETAILED INVESTIGATION")
        print(f"{'='*80}")
        
        # Analyze ERCs involved in mismatches
        all_mismatched = comparison.brute_force_only.union(comparison.hierarchy_aware_only)
        involved_ercs = set()
        for synergy in all_mismatched:
            involved_ercs.update([synergy.base1, synergy.base2, synergy.target])
        
        print(f"\nERCs involved in mismatches: {sorted(involved_ercs)}")
        
        # Print ERC details
        erc_dict = {erc.label: erc for erc in ercs}
        for erc_label in sorted(involved_ercs):
            if erc_label in erc_dict:
                erc = erc_dict[erc_label]
                closure = erc.get_closure(rn)
                from pyCOT.ERC_Hierarchy import species_list_to_names
                closure_names = species_list_to_names(closure)
                print(f"\n{erc_label}:")
                print(f"  Closure: {closure_names}")
                print(f"  Generators: {[species_list_to_names(gen) for gen in erc.min_generators]}")
        
        # Analyze specific mismatches
        if comparison.brute_force_only:
            print(f"\n🔸 ANALYZING BRUTE FORCE EXCLUSIVE SYNERGIES:")
            self.analyze_synergy_validity(comparison.brute_force_only, rn, erc_dict, "Brute Force")
        
        if comparison.hierarchy_aware_only:
            print(f"\n🔹 ANALYZING HIERARCHY AWARE EXCLUSIVE SYNERGIES:")
            self.analyze_synergy_validity(comparison.hierarchy_aware_only, rn, erc_dict, "Hierarchy Aware")
    
    def analyze_synergy_validity(self, synergies, rn, erc_dict, algorithm_name):
        """Analyze why certain synergies are found by one algorithm but not the other."""
        for synergy in sorted(synergies, key=str):
            print(f"\n   Analyzing: {synergy}")
            
            try:
                base1_erc = erc_dict[synergy.base1]
                base2_erc = erc_dict[synergy.base2]
                target_erc = erc_dict[synergy.target]
                
                # Check synergy validity with detailed output
                print(f"     Checking synergy validity...")
                
                # Use brute force validator as reference
                bf_calculator = BruteForceSynergyCalculator(rn)
                is_valid, details = bf_calculator.validates_synergy(base1_erc, base2_erc, target_erc, verbose=True)
                
                print(f"     Valid synergy: {is_valid}")
                
                if is_valid:
                    print(f"     Coverage ratio: {details.get('coverage_ratio', 'N/A'):.2%}")
                    print(f"     Synergy ratio: {details.get('synergy_ratio', 'N/A'):.2%}")
                
            except Exception as e:
                print(f"     Error analyzing synergy: {e}")
    
    def run_test_suite(self, test_cases):
        """
        Run comparison on multiple test cases.
        
        Parameters
        ----------
        test_cases : dict
            Dictionary mapping test names to reaction network sources
        """
        print(f"\n🧪 RUNNING SYNERGY ALGORITHM TEST SUITE")
        print(f"{'='*80}")
        print(f"Test cases: {len(test_cases)}")
        
        results_summary = []
        
        for test_name, rn_source in test_cases.items():
            print(f"\n{'─'*40}")
            print(f"🔬 Test Case: {test_name}")
            print(f"{'─'*40}")
            
            comparison = self.compare_algorithms(rn_source, test_name)
            
            if comparison:
                results_summary.append({
                    'test_name': test_name,
                    'bf_time': comparison.brute_force.execution_time,
                    'ha_time': comparison.hierarchy_aware.execution_time,
                    'speedup': comparison.speedup_factor,
                    'bf_count': comparison.brute_force.synergy_count,
                    'ha_count': comparison.hierarchy_aware.synergy_count,
                    'perfect_match': comparison.perfect_match,
                    'mismatches': len(comparison.brute_force_only) + len(comparison.hierarchy_aware_only)
                })
            else:
                results_summary.append({
                    'test_name': test_name,
                    'error': True
                })
        
        # Print summary
        self.print_test_suite_summary(results_summary)
        
        return results_summary
    
    def print_test_suite_summary(self, results_summary):
        """Print summary of test suite results."""
        print(f"\n📊 TEST SUITE SUMMARY")
        print(f"{'='*100}")
        
        valid_results = [r for r in results_summary if not r.get('error', False)]
        
        if not valid_results:
            print("❌ No valid test results")
            return
        
        print(f"{'Test Name':<20} {'BF Time':<10} {'HA Time':<10} {'Speedup':<10} {'BF Count':<10} {'HA Count':<10} {'Match':<8} {'Diffs'}")
        print("─" * 100)
        
        total_bf_time = 0
        total_ha_time = 0
        perfect_matches = 0
        
        for result in valid_results:
            total_bf_time += result['bf_time']
            total_ha_time += result['ha_time']
            if result['perfect_match']:
                perfect_matches += 1
            
            match_symbol = "✅" if result['perfect_match'] else "❌"
            
            print(f"{result['test_name']:<20} "
                  f"{result['bf_time']:<10.4f} "
                  f"{result['ha_time']:<10.4f} "
                  f"{result['speedup']:<10.2f} "
                  f"{result['bf_count']:<10} "
                  f"{result['ha_count']:<10} "
                  f"{match_symbol:<8} "
                  f"{result['mismatches']}")
        
        print("─" * 100)
        print(f"{'TOTALS':<20} "
              f"{total_bf_time:<10.4f} "
              f"{total_ha_time:<10.4f} "
              f"{total_bf_time/total_ha_time if total_ha_time > 0 else float('inf'):<10.2f} "
              f"{'─':<10} "
              f"{'─':<10} "
              f"{perfect_matches}/{len(valid_results):<8} "
              f"─")
        
        print(f"\n🎯 OVERALL STATISTICS:")
        print(f"   Perfect matches: {perfect_matches}/{len(valid_results)} ({perfect_matches/len(valid_results)*100:.1f}%)")
        if total_ha_time > 0:
            print(f"   Average speedup: {total_bf_time/total_ha_time:.2f}x")
        print(f"   Total execution time: {total_bf_time + total_ha_time:.4f}s")


def main():
    """Main function for running comparisons."""
    
    # Example test cases
    test_cases = {
        "Simple Test": """
        R1: a1+b1=>a2+b2;
        R2: a2+b1=>a3+b3;
        R3: a3+b3=>x3+a2;
        R4: a1+b4=>a5+b5;
        R5: a5+b5=>x5+b2;
        R6: a6+b1=>a7+b7;
        """,

        "File Test": "networks\RandomAlife\RN_Ns_20_Norg_6_id_137.txt"  # If the file exists
    }
    
    # Create comparator
    comparator = SynergyAlgorithmComparator(verbose=True)
    
    # Run single test
    print("🚀 Running single test comparison...")
    comparison = comparator.compare_algorithms(test_cases["File Test"], "File Test")

    # If you want detailed investigation of differences
    if comparison and not comparison.perfect_match:
        print("\n🔬 Running detailed investigation...")
        rn = comparator.load_reaction_network(test_cases["File Test"])
        ercs = ERC.ERCs(rn)
        comparator.print_detailed_investigation(comparison, rn, ercs)
    
    # Run test suite (comment out if not needed)
    # print("\n" + "="*80)
    # comparator.run_test_suite(test_cases)


if __name__ == "__main__":
    main()