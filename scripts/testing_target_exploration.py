#!/usr/bin/env python3
"""
ERC Exploration Testing Script

This script provides comprehensive testing of the ERCExplorer class,
demonstrating all functionality from basic target identification to
advanced synergy analysis and filtering.

Usage:
    python test_erc_exploration.py
    
Or with custom file:
    python test_erc_exploration.py --file your_network.txt
"""

import argparse
import sys
from pathlib import Path
import time

# Import necessary modules
from pyCOT.io.functions import read_txt
from pyCOT.ERC_Hierarchy import ERC, species_list_to_names
from pyCOT.ERC_Exploration import ERCExplorer


class ERCExplorationTester:
    """Comprehensive testing class for ERC exploration functionality."""
    
    def __init__(self, network_file):
        """Initialize the tester with a reaction network file."""
        self.network_file = network_file
        self.RN = None
        self.explorer = None
        
    def load_network(self):
        """Load the reaction network and create the explorer."""
        print(f"🔄 Loading reaction network from: {self.network_file}")
        try:
            self.RN = read_txt(self.network_file)
            print(f"✅ Network loaded successfully!")
            print(f"   Species: {len(self.RN.species())}")
            print(f"   Reactions: {len(self.RN.reactions())}")
            
            # Show some sample reactions
            print(f"\n📋 Sample reactions:")
            for i, reaction in enumerate(self.RN.reactions()[:5]):
                print(f"   {reaction}")
            if len(self.RN.reactions()) > 5:
                print(f"   ... and {len(self.RN.reactions()) - 5} more")
                
        except FileNotFoundError:
            print(f"❌ Error: File not found at {self.network_file}")
            sys.exit(1)
        except Exception as e:
            print(f"❌ Error loading network: {e}")
            sys.exit(1)
    
    def create_explorer(self):
        """Create the ERC explorer instance."""
        print(f"\n🔄 Creating ERC Explorer...")
        start_time = time.time()
        
        try:
            self.explorer = ERCExplorer(self.RN)
            creation_time = time.time() - start_time
            
            print(f"✅ Explorer created successfully! ({creation_time:.2f}s)")
            print(f"   ERCs generated: {len(self.explorer.ercs)}")
            print(f"   Hierarchy edges: {self.explorer.hierarchy.number_of_edges()}")
            
            # Show ERC details
            print(f"\n📋 ERC Overview:")
            for erc in self.explorer.ercs:
                closure = species_list_to_names(erc.get_closure(self.RN))
                min_gens = [species_list_to_names(gen) for gen in erc.min_generators]
                print(f"   {erc.label}: {len(closure)} species, {len(min_gens)} generators")
                
        except Exception as e:
            print(f"❌ Error creating explorer: {e}")
            import traceback
            traceback.print_exc()
            sys.exit(1)
    
    def test_hierarchy_analysis(self):
        """Test the hierarchy structural analysis."""
        print(f"\n" + "="*60)
        print(f"🏗️  TESTING HIERARCHY ANALYSIS")
        print(f"="*60)
        
        start_time = time.time()
        hierarchy_analysis = self.explorer.analyze_hierarchy_structure()
        analysis_time = time.time() - start_time
        
        print(f"✅ Hierarchy analysis completed! ({analysis_time:.2f}s)")
        
        # Display key metrics
        basic = hierarchy_analysis['basic_metrics']
        print(f"\n📊 Key Metrics:")
        print(f"   Hierarchy depth: {basic['hierarchy_depth']} levels")
        print(f"   Total containment relations: {basic['total_edges']}")
        print(f"   Level range: {basic['min_level']} to {basic['max_level']}")
        
        # Level distribution
        print(f"\n📈 Level Distribution:")
        for level in sorted(hierarchy_analysis['level_distribution'].keys()):
            count = hierarchy_analysis['level_distribution'][level]
            bar = "█" * min(count, 20)
            print(f"   Level {level:2}: {bar} {count}")
        
        # Root and leaf analysis
        root_leaf = hierarchy_analysis['root_leaf_analysis']
        print(f"\n🌳 Topology:")
        print(f"   Root ERCs: {root_leaf['num_roots']} {root_leaf['root_nodes']}")
        print(f"   Leaf ERCs: {root_leaf['num_leaves']} {root_leaf['leaf_nodes']}")
        
        return hierarchy_analysis
    
    def test_target_identification(self):
        """Test the target identification algorithm."""
        print(f"\n" + "="*60)
        print(f"🎯 TESTING TARGET IDENTIFICATION")
        print(f"="*60)
        
        start_time = time.time()
        potential_synergies = self.explorer.identify_potential_synergy_targets()
        identification_time = time.time() - start_time
        
        print(f"✅ Target identification completed! ({identification_time:.2f}s)")
        
        # Summary statistics
        base1_count = len(potential_synergies)
        total_pairs = sum(len(targets) for targets in potential_synergies.values())
        
        print(f"\n📊 Identification Results:")
        print(f"   Base1 ERCs with targets: {base1_count}")
        print(f"   Total potential synergy pairs: {total_pairs}")
        
        if base1_count > 0:
            avg_targets = total_pairs / base1_count
            print(f"   Average targets per base1: {avg_targets:.1f}")
            
            # Show distribution
            target_counts = [len(targets) for targets in potential_synergies.values()]
            print(f"   Target count range: {min(target_counts)} - {max(target_counts)}")
        
        # Show sample results
        if potential_synergies:
            print(f"\n📋 Sample Results:")
            for i, (base1, targets) in enumerate(list(potential_synergies.items())[:3]):
                print(f"   {base1.label} → {[t.label for t in targets]}")
            if len(potential_synergies) > 3:
                print(f"   ... and {len(potential_synergies) - 3} more base1 ERCs")
        
        return potential_synergies
    
    def test_target_centric_organization(self, potential_synergies):
        """Test the target-centric organization and filtering."""
        print(f"\n" + "="*60)
        print(f"🎛️  TESTING TARGET-CENTRIC ORGANIZATION")
        print(f"="*60)
        
        start_time = time.time()
        target_synergies = self.explorer.reorganize_to_target_centric(potential_synergies)
        organization_time = time.time() - start_time
        
        print(f"✅ Target-centric organization completed! ({organization_time:.2f}s)")
        
        # Summary statistics
        target_count = len(target_synergies)
        total_filtered_pairs = sum(len(base1_infos) for base1_infos in target_synergies.values())
        original_pairs = sum(len(targets) for targets in potential_synergies.values())
        
        print(f"\n📊 Organization Results:")
        print(f"   Targets with synergies: {target_count}")
        print(f"   Filtered synergy pairs: {total_filtered_pairs}")
        print(f"   Reduction from original: {original_pairs - total_filtered_pairs} pairs")
        if original_pairs > 0:
            reduction_pct = ((original_pairs - total_filtered_pairs) / original_pairs) * 100
            print(f"   Efficiency gain: {reduction_pct:.1f}% reduction")
        
        # Show base1 distribution per target
        if target_synergies:
            base1_counts = [len(base1_infos) for base1_infos in target_synergies.values()]
            print(f"   Base1s per target: {min(base1_counts)} - {max(base1_counts)} (avg: {sum(base1_counts)/len(base1_counts):.1f})")
        
        # Show sample filtered results
        if target_synergies:
            print(f"\n📋 Sample Filtered Results:")
            for i, (target, base1_infos) in enumerate(list(target_synergies.items())[:3]):
                base1_labels = [info.base1_erc.label for info in base1_infos]
                overlap_counts = [len(info.overlaps) for info in base1_infos]
                print(f"   {target.label} ← {base1_labels} (overlaps: {overlap_counts})")
            if len(target_synergies) > 3:
                print(f"   ... and {len(target_synergies) - 3} more targets")
        
        return target_synergies
    
    def test_detailed_overlap_analysis(self, target_synergies):
        """Test detailed overlap information analysis."""
        print(f"\n" + "="*60)
        print(f"🔍 TESTING DETAILED OVERLAP ANALYSIS")
        print(f"="*60)
        
        if not target_synergies:
            print(f"⚠️  No target synergies to analyze")
            return
        
        # Analyze all targets
        for sample_target, sample_base1_infos in target_synergies.items():
            print(f"\n{'='*50}")
        
            print(f"📋 Detailed Analysis for Target: {sample_target.label}")
            print(f"   Closure: {species_list_to_names(sample_target.get_closure(self.RN))}")
            print(f"   Generators: {[species_list_to_names(gen) for gen in sample_target.min_generators]}")
            
            for i, base1_info in enumerate(sample_base1_infos):
                print(f"\n   Base1 [{i+1}]: {base1_info.base1_erc.label} (Level {base1_info.hierarchy_level})")
                print(f"      Closure: {species_list_to_names(base1_info.base1_erc.get_closure(self.RN))}")
                print(f"      Overlap signature: {len(base1_info.overlap_signature)} elements")
                
                for j, overlap in enumerate(base1_info.overlaps):
                    print(f"      Overlap {j+1}:")
                    print(f"         Generator {overlap.generator_index + 1}: {sorted(overlap.generator_species)}")
                    print(f"         Overlapping: {sorted(overlap.overlapping_species)}")
                    print(f"         Base1 unique: {sorted(overlap.base1_unique) if overlap.base1_unique else 'None'}")
                    print(f"         Gen unique: {sorted(overlap.generator_unique) if overlap.generator_unique else 'None'}")
    
    def test_filtering_effectiveness(self, potential_synergies, target_synergies):
        """Test and analyze filtering effectiveness."""
        print(f"\n" + "="*60)
        print(f"⚡ TESTING FILTERING EFFECTIVENESS")
        print(f"="*60)
        
        # Analyze signature diversity
        all_signatures = set()
        signature_frequencies = {}
        
        for target, base1_infos in target_synergies.items():
            for base1_info in base1_infos:
                sig = base1_info.overlap_signature
                all_signatures.add(sig)
                signature_frequencies[sig] = signature_frequencies.get(sig, 0) + 1
        
        print(f"📊 Signature Analysis:")
        print(f"   Unique overlap signatures: {len(all_signatures)}")
        print(f"   Average signature frequency: {sum(signature_frequencies.values()) / len(signature_frequencies):.2f}")
        
        # Show most common signatures
        sorted_sigs = sorted(signature_frequencies.items(), key=lambda x: x[1], reverse=True)
        print(f"   Most frequent signatures:")
        for i, (sig, freq) in enumerate(sorted_sigs[:3]):
            sig_desc = f"({len(sig)} overlaps)" if sig else "(empty)"
            print(f"      {i+1}. {sig_desc}: {freq} times")
        
        # Level-based filtering analysis
        level_stats = {}
        for target, base1_infos in target_synergies.items():
            levels = [info.hierarchy_level for info in base1_infos]
            level_stats[target.label] = {
                'min_level': min(levels),
                'max_level': max(levels),
                'level_diversity': len(set(levels))
            }
        
        all_min_levels = [stats['min_level'] for stats in level_stats.values()]
        all_max_levels = [stats['max_level'] for stats in level_stats.values()]
        all_diversities = [stats['level_diversity'] for stats in level_stats.values()]
        
        print(f"\n📊 Level-based Filtering:")
        print(f"   Level range in results: {min(all_min_levels)} - {max(all_max_levels)}")
        print(f"   Average level diversity per target: {sum(all_diversities) / len(all_diversities):.1f}")
        print(f"   Targets with single level: {sum(1 for d in all_diversities if d == 1)}/{len(all_diversities)}")
    
    def run_comprehensive_test(self):
        """Run all tests in sequence."""
        print(f"\n" + "="*80)
        print(f"🚀 STARTING COMPREHENSIVE ERC EXPLORATION TEST")
        print(f"="*80)
        
        total_start_time = time.time()
        
        # Step 1: Load network and create explorer
        self.load_network()
        self.create_explorer()
        
        # Step 2: Test hierarchy analysis
        hierarchy_analysis = self.test_hierarchy_analysis()
        
        # Step 3: Test target identification
        potential_synergies = self.test_target_identification()
        
        # Step 4: Test target-centric organization
        target_synergies = self.test_target_centric_organization(potential_synergies)
        
        # Step 5: Test detailed overlap analysis
        self.test_detailed_overlap_analysis(target_synergies)
        
        # Step 6: Test filtering effectiveness
        self.test_filtering_effectiveness(potential_synergies, target_synergies)
        
        # Generate comprehensive report
        print(f"\n" + "="*60)
        print(f"📊 GENERATING COMPREHENSIVE REPORT")
        print(f"="*60)
        
        report_results = self.explorer.generate_comprehensive_report()
        
        # Final summary
        total_time = time.time() - total_start_time
        print(f"\n" + "="*80)
        print(f"✅ COMPREHENSIVE TEST COMPLETED SUCCESSFULLY!")
        print(f"="*80)
        print(f"Total execution time: {total_time:.2f} seconds")
        print(f"Network: {Path(self.network_file).name}")
        print(f"ERCs analyzed: {len(self.explorer.ercs)}")
        print(f"Potential synergy pairs found: {sum(len(targets) for targets in potential_synergies.values())}")
        print(f"Filtered synergy relationships: {sum(len(base1_infos) for base1_infos in target_synergies.values())}")
        
        return {
            'hierarchy_analysis': hierarchy_analysis,
            'potential_synergies': potential_synergies,
            'target_synergies': target_synergies,
            'execution_time': total_time
        }


def test_multiple_networks():
    """Test the explorer with multiple different networks."""
    test_networks = [
        'networks/testing/ERCs_test2.txt',
        'networks/testing/Farm.txt',
        # Add more test networks as available
    ]
    
    print(f"\n" + "="*80)
    print(f"🌐 TESTING MULTIPLE NETWORKS")
    print(f"="*80)
    
    results = {}
    
    for network_file in test_networks:
        if Path(network_file).exists():
            print(f"\n🔄 Testing network: {network_file}")
            try:
                tester = ERCExplorationTester(network_file)
                network_results = tester.run_comprehensive_test()
                results[network_file] = network_results
                print(f"✅ {network_file} completed successfully")
            except Exception as e:
                print(f"❌ {network_file} failed: {e}")
                results[network_file] = {'error': str(e)}
        else:
            print(f"⚠️  Network file not found: {network_file}")
    
    # Summary comparison
    if len(results) > 1:
        print(f"\n" + "="*60)
        print(f"📊 MULTI-NETWORK COMPARISON")
        print(f"="*60)
        
        for network, result in results.items():
            if 'error' not in result:
                potential_count = len(result['potential_synergies'])
                target_count = len(result['target_synergies'])
                time_taken = result['execution_time']
                print(f"{Path(network).name}: {potential_count} base1s → {target_count} targets ({time_taken:.1f}s)")
    
    return results


def main():
    """Main function with command line argument parsing."""
    parser = argparse.ArgumentParser(description='Test ERC Exploration functionality')
    parser.add_argument('--file', '-f', default='networks/RandomAlife/RN_Ns_40_Norg_20_id_396.txt',
                        help='Reaction network file to test (default: networks/testing/ERCs_test2.txt)')
    parser.add_argument('--multiple', '-m', action='store_true',
                        help='Test multiple networks')
    
    args = parser.parse_args()
    
    if args.multiple:
        return test_multiple_networks()
    else:
        tester = ERCExplorationTester(args.file)
        return tester.run_comprehensive_test()


if __name__ == "__main__":
    results = main()