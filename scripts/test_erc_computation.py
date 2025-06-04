#!/usr/bin/env python3
"""
Simplified ERC Analysis Script
Focuses on basic statistics and generator analysis only.
"""

import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from typing import List, Dict, Tuple, Set
import warnings
warnings.filterwarnings('ignore')

# Import the necessary classes
from pyCOT.ERC_Hierarchy import ERC, generators, closure, species_list_to_names
from pyCOT.io.functions import read_txt
from pyCOT.rn_rustworkx import ReactionNetwork

def levenshtein_distance(s1: List[str], s2: List[str]) -> int:
    """Calculate Levenshtein distance between two lists of strings."""
    if len(s1) < len(s2):
        return levenshtein_distance(s2, s1)
    
    if len(s2) == 0:
        return len(s1)
    
    previous_row = list(range(len(s2) + 1))
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row
    
    return previous_row[-1]

def jaccard_similarity(set1: Set[str], set2: Set[str]) -> float:
    """Calculate Jaccard similarity between two sets."""
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    return intersection / union if union > 0 else 0

class ERCAnalyzer:
    def __init__(self, reaction_network_file: str):
        """Initialize the ERC analyzer with a reaction network file."""
        self.rn = read_txt(reaction_network_file)
        self.ercs = None
        
    def compute_ercs(self):
        """Compute ERCs from the reaction network."""
        print("Computing ERCs...")
        self.ercs = ERC.ERCs(self.rn)
        print(f"Found {len(self.ercs)} ERCs")
        
    def get_basic_statistics(self) -> Dict:
        """Get basic statistics about the reaction network and ERCs."""
        if self.ercs is None:
            self.compute_ercs()
            
        stats = {}
        
        # Network-level statistics
        stats['total_species'] = len(self.rn.species())
        stats['total_reactions'] = len(self.rn.reactions())
        stats['total_ercs'] = len(self.ercs)
        
        # Minimal generators statistics
        min_gen_counts = [len(erc.min_generators) for erc in self.ercs]
        stats['avg_minimal_generators'] = np.mean(min_gen_counts) if min_gen_counts else 0
        stats['std_minimal_generators'] = np.std(min_gen_counts) if min_gen_counts else 0
        stats['min_minimal_generators'] = np.min(min_gen_counts) if min_gen_counts else 0
        stats['max_minimal_generators'] = np.max(min_gen_counts) if min_gen_counts else 0
        stats['total_minimal_generators'] = sum(min_gen_counts)
        
        # ERC closure sizes
        closure_sizes = [len(erc.get_closure(self.rn)) for erc in self.ercs]
        stats['avg_species_per_erc'] = np.mean(closure_sizes) if closure_sizes else 0
        stats['std_species_per_erc'] = np.std(closure_sizes) if closure_sizes else 0
        stats['min_species_per_erc'] = np.min(closure_sizes) if closure_sizes else 0
        stats['max_species_per_erc'] = np.max(closure_sizes) if closure_sizes else 0
        
        # Generator size statistics
        all_gen_sizes = [len(gen) for erc in self.ercs for gen in erc.min_generators]
        stats['avg_generator_size'] = np.mean(all_gen_sizes) if all_gen_sizes else 0
        stats['min_generator_size'] = np.min(all_gen_sizes) if all_gen_sizes else 0
        stats['max_generator_size'] = np.max(all_gen_sizes) if all_gen_sizes else 0
        
        return stats
    
    def calculate_generator_similarities(self) -> Dict:
        """Calculate similarities between minimal generators across all ERCs."""
        if self.ercs is None:
            self.compute_ercs()
            
        similarities = {
            'levenshtein_distances': [],
            'jaccard_similarities': [],
            'generator_overlaps': [],
            'pairwise_details': []
        }
        
        # Get all generators flattened
        all_generators = []
        for erc in self.ercs:
            for gen in erc.min_generators:
                gen_species = species_list_to_names(gen)
                all_generators.append({
                    'erc_label': erc.label,
                    'species': gen_species,
                    'species_set': set(gen_species)
                })
        
        # Compare all pairs of generators
        for i, gen1 in enumerate(all_generators):
            for j, gen2 in enumerate(all_generators[i+1:], i+1):
                
                # Calculate Levenshtein distance
                lev_dist = levenshtein_distance(gen1['species'], gen2['species'])
                similarities['levenshtein_distances'].append(lev_dist)
                
                # Calculate Jaccard similarity
                jaccard_sim = jaccard_similarity(gen1['species_set'], gen2['species_set'])
                similarities['jaccard_similarities'].append(jaccard_sim)
                
                # Calculate overlap count
                overlap = len(gen1['species_set'].intersection(gen2['species_set']))
                similarities['generator_overlaps'].append(overlap)
                
                # Store pairwise details
                similarities['pairwise_details'].append({
                    'erc1': gen1['erc_label'],
                    'erc2': gen2['erc_label'],
                    'gen1': gen1['species'],
                    'gen2': gen2['species'],
                    'levenshtein': lev_dist,
                    'jaccard': jaccard_sim,
                    'overlap': overlap
                })
        
        # Calculate summary statistics
        similarities['avg_levenshtein_distance'] = np.mean(similarities['levenshtein_distances']) if similarities['levenshtein_distances'] else 0
        similarities['avg_jaccard_similarity'] = np.mean(similarities['jaccard_similarities']) if similarities['jaccard_similarities'] else 0
        similarities['avg_generator_overlap'] = np.mean(similarities['generator_overlaps']) if similarities['generator_overlaps'] else 0
        similarities['max_jaccard_similarity'] = np.max(similarities['jaccard_similarities']) if similarities['jaccard_similarities'] else 0
        similarities['min_jaccard_similarity'] = np.min(similarities['jaccard_similarities']) if similarities['jaccard_similarities'] else 0
        
        return similarities
    
    def analyze_generator_patterns(self) -> Dict:
        """Analyze patterns in generator composition."""
        if self.ercs is None:
            self.compute_ercs()
            
        patterns = {
            'species_frequency': Counter(),
            'generator_size_distribution': Counter(),
            'unique_generators': set(),
            'repeated_generators': [],
            'erc_generator_details': []
        }
        
        # Analyze each ERC's generators
        for erc in self.ercs:
            erc_details = {
                'label': erc.label,
                'generators': [],
                'total_generators': len(erc.min_generators),
                'closure_size': len(erc.get_closure(self.rn))
            }
            
            for gen in erc.min_generators:
                gen_species = species_list_to_names(gen)
                gen_tuple = tuple(sorted(gen_species))
                
                # Track generator details
                erc_details['generators'].append(gen_species)
                
                # Count species frequency across all generators
                for species in gen_species:
                    patterns['species_frequency'][species] += 1
                
                # Track generator size distribution
                patterns['generator_size_distribution'][len(gen_species)] += 1
                
                # Track unique vs repeated generators
                if gen_tuple in patterns['unique_generators']:
                    patterns['repeated_generators'].append(gen_tuple)
                else:
                    patterns['unique_generators'].add(gen_tuple)
            
            patterns['erc_generator_details'].append(erc_details)
        
        return patterns
    
    def print_ascii_chart(self, data: List[Tuple[str, int]], title: str, max_width: int = 20):
        """Print an ASCII bar chart."""
        print(f"\n📈 {title.upper()}:")
        print("-" * 40)
        if not data:
            print("   No data available")
            return
            
        max_val = max(val for _, val in data)
        for label, value in data:
            bar_length = int((value / max_val) * max_width) if max_val > 0 else 0
            bar = "█" * bar_length + "░" * (max_width - bar_length)
            print(f"{label}: {bar} {value}")
    
    def analyze_reaction_coverage(self) -> Dict:
        """Analyze which reactions are covered by each ERC."""
        if self.ercs is None:
            self.compute_ercs()
            
        coverage = {}
        all_reactions = set(r.name() for r in self.rn.reactions())
        covered_reactions = set()
        
        for erc in self.ercs:
            reactions = erc.get_reacs(self.rn)
            erc_reactions = [r.name() for r in reactions]
            coverage[erc.label] = erc_reactions
            covered_reactions.update(erc_reactions)
        
        coverage['total_reactions'] = len(all_reactions)
        coverage['covered_reactions'] = len(covered_reactions)
        coverage['coverage_ratio'] = len(covered_reactions) / len(all_reactions) if all_reactions else 0
        coverage['uncovered_reactions'] = list(all_reactions - covered_reactions)
        
        return coverage
    
    def generate_comprehensive_report(self):
        """Generate a comprehensive analysis report focusing on basic stats and generators."""
        if self.ercs is None:
            self.compute_ercs()
            
        print("\n" + "="*80)
        print("🎯 ERC BASIC ANALYSIS REPORT")
        print("="*80)
        
        # === BASIC STATISTICS ===
        basic_stats = self.get_basic_statistics()
        print(f"\n📊 BASIC STATISTICS:")
        print("-" * 30)
        
        # Network overview
        print(f"Network Overview:")
        print(f"  Total Species:        {basic_stats['total_species']}")
        print(f"  Total Reactions:      {basic_stats['total_reactions']}")
        print(f"  Total ERCs:           {basic_stats['total_ercs']}")
        
        # Generator statistics
        print(f"\nMinimal Generators:")
        print(f"  Total Generators:     {basic_stats['total_minimal_generators']}")
        print(f"  Avg per ERC:          {basic_stats['avg_minimal_generators']:.2f} ± {basic_stats['std_minimal_generators']:.2f}")
        print(f"  Range per ERC:        {basic_stats['min_minimal_generators']} - {basic_stats['max_minimal_generators']}")
        print(f"  Avg Generator Size:   {basic_stats['avg_generator_size']:.2f}")
        print(f"  Generator Size Range: {basic_stats['min_generator_size']} - {basic_stats['max_generator_size']}")
        
        # Closure statistics
        print(f"\nERC Closures:")
        print(f"  Avg Species per ERC:  {basic_stats['avg_species_per_erc']:.2f} ± {basic_stats['std_species_per_erc']:.2f}")
        print(f"  Range:                {basic_stats['min_species_per_erc']} - {basic_stats['max_species_per_erc']}")
        
        # Reaction coverage
        coverage = self.analyze_reaction_coverage()
        print(f"\nReaction Coverage:")
        print(f"  Covered Reactions:    {coverage['covered_reactions']}/{coverage['total_reactions']} ({coverage['coverage_ratio']:.1%})")
        if coverage['uncovered_reactions']:
            print(f"  Uncovered:            {', '.join(coverage['uncovered_reactions'])}")
        
        # === ASCII VISUALIZATIONS ===
        
        # Species distribution per ERC
        species_data = [(erc.label, len(erc.get_closure(self.rn))) for erc in self.ercs]
        self.print_ascii_chart(species_data, "Species per ERC")
        
        # Generators distribution per ERC
        generator_data = [(erc.label, len(erc.min_generators)) for erc in self.ercs]
        self.print_ascii_chart(generator_data, "Generators per ERC")
        
        # Generator size distribution
        all_gen_sizes = [len(gen) for erc in self.ercs for gen in erc.min_generators]
        size_counter = Counter(all_gen_sizes)
        size_data = [(f"Size {size}", count) for size, count in sorted(size_counter.items())]
        self.print_ascii_chart(size_data, "Generator Size Distribution")
        
        # === GENERATOR ANALYSIS ===
        similarities = self.calculate_generator_similarities()
        patterns = self.analyze_generator_patterns()
        
        print(f"\n🔍 GENERATOR SIMILARITY ANALYSIS:")
        print("-" * 35)
        print(f"Pairwise Comparisons:     {len(similarities['pairwise_details'])}")
        print(f"Avg Levenshtein Distance: {similarities['avg_levenshtein_distance']:.2f}")
        print(f"Avg Jaccard Similarity:   {similarities['avg_jaccard_similarity']:.3f}")
        print(f"Jaccard Range:            {similarities['min_jaccard_similarity']:.3f} - {similarities['max_jaccard_similarity']:.3f}")
        print(f"Avg Generator Overlap:    {similarities['avg_generator_overlap']:.2f} species")
        
        # Most frequent species in generators
        print(f"\n🧬 GENERATOR COMPOSITION PATTERNS:")
        print("-" * 35)
        print(f"Unique Generator Types:   {len(patterns['unique_generators'])}")
        print(f"Repeated Generators:      {len(patterns['repeated_generators'])}")
        
        # Show most frequent species in generators
        most_frequent = patterns['species_frequency'].most_common(5)
        print(f"\nMost Frequent Generator Species:")
        for species, count in most_frequent:
            percentage = (count / basic_stats['total_minimal_generators']) * 100
            print(f"  {species}: {count} times ({percentage:.1f}%)")
        
        # Species frequency chart
        if most_frequent:
            freq_data = [(species, count) for species, count in most_frequent]
            self.print_ascii_chart(freq_data, "Species Frequency in Generators")
        
        # === DETAILED ERC BREAKDOWN ===
        print(f"\n📋 DETAILED ERC BREAKDOWN:")
        print("-" * 35)
        print("ERC | Generators | Species | Reactions | Generator Details")
        print("-" * 70)
        
        for erc in self.ercs:
            closure = erc.get_closure(self.rn)
            reactions = erc.get_reacs(self.rn)
            gen_details = [species_list_to_names(gen) for gen in erc.min_generators]
            gen_details_str = '; '.join([f"[{','.join(gen)}]" for gen in gen_details])
            
            print(f"{erc.label:3} | {len(erc.min_generators):10} | {len(closure):7} | {len(reactions):9} | {gen_details_str}")
        
        # Show detailed generator similarities
        if similarities['pairwise_details']:
            print(f"\n🔍 PAIRWISE GENERATOR SIMILARITIES:")
            print("-" * 45)
            print("ERC1-ERC2 | Jaccard | Overlap | Levenshtein | Details")
            print("-" * 65)
            
            # Sort by Jaccard similarity (descending)
            sorted_pairs = sorted(similarities['pairwise_details'], 
                                key=lambda x: x['jaccard'], reverse=True)
            
            # Show top 10 most similar pairs
            for pair in sorted_pairs[:10]:
                erc_pair = f"{pair['erc1']}-{pair['erc2']}"
                details = f"{pair['gen1']} vs {pair['gen2']}"
                print(f"{erc_pair:9} | {pair['jaccard']:7.3f} | {pair['overlap']:7} | {pair['levenshtein']:11} | {details}")
        
        # === DETAILED CLOSURES ===
        print(f"\n📋 DETAILED ERC CLOSURES:")
        print("-" * 35)
        for erc in self.ercs:
            closure = erc.get_closure(self.rn)
            closure_names = species_list_to_names(closure)
            min_gen_names = [species_list_to_names(gen) for gen in erc.min_generators]
            reactions = erc.get_reacs(self.rn)
            reaction_names = [r.name() for r in reactions]
            
            print(f"\n{erc.label}: {len(closure)} species, {len(erc.min_generators)} generators, {len(reactions)} reactions")
            print(f"   Generators: {min_gen_names}")
            print(f"   Closure:    {closure_names}")
            print(f"   Reactions:  {reaction_names}")
        
        print("\n" + "="*80)
        print("✅ ERC BASIC ANALYSIS COMPLETE")
        print("="*80)
        
        # Try to create matplotlib visualizations if available
        try:
            self.create_visualizations(basic_stats, similarities, patterns)
        except ImportError:
            print("\n📊 Matplotlib not available - ASCII charts shown above")
        except Exception as e:
            print(f"\n⚠️  Visualization error: {e}")
    
    def create_visualizations(self, basic_stats, similarities, patterns):
        """Create matplotlib visualizations."""
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        fig.suptitle('ERC Analysis Visualizations', fontsize=16, fontweight='bold')
        
        # 1. Species per ERC
        ax = axes[0, 0]
        erc_labels = [erc.label for erc in self.ercs]
        closure_sizes = [len(erc.get_closure(self.rn)) for erc in self.ercs]
        ax.bar(erc_labels, closure_sizes, color='skyblue', alpha=0.7, edgecolor='black')
        ax.set_title('Species per ERC')
        ax.set_xlabel('ERC')
        ax.set_ylabel('Number of Species')
        ax.grid(True, alpha=0.3)
        
        # 2. Generators per ERC
        ax = axes[0, 1]
        gen_counts = [len(erc.min_generators) for erc in self.ercs]
        ax.bar(erc_labels, gen_counts, color='lightcoral', alpha=0.7, edgecolor='black')
        ax.set_title('Generators per ERC')
        ax.set_xlabel('ERC')
        ax.set_ylabel('Number of Generators')
        ax.grid(True, alpha=0.3)
        
        # 3. Generator size distribution
        ax = axes[0, 2]
        all_gen_sizes = [len(gen) for erc in self.ercs for gen in erc.min_generators]
        if all_gen_sizes:
            ax.hist(all_gen_sizes, bins=range(min(all_gen_sizes), max(all_gen_sizes)+2), 
                   alpha=0.7, color='lightgreen', edgecolor='black')
        ax.set_title('Generator Size Distribution')
        ax.set_xlabel('Generator Size')
        ax.set_ylabel('Frequency')
        ax.grid(True, alpha=0.3)
        
        # 4. Jaccard similarities distribution
        ax = axes[1, 0]
        if similarities['jaccard_similarities']:
            ax.hist(similarities['jaccard_similarities'], bins=20, alpha=0.7, color='gold', edgecolor='black')
        ax.set_title('Jaccard Similarities Distribution')
        ax.set_xlabel('Jaccard Similarity')
        ax.set_ylabel('Frequency')
        ax.grid(True, alpha=0.3)
        
        # 5. Species frequency in generators
        ax = axes[1, 1]
        most_frequent = patterns['species_frequency'].most_common(10)
        if most_frequent:
            species, counts = zip(*most_frequent)
            ax.barh(range(len(species)), counts, color='mediumpurple', alpha=0.7, edgecolor='black')
            ax.set_yticks(range(len(species)))
            ax.set_yticklabels(species)
            ax.set_title('Most Frequent Species in Generators')
            ax.set_xlabel('Frequency')
        ax.grid(True, alpha=0.3)
        
        # 6. ERC size vs generators scatter
        ax = axes[1, 2]
        ax.scatter(closure_sizes, gen_counts, alpha=0.7, color='orange', s=100, edgecolor='black')
        for i, erc in enumerate(self.ercs):
            ax.annotate(erc.label, (closure_sizes[i], gen_counts[i]), 
                       xytext=(5, 5), textcoords='offset points', fontsize=10)
        ax.set_title('ERC Size vs Number of Generators')
        ax.set_xlabel('Number of Species in Closure')
        ax.set_ylabel('Number of Generators')
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.show()

def main():
    """Main function to run the simplified ERC analysis."""
    # Analyze the provided test file
    file_path = 'networks/testing/Farm.txt'
    file_path = 'networks/RandomAlife/RN_Ns_40_Norg_20_id_396.txt'
    file_path = 'networks/Navarino/RN_IN_05.txt'
    file_path = 'networks/biomodels_interesting/BIOMD0000000237_manyOrgs.txt'    
      # Adjust path as needed
    analyzer = ERCAnalyzer(file_path)
    analyzer.generate_comprehensive_report()

if __name__ == "__main__":
    main()