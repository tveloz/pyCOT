#!/usr/bin/env python3
"""
Relevant Closures Computation - Efficient Algorithm for Computing Productively Novel Closed Sets

This module implements an efficient algorithm for computing all relevant (productively novel)
closed sets of a reaction network by working at the ERC level and exploiting hierarchy structure.

@author: Efficient COT Algorithm
"""

import os
import sys
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from bitarray import bitarray, frozenbitarray
from collections import defaultdict, deque
from typing import List, Dict, Set, Tuple, Optional, FrozenSet
import networkx as nx
from pathlib import Path
import argparse

# Import pyCOT modules
from pyCOT.io.functions import read_txt, from_string
from pyCOT.rn_rustworkx import ReactionNetwork
from pyCOT.ERC_Hierarchy import ERC, species_list_to_names, closure
from pyCOT.Dynamical_Hierarchy_Synergy_Algorithm import efficient_fundamental_synergies


class ERCClosure:
    """Represents a closed set at the ERC level with its minimal covers."""
    
    def __init__(self, erc_bitarray, signature):
        self.erc_bitarray = erc_bitarray  # Which ERCs are in this closure
        self.minimal_covers = []  # List of minimal bitarrays generating this closure
        self.signature = signature  # Hash for fast comparison
        self.cover_construction = {}  # Maps cover -> (n_synergies, n_req_additions)
        
    def add_minimal_cover(self, cover_bitarray, n_synergies, n_req_additions):
        """Add a minimal cover with its construction history."""
        frozen_cover = frozenbitarray(cover_bitarray)
        self.minimal_covers.append(cover_bitarray.copy())
        self.cover_construction[frozen_cover] = (n_synergies, n_req_additions)
        
    def get_size(self):
        """Get number of ERCs in closure."""
        return self.erc_bitarray.count()
        
    def get_num_minimal_covers(self):
        """Get number of minimal covers."""
        return len(self.minimal_covers)


class RelevantClosuresComputation:
    """
    Efficient algorithm for computing all productively novel closed sets.
    Works at ERC level using hierarchy structure and bitarray operations.
    """
    
    def __init__(self, reaction_network: ReactionNetwork, verbose: bool = True):
        self.rn = reaction_network
        self.verbose = verbose
        
        # Compute ERCs and hierarchy
        if verbose:
            print("Computing ERCs and hierarchy...")
        self.ercs = ERC.ERCs(reaction_network)
        self.n_ercs = len(self.ercs)
        self.hierarchy_graph = ERC.build_hierarchy_graph(self.ercs, reaction_network)
        
        # ERC lookup structures
        self.erc_by_label = {erc.label: erc for erc in self.ercs}
        self.erc_index = {erc.label: i for i, erc in enumerate(self.ercs)}
        self.index_to_erc = {i: erc for i, erc in enumerate(self.ercs)}
        
        # Compute fundamental synergies
        if verbose:
            print("Computing fundamental synergies...")
        self.fundamental_synergies = efficient_fundamental_synergies(
            reaction_network, self.ercs, self.hierarchy_graph, verbose=False
        )
        
        # Pre-compute hierarchy relationships as bitarrays
        self._precompute_hierarchy()
        self._precompute_synergies()
        
        # Results
        self.closure_registry = {}
        self.processed_covers = set()
        
        # Statistics
        self.stats = {
            'n_ercs': self.n_ercs,
            'n_fundamental_synergies': len(self.fundamental_synergies),
            'n_closures_found': 0,
            'n_minimal_covers': 0,
            'n_pure_synergy_covers': 0,
            'n_complementarity_covers': 0,
            'n_mixed_covers': 0,
            'computation_time': 0,
            'max_closure_size': 0,
            'max_covers_per_closure': 0
        }
        
    def _precompute_hierarchy(self):
        """Pre-compute hierarchy relationships as bitarrays."""
        if self.verbose:
            print("Pre-computing hierarchy relationships...")
            
        self.descendants = {}
        self.ancestors = {}
        
        for i, erc in enumerate(self.ercs):
            # Compute descendants (ERCs contained in this ERC)
            desc_bitarray = bitarray(self.n_ercs)
            desc_bitarray[i] = 1  # Include self
            
            desc_labels = nx.descendants(self.hierarchy_graph, erc.label)
            for desc_label in desc_labels:
                desc_bitarray[self.erc_index[desc_label]] = 1
                
            self.descendants[i] = desc_bitarray
            
            # Compute ancestors (ERCs containing this ERC)
            anc_bitarray = bitarray(self.n_ercs)
            anc_bitarray[i] = 1  # Include self
            
            anc_labels = nx.ancestors(self.hierarchy_graph, erc.label)
            for anc_label in anc_labels:
                anc_bitarray[self.erc_index[anc_label]] = 1
                
            self.ancestors[i] = anc_bitarray
            
    def _precompute_synergies(self):
        """Pre-compute synergy relationships at ERC level."""
        if self.verbose:
            print(f"Pre-computing synergy targets from {len(self.fundamental_synergies)} fundamental synergies...")
            
        self.synergy_targets = {}
        self.synergy_partners = defaultdict(set)  # Maps ERC index to partners
        
        for syn in self.fundamental_synergies:
            i1 = self.erc_index[syn['base1'].label]
            i2 = self.erc_index[syn['base2'].label]
            target_idx = self.erc_index[syn['target'].label]
            
            # Create bitarray for target and its descendants
            target_bitarray = self.descendants[target_idx].copy()
            
            # Store both directions
            self.synergy_targets[(i1, i2)] = target_bitarray
            self.synergy_targets[(i2, i1)] = target_bitarray
            
            # Track partners
            self.synergy_partners[i1].add(i2)
            self.synergy_partners[i2].add(i1)
            
    def erc_closure(self, seed_ercs: bitarray) -> bitarray:
        """Compute closure at ERC level using hierarchy and synergies."""
        closure = seed_ercs.copy()
        
        # Add all descendants
        for erc_id in [i for i in range(self.n_ercs) if seed_ercs[i]]:
            closure |= self.descendants[erc_id]
            
        # Iteratively add synergy products
        changed = True
        while changed:
            old_count = closure.count()
            
            # Check all pairs in current closure
            current_ercs = [i for i in range(self.n_ercs) if closure[i]]
            for e1 in current_ercs:
                for e2 in current_ercs:
                    if (e1, e2) in self.synergy_targets:
                        closure |= self.synergy_targets[(e1, e2)]
                        
            changed = (closure.count() > old_count)
            
        return closure
        
    def is_minimal_cover(self, new_cover: bitarray, existing_covers: List[bitarray]) -> bool:
        """Check if a cover is minimal compared to existing covers."""
        new_count = new_cover.count()
        
        for existing in existing_covers:
            existing_count = existing.count()
            
            # If existing is smaller or equal and covers all of new, new is not minimal
            if existing_count <= new_count:
                if (new_cover & existing) == existing:  # existing is subset of new
                    return False
                    
        return True
        
    def compute_erc_req(self, closure_bitarray: bitarray) -> Set[str]:
        """Compute required species at the species level for a closure."""
        # Get all species in the closure
        closure_species = set()
        for i in range(self.n_ercs):
            if closure_bitarray[i]:
                erc = self.index_to_erc[i]
                closure_species.update(species_list_to_names(erc.get_closure(self.rn)))
                
        # Get all reactions in the closure
        closure_reactions = []
        for i in range(self.n_ercs):
            if closure_bitarray[i]:
                erc = self.index_to_erc[i]
                closure_reactions.extend(erc.get_reacs(self.rn))
                
        # Compute consumed and produced
        consumed = set()
        produced = set()
        for reaction in closure_reactions:
            consumed.update(reaction.support_names())
            produced.update(reaction.products_names())
            
        # Required = consumed - produced
        return consumed - produced
        
    def find_synergetic_extensions(self, current_cover: bitarray, current_closure: bitarray) -> List[int]:
        """Find all ERCs that have synergy with current cover."""
        candidates = set()
        
        # For each ERC in cover, get its synergy partners
        for erc_id in [i for i in range(self.n_ercs) if current_cover[i]]:
            candidates.update(self.synergy_partners[erc_id])
            
        # Filter candidates
        valid_candidates = []
        for cand_id in candidates:
            # Skip if already in closure
            if current_closure[cand_id]:
                continue
                
            # Skip if ancestor of any ERC in cover
            is_ancestor = False
            for erc_id in [i for i in range(self.n_ercs) if current_cover[i]]:
                if self.ancestors[cand_id][erc_id]:
                    is_ancestor = True
                    break
                    
            if not is_ancestor:
                valid_candidates.append(cand_id)
                
        return valid_candidates
        
    def find_complementarity_candidates(self, current_closure: bitarray, current_req: Set[str]) -> List[int]:
        """Find ERCs that could reduce required species."""
        if not current_req:
            return []
            
        candidates = []
        
        for cand_id in range(self.n_ercs):
            # Skip if in closure
            if current_closure[cand_id]:
                continue
                
            # Get what this ERC produces
            erc = self.index_to_erc[cand_id]
            erc_products = set()
            for reaction in erc.get_reacs(self.rn):
                erc_products.update(reaction.products_names())
                
            # Check if it produces anything we need
            if erc_products & current_req:
                candidates.append(cand_id)
                
        return candidates
        
    def compute_relevant_closures(self):
        """Main algorithm to compute all productively novel closed sets."""
        start_time = time.time()
        
        if self.verbose:
            print("\n" + "="*70)
            print("COMPUTING RELEVANT CLOSURES")
            print("="*70)
            
        # Phase 0: Initialize with single ERCs
        queue = deque()
        for i in range(self.n_ercs):
            cover = bitarray(self.n_ercs)
            cover[i] = 1
            queue.append((cover, i, 'initial'))
            
        # Phase 1: Synergetic construction
        if self.verbose:
            print("\nPhase 1: Synergetic construction...")
            
        while queue:
            current_cover, last_added, add_type = queue.popleft()
            
            # Skip if already processed
            frozen_cover = frozenbitarray(current_cover)
            if frozen_cover in self.processed_covers:
                continue
            self.processed_covers.add(frozen_cover)
            
            # Compute closure
            current_closure = self.erc_closure(current_cover)
            signature = hash(frozenbitarray(current_closure))
            
            # Count synergies in construction path
            n_synergies = current_cover.count() - 1 if add_type == 'synergy' else 0
            
            # Update registry
            if signature in self.closure_registry:
                erc_closure_obj = self.closure_registry[signature]
                if self.is_minimal_cover(current_cover, erc_closure_obj.minimal_covers):
                    erc_closure_obj.add_minimal_cover(current_cover, n_synergies, 0)
                    self.stats['n_minimal_covers'] += 1
                    if n_synergies == current_cover.count() - 1:
                        self.stats['n_pure_synergy_covers'] += 1
            else:
                # New closure
                erc_closure_obj = ERCClosure(current_closure, signature)
                erc_closure_obj.add_minimal_cover(current_cover, n_synergies, 0)
                self.closure_registry[signature] = erc_closure_obj
                self.stats['n_closures_found'] += 1
                self.stats['n_minimal_covers'] += 1
                if n_synergies == current_cover.count() - 1:
                    self.stats['n_pure_synergy_covers'] += 1
                    
            # Find synergetic extensions
            extensions = self.find_synergetic_extensions(current_cover, current_closure)
            for ext_id in extensions:
                new_cover = current_cover.copy()
                new_cover[ext_id] = 1
                queue.append((new_cover, ext_id, 'synergy'))
                
        # Phase 2: Complementarity construction
        if self.verbose:
            print("\nPhase 2: Complementarity construction...")
            
        # Process each closure for complementarity extensions
        for erc_closure_obj in list(self.closure_registry.values()):
            current_req = self.compute_erc_req(erc_closure_obj.erc_bitarray)
            
            if not current_req:
                continue
                
            # Try extending each minimal cover
            for minimal_cover in erc_closure_obj.minimal_covers:
                candidates = self.find_complementarity_candidates(
                    erc_closure_obj.erc_bitarray, current_req
                )
                
                for cand_id in candidates:
                    new_cover = minimal_cover.copy()
                    new_cover[cand_id] = 1
                    
                    frozen_new = frozenbitarray(new_cover)
                    if frozen_new in self.processed_covers:
                        continue
                    self.processed_covers.add(frozen_new)
                    
                    # Check complementarity criterion
                    new_closure = self.erc_closure(new_cover)
                    new_req = self.compute_erc_req(new_closure)
                    
                    cand_single = bitarray(self.n_ercs)
                    cand_single[cand_id] = 1
                    cand_req = self.compute_erc_req(cand_single)
                    
                    if len(new_req) < len(current_req) + len(cand_req):
                        # Valid complementarity
                        signature = hash(frozenbitarray(new_closure))
                        
                        frozen_minimal = frozenbitarray(minimal_cover)
                        n_syn = erc_closure_obj.cover_construction[frozen_minimal][0]
                        n_req = erc_closure_obj.cover_construction[frozen_minimal][1] + 1
                        
                        if signature in self.closure_registry:
                            closure_obj = self.closure_registry[signature]
                            if self.is_minimal_cover(new_cover, closure_obj.minimal_covers):
                                closure_obj.add_minimal_cover(new_cover, n_syn, n_req)
                                self.stats['n_minimal_covers'] += 1
                                self.stats['n_complementarity_covers'] += 1
                        else:
                            # New closure from complementarity
                            closure_obj = ERCClosure(new_closure, signature)
                            closure_obj.add_minimal_cover(new_cover, n_syn, n_req)
                            self.closure_registry[signature] = closure_obj
                            self.stats['n_closures_found'] += 1
                            self.stats['n_minimal_covers'] += 1
                            self.stats['n_complementarity_covers'] += 1
                            
        # Update final statistics
        self.stats['computation_time'] = time.time() - start_time
        self.stats['max_closure_size'] = max(
            (c.get_size() for c in self.closure_registry.values()), default=0
        )
        self.stats['max_covers_per_closure'] = max(
            (c.get_num_minimal_covers() for c in self.closure_registry.values()), default=0
        )
        
        # Count mixed covers
        for closure_obj in self.closure_registry.values():
            for cover, (n_syn, n_req) in closure_obj.cover_construction.items():
                if n_syn > 0 and n_req > 0:
                    self.stats['n_mixed_covers'] += 1
                    
        if self.verbose:
            print(f"\nComputation completed in {self.stats['computation_time']:.2f} seconds")
            print(f"Found {self.stats['n_closures_found']} relevant closures")
            print(f"Total minimal covers: {self.stats['n_minimal_covers']}")
            
    def get_statistics_summary(self) -> pd.DataFrame:
        """Get summary statistics as DataFrame."""
        # Calculate theoretical maximum
        theoretical_max = 2**self.n_ercs - 1  # All non-empty subsets
        
        stats_data = {
            'Metric': [
                'Number of ERCs',
                'Number of Fundamental Synergies',
                'Relevant Closures Found',
                'Theoretical Maximum (2^n - 1)',
                'Reduction Factor',
                'Total Minimal Covers',
                'Pure Synergy Covers',
                'Complementarity Covers',
                'Mixed Covers',
                'Max Closure Size (ERCs)',
                'Max Covers per Closure',
                'Computation Time (seconds)'
            ],
            'Value': [
                self.stats['n_ercs'],
                self.stats['n_fundamental_synergies'],
                self.stats['n_closures_found'],
                theoretical_max,
                f"{theoretical_max / max(self.stats['n_closures_found'], 1):.1f}x",
                self.stats['n_minimal_covers'],
                self.stats['n_pure_synergy_covers'],
                self.stats['n_complementarity_covers'],
                self.stats['n_mixed_covers'],
                self.stats['max_closure_size'],
                self.stats['max_covers_per_closure'],
                f"{self.stats['computation_time']:.3f}"
            ]
        }
        
        return pd.DataFrame(stats_data)
        
    def get_closure_details(self) -> pd.DataFrame:
        """Get detailed information about each closure."""
        closure_data = []
        
        for closure_obj in self.closure_registry.values():
            for cover, (n_syn, n_req) in closure_obj.cover_construction.items():
                # Get ERC labels in cover
                cover_ercs = [self.index_to_erc[i].label 
                             for i in range(self.n_ercs) if cover[i]]
                
                # Get all ERCs in closure
                closure_ercs = [self.index_to_erc[i].label 
                               for i in range(self.n_ercs) if closure_obj.erc_bitarray[i]]
                
                closure_data.append({
                    'Closure_Size': closure_obj.get_size(),
                    'Cover_Size': cover.count(),
                    'Cover_ERCs': ', '.join(cover_ercs),
                    'N_Synergies': n_syn,
                    'N_Req_Additions': n_req,
                    'Type': 'Pure Synergy' if n_req == 0 and n_syn > 0 
                           else 'Pure Complementarity' if n_syn == 0 
                           else 'Mixed' if n_syn > 0 and n_req > 0
                           else 'Single ERC'
                })
                
        return pd.DataFrame(closure_data).sort_values(['Closure_Size', 'Cover_Size'])
        
    def plot_results(self, save_path: Optional[str] = None, show_plots: bool = True):
        """Generate visualization plots."""
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        fig.suptitle('Relevant Closures Analysis', fontsize=16)
        
        # 1. Closure size distribution
        ax = axes[0, 0]
        closure_sizes = [c.get_size() for c in self.closure_registry.values()]
        ax.hist(closure_sizes, bins=min(20, max(closure_sizes)), 
                edgecolor='black', alpha=0.7)
        ax.set_xlabel('Closure Size (# ERCs)')
        ax.set_ylabel('Count')
        ax.set_title('Distribution of Closure Sizes')
        ax.grid(True, alpha=0.3)
        
        # 2. Cover types pie chart
        ax = axes[0, 1]
        cover_types = {
            'Pure Synergy': self.stats['n_pure_synergy_covers'],
            'Complementarity': self.stats['n_complementarity_covers'] - self.stats['n_mixed_covers'],
            'Mixed': self.stats['n_mixed_covers']
        }
        colors = ['#2ecc71', '#3498db', '#e74c3c']
        wedges, texts, autotexts = ax.pie(
            cover_types.values(), labels=cover_types.keys(), colors=colors,
            autopct='%1.1f%%', startangle=90
        )
        ax.set_title('Distribution of Cover Types')
        
        # 3. Covers per closure
        ax = axes[1, 0]
        covers_per_closure = [c.get_num_minimal_covers() 
                             for c in self.closure_registry.values()]
        unique_counts = sorted(set(covers_per_closure))
        count_freq = [covers_per_closure.count(c) for c in unique_counts]
        
        ax.bar(unique_counts, count_freq, edgecolor='black', alpha=0.7)
        ax.set_xlabel('Number of Minimal Covers')
        ax.set_ylabel('Number of Closures')
        ax.set_title('Distribution of Minimal Covers per Closure')
        ax.set_xticks(unique_counts)
        ax.grid(True, alpha=0.3, axis='y')
        
        # 4. Efficiency comparison
        ax = axes[1, 1]
        categories = ['Theoretical\nMaximum', 'Relevant\nClosures', 'Reduction\nFactor']
        theoretical = 2**self.stats['n_ercs'] - 1
        values = [theoretical, self.stats['n_closures_found'], 
                 theoretical / max(self.stats['n_closures_found'], 1)]
        
        bars = ax.bar(categories[:2], values[:2], color=['#e74c3c', '#2ecc71'], 
                      edgecolor='black', alpha=0.7)
        ax.set_ylabel('Count (log scale)')
        ax.set_yscale('log')
        ax.set_title('Efficiency: Relevant vs All Possible Closures')
        
        # Add reduction factor as text
        ax.text(0.5, theoretical/2, f"{values[2]:.1f}x\nreduction", 
               ha='center', va='center', fontsize=14, fontweight='bold',
               bbox=dict(boxstyle="round,pad=0.3", facecolor="yellow", alpha=0.7))
        
        ax.grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            if self.verbose:
                print(f"\nPlots saved to: {save_path}")
                
        if show_plots:
            plt.show()
        else:
            plt.close()
            
        return fig


def process_network_file(file_path: str, verbose: bool = True, 
                        save_plots: bool = True, show_plots: bool = True) -> Dict:
    """Process a single network file and return results."""
    if verbose:
        print(f"\nProcessing: {file_path}")
        print("-" * 70)
        
    # Load network
    try:
        rn = read_txt(file_path)
    except Exception as e:
        print(f"Error loading network: {e}")
        return None
        
    # Compute relevant closures
    rcc = RelevantClosuresComputation(rn, verbose=verbose)
    rcc.compute_relevant_closures()
    
    # Get results
    results = {
        'file_path': file_path,
        'stats': rcc.stats,
        'summary_df': rcc.get_statistics_summary(),
        'details_df': rcc.get_closure_details(),
        'computation': rcc
    }
    
    # Generate plots if requested
    if save_plots or show_plots:
        plot_path = None
        if save_plots:
            base_name = Path(file_path).stem
            plot_path = f"{base_name}_relevant_closures_analysis.png"
            
        rcc.plot_results(save_path=plot_path, show_plots=show_plots)
        
    return results


def process_folder(folder_path: str, verbose: bool = True,
                  save_plots: bool = True, show_plots: bool = False) -> List[Dict]:
    """Process all .txt files in a folder."""
    results = []
    txt_files = list(Path(folder_path).glob("*.txt"))
    
    if verbose:
        print(f"\nFound {len(txt_files)} network files in {folder_path}")
        
    for file_path in txt_files:
        result = process_network_file(
            str(file_path), verbose=verbose, 
            save_plots=save_plots, show_plots=show_plots
        )
        if result:
            results.append(result)
            
    return results


def compare_results(results: List[Dict], save_path: Optional[str] = None,
                   show_plot: bool = True):
    """Compare results across multiple networks."""
    if not results:
        print("No results to compare")
        return
        
    # Create comparison DataFrame
    comparison_data = []
    for result in results:
        stats = result['stats']
        comparison_data.append({
            'Network': Path(result['file_path']).stem,
            'ERCs': stats['n_ercs'],
            'Synergies': stats['n_fundamental_synergies'],
            'Closures': stats['n_closures_found'],
            'Reduction': 2**stats['n_ercs'] - 1 - stats['n_closures_found'],
            'Time (s)': stats['computation_time']
        })
        
    df = pd.DataFrame(comparison_data)
    
    # Create comparison plot
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    fig.suptitle('Network Comparison', fontsize=16)
    
    # 1. Network sizes
    ax = axes[0, 0]
    x = np.arange(len(df))
    width = 0.35
    ax.bar(x - width/2, df['ERCs'], width, label='ERCs', alpha=0.7)
    ax.bar(x + width/2, df['Synergies'], width, label='Synergies', alpha=0.7)
    ax.set_xlabel('Network')
    ax.set_ylabel('Count')
    ax.set_title('Network Characteristics')
    ax.set_xticks(x)
    ax.set_xticklabels(df['Network'], rotation=45, ha='right')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 2. Closures found
    ax = axes[0, 1]
    ax.bar(x, df['Closures'], alpha=0.7, color='green', edgecolor='black')
    ax.set_xlabel('Network')
    ax.set_ylabel('Relevant Closures')
    ax.set_title('Number of Relevant Closures')
    ax.set_xticks(x)
    ax.set_xticklabels(df['Network'], rotation=45, ha='right')
    ax.grid(True, alpha=0.3)
    
    # 3. Reduction efficiency
    ax = axes[1, 0]
    theoretical = [2**n - 1 for n in df['ERCs']]
    reduction_pct = [(t - c) / t * 100 for t, c in zip(theoretical, df['Closures'])]
    ax.bar(x, reduction_pct, alpha=0.7, color='blue', edgecolor='black')
    ax.set_xlabel('Network')
    ax.set_ylabel('Reduction %')
    ax.set_title('Computational Reduction')
    ax.set_xticks(x)
    ax.set_xticklabels(df['Network'], rotation=45, ha='right')
    ax.grid(True, alpha=0.3)
    
    # 4. Computation time
    ax = axes[1, 1]
    ax.bar(x, df['Time (s)'], alpha=0.7, color='orange', edgecolor='black')
    ax.set_xlabel('Network')
    ax.set_ylabel('Time (seconds)')
    ax.set_title('Computation Time')
    ax.set_xticks(x)
    ax.set_xticklabels(df['Network'], rotation=45, ha='right')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"\nComparison plot saved to: {save_path}")
        
    if show_plot:
        plt.show()
    else:
        plt.close()
        
    return df


def main():
    """Main function for command-line usage."""
    parser = argparse.ArgumentParser(
        description='Compute relevant closures (productively novel closed sets) for reaction networks'
    )
    parser.add_argument('input', help='Input file or folder path')
    parser.add_argument('--quiet', '-q', action='store_true', 
                       help='Suppress verbose output')
    parser.add_argument('--no-plots', action='store_true',
                       help='Skip generating plots')
    parser.add_argument('--save-plots', action='store_true',
                       help='Save plots to files')
    parser.add_argument('--export-csv', action='store_true',
                       help='Export results to CSV files')
    
    args = parser.parse_args()
    
    verbose = not args.quiet
    show_plots = not args.no_plots
    save_plots = args.save_plots or args.no_plots
    
    input_path = Path(args.input)
    
    if input_path.is_file():
        # Process single file
        result = process_network_file(
            str(input_path), verbose=verbose,
            save_plots=save_plots, show_plots=show_plots
        )
        
        if result and args.export_csv:
            base_name = input_path.stem
            result['summary_df'].to_csv(f"{base_name}_summary.csv", index=False)
            result['details_df'].to_csv(f"{base_name}_details.csv", index=False)
            if verbose:
                print(f"\nResults exported to CSV files")
                
    elif input_path.is_dir():
        # Process folder
        results = process_folder(
            str(input_path), verbose=verbose,
            save_plots=save_plots, show_plots=show_plots
        )
        
        if results:
            # Generate comparison
            comparison_df = compare_results(
                results, 
                save_path="network_comparison.png" if save_plots else None,
                show_plot=show_plots
            )
            
            if args.export_csv and comparison_df is not None:
                comparison_df.to_csv("network_comparison.csv", index=False)
                if verbose:
                    print(f"\nComparison results exported to network_comparison.csv")
    else:
        print(f"Error: {input_path} is neither a file nor a directory")
        return 1
        
    return 0


if __name__ == "__main__":
  
    
    # Demo with a small test network
    test_network = 'networks/Navarino/RN_IN_05.txt'
    test_network="networks/RandomAlife/RN_Ns_40_Norg_20_id_396.txt"
    test_network="networks/biomodels_interesting/central_ecoli.txt"

    rn = read_txt(test_network)
    rcc = RelevantClosuresComputation(rn, verbose=True)
    rcc.compute_relevant_closures()
    
    print("\n" + "="*70)
    print("RESULTS SUMMARY")
    print("="*70)
    print(rcc.get_statistics_summary().to_string(index=False))
    
    print("\n" + "="*70)
    print("CLOSURE DETAILS")
    print("="*70)
    print(rcc.get_closure_details().to_string(index=False))
    
    rcc.plot_results(show_plots=True)
