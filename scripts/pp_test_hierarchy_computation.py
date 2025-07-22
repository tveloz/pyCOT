#!/usr/bin/env python3
"""
Modular ERC Hierarchy Analysis Script

This script provides focused, separate analyses of ERC hierarchies with 
size filtering capabilities and performance-based clustering.

Usage:
    python erc_hierarchy_analysis.py <erc_hierarchy_benchmark.csv>
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import RobustScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.manifold import TSNE
from sklearn.metrics import silhouette_score
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Try to import UMAP, fall back to t-SNE if not available
try:
    from umap import UMAP
    HAS_UMAP = True
except ImportError:
    HAS_UMAP = False

class ModularERCAnalyzer:
    def __init__(self, data_file):
        """Initialize analyzer with ERC hierarchy benchmark results"""
        self.df = pd.read_csv(data_file)
        self.original_df = self.df.copy()  # Keep original for reference
        print(f"Loaded {len(self.df)} networks with {len(self.df.columns)} columns")
        self.inspect_data()
        self.prepare_features()
    
    def inspect_data(self):
        """Inspect data and show filtering options"""
        print("\n" + "="*60)
        print("DATASET OVERVIEW")
        print("="*60)
        
        # Show size ranges for filtering
        size_cols = ['n_species', 'n_reactions', 'n_ercs']
        for col in size_cols:
            if col in self.df.columns:
                print(f"{col}: {self.df[col].min():.0f} - {self.df[col].max():.0f}")
        
        # Show available algorithms
        perf_cols = ['original_time', 'optimized_time', 'ultra_optimized_time']
        available_algs = [col for col in perf_cols if col in self.df.columns]
        print(f"\nAvailable algorithms: {[col.replace('_time', '') for col in available_algs]}")
        
        # Success rates
        success_cols = ['success1', 'success2', 'success3'] 
        for i, col in enumerate(success_cols):
            if col in self.df.columns:
                alg_name = ['Original', 'Optimized', 'Ultra-optimized'][i]
                success_rate = self.df[col].mean() * 100
                print(f"{alg_name} success rate: {success_rate:.1f}%")
    
    def apply_size_filters(self, min_species=None, max_species=None,
                          min_reactions=None, max_reactions=None,
                          min_ercs=2, max_ercs=None):
        """Apply size filters to the dataset"""
        print("\n" + "="*60)
        print("APPLYING SIZE FILTERS")
        print("="*60)
        
        initial_size = len(self.df)
        
        # Apply filters
        if min_species is not None and 'n_species' in self.df.columns:
            self.df = self.df[self.df['n_species'] >= min_species]
            print(f"Min species >= {min_species}: {len(self.df)} networks remaining")
        
        if max_species is not None and 'n_species' in self.df.columns:
            self.df = self.df[self.df['n_species'] <= max_species]
            print(f"Max species <= {max_species}: {len(self.df)} networks remaining")
        
        if min_reactions is not None and 'n_reactions' in self.df.columns:
            self.df = self.df[self.df['n_reactions'] >= min_reactions]
            print(f"Min reactions >= {min_reactions}: {len(self.df)} networks remaining")
        
        if max_reactions is not None and 'n_reactions' in self.df.columns:
            self.df = self.df[self.df['n_reactions'] <= max_reactions]
            print(f"Max reactions <= {max_reactions}: {len(self.df)} networks remaining")
        
        if min_ercs is not None and 'n_ercs' in self.df.columns:
            self.df = self.df[self.df['n_ercs'] >= min_ercs]
            print(f"Min ERCs >= {min_ercs}: {len(self.df)} networks remaining")
        
        if max_ercs is not None and 'n_ercs' in self.df.columns:
            self.df = self.df[self.df['n_ercs'] <= max_ercs]
            print(f"Max ERCs <= {max_ercs}: {len(self.df)} networks remaining")
        
        filtered_out = initial_size - len(self.df)
        print(f"\nFiltered out {filtered_out} networks, {len(self.df)} remaining for analysis")
        
        # Update features after filtering
        self.prepare_features()
    
    def prepare_features(self):
        """Prepare feature groups for analysis"""
        # Create derived metrics
        if all(col in self.df.columns for col in ['n_levels', 'avg_branching']):
            self.df['hierarchy_complexity'] = self.df['n_levels'] * self.df['avg_branching']
        
        if all(col in self.df.columns for col in ['leaf_ratio', 'root_ratio']):
            self.df['structural_balance'] = self.df['leaf_ratio'] / (self.df['root_ratio'] + 1e-6)
        
        # Log transforms
        log_cols = ['n_species', 'n_reactions', 'n_ercs', 'ultra_optimized_time', 'optimized_time', 'original_time']
        for col in log_cols:
            if col in self.df.columns:
                self.df[f'log_{col}'] = np.log1p(self.df[col])
        
        # Define feature groups
        self.hierarchy_features = [
            'n_levels', 'max_chain_length', 'avg_branching', 'hierarchy_complexity',
            'avg_level_size', 'level_variance', 'hierarchy_density'
        ]
        
        self.connectivity_features = [
            'leaf_ratio', 'root_ratio', 'structural_balance', 'n_edges', 'avg_degree'
        ]
        
        self.performance_features = [
            'ultra_optimized_time', 'optimized_time', 'original_time',
            'log_ultra_optimized_time', 'log_optimized_time', 'log_original_time'
        ]
        
        # Filter to existing columns
        self.hierarchy_features = [f for f in self.hierarchy_features if f in self.df.columns]
        self.connectivity_features = [f for f in self.connectivity_features if f in self.df.columns]
        self.performance_features = [f for f in self.performance_features if f in self.df.columns]
        
        print(f"\nFeature groups prepared:")
        print(f"  Hierarchy: {len(self.hierarchy_features)} features")
        print(f"  Connectivity: {len(self.connectivity_features)} features")
        print(f"  Performance: {len(self.performance_features)} features")
    
    def analyze_hierarchy_structure(self):
        """Focused analysis of hierarchy structural properties"""
        print("\n" + "="*60)
        print("HIERARCHY STRUCTURE ANALYSIS")
        print("="*60)
        
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        fig.suptitle('ERC Hierarchy Structure Analysis', fontsize=16, fontweight='bold')
        
        # 1. Hierarchy depth distribution
        if 'n_levels' in self.df.columns:
            axes[0, 0].hist(self.df['n_levels'], bins=20, alpha=0.7, edgecolor='black', color='skyblue')
            axes[0, 0].set_xlabel('Number of Hierarchy Levels')
            axes[0, 0].set_ylabel('Frequency')
            axes[0, 0].set_title('Hierarchy Depth Distribution')
            axes[0, 0].grid(True, alpha=0.3)
        
        # 2. Chain length vs levels
        if 'max_chain_length' in self.df.columns and 'n_levels' in self.df.columns:
            axes[0, 1].scatter(self.df['n_levels'], self.df['max_chain_length'], alpha=0.6, color='coral')
            axes[0, 1].set_xlabel('Number of Levels')
            axes[0, 1].set_ylabel('Max Chain Length')
            axes[0, 1].set_title('Chain Length vs Hierarchy Depth')
            axes[0, 1].grid(True, alpha=0.3)
        
        # 3. Branching factor distribution
        if 'avg_branching' in self.df.columns:
            axes[0, 2].hist(self.df['avg_branching'], bins=20, alpha=0.7, edgecolor='black', color='lightgreen')
            axes[0, 2].set_xlabel('Average Branching Factor')
            axes[0, 2].set_ylabel('Frequency')
            axes[0, 2].set_title('Branching Distribution')
            axes[0, 2].grid(True, alpha=0.3)
        
        # 4. Structural balance
        if 'leaf_ratio' in self.df.columns and 'root_ratio' in self.df.columns:
            axes[1, 0].scatter(self.df['leaf_ratio'], self.df['root_ratio'], alpha=0.6, color='orchid')
            axes[1, 0].set_xlabel('Leaf Node Ratio')
            axes[1, 0].set_ylabel('Root Node Ratio')
            axes[1, 0].set_title('Structural Balance (Leaves vs Roots)')
            axes[1, 0].grid(True, alpha=0.3)
        
        # 5. Network size vs hierarchy complexity
        if 'log_n_ercs' in self.df.columns and 'hierarchy_complexity' in self.df.columns:
            axes[1, 1].scatter(self.df['log_n_ercs'], self.df['hierarchy_complexity'], alpha=0.6, color='gold')
            axes[1, 1].set_xlabel('Log(Number of ERCs)')
            axes[1, 1].set_ylabel('Hierarchy Complexity')
            axes[1, 1].set_title('Network Size vs Hierarchy Complexity')
            axes[1, 1].grid(True, alpha=0.3)
        
        # 6. Level size variance (hierarchy balance)
        if 'level_variance' in self.df.columns:
            axes[1, 2].hist(self.df['level_variance'], bins=20, alpha=0.7, edgecolor='black', color='lightcoral')
            axes[1, 2].set_xlabel('Level Size Variance')
            axes[1, 2].set_ylabel('Frequency')
            axes[1, 2].set_title('Hierarchy Balance')
            axes[1, 2].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('hierarchy_structure_analysis.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        # Print hierarchy statistics
        self.print_hierarchy_stats()
    
    def analyze_computational_performance(self):
        """Focused analysis of computational performance"""
        print("\n" + "="*60)
        print("COMPUTATIONAL PERFORMANCE ANALYSIS")
        print("="*60)
        
        # Find available algorithms
        alg_cols = {
            'Original': 'original_time',
            'Optimized': 'optimized_time', 
            'Ultra-optimized': 'ultra_optimized_time'
        }
        available_algs = {name: col for name, col in alg_cols.items() if col in self.df.columns}
        
        if not available_algs:
            print("No performance data available")
            return
        
        n_algs = len(available_algs)
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        fig.suptitle('Computational Performance Analysis', fontsize=16, fontweight='bold')
        
        # 1. Performance distribution comparison
        ax = axes[0, 0]
        perf_data = []
        labels = []
        for name, col in available_algs.items():
            data = self.df[col].dropna()
            if len(data) > 0:
                perf_data.append(data)
                labels.append(name)
        
        if perf_data:
            ax.boxplot(perf_data, labels=labels)
            ax.set_ylabel('Computation Time (seconds)')
            ax.set_title('Algorithm Performance Comparison')
            ax.set_yscale('log')
            ax.grid(True, alpha=0.3)
        
        # 2. Performance vs network size
        if 'log_n_ercs' in self.df.columns:
            ax = axes[0, 1]
            colors = ['red', 'orange', 'green']
            for i, (name, col) in enumerate(available_algs.items()):
                valid_data = self.df[self.df[col].notna()]
                if len(valid_data) > 0:
                    ax.scatter(valid_data['log_n_ercs'], valid_data[col], 
                              alpha=0.6, label=name, color=colors[i % len(colors)])
            ax.set_xlabel('Log(Number of ERCs)')
            ax.set_ylabel('Computation Time (seconds)')
            ax.set_title('Performance vs Network Size')
            ax.set_yscale('log')
            ax.legend()
            ax.grid(True, alpha=0.3)
        
        # 3. Speedup analysis (if multiple algorithms available)
        if len(available_algs) > 1:
            ax = axes[0, 2]
            # Find best and worst performing algorithms
            mean_times = {}
            for name, col in available_algs.items():
                mean_times[name] = self.df[col].dropna().mean()
            
            best_alg = min(mean_times.keys(), key=lambda x: mean_times[x])
            worst_alg = max(mean_times.keys(), key=lambda x: mean_times[x])
            
            best_col = available_algs[best_alg]
            worst_col = available_algs[worst_alg]
            
            # Calculate speedup with handling for infinite/very large values
            valid_mask = (self.df[best_col].notna()) & (self.df[worst_col].notna()) & (self.df[best_col] > 0)
            if valid_mask.any():
                speedup = self.df.loc[valid_mask, worst_col] / self.df.loc[valid_mask, best_col]
                # Filter out infinite or very large values
                finite_mask = speedup.abs() < 1e6  # Limit to reasonable speedup values
                speedup = speedup[finite_mask]
                
                if len(speedup) > 0:
                    ax.hist(speedup, bins=20, alpha=0.7, edgecolor='black', color='lightblue')
                    ax.set_xlabel('Speedup Factor')
                    ax.set_ylabel('Frequency')
                    ax.set_title(f'Speedup: {best_alg} vs {worst_alg}\n(filtered to reasonable values)')
                    ax.grid(True, alpha=0.3)
                else:
                    ax.text(0.5, 0.5, 'No valid speedup data\nafter filtering extremes',
                           ha='center', va='center', transform=ax.transAxes)
        
        # 4. Performance vs hierarchy complexity
        if 'hierarchy_complexity' in self.df.columns and available_algs:
            ax = axes[1, 0]
            best_alg_col = list(available_algs.values())[0]  # Use first available
            valid_data = self.df[self.df[best_alg_col].notna()]
            if len(valid_data) > 0:
                ax.scatter(valid_data['hierarchy_complexity'], valid_data[best_alg_col], 
                          alpha=0.6, color='purple')
                ax.set_xlabel('Hierarchy Complexity')
                ax.set_ylabel('Computation Time (seconds)')
                ax.set_title('Performance vs Structural Complexity')
                ax.set_yscale('log')
                ax.grid(True, alpha=0.3)
        
        # 5. Success rates
        success_cols = {'Original': 'success1', 'Optimized': 'success2', 'Ultra-optimized': 'success3'}
        available_success = {name: col for name, col in success_cols.items() if col in self.df.columns}
        
        if available_success:
            ax = axes[1, 1]
            success_rates = []
            names = []
            for name, col in available_success.items():
                rate = self.df[col].mean() * 100
                success_rates.append(rate)
                names.append(name)
            
            bars = ax.bar(names, success_rates, alpha=0.7, color=['red', 'orange', 'green'][:len(names)])
            ax.set_ylabel('Success Rate (%)')
            ax.set_title('Algorithm Success Rates')
            ax.set_ylim(0, 105)
            ax.grid(True, alpha=0.3)
            
            # Add value labels on bars
            for bar, rate in zip(bars, success_rates):
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2., height + 1,
                       f'{rate:.1f}%', ha='center', va='bottom')
        
        # 6. Performance scaling
        if 'log_n_ercs' in self.df.columns and available_algs:
            ax = axes[1, 2]
            best_alg_col = list(available_algs.values())[0]
            valid_data = self.df[self.df[best_alg_col].notna()]
            if len(valid_data) > 0:
                # Bin by network size and show performance distribution
                bins = pd.cut(valid_data['log_n_ercs'], bins=5)
                perf_by_size = [valid_data[valid_data.index.isin(bins[bins == bin_val].index)][best_alg_col].values 
                               for bin_val in bins.cat.categories]
                bin_labels = [f'{int(np.exp(cat.left))}-{int(np.exp(cat.right))}' for cat in bins.cat.categories]
                
                ax.boxplot([data for data in perf_by_size if len(data) > 0], 
                          labels=[label for i, label in enumerate(bin_labels) if len(perf_by_size[i]) > 0])
                ax.set_xlabel('ERC Size Range')
                ax.set_ylabel('Computation Time (seconds)')
                ax.set_title('Performance Scaling by Network Size')
                ax.set_yscale('log')
                plt.setp(ax.get_xticklabels(), rotation=45)
                ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('computational_performance_analysis.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        # Print performance statistics
        self.print_performance_stats(available_algs)
    
    def plot_correlation_matrix(self):
        """Create standalone correlation heatmap"""
        print("\n" + "="*60)
        print("CORRELATION MATRIX ANALYSIS")
        print("="*60)
        
        # Select variables for correlation
        all_vars = (self.hierarchy_features + self.connectivity_features + 
                   ['n_ercs', 'n_reactions', 'n_species'] +
                   [col for col in self.performance_features if 'log' not in col])
        available_vars = [v for v in all_vars if v in self.df.columns]
        
        if len(available_vars) < 3:
            print("Insufficient variables for correlation analysis")
            return
        
        plt.figure(figsize=(14, 12))
        
        corr_matrix = self.df[available_vars].corr()
        
        # Create heatmap
        sns.heatmap(corr_matrix, 
                   annot=True,
                   cmap='RdBu_r',
                   center=0,
                   square=True,
                   fmt='.2f',
                   cbar_kws={'shrink': 0.8},
                   xticklabels=True,
                   yticklabels=True)
        
        plt.title('ERC Hierarchy Variables - Correlation Matrix', fontsize=16, pad=20)
        plt.xticks(rotation=45, ha='right')
        plt.yticks(rotation=0)
        
        plt.tight_layout()
        plt.savefig('correlation_matrix.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        # Print strongest correlations
        self.print_strongest_correlations(corr_matrix)
    
    def analyze_performance_based_clustering(self):
        """Cluster networks based on best available performance algorithm"""
        print("\n" + "="*60)
        print("PERFORMANCE-BASED CLUSTERING ANALYSIS")
        print("="*60)
        
        # Find best performing algorithm
        perf_cols = ['ultra_optimized_time', 'optimized_time', 'original_time']
        best_perf_col = None
        for col in perf_cols:
            if col in self.df.columns and self.df[col].notna().sum() > 10:  # At least 10 data points
                best_perf_col = col
                break
        
        if not best_perf_col:
            print("No sufficient performance data for clustering")
            return None
        
        print(f"Using {best_perf_col.replace('_time', '')} algorithm for clustering")
        
        # Prepare clustering features
        cluster_features = self.hierarchy_features + self.connectivity_features
        cluster_features = [f for f in cluster_features if f in self.df.columns]
        
        if len(cluster_features) < 2:
            print("Insufficient features for clustering")
            return None
        
        # Filter to networks with performance data
        valid_data = self.df[self.df[best_perf_col].notna()].copy()
        
        # Prepare clustering data
        feature_data = valid_data[cluster_features].fillna(valid_data[cluster_features].median())
        scaler = RobustScaler()
        scaled_data = scaler.fit_transform(feature_data)
        
        # Find optimal number of clusters
        k_range = range(2, 5)
        silhouette_scores = []
        
        for k in k_range:
            kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
            labels = kmeans.fit_predict(scaled_data)
            score = silhouette_score(scaled_data, labels)
            silhouette_scores.append(score)
        
        optimal_k = k_range[np.argmax(silhouette_scores)]
        print(f"Optimal clusters: {optimal_k} (silhouette score: {max(silhouette_scores):.3f})")
        
        # Final clustering
        final_kmeans = KMeans(n_clusters=optimal_k, random_state=42, n_init=10)
        cluster_labels = final_kmeans.fit_predict(scaled_data)
        
        # Dimensionality reduction for visualization
        if HAS_UMAP:
            reducer = UMAP(n_components=2, random_state=42)
            viz_data = reducer.fit_transform(scaled_data)
            viz_method = "UMAP"
        else:
            reducer = TSNE(n_components=2, random_state=42, perplexity=min(30, len(scaled_data)//4))
            viz_data = reducer.fit_transform(scaled_data)
            viz_method = "t-SNE"
        
        # Create visualization
        fig, axes = plt.subplots(1, 2, figsize=(16, 6))
        fig.suptitle('Performance-Based Network Clustering', fontsize=16, fontweight='bold')
        
        # Cluster visualization
        scatter = axes[0].scatter(viz_data[:, 0], viz_data[:, 1], c=cluster_labels, 
                                 cmap='tab10', alpha=0.7, s=50)
        axes[0].set_xlabel(f'{viz_method} 1')
        axes[0].set_ylabel(f'{viz_method} 2')
        axes[0].set_title(f'Network Clusters ({viz_method})')
        axes[0].grid(True, alpha=0.3)
        
        # Performance by cluster
        cluster_performance = []
        cluster_names = []
        for cluster_id in range(optimal_k):
            cluster_mask = cluster_labels == cluster_id
            cluster_perf = valid_data.loc[valid_data.index[cluster_mask], best_perf_col]
            cluster_performance.append(cluster_perf)
            cluster_names.append(f'Cluster {cluster_id}\n(n={sum(cluster_mask)})')
        
        axes[1].boxplot(cluster_performance, labels=cluster_names)
        axes[1].set_ylabel('Computation Time (seconds)')
        axes[1].set_title('Performance Distribution by Cluster')
        axes[1].set_yscale('log')
        axes[1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('performance_clustering.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        # Store cluster results
        valid_data['cluster'] = cluster_labels
        self.cluster_results = {
            'data': valid_data,
            'labels': cluster_labels,
            'n_clusters': optimal_k,
            'performance_col': best_perf_col,
            'features': cluster_features
        }
        
        return self.cluster_results
    
    def analyze_performance_per_cluster(self, cluster_results=None):
        """Detailed analysis of performance characteristics per cluster"""
        if cluster_results is None:
            cluster_results = getattr(self, 'cluster_results', None)
        
        if cluster_results is None:
            print("No cluster results available. Run analyze_performance_based_clustering() first.")
            return
        
        print("\n" + "="*60)
        print("PERFORMANCE PER CLUSTER ANALYSIS")
        print("="*60)
        
        cluster_data = cluster_results['data']
        n_clusters = cluster_results['n_clusters']
        perf_col = cluster_results['performance_col']
        
        # Find all available performance columns
        all_perf_cols = ['original_time', 'optimized_time', 'ultra_optimized_time']
        available_perf_cols = [col for col in all_perf_cols if col in cluster_data.columns]
        
        if len(available_perf_cols) == 0:
            print("No performance data available")
            return
        
        n_perf_algs = len(available_perf_cols)
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        fig.suptitle('Performance Analysis by Cluster', fontsize=16, fontweight='bold')
        
        # 1. Performance comparison across clusters for each algorithm
        for i, perf_col in enumerate(available_perf_cols[:3]):
            ax = axes[0, i]
            cluster_perf = []
            cluster_names = []
            
            for cluster_id in range(n_clusters):
                cluster_mask = cluster_data['cluster'] == cluster_id
                perf_data = cluster_data.loc[cluster_mask, perf_col].dropna()
                if len(perf_data) > 0:
                    cluster_perf.append(perf_data)
                    cluster_names.append(f'C{cluster_id}\n(n={len(perf_data)})')
            
            if cluster_perf:
                ax.boxplot(cluster_perf, labels=cluster_names)
                ax.set_ylabel('Computation Time (seconds)')
                ax.set_title(f'{perf_col.replace("_time", "").title()} Performance')
                ax.set_yscale('log')
                ax.grid(True, alpha=0.3)
        
        # 2. Cluster characteristics comparison
        if len(cluster_results['features']) > 0:
            # Network size by cluster
            ax = axes[1, 0]
            if 'n_ercs' in cluster_data.columns:
                cluster_sizes = []
                cluster_names = []
                for cluster_id in range(n_clusters):
                    cluster_mask = cluster_data['cluster'] == cluster_id
                    size_data = cluster_data.loc[cluster_mask, 'n_ercs']
                    if len(size_data) > 0:
                        cluster_sizes.append(size_data)
                        cluster_names.append(f'Cluster {cluster_id}')
                
                ax.boxplot(cluster_sizes, labels=cluster_names)
                ax.set_ylabel('Number of ERCs')
                ax.set_title('Network Size by Cluster')
                ax.set_yscale('log')
                ax.grid(True, alpha=0.3)
            
            # Hierarchy complexity by cluster
            ax = axes[1, 1]
            if 'hierarchy_complexity' in cluster_data.columns:
                cluster_complexity = []
                for cluster_id in range(n_clusters):
                    cluster_mask = cluster_data['cluster'] == cluster_id
                    comp_data = cluster_data.loc[cluster_mask, 'hierarchy_complexity'].dropna()
                    if len(comp_data) > 0:
                        cluster_complexity.append(comp_data)
                
                if cluster_complexity:
                    ax.boxplot(cluster_complexity, labels=[f'C{i}' for i in range(len(cluster_complexity))])
                    ax.set_ylabel('Hierarchy Complexity')
                    ax.set_title('Structural Complexity by Cluster')
                    ax.grid(True, alpha=0.3)
            
            # Speedup analysis by cluster (if multiple algorithms available)
            if len(available_perf_cols) > 1:
                ax = axes[1, 2]
                best_col = available_perf_cols[0]  # Assume first is best
                worst_col = available_perf_cols[-1]  # Assume last is worst
                
                cluster_speedups = []
                for cluster_id in range(n_clusters):
                    cluster_mask = cluster_data['cluster'] == cluster_id
                    cluster_subset = cluster_data.loc[cluster_mask]
                    valid_mask = (cluster_subset[best_col].notna()) & (cluster_subset[worst_col].notna())
                    
                    if valid_mask.any():
                        speedup = cluster_subset.loc[valid_mask, worst_col] / cluster_subset.loc[valid_mask, best_col]
                        if len(speedup) > 0:
                            cluster_speedups.append(speedup)
                
                if cluster_speedups:
                    ax.boxplot(cluster_speedups, labels=[f'C{i}' for i in range(len(cluster_speedups))])
                    ax.set_ylabel('Speedup Factor')
                    ax.set_title('Algorithm Speedup by Cluster')
                    ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('performance_per_cluster.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        # Print cluster statistics
        self.print_cluster_statistics(cluster_data, n_clusters, available_perf_cols)
    
    def print_hierarchy_stats(self):
        """Print hierarchy structure statistics"""
        print("\nHierarchy Structure Statistics:")
        print("-" * 40)
        
        if 'n_levels' in self.df.columns:
            levels = self.df['n_levels'].dropna()
            print(f"Hierarchy Depth: {levels.mean():.2f} ± {levels.std():.2f} levels")
            print(f"Depth range: {levels.min():.0f} - {levels.max():.0f} levels")
        
        if 'max_chain_length' in self.df.columns:
            chains = self.df['max_chain_length'].dropna()
            print(f"Max Chain Length: {chains.mean():.2f} ± {chains.std():.2f}")
        
        if 'avg_branching' in self.df.columns:
            branching = self.df['avg_branching'].dropna()
            print(f"Average Branching: {branching.mean():.3f} ± {branching.std():.3f}")
    
    def print_performance_stats(self, available_algs):
        """Print performance statistics"""
        print("\nPerformance Statistics:")
        print("-" * 40)
        
        for name, col in available_algs.items():
            data = self.df[col].dropna()
            if len(data) > 0:
                print(f"{name}: {data.mean():.4f} ± {data.std():.4f} seconds")
                print(f"  Range: {data.min():.4f} - {data.max():.4f} seconds")
    
    def print_strongest_correlations(self, corr_matrix):
        """Print strongest correlations"""
        print("\nStrongest Correlations (|r| > 0.5):")
        print("-" * 50)
        
        correlations = []
        for i in range(len(corr_matrix.columns)):
            for j in range(i+1, len(corr_matrix.columns)):
                var1, var2 = corr_matrix.columns[i], corr_matrix.columns[j]
                corr_val = corr_matrix.iloc[i, j]
                if abs(corr_val) > 0.5:
                    correlations.append((abs(corr_val), corr_val, var1, var2))
        
        correlations.sort(reverse=True)
        
        if correlations:
            for abs_corr, corr, var1, var2 in correlations[:10]:  # Top 10
                direction = "+" if corr > 0 else "-"
                print(f"{var1} ↔ {var2}: {corr:.3f} ({direction})")
        else:
            print("No strong correlations found")
    
    def print_cluster_statistics(self, cluster_data, n_clusters, perf_cols):
        """Print detailed cluster statistics"""
        print("\nCluster Statistics:")
        print("=" * 50)
        
        for cluster_id in range(n_clusters):
            cluster_mask = cluster_data['cluster'] == cluster_id
            cluster_subset = cluster_data.loc[cluster_mask]
            
            print(f"\nCluster {cluster_id} (n={len(cluster_subset)}):")
            print("-" * 30)
            
            # Performance stats
            for col in perf_cols:
                perf_data = cluster_subset[col].dropna()
                if len(perf_data) > 0:
                    alg_name = col.replace('_time', '')
                    print(f"  {alg_name}: {perf_data.mean():.4f} ± {perf_data.std():.4f}s")
            
            # Network characteristics
            if 'n_ercs' in cluster_subset.columns:
                ercs = cluster_subset['n_ercs'].dropna()
                if len(ercs) > 0:
                    print(f"  ERCs: {ercs.mean():.1f} ± {ercs.std():.1f}")
            
            if 'n_levels' in cluster_subset.columns:
                levels = cluster_subset['n_levels'].dropna()
                if len(levels) > 0:
                    print(f"  Levels: {levels.mean():.1f} ± {levels.std():.1f}")
    
    def run_complete_analysis(self, **size_filters):
        """Run all analyses in sequence"""
        print("="*80)
        print("COMPLETE ERC HIERARCHY ANALYSIS")
        print("="*80)
        
        # Apply filters if provided
        if size_filters:
            self.apply_size_filters(**size_filters)
        
        # Run all analyses
        print("\n1. Hierarchy Structure Analysis...")
        self.analyze_hierarchy_structure()
        
        print("\n2. Computational Performance Analysis...")
        self.analyze_computational_performance()
        
        print("\n3. Correlation Matrix...")
        self.plot_correlation_matrix()
        
        print("\n4. Performance-Based Clustering...")
        cluster_results = self.analyze_performance_based_clustering()
        
        if cluster_results:
            print("\n5. Performance Per Cluster...")
            self.analyze_performance_per_cluster(cluster_results)
        
        print("\n" + "="*80)
        print("ANALYSIS COMPLETE - Generated visualizations:")
        print("  • hierarchy_structure_analysis.png")
        print("  • computational_performance_analysis.png") 
        print("  • correlation_matrix.png")
        print("  • performance_clustering.png")
        print("  • performance_per_cluster.png")
        print("="*80)

def main():
    """Main analysis function with interactive options"""
    import sys
    
    # Load data
    if len(sys.argv) == 2:
        data_file = sys.argv[1]
    else:
        default_file = './benchmark_results/erc_hierarchy_benchmark.csv'
        if Path(default_file).exists():
            data_file = default_file
            print(f"Using default file: {data_file}")
        else:
            print("Usage: python erc_hierarchy_analysis.py <erc_hierarchy_benchmark.csv>")
            return
    
    try:
        analyzer = ModularERCAnalyzer(data_file)
        
        # Example usage options:
        print("\n" + "="*60)
        print("ANALYSIS OPTIONS")
        print("="*60)
        print("1. Complete analysis (all networks):")
        print("   analyzer.run_complete_analysis()")
        print("\n2. Filtered analysis (e.g., small networks):")
        print("   analyzer.run_complete_analysis(max_ercs=100, max_reactions=200)")
        print("\n3. Individual analyses:")
        print("   analyzer.analyze_hierarchy_structure()")
        print("   analyzer.analyze_computational_performance()")
        print("   analyzer.plot_correlation_matrix()")
        print("   analyzer.analyze_performance_based_clustering()")
        print("   analyzer.analyze_performance_per_cluster()")
        
        # Run complete analysis by default
        print("\nRunning complete analysis...")
        analyzer.run_complete_analysis(min_ercs=4)

    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()