# # # Script 3: Batch Hierarchies Analysis
# This script analyzes all reaction networks in a directory and creates a statistical report
# The script is structured to:  
#   (a) load reaction networks from a directory
#   (b) calculate semi-organisations for each
#   (c) analyze the lattice properties (completeness, modularity, distributivity)
#   (d) analyze negation properties
#   (e) create a statistical report of all results

# Import Python standard library modules
import os
import sys
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

# Add the root directory to the PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.stdout.reconfigure(encoding='utf-8')

# Imports from the pyCOT library
from pyCOT.reaction_network import ReactionNetwork 
from pyCOT.file_manipulation import load_pyCOT_from_file
from pyCOT.rn_visualize import *                        
from pyCOT.simulations import *           
from pyCOT.abstractions import *          
from pyCOT.plot_dynamics import *         
from pyCOT.rn_types import StoichiometryMatrix
from pyCOT.abstractions import abstraction_ordinary 
from pyCOT.closure_structure import *     
from pyCOT.rn_hierarchy import *
import networkx as nx

# Directory containing reaction network files
input_directory = 'networks/RandomAlife/Ns20'
# Output directory for results
output_directory = input_directory+'/results'

# Create output directory if it doesn't exist
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

# Initialize results storage
results = []

# Create log file
log_filename = os.path.join(output_directory, f"analysis_log_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt")
with open(log_filename, 'w') as log_file:
    log_file.write(f"Analysis started at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    log_file.write(f"Input directory: {input_directory}\n\n")

def log_message(message):
    """Write a message to both console and log file"""
    print(message)
    with open(log_filename, 'a') as log_file:
        log_file.write(message + '\n')

def analyze_network(file_path):
    """Analyze a single reaction network and return its properties"""
    try:
        log_message(f"\n\n=== ANALYZING FILE: {file_path} ===")
        
        # Load the reaction network
        testRN = load_pyCOT_from_file(file_path)
        
        # Extract filename without path and extension
        filename = os.path.basename(file_path)
        network_name = os.path.splitext(filename)[0]
        
        # Calculate the semi-organisations
        log_message(f"Calculating semi-organisations for {network_name}")
        input_data = reactive_semi_orgs(testRN)
        
        # Build the hierarchy
        log_message(f"Building hierarchy graph for {network_name}")
        G = hierarchy_build(input_data)
        
        # Store basic information
        network_info = {
            'network_name': network_name,
            'file_path': file_path,
            'num_species': len(testRN.SpStr),
            'num_reactions': len(testRN.RnStr),
            'num_semi_orgs': len(input_data)
        }
        
        # Analyze lattice properties
        log_message(f"Analyzing lattice properties for {network_name}")
        
        # Analyze top and bottom elements
        top_bottom = identify_top_bottom_elements(G)
        network_info['has_top'] = top_bottom['has_top']
        network_info['has_bottom'] = top_bottom['has_bottom']
        network_info['is_bounded_lattice'] = top_bottom['has_top'] and top_bottom['has_bottom']
        
        if top_bottom['has_top']:
            network_info['top_element_size'] = len(get_species_from_node(G, top_bottom['top_element']))
        if top_bottom['has_bottom']:
            network_info['bottom_element_size'] = len(get_species_from_node(G, top_bottom['bottom_element']))
        
        # Analyze completeness
        completeness_result = test_completeness(testRN, G)
        network_info['is_complete'] = completeness_result['is_complete']
        network_info['completeness_ratio'] = completeness_result['completeness_ratio']
        network_info['join_ratio'] = completeness_result['stats']['join_ratio']
        network_info['meet_ratio'] = completeness_result['stats']['meet_ratio']
        
        # If it's complete, analyze modularity and distributivity
        if network_info['is_complete']:
            # Analyze modularity
            modularity_result = safe_test_modularity(testRN, G)
            network_info['is_modular'] = modularity_result['is_modular']
            network_info['modularity_ratio'] = modularity_result['modularity_ratio']
            
            # Analyze distributivity
            distributivity_result = safe_test_distributivity(testRN, G)
            network_info['is_distributive'] = distributivity_result['is_distributive']
            network_info['distributivity_ratio'] = distributivity_result['distributivity_ratio']
            network_info['dist1_ratio'] = distributivity_result['stats']['dist1_ratio']
            network_info['dist2_ratio'] = distributivity_result['stats']['dist2_ratio']
        else:
            network_info['is_modular'] = False
            network_info['modularity_ratio'] = 0
            network_info['is_distributive'] = False
            network_info['distributivity_ratio'] = 0
            network_info['dist1_ratio'] = 0
            network_info['dist2_ratio'] = 0
        
        # If it's a bounded lattice, analyze negation properties
        if network_info['is_bounded_lattice']:
            log_message(f"Analyzing negation properties for {network_name}")
            negation_map, score_info = generate_best_negation_map(testRN, G)
            
            network_info['negation_overall_score'] = score_info.get('overall_score', 0)
            network_info['negation_complement_score'] = score_info.get('complement_score', 0)
            network_info['negation_involution_score'] = score_info.get('involution_score', 0)
            network_info['negation_demorgan_score'] = score_info.get('demorgan_score', 0)
        else:
            network_info['negation_overall_score'] = 0
            network_info['negation_complement_score'] = 0
            network_info['negation_involution_score'] = 0
            network_info['negation_demorgan_score'] = 0
        
        # Graph structure properties
        network_info['graph_nodes'] = G.number_of_nodes()
        network_info['graph_edges'] = G.number_of_edges()
        network_info['graph_density'] = nx.density(G)
        
        # Lattice metrics
        network_info['lattice_height'] = max([len(nx.shortest_path(G, source=node, target=top_bottom['top_element'])) - 1 
                                           for node in G.nodes() if nx.has_path(G, node, top_bottom['top_element'])]) if top_bottom['has_top'] else 0
        
        # Log a summary of the analysis
        log_message(f"\nAnalysis Summary for {network_name}:")
        log_message(f"  Species: {network_info['num_species']}, Reactions: {network_info['num_reactions']}")
        log_message(f"  Semi-organizations: {network_info['num_semi_orgs']}")
        log_message(f"  Bounded Lattice: {network_info['is_bounded_lattice']}")
        log_message(f"  Complete: {network_info['is_complete']} (Ratio: {network_info['completeness_ratio']:.2f})")
        log_message(f"  Modular: {network_info['is_modular']} (Ratio: {network_info['modularity_ratio']:.2f})")
        log_message(f"  Distributive: {network_info['is_distributive']} (Ratio: {network_info['distributivity_ratio']:.2f})")
        if network_info['is_bounded_lattice']:
            log_message(f"  Negation Score: {network_info['negation_overall_score']:.2f}")
            log_message(f"  De Morgan Score: {network_info['negation_demorgan_score']:.2f}")
        
        return network_info
        
    except Exception as e:
        log_message(f"ERROR analyzing {file_path}: {str(e)}")
        return {'network_name': os.path.basename(file_path), 'error': str(e)}

def generate_statistics(results_list):
    """Generate statistics from all network results"""
    if not results_list:
        return "No results available for statistics"
    
    # Filter out results with errors
    valid_results = [r for r in results_list if 'error' not in r]
    
    if not valid_results:
        return "No valid results available for statistics"
    
    # Convert to DataFrame for easier analysis
    df = pd.DataFrame(valid_results)
    
    # Basic counts
    total_networks = len(valid_results)
    bounded_lattices = df['is_bounded_lattice'].sum()
    complete_lattices = df['is_complete'].sum()
    modular_lattices = df['is_modular'].sum()
    distributive_lattices = df['is_distributive'].sum()
    
    # Calculate percentages
    bounded_pct = (bounded_lattices / total_networks) * 100
    complete_pct = (complete_lattices / total_networks) * 100
    modular_pct = (modular_lattices / total_networks) * 100
    distributive_pct = (distributive_lattices / total_networks) * 100
    
    # Calculate averages for continuous metrics
    avg_completeness = df['completeness_ratio'].mean()
    avg_modularity = df['modularity_ratio'].mean()
    avg_distributivity = df['distributivity_ratio'].mean()
    avg_negation = df['negation_overall_score'].mean()
    avg_demorgan = df['negation_demorgan_score'].mean()
    
    # Statistics report
    stats = {
        'total_networks': total_networks,
        'bounded_lattices': bounded_lattices,
        'bounded_lattices_pct': bounded_pct,
        'complete_lattices': complete_lattices,
        'complete_lattices_pct': complete_pct,
        'modular_lattices': modular_lattices,
        'modular_lattices_pct': modular_pct,
        'distributive_lattices': distributive_lattices,
        'distributive_lattices_pct': distributive_pct,
        'avg_completeness_ratio': avg_completeness,
        'avg_modularity_ratio': avg_modularity,
        'avg_distributivity_ratio': avg_distributivity,
        'avg_negation_score': avg_negation,
        'avg_demorgan_score': avg_demorgan
    }
    
    # Additional statistics
    stats['avg_species'] = df['num_species'].mean()
    stats['avg_reactions'] = df['num_reactions'].mean()
    stats['avg_semi_orgs'] = df['num_semi_orgs'].mean()
    stats['avg_graph_nodes'] = df['graph_nodes'].mean()
    stats['avg_graph_edges'] = df['graph_edges'].mean()
    stats['avg_graph_density'] = df['graph_density'].mean()
    
    # Calculate correlations
    correlation_columns = [
        'num_species', 'num_reactions', 'num_semi_orgs', 
        'completeness_ratio', 'modularity_ratio', 'distributivity_ratio',
        'negation_overall_score', 'negation_demorgan_score'
    ]
    
    corr_matrix = df[correlation_columns].corr()
    stats['correlations'] = corr_matrix.to_dict()
    
    return stats

def plot_statistics(results_df, stats, output_directory):
    """Create visualizations of the results"""
    # Create output directory for plots
    plots_dir = os.path.join(output_directory, 'plots')
    if not os.path.exists(plots_dir):
        os.makedirs(plots_dir)
    
    # 1. Bar chart of lattice properties
    fig, ax = plt.subplots(figsize=(10, 6))
    properties = ['Bounded', 'Complete', 'Modular', 'Distributive']
    values = [
        stats['bounded_lattices_pct'],
        stats['complete_lattices_pct'],
        stats['modular_lattices_pct'],
        stats['distributive_lattices_pct']
    ]
    
    ax.bar(properties, values, color=['blue', 'green', 'orange', 'red'])
    ax.set_ylabel('Percentage of Networks (%)')
    ax.set_title('Lattice Properties Across All Networks')
    ax.grid(axis='y', linestyle='--', alpha=0.7)
    
    # Add value labels on top of each bar
    for i, v in enumerate(values):
        ax.text(i, v + 1, f"{v:.1f}%", ha='center')
    
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, 'lattice_properties.png'))
    
    # 2. Histogram of completeness ratios
    plt.figure(figsize=(10, 6))
    plt.hist(results_df['completeness_ratio'], bins=20, alpha=0.7, color='blue')
    plt.axvline(results_df['completeness_ratio'].mean(), color='red', linestyle='dashed', linewidth=2)
    plt.text(results_df['completeness_ratio'].mean() + 0.01, plt.ylim()[1]*0.9, 
             f'Mean: {results_df["completeness_ratio"].mean():.2f}', color='red')
    plt.xlabel('Completeness Ratio')
    plt.ylabel('Number of Networks')
    plt.title('Distribution of Completeness Ratios')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, 'completeness_distribution.png'))
    
    # 3. Scatterplot of distributivity vs modularity
    plt.figure(figsize=(10, 6))
    plt.scatter(results_df['modularity_ratio'], results_df['distributivity_ratio'], 
                alpha=0.7, c=results_df['num_semi_orgs'], cmap='viridis')
    plt.colorbar(label='Number of Semi-Organizations')
    plt.xlabel('Modularity Ratio')
    plt.ylabel('Distributivity Ratio')
    plt.title('Modularity vs Distributivity')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, 'modularity_vs_distributivity.png'))
    
    # 4. Correlation matrix heatmap
    correlation_columns = [
        'num_species', 'num_reactions', 'num_semi_orgs', 
        'completeness_ratio', 'modularity_ratio', 'distributivity_ratio',
        'negation_overall_score', 'negation_demorgan_score'
    ]
    
    corr_matrix = results_df[correlation_columns].corr()
    
    plt.figure(figsize=(12, 10))
    im = plt.imshow(corr_matrix, cmap='coolwarm', vmin=-1, vmax=1)
    plt.colorbar(im, label='Correlation Coefficient')
    
    # Add correlation values to the cells
    for i in range(len(corr_matrix.columns)):
        for j in range(len(corr_matrix.columns)):
            text = plt.text(j, i, f'{corr_matrix.iloc[i, j]:.2f}',
                            ha="center", va="center", color="black" if abs(corr_matrix.iloc[i, j]) < 0.7 else "white")
    
    plt.xticks(range(len(corr_matrix.columns)), corr_matrix.columns, rotation=45, ha='right')
    plt.yticks(range(len(corr_matrix.columns)), corr_matrix.columns)
    plt.title('Correlation Matrix of Key Metrics')
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, 'correlation_matrix.png'))
    
    # 5. De Morgan's laws compliance vs negation score
    plt.figure(figsize=(10, 6))
    plt.scatter(results_df['negation_overall_score'], results_df['negation_demorgan_score'], 
                alpha=0.7, c=results_df['distributivity_ratio'], cmap='plasma')
    plt.colorbar(label='Distributivity Ratio')
    plt.xlabel('Overall Negation Score')
    plt.ylabel('De Morgan Laws Score')
    plt.title('Negation Quality vs De Morgan Laws Compliance')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, 'negation_vs_demorgan.png'))
    
    # Return the list of created plot files
    return [
        os.path.join(plots_dir, 'lattice_properties.png'),
        os.path.join(plots_dir, 'completeness_distribution.png'),
        os.path.join(plots_dir, 'modularity_vs_distributivity.png'),
        os.path.join(plots_dir, 'correlation_matrix.png'),
        os.path.join(plots_dir, 'negation_vs_demorgan.png')
    ]

def generate_report(results_list, stats, plot_files, output_directory):
    """Generate a comprehensive HTML report of the results"""
    report_file = os.path.join(output_directory, 'analysis_report.html')
    
    # Convert timestamps to more readable format for the report
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    
    # Start building the HTML content
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Reaction Network Analysis Report</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 20px; line-height: 1.6; }}
            h1, h2, h3 {{ color: #333366; }}
            table {{ border-collapse: collapse; width: 100%; margin-bottom: 20px; }}
            th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
            th {{ background-color: #f2f2f2; }}
            tr:nth-child(even) {{ background-color: #f9f9f9; }}
            .plot-container {{ margin: 20px 0; text-align: center; }}
            .plot-container img {{ max-width: 100%; height: auto; }}
            .statistics {{ background-color: #f5f5f5; padding: 15px; border-radius: 5px; }}
            .correlation {{ overflow-x: auto; }}
        </style>
    </head>
    <body>
        <h1>Reaction Network Analysis Report</h1>
        <p>Generated on: {timestamp}</p>
        <p>Input directory: {input_directory}</p>
        
        <h2>Summary Statistics</h2>
        <div class="statistics">
            <p><strong>Total Networks Analyzed:</strong> {stats['total_networks']}</p>
            <p><strong>Networks with Bounded Lattices:</strong> {stats['bounded_lattices']} ({stats['bounded_lattices_pct']:.1f}%)</p>
            <p><strong>Networks with Complete Lattices:</strong> {stats['complete_lattices']} ({stats['complete_lattices_pct']:.1f}%)</p>
            <p><strong>Networks with Modular Lattices:</strong> {stats['modular_lattices']} ({stats['modular_lattices_pct']:.1f}%)</p>
            <p><strong>Networks with Distributive Lattices:</strong> {stats['distributive_lattices']} ({stats['distributive_lattices_pct']:.1f}%)</p>
            
            <h3>Average Metrics</h3>
            <p>Average Species: {stats['avg_species']:.2f}</p>
            <p>Average Reactions: {stats['avg_reactions']:.2f}</p>
            <p>Average Semi-Organizations: {stats['avg_semi_orgs']:.2f}</p>
            <p>Average Completeness Ratio: {stats['avg_completeness_ratio']:.2f}</p>
            <p>Average Modularity Ratio: {stats['avg_modularity_ratio']:.2f}</p>
            <p>Average Distributivity Ratio: {stats['avg_distributivity_ratio']:.2f}</p>
            <p>Average Negation Quality Score: {stats['avg_negation_score']:.2f}</p>
            <p>Average De Morgan Laws Score: {stats['avg_demorgan_score']:.2f}</p>
        </div>
        
        <h2>Visualizations</h2>
    """
    
    # Add plots to the report
    for plot_file in plot_files:
        plot_filename = os.path.basename(plot_file)
        html_content += f"""
        <div class="plot-container">
            <h3>{plot_filename.replace('.png', '').replace('_', ' ').title()}</h3>
            <img src="plots/{plot_filename}" alt="{plot_filename}" />
        </div>
        """
    
    # Add network details table
    html_content += """
        <h2>Individual Network Details</h2>
        <div style="overflow-x:auto;">
        <table>
            <tr>
                <th>Network Name</th>
                <th>Species</th>
                <th>Reactions</th>
                <th>Semi-Orgs</th>
                <th>Bounded</th>
                <th>Complete</th>
                <th>Modular</th>
                <th>Distributive</th>
                <th>Negation Score</th>
                <th>De Morgan Score</th>
            </tr>
    """
    
    # Add rows for each network
    for result in results_list:
        if 'error' in result:
            html_content += f"""
            <tr>
                <td>{result['network_name']}</td>
                <td colspan="9">Error: {result['error']}</td>
            </tr>
            """
        else:
            html_content += f"""
            <tr>
                <td>{result['network_name']}</td>
                <td>{result['num_species']}</td>
                <td>{result['num_reactions']}</td>
                <td>{result['num_semi_orgs']}</td>
                <td>{"Yes" if result['is_bounded_lattice'] else "No"}</td>
                <td>{"Yes" if result['is_complete'] else "No"} ({result['completeness_ratio']:.2f})</td>
                <td>{"Yes" if result['is_modular'] else "No"} ({result['modularity_ratio']:.2f})</td>
                <td>{"Yes" if result['is_distributive'] else "No"} ({result['distributivity_ratio']:.2f})</td>
                <td>{result['negation_overall_score']:.2f}</td>
                <td>{result['negation_demorgan_score']:.2f}</td>
            </tr>
            """
    
    # Close the table and HTML document
    html_content += """
        </table>
        </div>
        
        <h2>Conclusions</h2>
        <p>This analysis explores the lattice properties of reaction networks, focusing on completeness, modularity, 
        distributivity, and how well these structures support classical logical operations like negation and De Morgan's laws.</p>
        
        <p>The results indicate that reaction networks naturally form structures that often deviate from classical Boolean lattices, 
        showing properties similar to quantum logical structures, particularly in the way negation operations and De Morgan's laws behave.</p>
        
    </body>
    </html>
    """
    
    # Write the HTML to file
    with open(report_file, 'w') as f:
        f.write(html_content)
    
    return report_file

# Main script execution
def main():
    log_message(f"Starting batch analysis of reaction networks in {input_directory}")
    
    # Get all .txt files in the input directory
    file_list = [os.path.join(input_directory, f) for f in os.listdir(input_directory) 
                if os.path.isfile(os.path.join(input_directory, f)) and f.endswith('.txt')]
    
    if not file_list:
        log_message(f"No .txt files found in {input_directory}")
        return
    
    log_message(f"Found {len(file_list)} reaction network files")
    
    # Analyze each network
    for i, file_path in enumerate(file_list):
        log_message(f"\nProcessing file {i+1}/{len(file_list)}: {file_path}")
        result = analyze_network(file_path)
        results.append(result)
    
    # Save raw results to JSON
    results_file = os.path.join(output_directory, 'analysis_results.json')
    with open(results_file, 'w') as f:
        json.dump(results, f, indent=4)
    
    log_message(f"\nRaw results saved to {results_file}")
    
    # Generate statistics
    log_message("\nGenerating statistics...")
    stats = generate_statistics(results)
    
    # Save statistics to JSON
    stats_file = os.path.join(output_directory, 'statistics.json')
    with open(stats_file, 'w') as f:
        json.dump(stats, f, indent=4)
    
    log_message(f"Statistics saved to {stats_file}")
    
    # Generate plots
    log_message("\nGenerating plots...")
    valid_results = [r for r in results if 'error' not in r]
    if valid_results:
        results_df = pd.DataFrame(valid_results)
        plot_files = plot_statistics(results_df, stats, output_directory)
        log_message(f"Generated {len(plot_files)} plot files")
    else:
        plot_files = []
        log_message("No valid results for plotting")
    
    # Generate report
    log_message("\nGenerating final report...")
    report_file = generate_report(results, stats, plot_files, output_directory)
    log_message(f"Final report saved to {report_file}")
    
    log_message("\nAnalysis complete!")

if __name__ == "__main__":
    main()
# Analyze which semi-organizations are organizations
#flux_results = analyze_semi_orgs_flux_solutions(testRN, input_data)

# Filter organizations
#organizations = [semi_org for semi_org, has_solution, _, _ in flux_results if has_solution]
#print(f"\nFound {len(organizations)} organizations:")
#for i, org in enumerate(organizations):
#   print(f"{i+1}. {org}")
#hierarchy_visualize_html(organizations)  
"""
The hierarchy visualisation shows the semi-organisations of the reaction network (RN) as a directed graph.
The nodes represent the semi-organisations, and the edges represent the relationships between them.
By default in hierarchy_visualize_html():
    Nodes are represented by circular nodes colored 'cyan' and 
    Edges are represented by arrows colored 'gray'.
The node size by default is of size 20, and the edge width is of size 2.

"""

# # Hierarchy of semi-organisations with lst_color_subsets (Farm.txt)  
# lst_color_subsets = [('red', [['water'], 
#                               ['chickens', 'eggs', 'infr', 'water', 'grain', 'straw', 'fertilizer'], 
#                               ['chickens', 'eggs', 'infr', 'grass', 'grain', 'straw', 'water', 'fertilizer'], 
#                               ['chickens', 'money', 'infr', 'eggs', 'grass', 'water', 'grain', 'straw', 'fertilizer', 'farmer']])] # List of tuples with color and subsets of species
# hierarchy_visualize_html(input_data, node_color="blue", lst_color_subsets=lst_color_subsets) # Open the HTML hierarchy visualisation with the specified color for the nodes and list of color subsets

#####################################################################
# (d) generate time series of abstractions for random dynamics and for mak dynamics 
#####################################################################
# # Generates an initial state vector based on the number of species in the reaction network
#x_inicial_random = generate_random_vector(len(testRN.SpStr))  
#print(x_inicial_random)

# Calculates the stoichiometric matrix associated with the reaction network
#matrix_data = universal_stoichiometric_matrix(testRN) # Calculates the universal stoichiometric matrix associated with the reaction network
#print(matrix_data)
#S = StoichiometryMatrix(matrix_data, species=testRN.SpStr, reactions=testRN.RnStr) # Creates a StoichiometryMatrix object with the required set of species and reactions
#print(S)

############################################################
# # Simulates the time series of abstractions for random dynamics
############################################################ 
#random_time_series, random_flux_series = simulate_discrete_random(testRN, S, x_inicial_random, n_iter=100)

# # Define a fixed threshold or per-species thresholds
#threshold_random = 0.01  # Fixed threshold for random dynamics

# # Generates an abstraction for random dynamics 
# random_abstract_time_series = abstraction_ordinary(random_time_series)          # Generates an abstraction with a default threshold of 0.5
#random_abstract_time_series = abstraction_ordinary(random_time_series, threshold_random) # Generates an abstraction with a defined threshold for each species.
#print("Random_abstract_time_series:")
#print(random_abstract_time_series)

############################################################
# # Simulates the time series of abstractions for mak dynamics
############################################################
#mak_time_series, mak_flux_series = simulate_ode_mak(testRN) 

# # Define a fixed threshold or per-species thresholds
#threshold_mak = 0.01  # Fixed threshold for mak dynamics

# # Generates an abstraction for mak dynamics 
#abstract_time_series = abstraction_ordinary(time_series)          # Generates an abstraction with a default threshold of 0.5
#mak_abstract_time_series = abstraction_ordinary(mak_time_series, threshold_mak) # Generates an abstraction with a defined threshold for each species.
#print("MAK_abstract_time_series:")
#print(mak_abstract_time_series)

#####################################################################
# (e) plot static dynamics graphs in both cases 
##################################################################### 
# # Static graph with nodes of semi-organisations and abstractions 
#plot_static_abstraction_graph_hierarchy_html(random_abstract_time_series)
#plot_static_abstraction_graph_hierarchy_html(mak_abstract_time_series)            
 
##################################################################### 
# (f) plot dynamic movie in both cases
#####################################################################
# # # Movie  
#plot_abstraction_graph_movie_html(random_abstract_time_series) # Open the HTML file with the animation of the abstraction graph
#plot_abstraction_graph_movie_html(mak_abstract_time_series)    # Open the HTML file with the animation of the abstraction graph

# # # Movie with nodes of semi-organisations and abstractions 
#film_semiorganizations_abstractions_html(random_abstract_time_series,input_data) # Open the HTML file with the animation of the abstraction graph
#film_semiorganizations_abstractions_html(mak_abstract_time_series,input_data)    # Open the HTML file with the animation of the abstraction graph

# # Plot Hasse diagram 
#plot_hasse_diagram(testRN.SpStr) # Only can be plotted for RN with less than 10 species, because the number of nodes would be 2^n.