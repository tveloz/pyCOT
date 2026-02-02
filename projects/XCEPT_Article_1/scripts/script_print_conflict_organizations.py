import os
import sys
# Add the project root directory to the PYTHONPATH (go up 4 levels from scripts/XCEPT_project/Article_1/)
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))

from pyCOT.io.functions import read_txt
from pyCOT.analysis.Persistent_Modules_Generator import compute_all_organizations
from pyCOT.visualization.rn_visualize import hierarchy_visualize_html
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from collections import defaultdict
import math

def analyze_organizations_with_conflict(file_path):
    """
    Analyze organizations and visualize them colored by conflict species presence.
    
    Parameters:
    -----------
    file_path : str
        Path to the reaction network file
    """
    print(f"Analyzing network: {file_path}")
    print("=" * 60)
    
    # Check if file exists
    if not os.path.exists(file_path):
        print(f"ERROR: File not found: {file_path}")
        print("Please provide a valid network file path.")
        return None, None, None
    
    # Load the reaction network
    print(f"Loading network from: {file_path}")
    rn = read_txt(file_path)
    
    print(f"Network loaded: {len(rn.species())} species, {len(rn.reactions())} reactions")
    
    # Compute all organizations
    print("Computing organizations...")
    results = compute_all_organizations(
        rn,
        max_generator_size=8,
        max_organization_size=5,
        verbose=False
    )
    
    # Get organizations
    elementary_sos = results['elementary_sos']
    all_orgs = results['all_organizations']
    
    # Define conflict species (customize this based on your model)
    conflict_species = {
        'V', 'AG_SL', 'AG_RL', 'WR_SL', 'WR_RL',  # Violence and armed groups
    }
    
    # Track unique organizations
    organizations = []
    seen_closures = set()
    
    # Process organizations
    for idx, org in enumerate(elementary_sos + all_orgs):
        # Get species list
        if hasattr(org, "combined_closure"):  # Non-elementary
            species_list = [sp.name for sp in org.combined_closure]
        elif hasattr(org, "closure_species"):  # Elementary
            species_list = [sp.name for sp in org.closure_species]
        else:
            continue
            
        species_set = frozenset(species_list)
        
        # Skip duplicates
        if species_set in seen_closures:
            continue
        seen_closures.add(species_set)
        
        # Count conflict species
        conflict_count = sum(1 for sp in species_list if sp in conflict_species)
        total_species = len(species_list)
        
        # Calculate conflict level (0 to 1)
        conflict_level = conflict_count / total_species if total_species > 0 else 0
        
        organizations.append({
            'index': len(organizations),
            'species_list': species_list,
            'species_set': set(species_list),
            'size': total_species,
            'conflict_count': conflict_count,
            'conflict_level': conflict_level,
            'has_conflict': conflict_count > 0,
            'is_elementary': idx < len(elementary_sos)
        })
    
    print(f"\nFound {len(organizations)} unique organizations")
    return organizations, conflict_species, rn

def create_conflict_colored_visualization(organizations, conflict_species, filename="conflict_organizations.html"):
    """
    Create a hierarchy visualization with organizations colored by conflict level
    and sized by number of species.
    
    Parameters:
    -----------
    organizations : list
        List of organization dictionaries
    conflict_species : set
        Set of conflict species names
    filename : str
        Output HTML filename
    """
    print(f"\nCreating visualization: {filename}")
    
    # Prepare data for hierarchy visualization
    org_sets = [org['species_set'] for org in organizations]
    
    # Calculate size range for scaling
    sizes = [org['size'] for org in organizations]
    min_size = min(sizes) if sizes else 1
    max_size = max(sizes) if sizes else 1
    
    # Create colormap from green (no conflict) to red (full conflict)
    cmap = plt.cm.RdYlGn_r  # Red-Yellow-Green, reversed (red=high conflict)
    
    # Create visualization manually using pyvis (not using hierarchy_visualize_html)
    # because we need individual node sizes
    from pyvis.network import Network
    
    # Initialize the network
    net = Network(height="750px", width="100%", directed=True, notebook=False)
    
    # Configure layout for hierarchy
    net.set_options(f"""
    {{
      "nodes": {{
        "font": {{
          "size": 14,
          "align": "center"
        }},
        "borderWidth": 2,
        "borderColor": "black"
      }},
      "edges": {{
        "smooth": false,
        "color": "gray",
        "width": 2
      }},
      "physics": {{
        "enabled": false,
        "stabilization": {{
          "enabled": false
        }},
        "hierarchicalRepulsion": {{
          "nodeDistance": 150
        }}
      }},
      "layout": {{
        "hierarchical": {{
          "enabled": true,
          "direction": "DU",
          "sortMethod": "directed"
        }}
      }}
    }}
    """)
    
    # Add nodes with individual sizes and colors
    for i, org in enumerate(organizations):
        node_id = f"Org{org['index']}"
        
        # Calculate node size (scaled between 20 and 50)
        base_size = 20
        if max_size > min_size:
            # Logarithmic scaling for better visualization
            size_ratio = (org['size'] - min_size) / (max_size - min_size)
            node_size = base_size + size_ratio * 30
        else:
            node_size = base_size + 15
        
        # Calculate node color based on conflict level
        rgba = cmap(org['conflict_level'])
        color_hex = mcolors.to_hex(rgba)
        
        # Create hover text with information
        hover_text = (
            f"Organization {org['index']}\n"
            f"Size: {org['size']} species\n"
            f"Conflict level: {org['conflict_level']*100:.1f}%\n"
            f"Conflict species: {org['conflict_count']}\n"
            f"Species: {', '.join(sorted(org['species_list']))}"
        )
        
        # Add node with custom properties
        net.add_node(
            node_id,
            label=node_id,
            title=hover_text,
            color=color_hex,
            size=node_size,
            shape='dot',
            font={"size": 14}
        )
    
    # Add edges based on subset relationships (containment hierarchy)
    for i, org1 in enumerate(organizations):
        for j, org2 in enumerate(organizations):
            if i == j:
                continue
                
            # Check if org1 is subset of org2
            if org1['species_set'].issubset(org2['species_set']):
                # Check if it's a direct subset (no intermediate)
                is_direct = True
                for k, org3 in enumerate(organizations):
                    if (i != k and j != k and 
                        org1['species_set'].issubset(org3['species_set']) and 
                        org3['species_set'].issubset(org2['species_set'])):
                        is_direct = False
                        break
                
                if is_direct:
                    net.add_edge(f"Org{org1['index']}", f"Org{org2['index']}")
    
    # Save the visualization
    net.write_html(filename)
    
    print(f"Visualization created successfully: {filename}")
    
    # Create a legend
    create_color_legend(organizations, conflict_species, cmap, min_size, max_size)

def create_color_legend(organizations, conflict_species, cmap, min_size, max_size):
    """
    Create a text legend explaining the color and size coding.
    """
    print("\n" + "=" * 60)
    print("VISUALIZATION LEGEND")
    print("=" * 60)
    print("NODE COLORS:")
    print("  • Green → Red = No conflict → High conflict")
    print("  • Color intensity indicates conflict species presence")
    print("\nNODE SIZES:")
    print(f"  • Smallest: {min_size} species")
    print(f"  • Largest: {max_size} species")
    print("  • Size increases with number of species in organization")
    print("\nNODE LABELS:")
    print("  • OrgX: Organization number X")
    print("  • Hover over nodes to see details")
    
    print("\nCONFLICT SPECIES CONSIDERED:")
    for species in sorted(conflict_species):
        print(f"  - {species}")
    
    # Show statistics
    print("\nORGANIZATION STATISTICS:")
    print(f"  Total organizations: {len(organizations)}")
    
    if organizations:
        avg_size = sum(org['size'] for org in organizations) / len(organizations)
        avg_conflict = sum(org['conflict_level'] for org in organizations) / len(organizations) * 100
        
        print(f"  Average size: {avg_size:.1f} species")
        print(f"  Average conflict level: {avg_conflict:.1f}%")
        
        # Size distribution
        size_bins = defaultdict(int)
        for org in organizations:
            if org['size'] <= 5:
                size_bins["1-5"] += 1
            elif org['size'] <= 10:
                size_bins["6-10"] += 1
            elif org['size'] <= 20:
                size_bins["11-20"] += 1
            else:
                size_bins["21+"] += 1
        
        print(f"\n  Size distribution:")
        for bin_name, count in sorted(size_bins.items()):
            percentage = (count / len(organizations)) * 100
            print(f"    {bin_name} species: {count} orgs ({percentage:.1f}%)")
    
    # Conflict species frequency
    print(f"\nMOST COMMON CONFLICT SPECIES:")
    conflict_counts = defaultdict(int)
    for org in organizations:
        for species in org['species_list']:
            if species in conflict_species:
                conflict_counts[species] += 1
    
    if conflict_counts:
        for species, count in sorted(conflict_counts.items(), key=lambda x: x[1], reverse=True)[:10]:
            percentage = (count / len(organizations)) * 100
            print(f"  {species}: {count} organizations ({percentage:.1f}%)")
    else:
        print("  No conflict species found in any organization")

def save_text_report(organizations, conflict_species, filename="conflict_analysis_report.txt"):
    """
    Save a detailed text report of the analysis.
    """
    with open(filename, 'w') as f:
        f.write("CONFLICT ORGANIZATION ANALYSIS REPORT\n")
        f.write("=" * 60 + "\n\n")
        
        f.write(f"Total organizations: {len(organizations)}\n")
        f.write(f"Conflict species considered: {len(conflict_species)}\n")
        f.write("\n")
        
        # List conflict species
        f.write("Conflict species:\n")
        for species in sorted(conflict_species):
            f.write(f"  - {species}\n")
        f.write("\n")
        
        # Organizations by conflict level
        f.write("ORGANIZATIONS BY CONFLICT LEVEL (sorted high to low):\n")
        f.write("-" * 60 + "\n")
        
        for org in sorted(organizations, key=lambda x: x['conflict_level'], reverse=True):
            conflict_percent = org['conflict_level'] * 100
            type_str = "Elementary" if org['is_elementary'] else "Non-elementary"
            
            f.write(f"\nOrganization {org['index']} ({type_str}):\n")
            f.write(f"  Size: {org['size']} species\n")
            f.write(f"  Conflict level: {conflict_percent:.1f}% ({org['conflict_count']} conflict species)\n")
            f.write(f"  Species: {sorted(org['species_list'])}\n")
            
            # List conflict species in this organization
            conflict_in_org = [sp for sp in org['species_list'] if sp in conflict_species]
            if conflict_in_org:
                f.write(f"  Conflict species present: {sorted(conflict_in_org)}\n")
        
        # Statistics
        f.write("\n" + "=" * 60 + "\n")
        f.write("STATISTICS\n")
        f.write("-" * 60 + "\n")
        
        total_conflict_orgs = sum(1 for org in organizations if org['conflict_count'] > 0)
        f.write(f"Organizations with conflict species: {total_conflict_orgs} ({total_conflict_orgs/len(organizations)*100:.1f}%)\n")
        
        if organizations:
            avg_conflict = sum(org['conflict_level'] for org in organizations) / len(organizations) * 100
            f.write(f"Average conflict level: {avg_conflict:.1f}%\n")
        
        # Conflict species frequency
        f.write("\nCONFLICT SPECIES FREQUENCY:\n")
        conflict_counts = defaultdict(int)
        for org in organizations:
            for species in org['species_list']:
                if species in conflict_species:
                    conflict_counts[species] += 1
        
        for species, count in sorted(conflict_counts.items(), key=lambda x: x[1], reverse=True):
            percentage = (count / len(organizations)) * 100
            f.write(f"  {species}: {count} organizations ({percentage:.1f}%)\n")
    
    print(f"\nDetailed report saved to: {filename}")

# Main execution
if __name__ == "__main__":
    print("CONFLICT ORGANIZATION ANALYZER")
    print("=" * 60)
    
    # Choose your network file
    # Get the project root directory (where pyCOT is located)
    project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
    
    file_path = os.path.join(project_root, 'networks/Conflict_Theory/Resource_Community_Insurgency_Loops_model3.txt')
    # file_path = os.path.join(project_root, 'networks/Conflict_Theory/conflict_toy_model0.txt')
    # file_path = os.path.join(project_root, 'networks/Conflict_Theory/cause_driven_conflict_gov.txt')
    # file_path = os.path.join(project_root, 'networks/Conflict_Theory/Resource_Scarcity_Toy_Model2.txt')
    # file_path = 'Txt/Farm.txt'  # For testing
    
    # Check if file exists
    if not os.path.exists(file_path):
        print(f"ERROR: File not found: {file_path}")
        print("\nPlease provide a valid network file path.")
        print("\nAvailable options:")
        print("1. networks/Conflict_Theory/conflict_toy_model0.txt")
        print("2. networks/Conflict_Theory/cause_driven_conflict_gov.txt")
        print("3. networks/Conflict_Theory/Resource_Community_Insurgency_Loops_model1.txt")
        print("4. networks/Conflict_Theory/Resource_Scarcity_Toy_Model2.txt")
        print("5. Txt/Farm.txt (simple test network)")
        sys.exit(1)
    
    # Analyze organizations
    organizations, conflict_species, rn = analyze_organizations_with_conflict(file_path)
    
    if organizations is None:
        sys.exit(1)
    
    # Create visualization with variable node sizes
    base_filename = os.path.basename(file_path).replace('.txt', '')
    os.makedirs("reports", exist_ok=True)
    create_conflict_colored_visualization(
        organizations, 
        conflict_species,
        filename=f"visualizations/conflict_orgs_{base_filename}.html"
    )
    
    # Save detailed report
    save_text_report(organizations, conflict_species, 
                    filename=f"reports/conflict_report_{base_filename}.txt")
    
    print("\n" + "=" * 60)
    print("ANALYSIS COMPLETE")
    print("=" * 60)
    print(f"\nNetwork analyzed: {file_path}")
    print(f"Total organizations found: {len(organizations)}")
    print(f"\nOutputs generated:")
    print(f"  1. Interactive HTML visualization: conflict_orgs_{base_filename}.html")
    print(f"  2. Detailed text report: conflict_report_{base_filename}.txt")
    print(f"\nOpen the HTML file in your browser to see the visualization.")
    print("\nIn the visualization:")
    print("  • Node color: Green → Red (no conflict → high conflict)")
    print("  • Node size: Larger = more species in the organization")
    print("  • Hover over nodes to see details")
    print("  • Click and drag to rearrange if needed")