from pyCOT.rn_rustworkx import ReactionNetwork
from pyCOT.ERC_Hierarchy import ERC, species_list_to_names
from pyCOT.io.functions import read_txt
from pyCOT.synergy_calculator import SynergyCalculator
import matplotlib.pyplot as plt
import time
import networkx as nx
import numpy as np


# Load network and compute ERCs
file_path= 'networks/testing/ERCs_test2.txt'  # Adjust path as needed
print("Loading network and computing ERCs...")
RN = read_txt(file_path)
ercs = ERC.ERCs(RN)
hierarchy = ERC.build_hierarchy_graph(ercs, RN)

print(f"\nNetwork statistics:")
print(f"- Reactions: {len(RN.reactions())}")
print(f"- Species: {len(RN.species())}")
print(f"- ERCs found: {len(ercs)}")

# Plot original hierarchy
print("ERC Hierarchy")
plt.figure(figsize=(10, 8))
ERC.plot_hierarchy_synergies(ercs, RN, hierarchy, title="Original ERC Hierarchy")
print("Displaying list of ERCs including generators and closures:")
for erc in ercs:
    print(f"\nERC {erc.label}:")
    print(f"  Minimal Generators: {species_list_to_names(erc.min_generators)}")
    print(f"  Closure: {species_list_to_names(erc.get_closure(RN))}")
    print(f"  All Generators: {species_list_to_names(erc.all_generators)}")
#plt.savefig(f"hierarchy.png")


# Initialize and run synergy calculator
print("Computing Synergies")
calculator = SynergyCalculator(RN, ercs, hierarchy)

start_time = time.time()
minimal_synergies = calculator.find_minimal_synergies()
computation_time = time.time() - start_time

print(f"\nComputation time: {computation_time:.2f} seconds")
print(f"Found {len(minimal_synergies)} minimal synergies")

# Display synergies
print("\nSynergy List with Contributions")
print("===============================")
synergy_count = 0
for base1, base2, target in minimal_synergies:
    synergy_count += 1
    print(f"\n{synergy_count}. {base1.label}+{base2.label}→{target.label}")
    
    # Get closures
    base1_species = set(species_list_to_names(base1.get_closure(RN)))
    base2_species = set(species_list_to_names(base2.get_closure(RN)))
    target_species = set(species_list_to_names(target.get_closure(RN)))
    
    # Show closures
    print(f"   Base1 closure: {base1_species}")
    print(f"   Base2 closure: {base2_species}")
    print(f"   Target closure: {target_species}")
    
    # Analyze contributions for each target generator
    print("   Generator Analysis:")
    for gen in target.min_generators:
        gen_species = set(species_list_to_names(gen))
        if gen_species.issubset(base1_species | base2_species):
            # Calculate unique contributions
            base1_contrib = gen_species & base1_species - base2_species
            base2_contrib = gen_species & base2_species - base1_species
            shared_contrib = gen_species & base1_species & base2_species
            
            if base1_contrib or base2_contrib:
                print(f"   - Generator {gen_species}:")
                if base1_contrib:
                    print(f"     {base1.label} contributes: {base1_contrib}")
                if base2_contrib:
                    print(f"     {base2.label} contributes: {base2_contrib}")
                if shared_contrib:
                    print(f"     Shared species: {shared_contrib}")
    print()

# Convert to visualization format and plot
print("Synergy Hierarchy Visualization")
synergy_details = calculator.get_synergy_details(minimal_synergies)

plt.figure(figsize=(12, 8))
ERC.plot_hierarchy_synergies(ercs, RN, hierarchy, synergy_details, title="ERC Hierarchy with Synergies")
#plt.title(f"ERC Hierarchy with Synergies\n({len(minimal_synergies)} synergies found)")
#plt.savefig(f"synergy_hierarchy_{name.replace(' ', '_')}.png")
plt.show()
