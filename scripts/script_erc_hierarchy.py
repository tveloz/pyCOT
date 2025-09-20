from pyCOT.rn_rustworkx import ReactionNetwork
from pyCOT.ERC_Hierarchy import ERC, species_list_to_names
from pyCOT.io.functions import read_txt
#from pyCOT.synergy_calculator import SynergyCalculator
import matplotlib.pyplot as plt
import time
import networkx as nx
import numpy as np


# Load network and compute ERCs
file_path= 'networks/testing/ERCs_test2.txt'
#file_path= 'networks/testing/Farm.txt'
#file_path= 'networks/Navarino/RN_IN_05.txt'
#file_path= 'networks/biomodels_interesting/bigg_iAF692.txt'  # Adjust path as needed
print("Loading network and computing ERCs...")
RN = read_txt(file_path)
ercs = ERC.ERCs(RN)
hierarchy= ERC.build_hierarchy_graph(ercs, RN)

print(f"\nNetwork statistics:")
print(f"- Reactions: {len(RN.reactions())}")
print(f"- Species: {len(RN.species())}")
print(f"- ERCs found: {len(ercs)}")

# Plot original hierarchy

print("Displaying list of ERCs including generators and closures:")
for erc in ercs:
    print(f"\nERC {erc.label}:")
    print(f"  Minimal Generators: {species_list_to_names(erc.min_generators)}")
    print(f"  Closure: {species_list_to_names(erc.get_closure(RN))}")
    print(f"  All Generators: {species_list_to_names(erc.all_generators)}")
#plt.savefig(f"hierarchy.png")

print("ERC Hierarchy")
ERC.plot_hierarchy(ercs, RN, hierarchy, title="Original ERC Hierarchy")

