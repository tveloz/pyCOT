#!/usr/bin/env python
"""
Script to load a reaction network from a file using rn_rustworkx
"""
import sys

import argparse
import os
from pyCOT.rn_rustworkx import ReactionNetwork
from pyCOT.io.functions import read_txt

file_path = "Txt/autopoietic.txt"
rn = read_txt(file_path)


# Optional: Display the stoichiometry matrix
try:
    stoich_matrix = rn.stoichiometry_matrix()
    print("\nStoichiometry Matrix:")
    print(stoich_matrix)
except Exception as e:
    print(f"\nCould not display stoichiometry matrix: {str(e)}")


inflow = rn.inflow_species()
print(f"Inflow species: {[s.name for s in inflow]}")

# For each species, find what can be produced from it
for species in rn.species():
    products = rn.get_prod_from_species(species.name)
    if products:
        print(f"From {species.name}, can produce: {[p.name for p in products]}")
