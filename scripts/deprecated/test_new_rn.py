# Import necessary libraries
import os
import sys
# Add the root directory to the PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import the necessary modules from pyCOT  
from pyCOT.io.functions import read_txt  
from pyCOT.simulations import build_reaction_dict

###################################################################################
# # # Load the reaction network from a file
# file_path = 'Txt/autopoietic.txt'
file_path = 'Txt_sub_network/sub_network_autopoietic.txt'

# Creates the ReactionNetwork object: rn
rn = read_txt(file_path) 

species = [specie.name for specie in rn.species()]
print("Species Set =",species)

reactions = [reaction.name() for reaction in rn.reactions()]
print("\nReactions Set =",reactions)

rn.add_species("p1")  # Agrega una nueva especie llamada p1
rn.add_species("p2")  # Agrega una nueva especie llamada p2
rn.add_reaction("R1", ["p1", "s1"], ["p2"], [1, 1], [2])  # Agrega una nueva reacción R4
rn.add_reaction("R6", ["p2", "s2"], ["p2"], [1, 1], [3])
import inspect
print(inspect.signature(rn.add_reaction))

# print("Object Species:\n",rn.species())
species = [specie.name for specie in rn.species()]
print("Species Set =",species)

# print("\nObject Reactions:\n",rn.reactions())
reactions = [reaction.name() for reaction in rn.reactions()]
print("\nReactions Set =",reactions)

# Print the stoichiometry matrix
S = rn.stoichiometry_matrix()
print("\nStoichiometry matrx:\n",S,"\n")

# # Print the reactions and their reactants and products
for i, reaction in enumerate(rn.reactions()):
      print(reaction)

# Función para construir un diccionario de reacciones del módulo simulations
rn_dict = build_reaction_dict(rn) 
print("\nrn_dict:",rn_dict)

# Print the reactions and their reactants and products
print("\nReactions and their reactants and products of rn_dict:")
for reaction_name, (reactants, products) in rn_dict.items():
    print(f"{reaction_name}: {reactants} -> {products}")

# Get the attributes of a specie
s = rn.get_species("p1")
print("name",s.name)
print("index",s.index)
print("quantily",s.quantity)

# Check if a species or reaction exists in the network
print(rn.has_species("Glucose"))  # True o False
print(rn.has_reaction("R1"))  # True o False

# Get the index of a species and a reaction
first_species = rn.get_species_by_index(6)  # Get the specie
print("Get the specie by index:",first_species.name)
first_reaction = rn.get_reaction_by_index(0)  # Get the reaction
print("Get the reaction by index:",first_reaction.name)

# Convierte distintas formas de entrada (string, objeto o lista) a una lista de objetos Species.
species_list = rn._parse_species_input(["p1", "p2"])
print([s.name for s in species_list])

# Get the edges of a reaction
rxn = rn.get_reaction("R1")
print(rxn.name())
print(rxn.edges)
edges = rn.get_reaction_edges_by_index(0) 
print(edges)

# Check if a reaction exists in the network
print(rn.has_reaction("R1"))  # True o False
print(rn.has_reaction("R7"))  # True o False

# Add a new reaction and species
rn.add_from_reaction_string("R7: A + B => C")
rn.add_species("l")  # si no está previamente
rn.add_from_reaction_string("R8: p1 + l => 2p2")

# Print the stoichiometry matrix
S = rn.stoichiometry_matrix()
print("\nStoichiometry matrx:\n",S,"\n")