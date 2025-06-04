###############################################################################
"""
Grafo Dirigido Bipartito de la Red de Reacciones:
	•Red de Reacciones Químicas: 
       R1: A + B → C
    •Nodos:            
       -Conjunto 1 (Especies):
        A, B, C
       -Conjunto 2 (Reacciones):
        R1
	•Aristas (Las aristas solo conectan nodos de diferentes conjuntos):
       A  ─▶ R1
       B  ─▶ R1
       R1 ─▶ C

Equema visual del flujo de clases:
          ReactionNetwork
               │
     ┌─────────┴─────────┐
     │                   │
 PyDiGraph          Métodos propios
     │               (add_reaction,
     │               get_species, etc.)
     ▼
[Nodos con info]
(Especies, Reacciones)
     │
     ▼
[Aristas dirigidas]
(Reactantes ➜ Reacción ➜ Productos)       
"""
# Import the necessary modules from rustworkx
from rustworkx import PyDiGraph

G = PyDiGraph()
A = G.add_node({"type": "species", "name": "A"})
B = G.add_node({"type": "species", "name": "B"})
C = G.add_node({"type": "species", "name": "C"})
R1 = G.add_node({"type": "reaction", "name": "R1"})

G.add_edge(A, R1, {"stoich": 1})
G.add_edge(B, R1, {"stoich": 1})
G.add_edge(R1, C, {"stoich": 1})

# Imprimir nodos
print("\n=== NODOS ===") 
print("G.nodes() =",G.nodes())
# Imprimir nodos de tipo 'species' y 'reaction' 
species_names = []
reaction_names = []
for node_data in G.nodes():
    if node_data['type'] == 'species':
        species_names.append(node_data['name'])
    if node_data['type'] == 'reaction':
        reaction_names.append(node_data['name'])
print(species_names)
print(reaction_names)        


# Imprimir aristas
print("\n=== ARISTAS + DATOS ===")
print("G.edges() =",G.edges())
print("G.edge_list() =",G.edge_list())
for u, v in G.edge_list():
    edge_data = G.get_edge_data(u, v)
    print(f"Edge from {u} to {v}: {edge_data}")