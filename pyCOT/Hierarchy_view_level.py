import networkx as nx
import matplotlib.pyplot as plt
from mplcursors import cursor
import numpy as np

#####################################################################################
# Función para crear y visualizar el grafo jerárquico
def Hierarchy_view_level(input_data):
    #Ordenar los conjuntos por la cantidad de elementos de la segunda coordenada 
    # input_data.sort(key=lambda x: len(x[1]))
    input_data.sort(key=lambda x: x[2])# Ordenar por el tercer elemento

    # Crear un grafo dirigido
    G = nx.DiGraph()

    # Añadir nodos y sus elementos, asignando su nivel
    for name, elements, level in input_data:
        G.add_node(name, elements=elements, level=level)

    # Definir automáticamente las relaciones de inclusión entre niveles
    for i in range(len(input_data)):
        for j in range(len(input_data)):
            if i != j:  # Asegurarse de que no se compare el nodo consigo mismo
                set_i = set(input_data[i][1])  # Convertir la segunda coordenada de la entrada 'i' en un conjunto
                set_j = set(input_data[j][1])  # Convertir la segunda coordenada de la entrada 'j' en un conjunto

                # Calcular la diferencia de niveles
                k = input_data[j][2] - input_data[i][2]
                
                # Agregar flechas desde el nivel 0 al nivel 1
                for node in input_data:
                    if node[2] == k and k>3:  # Nivel 1
                        G.add_edge(input_data[0][0], node[0])  # Añadir arista desde "S0" a cada nodo en el nivel j
                                
                # Comprobar si 'set_i' está contenido en 'set_j' y si k >= 1
                if set_i.issubset(set_j) and k >= 1:
                    # Comprobar si ya se ha añadido una flecha indirecta
                    indirect_edge_exists = False
                    for neighbor in G.successors(input_data[i][0]):
                        if set(input_data[G.nodes[neighbor]['level']][1]).issubset(set_j):
                            indirect_edge_exists = True
                            break 
                    if not indirect_edge_exists:
                        G.add_edge(input_data[i][0], input_data[j][0])  # Añadir una arista dirigida entre 'i' y 'j'
                        break 
 
   # Agregar flechas desde el nivel 0 al nivel 1
    for node in input_data:
        if node[2] == 1:  # Nivel 1
            G.add_edge(input_data[0][0], node[0])  # Añadir arista desde "S0" a cada nodo en el nivel 1

    #####################################################################################
    # Función para generar posiciones automáticas del grafo
    def generate_positions(graph):
        pos = {}
        levels = {}

        for node in graph.nodes:
            level = G.nodes[node]['level']
            if level not in levels:
                levels[level] = []
            levels[level].append(node)

        for level, nodes in levels.items():
            y_offset = level

            if len(nodes) == 1:
                x_offset = [0]
            else:
                x_offset = np.linspace(-1, 1, len(nodes))

            for i, node in enumerate(nodes):
                pos[node] = (x_offset[i], y_offset)

        return pos

    # Generar las posiciones
    pos_hierarchy = generate_positions(G)

    # Dibujar el grafo
    plt.figure(figsize=(8, 6))

    # Dibujar nodos
    nx.draw_networkx_nodes(G, pos=pos_hierarchy, node_shape='s', node_size=350, node_color='lightblue')

    # Dibujar las aristas y etiquetas con flechas rectas
    nx.draw_networkx_edges(G, pos=pos_hierarchy, arrows=True, arrowstyle='-|>', arrowsize=20, edge_color='gray', width=2,
                         connectionstyle='arc3,rad=0.') #rad=0.2:cambia la curvatura de la flecha

    nx.draw_networkx_labels(G, pos=pos_hierarchy, font_size=10, font_weight='bold')

    #####################################################################################
    # Crear un cursor que muestre solo los elementos de los conjuntos al pasar el mouse
    cursor_labels = [
        ', '.join(elements) if elements else '∅'
        for name, elements, _ in input_data
    ]

    cursor_obj = cursor(plt.gca(), hover=True)

    @cursor_obj.connect("add")
    def on_add(sel):
        index = sel.index
        if index < len(cursor_labels):
            sel.annotation.set_text(cursor_labels[index])
            sel.annotation.set_visible(True)

    @cursor_obj.connect("remove")
    def on_remove(sel):
        sel.annotation.set_visible(False)

    plt.title("Jerarquías por niveles")
    plt.gca().set_frame_on(False)

    plt.show()

#####################################################################################
# Llamar a la función con los datos de entrada
# input_data = [
#     ["S0", [], 0],
#     ["S1", ["s1", "s2"], 1],
#     ["S2", ["s1", "s3"], 1],
#     ["S3", ["s1", "s2", "s3"], 3]
# ]

# input_data = [
#     ["S0", [], 0],
#     ["S1", ["s1"], 1],
#     ["S2", ["s2"], 1],
#     ["S3", ["s3"], 1],
#     ["S4", ["s1", "s3"], 3],
#     ["S5", ["s1", "s2", "s3"], 4]
# ]

input_data = [
    ["S0", [], 0],  # Nodo con su nivel
    ["S1", ["s1"], 1],
    ["S2", ["s2"], 1],
    ["S3", ["s3"], 1],
    ["S4", ["s4"], 1],
    ["S5", ["s1", "s3"], 2],
    ["S6", ["s2", "s4"], 2],
    ["S7", ["s1", "s2", "s3", "s4"], 3]
]

# input_data = [
#     ["S0", [], 0],  
#     ["S1", ["s1", "s2"], 1],
#     ["S2", ["s1", "s3"], 1],
#     ["S3", ["s2", "s3"], 1],
#     ["S4", ["s1", "s2", "s4"], 2],
#     ["S5", ["s1", "s3", "s4"], 2],
#     ["S6", ["s2", "s3", "s4"], 2],
#     ["S7", ["s1", "s2", "s3", "s4"], 3]
# ]
 
 
# input_data = [
#     ["S0", [],0],
#     ["S1", ["s4"],1],     
#     ["S2", ["s1", "s2"],2],
#     ["S3", ["s1", "s2", "s3"],3],   
#     ["S4", ["s1", "s2", "s3", "s4"],4]
# ]

# input_data = [
#     ["A", [],0],
#     ["B", ["s1"],1],
#     ["K", ["s3"],1],    
#     ["C", ["s1", "s2"],2],  
#     ["D", ["s1", "s3"],2],
#     ["E", ["s2", "s3"],2],
#     ["F", ["s1", "s2", "s4"],3],
#     ["G", ["s1", "s3", "s4"],3],
#     ["H", ["s2", "s3", "s4"],3],
#     ["I", ["s1", "s2", "s3", "s4"],4],      
#     ["J", ["s1", "s2", "s3", "s4", "s5"],5]
# ]
 

# input_data = [
#     ["∅", [],0],
#     ["X", ["s1", "s2"],2],
#     ["Y", ["s1", "s2", "s3"],3],
#     ["Z", ["s4"],1],    
#     ["M", ["s1", "s2", "s3", "s4"],4]
# ]

# input_data = [
#     ["∅", [],0],
#     ["X", ["s1", "s2"],1],
#     ["Y", ["s1", "s2", "s3"],2],
#     ["Z", ["s4"],2],    
#     ["M", ["s1", "s2", "s3", "s4"],3]
# ]

Hierarchy_view_level(input_data)
