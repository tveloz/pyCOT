from pyvis.network import Network
import re

######################################################################################
# Paso 1) Función para graficar una red de reacciones en HTML:
def RN_view_html(RN, list_color_spcs, list_color_RnStr, filename="network.html"):
    net = Network(height='750px', width='100%', notebook=True, directed=True, cdn_resources='in_line')
    
    # Inicializar el conjunto de conexiones para verificar aristas opuestas
    connections = set()    
    species_color_map = {}  # Mapa para almacenar colores de especies
    
    ######################################
    # Procesar list_color_spcs para asignar colores a las especies
    for color, lst_spcs in list_color_spcs:
        for spcs in lst_spcs:
            species_color_map[spcs] = color
    
    # Añadir nodos de especies con colores definidos
    for species in species_color_map:
        net.add_node(species, shape='circle', label=species, color=species_color_map[species], size=15)

    reaction_color_map = {}  # Mapa para almacenar colores de reacciones
    
    ######################################
    # Procesar list_color_RnStr para asignar colores a las reacciones
    for color, lst_RnStr in list_color_RnStr:
        for RnStr in lst_RnStr:
            reaction_color_map[RnStr] = color  # Corregido para que use las reacciones en lugar de especies
    
    # Añadir nodos de reacciones con colores definidos
    for reaction in RN:
        net.add_node(reaction, shape='box', label=reaction, color=reaction_color_map.get(reaction, 'lightgray'), size=10)
    
    ######################################
     # Añadir aristas y etiquetas
    for reaction, (inputs, outputs) in RN.items():
        for species, coef in inputs:
            edge_id = (species, reaction)  # Identificador de la arista

            if (reaction, species) in connections:
                # Si ya existe una arista en dirección opuesta, curva ambas aristas
                if coef > 1:
                    net.add_edge(species, reaction, arrows='to', label=str(coef), color='red', font_color='red', smooth={'type': 'curvedCCW'})
                else:
                    net.add_edge(species, reaction, arrows='to', color='red', font_color='red', smooth={'type': 'curvedCCW'})
            else:
                if coef > 1:
                    net.add_edge(species, reaction, arrows='to', label=str(coef), color='red', font_color='red')
                else:
                    net.add_edge(species, reaction, arrows='to', color='red', font_color='red')

            # Guardar la conexión en el conjunto
            connections.add(edge_id)

        ######################################
        for species, coef in outputs:
            edge_id = (reaction, species)  # Identificador de la arista

            if (species, reaction) in connections:
                # Si ya existe una arista en dirección opuesta, curva ambas aristas
                if coef > 1:
                    net.add_edge(reaction, species, arrows='to', label=str(coef), color='green', font_color='green', smooth={'type': 'curvedCW'})
                else:
                    net.add_edge(reaction, species, arrows='to', color='green', font_color='green', smooth={'type': 'curvedCW'})
            else:
                if coef > 1:
                    net.add_edge(reaction, species, arrows='to', label=str(coef), color='green', font_color='green')
                else:
                    net.add_edge(reaction, species, arrows='to', color='green', font_color='green')

            # Guardar la conexión en el conjunto
            connections.add(edge_id)
            
    ######################################
    # Establecer opciones de visualización
    net.set_options('''{
        "physics": {"enabled": false},
        "edges": {"label": {"enabled": true, "font": {"size": 14}}}
    }''')
    
    # Guardar la visualización en un archivo HTML
    net.show(filename)

######################################################################################
# Paso 2) Función para interpretar las reacciones desde un archivo.
def interpret_reactions_from_archive(filename):
    RN = {}  # Inicializa un diccionario vacío para almacenar las reacciones.

    with open(filename, 'r') as file:
        reactions = file.readlines()  # Lee todas las líneas del archivo.

    for i, reaction in enumerate(reactions):
        reaction = reaction.strip().replace(";", "")  # Limpieza de la línea.

        if ":" in reaction:
            reaction = reaction.split(":")[-1].strip()  # Ignora las etiquetas.

        if "=>" not in reaction:
            continue  # Ignora líneas que no tengan la flecha de productos.

        reactants, products = reaction.split("=>")  # Divide en reactivos y productos.

        # Procesar reactantes.
        reactants_list = [item.strip() for item in reactants.split("+")]
        reactants_coeff = []
        for reactant in reactants_list:
            if reactant:
                match = re.match(r'(\d*)(\w+)', reactant.strip())
                if match:
                    coeff_str, var = match.groups()
                    coeff = int(coeff_str) if coeff_str else 1
                    reactants_coeff.append((var.strip(), coeff))

        # Procesar productos.
        products_list = [item.strip() for item in products.split("+")]
        products_coeff = []
        for product in products_list:
            if product:
                match = re.match(r'(\d*)(\w+)', product.strip())
                if match:
                    coeff_str, var = match.groups()
                    coeff = int(coeff_str) if coeff_str else 1
                    products_coeff.append((var.strip(), coeff))

        # Agregar la reacción al diccionario RN.
        RN[f'r{i}'] = (reactants_coeff, products_coeff)

    return RN


######################################################################################
# Paso 3) Especificar el archivo.txt que se va a leer.
filename = '/Users/yvanomarbalderamoreno/Desktop/Projects/pyCOT/pyCOT/autopoietic.txt'

# Interpretar el archivo y obtener la red de reacciones.
RN = interpret_reactions_from_archive(filename)

######################################################################################
# Paso 4) Obtener el conjunto de especies y reacciones.
species_set = set()
for inputs, outputs in RN.values():
    species_set.update(species for species, _ in inputs)
    species_set.update(species for species, _ in outputs)

# Vector de especies y reacciones.
species_vector  = list(species_set)
reaction_vector = list(RN.keys())

######################################################################################
# Paso 5) Definir colores de los nodos de especies y reacciones
color_species = 'cyan'       # Color del nodo de especies
list_color_species = [(color_species, species_vector)] #Color del nodo de especies

color_reaction = 'lightgray' #Color del nodo de reacciones
list_color_reaction = [(color_reaction, reaction_vector)]

# Imprimir la red de reacciones, especies y reacciones.
print("Network Reaction:", RN)
print("Species:", list_color_species)
print("Reactions:", list_color_reaction)
######################################################################################
# Paso 5) Llamar a la función para graficar la red de reacciones y generar archivo HTML.
RN_view_html(RN, list_color_species, list_color_reaction, filename="RN_autopoietic.html")
