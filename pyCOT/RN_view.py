from pyvis.network import Network
import re
import os
import sys

# Añadir el directorio raíz al PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

######################################################################################
# Función para graficar una red de reacciones en HTML:    
def RN_view(RN, lst_color_spcs=None, lst_color_reacs=None, 
            global_species_color=None, global_reaction_color=None,
            global_input_edge_color=None, global_output_edge_color=None, filename="network.html"):
    
    net = Network(height='750px', width='100%', notebook=True, directed=True, cdn_resources='in_line')
    
    # Colores predeterminados
    default_species_color = global_species_color or 'cyan'
    default_reaction_color = global_reaction_color or 'lightgray'
    input_edge_color = global_input_edge_color or 'red'     # Color de flechas de entrada
    output_edge_color = global_output_edge_color or 'green'  # Color de flechas de salida
    
    # Conjunto de conexiones para verificar aristas opuestas
    connections = set()
    
    # Diccionarios para mapear colores específicos
    species_colors = {species: color for color, species_list in (lst_color_spcs or []) for species in species_list}
    reaction_colors = {reaction: color for color, reaction_list in (lst_color_reacs or []) for reaction in reaction_list}

    # Identificar todas las especies en la red para verificar su presencia en lst_color_spcs
    species_vector = sorted(set([spcs for reactants, products in RN.values() for spcs, _ in reactants + products]))
    species_set = set(species_vector)  # Convertimos a un conjunto para facilitar la comparación
    
    # Verificación de especies en lst_color_spcs que no pertenecen a species_vector
    if lst_color_spcs:
        for color, species_list in lst_color_spcs:
            for species in species_list:
                if species not in species_set:
                    print(f"Warning: The species '{species}' specified in lst_color_spcs does not belong to the species of the network.")
    
    # Identificar todas las reacciones en la red para verificar su presencia en lst_color_reacs
    reaction_vector = list(RN.keys())
    reaction_set = set(reaction_vector)  # Convertimos a un conjunto para facilitar la comparación

    # Verificación de reacciones en lst_color_reacs que no pertenecen a reaction_vector
    if lst_color_reacs:
        for color, reaction_list in lst_color_reacs:
            for reaction in reaction_list:
                if reaction not in reaction_set:
                    print(f"Warning: The reaction '{reaction}' specified in lst_color_reacs does not belong to the network reactions.")
    
    # Añadir nodos de especies con color definido en `lst_color_spcs` o color global/predeterminado
    for species in species_vector:
        color = species_colors.get(species, default_species_color)
        net.add_node(species, shape='circle', label=species, color=color, size=15)

    # Añadir nodos de reacciones con color definido en `lst_color_reacs` o color global/predeterminado
    for reaction in reaction_vector:
        color = reaction_colors.get(reaction, default_reaction_color)
        net.add_node(reaction, shape='box', label=reaction, color=color, size=10)
    
    ######################################
    # Añadir aristas y etiquetas
    for reaction, (inputs, outputs) in RN.items():
        for species, coef in inputs:
            edge_id = (species, reaction)

            if (reaction, species) in connections:
                if coef > 1:
                    net.add_edge(species, reaction, arrows='to', label=str(coef), color=input_edge_color, font_color=input_edge_color, smooth={'type': 'curvedCCW'})
                else:
                    net.add_edge(species, reaction, arrows='to', color=input_edge_color, font_color=input_edge_color, smooth={'type': 'curvedCCW'})
            else:
                if coef > 1:
                    net.add_edge(species, reaction, arrows='to', label=str(coef), color=input_edge_color, font_color=input_edge_color)
                else:
                    net.add_edge(species, reaction, arrows='to', color=input_edge_color, font_color=input_edge_color)

            connections.add(edge_id)

        for species, coef in outputs:
            edge_id = (reaction, species)

            if (species, reaction) in connections:
                if coef > 1:
                    net.add_edge(reaction, species, arrows='to', label=str(coef), color=output_edge_color, font_color=output_edge_color, smooth={'type': 'curvedCW'})
                else:
                    net.add_edge(reaction, species, arrows='to', color=output_edge_color, font_color=output_edge_color, smooth={'type': 'curvedCW'})
            else:
                if coef > 1:
                    net.add_edge(reaction, species, arrows='to', label=str(coef), color=output_edge_color, font_color=output_edge_color)
                else:
                    net.add_edge(reaction, species, arrows='to', color=output_edge_color, font_color=output_edge_color)

            connections.add(edge_id)
            
    ######################################
    # # Establecer opciones de visualización
    # net.set_options('''{
    #     "physics": {"enabled": false},
    #     "edges": {"label": {"enabled": true, "font": {"size": 14}}}
    # }''')
    
    # # Guardar la visualización en un archivo HTML
    # net.show(filename)
    # Generate HTML content
    html_content = net.generate_html()

    # Save the HTML content with UTF-8 encoding
    with open(filename, 'w', encoding='utf-8') as file:
        file.write(html_content)



######################################################################################
# Función para interpretar las reacciones desde un archivo.
def interpret_reactions_from_archive(filename):
    RN = {}  # Diccionario para almacenar las reacciones.
    species_set = set()  # Conjunto para almacenar especies únicas.
    reaction_list = []  # Lista para almacenar nombres de las reacciones.

    # Verificación de existencia del archivo
    if not os.path.isfile(filename):
        print(f"The path to the file named {filename} is not found.")
        sys.exit(1)

    # Lectura y procesamiento del archivo
    with open(filename) as file:
        reactions = file.readlines()
        # print(f"El código se ha ejecutado con éxito.")

    for i, reaction in enumerate(reactions):
        reaction = reaction.strip().replace(";", "")
        
        if ":" in reaction:
            reaction = reaction.split(":")[-1].strip()
            
        if "=>" not in reaction:
            continue

        reactants, products = reaction.split("=>")
        
        reactants_list = [item.strip() for item in reactants.split("+")]
        reactants_coeff = []
        for reactant in reactants_list:
            if reactant:
                match = re.match(r'(\d*)(\w+)', reactant.strip())
                if match:
                    coeff_str, var = match.groups()
                    coeff = int(coeff_str) if coeff_str else 1
                    reactants_coeff.append((var.strip(), coeff))
                    species_set.add(var.strip())  # Agrega la especie al conjunto de especies

        products_list = [item.strip() for item in products.split("+")]
        products_coeff = []
        for product in products_list:
            if product:
                match = re.match(r'(\d*)(\w+)', product.strip())
                if match:
                    coeff_str, var = match.groups()
                    coeff = int(coeff_str) if coeff_str else 1
                    products_coeff.append((var.strip(), coeff))
                    species_set.add(var.strip())  # Agrega la especie al conjunto de especies

        reaction_name = f'R{i+1}'
        RN[reaction_name] = (reactants_coeff, products_coeff)
        reaction_list.append(reaction_name)  # Agrega la reacción a la lista de reacciones

    # print("Red de Reacciones:", RN)
    print("Species:", sorted(species_set))
    print("Reactions:", reaction_list)
    # Retorna el diccionario de la red, las especies y las reacciones
    return RN, sorted(species_set), reaction_list

