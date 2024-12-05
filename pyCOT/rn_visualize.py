from pyvis.network import Network 


######################################################################################
# Function to plot a reaction network in HTML:
def rn_get_visualization(RN, lst_color_spcs=None, lst_color_reacs=None, 
                 global_species_color=None, global_reaction_color=None,
                 global_input_edge_color=None, global_output_edge_color=None, 
                 node_size=1000, curvature=None, physics_enabled=False, 
                 filename="reaction_network.html"):
    """
    Visualizes a reaction network in an interactive HTML file. 

    Parameters:
    - RN: ReactionNetwork object 
    - lst_color_spcs (list, optional): List of tuples that associates specific colors with species (e.g., [('cyan', ['A', 'B'])]).
    - lst_color_reacs (list, optional): List of tuples that associates specific colors with reactions (e.g., [('gray', ['R1', 'R2'])]).
    - global_species_color (str, optional): Default color for species, if no other color is specified.
    - global_reaction_color (str, optional): Default color for reactions, if no other color is specified.
    - global_input_edge_color (str, optional): Color of input arrows (consumption) to reactions.
    - global_output_edge_color (str, optional): Color of output arrows (production) from reactions.
    - node_size (int, optional): Base size for network nodes.
    - curvature (str, optional): None=False (Without curvature), True=curvedCCW (Counter Clockwise)
    - filename (str, optional): Name of the output HTML file.    

    Returns:
    - None. The function saves the visualization to an HTML file specified by `filename`.

    Example of usage: 
    # Define the file path where the reaction network model is located.
    file_path = 'Txt/Farm.txt'  

    # Load the reaction network from the specified file in `file_path`.
    testRN = load_pyCOT_from_file(file_path)  

    # Visualize the reaction network `testRN`, assigning the color yellow to species 's1' and saving the result in an HTML file.
    rn_visualize(testRN.RN, lst_color_spcs=[('yellow', ['s1'])], filename="Reaction_network.html")  
   
    """
 
    # Initialize the network visualization with specific size, settings, and directionality
    net = Network(height='750px', width='100%', notebook=True, directed=True, cdn_resources='in_line')

    ######################################     
    # Set physics options depending on ‘physics_enabled’.
    options = f"""
    var options = {{
        "physics": {{
            "enabled": {str(physics_enabled).lower()}
        }}
    }}
    """
    net.set_options(options)

    # Default colors: Define fallback colors for species, reactions, and edges if not provided
    default_species_color = global_species_color or 'cyan'
    default_reaction_color = global_reaction_color or 'lightgray'
    input_edge_color = global_input_edge_color or 'red'     # Color of input arrows to the reaction (Consume)
    output_edge_color = global_output_edge_color or 'green'  # Color of output arrows from the reaction (Produce)
    
    # Map specific colors to certain species or reactions if lists provided, otherwise use defaults
    species_colors = {species: color for color, species_list in (lst_color_spcs or []) for species in species_list}
    reaction_colors = {reaction: color for color, reaction_list in (lst_color_reacs or []) for reaction in reaction_list}
    
    #Construct RN dictionary from reactionnetwork object 
    RN_dict = {}  
    
    # Lista para almacenar los nombres de las reacciones
    reaction_names = []  

    # Iterar sobre las reacciones
    for i in range(len(RN.RnStr)):  # Usamos RnStr directamente, ya que RN aún no está completamente construido
        reactants_coeff = []
        products_coeff = []

        # Procesar reactivos
        for j in range(len(RN.SpStr)):  # Usamos SpStr, ya que SpStr contiene las especies
            coef_r = RN.RnMsupp[i][j]  # Coeficiente del reactivo en la reacción i para la especie j
            if coef_r > 0:
                reactants_coeff.append((RN.SpStr[j], coef_r))
        
        # Procesar productos
        for j in range(len(RN.SpStr)):  # Lo mismo para los productos
            coef_r = RN.RnMprod[i][j]  # Coeficiente del producto en la reacción i para la especie j
            if coef_r > 0:
                products_coeff.append((RN.SpStr[j], coef_r))

        # Crear el nombre de la reacción
        reaction_name = RN.RnStr[i] #f'R{i+1}'

        # Guardar la reacción en el diccionario
        RN_dict[reaction_name] = (reactants_coeff, products_coeff)

        # Agregar el nombre de la reacción a la lista de nombres de reacciones
        reaction_names.append(reaction_name)



    # Identify all species in the network to check their presence in lst_color_spcs
    species_vector = sorted(set([spcs for reactants, products in RN_dict.values() for spcs, _ in reactants + products]))
    species_set = set(species_vector)  # Convert to a set for easier comparison
    
    # Validate species in lst_color_spcs to ensure they belong to the network's species
    if lst_color_spcs:
        for color, species_list in lst_color_spcs:
            for species in species_list:
                if species not in species_set:
                    print(f"Warning: The species '{species}' specified in lst_color_spcs does not belong to the species of the network.")

    # Identify all reactions in the network to check their presence in lst_color_reacs
    reaction_vector = list(RN_dict.keys())
    reaction_set = set(reaction_vector)  # Convert to a set for easier comparison

    # Validate reactions in lst_color_reacs to ensure they belong to the network's reactions
    if lst_color_reacs:
        for color, reaction_list in lst_color_reacs:
            for reaction in reaction_list:
                if reaction not in reaction_set:
                    print(f"Warning: The reaction '{reaction}' specified in lst_color_reacs does not belong to the network reactions.")
    
    # Calculate node size based on character length of species name for better visualization 
    longest_species = max(RN.SpStr, key=len)
    node_size = len(longest_species) * 10   
    # Add species nodes with color defined in `lst_color_spcs` or global/default color
    for species in species_vector:
        color = species_colors.get(species, default_species_color)
        net.add_node(species, shape='circle', label=species, color=color, 
                     size=node_size, font={'size': 14, 'color': 'black'})

    # y_offset = 5  # Desplazamiento hacia abajo para la etiqueta

    # for species in species_vector:
    #     color = species_colors.get(species, default_species_color) 

    #     # Nodo principal
    #     net.add_node(species, shape='circle', color=color, 
    #                 size=node_size, font={'size': 14, 'color': 'black'})

    #     # Nodo para la etiqueta (debajo del nodo principal)
    #     label_node_id = f"{species}_label"
    #     net.add_node(label_node_id, label=species, shape='text', 
    #                 font={'size': 14, 'color': 'black'})
        
    #     # Agregar una arista invisible para posicionar la etiqueta debajo del nodo principal
    #     net.add_edge(species, label_node_id, color='rgba(0,0,0,0)', 
    #                 physics=False, smooth={'type': 'cubicBezier'},
    #                 length=y_offset)
 
################################################################################################ 

    # Add reaction nodes with color defined in `lst_color_reacs` or global/default color
    for reaction in reaction_vector:
        color = reaction_colors.get(reaction, default_reaction_color)
        net.add_node(reaction, shape='box', label=reaction, color=color, size=10, font={'size': 14, 'color': 'black'})

    ######################################
    # Add edges and labels
    connections = set() # Set of connections to track opposite edges and avoid redundant edges
    for reaction, (inputs, outputs) in RN_dict.items():
        for species, coef in inputs:
            edge_id = (species, reaction)
            default_curve='curvedCCW'
            smooth_type = {'type': default_curve} if curvature else None
            # Check for opposite edges and set edge direction and label accordingly
            if (reaction, species) in connections:
                if coef > 1:
                    net.add_edge(species, reaction, arrows='to', label=str(coef), color=input_edge_color, 
                                 font_color=input_edge_color, smooth=smooth_type)
                else:
                    net.add_edge(species, reaction, arrows='to', color=input_edge_color, 
                                 font_color=input_edge_color, smooth=smooth_type)
            else:
                if coef > 1:
                    net.add_edge(species, reaction, arrows='to', label=str(coef), color=input_edge_color, 
                                 font_color=input_edge_color)
                else:
                    net.add_edge(species, reaction, arrows='to', color=input_edge_color, 
                                 font_color=input_edge_color)

            # Add edge to connections set to avoid redundant edges
            connections.add(edge_id)

        for species, coef in outputs:
            edge_id = (reaction, species)
            default_curve='curvedCCW'
            smooth_type = {'type': default_curve} if curvature else None
            # Check for opposite edges and set edge direction and label accordingly
            if (species, reaction) in connections:
                if coef > 1:
                    net.add_edge(reaction, species, arrows='to', label=str(coef), color=output_edge_color, 
                                 font_color=output_edge_color, smooth=smooth_type)
                else:
                    net.add_edge(reaction, species, arrows='to', color=output_edge_color, 
                                 font_color=output_edge_color, smooth=smooth_type)
            else:
                if coef > 1:
                    net.add_edge(reaction, species, arrows='to', label=str(coef), color=output_edge_color, 
                                 font_color=output_edge_color)
                else:
                    net.add_edge(reaction, species, arrows='to', color=output_edge_color, 
                                 font_color=output_edge_color)

            # Add edge to connections set to avoid redundant edges
            connections.add(edge_id)

    # Save the visualization to an HTML file
    # net.show(filename)
    net.save_graph(filename) 
    return filename

 

def rn_visualize(RN, lst_color_spcs=None, lst_color_reacs=None, 
                 global_species_color=None, global_reaction_color=None,
                 global_input_edge_color=None, global_output_edge_color=None, 
                 node_size=1000, curvature=None, physics_enabled=False, 
                 filename="reaction_network.html"):
    
    # Llamar a la función `rn_get_visualization` para obtener la visualización
    visualization = rn_get_visualization(
        RN, lst_color_spcs, lst_color_reacs, 
        global_species_color, global_reaction_color,
        global_input_edge_color, global_output_edge_color, 
        node_size=node_size, curvature=curvature, 
        physics_enabled=physics_enabled, filename=filename
    )
    
    # Imprimir el valor de visualization para depuración
    print(f"The reaction network is: {visualization}")
    
    # Verificar si la visualización es válida
    if visualization is None or not isinstance(visualization, str):
        raise ValueError("The display was not generated correctly. A filename of type string was expected.")
    



