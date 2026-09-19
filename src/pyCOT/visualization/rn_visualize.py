from collections import Counter
from collections import defaultdict
from graphviz import Digraph
from IPython.display import Image
# from PIL import Image
from pyvis.network import Network
import networkx as nx
import rustworkx as rx
import matplotlib.pyplot as plt
import mplcursors 
import os
import re
import base64
import webbrowser
import sys
sys.stdout.reconfigure(encoding='utf-8')
# import codecs  # Agregado para manejar codificación UTF-8
    
import tempfile

from pyCOT.simulations.core import build_reaction_dict

##################################################################
# Función de apoyo para renderizar "LaTeX" visual con Unicode
##################################################################
def _tokenize_body(remainder, force_literal=False):
    """Tokeniza la parte sin carga de un nombre de especie.

    - force_literal=True (viene del prefijo '=', ver _tokenize_chemical_name):
      fuerza texto literal sin importar dígitos/guiones. Para nombres-etiqueta
      como P680/P700 que NO son fórmulas químicas y no deben subindexarse.
    - Si contiene guion, espacio o coma: se trata de un nombre trivial/IUPAC
      con locantes (posiciones), NO de una fórmula química. Se renderiza
      literal, sin convertir dígitos en subíndice (ej. "Glucose-6-Phosphate",
      "1,3-Bisphosphoglycerate", "S-Adenosyl Methionine").
    - Si NO contiene esos separadores y tiene dígitos embebidos: se asume
      fórmula química compacta y los dígitos se tokenizan como subíndice
      (ej. "H2O" -> H,2(sub),O ; "FADH2" -> FADH,2(sub)).
    - En cualquier otro caso (abreviaturas sin dígito: "CoA", "Pi", "PPi",
      "ATP", "FAD", ...): segmento único literal, sin heurística adicional.
    """
    if not remainder:
        return []

    if force_literal:
        return [(remainder, 'normal')]

    if re.search(r'[-,\s]', remainder):
        return [(remainder, 'normal')]

    if re.search(r'\d', remainder):
        segments = []
        for chunk in re.findall(r'[A-Za-z]+|\d+', remainder):
            segments.append((chunk, 'sub' if chunk.isdigit() else 'normal'))
        return segments

    return [(remainder, 'normal')]


def _tokenize_chemical_name(name):
    """
    Descompone un nombre de especie en segmentos ordenados (texto, tipo),
    tipo en {'normal', 'sub', 'super'}, con la siguiente prioridad:

    0. Prefijo '=' explícito: fuerza el cuerpo a texto literal (sin
       auto-subíndice de dígitos), preservando la detección de carga
       posterior. Uso: nombres-etiqueta con dígitos que no son fórmulas
       (ej. "=P680", "=P680+", "=P680^*").
    1. Guion bajo explícito '_': todo lo posterior es subíndice (ej. P_act).
    2. Circunflejo explícito '^': escape hatch para forzar la lectura
       correcta cuando la convención por defecto (punto 3) no aplica, o
       para cargas de magnitud >1 sin conteo atómico (ej. "Ca^2+" -> Ca²⁺).
    3. Carga iónica final detectada automáticamente: uno o más '+'/'-' al
       final de la cadena. CONVENCIÓN POR DEFECTO: cualquier dígito
       inmediatamente anterior al signo se interpreta como parte del
       cuerpo (conteo atómico -> subíndice), y el signo como carga unitaria
       en superíndice. Resuelve correctamente "NH4+" -> NH₄⁺,
       "H3O+" -> H₃O⁺, "HCO3-" -> HCO₃⁻ sin necesidad de escape.
       LIMITACIÓN: iones con carga de magnitud >1 y SIN conteo atómico
       adyacente (ej. Ca2+, Fe3+) requieren el escape '^' ("Ca^2+"),
       porque el mismo patrón sintáctico (dígito+signo) es ambiguo entre
       "conteo+carga unitaria" y "carga de magnitud N"; no es resoluble
       por regex sin conocimiento semántico del ion.
    4. Nombre trivial/IUPAC con guion, espacio o coma: se renderiza
       literal (ver _tokenize_body).
    5. Fórmula compacta sin separadores: dígitos embebidos -> subíndice.
    """
    force_literal = name.startswith('=')
    if force_literal:
        name = name[1:]

    if '_' in name:
        base, sub = name.split('_', 1)
        segments = [(base, 'normal')]
        if sub:
            segments.append((sub, 'sub'))
        return segments

    if '^' in name:
        body, charge = name.split('^', 1)
        segments = _tokenize_body(body, force_literal)
        if charge:
            segments.append((charge, 'super'))
        return segments

    charge_match = re.match(r'^(.+?)([+-]+)$', name)
    if charge_match:
        remainder, charge_sign = charge_match.groups()
        segments = _tokenize_body(remainder, force_literal)
        segments.append((charge_sign, 'super'))
        return segments

    return _tokenize_body(name, force_literal)


def _segments_to_svg_text(segments, font_size, sub_font_size, x, y):
    """
    Convierte una lista de segmentos (texto, tipo) en un elemento <text> SVG
    con <tspan> anidados para subíndices y superíndices, ajustando la posición vertical según el tipo. 
    Los segmentos se renderizan en orden, y 
    los superíndices se apilan sobre los subíndices inmediatamente anteriores si están adyacentes.
    """
    CHAR_WIDTH_FACTOR = 0.75
    SUB_DY_FACTOR = 0.25
    SUPER_DY_FACTOR = 0.35
    STACK_OFFSET_RATIO = 0.7   # <-- NUEVO: 1.0 = justo encima, 0.0 = secuencial a la derecha.
                               #     Baja este valor para correr el superíndice más a la derecha.

    def shift_for(kind):
        if kind == 'sub':
            return font_size * SUB_DY_FACTOR
        if kind == 'super':
            return -font_size * SUPER_DY_FACTOR
        return 0.0

    parts = []
    current_shift = 0.0
    prev_kind, prev_text = None, ""
    for text, kind in segments:
        size = font_size if kind == 'normal' else sub_font_size
        target_shift = shift_for(kind)
        dy = target_shift - current_shift

        dx_attr = ""
        if kind == 'super' and prev_kind == 'sub':
            back = len(prev_text) * sub_font_size * CHAR_WIDTH_FACTOR * STACK_OFFSET_RATIO
            dx_attr = f' dx="-{back:.2f}"'

        dy_attr = f' dy="{dy:.2f}"' if dy != 0 else ""
        parts.append(f'<tspan{dx_attr}{dy_attr} font-size="{size}">{text}</tspan>')

        current_shift = target_shift
        prev_kind, prev_text = kind, text

    inner = ''.join(parts)
    return f'<text x="{x}" y="{y}" font-family="Arial, sans-serif" text-anchor="middle" fill="black">{inner}</text>'


def _segments_width(segments, normal_w, sub_w):
    """Ancho efectivo ponderado. Cuando un 'super' se apila sobre un 'sub'
    inmediatamente anterior (ver _segments_to_svg_text), no se suman por
    separado: ocupan la misma columna, así que se toma el máximo de los
    dos anchos en vez de la suma."""
    total = 0.0
    i, n = 0, len(segments)
    while i < n:
        text, kind = segments[i]
        w = normal_w if kind == 'normal' else sub_w
        this_width = len(text) * w
        if kind == 'sub' and i + 1 < n and segments[i + 1][1] == 'super':
            next_text = segments[i + 1][0]
            this_width = max(this_width, len(next_text) * sub_w)
            total += this_width
            i += 2
            continue
        total += this_width
        i += 1
    return total


def generate_svg_data_uri(name, bg_color, shape_type, tipo='specie', use_latex_style=False):
    """
    Genera una imagen SVG. Soporta subíndices y superíndices (carga iónica),
    y ajusta la posición del texto dependiendo de si la forma es 'dot'
    (texto abajo) o encapsulada (texto adentro) en el nodo.
    """
    if not use_latex_style:
        return None

    segments = _tokenize_chemical_name(name)

    # --- AJUSTE DINÁMICO DE TAMAÑO Y POSICIÓN ---
    if tipo == 'reaction':
        # LAS REACCIONES MANTIENEN SU FORMATO COMPACTO CON TEXTO ADENTRO
        font_size = 6          
        sub_font_size = 6
        width = max(20, 5 + _segments_width(segments, normal_w=4, sub_w=2))
        height = 15             
        cx, cy = width / 2, height / 3
        
        shape_svg = f'<rect x="2" y="2" width="{width-4}" height="{height-4}" rx="4" ry="4" fill="{bg_color}" stroke="{bg_color}" stroke-width="2"/>'
        text_svg = _segments_to_svg_text(segments, font_size, sub_font_size, cx, cy + 4)

    else:
        # LAS ESPECIES EVALÚAN SI DEBEN PONER EL TEXTO AFUERA O ADENTRO
        font_size = 16          
        sub_font_size = 12      
        width = max(45, 20 + _segments_width(segments, normal_w=12, sub_w=9))
        cx = width / 2

        if shape_type == 'dot':
            # CASO 'DOT': Círculo en la parte superior, texto en la parte inferior
            radio = 16
            gap = 6  # Espacio entre el círculo y el texto
            y_text = (radio * 2) + gap + 10 # Posición Y del texto
            height = y_text + 10             # Alto total del lienzo aumentado
            
            # Dibujamos el círculo en la parte superior (cy = radio + 2 para el borde)
            shape_svg = f'<circle cx="{cx}" cy="{radio + 2}" r="{radio}" fill="{bg_color}" stroke="{bg_color}" stroke-width="2"/>'
            text_svg = _segments_to_svg_text(segments, font_size, sub_font_size, cx, y_text)
        
        else:
            # OTROS CASOS ('circle', 'box'): Texto centrado adentro de la forma
            height = 32             
            cy = height / 2
            
            if shape_type in ['circle', 'ellipse']:
                shape_svg = f'<circle cx="{cx}" cy="{cy}" r="{height/2 - 2}" fill="{bg_color}" stroke="{bg_color}" stroke-width="2"/>'
            else:
                shape_svg = f'<rect x="2" y="2" width="{width-4}" height="{height-4}" rx="4" ry="4" fill="{bg_color}" stroke="{bg_color}" stroke-width="2"/>'
                
            text_svg = _segments_to_svg_text(segments, font_size, sub_font_size, cx, cy + 5)

    # Empaquetar y exportar
    svg = f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">{shape_svg}{text_svg}</svg>'
    b64_encoded = base64.b64encode(svg.encode('utf-8')).decode('utf-8')
    return f"data:image/svg+xml;base64,{b64_encoded}"


# ##################################################################
# # Plot a reaction network in HTML
# ##################################################################
def rn_get_visualization(rn, lst_color_spcs=None, lst_color_reacs=None, 
                         global_species_color=None, global_reaction_color=None,
                         global_input_edge_color=None, global_output_edge_color=None, 
                         node_size=20, shape_species_node='dot', shape_reactions_node='box', 
                         curvature=None, physics_enabled=False, 
                         use_latex_style=False, # NUEVO PARÁMETRO AQUÍ
                         species_display_names=None, # NUEVO: nombre -> etiqueta LaTeX-style para el SVG (no afecta el id real del nodo)
                         filename="reaction_network.html"):
    """
    Visualizes a reaction network as an interactive HTML file.
    """
    net = Network(height='100vh', width='100%', notebook=True, directed=True, cdn_resources='in_line') # Antes height='750px'
 
    options = f"""
    var options = {{
        "physics": {{
            "enabled": {str(physics_enabled).lower()}
        }}
    }}
    """
    net.set_options(options)

    default_species_color = global_species_color or 'cyan'
    default_reaction_color = global_reaction_color or 'lightgray'
    input_edge_color = global_input_edge_color or 'red'
    output_edge_color = global_output_edge_color or 'green' 
    
    species_colors = {species: color for color, species_list in (lst_color_spcs or []) for species in species_list}
    reaction_colors = {reaction: color for color, reaction_list in (lst_color_reacs or []) for reaction in reaction_list}
    
    RN_dict = build_reaction_dict(rn)

    species_vector = sorted(set([spcs for reactants, products in RN_dict.values() for spcs, _ in reactants + products]))
    species_set = set(species_vector) 
    
    if lst_color_spcs:
        for color, species_list in lst_color_spcs:
            for species in species_list:
                if species not in species_set:
                    print(f"Warning: The species '{species}' specified in lst_color_spcs does not belong to the species of the network.")

    reaction_vector = list(RN_dict.keys())
    reaction_set = set(reaction_vector) 

    if lst_color_reacs:
        for color, reaction_list in lst_color_reacs:
            for reaction in reaction_list:
                if reaction not in reaction_set:
                    print(f"Warning: The reaction '{reaction}' specified in lst_color_reacs does not belong to the network reactions.")

    ######################################
    # AGREGAR NODOS DE ESPECIES
    ######################################
    display_names = species_display_names or {}
    for species in species_vector:
        color = species_colors.get(species, default_species_color)
        # render_size = node_size + 10 if use_latex_style else node_size
        # Quitamos o reducimos el +10 para que PyVis no agrande la imagen
        render_size = node_size if use_latex_style else node_size
        
        if use_latex_style:
            # Pasamos tipo='specie' a la función generadora.
            # display_names solo cambia la ETIQUETA renderizada; el id real
            # del nodo (species) no se toca en ningún lado.
            svg_uri = generate_svg_data_uri(display_names.get(species, species), color, shape_species_node, tipo='specie', use_latex_style=True)
            net.add_node(species, shape='image', image=svg_uri, label=" ", size=render_size)
        else:
            net.add_node(species, shape=shape_species_node, label=species, color=color, 
                         size=render_size, font={'size': 14, 'color': 'black'})
 
    ######################################
    # AGREGAR NODOS DE REACCIÓN
    ######################################
    for reaction in reaction_vector:
        color = reaction_colors.get(reaction, default_reaction_color)
        # Reducimos el tamaño para las reacciones para que se vean más discretas
        render_size = max(10, node_size - 5) if use_latex_style else max(5, node_size - 10)
        
        if use_latex_style:
            # Pasamos tipo='reaction' a la función generadora
            svg_uri = generate_svg_data_uri(reaction, color, shape_reactions_node, tipo='reaction', use_latex_style=True)
            net.add_node(reaction, shape='image', image=svg_uri, label=" ", size=render_size)
        else:
            net.add_node(reaction, shape=shape_reactions_node, label=reaction, color=color, 
                         size=render_size, font={'size': 14, 'color': 'black'})    

    ######################################
    # ARISTAS (Se mantiene igual)
    ######################################
    connections = set() 
    for reaction, (inputs, outputs) in RN_dict.items():
        for species, coef in inputs:
            if coef.is_integer():
                coef = int(coef)
            else:
                coef = float(coef)
                        
            edge_id = (species, reaction)
            default_curve='curvedCCW'
            smooth_type = {'type': default_curve} if curvature else None
            
            if (reaction, species) in connections:
                if coef != 1:
                    net.add_edge(species, reaction, arrows='to', label=str(coef), color=input_edge_color, 
                                 font_color=input_edge_color, smooth=smooth_type)
                else:
                    net.add_edge(species, reaction, arrows='to', color=input_edge_color, 
                                 font_color=input_edge_color, smooth=smooth_type)
            else:
                if coef != 1:
                    net.add_edge(species, reaction, arrows='to', label=str(coef), color=input_edge_color, 
                                 font_color=input_edge_color)
                else:
                    net.add_edge(species, reaction, arrows='to', color=input_edge_color, 
                                 font_color=input_edge_color)

            connections.add(edge_id)

        for species, coef in outputs:
            if coef.is_integer():
                coef = int(coef)
            else:
                coef = float(coef)

            edge_id = (reaction, species)
            default_curve='curvedCCW'
            smooth_type = {'type': default_curve} if curvature else None
            
            if (species, reaction) in connections:
                if coef != 1:
                    net.add_edge(reaction, species, arrows='to', label=str(coef), color=output_edge_color, 
                                 font_color=output_edge_color, smooth=smooth_type)
                else:
                    net.add_edge(reaction, species, arrows='to', color=output_edge_color, 
                                 font_color=output_edge_color, smooth=smooth_type)
            else:
                if coef != 1:
                    net.add_edge(reaction, species, arrows='to', label=str(coef), color=output_edge_color, 
                                 font_color=output_edge_color)
                else:
                    net.add_edge(reaction, species, arrows='to', color=output_edge_color, 
                                 font_color=output_edge_color)

            connections.add(edge_id)
    
    net.html = net.generate_html()
    with open(filename, "w", encoding="utf-8") as f:
        f.write(net.html)
    return filename


###########################################################################################
# Function to visualize the reaction network and open the HTML file in a browser
###########################################################################################
def rn_visualize_html(rn, lst_color_spcs=None, lst_color_reacs=None, 
                 global_species_color=None, global_reaction_color=None,
                 global_input_edge_color=None, global_output_edge_color=None, 
                 node_size=20, shape_species_node='dot', shape_reactions_node='box', 
                 curvature=None, physics_enabled=False, 
                 use_latex_style=False, # NUEVO PARÁMETRO AQUÍ
                 species_display_names=None, # NUEVO: ver rn_get_visualization
                 filename="reaction_network.html"):
    
    visualizations_dir = "visualizations/rn_visualize_html"
    if not os.path.exists(visualizations_dir):
        os.makedirs(visualizations_dir)
    
    full_path = os.path.join(visualizations_dir, filename)
    
    # Pasamos el nuevo parámetro a rn_get_visualization
    rn_get_visualization(
        rn, lst_color_spcs, lst_color_reacs, 
        global_species_color, global_reaction_color,
        global_input_edge_color, global_output_edge_color, 
        node_size=node_size, shape_species_node=shape_species_node, shape_reactions_node=shape_reactions_node,
        curvature=curvature, physics_enabled=physics_enabled, 
        use_latex_style=use_latex_style, # SE PASA AQUÍ
        species_display_names=species_display_names, # SE PASA AQUÍ
        filename=full_path
    )
    
    abs_path = os.path.abspath(full_path)  
    if not os.path.isfile(abs_path):
        print(f"\nFile not found at {abs_path}") 
        raise FileNotFoundError(f"\nThe file {abs_path} was not found. Check if rn_get_visualization generated the file correctly.")
    
    print(f"\nThe visualization in HTML format of the reaction network was saved in:\n{abs_path}\n")
    webbrowser.open(f"file://{abs_path}")


###########################################################################################
# # Plot the Hierarchy
###########################################################################################
# Function to access species sets from nodes
def get_species_from_node(G, node_name):
    """
    Get the species set stored in a node.
    
    Args:
        G (nx.DiGraph): The hierarchy graph
        node_name (str): The name of the node (e.g., "S1")
        
    Returns:
        set: The set of species stored in the node
    """
    if node_name in G:
        return G.nodes[node_name]['species_set']
    else:
        raise ValueError(f"Node {node_name} not found in the graph")

##################################################################
# Function to visualize the hierarchy of sets and save it as an HTML file
def hierarchy_get_visualization_html(
    input_data, 
    node_size=20, 
    node_color="cyan", 
    edge_color="gray", 
    shape_node='dot',
    lst_color_subsets=None,
    node_font_size=14, 
    edge_width=2,  
    use_latex_style=False,  # NUEVO PARÁMETRO
    filename="hierarchy_visualization.html"
):
    """
    Visualizes the containment hierarchy among sets with automatic positions (inverted hierarchy).
    """
    # Convert input data to unique sets
    unique_subsets = []
    for sublist in input_data:
        if set(sublist) not in [set(x) for x in unique_subsets]:
            unique_subsets.append(sublist)

    # Create a list of unique sets for reference
    Set_of_sets = [set(s) for s in unique_subsets]

    # Sort the sets by size (from smallest to largest)
    Set_of_sets.sort(key=lambda x: len(x))
    
    # Create set names based on their level
    Set_names = [f"X{i+1}" for i in range(len(Set_of_sets))]
    
    # Create a dictionary of labels for the nodes
    labels = {f"X{i+1}": f"{', '.join(sorted(s))}" for i, s in enumerate(Set_of_sets)}

    # Create set labels to display when hovering over nodes
    cursor_labels = [
        ', '.join(sorted(list(s))) if s else '∅'  # Use ∅ for empty sets
        for s in Set_of_sets
    ]

    # Initialize the PyVis network
    net = Network(height='100vh', width='100%', notebook=True, directed=True, cdn_resources='in_line') 
    # net = Network(height="750px", width="100%", directed=True, notebook=False)
    net.set_options(f"""
    {{
      "nodes": {{
        "shape": "{shape_node}",  
        "font": {{
        "size": {node_font_size},  
        "align": "center"  
        }},
        "borderWidth": 2,  
        "borderColor": "black"  
      }},
      "edges": {{
        "smooth": false,
        "color": "{edge_color}",
        "width": {edge_width}  
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

    # Assign colors to nodes based on lst_color_subsets
    color_map = {}
    if lst_color_subsets:
        for color, subsets in lst_color_subsets:
            for subset in subsets:
                for i, s in enumerate(Set_of_sets):
                    if s == set(subset):
                        color_map[Set_names[i]] = color

    ######################################
    # BUCLE DE NODOS ACTUALIZADO
    ######################################
    for name, hover_text in zip(Set_names, cursor_labels):
        color = color_map.get(name, node_color)  
        
        if use_latex_style:
            # Reutilizamos la función generadora. Usamos 'specie' para que tengan el tamaño estándar
            svg_uri = generate_svg_data_uri(name, color, shape_node, tipo='specie', use_latex_style=True)
            
            net.add_node(
                name,  
                label=" ",                    # Ocultamos la etiqueta por defecto
                title=hover_text,  
                color=color,  
                size=node_size + 10,          # Aumentamos un poco para la imagen SVG
                shape='image',                # Cambiamos la forma a imagen
                image=svg_uri                 # Inyectamos el SVG con el subíndice
            )
        else:
            net.add_node(
                name,  
                label=name,  
                title=hover_text,  
                color=color,  
                size=node_size,  
                font={"size": node_font_size},  
                shape=shape_node  
            )

    # Add edges based on containment relationships
    for i, child_set in enumerate(Set_of_sets):
        for j, parent_set in enumerate(Set_of_sets):
            if i < j and child_set.issubset(parent_set):  
                is_direct = True
                for k, intermediate_set in enumerate(Set_of_sets):
                    if i < k < j and child_set.issubset(intermediate_set) and intermediate_set.issubset(parent_set):
                        is_direct = False
                        break
                if is_direct:
                    net.add_edge(Set_names[i], Set_names[j])

    # Print tuple with labels and subsets
    label_set_pairs = [(label, s) for label, s in zip(Set_names, Set_of_sets)]
    print(f"\nTuples of Labels and Subsets:\n{label_set_pairs}")

    # Save the visualization to an HTML file 
    net.html = net.generate_html() 
    with open(filename, "w", encoding="utf-8") as f:
        f.write(net.html)
    return net, label_set_pairs

##################################################################
# Function to visualize the hierarchy of sets and open the HTML file in a browser   
import os
import webbrowser

def hierarchy_visualize_html(
    input_data, 
    node_size=20, 
    node_color="cyan", 
    edge_color="gray", 
    shape_node='dot',    
    lst_color_subsets=None, 
    node_font_size=14, 
    edge_width=2,  
    use_latex_style=False,  # NUEVO PARÁMETRO AQUÍ
    filename="hierarchy_visualization.html"
):        
    """
    Wrapper function to generate and visualize the containment hierarchy among sets as an HTML file.
    """
    # Create the directory structure if it doesn't exist
    target_dir = "visualizations/hierarchy_visualize_html"
    if not os.path.exists(target_dir):
        os.makedirs(target_dir, exist_ok=True)
    
    # Create the full path including the directory structure
    full_path = os.path.join(target_dir, filename)
    
    # Call the hierarchy_get_visualization_html function to generate the HTML file with the additional parameters
    hierarchy_get_visualization_html(
        input_data,  
        node_size=node_size, 
        node_color=node_color, 
        edge_color=edge_color, 
        shape_node=shape_node,
        lst_color_subsets=lst_color_subsets,  
        node_font_size=node_font_size, 
        edge_width=edge_width,  
        use_latex_style=use_latex_style,  # PASAMOS EL PARÁMETRO A LA FUNCIÓN INTERNA
        filename=full_path
    )
    
    # Convert to an absolute path
    abs_path = os.path.abspath(full_path) 
    
    # Check if the file was created correctly
    if not os.path.isfile(abs_path):
        print(f"File not found at {abs_path}")  # Additional message for debugging
        raise FileNotFoundError(f"The file {abs_path} was not found. Check if hierarchy_get_visualization_html generated the file correctly.")
    
    # Inform the user about the file's location
    print(f"\nThe hierarchy visualization was saved to:\n{abs_path}\n")
    
    # Open the HTML file in the default browser
    webbrowser.open(f"file://{abs_path}")

def rn_get_string(rn):
    # Print basic information about the loaded network
    print(f"Loaded reaction network with:")
    print(f"  - {len(rn.species())} species")
    print(f"  - {len(rn.reactions())} reactions")

    # Print the species
    print("\nSpecies:")
    for species in rn.species():
        print(f"  - {species.name}")

    # Print the reactions
    print("\nReactions:")
    for reaction in rn.reactions():
        support = " + ".join([f"{edge.coefficient}*{edge.species_name}" if edge.coefficient != 1 else edge.species_name 
                                for edge in reaction.support_edges()])
        products = " + ".join([f"{edge.coefficient}*{edge.species_name}" if edge.coefficient != 1 else edge.species_name 
                                for edge in reaction.products_edges()])
        
        if not support:
            support = "∅"  # Empty set symbol for inflow reactions
        if not products:
            products = "∅"  # Empty set symbol for outflow reactions
            
        print(f"  - {reaction.name()}: {support} => {products}")
    def rn_get_reaction_string(rn,reaction):
        support = " + ".join([f"{edge.coefficient}*{edge.species_name}" if edge.coefficient != 1 else edge.species_name 
                                for edge in reaction.support_edges()])
        products = " + ".join([f"{edge.coefficient}*{edge.species_name}" if edge.coefficient != 1 else edge.species_name 
                                for edge in reaction.products_edges()])
        
        if not support:
            support = "∅"  # Empty set symbol for inflow reactions
        if not products:
            products = "∅"  # Empty set symbol for outflow reactions
            
        print(f"  - {reaction.name()}: {support} => {products}")

##################################################################
# # Plot a bipartite metabolic network graph from a ReactionNetwork object
##################################################################
# Function to create a bipartite graph from a ReactionNetwork object
def create_bipartite_graph_from_rn(rn):
    graph = rx.PyDiGraph()
    metabolite_nodes = {}
    reaction_nodes = {} 

    for specie in rn.species():
        idx = graph.add_node(('specie', specie.name))
        metabolite_nodes[specie.name] = idx

    for reaction in rn.reactions():
        idx = graph.add_node(('reaction', reaction.name()))
        reaction_nodes[reaction.name()] = idx 

    # Agregar aristas según edges de cada reacción
    for reaction in rn.reactions():
        rxn_name = reaction.name()
        rxn_node = reaction_nodes[rxn_name]

        for edge in reaction.edges:
            met_name = edge.species_name
            coeff = edge.coefficient
            if edge.type == 'reactant':
                graph.add_edge(metabolite_nodes[met_name], rxn_node, coeff)
            elif edge.type == 'product':
                graph.add_edge(rxn_node, metabolite_nodes[met_name], coeff)

    return graph, metabolite_nodes, reaction_nodes

# Function to plot the bipartite graph using Graphviz library
import os
from graphviz import Digraph
from IPython.display import Image

def rn_visualize_png_in_out(
    graph,
    lst_color_spcs=None,
    lst_color_reacs=None,
    global_species_color='cyan',
    global_reaction_color='lightgray',
    global_input_edge_color='red',
    global_output_edge_color='green',
    node_size=20,
    shape_species_node='circle',
    shape_reactions_node='box',
    filename="rn_visualize_png_in_out"  # Quita la extensión .png aquí
):
    # Crear el directorio si no existe
    output_dir = "visualizations/rn_visualize_png_in_out"
    os.makedirs(output_dir, exist_ok=True)
    
    # Construir la ruta completa del archivo (sin extensión)
    filepath = os.path.join(output_dir, filename)
    
    dot = Digraph(comment="Bipartite Network")

    # Crear diccionarios para color específico por nombre
    species_colors = {species: color for color, species_list in (lst_color_spcs or []) for species in species_list}
    reaction_colors = {reaction: color for color, reaction_list in (lst_color_reacs or []) for reaction in reaction_list}

    for idx, (tipo, nombre) in enumerate(graph.nodes()):
        if tipo == 'specie':
            color = species_colors.get(nombre, global_species_color)
            shape = shape_species_node
        else:
            color = reaction_colors.get(nombre, global_reaction_color)
            shape = shape_reactions_node

        dot.node(
            str(idx),
            nombre,
            shape=shape,
            style='filled',
            fillcolor=color,
            width=str(node_size / 72),  # Aproximación para escalar (Graphviz usa pulgadas)
            fontsize='10'
        )

    for src, dst in graph.edge_list():
        data = graph.get_edge_data(src, dst)
        label = str(data)

        # Determinar tipo de arista
        src_tipo = graph.nodes()[src][0]
        color = global_input_edge_color if src_tipo == 'specie' else global_output_edge_color

        # No mostrar etiqueta si peso = 1
        if label == '1':
            dot.edge(str(src), str(dst), color=color)
        else:
            dot.edge(str(src), str(dst), label=label, color=color)

    # Renderizar en la ruta especificada (graphviz añadirá automáticamente .png)
    dot.render(filepath, format='png', cleanup=True)
    full_path = os.path.abspath(f"{filepath}.png")
    print(f"Reaction network saved as: {full_path}")

    return Image(f"{filepath}.png")

# ##################################################################
# # get_rn_visualize_html_in_out
# ################################################################## 
# Install necessary libraries for visualization
# pip install pyvis
# pip install rustworkx
# pip install networkx
# pip install pydot
def get_rn_visualize_html_in_out(graph, lst_color_spcs=None, lst_color_reacs=None, 
                         global_species_color=None, global_reaction_color=None,
                         global_input_edge_color=None, global_output_edge_color=None, 
                         node_size=20, shape_species_node='dot', shape_reactions_node='box', 
                         curvature=None, physics_enabled=False, 
                         use_latex_style=False, # NUEVO PARÁMETRO
                         species_display_names=None, # NUEVO: ver rn_get_visualization
                         filename="metabolic_network.html"):
    """
    (Docstring original...)
    """
    net = Network(height='100vh', width='100%', notebook=True, directed=True, cdn_resources='in_line') 
    
    if physics_enabled:
        net.barnes_hut()
    else:
        net.toggle_physics(False)
    
    default_species_color = global_species_color or 'cyan'
    default_reaction_color = global_reaction_color or 'lightgray'
    input_edge_color = global_input_edge_color or 'red'
    output_edge_color = global_output_edge_color or 'green'
    
    species_colors = {species: color for color, species_list in (lst_color_spcs or []) for species in species_list}
    reaction_colors = {reaction: color for color, reaction_list in (lst_color_reacs or []) for reaction in reaction_list}
    
    g_nx = nx.DiGraph()
    node_mapping = {}  
    species_set = set()
    reaction_set = set()
    
    for idx, (tipo, nombre) in enumerate(graph.nodes()):
        g_nx.add_node(idx, tipo=tipo, nombre=nombre)
        node_mapping[idx] = nombre
        if tipo == 'specie':
            species_set.add(nombre)
        else:
            reaction_set.add(nombre)

    if lst_color_spcs:
        for color, species_list in lst_color_spcs:
            for species in species_list:
                if species not in species_set:
                    print(f"Warning: The species '{species}' specified in lst_color_spcs does not belong to the species of the network.")

    if lst_color_reacs:
        for color, reaction_list in lst_color_reacs:
            for reaction in reaction_list:
                if reaction not in reaction_set:
                    print(f"Warning: The reaction '{reaction}' specified in lst_color_reacs does not belong to the network reactions.")

    edge_counts = Counter((src, dst) for src, dst in graph.edge_list())
    
    for src, dst in graph.edge_list():
        g_nx.add_edge(src, dst, weight=graph.get_edge_data(src, dst))

    if not physics_enabled:
        try:
            pos = nx.nx_pydot.graphviz_layout(g_nx, prog="dot")
        except:
            pos = nx.spring_layout(g_nx)
    else:
        pos = {}

    ######################################
    # BUCLE DE NODOS ACTUALIZADO
    ######################################
    display_names = species_display_names or {}
    for idx in g_nx.nodes():
        tipo = g_nx.nodes[idx]['tipo']
        nombre = g_nx.nodes[idx]['nombre']
        
        # Ajustamos el tamaño base de renderizado según si es especie o reacción
        if tipo == 'specie':
            color = species_colors.get(nombre, default_species_color)
            shape = shape_species_node
            # Como la imagen SVG ya tiene la proporción correcta para el texto abajo, 
            # usamos el tamaño original del nodo sin aumentarlo.
            render_size = node_size 
        else:
            color = reaction_colors.get(nombre, default_reaction_color)
            shape = shape_reactions_node
            # Reducimos el tamaño de las reacciones en PyVis
            render_size = max(10, node_size - 5) if use_latex_style else max(5, node_size - 10)
        
        # display_names solo aplica a especies; el id real del grafo (nombre)
        # no se toca en ningún lado, solo la etiqueta pasada al SVG.
        render_name = display_names.get(nombre, nombre) if tipo == 'specie' else nombre

        if use_latex_style:
            # Pasamos la variable 'tipo' explícitamente a la función SVG
            svg_uri = generate_svg_data_uri(render_name, color, shape, tipo=tipo, use_latex_style=True)
            node_shape = 'image'
            display_label = " " 
        else:
            svg_uri = None
            node_shape = shape
            display_label = nombre
        
        if not physics_enabled and idx in pos:
            x, y = pos[idx]
            if use_latex_style:
                net.add_node(n_id=idx, shape=node_shape, image=svg_uri, label=display_label, 
                             size=render_size, x=x, y=-y, physics=False)
            else:
                net.add_node(n_id=idx, label=display_label, shape=node_shape, color=color, 
                             size=render_size, x=x, y=-y, physics=False, 
                             font={'size': 14, 'color': 'black'})
        else:
            if use_latex_style:
                net.add_node(n_id=idx, shape=node_shape, image=svg_uri, label=display_label, 
                             size=render_size)
            else:
                net.add_node(n_id=idx, label=display_label, shape=node_shape, color=color, 
                             size=render_size, font={'size': 14, 'color': 'black'})

    ######################################
    # LÓGICA DE ARISTAS (Se mantiene igual)
    ######################################
    edge_usage = Counter()
    connections = set()
    
    for src, dst in g_nx.edges():
        weight = g_nx.edges[src, dst]['weight']
        count = edge_counts[(src, dst)]
        edge_usage[(src, dst)] += 1
        
        src_tipo = g_nx.nodes[src]['tipo']
        dst_tipo = g_nx.nodes[dst]['tipo']
        
        if src_tipo == 'specie' and dst_tipo == 'reaction':
            edge_color = input_edge_color
        elif src_tipo == 'reaction' and dst_tipo == 'specie':
            edge_color = output_edge_color
        else:
            edge_color = 'gray' 
        
        smooth_config = {}

        if curvature:
            if count != 1 or (dst, src) in connections:
                curve_type = "cubicBezier" if edge_usage[(src, dst)] % 2 == 0 else "curvedCCW"
                smooth_config = {
                    "smooth": {
                        "type": curve_type,
                        "roundness": 0.3
                    }
                }
            else:
                smooth_config = {
                    "smooth": {
                        "type": "cubicBezier",
                        "forceDirection": "vertical",
                        "roundness": 0.4
                    }
                }
        elif count != 1 or (dst, src) in g_nx.edges():
            curve_type = "curvedCW" if edge_usage[(src, dst)] % 2 == 0 else "curvedCCW"
            smooth_config = {
                "smooth": {
                    "type": curve_type,
                    "forceDirection": "vertical",
                    "roundness": 0.4
                }
            }
        
        if weight != 1: 
            net.add_edge(src, dst, label=str(weight), arrows='to', 
                         color=edge_color, font_color=edge_color, **smooth_config)
        else: 
            net.add_edge(src, dst, arrows='to', color=edge_color, 
                         font_color=edge_color, **smooth_config)

        connections.add((src, dst)) 
          
    net.html = net.generate_html()
    with open(filename, "w", encoding="utf-8") as f:
        f.write(net.html)
    return filename 

##################################################################
# Modificación en rn_visualize_html_in_out
##################################################################
def rn_visualize_html_in_out(graph, lst_color_spcs=None, lst_color_reacs=None, 
                         global_species_color=None, global_reaction_color=None,
                         global_input_edge_color=None, global_output_edge_color=None, 
                         node_size=20, shape_species_node='dot', shape_reactions_node='box', 
                         curvature=None, physics_enabled=False, 
                         use_latex_style=False, # NUEVO PARÁMETRO
                         species_display_names=None, # NUEVO: ver rn_get_visualization
                         filename="metabolic_network.html"):
    """
    Visualizes an interactive graph object and saves it in the specified directory.
    """
    output_dir = "visualizations/rn_visualize_html_in_out"
    os.makedirs(output_dir, exist_ok=True)
    
    filepath = os.path.join(output_dir, filename)
    
    get_rn_visualize_html_in_out(
        graph, lst_color_spcs, lst_color_reacs, 
        global_species_color, global_reaction_color,
        global_input_edge_color, global_output_edge_color, 
        node_size=node_size, shape_species_node=shape_species_node, shape_reactions_node=shape_reactions_node,
        curvature=curvature, physics_enabled=physics_enabled, 
        use_latex_style=use_latex_style, # PASANDO EL PARÁMETRO
        species_display_names=species_display_names, # PASANDO EL PARÁMETRO
        filename=filepath  
    )
    
    abs_path = os.path.abspath(filepath)  
    
    if not os.path.isfile(abs_path):
        print(f"\nFile not found at {abs_path}") 
        raise FileNotFoundError(f"\nThe file {abs_path} was not found. Check if get_rn_visualize_html_in_out created the file correctly.")
    
    print(f"\nThe visualization of the bipartite graph of the metabolic network was saved in HTML format:\n{abs_path}\n")
    
    webbrowser.open(f"file://{abs_path}")



##################################################################
# Weighted hierarchy: Hasse diagram with directional reachability
# edge labels.
##################################################################

def _reachability_to_color(r):
    """Map r in [0,1] to a hex color: low=blue, mid=yellow, high=green."""
    r = max(0.0, min(1.0, r))
    if r < 0.5:
        t = r / 0.5
        red   = int(68  + t * (220 - 68))
        green = int(68  + t * (200 - 68))
        blue  = int(200 + t * (60  - 200))
    else:
        t = (r - 0.5) / 0.5
        red   = int(220 - t * (220 - 40))
        green = int(200 - t * (200 - 160))
        blue  = int(60  - t * (60  - 40))
    return f'#{red:02X}{green:02X}{blue:02X}'


def hierarchy_get_visualization_html_weighted(
    input_data,
    edge_weights=None,
    node_size=20,
    node_color="cyan",
    edge_color="gray",
    shape_node='dot',
    lst_color_subsets=None,
    node_font_size=14,
    edge_width=2,
    filename="hierarchy_weighted.html"
):
    """
    Hasse-diagram visualisation identical to hierarchy_get_visualization_html
    but with directed reachability labels on every cover edge.

    Parameters
    ----------
    input_data : list of sets / frozensets
        The organisations to display.
    edge_weights : dict, optional
        Maps ``(frozenset_A, frozenset_B)`` — where A ⊂ B — to a tuple
        ``(r_AtoB, r_BtoA)`` of reachability coefficients in [0, 1].
        A is the smaller set (lower node), B the larger (upper node).
        Edge label: "↑{r_AtoB:.2f} ↓{r_BtoA:.2f}".
        Edge width and colour are scaled by ``max(r_AtoB, r_BtoA)``.
    All other parameters identical to hierarchy_get_visualization_html.

    Returns
    -------
    (net, label_set_pairs) — same as the unweighted version.
    """
    unique_subsets = []
    for sublist in input_data:
        if set(sublist) not in [set(x) for x in unique_subsets]:
            unique_subsets.append(sublist)

    Set_of_sets = [set(s) for s in unique_subsets]
    Set_of_sets.sort(key=lambda x: len(x))
    Set_names = [f"O{i+1}" for i in range(len(Set_of_sets))]
    labels = {f"O{i+1}": f"{', '.join(sorted(s))}" for i, s in enumerate(Set_of_sets)}
    cursor_labels = [
        ', '.join(sorted(list(s))) if s else '∅'
        for s in Set_of_sets
    ]

    net = Network(height="750px", width="100%", directed=True, notebook=False)
    net.set_options(f"""
    {{
      "nodes": {{
        "shape": "{shape_node}",
        "font": {{"size": {node_font_size}, "align": "center"}},
        "borderWidth": 2
      }},
      "edges": {{
        "smooth": false,
        "font": {{"size": 11, "align": "middle"}},
        "arrows": {{"to": {{"enabled": true, "scaleFactor": 0.6}}}}
      }},
      "physics": {{
        "enabled": false,
        "stabilization": {{"enabled": false}}
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

    color_map = {}
    if lst_color_subsets:
        for color, subsets in lst_color_subsets:
            for subset in subsets:
                for i, s in enumerate(Set_of_sets):
                    if s == set(subset):
                        color_map[Set_names[i]] = color

    for name, hover_text in zip(Set_names, cursor_labels):
        color = color_map.get(name, node_color)
        net.add_node(
            name, label=name, title=hover_text,
            color=color, size=node_size,
            font={"size": node_font_size}, shape=shape_node
        )

    for i, child_set in enumerate(Set_of_sets):
        for j, parent_set in enumerate(Set_of_sets):
            if i >= j or not child_set.issubset(parent_set):
                continue
            is_direct = not any(
                i < k < j
                and child_set.issubset(Set_of_sets[k])
                and Set_of_sets[k].issubset(parent_set)
                for k in range(len(Set_of_sets))
            )
            if not is_direct:
                continue

            key = (frozenset(child_set), frozenset(parent_set))
            if edge_weights and key in edge_weights:
                w = edge_weights[key]
                if len(w) == 6:
                    r_fwd, r_min_fwd, r_max_fwd, r_bwd, r_min_bwd, r_max_bwd = w
                    r_avg  = (r_fwd + r_bwd) / 2.0
                    label  = (f"↑{r_fwd:.2f}[{r_min_fwd:.2f}–{r_max_fwd:.2f}]\n"
                              f"↓{r_bwd:.2f}[{r_min_bwd:.2f}–{r_max_bwd:.2f}]")
                    title  = (f"{Set_names[i]}→{Set_names[j]}: avg={r_fwd:.3f} "
                              f"[{r_min_fwd:.3f}–{r_max_fwd:.3f}]\n"
                              f"{Set_names[j]}→{Set_names[i]}: avg={r_bwd:.3f} "
                              f"[{r_min_bwd:.3f}–{r_max_bwd:.3f}]")
                else:
                    r_fwd, r_bwd = w
                    r_avg  = (r_fwd + r_bwd) / 2.0
                    label  = f"↑{r_fwd:.2f} ↓{r_bwd:.2f}"
                    title  = (f"{Set_names[i]}→{Set_names[j]}: {r_fwd:.3f}\n"
                              f"{Set_names[j]}→{Set_names[i]}: {r_bwd:.3f}")
                ecolor = _reachability_to_color(r_avg)
                ewidth = 1.0 + 5.0 * r_avg
            else:
                label  = ""
                title  = ""
                ecolor = edge_color
                ewidth = float(edge_width)

            net.add_edge(
                Set_names[i], Set_names[j],
                label=label, title=title,
                color=ecolor, width=ewidth
            )

    label_set_pairs = [(label, s) for label, s in zip(Set_names, Set_of_sets)]
    print(f"\nTuples of Labels and Subsets:\n{label_set_pairs}")

    net.html = net.generate_html()
    with open(filename, "w", encoding="utf-8") as f:
        f.write(net.html)
    return net, label_set_pairs


def hierarchy_visualize_html_weighted(
    input_data,
    edge_weights=None,
    node_size=20,
    node_color="cyan",
    edge_color="gray",
    shape_node='dot',
    lst_color_subsets=None,
    node_font_size=14,
    edge_width=2,
    filename="hierarchy_weighted.html"
):
    """
    Wrapper: generate weighted Hasse diagram and open it in the default browser.

    Parameters identical to hierarchy_get_visualization_html_weighted.
    File is saved under visualizations/hierarchy_visualize_html/<filename>.
    """
    target_dir = "visualizations/hierarchy_visualize_html"
    os.makedirs(target_dir, exist_ok=True)
    full_path = os.path.join(target_dir, filename)

    hierarchy_get_visualization_html_weighted(
        input_data,
        edge_weights=edge_weights,
        node_size=node_size,
        node_color=node_color,
        edge_color=edge_color,
        shape_node=shape_node,
        lst_color_subsets=lst_color_subsets,
        node_font_size=node_font_size,
        edge_width=edge_width,
        filename=full_path,
    )

    abs_path = os.path.abspath(full_path)
    if not os.path.isfile(abs_path):
        raise FileNotFoundError(
            f"Weighted hierarchy file not found at {abs_path}."
        )
    print(f"\nWeighted hierarchy visualization saved to:\n{abs_path}\n")
    webbrowser.open(f"file://{abs_path}")

    # Open the HTML file in the default browser
    webbrowser.open(f"file://{abs_path}")