######################################################################################
import os
import sys
# Añadir el directorio raíz al PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
 

from pyCOT.RN_view import *

#####################################################################################
# Ejemplos

# Ruta del archivo
# file_path = 'Txt/autopoietic.txt'
# file_path = 'Txt/2019fig1.txt'
# file_path = 'Txt/2019fig2.txt'
# file_path = 'Txt/non_connected_example.txt' 
file_path = 'Txt/Farm.txt'

# Listas de conjuntos de especies y reaciones con colores 
lst_color_specs = [("blue", ["farmer","fertilizer"]), ("yellow", ["eggs","chickens"])]
lst_color_reacs = [("purple", ["R1"]), ("orange", ["R2","R3","R4"]), ("green", ["R12","R13","R14"])]

RN, species, reactions = interpret_reactions_from_archive(file_path)


# Visualizar la red con lst_color_spcs y lst_color_reacs aplicados 
RN_view(RN, filename='RN_view.html') 
# RN_view(RN, lst_color_spcs=lst_color_specs, filename="RN_view.html")
# RN_view(RN, lst_color_reacs=lst_color_reacs, filename="RN_view.html")
# RN_view(RN, lst_color_spcs=lst_color_specs, lst_color_reacs=lst_color_reacs, filename="RN_view.html")
# RN_view(RN, global_species_color='blue', global_reaction_color='orange', filename="RN_view.html")
# RN_view(RN,global_species_color='blue', global_reaction_color='orange',global_input_edge_color='blue', global_output_edge_color='gray', filename="RN_view.html")
# RN_view(RN,lst_color_spcs=lst_color_specs, lst_color_reacs=lst_color_reacs, global_input_edge_color='blue', global_output_edge_color='gray', filename="RN_view.html")