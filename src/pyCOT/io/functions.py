from pyCOT.core.rn_rustworkx import ReactionNetwork
from pyCOT.io._utils import separate_string, remove_comments
#from pyCOT.Persistent_Modules import is_self_maintaining

# Import libraries 
import os

def from_string(string: str) -> ReactionNetwork:
    """
    Creates a ReactionNetwork from a string

    Parameters
    ----------
    string : str
        The string representation of the ReactionNetwork

    Returns
    -------
    ReactionNetwork
        The created ReactionNetwork
    """
    reaction_string_error = ValueError("Reaction equation must be in the format 'reaction_name: reactant1 + reactant2 + ... => product1 + product2 + ...'")

    reactions = separate_string(remove_comments(string))
    reactions = [reaction for reaction in reactions if reaction]
    if len(reactions) == 0:
        raise reaction_string_error

    rn = ReactionNetwork()
    for reaction in reactions:
        reaction_splited = reaction.split(':', 1)
        if len(reaction_splited) != 2:
            raise reaction_string_error
        reaction_name = reaction_splited[0].strip()
        if rn.has_reaction(reaction_name):
            raise ValueError(f"Reaction '{reaction_name}' already exists in the ReactionNetwork")
        rn.add_from_reaction_string(reaction)

    return rn

def read_txt(file: str) -> ReactionNetwork:
    """
    Loads a ReactionNetwork from a .txt file with optional comments after ';'.

    Parameters
    ----------
    file : str
        The path to the .txt file.

    Returns
    -------
    ReactionNetwork
        The loaded ReactionNetwork

    Details
    -------
    The format of the .txt file is as the following example:

    R1: => l; Inflow (Entrada de la sustancia l al sistema)
    R2: s1 + l => 2s1; Autocatálisis (s1 se duplica usando l como reactivo)
    R3: s1 => s2;
    R4: s2 + l => s1; Retroconversión (s2 se convierte en s1 usando l)
    R5: s2 => 
    """
    rn = ReactionNetwork()
    reaction_comments = {}  # opcional: guardar comentarios asociados a cada reacción

    # Read reactions line by line and add them to the ReactionNetwork
    with open(file, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                # Separar reacción del comentario (si existe)
                if ";" in line:
                    reaction_part, comment_part = line.split(";", 1)
                    comment_part = comment_part.strip()
                else:
                    reaction_part, comment_part = line, None

                reaction_part = reaction_part.strip()

                # Procesar la reacción
                reactions = separate_string(remove_comments(reaction_part))
                reactions = [reaction for reaction in reactions if reaction]

                for reaction in reactions:
                    reaction_name = reaction.split(':', 1)[0]
                    if rn.has_reaction(reaction_name):
                        raise ValueError(f"Reaction '{reaction_name}' already exists in the ReactionNetwork")
                    
                    rn.add_from_reaction_string(reaction)

                    # Guardar comentario si existe
                    if comment_part:
                        reaction_comments[reaction_name] = comment_part

    # Si quieres asociar los comentarios al ReactionNetwork
    rn.reaction_comments = reaction_comments  

    return rn

#############################################################
# Function for printing the reaction network
def print_reaction_network(rn):
    reactions = rn.reactions()
    for i, reaction in enumerate(reactions):
        name = reaction.name()
        reactants = []
        products = []
        for edge in reaction.edges:
            species = edge.species_name
            coef = edge.coefficient
            entry = f"{coef if coef != 1 else ''}{species}"
            if edge.type == "reactant":
                reactants.append(entry)
            elif edge.type == "product":
                products.append(entry)
        reactant_str = " + ".join(reactants)
        product_str = " + ".join(products)
        ending = ";" if i < len(reactions) - 1 else ""
        print(f"{name}: {reactant_str} => {product_str}{ending}")

#############################################################
# Function to generate a .txt file with the subnetwork 
def generate_subnetwork_txt(species, reactions, rn, folder_name="Txt_sub_network", file_name="sub_network.txt"):
    """
    Generates a .txt file with the species, reactions, rn,
    saved the file inside a specified folder.
    
    Args:
        species (list): list of species names (e.g., ['s1', 's2', 's6'])
        reactions (list): list of reaction labels (e.g., ['R1', 'R3', 'R7'])
        rn: reaction network object
        folder_name (str): name of the folder where the file will be saved
        file_name (str): name of the output file (e.g., 'sub_network.txt')
        
    Returns:
        str: absolute path to the saved file
    """
    # Get the stoichiometry matrix (for access to species and reactions lists)
    stoich_matrix = rn.stoichiometry_matrix()
    
    # Get the lists of all species and reactions
    all_species = stoich_matrix.species
    all_reactions = stoich_matrix.reactions
    
    # For better coefficient extraction, get reactants and products matrices
    reactants_matrix = rn.reactants_matrix()
    products_matrix = rn.products_matrix()
    
    # Create folder if it doesn't exist
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
    
    # Full file path
    file_path = os.path.join(folder_name, file_name)
    
    def format_reaction(reactants, products):
        left = '+'.join(reactants) if reactants else ''
        right = '+'.join(products) if products else ''
        return f"{left}=>{right}"
    
    with open(file_path, 'w') as f:
        for reaction_name in reactions:
            if reaction_name in all_reactions:
                col_idx = all_reactions.index(reaction_name)
                
                reactants = []
                products = []
                
                # Check each species' coefficients for this reaction
                for species_name in species:
                    if species_name in all_species:
                        species_idx = all_species.index(species_name)
                        
                        # Extract reactant coefficient (if any)
                        reactant_coef = int(reactants_matrix.data[species_idx, col_idx])
                        if reactant_coef > 0:
                            reactants.append(f"{reactant_coef if reactant_coef!=1 else ''}{species_name}")
                        
                        # Extract product coefficient (if any)
                        product_coef = int(products_matrix.data[species_idx, col_idx])
                        if product_coef > 0:
                            products.append(f"{product_coef if product_coef!=1 else ''}{species_name}")
                
                reaction_str = format_reaction(reactants, products)
                f.write(f"{reaction_name}:\t{reaction_str};\n")
    
    # Print absolute file path
    abs_path = os.path.abspath(file_path)
    print(f"\nFile saved at:\n{abs_path}")
    return abs_path