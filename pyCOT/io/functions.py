from pyCOT.rn_rustworkx import ReactionNetwork
from pyCOT.io._utils import separate_string, remove_comments
from pyCOT.closure_structure import is_self_maintaining

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
    Loads a ReactionNetwork from a .txt file

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

    R1:	=> l;
    R2:	s1+l => 2s1;
    R3:	s1 => s2;
    R4:	s2 + l=>s1;
    R5:	s2=>
    """
    rn = ReactionNetwork()

    # Read reactions line by line and add them to the ReactionNetwork
    with open(file, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                reactions = separate_string(remove_comments(line))
                reactions = [reaction for reaction in reactions if reaction]
                for reaction in reactions:
                    reaction_name = reaction.split(':', 1)[0]
                    if rn.has_reaction(reaction_name):
                        raise ValueError(f"Reaction '{reaction_name}' already exists in the ReactionNetwork")
                    rn.add_from_reaction_string(reaction)
    return rn

#############################################################
def remove_duplicates(semi_org):
    """
    Remove duplicate sublists from a list of lists.

    Args:
        semi_org (list of list): A list of subsets (lists) that may contain duplicates.

    Returns:
        list of list: The input list with duplicate subsets removed, preserving original order.
    """
    unique_subsets = []
    seen_sets = []  # List to store already seen sets

    for sublist in semi_org:
        current_set = frozenset(sublist)  # Use frozenset to compare contents
        if current_set not in seen_sets:
            seen_sets.append(current_set)
            unique_subsets.append(sublist)  # Keep the original list
 
    semi_organisations = []
    for i, semi_org_new in enumerate(unique_subsets):
        semi_organisations.append(semi_org_new)

    return semi_organisations

#############################################################
# Function to find organisations in a reaction network
# based on the semi-organisations and check if they are self-maintaining
# and return the organisations, vector_process and vector_x
def find_organisations(rn, semi_orgs):
    """
    Find self-maintaining organisations from a list of candidate subsets.

    Args:
        rn: A reaction network object or structure required by `is_self_maintaining`.
        semi_orgs (list of list): A list of candidate semi-organisations.

    Returns:
        tuple:
            organisations (list of list): The self-maintaining organisations.
            vector_process (list of list): The corresponding process vectors.
            vector_x (list of list): The corresponding state vectors.
    """
    semi_orgs = remove_duplicates(semi_orgs)  # Remove duplicates
    semi_orgs.sort(key=lambda x: len(x))  # Sort by size
    print("Semi-organisations without duplicates and sorted by size:")
    for i, semi_org_new in enumerate(semi_orgs):
        print(f"S{i+1}:", semi_org_new)
    print("-" * 70)

    print("find_organisations starting")
    organisations = []
    vector_process = []
    vector_x = []
    for i, semi_organisation in enumerate(semi_orgs):
        if len(semi_organisation) == 0:
            print("Semi-organisation = []\n")
            print("     Empty semi-organisation.")
            continue

        print(f"\nSemi-organisation_{i+1} = {semi_organisation}")
        print(f"Semi-organisation_{i+1} size =", len(semi_organisation))
        res = is_self_maintaining(rn, X=semi_organisation)
        if res[0] == True:
            print("     Is self-maintaining:", res[0])
            organisations.append(semi_organisation)
            vector_process.append([x.tolist() for x in res[1]])
            vector_x.append([x.tolist() for x in res[2]])
        else:
            print("     Is self-maintaining: False")

    return organisations, vector_process, vector_x

#############################################################
# Function to generate a .txt file with the subnetwork
# corresponding to a set of species and reactions
# and save it in a specified folder
def generate_subnetwork_txt(species, reactions, rn, folder_name="Txt_sub_network", file_name="sub_network.txt"):
    """
    Generates a .txt file with the species, reactions, and stoichiometric matrix corresponding to a subnetwork,
    saved inside a specified folder.
    
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