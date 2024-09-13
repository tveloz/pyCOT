from pyCOT.rn_rustworkx import ReactionNetwork
from pyCOT.io._utils import separate_string, remove_comments


def from_string(string: str) -> ReactionNetworK
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

def read_txt(file: str) -> ReactionNetworK
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
