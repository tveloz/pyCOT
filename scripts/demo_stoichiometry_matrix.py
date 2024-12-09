import sys
import os

# Añadir el directorio raíz al PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pyCOT.rn_rustworkx import ReactionNetwork
from pyCOT.rn_types import StoichiometryMatrix

def main():
    # Initialize ReactionNetwork
    rn = ReactionNetwork()

    # Define species
    rn.add_species('H2')
    rn.add_species('O2')
    rn.add_species('H2O')

    # Define reactions
    rn.add_reaction("R1", ["H2", "O2"], ["H2O"], [2, 1], [2])

    # Generate stoichiometry matrix
    stoich_matrix = rn.stoichiometry_matrix()

    # Print the stoichiometry matrix
    print(stoich_matrix)

if __name__ == "__main__":
    main()