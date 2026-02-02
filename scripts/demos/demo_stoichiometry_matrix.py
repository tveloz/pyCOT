import sys
import os

import numpy as np

# Añadir el directorio raíz al PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pyCOT.core.rn_rustworkx import ReactionNetwork
from pyCOT.core.rn_types import StoichiometryMatrix

def main():
    # Initialize ReactionNetwork
    rn = ReactionNetwork()

    # Define species
    rn.add_species('H2')
    rn.add_species('O2')
    rn.add_species('H2O')

    # Define reactions
    rn.add_reaction("R1", ["H2", "O2"], ["H2O"], [2, 1], [2])
    rn.add_reaction("R2", ["H2O", "O2"], ["H2O"], [1, 1], [2])
    rn.add_reaction("R3", ["O2", "H2O"], ["H2"], [1, 1], [2])

    # Generate stoichiometry matrix
    stoich_matrix = rn.stoichiometry_matrix()

    # Print the stoichiometry matrix
    print("Stoichiometry matrix:")
    print(stoich_matrix)

    # vector = np.array([2,])

    # # Multiply the stoichiometry matrix by a vector
    # result = stoich_matrix @ vector

    # print("Print the result for H2: ")
    # print(result["H2"])

if __name__ == "__main__":
    main()