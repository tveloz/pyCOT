import pytest
import numpy as np
from pyCOT.core.rn_rustworkx import ReactionNetwork

@pytest.fixture
def reaction_network():
    rn = ReactionNetwork()
    rn.add_species("A", 10)
    rn.add_species("B", 5)
    rn.add_species("C", 0)
    rn.add_reaction("R1", ["A"], ["B"], [2], [1], rate=1.0)
    rn.add_reaction("R2", ["B"], ["C"], [1], [3], rate=0.5)
    return rn

def test_stoichiometry_matrix_shape(reaction_network):
    matrix = reaction_network.stoichiometry_matrix()
    assert matrix.shape == (3, 2)

def test_stoichiometry_matrix_values(reaction_network):
    matrix = reaction_network.stoichiometry_matrix()
    expected = np.array([
        [-2,  0],
        [ 1, -1],
        [ 0,  3]
    ])
    np.testing.assert_array_equal(matrix, expected)

def test_stoichiometry_matrix_species_order(reaction_network):
    stoich_matrix = reaction_network.stoichiometry_matrix()
    species = stoich_matrix.species
    assert species == ["A", "B", "C"]

def test_stoichiometry_matrix_reaction_order(reaction_network):
    stoich_matrix = reaction_network.stoichiometry_matrix()
    reactions = stoich_matrix.reactions
    assert reactions == ["R1", "R2"]

def test_empty_reaction_network():
    rn = ReactionNetwork()
    matrix = rn.stoichiometry_matrix()
    assert matrix.size == 0
