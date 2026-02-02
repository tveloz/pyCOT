import pytest
import numpy as np
from pyCOT.core.rn_types import Species, Reaction, ReactionEdge, ReactionNode, StoichiometryMatrix, NamedVector

@pytest.fixture
def species_list():
    return [
        Species(index=0, name='A'),
        Species(index=1, name='B'),
        Species(index=2, name='C')
    ]

@pytest.fixture
def reaction_nodes():
    return [
        ReactionNode(index=0, name='R1', rate=1.0),
        ReactionNode(index=1, name='R2', rate=2.0)
    ]

@pytest.fixture
def reaction_edges():
    return [
        ReactionEdge(index=0, source_index=0, target_index=1, species_name='A', type='reactant', coefficient=2.0),
        ReactionEdge(index=1, source_index=1, target_index=2, species_name='B', type='product', coefficient=1.0),
        ReactionEdge(index=2, source_index=2, target_index=0, species_name='C', type='reactant', coefficient=3.0),
        ReactionEdge(index=3, source_index=0, target_index=2, species_name='A', type='product', coefficient=1.5)
    ]

@pytest.fixture
def reactions(reaction_nodes, reaction_edges):
    return [
        Reaction(node=reaction_nodes[0], edges=[reaction_edges[0], reaction_edges[1]]),
        Reaction(node=reaction_nodes[1], edges=[reaction_edges[2], reaction_edges[3]])
    ]

@pytest.fixture
def stoichiometry_matrix(species_list, reactions):
    matrix = np.array([
        [-2.0, 1.5],  # A: -2 from R1, +1.5 from R2
        [1.0, 0.0],   # B: +1 from R1, 0 from R2
        [0.0, -3.0]    # C: 0 from R1, -3 from R2
    ])
    species_names = [species.name for species in species_list]
    reaction_names = [reaction.name() for reaction in reactions]
    return StoichiometryMatrix(matrix, species_names, reaction_names)  # Ensure only three arguments are passed

def test_stoichiometry_matrix_initialization(stoichiometry_matrix, species_list, reactions):
    assert isinstance(stoichiometry_matrix, StoichiometryMatrix)
    assert stoichiometry_matrix.shape == (3, 2)
    assert stoichiometry_matrix.species == [species.name for species in species_list]
    assert stoichiometry_matrix.reactions == [reaction.name() for reaction in reactions]  # Changed to call the method

def test_stoichiometry_matrix_repr(stoichiometry_matrix):
    repr_str = repr(stoichiometry_matrix)
    assert "StoichiometryMatrix" in repr_str

def test_stoichiometry_matrix_str(stoichiometry_matrix):
    str_output = str(stoichiometry_matrix)
    assert isinstance(str_output, str)
    assert "A" in str_output
    assert "R1" in str_output

def test_stoichiometric_coefficients(stoichiometry_matrix):
    coefficients = stoichiometry_matrix.sum(axis=1)
    expected = [-0.5, 1, -3]
    np.testing.assert_array_almost_equal(coefficients, expected)

def test_stoichiometry_matrix_multiplication(stoichiometry_matrix):
    vector = np.array([1.0, 2.0])
    result = stoichiometry_matrix @ vector
    assert isinstance(result, NamedVector)
    expected_values = stoichiometry_matrix @ vector
    np.testing.assert_array_almost_equal(result.vector, expected_values)
    expected_names = stoichiometry_matrix.species
    assert result.names == expected_names