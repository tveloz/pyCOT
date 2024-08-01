import pytest
from pyCOT.rn_rustworkx import ReactionEdge, ReactionNetwork

@pytest.fixture
def rn():
    rn = ReactionNetwork()
    rn.add_species("A", 10)
    rn.add_species("B", 20)
    rn.add_species("C", 30)
    rn.add_reaction("R1", ["A", "B"], ["C"], [1, 1], [1])
    return rn

def test_get_reaction_edges_by_index_valid_index(rn: ReactionNetwork):
    reaction_index = 3
    edges = rn.get_reaction_edges_by_index(reaction_index)
    assert len(edges) == 3
    for edge in edges:
        assert isinstance(edge, ReactionEdge)
    

def test_get_reaction_edges_by_index_invalid_index(rn: ReactionNetwork):
    no_edges = rn.get_reaction_edges_by_index(4)
    assert len(no_edges) == 0