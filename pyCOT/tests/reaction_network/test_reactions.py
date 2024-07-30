import pytest
from pyCOT.rn_rustworkx import ReactionNetwork, Reaction

@pytest.fixture
def rn():
    rn = ReactionNetwork()
    rn.add_species("A", 10)
    rn.add_species("B", 20)
    rn.add_species("C", 30)
    rn.add_reaction("R1", ["A"], ["B"], [1], [1], 1)
    rn.add_reaction("R2", ["B"], ["C"], [1], [1], 1)
    return rn

def test_reactions(rn: ReactionNetwork):
    reactions = rn.reactions()
    assert len(reactions) == 2
    assert reactions[0].node.name == "R1"
    assert reactions[1].node.name == "R2"

def test_reactions_empty():
    rn = ReactionNetwork()
    reactions = rn.reactions()
    assert len(reactions) == 0