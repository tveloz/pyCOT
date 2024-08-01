import pytest
from pyCOT.rn_rustworkx import ReactionNetwork

@pytest.fixture
def rn():
    rn = ReactionNetwork()
    rn.add_species("A", 10)
    rn.add_species("B", 20)
    rn.add_reaction("R1", ["A"], ["B"], [1], [1], 1)
    return rn

def test_has_reaction_true(rn: ReactionNetwork):
    assert rn.has_reaction("R1") is True

def test_has_reaction_false(rn: ReactionNetwork):
    assert rn.has_reaction("R2") is False

def test_has_reaction_empty_reaction_network():
    rn = ReactionNetwork()
    assert rn.has_reaction("R1") is False