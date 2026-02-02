import pytest
from rustworkx import InvalidNode
from pyCOT.core.rn_rustworkx import ReactionNetwork

@pytest.fixture
def rn():
    rn = ReactionNetwork()
    rn.add_species("A", 2)
    rn.add_species("B", 0)
    rn.add_reaction("R1", ["A"], ["B"], [2], [1], 1)
    rn.add_reaction("R2", ["B"], ["A"], [1], [1], 1)
    return rn

def test_is_active_reaction_true(rn: ReactionNetwork):
    assert rn.is_active_reaction("R1") is True

def test_is_active_reaction_false(rn: ReactionNetwork):
    assert rn.is_active_reaction("R2") is False

def test_is_active_reaction_empty_reaction_network(rn: ReactionNetwork):
    with pytest.raises(InvalidNode):
        rn.is_active_reaction("R3")

def test_is_active_reaction_unknown_quantity(rn: ReactionNetwork):
    rn.add_species("C", None)
    rn.add_reaction("R3", ["A", "C"], ["B"], [1, 1], [1], 1)
    with pytest.raises(ValueError):
        rn.is_active_reaction("R3")