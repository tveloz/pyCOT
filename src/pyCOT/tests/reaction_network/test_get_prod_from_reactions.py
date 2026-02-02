import pytest
from pyCOT.core.rn_rustworkx import ReactionNetwork, InvalidNode

@pytest.fixture
def rn():
    rn = ReactionNetwork()
    rn.add_species("s1", 10)
    rn.add_species("s2", 10)
    rn.add_reaction("r1", ["s1"], ["s2"], [1], [1])
    rn.add_reaction("r2", ["s1"], ["s1", "s2"], [1], [1, 1])
    return rn

def test_get_prod_from_reactions_single_reaction(rn: ReactionNetwork):
    products = rn.get_prod_from_reactions("r1")
    assert len(products) == 1
    assert products[0] == rn[1] # index 1 is s2

def test_get_prod_from_reactions_multiple_reactions(rn: ReactionNetwork):
    products = rn.get_prod_from_reactions(["r1", "r2"])
    assert len(products) == 2
    assert rn[0] in products
    assert rn[1] in products

def test_get_prod_from_reactions_invalid_reaction(rn: ReactionNetwork):
    with pytest.raises(InvalidNode):
        rn.get_prod_from_reactions("invalid_reaction")