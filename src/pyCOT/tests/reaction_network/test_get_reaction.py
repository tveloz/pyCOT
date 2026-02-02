import pytest
from pyCOT.core.rn_rustworkx import ReactionNetwork, Reaction, InvalidNode

@pytest.fixture
def rn():
    rn = ReactionNetwork()
    rn.add_species("A", 10)
    rn.add_species("B", 20)
    rn.add_reaction("R1", ["A"], ["B"], [1], [1], 1)
    return rn

def test_get_reaction_existing(rn: ReactionNetwork):
    reaction = rn.get_reaction("R1")

    assert isinstance(reaction, Reaction)
    assert reaction.name() == "R1"

    reactant_edge = reaction.support_edges()[0]
    assert reactant_edge.type == "reactant"
    assert reactant_edge.source_index == 0
    assert reactant_edge.target_index == 2
    
    product_edge = reaction.products_edges()[0]
    assert product_edge.type == "product"
    assert product_edge.source_index == 2
    assert product_edge.target_index == 1

def test_get_reaction_not_existing(rn: ReactionNetwork):
    with pytest.raises(InvalidNode):
        rn.get_reaction("R2")