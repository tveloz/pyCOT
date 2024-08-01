import pytest
from pyCOT.rn_rustworkx import ReactionEdge, ReactionNetwork, ReactionNode

@pytest.fixture
def rn():
    rn = ReactionNetwork()
    rn.add_species("A", 10)
    rn.add_species("B", 20)
    return rn

def test_add_reaction_valid_input(rn: ReactionNetwork):
    reaction_node_index = rn.add_reaction("R", ["A"], ["B"], [2], [1], 1)

    assert reaction_node_index == 2

    reaction_node = rn[reaction_node_index]
    assert isinstance(reaction_node, ReactionNode)
    assert reaction_node.name == "R"
    assert reaction_node.rate == 1

    assert rn.get_edge_endpoints_by_index(0) == (0, 2)
    reaction_edge = rn.get_edge_data_by_index(0)
    assert isinstance(reaction_edge, ReactionEdge)
    assert reaction_edge.source == 0
    assert reaction_edge.target == 2
    assert reaction_edge.type == "reactant"
    assert reaction_edge.coefficient == 2

    assert rn.get_edge_endpoints_by_index(1) == (2, 1)
    reaction_edge = rn.get_edge_data_by_index(1)
    assert isinstance(reaction_edge, ReactionEdge)
    assert reaction_edge.source == 2
    assert reaction_edge.target == 1
    assert reaction_edge.type == "product"
    assert reaction_edge.coefficient == 1
    

# def test_add_reaction_invalid_name(rn):
#     with pytest.raises(ValueError):
#         rn.add_reaction("R", ["A"], ["B"], [1], [1], 1)
#         rn.add_reaction("R", ["A"], ["B"], [1], [1], 1)

# def test_add_reaction_invalid_reactant(rn):
#     with pytest.raises(ValueError):
#         rn.add_reaction("R", ["C"], ["B"], [1], [1], 1)

# def test_add_reaction_invalid_product(rn):
#     with pytest.raises(ValueError):
#         rn.add_reaction("R", ["A"], ["C"], [1], [1], 1)