import pytest
from pyCOT.rn_rustworkx import ReactionNetwork, InvalidNode


@pytest.fixture
def rn():
    return ReactionNetwork()

def test_add_from_reaction_string_valid_string(rn: ReactionNetwork):
    string = "R1: 2*A + 3*B => 4*C"
    index = rn.add_from_reaction_string(string)
    assert rn.has_reaction("R1")
    assert index == 3

def test_add_from_reaction_string_mixed_coefficients(rn: ReactionNetwork):
    string = "R1: 0.2*A + .3*B => C"
    index = rn.add_from_reaction_string(string)
    assert rn.has_reaction("R1")
    assert index == 3
    r1 = rn.get_reaction("R1")
    print(r1)
    print(r1.edges)
    assert {edge.coefficient for edge in r1.support_edges()} == {0.2, 0.3}
    assert [edge.coefficient for edge in r1.products_edges()] == [1]

def test_add_from_reaction_string_invalid_string_format(rn: ReactionNetwork):
    string = "R1: 2*A + 3*B => 4*C => 5*D"
    with pytest.raises(ValueError):
        rn.add_from_reaction_string(string)

def test_add_from_reaction_string_reaction_name_already_exists(rn: ReactionNetwork):
    string1 = "R1: 2*A + 3*B => 4*C"
    string2 = "R1: 2*A + 3*B => 4*C"
    rn.add_from_reaction_string(string1)
    with pytest.raises(ValueError):
        rn.add_from_reaction_string(string2)

def test_add_from_reaction_string_empty_string(rn: ReactionNetwork):
    string = ""
    with pytest.raises(ValueError):
        rn.add_from_reaction_string(string)

def test_add_from_reaction_string_simplifiable(rn: ReactionNetwork):
    string = "R1: A + A => B"
    index = rn.add_from_reaction_string(string)
    r1 = rn.get_reaction("R1")
    assert r1.name() == "R1"
    assert index == 2
    assert r1.support_names() == ["A"]
    assert r1.products_names() == ["B"]
    assert r1.support_edges()[0].coefficient == 2
    assert {edge.coefficient for edge in r1.products_edges()} == {1}

def test_add_from_reaction_string_inflow(rn: ReactionNetwork):
    string = "R1: => A"
    index = rn.add_from_reaction_string(string)
    r1 = rn.get_reaction("R1")
    assert r1.name() == "R1"
    assert index == 1
    assert r1.support_names() == []
    assert r1.products_names() == ["A"]
    assert [edge.coefficient for edge in r1.support_edges()] == []
    assert [edge.coefficient for edge in r1.products_edges()] == [1]
