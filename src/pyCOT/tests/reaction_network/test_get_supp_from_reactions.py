import pytest
from pyCOT.core.rn_rustworkx import ReactionNetwork, InvalidNode, Species

@pytest.fixture
def rn():
    # create a ReactionNetwork object with some reactions and species
    rn = ReactionNetwork()
    rn.add_species("A", 10)
    rn.add_species("B", 10)
    rn.add_species("C", 10)
    rn.add_species("D", 10)
    rn.add_reaction("R1", ["A", "B"], ["C"], [1, 1], [1])
    rn.add_reaction("R2", ["B", "C"], ["D"], [1, 1], [1])
    return rn

def test_get_supp_from_reactions_one_reaction(rn: ReactionNetwork):
    # test with one reaction
    supp = rn.get_supp_from_reactions("R1")
    assert len(supp) == 2
    assert all(isinstance(species, Species) for species in supp)
    assert {supp[0].name, supp[1].name} == {"A", "B"}

def test_get_supp_from_reactions_multiple_reactions(rn: ReactionNetwork):
    # test with multiple reactions
    supp = rn.get_supp_from_reactions(["R1", "R2"])
    assert len(supp) == 3
    assert all(isinstance(species, Species) for species in supp)
    assert {supp[0].name, supp[1].name, supp[2].name} == {"A", "B", "C"}

def test_get_supp_from_reactions_empty_set(rn: ReactionNetwork):
    # test with an empty set of reactions
    supp = rn.get_supp_from_reactions([])
    assert len(supp) == 0

def test_get_supp_from_reactions_invalid_reaction(rn: ReactionNetwork):
    # test with an invalid reaction
    with pytest.raises(InvalidNode):
        rn.get_supp_from_reactions("R3")