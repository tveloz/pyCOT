import pytest
from rustworkx import InvalidNode
from pyCOT.rn_rustworkx import ReactionNetwork

@pytest.fixture
def rn():
    rn = ReactionNetwork()
    rn.add_species("A", 10)
    rn.add_species("B", 10)
    rn.add_species("C", 10)
    rn.add_reaction("R1", ["A"], ["B"], [1], [1], 1)
    rn.add_reaction("R2", ["B"], ["A"], [1], [1], 1)
    rn.add_reaction("R3", ["A", "C"], ["B"], [1, 1], [1], 1)
    return rn

def test_get_reactions_from_species_single_species(rn: ReactionNetwork):
    reactions = rn.get_reactions_from_species("A")
    assert len(reactions) == 1
    assert reactions[0].name() == "R1"

def test_get_reactions_from_species_multiple_species(rn: ReactionNetwork):
    reactions = rn.get_reactions_from_species(["A", "B"])
    assert len(reactions) == 2
    assert {"R1", "R2"} == set(reaction.name() for reaction in reactions)

def test_get_reactions_from_species_incomplete_support(rn: ReactionNetwork):
    reactions = rn.get_reactions_from_species(["B", "C"])
    assert len(reactions) == 1
    assert reactions[0].name() == "R2"

def test_get_reactions_from_species_nonexistent_species(rn: ReactionNetwork):
    with pytest.raises(InvalidNode):
        rn.get_reactions_from_species("D")

def test_get_reactions_from_species_species(rn: ReactionNetwork):
    species = rn[0]
    reactions = rn.get_reactions_from_species(species)
    assert len(reactions) == 1
    assert reactions[0].name() == "R1"

def test_get_reactions_from_species_multiple_species_instances(rn: ReactionNetwork):
    species = [rn[0], rn[1]]
    reactions = rn.get_reactions_from_species(species)
    assert len(reactions) == 2
    assert {"R1", "R2"} == set(reaction.name() for reaction in reactions)