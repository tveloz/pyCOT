import pytest
from pyCOT.rn_rustworkx import ReactionNetwork, Species, InvalidNode

# Fixture to create a ReactionNetwork instance
@pytest.fixture
def reaction_network():
    rn = ReactionNetwork()
    rn.add_species("A", 10)
    rn.add_species("B", 20)
    rn.add_species("C", 30)
    rn.add_reaction("R1", ["A"], ["B"], [1], [1], 1)
    rn.add_reaction("R2", ["B"], ["C"], [1], [1], 1)
    return rn

# Test get_species_from_reactions with a single reaction name
def test_get_species_from_reactions_single_reaction(reaction_network: ReactionNetwork):
    print(reaction_network.get_reaction("R1"))
    species = reaction_network.get_species_from_reactions("R1")
    assert len(species) == 2
    assert all(isinstance(sp, Species) for sp in species)
    assert {sp.name for sp in species} == {"A", "B"}

# Test get_species_from_reactions with multiple reaction names
def test_get_species_from_reactions_multiple_reactions(reaction_network: ReactionNetwork):
    species = reaction_network.get_species_from_reactions(["R1", "R2"])
    assert len(species) == 3
    assert all(isinstance(sp, Species) for sp in species)
    assert {sp.name for sp in species} == {"A", "B", "C"}

# Test get_species_from_reactions with a non-existent reaction name
def test_get_species_from_reactions_non_existent_reaction(reaction_network: ReactionNetwork):
    with pytest.raises(InvalidNode):
        reaction_network.get_species_from_reactions("R3")

# Test get_species_from_reactions with an empty reaction name list
def test_get_species_from_reactions_empty_reaction_list(reaction_network: ReactionNetwork):
    species = reaction_network.get_species_from_reactions([])
    assert len(species) == 0

# Test get_species_from_reactions with a reaction name that is not a string
def test_get_species_from_reactions_invalid_reaction_name(reaction_network: ReactionNetwork):
    with pytest.raises(TypeError):
        reaction_network.get_species_from_reactions(123)