import pytest
from pyCOT.core.rn_rustworkx import ReactionNetwork, Species

@pytest.fixture
def rn():
    rn = ReactionNetwork()
    rn.add_species("A", 10)
    rn.add_species("B", 20)
    rn.add_reaction("R", ["A"], ["B"], [2], [1], 1)
    return rn

def test_get_species_by_index_valid_index(rn: ReactionNetwork):
    species = rn.get_species_by_index(0)
    assert isinstance(species, Species)
    assert species.index == 0
    assert species.name == "A"
    assert species.quantity == 10

def test_get_species_by_index_invalid_species(rn: ReactionNetwork):
    with pytest.raises(ValueError):
        rn.get_species_by_index(2)

def test_get_species_by_index_invalid_index(rn: ReactionNetwork):
    with pytest.raises(IndexError):
        rn.get_species_by_index(3)