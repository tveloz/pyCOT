import pytest
from pyCOT.rn_rustworkx import ReactionNetwork, Species, InvalidNode

@pytest.fixture
def rn():
    rn = ReactionNetwork()
    rn.add_species("A", 10)
    rn.add_species("B", 20)
    return rn

def test_get_species_existing_species(rn: ReactionNetwork):
    species = rn.get_species("A")
    assert isinstance(species, Species)
    assert species.name == "A"
    assert species.quantity == 10

def test_get_species_nonexistent_species(rn: ReactionNetwork):
    with pytest.raises(InvalidNode):
        rn.get_species("C")