import pytest
from pyCOT.rn_rustworkx import ReactionNetwork, Species

@pytest.fixture
def rn() -> ReactionNetwork:
    rn = ReactionNetwork()
    rn.add_species("A", 10)
    rn.add_species("B", 20)
    return rn

def test_species_empty_rn():
    rn = ReactionNetwork()
    assert rn.species() == []

def test_species_one_species(rn: ReactionNetwork):
    species_list = rn.species()
    assert len(species_list) == 2
    assert species_list[0].name == "A"
    assert species_list[0].quantity == 10
    assert species_list[1].name == "B"
    assert species_list[1].quantity == 20