import pytest

from pyCOT.rn_rustworkx import ReactionNetwork, Species

@pytest.fixture
def rn():
    return ReactionNetwork()

def test_add_species_new_species(rn: ReactionNetwork):
    new_index = rn.add_species("A", 10)
    new_species = rn[new_index]
    assert new_index == 0
    assert isinstance(new_species, Species)
    assert new_species.index == 0
    assert new_species.name == "A"
    assert new_species.quantity == 10

def test_add_species_existing_species(rn: ReactionNetwork):
    rn.add_species("A", 10)
    with pytest.raises(ValueError):
        rn.add_species("A", 20)

# def test_add_species_invalid_name(rn: ReactionNetwork):
#     with pytest.raises(ValueError):
#         rn.add_species(10, 10)

# def test_add_species_invalid_quantity(rn: ReactionNetwork):
#     with pytest.raises(ValueError):
#         rn.add_species("A", "10")