import pytest
from pyCOT.core.rn_rustworkx import ReactionNetwork

@pytest.fixture
def rn_A_species():
    rn = ReactionNetwork()
    rn.add_species("A", 10)
    return rn

def test_has_species_true(rn_A_species: ReactionNetwork):
    assert rn_A_species.has_species("A") == True

def test_has_species_false(rn_A_species: ReactionNetwork):
    assert rn_A_species.has_species("B") == False

# def test_has_species_invalid_name(rn_A_species: ReactionNetwork):
#     with pytest.raises(TypeError):
#         rn_A_species.has_species(10)