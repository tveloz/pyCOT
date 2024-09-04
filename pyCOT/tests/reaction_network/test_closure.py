import pytest

from pyCOT.rn_rustworkx import ReactionNetwork, Species


@pytest.fixture
def rn():
    rn = ReactionNetwork()
    rn.add_species("S1")
    rn.add_species("S2")
    rn.add_species("S3")
    rn.add_species("S4")
    rn.add_species("S5")
    rn.add_species("S6")
    rn.add_species("S7") # Non reactive species
    rn.add_species("S8")
    rn.add_species("S9")

    rn.add_reaction("R1", ["S1"], ["S2", "S3"], [1], [1, 1], 1)
    rn.add_reaction("R2", ["S2"], ["S3"], [1], [1], 1)
    rn.add_reaction("R3", ["S3", "S4"], ["S5"], [1, 1], [1], 1)
    rn.add_reaction("R4", None, ["S4", "S6"], None, [1, 1], 1)
    rn.add_reaction("R5", ["S4", "S8"], ["S9"], [1, 1], [1], 1)

    return rn

def test_closure_bigger(rn: ReactionNetwork):
    """Test a species set that generates a clousure with more species."""
    closure = rn.clousure(["S1", "S5", "S6"])

    assert {s.name for s in closure} == {"S1", "S2", "S3", "S4", "S5", "S6"}

def test_non_reactive(rn: ReactionNetwork):
    """Test closure for a non reactive species."""
    closure = rn.clousure(["S7"])

    assert {s.name for s in closure} == {"S4", "S6", "S7"}