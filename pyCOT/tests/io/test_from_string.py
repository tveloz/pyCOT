import pytest
from pyCOT.io.functions import from_string
from pyCOT.rn_rustworkx import ReactionNetwork

def test_empty_string():
    with pytest.raises(ValueError):
        rn = from_string("")
        print(rn)

def test_single_reaction():
    string = "R1: A => B"
    rn = from_string(string)
    assert isinstance(rn, ReactionNetwork)
    assert rn.has_reaction("R1")

def test_multiple_reactions():
    string = "R1: A => B\nR2: B => C"
    rn = from_string(string)
    assert isinstance(rn, ReactionNetwork)
    assert rn.has_reaction("R1")
    assert rn.has_reaction("R2")

def test_reaction_already_exists():
    string = "R1: A => B\nR1: C => D"
    with pytest.raises(ValueError):
        from_string(string)

def test_invalid_reaction_string():
    string = "Invalid reaction string"
    with pytest.raises(ValueError):
        from_string(string)

def test_inflow_and_outflow():
    string = """
    R1:	=> l;
    R2:	l=>
    """

    rn = from_string(string)
    assert rn.has_reaction("R1")
    assert rn.has_reaction("R2")