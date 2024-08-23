import pytest
import tempfile
import os
from pyCOT.io.functions import read_txt

@pytest.fixture
def temp_txt_file():
    """Fixture to create a temporary .txt file."""
    temp_file = tempfile.NamedTemporaryFile(delete=False, suffix='.txt')
    yield temp_file
    os.remove(temp_file.name)

def test_read_txt_basic(temp_txt_file):
    """Test basic reaction network parsing."""
    content = """
    R1:    => l;
    R2:    s1+l => 2s1;
    R3:    s1 => s2;
    R4:    s2 + l=>s1;
    R5:    s2=>
    """
    temp_txt_file.write(content.encode('utf-8'))
    temp_txt_file.close()

    rn = read_txt(temp_txt_file.name)
    assert rn.has_species('l')
    assert rn.has_species('s1')
    assert rn.has_species('s2')
    assert rn.has_reaction('R1')
    assert rn.has_reaction('R2')
    assert rn.has_reaction('R3')
    assert rn.has_reaction('R4')
    assert rn.has_reaction('R5')

# def test_read_txt_empty_reaction(temp_txt_file):
#     """Test handling of empty reactions."""
#     content = "R1:    => ;"
#     temp_txt_file.write(content.encode('utf-8'))
#     temp_txt_file.close()

#     with pytest.raises(IndexError): # TODO: define behavior for empty reaction
#         rn = read_txt(temp_txt_file.name)


def test_read_txt_duplicate_reaction(temp_txt_file):
    """Test handling of duplicate reaction names."""
    content = """
    R1:    => l;
    R1:    s1+l => 2s1;
    """
    temp_txt_file.write(content.encode('utf-8'))
    temp_txt_file.close()

    with pytest.raises(ValueError, match="Reaction 'R1' already exists in the ReactionNetwork"):
        read_txt(temp_txt_file.name)

def test_read_txt_comments_and_whitespace(temp_txt_file):
    """Test handling of comments and various whitespace scenarios."""
    content = """
    # This is a comment
    R1:    => l;  # Inline comment
    R2: s1+l => 2s1;
        
    # Another comment
    R3:    s1 => s2;
    """
    temp_txt_file.write(content.encode('utf-8'))
    temp_txt_file.close()

    rn = read_txt(temp_txt_file.name)

    assert rn.has_reaction('R1')
    assert rn.has_reaction('R2')
    assert rn.has_reaction('R3')

def test_read_txt_complex_species_names(temp_txt_file):
    """Test handling of complex species names."""
    content = """
    R1:    s_1 + l => s_2_3;
    R2:    s_2_3 => s_4;
    """
    temp_txt_file.write(content.encode('utf-8'))
    temp_txt_file.close()

    rn = read_txt(temp_txt_file.name)

    assert rn.has_species('s_1')
    assert rn.has_species('l')
    assert rn.has_species('s_2_3')
    assert rn.has_species('s_4')

def test_read_txt_invalid_format(temp_txt_file):
    """Test handling of invalid file format."""
    content = """
    This is an invalid format
    R1:    => l;
    """
    temp_txt_file.write(content.encode('utf-8'))
    temp_txt_file.close()

    with pytest.raises(Exception): # TODO: define a more specific exception
        read_txt(temp_txt_file.name)
