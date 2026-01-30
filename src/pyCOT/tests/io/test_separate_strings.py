import pytest
from pyCOT.io._utils import separate_string

def test_separate_string_empty():
    """Test that an empty string returns an empty list."""
    assert separate_string("") == []

def test_separate_string_single_line():
    """Test that a single line of text returns a list with one element."""
    assert separate_string("Hello world") == ["Hello world"]

def test_separate_string_multiple_lines():
    """Test that multiple lines of text return a list with multiple elements."""
    assert separate_string("Hello world\nThis is a test") == ["Hello world", "This is a test"]

def test_separate_string_multiple_lines_with_semicolons():
    """Test that multiple lines of text with semicolons return a list with multiple elements."""
    assert separate_string("Hello world; This is a test\nAnother line") == ["Hello world", "This is a test", "Another line"]

def test_separate_string_single_line_with_semicolons():
    """Test that a single line of text with semicolons return a list with multiple elements."""
    assert separate_string("Hello world; This is a test; Another line") == ["Hello world", "This is a test", "Another line"]

def test_separate_string_leading_trailing_whitespace():
    """Test that leading and trailing whitespace is ignored."""
    assert separate_string("   Hello world   ") == ["Hello world"]

def test_separate_string_empty_lines():
    """Test that empty lines are ignored."""
    assert separate_string("Hello world\n\nThis is a test") == ["Hello world", "This is a test"]

def test_separate_string_only_semicolons():
    """Test that a string with only semicolons returns an empty list."""
    assert separate_string(";;;") == []

def test_separate_string_only_newlines():
    """Test that a string with only newlines returns an empty list."""
    assert separate_string("\n\n\n") == []