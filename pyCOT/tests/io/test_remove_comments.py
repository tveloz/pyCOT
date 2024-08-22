import pytest
from pyCOT.io._utils import remove_comments


def test_remove_inline_comment():
    input_string = "1+1#sum\n2+2"
    expected_output = "1+1\n2+2"
    assert remove_comments(input_string) == expected_output

def test_remove_full_line_comment():
    input_string = "# This is a comment\n1+1"
    expected_output = "\n1+1"
    assert remove_comments(input_string) == expected_output

def test_remove_comment_with_semicolon():
    input_string = "1+1;2+2#sum\n3+3"
    expected_output = "1+1;2+2\n3+3"
    assert remove_comments(input_string) == expected_output

def test_comment_at_end_of_line():
    input_string = "a=1 # initialize a\nb=2 # initialize b"
    expected_output = "a=1 \nb=2 "
    assert remove_comments(input_string) == expected_output

def test_multiple_comments():
    input_string = "# Comment 1\n1+1 # Comment 2\n# Comment 3\n2+2"
    expected_output = "\n1+1 \n\n2+2"
    assert remove_comments(input_string) == expected_output

def test_no_comments():
    input_string = "1+1\n2+2"
    expected_output = "1+1\n2+2"
    assert remove_comments(input_string) == expected_output

def test_comment_with_semicolon_and_newline():
    input_string = "x=5;y=10 # Set x and y;\nz=15 # Set z\n"
    expected_output = "x=5;y=10 \nz=15 \n"
    assert remove_comments(input_string) == expected_output

def test_empty_string():
    input_string = ""
    expected_output = ""
    assert remove_comments(input_string) == expected_output
