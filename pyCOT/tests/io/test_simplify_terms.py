import pytest
from pyCOT.io._utils import simplify_terms

def test_simplify_terms_basic_sum():
    # Prueba la suma de coeficientes para el mismo término
    terms = [(2, 'H2O'), (3, 'H2O')]
    expected = [(5, 'H2O')]
    assert simplify_terms(terms) == expected

def test_simplify_terms_remove_zero_coefficient():
    # Prueba la eliminación de términos con coeficiente 0
    terms = [(2, 'H2O'), (0, 'CO2'), (3, 'O2')]
    expected = [(2, 'H2O'), (3, 'O2')]
    assert simplify_terms(terms) == expected

def test_simplify_terms_negative_and_positive_sum():
    # Prueba la cancelación de términos cuando la suma da 0
    terms = [(2, 'H2O'), (-2, 'H2O')]
    expected = []
    assert simplify_terms(terms) == expected

def test_simplify_terms_multiple_species():
    # Prueba la función con múltiples especies
    terms = [(1, 'H2O'), (3, 'O2'), (2, 'CO2')]
    expected = [(1, 'H2O'), (3, 'O2'), (2, 'CO2')]
    assert simplify_terms(terms) == expected

def test_simplify_terms_combination():
    # Prueba una combinación de los casos anteriores
    terms = [(2, 'H2O'), (3, 'O2'), (2, 'CO2'), (-2, 'H2O'), (1, 'O2')]
    expected = [(4, 'O2'), (2, 'CO2')]
    assert simplify_terms(terms) == expected

def test_simplify_terms_all_zero_coefficients():
    # Prueba con todos los coeficientes a 0
    terms = [(0, 'H2O'), (0, 'O2'), (0, 'CO2')]
    expected = []
    assert simplify_terms(terms) == expected

def test_simplify_terms_empty_input():
    # Prueba con una entrada vacía
    terms = []
    expected = []
    assert simplify_terms(terms) == expected
