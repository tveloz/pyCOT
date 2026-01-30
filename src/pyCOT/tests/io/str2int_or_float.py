import pytest
from pyCOT.io._utils import str2int_or_float

def test_str2int_with_integer_string():
    # Prueba con un string que representa un entero
    assert str2int_or_float("123") == 123

def test_str2int_with_float_string():
    # Prueba con un string que representa un float
    assert str2int_or_float("123.45") == 123.45

def test_str2int_with_negative_integer_string():
    # Prueba con un string que representa un entero negativo
    assert str2int_or_float("-123") == -123

def test_str2int_with_negative_float_string():
    # Prueba con un string que representa un float negativo
    assert str2int_or_float("-123.45") == -123.45

def test_str2int_with_zero_integer_string():
    # Prueba con un string que representa un cero entero
    assert str2int_or_float("0") == 0

def test_str2int_with_zero_float_string():
    # Prueba con un string que representa un cero flotante
    assert str2int_or_float("0.0") == 0.0

def test_str2int_or_float_dot_decimal():
    # Prueba con un string que representa un punto decimal
    assert str2int_or_float(".456") == 0.456

def test_str2int_with_invalid_string():
    # Prueba con un string que no puede ser convertido a int o float, debería lanzar un ValueError.
    with pytest.raises(ValueError):
        str2int_or_float("abc")

def test_str2int_with_large_number_string():
    # Prueba con un número muy grande
    assert str2int_or_float("12345678901234567890") == 12345678901234567890

def test_str2int_with_large_float_string():
    # Prueba con un número flotante muy grande
    assert str2int_or_float("1.2345678901234567890") == 1.2345678901234567  # Floats pueden perder precisión
