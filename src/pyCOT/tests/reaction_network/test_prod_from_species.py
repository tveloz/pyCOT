import pytest
from rustworkx import InvalidNode

from pyCOT.core.rn_rustworkx import ReactionNetwork

@pytest.fixture
def rn():
    network = ReactionNetwork()

    # Crear especies
    network.add_species("A")
    network.add_species("B")
    network.add_species("C")
    network.add_species("D")
    network.add_species("E")

    # Crear reacciones
    network.add_reaction("R1", ["A"], ["B"], [1], [1])
    network.add_reaction("R2", ["B"], ["C"], [1], [1])
    network.add_reaction("R3", ["A", "C"], ["D"], [1, 1], [1])
    network.add_reaction("R4", ["E"], ["A"], [1], [1])
    network.add_reaction("R5", ["B"], ["D"], [1], [1])

    return network

def test_single_reactant_produces_single_product(rn: ReactionNetwork):
    products = rn.get_prod_from_species("A")
    product_names = [species.name for species in products]
    expected = {"B"}
    assert set(product_names) == expected, f"Esperado: {expected}, Obtenido: {product_names}"

def test_multiple_reactants_produce_multiple_products(rn: ReactionNetwork):
    products = rn.get_prod_from_species(["A", "C"])
    product_names = [species.name for species in products]
    expected = {"B", "D"}
    assert set(product_names) == expected, f"Esperado: {expected}, Obtenido: {product_names}"

def test_reactant_not_involved_in_any_reaction(rn: ReactionNetwork):
    products = rn.get_prod_from_species("D")
    product_names = [species.name for species in products]
    expected = []
    assert product_names == expected, f"Esperado: {expected}, Obtenido: {product_names}"

def test_empty_reactant_list(rn: ReactionNetwork):
    products = rn.get_prod_from_species([])
    assert products == [], "Esperado: [], Obtenido: productos"

def test_reactant_in_multiple_reactions(rn: ReactionNetwork):
    products = rn.get_prod_from_species("B")
    product_names = [species.name for species in products]
    expected = {"C", "D"}
    assert set(product_names) == expected, f"Esperado: {expected}, Obtenido: {product_names}"

def test_non_existent_reactant(rn: ReactionNetwork):
    with pytest.raises(InvalidNode):
        rn.get_prod_from_species("Z")
