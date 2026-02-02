import pytest
from pyCOT.core.rn_rustworkx import ReactionNetwork
from pyCOT.core.rn_types import ReactionEdge, Species, ReactionNode

def test_sub_network_single_species():
    rn = ReactionNetwork()
    rn.add_species("A")
    sub_rn = rn.sub_reaction_network("A")
    # Verifica que el subgrafo contiene solo la especie "A"
    species_names = [node.name for node in sub_rn.nodes() if isinstance(node, Species)]
    assert species_names == ["A"]
    # No debe haber reacciones
    reaction_nodes = [node for node in sub_rn.nodes() if isinstance(node, ReactionNode)]
    assert reaction_nodes == []

def test_sub_network_with_reaction():
    rn = ReactionNetwork()
    rn.add_species("A")
    rn.add_species("B")
    # Agrega la reacción: A -> B
    rn.add_reaction("r1", support=["A"], products=["B"], support_coefficients=[1], products_coefficients=[1])
    sub_rn = rn.sub_reaction_network("A")
    # La clausura debe incluir ambas especies: A y B
    species_names = sorted([node.name for node in sub_rn.nodes() if isinstance(node, Species)])
    assert species_names == ["A", "B"]
    # Verifica que la reacción "r1" se incluya en el subgrafo
    reaction_names = [node.name for node in sub_rn.nodes() if isinstance(node, ReactionNode)]
    assert "r1" in reaction_names
    # Verifica que los reaction edgses estén presentes
    reaction_edges = [edge for edge in sub_rn.edges() if isinstance(edge, ReactionEdge)]
    assert len(reaction_edges) == 2
    reactant_edges = [
        edge 
        for reaction in sub_rn.reactions() 
        for edge in reaction.support_edges() 
        if isinstance(edge, ReactionEdge) 
    ]
    assert len(reactant_edges) == 1  
    assert reactant_edges[0].coefficient == 1 # Verifica el coeficiente de la reacción
    assert reactant_edges[0].species_name == "A" # Verifica el nombre de la especie reactante
    product_edges = [
        edge 
        for reaction in sub_rn.reactions() 
        for edge in reaction.products_edges() 
        if isinstance(edge, ReactionEdge) 
    ]
    assert len(product_edges) == 1
    assert product_edges[0].coefficient == 1 # Verifica el coeficiente de la reacción
    assert product_edges[0].species_name == "B" # Verifica el nombre de la especie producto

def test_sub_network_inflow_reaction():
    rn = ReactionNetwork()
    rn.add_species("C")
    # Agrega una reacción de inflow (sin soportes) que produce "C"
    rn.add_reaction("inflow", support=None, products=["C"], support_coefficients=[], products_coefficients=[1])
    sub_rn = rn.sub_reaction_network("C")
    species_names = [node.name for node in sub_rn.nodes() if isinstance(node, Species)]
    assert "C" in species_names
    # Verifica que la reacción de inflow esté presente
    reaction_names = [node.name for node in sub_rn.nodes() if isinstance(node, ReactionNode)]
    assert "inflow" in reaction_names
