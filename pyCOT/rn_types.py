from dataclasses import dataclass
from numbers import Real
from typing import Literal


@dataclass(slots=True)
class Species:
    index: int
    name: str
    quantity: Real | None = None
    
    def __hash__(self):
        return hash(self.index)

@dataclass(slots=True)
class ReactionNode:
    index: int
    name: str
    rate: Real | None = None

@dataclass(slots=True)
class ReactionEdge:
    index: int
    source_index: int
    target_index: int
    species_name: str
    type: Literal["reactant", "product"]
    coefficient: Real

    def print_term(self, precision: int = 4, decimals: int = 1) -> str:
        coeff_precision = round(self.coefficient, precision)

        if coeff_precision == 1:
            coeff_str = f"{self.species_name}"
        else:
            coeff_str = f"{coeff_precision:.{decimals}g}*{self.species_name}"
        return coeff_str


@dataclass(slots=True)
class Reaction:
    node: ReactionNode 
    edges: list[ReactionEdge]

    def name(self) -> str:
        return self.node.name


    def support_edges(self) -> list[ReactionEdge]:
        return [edge for edge in self.edges if edge.type == "reactant"]
    

    def products_edges(self) -> list[ReactionEdge]:
        return [edge for edge in self.edges if edge.type == "product"]
    

    def support_indices(self) -> list[int]:
        return [edge.source_index for edge in self.support_edges()]


    def support_names(self) -> list[str]:
        return [edge.species_name for edge in self.support_edges()]
    

    def products_indices(self) -> list[int]:
        return [edge.target_index for edge in self.products_edges()]
    

    def products_names(self) -> list[str]:
        return [edge.species_name for edge in self.products_edges()]
    
    def species_indices(self) -> list[int]:
        return [edge.source_index if edge.type == "reactant" else edge.target_index for edge in self.edges]
    

    def species_names(self) -> list[str]:
        return [edge.species_name for edge in self.edges]
    
    # TODO: add *_coefficients getter methods
    def stoichiometric_coefficients(self) -> list[Real]:
        """Get the stoichiometric coefficients of the reaction ordered by species index."""
        indices = [edge.source_index if edge.type == "reactant" else edge.target_index for edge in self.edges]
        coefficients = [edge.coefficient for edge in self.edges]
        sorted_indices_coefficients = sorted(zip(indices, coefficients), key=lambda x: x[0])
        sorted_coefficients = [coeff for _, coeff in sorted_indices_coefficients]
        return sorted_coefficients


    def is_inflow(self) -> bool:
        """Checks if the reaction is an inflow reaction."""
        return all(edge.type == "product" for edge in self.edges)
    
    
    def is_outflow(self) -> bool:
        """Checks if the reaction is an outflow reaction."""
        return all(edge.type == "reactant" for edge in self.edges)


    def __str__(self) -> str:
        support_str = " + ".join([edge.print_term() for edge in self.support_edges()])
        products_str = " + ".join([edge.print_term() for edge in self.products_edges()])
        return f"Reaction {self.name()} (rate = {self.node.rate}): {support_str} -> {products_str}"