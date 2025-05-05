from dataclasses import dataclass
from numbers import Real
from typing import Literal

import numpy as np


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
        return self.node.name  # Changed from self.node.name() to self.node.name


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
    
    def support_coefficients(self) -> list[Real]:
        return [edge.coefficient for edge in self.support_edges()]

    def products_coefficients(self) -> list[Real]:
        return [edge.coefficient for edge in self.products_edges()]

    def coefficients(self) -> list[Real]:
        return [edge.coefficient for edge in self.edges]

    def stoichiometric_coefficients(self) -> list[Real]:
        """Get the stoichiometric coefficients of the reaction ordered by species index."""
        indices = self.species_indices()
        coefficients = self.coefficients()
        sorted_pairs = sorted(zip(indices, coefficients), key=lambda x: x[0])
        return [coeff for _, coeff in sorted_pairs]


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
    


class NamedVector(np.ndarray):
    vector: np.ndarray
    names: list[str]

    def __new__(cls, vector: np.ndarray, names: list[str]):
        obj = np.asarray(vector).view(cls)
        obj.vector = vector
        obj.names = names
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self.vector = getattr(obj, 'vector', None)  # Assigns 'vector' if exists
        self.names = getattr(obj, 'names', [])

    def __getitem__(self, key: str) -> Real:
        if isinstance(key, str):
            return super().__getitem__(self.names.index(key))  # Llama al método de la clase base
        else:
            return super().__getitem__(key)  # Maneja claves no string


class StoichiometryMatrix(np.ndarray):
    species: list[str]
    reactions: list[str]

    def __new__(cls, input_array: np.ndarray, species: list[str], reactions: list[str]):
        obj = np.asarray(input_array).view(cls)
        obj.species = species
        obj.reactions = reactions
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self.species = getattr(obj, 'species', [])
        self.reactions = getattr(obj, 'reactions', [])

    def __repr__(self):
        return f"StoichiometryMatrix({super().__repr__()}, species={self.species}, reactions={self.reactions})"

    def __matmul__(self, other: np.ndarray) -> NamedVector:
        result = super().__matmul__(other)
        return NamedVector(result, names=self.species)

    def __str__(self) -> str:
        def truncate_list(lst):
            return lst[:8] + ["..."] if len(lst) > 10 else lst

        display_species = truncate_list(self.species)
        display_reactions = truncate_list(self.reactions)

        matrix = self  # Changed from self.matrix to self
        if len(display_species) > 9:
            matrix = matrix[:8, :]
        if len(display_reactions) > 9:
            matrix = matrix[:, :8]

        species_header = "\t" + "\t".join(display_reactions)
        rows = []
        for name, row in zip(display_species, matrix):
            row_values = "\t".join(f"{val:.2f}" for val in row)
            rows.append(f"{name}\t{row_values}")

        return f"{species_header}\n" + "\n".join(rows)







