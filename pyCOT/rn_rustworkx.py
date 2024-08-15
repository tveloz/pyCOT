from collections.abc import Collection
from dataclasses import dataclass
from numbers import Real
from typing import Literal

from rustworkx import PyDiGraph, InvalidNode

@dataclass(slots=True)
class Species:
    index: int
    name: str
    quantity: Real | None = None

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
    

    def __str__(self) -> str:
        support_str = " + ".join([edge.print_term() for edge in self.support_edges()])
        products_str = " + ".join([edge.print_term() for edge in self.products_edges()])
        return f"Reaction {self.name()} (rate = {self.node.rate}): {support_str} -> {products_str}"


class ReactionNetwork(PyDiGraph):

    ########################################################################################################
    #### Basic methods for species #######################################################################
    ########################################################################################################
    def add_species(self, name: str, quantity: Real | None = None) -> int:
        """Add a new species to the reaction network."""
        if self.has_species(name):
            raise ValueError(f"Species '{name}' already exists.")

        new_index = self.add_node(None)
        self[new_index] = Species(new_index, name, quantity)

        return new_index
    

    def species(self) -> list[Species]:
        """Get the list of species in the reaction network."""
        nodes_data = self.nodes()
        return [node_data for node_data in nodes_data if isinstance(node_data, Species)]


    def get_species_by_index(self, index: int) -> Species:
        species = self[index]
        
        if not isinstance(species, Species):
            raise ValueError(f"Node at index {index} is not a species.")
        else:
            return species


    def get_species(self, name: str) -> Species:
        def filter_species_by_name(payload: Species | ReactionNode) -> bool:
            return isinstance(payload, Species) and payload.name == name
        
        indices = self.filter_nodes(filter_species_by_name)

        if len(indices) == 0:
            raise InvalidNode(f"Species '{name}' does not exist.")

        if len(indices) > 1:  # TODO: Estudiar cómo evitar llegar a este punto
            raise ValueError(f"Malformed reaction network. Species '{name}' is not unique.")

        return self.get_species_by_index(indices[0])
    

    def has_species(self, name: str) -> bool:
        try:
            self.get_species(name)
        except InvalidNode:
            is_present = False
        else:
            is_present = True

        return is_present


    ########################################################################################################
    #### Basic methods for reactions #######################################################################
    ########################################################################################################  
    def add_reaction(
            self, name: str, support: list[str] | None, products: list[str] | None, 
            support_coefficients: list[Real], products_coefficients: list[Real], rate: Real = None
    ) -> int:
        """Add a new reaction to the reaction network."""
        if self.has_reaction(name):
            raise ValueError(f"Reaction '{name}' already exists.")
        reaction_node_index = self.add_node(None)
        self[reaction_node_index] = ReactionNode(reaction_node_index, name, rate)

        if support is not None:
            for reactant in support:
                if not self.has_species(reactant):
                    raise ValueError(f"Reactant '{reactant}' must be a declared species.")
        
            support_indices = [self.get_species(reactant).index for reactant in support]
            n_reactants = len(support)

            support_edges_indices = self.add_edges_from(list(zip(support_indices, [reaction_node_index] * n_reactants, [None] * n_reactants)))
            for i, edge_index in enumerate(support_edges_indices):
                edge_data = ReactionEdge(
                    index = edge_index, 
                    source_index = support_indices[i], 
                    target_index = reaction_node_index,
                    species_name = support[i],
                    type = "reactant", 
                    coefficient = support_coefficients[i]
                )
                self.update_edge_by_index(edge_index, edge_data)

        if products is not None:
            for product in products:
                if not None and not self.has_species(product):
                    raise ValueError(f"Product '{product}' must be a declared species.")

            products_indices = [self.get_species(product).index for product in products]
            n_products = len(products)


            products_edges_indices = self.add_edges_from(list(zip([reaction_node_index] * n_products, products_indices, [None] * n_products)))
            for i, edge_index in enumerate(products_edges_indices):
                edge_data = ReactionEdge(
                    index = edge_index, 
                    source_index = reaction_node_index, 
                    target_index = products_indices[i],
                    species_name = products[i], 
                    type = "product", 
                    coefficient = products_coefficients[i]
                )
                self.update_edge_by_index(edge_index, edge_data)

        return reaction_node_index
    

    def get_reaction_edges_by_index(self, reaction_index: int) -> list[ReactionEdge]:
        """Get the list of reaction edges incident to the reaction with the given index."""
        return [edge[2] for edge in self.incident_edge_index_map(reaction_index, all_edges = True).values()]
    

    def reactions(self) -> list[Reaction]:
        """Get the list of reactions in the reaction network."""
        reaction_nodes = [node_data for node_data in self.nodes() if isinstance(node_data, ReactionNode)]
        return [
            Reaction(reaction_node, self.get_reaction_edges_by_index(reaction_node.index)) 
            for reaction_node in reaction_nodes
        ]


    def get_reaction(self, name: str) -> Reaction:
        def filter_reactions_by_name(payload: Species | ReactionNode) -> bool:
            return isinstance(payload, ReactionNode) and payload.name == name
        
        index = self.filter_nodes(filter_reactions_by_name)

        if len(index) == 0:
            raise InvalidNode(f"Reaction '{name}' does not exist.")

        if len(index) > 1:  # TODO: Estudiar qué evitar llegar a este punto
            raise ValueError(f"Malformed reaction network. Reaction '{name}' is not unique.")

        index = index[0]
        return Reaction(self[index], self.get_reaction_edges_by_index(index))
    

    def get_reaction_by_index(self, index: int) -> Reaction:
        return Reaction(self[index], self.get_reaction_edges_by_index(index))


    def has_reaction(self, name: str) -> bool:
        try:
            self.get_reaction(name)
        except InvalidNode:
            is_present = False
        else:
            is_present = True

        return is_present
    
    

    def is_active_reaction(self, reaction_name: str) -> bool: # TODO: Considerar overloads
        """
        Check if a reaction is active.

        Parameters
        ----------
        reaction : str
            The reaction name.

        Returns
        -------
        bool
            True if the reaction is active, False otherwise.
        """
        reaction = self.get_reaction(reaction_name)

        for edge in reaction.support_edges():
            reactant = self.get_species_by_index(edge.source_index)

            if reactant.quantity is None:
                raise ValueError(f"Reactant '{reactant.name}' does not have a quantity.")
            if reactant.quantity < edge.coefficient:
                active = False
                break
        else:
            active = True

        return active


    ########################################################################################################
    #### Connectivity queries ##############################################################################
    ########################################################################################################  

    def get_reactions_from_species(self, species_names: str | Collection[str]) -> list[Reaction]:
        """
        Obtain the reactions potentially activated by a given species set.

        Parameters
        ----------
        species : str | Collection[str]
            The species set.

        Returns
        -------
        list[Reaction]
            Reactions potentially activated by the species set (i.e. the amount of species and the stoichiometry of the reaction is not considered).
        """
        if isinstance(species_names, str):
            species_names = [species_names]
        
        species_indices = [self.get_species(species_name).index for species_name in species_names]
        candidates_indices = set(
            reaction_index 
            for species_index in species_indices
            for reaction_index in self.successor_indices(species_index)
        )
        candidates = [self.get_reaction_by_index(reaction_index) for reaction_index in candidates_indices]
        
        def is_support_present(reaction: Reaction) -> bool:
            return all(self.get_species_by_index(edge.source_index).name in species_names for edge in reaction.support_edges())
        
        return list(filter(is_support_present, candidates))
    

    def get_supp_from_reactions(self, reaction_names: str | Collection[str]) -> list[Species]:
        """
        Obtain the species in the support of a given set of reactions.

        Parameters
        ----------
        reaction_names : str | Collection[str]
            The reaction set.

        Returns
        -------
        list[Species]
        """
        if isinstance(reaction_names, str):
            reaction_names = [reaction_names]
        
        reactions = (self.get_reaction(reaction_name) for reaction_name in reaction_names)
        reactants_indices = {
            reactant_index
            for reaction in reactions
            for reactant_index in reaction.support_indices()
        }
        return [self[reactant_index] for reactant_index in reactants_indices]


    def get_prod_from_reactions(self, reaction_names: str | Collection[str]) -> list[Species]:
        """
        Obtain the species in the product sets of a given set of reactions.

        Parameters
        ----------
        reaction_names : str | Collection[str]
            The reaction set.

        Returns
        -------
        list[Species]
        """
        if isinstance(reaction_names, str):
            reaction_names = [reaction_names]
        
        reactions = (self.get_reaction(reaction_name) for reaction_name in reaction_names)
        products_indices = {
            product_index
            for reaction in reactions
            for product_index in reaction.products_indices()
        }
        return [self[product_index] for product_index in products_indices]


    def get_species_from_reactions(self, reaction_names: str | Collection[str]) -> list[Species]:
        if isinstance(reaction_names, str):
            reaction_names = [reaction_names]
        
        reactions = (self.get_reaction(reaction_name) for reaction_name in reaction_names)

        species_indices = {
            reactant_index
            for reaction in reactions
            for reactant_index in reaction.species_indices()
        }

        return [self[species_index] for species_index in species_indices]    

    # def get_prod_from_species(self):
    #     ...
    

    # def get_reactions_consuming_species(self):
    #     ...


    # def get_reactions_producing_species(self):
    #     ...