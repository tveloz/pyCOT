from __future__ import annotations
from collections.abc import Iterable
from itertools import starmap
from numbers import Real
from re import findall, compile
from typing import Literal


import numpy as np
from rustworkx import PyDiGraph, InvalidNode
from pyCOT.io._utils import simplify_terms, str2int_or_float
from pyCOT.rn_types import Species, Reaction, ReactionEdge, ReactionNode, StoichiometryMatrix



class ReactionNetwork(PyDiGraph):

    def __init__(self):
        super().__init__()
        self._species_map = {} 
        self._reaction_map = {}
        
    ########################################################################################################
    #### Basic methods for species #######################################################################
    ########################################################################################################
    def add_species(self, name: str, quantity: Real | None = None) -> int:
        """Add a new species to the reaction network."""
        if self.has_species(name):
            raise ValueError(f"Species '{name}' already exists.")

        new_index = self.add_node(None)
        self[new_index] = Species(new_index, name, quantity)

        self._species_map[name] = new_index

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
        """Obtain the species with the given name."""
        if name not in self._species_map:
            raise InvalidNode(f"Species '{name}' does not exist.")

        index = self._species_map[name]
        return self[index]
    

    def has_species(self, name: str) -> bool:
        try:
            self._species_map[name]
        except KeyError:
            is_present = False
        else:
            is_present = True

        return is_present
    

    def _parse_species_input(self, species: str | Species | Iterable[str] | Iterable[Species]) -> list[Species]:
        if isinstance(species, str):
            species = [species]

        if isinstance(species, Species):
            species = [species]

        if all(isinstance(species_item, str) for species_item in species):
            species = [self.get_species(species_name) for species_name in species]
        return list(species)


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

        self._reaction_map[name] = reaction_node_index

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
        """Obtain the reaction with the given name."""
        if name not in self._reaction_map:
            raise InvalidNode(f"Reaction '{name}' does not exist.")

        index = self._reaction_map[name]
        return Reaction(self[index], self.get_reaction_edges_by_index(index)) # TODO: Es get_reaction_edges_by_index
    

    def get_reaction_by_index(self, index: int) -> Reaction:
        """Obtain the reaction with the given index."""
        return Reaction(self[index], self.get_reaction_edges_by_index(index))


    def has_reaction(self, name: str) -> bool:
        """Checks if the reaction network has a reaction with the given name."""
        try:
            self._reaction_map[name]
        except KeyError:
            is_present = False
        else:
            is_present = True

        return is_present
    

    def add_from_reaction_string(self, string: str) -> int:
        """
        Add a new species and reaction from a string.
        
        Parameters
        ----------
        string : str
            The string to parse.
        
        Returns
        -------
        int
            The index of the added reaction.

        Details
        -------
        The string must be in the following format:
            'reaction_name: coef1*reactant1 + coef2*reactant2 + ... => coef3*product1 + coef4*product2 + ...'
        """
        reaction_string_error = ValueError("Reaction equation must be in the format 'reaction_name: reactant1 + reactant2 + ... => product1 + product2 + ...'")

        parts = string.split(':')

        if len(parts) != 2:
            raise reaction_string_error

        reaction_name = parts[0].strip()
        reaction_equation = parts[1].split("=>")

        if len(reaction_equation) != 2:
            raise reaction_string_error

        if self.has_reaction(reaction_name):
            raise ValueError(f"Reaction '{reaction_name}' already exists in the ReactionNetwork.")
        
        term_regex = compile(r'(\d*(?:\.\d+)?)?\*?([a-zA-Z_]\w*)')
        support_terms = findall(term_regex, reaction_equation[0]) # TODO: Evaluar en soportes o productos vacíos
        support_terms = starmap(lambda coef, term: (str2int_or_float(coef) if coef != '' else 1, term), support_terms)
        support_terms = simplify_terms(support_terms)

        products_terms = findall(term_regex, reaction_equation[1])
        products_terms = starmap(lambda coef, term: (str2int_or_float(coef) if coef != '' else 1, term), products_terms)
        products_terms = simplify_terms(products_terms)

        # Add species to the ReactionNetwork if they don't already exist
        reaction_terms = support_terms + products_terms
        for term in reaction_terms:
            if not self.has_species(term[1]):
                self.add_species(term[1], None)


        # Add reaction to the ReactionNetwork
        support_coefficients = [term[0] for term in support_terms]
        support_species = [term[1] for term in support_terms]

        products_coefficients = [term[0] for term in products_terms]
        products_species = [term[1] for term in products_terms]

        return self.add_reaction(reaction_name, support_species, products_species, support_coefficients, products_coefficients, None)

    def is_active_reaction(self, reaction: str | Reaction) -> bool:
        """
        Check if a reaction is active.

        Parameters
        ----------
        reaction : str | Reaction
            The reaction to check.

        Returns
        -------
        bool
            True if the reaction is active, False otherwise.
        """
        if isinstance(reaction, str):
            reaction = self.get_reaction(reaction)

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
    

    def _parse_reactions_input(self, reactions: str | Reaction | Iterable[str] | Iterable[Reaction]) -> list[Reaction]:
        if isinstance(reactions, str):
            reactions = [reactions]

        if isinstance(reactions, Reaction):
            reactions = [reactions]
        
        if all(isinstance(r, str) for r in reactions):
            reactions = [self.get_reaction(r) for r in reactions]

        return list(reactions)


    ########################################################################################################
    #### Other representations #############################################################################
    ########################################################################################################  

    def _build_matrix(self, include_edge: Literal["both", "reactant", "product"]) -> StoichiometryMatrix:
        """Private helper to build stoichiometry matrices.

        include_edge: "reactant", "product", or "both"
        """
        species = sorted(self.species(), key=lambda s: s.index)
        species_names = [s.name for s in species]
        reactions = sorted(self.reactions(), key=lambda r: r.node.index)
        reactions_names = [r.name() for r in reactions]
        stoich = np.zeros((len(species), len(reactions)), dtype=float)
        species_index_map = {s.index: i for i, s in enumerate(species)}
        for j, reaction in enumerate(reactions):
            for edge in reaction.edges:
                if include_edge == "both" or edge.type == include_edge:
                    idx = edge.source_index if edge.type == "reactant" else edge.target_index
                    val = edge.coefficient
                    if include_edge == "both" and edge.type == "reactant":
                        val = -val
                    stoich[species_index_map[idx], j] += val
        return StoichiometryMatrix(stoich, species_names, reactions_names)
    
    def reactants_matrix(self) -> StoichiometryMatrix:
        """Obtain the stoichiometry matrix of reactants."""
        return self._build_matrix("reactant")
    
    def products_matrix(self) -> StoichiometryMatrix:
        """Obtain the stoichiometry matrix of products."""
        return self._build_matrix("product")    
    
    def stoichiometry_matrix(self) -> StoichiometryMatrix:
        """
        Obtain the stoichiometry matrix of the reaction network.

        Returns
        -------
        StoichiometryMatrix
            The matrix of stoichiometric coefficients. Each row corresponds to a species and each column corresponds to a reaction.
        """
        return self._build_matrix("both")        



    ########################################################################################################
    #### Connectivity queries ##############################################################################
    ########################################################################################################  


    def get_reactions_from_species(self, species: str | Species | Iterable[str] | Iterable[Species]) -> list[Reaction]:
        """
        Obtain the reactions activated by a given species set.

        Parameters
        ----------
        species : str | Species | Iterable[str] | Iterable[Species]
            The species set.

        Returns
        -------
        list[Reaction]
            Reactions potentially activated by the species set (including inflow reactions).
        """
        species = self._parse_species_input(species)
        species_indices = [sp.index for sp in species]
        candidates_indices = set(
            reaction_index 
            for species_index in species_indices
            for reaction_index in self.successor_indices(species_index)
        )
        candidates = [self.get_reaction_by_index(reaction_index) for reaction_index in candidates_indices]
        
        def is_support_present(reaction: Reaction) -> bool:
            return all(reactant_index in species_indices for reactant_index in reaction.support_indices())
        
        # Get reactions with support/reactants present in the given species
        activated_reactions = list(filter(is_support_present, candidates))
        
        # Add all inflow reactions (reactions with no reactants/support)
        inflow_reactions = [reaction for reaction in self.reactions() if not reaction.support_indices()]
        
        # Combine and return unique reactions
        all_reactions = activated_reactions + inflow_reactions
        # Remove duplicates by reaction name
        unique_reactions = []
        seen_names = set()
        for reaction in all_reactions:
            if reaction.name() not in seen_names:
                unique_reactions.append(reaction)
                seen_names.add(reaction.name())
        
        return unique_reactions

    
    def get_reactions_partially_intersecting_support(self, species: str | Species | Iterable[str] | Iterable[Species]) -> list[Reaction]:
         """
        Obtain the reactions partially but not totally activated by a given species set.

        Parameters
        ----------
        species : str | Species | Iterable[str] | Iterable[Species]
            The species set.

        Returns
        -------
        list[Reaction]
            Reactions partially but not totally activated by the species set (excluding inflow and reactions that can be fully triggered).
        """
         species = self._parse_species_input(species)
         species_indices = [sp.index for sp in species]
         candidates_indices = set(
             reaction_index 
             for species_index in species_indices
             for reaction_index in self.successor_indices(species_index))
         candidates = [self.get_reaction_by_index(reaction_index) for reaction_index in candidates_indices]

         def is_support_partially_present(reaction: Reaction) -> bool:
             return any(reactant_index in species_indices for reactant_index in reaction.support_indices()) and not all(reactant_index in species_indices for reactant_index in reaction.support_indices())
        
        # Get reactions with support/reactants partially present in the given species
         activated_reactions = list(filter(is_support_partially_present, candidates))
        
        # Remove inflow reactions (reactions with no reactants/support)
         activated_reactions = [reaction for reaction in activated_reactions if reaction.support_indices()]
        
         return activated_reactions
    def get_supp_from_reactions(self, reactions: str | Reaction | Iterable[str] | Iterable[Reaction]) -> list[Species]:
        """
        Get the list of species produced by a given set of species.

        Parameters
        ----------
        species : str | Species | Iterable[str] | Iterable[Species]
            The species set.

        Returns
        -------
        list[Species]
            A list of species produced by the given set of species.

        Raises
        ------
        ValueError
            If the 'species' parameter is empty or None.
        """
       
        if not reactions:
            raise ValueError("The 'reactions' parameter cannot be empty or None.")
        if not reactions:
            return []  # Return an empty list if reactions is empty or None
        try:
            reactions = self._parse_reactions_input(reactions)
        except Exception as e:
            raise ValueError(f"Invalid reactions input: {e}")
        return [self.get_species_by_index(edge.source_index) for reaction in reactions for edge in reaction.support_edges()]    
  
    def get_prod_from_species(self, species):
        """Get all species produced by reactions that use any of the given species."""
        products = []
        for sp in species:
            for target_idx in self.successor_indices(sp.index):
                try:
                    product = self.get_species_by_index(target_idx)
                    if product not in products:  # avoid duplicates
                        products.append(product)
                except ValueError:
                    continue  # Skip if the node is not a species
        return products

    # def get_reactions_consuming_species(self):
    #     ...


    # def get_reactions_producing_species(self):
    #     ...


    def inflow_species(self) -> list[Species]: # TODO: Add tests
        """Get the list of species produced by the inflow reactions in the reaction network."""
        inflow_reactions = filter(lambda r: r.is_inflow(), self.reactions())
        return self.get_prod_from_reactions([r.name() for r in inflow_reactions])


    def generated_closure(self, species: str | Species | Iterable[str] | Iterable[Species]) -> list[Species]:
        """Obtain the smallest closure set containing a given set of species."""
        species = self._parse_species_input(species)

        species = set(species).union(self.inflow_species())

        while True:
            products = set(self.get_prod_from_species([species.name for species in species]))
            if products.issubset(species):
                break
            species = species.union(products)

        return species
    
    def sub_reaction_network(self, species: str | Species | Iterable[str] | Iterable[Species]) -> ReactionNetwork:
        """Generate a sub-reaction network from the closure of a given set of species."""
        species = self._parse_species_input(species)
        closure = self.generated_closure(species)
        reactions = self.get_reactions_from_species(closure)
        nodes = [sp.index for sp in closure] + [r.node.index for r in reactions]
        subnet = self.subgraph(nodes)
        return ReactionNetwork.from_pydigraph_unsafe(subnet)

    @classmethod
    def from_pydigraph_unsafe(cls, graph: PyDiGraph) -> ReactionNetwork:
        new_net = cls()
        # Add all species using add_species
        for old_idx in graph.node_indices():
            node = graph[old_idx]
            if isinstance(node, Species):
                new_net.add_species(node.name, node.quantity)
        # Collect reaction nodes
        reaction_nodes = {old_idx: graph[old_idx] for old_idx in graph.node_indices() if isinstance(graph[old_idx], ReactionNode)}
        # Map reactants and products for each reaction
        reaction_edges = {old_idx: {'reactants': [], 'products': []} for old_idx in reaction_nodes}
        for u, v, edge in graph.weighted_edge_list():
            if v in reaction_edges:
                reaction_edges[v]['reactants'].append((edge.species_name, edge.coefficient))
            if u in reaction_edges:
                reaction_edges[u]['products'].append((edge.species_name, edge.coefficient))
        # Add reactions using add_reaction
        for old_idx, node in reaction_nodes.items():
            react_list = reaction_edges[old_idx]['reactants']
            prod_list = reaction_edges[old_idx]['products']
            support_species = [name for name, _ in react_list] or None
            support_coeffs = [coeff for _, coeff in react_list]
            products_species = [name for name, _ in prod_list] or None
            products_coeffs = [coeff for _, coeff in prod_list]
            new_net.add_reaction(node.name, support_species, products_species, support_coeffs, products_coeffs, node.rate)
        return new_net