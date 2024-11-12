from collections.abc import Iterable
from itertools import starmap
from numbers import Real
from re import findall, compile

from rustworkx import PyDiGraph, InvalidNode
from pyCOT.io._utils import simplify_terms, str2int_or_float
from pyCOT.rn_types import Species, Reaction, ReactionEdge, ReactionNode



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
    #### Connectivity queries ##############################################################################
    ########################################################################################################  

    def get_reactions_from_species(self, species: str | Species | Iterable[str] | Iterable[Species]) -> list[Reaction]:
        """
        Obtain the reactions potentially activated by a given species set.

        Parameters
        ----------
        species : str | Species | Iterable[str] | Iterable[Species]
            The species set.

        Returns
        -------
        list[Reaction]
            Reactions potentially activated by the species set (i.e. the amount of species and the stoichiometry of the reaction is not considered).
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
        
        return list(filter(is_support_present, candidates))

    
    def get_supp_from_reactions(self, reaction: str | Reaction | Iterable[str] | Iterable[Reaction]) -> list[Species]:
        """
        Obtain the species in the support of a given set of reactions.

        Parameters
        ----------
        reaction : str | Reaction | Iterable[str] | Iterable[Reaction]
            The reaction set.

        Returns
        -------
        list[Species]
        """
        reactions = self._parse_reactions_input(reaction)
        reactants_indices = {
            reactant_index
            for reaction in reactions
            for reactant_index in reaction.support_indices()
        }
        return [self[reactant_index] for reactant_index in reactants_indices]


    def get_prod_from_reactions(self, reaction: str | Reaction | Iterable[str] | Iterable[Reaction]) -> list[Species]:
        """
        Obtain the species in the product sets of a given set of reactions.

        Parameters
        ----------
        reaction : str | Reaction | Iterable[str] | Iterable[Reaction]
            The reaction set.

        Returns
        -------
        list[Species]
        """
        reactions = self._parse_reactions_input(reaction)
        products_indices = {
            product_index
            for reaction in reactions
            for product_index in reaction.products_indices()
        }
        return [self[product_index] for product_index in products_indices]


    def get_species_from_reactions(self, reaction: str | Reaction | Iterable[str] | Iterable[Reaction]) -> list[Species]:
        """Obtain the species in the support and product sets of a given set of reactions."""
        reactions = self._parse_reactions_input(reaction)

        species_indices = {
            reactant_index
            for reaction in reactions
            for reactant_index in reaction.species_indices()
        }

        return [self[species_index] for species_index in species_indices]    

    def get_prod_from_species(self, species: str | Species | Iterable[str] | Iterable[Species]) -> list[Species]:
        """Obtain the species produced by a given set of reactants."""
        potencially_active_reactions = self.get_reactions_from_species(species)

        products_indices = {
            product_index
            for reaction in potencially_active_reactions
            for product_index in reaction.products_indices()
        }

        return [self[product_index] for product_index in products_indices]
    

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