"""
types.py — Result types for the decomposition module.
"""
from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class FragileCircuit:
    """
    One fragile circuit Di: an equivalence class of C = X \\ (E ∪ F)
    under dynamical connection.

    species_mask : bitmask of Di's species (cot_gen species-index space).
    reaction_ids : sorted reaction indices of R*_i — reactions of R_X
                   that consume (have in their support) a species of Di.
                   Indices are into `rn.reactions()` / RN.stoichiometry_matrix()
                   columns for the ORIGINAL full network, not a sub-index.
    is_self_maintaining : result of Theorem 2.16's local LP check.
    flux : optimal flux vector over reaction_ids (only when self-maintaining).
    """
    species_mask: int
    reaction_ids: tuple[int, ...]
    is_self_maintaining: bool
    flux: tuple[float, ...] | None = None

    def size(self) -> int:
        return bin(self.species_mask).count('1')


@dataclass
class DecompositionResult:
    """
    Full decomposition of one semi-organization X.

    sp_mask   : X's species mask, in cot_gen's quotiented species-index
                space (E0/native-food species excluded — see bridge.py).
    E_mask    : catalysts (zero net effect in every reaction of R_X).
    F_mask    : maximal simultaneously overproducible set, INCLUDING
                native food species (Prop 5.1a: F0 ⊆ F trivially).
    circuits  : fragile circuits D1..Dm (C = X_full \\ (E ∪ F) partitioned
                by dynamical connection). Empty list ⟺ X is fully
                overproduced/catalytic (organization iff C is empty).
    R_X       : sorted reaction indices with full support ⊆ X_full.
    is_organization : True iff every circuit is self-maintaining
                (Theorem 2.16). Assumes X is already known closed
                (guaranteed for any SO coming out of cot_gen).
    """
    sp_mask: int
    E_mask: int
    F_mask: int
    circuits: list[FragileCircuit]
    R_X: tuple[int, ...]
    is_organization: bool

    def n_species_full(self) -> int:
        return bin(self.sp_mask | self.F_mask | self.E_mask).count('1')

    def summary(self) -> str:
        n_e = bin(self.E_mask).count('1')
        n_f = bin(self.F_mask).count('1')
        sizes = sorted((c.size() for c in self.circuits), reverse=True)
        return (f"|E|={n_e} |F|={n_f} circuits={len(self.circuits)} "
                f"sizes={sizes} organization={self.is_organization}")
