"""
types.py — Core data types for cot_gen.

BitSet: immutable, hashable wrapper around Python int (big-int bitset).
        Bit position i is set  ⟺  species i is in the set.

RNData: compiled, bitset-indexed representation of a reaction network.
        Built once (via io_pyCOT) and then treated as read-only.
        All per-species and per-reaction sets are plain Python ints
        (the BitSet type is available for ergonomic access when needed).

ERCData: one Elementary Reaction Closure.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Iterator


# ---------------------------------------------------------------------------
# BitSet
# ---------------------------------------------------------------------------

class BitSet:
    """
    Immutable bitset backed by a Python int.

    All set operations delegate to Python int bit-ops, which use
    arbitrary-precision arithmetic and map to hardware integers for
    small networks.  The wrapper exists so the representation can be
    swapped (e.g. to numpy uint64 arrays) without changing call sites.
    """
    __slots__ = ('_v',)

    def __init__(self, v: int = 0):
        self._v = int(v)

    # ---- constructors -------------------------------------------------------

    @classmethod
    def from_indices(cls, indices) -> 'BitSet':
        """Build a BitSet from an iterable of bit positions (species indices)."""
        v = 0
        for i in indices:
            v |= 1 << i
        return cls(v)

    @classmethod
    def full(cls, n: int) -> 'BitSet':
        """BitSet with bits 0..n-1 all set."""
        return cls((1 << n) - 1)

    # ---- raw access (for hot paths) ----------------------------------------

    @property
    def value(self) -> int:
        return self._v

    def __int__(self) -> int:
        return self._v

    # ---- membership --------------------------------------------------------

    def __contains__(self, i: int) -> bool:
        return bool(self._v & (1 << i))

    # ---- set operations (return new BitSet) --------------------------------

    def __or__(self, other: 'BitSet') -> 'BitSet':
        return BitSet(self._v | other._v)

    def __and__(self, other: 'BitSet') -> 'BitSet':
        return BitSet(self._v & other._v)

    def __xor__(self, other: 'BitSet') -> 'BitSet':
        return BitSet(self._v ^ other._v)

    def andnot(self, other: 'BitSet') -> 'BitSet':
        """self & ~other  (avoids unbounded ~other in Python)."""
        return BitSet(self._v & ~other._v)

    # ---- predicates --------------------------------------------------------

    def issubset(self, other: 'BitSet') -> bool:
        return (self._v & other._v) == self._v

    def issuperset(self, other: 'BitSet') -> bool:
        return (other._v & self._v) == other._v

    def isdisjoint(self, other: 'BitSet') -> bool:
        return (self._v & other._v) == 0

    def __eq__(self, other) -> bool:
        if isinstance(other, BitSet):
            return self._v == other._v
        return NotImplemented

    def __lt__(self, other: 'BitSet') -> bool:
        """Strict subset."""
        return self._v != other._v and (self._v & other._v) == self._v

    def __le__(self, other: 'BitSet') -> bool:
        return (self._v & other._v) == self._v

    # ---- collection interface ----------------------------------------------

    def __bool__(self) -> bool:
        return bool(self._v)

    def __len__(self) -> int:
        return bin(self._v).count('1')

    def __iter__(self) -> Iterator[int]:
        """Yield bit positions in ascending order."""
        v = self._v
        while v:
            lsb = v & (-v)
            yield lsb.bit_length() - 1
            v &= v - 1

    def __hash__(self) -> int:
        return hash(self._v)

    def __repr__(self) -> str:
        bits = list(self)
        return f"BitSet({bits})"


# ---------------------------------------------------------------------------
# RNData — compiled reaction network (read-only after construction)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class RNData:
    """
    Compiled, bitset-indexed reaction network.

    Species are indexed 0..n_species-1.
    Reactions are indexed 0..n_reactions-1.
    All sets are plain Python ints (bitsets); use BitSet() to wrap if needed.

    Stage-0 invariants (asserted in io_pyCOT.build_rndata):
      • E0_mask = closure_oracle(supp_raw, prod_raw, inflow_prod_mask)
      • For every r: supp_q[r] = supp_raw[r] & ~E0_mask
      • For every r: prod_q[r] = prod_raw[r] & ~E0_mask
      • If supp_q[r] == 0 then prod_q[r] == 0  (provable from E0 correctness)
      • species_to_reactions[s] = sorted list of reaction indices r with
        bit s set in supp_q[r]
    """
    # ---- species universe --------------------------------------------------
    n_species: int
    species_names: tuple[str, ...]       # index → name
    species_index: tuple[tuple[str, int], ...]  # stored as items of dict for freeze

    # ---- reactions (raw) ---------------------------------------------------
    n_reactions: int
    reaction_names: tuple[str, ...]
    supp_raw: tuple[int, ...]            # support bitmasks (pre-E0)
    prod_raw: tuple[int, ...]            # product bitmasks (pre-E0)

    # ---- E0 ---------------------------------------------------------------
    E0_mask: int                         # bitmask of E0 species

    # ---- quotiented (E0 stripped) -----------------------------------------
    supp_q: tuple[int, ...]
    prod_q: tuple[int, ...]

    # ---- inverted index (for Horn propagation) ----------------------------
    # species_to_reactions[i] = tuple of reaction indices r with i in supp_q[r]
    species_to_reactions: tuple[tuple[int, ...], ...]

    # ---- convenience -------------------------------------------------------

    def _sp_index(self) -> dict[str, int]:
        """Reconstruct the species-name→index dict (not stored frozen)."""
        return dict(self.species_index)

    def species_name(self, i: int) -> str:
        return self.species_names[i]

    def species_idx(self, name: str) -> int:
        return dict(self.species_index)[name]

    def bitset_to_names(self, mask: int) -> list[str]:
        result = []
        m = mask
        i = 0
        while m:
            if m & 1:
                result.append(self.species_names[i])
            m >>= 1
            i += 1
        return sorted(result)

    def names_to_bitset(self, names) -> int:
        idx = dict(self.species_index)
        v = 0
        for name in names:
            v |= 1 << idx[name]
        return v

    def reaction_name(self, r: int) -> str:
        return self.reaction_names[r]

    # ---- Stage-0 assertions ------------------------------------------------

    def assert_stage0_invariants(self) -> None:
        """Raise AssertionError if any Stage-0 invariant is violated."""
        idx = dict(self.species_index)

        # species_index must be a bijection with species_names
        assert len(set(self.species_names)) == self.n_species, "duplicate species names"
        assert len(idx) == self.n_species, "species_index length mismatch"
        for name, i in idx.items():
            assert self.species_names[i] == name, f"species_index/names mismatch at {name}"

        # Quotiented masks must equal raw & ~E0
        for r in range(self.n_reactions):
            assert self.supp_q[r] == (self.supp_raw[r] & ~self.E0_mask), \
                f"supp_q mismatch at reaction {r}"
            assert self.prod_q[r] == (self.prod_raw[r] & ~self.E0_mask), \
                f"prod_q mismatch at reaction {r}"
            # Correctness of E0: trivial reactions must produce nothing new
            if self.supp_q[r] == 0:
                assert self.prod_q[r] == 0, \
                    f"reaction {r}: supp_q=0 but prod_q≠0 — E0 computation error"

        # Inverted index consistency
        assert len(self.species_to_reactions) == self.n_species, "inv_idx length"
        for s, rx_list in enumerate(self.species_to_reactions):
            for r in rx_list:
                assert (self.supp_q[r] >> s) & 1, \
                    f"inv_idx[{s}] lists reaction {r} but bit {s} not in supp_q[{r}]"
        # reverse check: every (r,s) in supp_q appears in inv_idx
        for r in range(self.n_reactions):
            m = self.supp_q[r]
            while m:
                lsb = m & (-m)
                s = lsb.bit_length() - 1
                assert r in self.species_to_reactions[s], \
                    f"supp_q[{r}] has bit {s} but inv_idx[{s}] missing reaction {r}"
                m &= m - 1


# ---------------------------------------------------------------------------
# ERCData — one ERC
# ---------------------------------------------------------------------------

@dataclass
class ERCData:
    """
    One Elementary Reaction Closure.

    species_mask: the set of species in this ERC (quotiented, i.e. E0 excluded).
    reaction_indices: all r with closure(supp_q[r]) == species_mask.
    min_bases: ⊆-minimal elements of {supp_q[r] : r in reaction_indices}.
               These are the minimal ways to "ignite" this ERC.
    req_mask: supp(R_E) & ~prod(R_E) — species consumed but not produced.
    prod_mask: prod(R_E) — species produced by any active reaction in E.
    """
    erc_id: int          # position in the sorted ERC list (by species_mask popcount)
    species_mask: int    # bitset of species
    reaction_indices: list[int]
    min_bases: list[int]
    req_mask: int = 0
    prod_mask: int = 0

    # ---- convenience -------------------------------------------------------

    def is_persistent(self) -> bool:
        """True if this ERC is semi-self-maintaining: req_mask == 0."""
        return self.req_mask == 0

    def size(self) -> int:
        return bin(self.species_mask).count('1')

    def __hash__(self):
        return hash(self.species_mask)

    def __eq__(self, other):
        if isinstance(other, ERCData):
            return self.species_mask == other.species_mask
        return NotImplemented

    def __repr__(self) -> str:
        return (f"ERC(id={self.erc_id}, size={self.size()}, "
                f"persistent={self.is_persistent()}, "
                f"n_reactions={len(self.reaction_indices)}, "
                f"n_min_bases={len(self.min_bases)})")
