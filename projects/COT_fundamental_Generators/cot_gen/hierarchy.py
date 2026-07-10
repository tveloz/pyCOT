"""
hierarchy.py — Bitset-based ERC containment order and Hasse diagram.

Public API
----------
build_hierarchy(ercs) -> HierarchyData
    Compute the Hasse diagram of the ERC containment order.

    E1 ≤ E2  ⟺  E1.species_mask ⊆ E2.species_mask
    (E1 is "contained in" / "below" E2)

    The Hasse diagram stores only the *direct* cover relations:
      E2 directly covers E1  ⟺  E1 < E2 and ∄ E3: E1 < E3 < E2

Returns HierarchyData with:
  • parents[i]   : list of ERC indices j where ercs[j] directly covers ercs[i]
  • children[i]  : list of ERC indices j directly covered by ercs[i]
  • ancestors[i] : frozenset of all ERC indices strictly above ercs[i]
  • descendants[i]: frozenset of all ERC indices strictly below ercs[i]

Indexing: the index here is the position in the `ercs` list, NOT erc_id.

Complexity: O(|E|²) for Hasse edges, O(|E|²) for ancestor/descendant closure.
This is acceptable because |E| ≪ 2^|M| and the hierarchy is small.

See also: cot_gen/synergy.py which uses HierarchyData for fast pair filtering.
"""
from __future__ import annotations

from dataclasses import dataclass, field


# ---------------------------------------------------------------------------
# HierarchyData
# ---------------------------------------------------------------------------

@dataclass
class HierarchyData:
    """
    Containment Hasse diagram for a list of ERCData objects.

    All indices refer to positions in the `ercs` list passed to build_hierarchy().
    i < j in the partial order  ⟺  ercs[i].species_mask ⊊ ercs[j].species_mask.
    """
    n: int                              # number of ERCs
    parents:    list[list[int]]         # parents[i]  = direct covers of i
    children:   list[list[int]]         # children[i] = ERCs directly covered by i
    ancestors:  list[frozenset[int]]    # all j with ercs[j] ⊋ ercs[i]
    descendants: list[frozenset[int]]   # all j with ercs[j] ⊊ ercs[i]

    def is_ancestor(self, i: int, j: int) -> bool:
        """Return True if ercs[j] strictly contains ercs[i] (j is above i)."""
        return j in self.ancestors[i]

    def is_comparable(self, i: int, j: int) -> bool:
        """Return True if i and j are in containment relation (either direction)."""
        return j in self.ancestors[i] or i in self.ancestors[j]

    def can_interact(self, i: int, j: int) -> bool:
        """Return True if i and j are incomparable (neither contains the other)."""
        return not self.is_comparable(i, j)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def build_hierarchy(ercs) -> HierarchyData:
    """
    Build the ERC containment Hasse diagram from a list of ERCData.

    Algorithm
    ---------
    1. For each pair (i, j): check containment via (mask_i & mask_j) == mask_i.
       Record all strict-subset pairs.
    2. Reduce to Hasse diagram: (i, j) is a direct edge iff no k with i < k < j.
    3. Compute transitive closure for ancestors/descendants.

    Parameters
    ----------
    ercs : list[ERCData]  (from compute_ercs)

    Returns
    -------
    HierarchyData
    """
    n = len(ercs)
    masks = [e.species_mask for e in ercs]

    # --- Step 1: all pairs in strict containment order ----------------------
    # below[i] = set of j with mask_i ⊊ mask_j (j is strictly above i)
    below: list[set[int]] = [set() for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            mi, mj = masks[i], masks[j]
            if mi == mj:
                continue
            if (mi & mj) == mi:   # mi ⊊ mj  →  j is strictly above i
                below[i].add(j)

    # --- Step 2: Hasse reduction (remove transitive edges) ------------------
    # Direct cover: (i, j) is direct iff no k with j in below[k] and k in below[i]
    parents:  list[list[int]] = [[] for _ in range(n)]
    children: list[list[int]] = [[] for _ in range(n)]

    for i in range(n):
        for j in below[i]:
            # j directly covers i iff no k in below[i] with i < k and k < j
            # i.e., j NOT in below[k] for any k that is between i and j
            is_direct = True
            for k in below[i]:
                if k != j and j in below[k]:
                    # k is above i and j is above k → j is not a direct cover of i
                    is_direct = False
                    break
            if is_direct:
                parents[i].append(j)
                children[j].append(i)

    # --- Step 3: transitive closure (ancestors and descendants) -------------
    # ancestors[i] = all j with mask_i ⊊ mask_j (entire upset of i, excl. i)
    ancestors:   list[frozenset[int]] = [frozenset(below[i]) for i in range(n)]

    # descendants[i] = all j with mask_j ⊊ mask_i (entire downset of i, excl. i)
    above: list[set[int]] = [set() for _ in range(n)]
    for i in range(n):
        for j in below[i]:
            above[j].add(i)
    descendants: list[frozenset[int]] = [frozenset(above[i]) for i in range(n)]

    return HierarchyData(
        n=n,
        parents=parents,
        children=children,
        ancestors=ancestors,
        descendants=descendants,
    )
