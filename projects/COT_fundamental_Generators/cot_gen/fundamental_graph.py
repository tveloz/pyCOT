"""
fundamental_graph.py — Unified ERC adjacency graph for EPM/ESPM traversal.

All traversal state lives in ERC-index space.  No species-level closure
computation is performed during the DFS.

Mathematical basis
------------------
For a set of ERCs S, define:

  req(S)  = (∪ req_mask[i]) & ~(∪ prod_mask[i])   for i ∈ S
  prod(S) = ∪ prod_mask[i]                          for i ∈ S
  sp(S)   = ∪ species_mask[i]                       for i ∈ S

These formulae are EXACT because:
  • req_mask[i] = supp(R_{E_i}) & ~prod(R_{E_i})
    incorporates every reaction inside E_i, including those of all sub-ERCs
    E_j ⊊ E_i (since R_{E_j} ⊆ R_{E_i} by containment in the hierarchy).
  • When new ERCs join S via a fundamental synergy (E_a, E_b) → E_k,
    E_k captures all cross-boundary reactions activated by that combination.
    Propagating these implications (Horn-style) via erc_syn_close() ensures
    no active reaction is missed.

Therefore: SSM(S) iff req(S) == 0 — no reaction scan required.

Public API
----------
FundamentalGraph(ercs, hier, syn_result, comp_result)
    Build the unified graph from Stage-1/2/3 outputs.

DFSState  — immutable traversal state: (erc_set, sp, req, prod)

FundamentalGraph.make_seed_state(i) -> DFSState
    Initial state for a single-ERC seed.

FundamentalGraph.extend_state(state, new_idx) -> DFSState
    Add one ERC (plus all synergy implications) and return the new state.

FundamentalGraph.erc_syn_close(current_set, new_indices) -> frozenset
    Horn propagation: find all ERCs implied by adding new_indices to current_set.
"""
from __future__ import annotations
from dataclasses import dataclass


# ---------------------------------------------------------------------------
# DFSState — immutable traversal state
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class DFSState:
    """
    One node in the Mode-1 DFS, expressed in ERC-index space.

    Fields
    ------
    erc_set : frozenset[int]
        Indices of ERCs currently active.  Maintained synergy-closed:
        if (E_a, E_b) → E_k is a fundamental synergy and E_a, E_b ∈ erc_set
        then E_k ∈ erc_set.

    sp : int
        Union of species_mask[i] for i in erc_set.
        The full set of species present in this state.

    req : int
        (∪ req_mask[i]) & ~(∪ prod_mask[i]) for i in erc_set.
        Species still required from outside.  Zero iff the state is SSM.

    prod : int
        ∪ prod_mask[i] for i in erc_set.
        All species produced by reactions active in this state.

    min_ext : int
        Canonical ordering gate: the next *explicit* extension (comp or syn
        partner choice) must have ERC index ≥ min_ext.  This enforces that
        every SSM {e_1 < e_2 < … < e_k} is assembled in strictly increasing
        index order, eliminating the k! redundant orderings.
        Default 0 = no constraint (used by Mode-2 / external callers).

    All fields except min_ext are maintained incrementally by extend_state();
    no reaction scan is ever performed after the graph is built.
    """
    erc_set: frozenset
    sp:      int
    req:     int
    prod:    int
    min_ext: int = 0  # canonical ordering gate (0 = unconstrained)

    @property
    def is_ssm(self) -> bool:
        """True iff this state is semi-self-maintaining (req == 0)."""
        return self.req == 0


# ---------------------------------------------------------------------------
# FundamentalGraph
# ---------------------------------------------------------------------------

class FundamentalGraph:
    """
    Unified ERC adjacency graph combining hierarchy + synergy + complementarity.

    Nodes     : ERC indices 0..n-1
    Syn edges : fundamental synergy hyperedges (i,j)→k, stored bidirectionally
    Comp edges: fundamental complementarity, indexed by species bit
    Node data : req_mask, prod_mask, species_mask, is_persistent per ERC

    Usage: build once, then call make_seed_state / extend_state for DFS.
    """

    def __init__(self, ercs, hier, syn_result, comp_result):
        """
        Build the graph from Stage outputs.

        Parameters
        ----------
        ercs        : list[ERCData]       — all ERCs (sorted by species_mask)
        hier        : HierarchyData       — ERC containment hierarchy
        syn_result  : SynergyResult       — fundamental synergies (st.i, st.j, st.k)
        comp_result : CompResult          — fundamental complementarities (fc.species, fc.prod_idx)
        """
        n = len(ercs)
        self.n = n

        # ── Per-ERC masks ─────────────────────────────────────────────────
        self.req_mask      = [e.req_mask      for e in ercs]
        self.prod_mask     = [e.prod_mask     for e in ercs]
        self.species_mask  = [e.species_mask  for e in ercs]
        self.is_persistent = [e.is_persistent() for e in ercs]

        # ── Hierarchy (reuse frozensets from HierarchyData) ───────────────
        self.ancestors   = hier.ancestors    # ancestors[i]   = frozenset of strict supersets
        self.descendants = hier.descendants  # descendants[i] = frozenset of strict subsets
        self.parents     = hier.parents      # parents[i]     = direct covers of i (Hasse edges)

        # ── Fundamental synergy edges ──────────────────────────────────────
        # syn_from[i] = [(other, target), ...]
        # Meaning: if BOTH i and other are in erc_set, then target is implied.
        # Stored symmetrically: (i,j)→k gives entries in both syn_from[i] and syn_from[j].
        syn_from: dict[int, list[tuple[int, int]]] = {}
        for st in syn_result.fundamental:
            syn_from.setdefault(st.i, []).append((st.j, st.k))
            syn_from.setdefault(st.j, []).append((st.i, st.k))
        self.syn_from = syn_from

        # ── Fundamental complementarity edges ─────────────────────────────
        # comp_by_species[species_bit] = [prod_erc_idx, ...]
        # Meaning: for species s missing from the current state, these are the
        # ⊆-minimal ERCs that produce s.
        #
        # comp_consumers_by_species[species_bit] = [cons_erc_idx, ...]
        # Meaning: for species s already produced by the current state, these
        # are the ⊆-minimal ERCs that require s externally (Def 26, consumer
        # side).  Used by Mode-2 to grow an already-SSM module outward via
        # complementarity: attaching cons_erc_idx does not close any open
        # requirement (there is none — the module is already SSM) but adds
        # productive novelty by discharging cons_erc_idx's own requirement
        # for s "for free".
        comp_by_species: dict[int, list[int]] = {}
        comp_consumers_by_species: dict[int, list[int]] = {}
        seen: set[tuple[int, int]] = set()
        seen_cons: set[tuple[int, int]] = set()
        for fc in comp_result.fundamental:
            key = (fc.species, fc.prod_idx)
            if key not in seen:
                seen.add(key)
                comp_by_species.setdefault(fc.species, []).append(fc.prod_idx)
            ckey = (fc.species, fc.cons_idx)
            if ckey not in seen_cons:
                seen_cons.add(ckey)
                comp_consumers_by_species.setdefault(fc.species, []).append(fc.cons_idx)
        self.comp_by_species = comp_by_species
        self.comp_consumers_by_species = comp_consumers_by_species

    # -----------------------------------------------------------------------
    # Synergy Horn propagation
    # -----------------------------------------------------------------------

    def erc_syn_close(self, current_sp: int, new_indices) -> tuple[frozenset, int]:
        """
        Propagate synergy implications when new_indices join a state with
        species coverage current_sp.

        Firing test (Lemma "Ignition", companion paper): a fundamental synergy
        (i, j) → k fires once species_mask[i] AND species_mask[j] are both
        ⊆ the state's species coverage — NOT once i and j are literally ERC
        *indices* present in erc_set.  These differ whenever a bigger ERC is
        in the state (e.g. via a vertical lift) that species-wise dominates
        the minimal reactant ERC i without i itself ever being an explicit
        member: species_mask[i] ⊆ species_mask[bigger] ⊆ sp(state) is enough
        to ignite i, even though i ∉ erc_set.  Testing ERC-index membership
        instead of a species-subset test under-fires in exactly that case,
        silently leaving the state's species/req/prod out of sync with its
        true (reaction-level) closure — a reaction whose support becomes
        covered would then never contribute its product.

        This is Horn propagation on the synergy hyperedge graph:
          clause: (species_mask[i] ⊆ sp) ∧ (species_mask[j] ⊆ sp) → (k ∈ S)

        Lookup keys are seeded not just from new_indices but also from their
        hierarchy descendants: a descendant d of a newly added ERC e has
        species_mask[d] ⊆ species_mask[e] ⊆ sp by construction, so d is
        ignited too and may itself be the minimal reactant half of some
        registered fundamental synergy (indexed under d, not under e).

        Parameters
        ----------
        current_sp  : int — species bitmask already covered before this addition
        new_indices : iterable[int] — ERCs being added now

        Returns
        -------
        (implied, sp)
          implied : frozenset — new_indices ∪ all transitively implied targets
          sp      : int — species coverage after adding implied's species
        """
        implied = set(new_indices)
        sp = current_sp
        for i in new_indices:
            sp |= self.species_mask[i]

        queue: list[int] = []
        for e in new_indices:
            queue.append(e)
            queue.extend(self.descendants[e])

        seen_lookup: set[int] = set()
        while queue:
            j = queue.pop()
            if j in seen_lookup:
                continue
            seen_lookup.add(j)
            for (other, target) in self.syn_from.get(j, []):
                if target in implied:
                    continue
                if (self.species_mask[other] & sp) == self.species_mask[other]:
                    implied.add(target)
                    sp |= self.species_mask[target]
                    queue.append(target)
                    queue.extend(self.descendants[target])
        return frozenset(implied), sp

    # -----------------------------------------------------------------------
    # State construction
    # -----------------------------------------------------------------------

    def make_seed_state(self, i: int) -> DFSState:
        """
        Build the initial DFSState for a single-ERC seed i.

        A single ERC cannot form a synergy pair with itself, so no synergy
        propagation is needed.  req and prod are taken directly from ERC i's
        masks, which already incorporate all sub-ERC reactions.

        min_ext is set to i+1: the canonical-ordering rule requires that every
        subsequent explicit extension adds an ERC with index strictly greater
        than i.  This ensures each SSM {e_1 < … < e_k} is built in increasing
        order and found only from seed e_1.

        Parameters
        ----------
        i : int — index of the seed ERC in the ercs list

        Returns
        -------
        DFSState with erc_set = {i}, min_ext = i+1
        """
        return DFSState(
            erc_set=frozenset([i]),
            sp=self.species_mask[i],
            req=self.req_mask[i],
            prod=self.prod_mask[i],
            min_ext=i + 1,
        )

    def extend_state(self, state: DFSState, new_idx: int) -> DFSState:
        """
        Add ERC new_idx (and all synergy implications) to state.

        Algorithm
        ---------
        1. Run erc_syn_close to find all ERCs implied by adding new_idx
           given what is already in state.erc_set.
        2. Accumulate req_mask, prod_mask, species_mask for all implied ERCs.
        3. Update req incrementally:
             new_req = (state.req | ∪ req_mask[implied]) & ~(state.prod | ∪ prod_mask[implied])
           This is correct because species internally recycled within existing
           ERCs (in state.prod but not state.req) cannot become newly needed.

        Parameters
        ----------
        state   : DFSState — current state (new_idx must NOT be in state.erc_set)
        new_idx : int — ERC index to add

        Returns
        -------
        New DFSState with all four fields updated.
        """
        implied, _ = self.erc_syn_close(state.sp, [new_idx])

        add_req  = 0
        add_prod = 0
        add_sp   = 0
        for i in implied:
            add_req  |= self.req_mask[i]
            add_prod |= self.prod_mask[i]
            add_sp   |= self.species_mask[i]

        new_prod = state.prod | add_prod
        new_req  = (state.req | add_req) & ~new_prod

        return DFSState(
            erc_set=state.erc_set | implied,
            sp=state.sp | add_sp,
            req=new_req,
            prod=new_prod,
            min_ext=new_idx + 1,  # next explicit choice must exceed this extension
        )
