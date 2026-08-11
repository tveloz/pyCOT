"""
explorer.py — Interactive manual exploration of the fundamental graph
(ERC hierarchy + fundamental synergy + fundamental complementarity), and of
how a hand-picked generator builds up towards a semi-organization (SSM).

This is a research/understanding tool, not part of the validated EPM/ESPM
search pipeline (cot_gen/epm.py).  It reuses the same precomputed
structures (ERCData, HierarchyData, SynergyResult, CompResult) and the
same state-tracking machinery (FundamentalGraph.extend_state) that the
search uses, so "what the generator has reached so far" is always computed
the validated way — but *which* ERCs get added, in *what* order, is
entirely up to you: nothing here enforces canonical ordering, minimality,
or fundamentality.  It's meant to let you build (and un-build) a partial
module by hand and watch req/prod/sp evolve.

Typical use — from the companion script's embedded console, or directly
from a Jupyter notebook / IPython session:

    from cot_gen.explorer import GraphExplorer
    gx = GraphExplorer.from_network("e_coli_core")

    gx.describe(0)                  # inspect one ERC
    gx.neighbors(0, radius=2)       # what's near it (containment/syn/comp)
    gx.local_view(0, radius=2)      # ... and draw that neighborhood

    gx.add(0)                       # start a generator
    gx.candidates()                 # what could be added next, and why
    gx.add(7)
    gx.summary()                    # species reached, req, prod, is_ssm
    gx.degree_stats()               # connections per node + distribution
    gx.history()                    # how req/sp evolved as ERCs were added

    gx.remove(7)                    # roll back
    gx.checkpoint("branch_a")
    gx.add(12)                      # try an alternative extension
    gx.restore("branch_a")          # ... or go back and try another

    gx.generator_view()             # exactly the generator's own subgraph + candidates
    gx.context_view()               # whole-network hierarchy, generator lit up
    gx.local_view()                 # generator + its radius-hop neighborhood

Public API
----------
GraphExplorer.from_network(name_or_path, *, data_root=None) -> GraphExplorer
GraphExplorer(rn, ercs, hier, syn, comp)

  Generator manipulation
    .add(idx)  .remove(idx)  .reset()
    .undo(n=1)
    .checkpoint(label)  .restore(label)  .checkpoints()
    .generator            -- current ordered list[int] of ERC indices
    .state                -- current DFSState (sp/req/prod), or None if empty

  Inspection (prints to console; also returns data for programmatic use)
    .describe(idx)
    .neighbors(idx, radius=1, relations=('containment','synergy','complementarity')) -> set[int]
    .candidates() -> dict[str, list[tuple[int, str]]]
    .summary()
    .degree_stats()
    .synergy_level_distribution() -> dict[int, int]
    .history()

  Visualization (writes an HTML file and opens it in the browser)
    .generator_view(show_candidates=True, relations=..., synergy_layout="below", filename=None)
    .local_view(idx=None, radius=2, relations=..., synergy_layout="below", filename=None)
    .context_view(relations=..., max_nodes=150, synergy_layout="below", filename=None)

  synergy_layout="interleaved" (any view) puts each fundamental synergy's
  diamond in the inter-level gap just above its lower reactant's own
  containment level, instead of one shared band below everything -- so
  the band with the most diamonds shows which hierarchy depth generates
  the most fundamental synergies.
"""
from __future__ import annotations

import os
import webbrowser
from dataclasses import dataclass, field

from pyCOT.analysis.organizations.erc import compute_ercs
from pyCOT.analysis.organizations.hierarchy import build_hierarchy
from pyCOT.analysis.organizations.synergy import compute_synergies_basis_first
from pyCOT.analysis.organizations.complementarity import compute_complementarities
from pyCOT.analysis.organizations.fundamental_graph import FundamentalGraph, DFSState

RELATIONS = ("containment", "synergy", "complementarity")

_REL_COLOR = {
    "containment":    "#7f8c8d",   # gray
    "synergy":        "#e67e22",   # orange
    "complementarity": "#2980b9",  # blue
}


def _bits(mask: int):
    m = mask
    while m:
        lsb = m & (-m)
        yield lsb.bit_length() - 1
        m &= m - 1


# ---------------------------------------------------------------------------
# Loading a network
# ---------------------------------------------------------------------------

def _default_data_root() -> str:
    _here = os.path.dirname(os.path.abspath(__file__))
    return os.path.normpath(os.path.join(_here, "..", "..", "..", "data", "biomodels"))


def _discover_networks(root: str) -> dict[str, str]:
    seen: dict[str, tuple[str, int]] = {}
    for dirpath, _dirs, files in os.walk(root):
        for fname in files:
            if not fname.endswith(".txt"):
                continue
            full = os.path.join(dirpath, fname)
            name = os.path.splitext(fname)[0]
            if name.startswith("bigg_"):
                name = name[5:]
            if name not in seen or len(full) < seen[name][1]:
                seen[name] = (full, len(full))
    return {name: path for name, (path, _) in seen.items()}


def load_network(name_or_path: str, *, data_root: str | None = None):
    """
    Resolve a network name or file path to (RNData, list[ERCData], net_id).

    Mirrors run_network.py's resolution logic: exact file path, exact
    catalogue name, or unique case-insensitive prefix match.
    """
    from pyCOT.io.functions import read_txt
    from pyCOT.analysis.organizations.io_pyCOT import build_rndata

    if os.path.isfile(name_or_path):
        path = os.path.abspath(name_or_path)
        net_id = os.path.splitext(os.path.basename(path))[0]
    else:
        root = data_root or _default_data_root()
        catalogue = _discover_networks(root)
        if name_or_path in catalogue:
            net_id, path = name_or_path, catalogue[name_or_path]
        else:
            lower = name_or_path.lower()
            matches = [(k, v) for k, v in catalogue.items() if k.lower().startswith(lower)]
            if len(matches) == 1:
                net_id, path = matches[0]
            elif len(matches) > 1:
                raise ValueError(
                    f"Ambiguous network name '{name_or_path}'. Matches: "
                    f"{[k for k, _ in matches[:10]]}"
                )
            else:
                raise FileNotFoundError(
                    f"Network '{name_or_path}' not found under {root}"
                )

    rn_pycot = read_txt(path)
    rn = build_rndata(rn_pycot, network_id=net_id)
    ercs = compute_ercs(rn)
    return rn, ercs, net_id


# ---------------------------------------------------------------------------
# GraphExplorer
# ---------------------------------------------------------------------------

@dataclass
class GraphExplorer:
    rn: object
    ercs: list
    hier: object
    syn: object
    comp: object

    generator: list = field(default_factory=list)   # ordered ERC indices
    state: DFSState | None = field(default=None, repr=False)

    _undo_stack: list = field(default_factory=list, repr=False)
    _checkpoints: dict = field(default_factory=dict, repr=False)

    def __post_init__(self):
        self._g = FundamentalGraph(self.ercs, self.hier, self.syn, self.comp)

        # Display-only aggregation: a single reactant pair (i,j) can have
        # several simultaneous fundamental synergy targets (Emax(E,E') is an
        # antichain in general, not a singleton).  FundamentalGraph.syn_from
        # stays unaggregated (one entry per (partner,target) pair) because
        # aggregating it measurably slowed the search on genome-scale
        # networks (see fundamental_graph.py); for display purposes here,
        # where call volume is tiny by comparison, aggregating by pair is
        # unambiguously nicer to read and costs nothing noticeable.
        self._syn_by_pair: dict[tuple[int, int], set[int]] = {}
        self._syn_partners: dict[int, set[int]] = {}
        for st in self.syn.fundamental:
            key = (min(st.i, st.j), max(st.i, st.j))
            self._syn_by_pair.setdefault(key, set()).add(st.k)
            self._syn_partners.setdefault(st.i, set()).add(st.j)
            self._syn_partners.setdefault(st.j, set()).add(st.i)

        self._containment_level = self._compute_containment_levels()

    def _compute_containment_levels(self) -> list[int]:
        """
        Level 0 = ⊆-minimal ERCs (no children in the containment hierarchy);
        level increases by 1 per containment step upward, so a bigger,
        containing ERC always sits strictly above every ERC it contains.
        Computed purely from hier.children — never from synergy or
        complementarity — so plotting stays a legible tree/DAG regardless
        of which other relations get drawn on top of it (see _plot).
        """
        n = len(self.ercs)
        level = [0] * n
        # Process nodes in increasing order of descendant-set size: a child
        # always has a strictly smaller descendant set than its parent, so
        # this guarantees every child is leveled before its parent.
        order = sorted(range(n), key=lambda i: len(self.hier.descendants[i]))
        for i in order:
            children = self.hier.children[i]
            if children:
                level[i] = 1 + max(level[c] for c in children)
        return level

    # -----------------------------------------------------------------------
    # Construction
    # -----------------------------------------------------------------------

    @classmethod
    def from_network(cls, name_or_path: str, *, data_root: str | None = None) -> "GraphExplorer":
        rn, ercs, net_id = load_network(name_or_path, data_root=data_root)
        hier = build_hierarchy(ercs)
        syn = compute_synergies_basis_first(ercs, hier)
        comp = compute_complementarities(ercs, hier, syn_result=syn)
        print(f"Loaded {net_id}: "
              f"{rn.n_species} species, {rn.n_reactions} reactions, {len(ercs)} ERCs "
              f"({sum(1 for e in ercs if e.is_persistent())} persistent)")
        print(f"Fundamental synergies: {len(syn.fundamental)}   "
              f"Fundamental complementarities: {len(comp.fundamental)}")
        return cls(rn=rn, ercs=ercs, hier=hier, syn=syn, comp=comp)

    # -----------------------------------------------------------------------
    # Naming helpers
    # -----------------------------------------------------------------------

    def erc_label(self, idx: int) -> str:
        return f"E{idx}"

    def erc_species_names(self, idx: int) -> list[str]:
        return self.rn.bitset_to_names(self.ercs[idx].species_mask)

    def species_names(self, mask: int) -> list[str]:
        return self.rn.bitset_to_names(mask)

    # -----------------------------------------------------------------------
    # Generator manipulation
    # -----------------------------------------------------------------------

    def _snapshot(self):
        self._undo_stack.append(list(self.generator))

    def _recompute_state(self):
        if not self.generator:
            self.state = None
            return
        state = self._g.make_seed_state(self.generator[0])
        for idx in self.generator[1:]:
            if idx not in state.erc_set:
                state = self._g.extend_state(state, idx)
        self.state = state

    def add(self, idx: int) -> None:
        """Add ERC idx to the generator (no-op with a warning if already
        present, directly or as a synergy implication of what's there)."""
        self._require_valid_idx(idx)
        if self.state is not None and idx in self.state.erc_set:
            print(f"  {self.erc_label(idx)} is already reached by the current "
                  f"generator (directly or via synergy closure) — nothing to add.")
            return
        self._snapshot()
        self.generator.append(idx)
        self._recompute_state()
        self._print_add_summary(idx)

    def remove(self, idx: int) -> None:
        """Remove ERC idx from the generator and recompute state from the
        remaining ERCs (this is the "roll back" / "eliminate ERC" move)."""
        if idx not in self.generator:
            print(f"  {self.erc_label(idx)} is not an explicit member of the "
                  f"generator (it may only be present via synergy closure).")
            return
        self._snapshot()
        self.generator.remove(idx)
        self._recompute_state()
        print(f"  Removed {self.erc_label(idx)}. Generator now: "
              f"{[self.erc_label(i) for i in self.generator]}")

    def reset(self) -> None:
        self._snapshot()
        self.generator = []
        self.state = None
        print("  Generator cleared.")

    def undo(self, n: int = 1) -> None:
        for _ in range(n):
            if not self._undo_stack:
                print("  Nothing to undo.")
                return
            self.generator = self._undo_stack.pop()
            self._recompute_state()
        print(f"  Generator now: {[self.erc_label(i) for i in self.generator]}")

    def checkpoint(self, label: str) -> None:
        """Save the current generator under a name, to branch from later."""
        self._checkpoints[label] = list(self.generator)
        print(f"  Checkpoint '{label}' saved: {[self.erc_label(i) for i in self.generator]}")

    def restore(self, label: str) -> None:
        if label not in self._checkpoints:
            print(f"  No checkpoint named '{label}'. Known: {list(self._checkpoints)}")
            return
        self._snapshot()
        self.generator = list(self._checkpoints[label])
        self._recompute_state()
        print(f"  Restored '{label}': {[self.erc_label(i) for i in self.generator]}")

    def checkpoints(self) -> None:
        if not self._checkpoints:
            print("  No checkpoints saved yet. Use gx.checkpoint('label').")
            return
        for label, gen in self._checkpoints.items():
            print(f"  {label:<20} {[self.erc_label(i) for i in gen]}")

    def _require_valid_idx(self, idx: int) -> None:
        if not (0 <= idx < len(self.ercs)):
            raise IndexError(f"ERC index {idx} out of range 0..{len(self.ercs)-1}")

    def _print_add_summary(self, idx: int) -> None:
        s = self.state
        req_names = self.species_names(s.req)
        print(f"  Added {self.erc_label(idx)}. Generator now: "
              f"{[self.erc_label(i) for i in self.generator]}"
              f"  (erc_set closure: {sorted(self.erc_label(i) for i in s.erc_set)})")
        print(f"    species reached: {bin(s.sp).count('1')}  "
              f"| req: {req_names if req_names else '(none — SSM!)'}")

    # -----------------------------------------------------------------------
    # Inspection
    # -----------------------------------------------------------------------

    def describe(self, idx: int) -> None:
        """Print one ERC's species/req/prod and its direct relations."""
        self._require_valid_idx(idx)
        e = self.ercs[idx]
        print(f"{self.erc_label(idx)}  (erc_id={e.erc_id}, size={e.size()}, "
              f"persistent={e.is_persistent()})")
        print(f"  species : {self.erc_species_names(idx)}")
        print(f"  req     : {self.species_names(e.req_mask) or '(none)'}")
        print(f"  prod    : {self.species_names(e.prod_mask)}")
        print(f"  reactions: {[self.rn.reaction_name(r) for r in e.reaction_indices]}")

        parents = self.hier.parents[idx]
        children = self.hier.children[idx]
        print(f"  containment: parents={[self.erc_label(i) for i in parents]}  "
              f"children={[self.erc_label(i) for i in children]}")

        partners = sorted(self._syn_partners.get(idx, ()))
        n_targets_as_reactant = sum(
            1 for st in self.syn.fundamental if st.i == idx or st.j == idx
        )
        n_produced_by = sum(1 for st in self.syn.fundamental if st.k == idx)
        pair_view = [
            (self.erc_label(j), [self.erc_label(k) for k in sorted(
                self._syn_by_pair[(min(idx, j), max(idx, j))])])
            for j in partners
        ]
        print(f"  synergy: {len(partners)} partner(s), "
              f"{n_targets_as_reactant} fundamental target(s) as reactant -> {pair_view}")
        print(f"  synergy: produced by {n_produced_by} fundamental synergy(-ies) as target")

        comp_out, comp_in = self._comp_edges_for(idx)
        print(f"  complementarity (as producer): "
              f"{[(self.erc_label(c), self.species_names(1<<s)[0]) for c, s in comp_out]}")
        print(f"  complementarity (as consumer): "
              f"{[(self.erc_label(p), self.species_names(1<<s)[0]) for p, s in comp_in]}")

    def _comp_edges_for(self, idx: int):
        """(as-producer list of (consumer_idx, species_bit), as-consumer list of (producer_idx, species_bit))"""
        out_edges, in_edges = [], []
        for fc in self.comp.fundamental:
            if fc.prod_idx == idx:
                out_edges.append((fc.cons_idx, fc.species))
            if fc.cons_idx == idx:
                in_edges.append((fc.prod_idx, fc.species))
        return out_edges, in_edges

    def neighbors(
        self,
        idx: int,
        radius: int = 1,
        relations: tuple = RELATIONS,
    ) -> set[int]:
        """
        ERCs reachable from `idx` within `radius` hops along the selected
        relation types (any mix of 'containment', 'synergy',
        'complementarity'). Does not modify the generator.
        """
        self._require_valid_idx(idx)
        frontier = {idx}
        seen = {idx}
        for _ in range(radius):
            nxt = set()
            for i in frontier:
                if "containment" in relations:
                    nxt |= set(self.hier.parents[i]) | set(self.hier.children[i])
                if "synergy" in relations:
                    for j, k in self._g.syn_from.get(i, []):
                        nxt.add(j)
                        nxt.add(k)
                if "complementarity" in relations:
                    out_e, in_e = self._comp_edges_for(i)
                    nxt |= {c for c, _ in out_e} | {p for p, _ in in_e}
            nxt -= seen
            if not nxt:
                break
            seen |= nxt
            frontier = nxt
        return seen

    def _compute_candidates(self) -> dict:
        """
        Silent version of candidates(): same return value, no printing.
        Used internally by candidates() itself and by the plotting methods
        (which need the candidate set for highlighting without triggering a
        console dump every time a view is drawn).
        """
        result: dict[str, list[tuple[int, str]]] = {
            "producer": [], "consumer": [], "synergy_partner": [], "vertical_lift": [],
        }
        if self.state is None:
            return result
        erc_set = self.state.erc_set

        for s_bit in _bits(self.state.req):
            sp_name = self.species_names(1 << s_bit)[0]
            for prod_idx in self._g.comp_by_species.get(s_bit, []):
                if prod_idx not in erc_set:
                    result["producer"].append((prod_idx, f"supplies required '{sp_name}'"))
        for s_bit in _bits(self.state.prod):
            sp_name = self.species_names(1 << s_bit)[0]
            for cons_idx in self._g.comp_consumers_by_species.get(s_bit, []):
                if cons_idx not in erc_set:
                    result["consumer"].append((cons_idx, f"would consume already-produced '{sp_name}'"))
        for i in erc_set:
            for j in self._syn_partners.get(i, ()):
                if j not in erc_set:
                    targets = self._syn_by_pair[(min(i, j), max(i, j))]
                    target_lbls = ", ".join(self.erc_label(k) for k in sorted(targets))
                    result["synergy_partner"].append(
                        (j, f"synergy with {self.erc_label(i)} -> {target_lbls}"))
            for a in self.hier.parents[i]:
                if a not in erc_set:
                    result["vertical_lift"].append((a, f"ancestor of {self.erc_label(i)}"))
        return result

    def candidates(self) -> dict:
        """
        ERCs connected to the CURRENT generator that could be added next,
        grouped by the move that connects them (mirrors the search
        algorithm's own extension moves — see cot_gen/epm.py):

          'producer'        : minimal producer of a still-required species
          'consumer'        : minimal consumer of an already-produced species
          'synergy_partner' : fundamental-synergy partner of a member
          'vertical_lift'   : direct hierarchy ancestor of a member

        Prints a readable table and returns the same data as a dict.
        """
        if self.state is None:
            print("  Generator is empty — add a first ERC, e.g. gx.add(0).")
            return self._compute_candidates()

        result = self._compute_candidates()
        print(f"Candidates for extending {[self.erc_label(i) for i in self.generator]}:")
        for move, items in result.items():
            uniq = sorted(set(items), key=lambda t: t[0])
            if not uniq:
                continue
            print(f"  [{move}]")
            for idx, why in uniq:
                print(f"    {self.erc_label(idx):<6} {why}")
        if not any(result.values()):
            print("  (none — this generator is a dead end, or already SSM with no "
                  "further fundamental relation to explore)")
        return result

    def summary(self) -> None:
        """Print the full state of the current generator."""
        if self.state is None:
            print("Generator is empty.")
            return
        s = self.state
        sp_names = self.species_names(s.sp)
        req_names = self.species_names(s.req)
        prod_names = self.species_names(s.prod)
        rxns = sorted({r for i in s.erc_set for r in self.ercs[i].reaction_indices})

        print(f"Generator: {[self.erc_label(i) for i in self.generator]}")
        print(f"Closure (erc_set, incl. synergy implications): "
              f"{sorted(self.erc_label(i) for i in s.erc_set)}")
        print(f"Species reached ({len(sp_names)}): {sp_names}")
        print(f"Reactions reached ({len(rxns)}): {[self.rn.reaction_name(r) for r in rxns]}")
        print(f"Produced ({len(prod_names)}): {prod_names}")
        print(f"Required — still missing ({len(req_names)}): {req_names or '(none)'}")
        print(f"Semi-self-maintaining (SSM): {s.is_ssm}"
              f"{'  <-- this generator is already a persistent module!' if s.is_ssm else ''}")

    def degree_stats(self) -> None:
        """
        Per-node connection counts (containment/synergy/complementarity) for
        every ERC currently in the generator, plus the degree distribution
        across them.

        Synergy is reported two ways, since they answer different
        questions: 'syn_reactant' is how many DISTINCT partner ERCs this
        one can combine with as a reactant (syn_from is already aggregated
        by partner, so a partner that simultaneously triggers several
        fundamental targets still counts once); 'syn_produced_by' is how
        many different fundamental synergies produce this ERC as their
        target — i.e. how many distinct reactant pairs can generate it.
        """
        if not self.generator:
            print("Generator is empty.")
            return
        produced_by: dict[int, int] = {}
        for st in self.syn.fundamental:
            produced_by[st.k] = produced_by.get(st.k, 0) + 1

        rows = []
        for i in self.generator:
            n_cont = len(self.hier.parents[i]) + len(self.hier.children[i])
            n_syn_partner = len(self._syn_partners.get(i, ()))
            n_syn_produced_by = produced_by.get(i, 0)
            out_e, in_e = self._comp_edges_for(i)
            n_comp = len(out_e) + len(in_e)
            total = n_cont + n_syn_partner + n_syn_produced_by + n_comp
            rows.append((i, n_cont, n_syn_partner, n_syn_produced_by, n_comp, total))

        print(f"{'ERC':<6} {'contain':>8} {'syn_react':>10} {'syn_prodby':>11} "
              f"{'comp':>8} {'total':>8}")
        for i, c, sr, spb, cp, t in rows:
            print(f"{self.erc_label(i):<6} {c:>8} {sr:>10} {spb:>11} {cp:>8} {t:>8}")

        totals = [t for *_, t in rows]
        if totals:
            import statistics as _st
            print(f"\nTotal-degree distribution over generator "
                  f"(n={len(totals)}): min={min(totals)} "
                  f"median={_st.median(totals):.1f} max={max(totals)} "
                  f"mean={_st.mean(totals):.2f}")
            buckets: dict[int, int] = {}
            for t in totals:
                buckets[t] = buckets.get(t, 0) + 1
            for deg in sorted(buckets):
                print(f"    degree {deg:>3}: {'#' * buckets[deg]} ({buckets[deg]})")

    def synergy_level_distribution(self) -> dict[int, int]:
        """
        Count fundamental synergies (network-wide, not just the current
        generator) by the containment level of their LOWER reactant --
        i.e. the same band each one is placed in by
        local_view(synergy_layout="interleaved") etc.  Prints a bar chart
        and returns {level: count}; the level(s) with the largest count
        are where fundamental synergies fire most often, hierarchy-depth
        -wise.
        """
        buckets: dict[int, int] = {}
        seen_pairs: set[tuple[int, int]] = set()
        for st in self.syn.fundamental:
            pair_key = (min(st.i, st.j), max(st.i, st.j))
            if pair_key in seen_pairs:
                continue
            seen_pairs.add(pair_key)
            lvl = min(self._containment_level[st.i], self._containment_level[st.j])
            buckets[lvl] = buckets.get(lvl, 0) + 1

        if not buckets:
            print("  No fundamental synergies in this network.")
            return buckets
        print(f"Fundamental synergy pairs by lower-reactant containment level "
              f"(n={sum(buckets.values())} pairs):")
        for lvl in sorted(buckets):
            print(f"  level {lvl}-{lvl+1}: {'#' * buckets[lvl]} ({buckets[lvl]})")
        return buckets

    def history(self) -> None:
        """
        Replay the generator's construction in insertion order, showing how
        req/sp evolved and exactly when (if ever) it became SSM -- i.e.
        when the growing set first formed a semi-organization.
        """
        if not self.generator:
            print("Generator is empty.")
            return
        state = self._g.make_seed_state(self.generator[0])
        print(f"{'step':<5} {'add':<6} {'#species':>9} {'#req':>6} {'SSM?':>6}")
        print(f"{0:<5} {self.erc_label(self.generator[0]):<6} "
              f"{bin(state.sp).count('1'):>9} {bin(state.req).count('1'):>6} "
              f"{'YES' if state.is_ssm else 'no':>6}"
              f"{'  <- first SSM (semi-org formed)' if state.is_ssm else ''}")
        first_ssm_step = 0 if state.is_ssm else None
        for step, idx in enumerate(self.generator[1:], start=1):
            if idx in state.erc_set:
                print(f"{step:<5} {self.erc_label(idx):<6}  (already implied — no change)")
                continue
            state = self._g.extend_state(state, idx)
            became_ssm = state.is_ssm and first_ssm_step is None
            if became_ssm:
                first_ssm_step = step
            print(f"{step:<5} {self.erc_label(idx):<6} "
                  f"{bin(state.sp).count('1'):>9} {bin(state.req).count('1'):>6} "
                  f"{'YES' if state.is_ssm else 'no':>6}"
                  f"{'  <- first SSM (semi-org formed)' if became_ssm else ''}")
        if first_ssm_step is None:
            print("\nThis generator never reached SSM (still an open, non-persistent partial module).")

    # -----------------------------------------------------------------------
    # Visualization
    # -----------------------------------------------------------------------

    def generator_view(
        self,
        show_candidates: bool = True,
        relations: tuple = RELATIONS,
        synergy_layout: str = "below",
        filename: str | None = None,
    ) -> str:
        """
        Draw exactly what the generator's closure actually contains, plus
        (optionally) the ERCs that are one genuine fundamental relation
        away from it. Two-step construction, both grounded directly in
        what's definitely already part of the closure -- never in
        second-order relations of a relation:

          1. "Contained" = self.state.erc_set (the generator's true
             synergy-closure, exactly as computed by extend_state /
             erc_syn_close -- so any fundamental synergy whose reactants
             are already covered has already fired and its target is
             already in here) UNION each of those ERCs' hierarchy
             descendants (a descendant's species are a subset of its
             containing ERC's, so it is trivially part of the closure too,
             even though erc_syn_close never lists it explicitly -- it only
             uses descendants as lookup keys internally).
          2. Candidates (only if show_candidates=True) = every ERC that has
             a *direct* fundamental relation to some member of "contained"
             -- a fundamental complementarity in EITHER direction (no
             producer/consumer distinction), or a fundamental synergy
             partnership. Nothing is shown merely because it relates to
             another candidate; every node's hover text says exactly which
             contained ERC it connects to and how.

        Vertical-lift (hierarchy ancestor) candidates are deliberately not
        included here -- this view is about the closure and its direct
        complementarity/synergy relations only, not every possible
        extension move (see .candidates() / local_view for that).

        synergy_layout: see _plot -- "below" (default) or "interleaved".
        """
        if not self.generator:
            print("  Generator is empty — add a first ERC, e.g. gx.add(0).")
            return ""

        closure = set(self.state.erc_set) if self.state is not None else set(self.generator)
        reasons: dict[int, str] = {i: "explicit generator member" for i in self.generator}
        for m in closure - set(self.generator):
            reasons.setdefault(m, "fundamental-synergy closure member (reachable once its reactants were covered)")

        contained = set(closure)
        for m in list(closure):
            for d in self.hier.descendants[m]:
                if d not in contained:
                    contained.add(d)
                    reasons.setdefault(d, f"descendant of {self.erc_label(m)} (already inside its closure)")

        node_set = set(contained)
        if show_candidates:
            for fc in self.comp.fundamental:
                sp_name = self.species_names(1 << fc.species)[0]
                if fc.prod_idx in contained and fc.cons_idx not in contained:
                    node_set.add(fc.cons_idx)
                    reasons.setdefault(
                        fc.cons_idx,
                        f"requires '{sp_name}', which {self.erc_label(fc.prod_idx)} (in the closure) produces")
                elif fc.cons_idx in contained and fc.prod_idx not in contained:
                    node_set.add(fc.prod_idx)
                    reasons.setdefault(
                        fc.prod_idx,
                        f"produces '{sp_name}', which {self.erc_label(fc.cons_idx)} (in the closure) requires")
            for i in contained:
                for j in self._syn_partners.get(i, ()):
                    if j not in contained:
                        node_set.add(j)
                        targets = self._syn_by_pair[(min(i, j), max(i, j))]
                        target_lbls = ", ".join(self.erc_label(k) for k in sorted(targets))
                        reasons.setdefault(
                            j, f"fundamental synergy with {self.erc_label(i)} (in the closure) -> {target_lbls}")
                        # The diamond needs its target ERC(s) present too, or
                        # _plot has nothing to draw the third edge to and
                        # silently skips the whole diamond, leaving the
                        # candidate with no visible connection at all.
                        for k in targets:
                            if k not in node_set:
                                node_set.add(k)
                                reasons.setdefault(
                                    k, f"would become reachable via the fundamental synergy "
                                       f"{self.erc_label(i)}+{self.erc_label(j)}")

        out = filename or "generator_view.html"
        return self._plot(sorted(node_set), relations, out, synergy_layout=synergy_layout,
                           anchor_set=contained, reasons=reasons)

    def local_view(
        self,
        idx: int | None = None,
        radius: int = 2,
        relations: tuple = RELATIONS,
        synergy_layout: str = "below",
        filename: str | None = None,
    ) -> str:
        """
        Draw ONLY the neighborhood of `idx` (default: the whole current
        generator) within `radius` hops along `relations` -- small and
        focused, safe on any network size. Opens in the browser.

        synergy_layout: see _plot -- "below" (default) or "interleaved".
        """
        centers = {idx} if idx is not None else set(self.generator)
        if not centers:
            print("  Nothing to center on: pass an ERC index, or gx.add(...) first.")
            return ""
        node_set: set[int] = set(centers)
        for c in centers:
            node_set |= self.neighbors(c, radius=radius, relations=relations)
        out = filename or f"local_view_{'_'.join(map(str, sorted(centers)))}.html"
        return self._plot(sorted(node_set), relations, out, synergy_layout=synergy_layout)

    def context_view(
        self,
        relations: tuple = RELATIONS,
        max_nodes: int = 150,
        synergy_layout: str = "below",
        filename: str | None = None,
    ) -> str:
        """
        Draw the WHOLE network's ERC graph with the current generator and
        its candidates highlighted, so you can see where it sits in
        context. Guarded by `max_nodes` since the full hierarchy can be
        very large; raise the cap explicitly if you really want everything.

        synergy_layout: see _plot -- "below" (default) or "interleaved".
        """
        n = len(self.ercs)
        if n > max_nodes:
            print(f"  Network has {n} ERCs > max_nodes={max_nodes}. "
                  f"Showing the {max_nodes} ERCs closest to the current generator "
                  f"instead (raise max_nodes to force the full network).")
            if self.generator:
                node_set: set[int] = set()
                radius = 1
                while len(node_set) < max_nodes and radius < n:
                    node_set = set(self.generator)
                    for c in self.generator:
                        node_set |= self.neighbors(c, radius=radius, relations=relations)
                    radius += 1
                node_ids = sorted(node_set)[:max_nodes]
            else:
                node_ids = list(range(min(max_nodes, n)))
        else:
            node_ids = list(range(n))
        out = filename or "context_view.html"
        return self._plot(node_ids, relations, out, synergy_layout=synergy_layout)

    def _plot(
        self,
        node_ids: list[int],
        relations: tuple,
        filename: str,
        synergy_layout: str = "below",
        anchor_set: set[int] | None = None,
        reasons: dict[int, str] | None = None,
    ) -> str:
        """
        Always laid out hierarchically by ERC containment level (see
        _compute_containment_levels): ERCs are never scattered by a
        force-directed simulation, so the containment structure the network
        actually has is always what you see. Fundamental synergies are not
        drawn as edges directly between their two reactant ERCs (which
        would cut straight across the hierarchy and land visually on top of
        whatever ERC nodes happen to be in between); each is its own
        diamond "reaction" node.

        synergy_layout controls where that diamond row sits:
          "below"       (default) -- every diamond pinned to one dedicated
                         level below the lowest ERC level shown, so all
                         synergies show up organized in a single row
                         beneath the hierarchy rather than scattered over
                         it.
          "interleaved" -- each diamond sits in the gap directly above the
                         LOWER of its two reactants' containment levels
                         (halfway between that level and the next one up),
                         so synergies are spread across one inter-level
                         band per hierarchy depth instead of one shared
                         band under everything. The band with the most
                         diamonds is where fundamental synergies are firing
                         most often, level-wise.

        anchor_set, if given (generator_view uses this; local_view and
        context_view leave it None), restricts every synergy diamond and
        every complementarity edge to ones touching at least one anchor
        node -- so a relation between two peripheral nodes that both merely
        happen to be shown is never drawn, only relations that actually
        connect to the anchored closure. Nodes in anchor_set but not in
        self.generator are colored as "implied" (teal), distinct from both
        explicit generator members (green) and further candidates (yellow).

        reasons, if given, is a {erc_idx: text} map merged into each node's
        hover tooltip, so every non-generator node's title says exactly
        which relation put it on the graph.
        """
        from pyvis.network import Network

        net = Network(height="800px", width="100%", directed=True, notebook=False)
        net.set_options("""
        { "physics": {"enabled": false},
          "layout": {"hierarchical": {"enabled": true, "direction": "DU",
                                       "sortMethod": "hubsize",
                                       "levelSeparation": 140,
                                       "nodeSpacing": 110} } }
        """)

        if synergy_layout not in ("below", "interleaved"):
            raise ValueError(
                f"synergy_layout must be 'below' or 'interleaved', got {synergy_layout!r}")

        node_set = set(node_ids)
        min_erc_level = min((self._containment_level[i] for i in node_ids), default=0)
        syn_level = min_erc_level - 1  # "below" mode: one dedicated band under everything
        gen_set = set(self.generator)
        if anchor_set is not None:
            implied_set = anchor_set - gen_set   # e.g. descendants/closure members, teal
            cand_set = node_set - anchor_set      # genuine one-hop extension candidates, yellow
        else:
            implied_set = set()
            cand_set = {i for items in self._compute_candidates().values() for i, _ in items}

        for i in node_ids:
            e = self.ercs[i]
            if i in gen_set:
                color = "#27ae60"   # green — explicit generator member
            elif i in implied_set:
                color = "#48c9b0"   # teal — already inside the closure (descendant/synergy-cascade)
            elif i in cand_set:
                color = "#f1c40f"   # yellow — a one-hop candidate, not yet in the closure
            elif e.is_persistent():
                color = "#5dade2"   # light blue — persistent ERC
            else:
                color = "#d0d3d4"   # gray — everything else
            why = f"<br>why shown: {reasons[i]}" if reasons and i in reasons else ""
            title = (f"{self.erc_label(i)}  size={e.size()}  "
                     f"persistent={e.is_persistent()}<br>"
                     f"species: {self.erc_species_names(i)}<br>"
                     f"req: {self.species_names(e.req_mask)}{why}")
            net.add_node(self.erc_label(i), label=self.erc_label(i), title=title,
                         color=color, size=18 + 2 * e.size(),
                         level=self._containment_level[i])

        if "containment" in relations:
            # Direction is child -> parent (the smaller ERC's arrow points at
            # the bigger one it's a subset of); the "⊂" label makes that
            # explicit on the graph itself, not just on hover, since
            # "which end is contained in which" is otherwise easy to guess
            # backwards even with the arrowhead visible.
            for i in node_ids:
                for j in self.hier.parents[i]:
                    if j in node_set:
                        net.add_edge(self.erc_label(i), self.erc_label(j),
                                     color=_REL_COLOR["containment"], width=2,
                                     label="⊂", font={"color": _REL_COLOR["containment"], "size": 12},
                                     title=f"{self.erc_label(i)} ⊂ {self.erc_label(j)}  "
                                           f"({self.erc_label(i)} is contained in {self.erc_label(j)})")
        if "synergy" in relations:
            # Uses self._syn_by_pair (built from self.syn.fundamental only —
            # never .basic/.maximal), which groups by reactant PAIR: a pair
            # can have several simultaneous fundamental targets (an
            # antichain, per the companion theory), and they are drawn as
            # ONE diamond "reaction" node per pair, with edges out to every
            # one of its targets — not one diamond per (pair, target)
            # triple, so a pair that jointly ignites several ERCs at once
            # reads as the single relation it is.  (FundamentalGraph.syn_from
            # itself stays unaggregated for search-performance reasons — see
            # fundamental_graph.py — this grouping is display-only.)
            seen_pairs: set[tuple[int, int]] = set()
            for i in node_ids:
                for j in self._syn_partners.get(i, ()):
                    pair_key = (min(i, j), max(i, j))
                    if pair_key in seen_pairs or j not in node_set:
                        continue
                    if anchor_set is not None and i not in anchor_set and j not in anchor_set:
                        continue  # neither reactant touches the anchored closure -- skip
                    seen_pairs.add(pair_key)
                    targets_in_set = sorted(
                        k for k in self._syn_by_pair[pair_key] if k in node_set)
                    if not targets_in_set:
                        continue
                    target_lbls = ", ".join(self.erc_label(k) for k in targets_in_set)
                    syn_node = f"syn:{i}+{j}"
                    if synergy_layout == "interleaved":
                        # In the gap directly above the LOWER reactant's own
                        # containment level -- i.e. between that level and
                        # the next one up -- so the diamond's row reflects
                        # where in the hierarchy this particular synergy
                        # fires, and a crowded band shows a depth at which
                        # fundamental synergies are especially common.
                        pair_level = min(self._containment_level[i], self._containment_level[j])
                        node_level = pair_level + 0.5
                    else:
                        node_level = syn_level
                    net.add_node(syn_node, label="+", shape="diamond", size=12,
                                 color=_REL_COLOR["synergy"], level=node_level,
                                 title=f"fundamental synergy: {self.erc_label(i)} + "
                                       f"{self.erc_label(j)} -> {target_lbls}"
                                       f"{f'  (level {pair_level}–{pair_level+1})' if synergy_layout == 'interleaved' else ''}")
                    net.add_edge(self.erc_label(i), syn_node,
                                 color=_REL_COLOR["synergy"], width=2)
                    net.add_edge(self.erc_label(j), syn_node,
                                 color=_REL_COLOR["synergy"], width=2)
                    for k in targets_in_set:
                        net.add_edge(syn_node, self.erc_label(k),
                                     color=_REL_COLOR["synergy"], width=2)
        if "complementarity" in relations:
            # Each FundComp(prod_idx, cons_idx, species) is one direction only
            # (prod_idx supplies species to cons_idx); if a pair is
            # complementary in both directions, two separate FundComp records
            # exist and are drawn as two separately curved arrows so neither
            # hides the other.
            for fc in self.comp.fundamental:
                if fc.prod_idx in node_set and fc.cons_idx in node_set:
                    if (anchor_set is not None
                            and fc.prod_idx not in anchor_set and fc.cons_idx not in anchor_set):
                        continue  # neither end touches the anchored closure -- skip
                    sp_name = self.species_names(1 << fc.species)[0]
                    net.add_edge(
                        self.erc_label(fc.prod_idx), self.erc_label(fc.cons_idx),
                        color=_REL_COLOR["complementarity"], width=1.5, dashes=True,
                        label=f"provides {sp_name}",
                        title=f"{self.erc_label(fc.prod_idx)} provides '{sp_name}' "
                              f"to {self.erc_label(fc.cons_idx)}",
                        smooth={"type": "curvedCW", "roundness": 0.2},
                    )

        out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "outputs", "explorer_html")
        os.makedirs(out_dir, exist_ok=True)
        full_path = os.path.abspath(os.path.join(out_dir, filename))
        net.html = net.generate_html()
        with open(full_path, "w", encoding="utf-8") as f:
            f.write(net.html)
        print(f"  Saved: {full_path}")
        webbrowser.open(f"file://{full_path}")
        return full_path
