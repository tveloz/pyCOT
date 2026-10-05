"""
metanetwork.py — Combined structure: hierarchy, synergy, complementarity,
elementary semi-organizations (Stage 4).

Bundles all computed structures into one queryable object with summary statistics
and export helpers for visualization or downstream analysis.

Public API
----------
build_metanetwork(ercs, hier, syn, comp, gen, elementary=None) -> MetaNetwork

MetaNetwork:
  .stats()        -> dict   — flat dict of scalar statistics
  .to_node_list() -> list[dict]
  .to_edge_list() -> list[dict]
  .print_summary()
"""
from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class MetaNetwork:
    """
    Aggregate of all Stage 1–4 outputs for a reaction network.

    Attributes
    ----------
    ercs       : list[ERCData]
    hier       : HierarchyData
    syn        : SynergyResult  (None when skipped)
    comp       : CompResult
    gen        : GeneratorResult (None when skipped)
    elementary : ElementarySOResult (None when not computed)
    """
    ercs: list
    hier: object
    syn:  object
    comp: object
    gen:  object
    elementary: object = None

    # ------------------------------------------------------------------ stats

    def stats(self) -> dict:
        n = len(self.ercs)
        n_persistent = sum(1 for e in self.ercs if e.is_persistent())
        n_pairs = n * (n - 1) // 2
        n_comparable = sum(len(self.hier.ancestors[i]) for i in range(self.hier.n))
        n_incomparable = n_pairs - n_comparable
        hasse_edges = sum(len(self.hier.parents[i]) for i in range(self.hier.n))

        if self.syn is not None:
            n_basic = len(self.syn.basic)
            n_maxim = len(self.syn.maximal)
            n_fund  = len(self.syn.fundamental)
        else:
            n_basic = n_maxim = n_fund = None

        if self.gen is not None:
            n_prim   = len(self.gen.primitive_indices)
            n_reach  = len(self.gen.basis_reach)
            coverage = self.gen.coverage
            complete = self.gen.is_complete()
        else:
            n_prim = n_reach = coverage = complete = None

        elem_single = elem_multi = elem_total = None
        if self.elementary is not None:
            elem_single = len(self.elementary.single_erc_indices)
            elem_multi  = len(self.elementary.multi_erc_masks)
            elem_total  = len(self.elementary.all_elementary_masks)

        return {
            "n_ercs":              n,
            "n_persistent_ercs":   n_persistent,
            "n_non_persistent":    n - n_persistent,
            "n_pairs":             n_pairs,
            "n_comparable":        n_comparable,
            "n_incomparable":      n_incomparable,
            "n_hasse_edges":       hasse_edges,
            "n_basic_syn":         n_basic,
            "n_maximal_syn":       n_maxim,
            "n_fundamental_syn":   n_fund,
            "n_comp_basic":        len(self.comp.basic),
            "n_comp_pure":         len(self.comp.pure),
            "n_comp_fundamental":  len(self.comp.fundamental),
            "n_primitives":        n_prim,
            "n_reachable":         n_reach,
            "coverage":            coverage,
            "basis_complete":      complete,
            "n_elementary_single": elem_single,
            "n_elementary_multi":  elem_multi,
            "n_elementary_total":  elem_total,
        }

    # ----------------------------------------------------------------- nodes

    def to_node_list(self) -> list[dict]:
        primitives  = set(self.gen.primitive_indices) if self.gen is not None else set()
        basis_reach = self.gen.basis_reach             if self.gen is not None else frozenset()
        elem_set    = set(self.elementary.single_erc_indices) if self.elementary is not None else set()
        nodes = []
        for i, e in enumerate(self.ercs):
            nodes.append({
                "idx":                  i,
                "erc_id":               e.erc_id,
                "size":                 e.size(),
                "is_persistent":        e.is_persistent(),
                "is_primitive":         i in primitives,
                "in_basis_reach":       i in basis_reach,
                "is_single_elementary": i in elem_set,
                "n_reactions":          len(e.reaction_indices),
                "n_min_bases":          len(e.min_bases),
                "req_popcount":         bin(e.req_mask).count('1'),
                "prod_popcount":        bin(e.prod_mask).count('1'),
            })
        return nodes

    # ----------------------------------------------------------------- edges

    def to_edge_list(self) -> list[dict]:
        edges = []

        # Hasse edges
        for i in range(self.hier.n):
            for j in self.hier.parents[i]:
                edges.append({"src": i, "tgt": j, "rel_type": "hasse"})

        # Synergy
        if self.syn is not None:
            for s in self.syn.fundamental:
                edges.append({"src": s.i, "tgt": s.j,
                              "rel_type": "syn_fundamental", "via": s.k})
            for s in self.syn.maximal:
                edges.append({"src": s.i, "tgt": s.j,
                              "rel_type": "syn_maximal", "via": s.k})
            for s in self.syn.basic:
                edges.append({"src": s.i, "tgt": s.j,
                              "rel_type": "syn_basic", "via": s.k})

        # Complementarity (paper definitions)
        for c in self.comp.fundamental:
            edges.append({"src": c.prod_idx, "tgt": c.cons_idx,
                          "rel_type": "comp_fundamental", "species": c.species})
        for c in self.comp.basic:
            edges.append({"src": c.i, "tgt": c.j,
                          "rel_type": "comp_basic",
                          "fwd_supply": c.fwd_supply, "bwd_supply": c.bwd_supply})

        return edges

    # --------------------------------------------------------- summary print

    def print_summary(self, prefix: str = "") -> None:
        s = self.stats()
        p = prefix
        print(f"{p}ERCs:              {s['n_ercs']}  "
              f"({s['n_persistent_ercs']} persistent, "
              f"{s['n_non_persistent']} non-persistent)")
        print(f"{p}Hierarchy:         {s['n_hasse_edges']} Hasse edges  "
              f"({s['n_comparable']} comparable pairs)")
        if s['n_fundamental_syn'] is None:
            print(f"{p}Synergy:           (skipped — network too large)")
        else:
            print(f"{p}Synergy:           {s['n_fundamental_syn']} fundamental  "
                  f"/ {s['n_maximal_syn']} maximal  "
                  f"/ {s['n_basic_syn']} basic")
        print(f"{p}Complementarity:   {s['n_comp_basic']} basic  "
              f"/ {s['n_comp_pure']} pure  "
              f"/ {s['n_comp_fundamental']} fundamental")
        if s['n_primitives'] is None:
            print(f"{p}Generators:        (skipped)")
        else:
            print(f"{p}Generators:        {s['n_primitives']} primitive  "
                  f"coverage={s['coverage']:.1%}  "
                  f"complete={s['basis_complete']}")
        if s['n_elementary_total'] is not None:
            print(f"{p}Elementary SOs:    {s['n_elementary_total']} total  "
                  f"({s['n_elementary_single']} single-ERC, "
                  f"{s['n_elementary_multi']} multi-ERC)")


# ---------------------------------------------------------------------------
# Factory
# ---------------------------------------------------------------------------

def build_metanetwork(ercs, hier, syn, comp, gen, elementary=None) -> MetaNetwork:
    return MetaNetwork(ercs=ercs, hier=hier, syn=syn, comp=comp, gen=gen, elementary=elementary)
