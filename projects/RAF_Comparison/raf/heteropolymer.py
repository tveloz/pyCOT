"""
heteropolymer.py -- binary-string heteropolymer network generator (the
classic Kauffman/Farmer/Bagley "binary polymer model" used throughout
Hordijk & Steel's RAF papers as their canonical synthetic test case),
rewritten from a first draft that had two structural problems:

1. LIGATION-ONLY NETWORKS CANNOT CONTAIN A FRAGILE CIRCUIT.
   Ligation A+B->AB is strictly length-increasing, so the "produced by"
   relation among species is a DAG ordered by chain length: nothing
   downstream of food can ever come back around to help produce something
   it depends on. Every closed set is therefore either not self-maintaining
   at all, or self-maintaining with F = everything and C = empty -- the
   decomposition (Def:overprod / Def:C in decomposing_RAF_v5.tex) has
   nothing to say. Cyclic/mutual structure (the whole point of fragile
   circuits) REQUIRES a reverse reaction that can consume a longer species
   to help produce a shorter one that in turn helps produce it. This module
   therefore adds the reverse cleavage AB->A+B for every ligation A+B->AB
   by default (`include_cleavage=True`) -- not a stylistic choice, a
   structural necessity if the network is meant to exercise the
   decomposition theorem at all.

2. ONE REACTION PER CATALYST WASTES REACTIONS AND MISREPRESENTS RAF
   SEMANTICS. RAF's C: R -> 2^M assigns a SET of alternative catalysts to
   ONE reaction (Def 2.1(b) in decomposing_RAF_v5.tex: "there is A catalyst
   c in C(r)", existential -- any one suffices). A first-draft version that
   emits a separate reaction line per candidate catalyst is semantically
   equivalent for RAF purposes but multiplies |R| by the average catalyst
   count for no benefit, and turns what should be one stoichiometric
   transformation into many duplicate columns of the self-maintenance LP.
   This module instead writes each reaction ONCE, listing every assigned
   catalyst as a net-zero species on both sides of that single line
   (`... + 1 c1 + 1 c2 => ... + 1 c1 + 1 c2`) and lets the SAME net-zero
   catalyst induction already used for kinetic/BIOMD networks
   (`biomodel_crs.crs_from_biomodel`) recover the catalyst set -- each
   catalyst independently satisfies `r.catalysts & W` in raf_algo.is_raf,
   giving correct OR semantics with no new induction code needed.

CATALYSIS MODES
----------------
  'uniform'      : classical Kauffman null model -- each (reaction,
                   candidate species) pair is independently catalyzed with
                   probability `p_catalyst`, no correlation with sequence
                   identity. Some reactions end up genuinely uncatalyzed
                   (no fallback bootstrap catalyst is forced in) -- keeping
                   catalysis rare and structurally meaningful, rather than
                   trivializing maxRAF into "almost the whole network".
  'template'     : deterministic structural (motif) catalysis -- species s
                   catalyzes A+B->AB (and its reverse) iff s occurs as a
                   contiguous substring of AB that straddles the ligation
                   junction (i.e. s starts at some position <= len(A)-1 and
                   ends at some position >= len(A)), modelling a
                   template/ribozyme that recognizes the splice site rather
                   than an arbitrary internal motif. No randomness at all;
                   reproducible from the sequences alone.
  'preferential' : degree-weighted random catalysis -- like 'uniform', but
                   each species' per-reaction selection probability is
                   scaled by its precomputed "structural degree" (how many
                   OTHER species/reaction-products it occurs in as a
                   substring), a static proxy for Barabasi-Albert-style
                   preferential attachment: already-central sequences are
                   more likely to pick up new catalytic roles, producing a
                   heavy-tailed catalyst-degree distribution instead of the
                   uniform model's flat one.

COMBINATORIAL CONTROL
----------------------
Species grow as 2^(L+1)-2 and ligation pairs (ordered, since AB != BA in
general) as roughly O(4^L) before the length filter, so keep `max_length`
small (default 4). `reaction_keep_prob` additionally subsamples WHICH
chemically-possible ligations are actually instantiated (independent of
catalysis), giving a second, more standard RAF-literature knob: network
density, separate from mean catalyst count.
"""
from __future__ import annotations

import random
from dataclasses import dataclass, field


@dataclass
class HeteropolymerNetwork:
    food: list[str]
    species: list[str]
    txt: str
    stats: dict = field(default_factory=dict)


def _all_species(max_length: int) -> list[str]:
    species = []
    for length in range(1, max_length + 1):
        for i in range(2 ** length):
            species.append(format(i, f"0{length}b"))
    return species


def _structural_degree(species: list[str]) -> dict[str, int]:
    """How many times each species occurs as a substring of some OTHER
    species in the pool -- a cheap static proxy for 'already central'."""
    deg = {s: 0 for s in species}
    for s in species:
        for t in species:
            if s != t and s in t:
                deg[s] += 1
    return deg


def _template_catalysts(A: str, B: str, product: str, species_set: set[str]) -> list[str]:
    """Species straddling the A|B junction of `product`, restricted to
    length >= 2 (a length-1 'catalyst' would just be an arbitrary monomer
    present in the product, not a real motif match)."""
    j = len(A)
    n = len(product)
    cats = []
    for start in range(0, j):          # start <= j-1: reaches into A
        for end in range(j + 1, n + 1):  # end >= j+1: reaches into B
            motif = product[start:end]
            if motif in species_set:
                cats.append(motif)
    return cats


def generate_heteropolymer_network(
    max_length: int = 4,
    *,
    catalysis_mode: str = "uniform",   # 'uniform' | 'template' | 'preferential'
    p_catalyst: float = 0.05,
    include_cleavage: bool = True,
    cleavage_prob: float = 1.0,
    reaction_keep_prob: float = 1.0,
    seed: int | None = 0,
) -> HeteropolymerNetwork:
    """Build a binary-string heteropolymer network and its pyCOT .txt
    encoding. See module docstring for the two structural fixes relative
    to the first-draft generator (cleavage reactions, single-line
    multi-catalyst reactions) and for the three catalysis modes.

    `cleavage_prob` (only relevant when `include_cleavage=True`) is an
    INDEPENDENT per-ligation coin flip deciding whether that reaction's
    reverse is included, distinct from `reaction_keep_prob` (which decides
    whether the ligation itself exists at all). Rmk:ligation-acyclic only
    requires SOME length-decreasing reaction to exist, not a matched
    reverse for every ligation -- and universal 1:1 reversibility is a
    modelling liability, not just an idealization: whenever a reaction and
    its exact reverse share identical stoichiometry, setting equal flux on
    both gives Sv=0 for that pair for free, in ANY closed set containing
    it, which trivializes self-maintenance network-wide (empirically,
    `cleavage_prob=1.0` at L=3/uniform/p=0.2 makes every one of the 21
    discovered semi-organizations already an organization -- the
    semi-organization/organization gap the theory is built to detect never
    shows up). `cleavage_prob<1.0` keeps the necessity argument intact
    while no longer handing out that escape hatch everywhere.
    """
    if catalysis_mode not in ("uniform", "template", "preferential"):
        raise ValueError(f"unknown catalysis_mode: {catalysis_mode!r}")

    rng = random.Random(seed)
    species = _all_species(max_length)
    species_set = set(species)
    food = [s for s in species if len(s) == 1]

    degree = _structural_degree(species) if catalysis_mode == "preferential" else None
    max_deg = max(degree.values()) if degree and max(degree.values()) > 0 else 1

    def choose_catalysts(A: str, B: str, product: str) -> list[str]:
        if catalysis_mode == "template":
            return sorted(set(_template_catalysts(A, B, product, species_set)))
        if catalysis_mode == "uniform":
            return [s for s in species if rng.random() < p_catalyst]
        # preferential: weight by structural degree, normalized so the
        # AVERAGE selection probability across species still equals
        # p_catalyst (keeps modes comparable at fixed p_catalyst).
        mean_deg = sum(degree.values()) / len(degree)
        mean_deg = mean_deg if mean_deg > 0 else 1
        chosen = []
        for s in species:
            weight = (degree[s] + 1) / (mean_deg + 1)
            if rng.random() < min(1.0, p_catalyst * weight):
                chosen.append(s)
        return chosen

    lines = ["# Heteropolymer network (binary strings, length 1.." + str(max_length) + ")",
             f"# catalysis_mode={catalysis_mode}  p_catalyst={p_catalyst}  "
             f"include_cleavage={include_cleavage}  cleavage_prob={cleavage_prob}  "
             f"reaction_keep_prob={reaction_keep_prob}  seed={seed}"]

    for i, tok in enumerate(food):
        lines.append(f"INFLOW_{i}:  => 1 {tok};")

    r_id = 0
    n_ligations_possible = 0
    n_ligations_kept = 0
    n_with_catalyst = 0
    catalyst_counts = []
    for A in species:
        for B in species:
            product = A + B
            if len(product) > max_length:
                continue
            n_ligations_possible += 1
            if reaction_keep_prob < 1.0 and rng.random() >= reaction_keep_prob:
                continue
            n_ligations_kept += 1

            cats = choose_catalysts(A, B, product)
            catalyst_counts.append(len(cats))
            if cats:
                n_with_catalyst += 1
            cat_terms = "".join(f" + 1 {c}" for c in cats)

            r_id += 1
            fwd = f"r{r_id}_lig: 1 {A} + 1 {B}{cat_terms} => 1 {product}{cat_terms};"
            lines.append(fwd)
            if include_cleavage and (cleavage_prob >= 1.0 or rng.random() < cleavage_prob):
                r_id += 1
                rev = f"r{r_id}_cle: 1 {product}{cat_terms} => 1 {A} + 1 {B}{cat_terms};"
                lines.append(rev)

    txt = "\n".join(lines) + "\n"
    stats = dict(
        n_species=len(species),
        n_food=len(food),
        n_ligations_possible=n_ligations_possible,
        n_ligations_instantiated=n_ligations_kept,
        n_reactions_total=r_id,
        n_ligations_with_catalyst=n_with_catalyst,
        fraction_ligations_catalyzed=(n_with_catalyst / n_ligations_kept) if n_ligations_kept else 0.0,
        mean_catalysts_per_ligation=(sum(catalyst_counts) / len(catalyst_counts)) if catalyst_counts else 0.0,
    )
    return HeteropolymerNetwork(food=food, species=species, txt=txt, stats=stats)
