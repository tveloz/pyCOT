# Endosymbiotic Fusion of Elementary Persistent Modules — a modelling framework

Status: working theory document, written before implementation. Updated as the
project progresses; see `PROGRESS.md` for the current state of the actual
computational work.

## 1. Motivation

The companion paper (`projects/COT_Fundamental_Generators_Exploration/paper_ecoli_comparative`)
found that supplying an *E. coli* genome-scale reconstruction with its full native
medium — trace metals and cofactor precursors absent from a minimal core medium —
does not add new, independent elementary persistent modules (EPMs). Instead it
**fuses** existing ones: EPM count falls while mean EPM size roughly triples. The
mechanistic reading was that Mo/W/Ni/Se-cofactor-dependent respiratory enzymes sit
at junctions between what were, under the core medium, two separately
self-sustaining modules; supplying the cofactor lets those modules act as one.

That finding raises two questions this project takes up:

1. **Is fusion a generic, searchable property of certain inflows, or was it
   specific to E. coli's trace metals?** Given any reaction network and a
   baseline food set, can we systematically search over candidate additional
   inflow species and rank them by how much they fuse (rather than merely
   extend) the network's elementary persistent structure?

2. **What happens if the "inflow" is not a single abiotic species but an
   entire second living system's worth of continuously-regenerated metabolic
   output?** This is a structural description of endosymbiosis: a host and a
   symbiont, each independently self-maintaining under their own native food
   sets, become physically coupled so that products of one are available to
   the other. Does this produce the same fusion signature as an abiotic
   cofactor addition, but at a qualitatively larger scale — and can the
   computable signature of that jump serve as a structural, falsifiable model
   of "major evolutionary transition"-style complexification?

## 2. Framework

### 2.1 Fusigenic inflow search (Question 1)

Given a compiled network `N` with elementary reaction closures (ERCs), a
containment hierarchy, fundamental synergies/complementarities, and EPMs
computed under baseline food set `F`:

For each candidate species `c` not already in `F` (or, for large networks, a
biologically filtered/sampled subset — e.g. species that are the *req_mask* of
at least two distinct persistent ERCs, since those are structurally the
candidates capable of bridging two closures at all), recompute EPMs under
`F ∪ {c}` and record:

- `Δcount = |EPMs(F∪{c})| − |EPMs(F)|`
- `Δmean_size`, `Δmax_size` (species-count basis, EPM-species-mask ∪ E0)
- a **fusion score**, operationalised as `−Δcount × Δmean_size` when
  `Δcount ≤ 0` and `Δmean_size > 0` (proliferation — `Δcount > 0` — is scored
  separately and is not fusion by this definition, even if some individual
  EPM grew)

Species with high fusion scores are candidates for the same kind of role the
E. coli paper's trace metals played: a bottleneck cofactor/precursor that
bridges two previously-separate closures without itself producing a large
number of new independent products. This is cheap: EPM computation on
networks in the hundreds-to-low-thousands of reactions completes in seconds,
so an exhaustive per-species sweep is tractable up to genome scale in minutes,
and importantly needs **no LP verification** — the same scope decision as the
main paper.

### 2.2 Endosymbiotic merger (Question 2)

Two independently self-maintaining networks:

- **Host** `H`, species `X_H`, reactions `R_H`, native food set `F_H`.
- **Symbiont** `S`, species `X_S`, reactions `R_S`, native food set `F_S`.

computed independently first: `EPMs(H, F_H)` and `EPMs(S, F_S)` are the
pre-fusion baseline.

**Compartmentalised union.** Naively concatenating `H` and `S`'s reaction
lists is *wrong*: BiGG-style species tokens like `atp_c` refer to a specific
cell's cytoplasm, and pre-endosymbiosis, host-cytoplasm ATP and
symbiont-cytoplasm ATP are physically different pools in different cells. We
therefore relabel every species in `S` with a distinguishing tag (e.g.
`atp_c` → `atp_c__endo`) before union, so `H`'s and `S`'s internal chemistries
stay physically distinct except where we deliberately connect them. This is
the modelling step that makes the framework a model of *endosymbiosis*
specifically, rather than of two networks that happen to share metabolite
names.

**Interface.** A small, explicit, curated set of transport/exchange reactions
connecting specific `H` and `S` species pairs, e.g. `pyr_c <=> pyr_c__endo`,
representing a real or hypothesised transporter. This is the framework's main
free parameter: it defines exactly what crosses the new compartment boundary,
in which direction(s), and can be swept from a minimal one-metabolite
interface to a rich multi-metabolite one to see how the fusion effect scales
with interface richness — directly analogous to the real evolutionary
trajectory from a loosely-associated symbiont (few transporters) to a fully
integrated organelle (many).

**Structural synergy reactions (optional, curated, explicitly flagged as
hypotheses).** Reactions present in *neither* `R_H` nor `R_S`, representing a
catalytic capability that becomes chemically possible only once host and
symbiont are co-localised — e.g. a host enzyme whose substrate was simply
never available in the free-living host now receiving it from the symbiont,
or a genuinely novel combination requiring real biochemical judgement about
enzyme promiscuity/substrate compatibility. Every such reaction used in the
realistic (non-toy) example must carry a stated biochemical justification; we
do not derive these automatically, in line with the user's framing that this
step "requires biochemical knowledge beyond the models themselves."

**Merged network** `M = H ⊎ S_relabelled ⊎ Interface ⊎ Synergy`, food set
`F_M = F_H ∪ F_S_relabelled` (the merger does not remove either partner's own
abiotic requirements — it adds a new internal supply channel between them, on
top of, not instead of, what each already needed from the environment).

### 2.3 Complexification metrics

Computed by comparing `EPMs(M, F_M)` against the pre-fusion baseline
`EPMs(H, F_H) ⊔ EPMs(S, F_S)`:

- **`ΔEPM_max`**: does the merged network's largest EPM exceed what either
  partner could achieve alone?
- **Hybrid EPMs**: EPMs of `M` whose species set has a non-empty intersection
  with *both* `X_H` and the relabelled `X_S` — i.e. modules that are only
  self-maintaining because they draw on machinery from both partners. This is
  the cleanest operational signature of genuine integration, as opposed to
  mere juxtaposition: two systems placed side by side without any interface
  reactions would simply reproduce `EPMs(H) ⊔ EPMs(S)` unchanged, with zero
  hybrid EPMs. Hybrid-EPM count and size are the project's primary
  "endosymbiotic complexification" readout.
- **EPM count change**: fusion (count falls, sizes grow) vs. mere
  proliferation (count grows via a few new small opportunistic modules using
  the interface, without deep integration) are structurally distinguishable
  and both biologically meaningful (the latter looks more like a loose,
  facultative/food-sharing association; the former looks more like the
  beginning of obligate integration).
- **Interface fragility**: removing one interface reaction at a time and
  recomputing identifies which specific cross-feeding link is load-bearing
  for which hybrid EPMs — a computable proxy for the real evolutionary
  question of which transporter, if lost, would revert or break the
  symbiosis, and which genes/functions are consequently under the strongest
  selective pressure to be retained (or transferred to the host genome) during
  organelle-genome reduction.

## 3. Worked plan

1. **Toy model**, deliberately modelled on the textbook mitochondrial
   endosymbiosis story so every step can be hand-verified: a small host
   network capturing glycolysis (substrate-level ATP production, ends at
   pyruvate, cannot itself extract much more energy) and a small symbiont
   network capturing the TCA cycle + a minimal electron-transport/oxidative
   phosphorylation stand-in (consumes pyruvate-derived carbon and O2,
   produces much more ATP). Interface: pyruvate host→symbiont, ATP (and/or a
   reducing-equivalent proxy) symbiont→host. This lets us confirm the whole
   pipeline — relabelling, interface construction, EPM comparison, hybrid-EPM
   detection — behaves as expected on a system small enough to check by hand
   before trusting it on real data.

2. **Fusigenic inflow search** run on the toy model and (time permitting) on
   one or two of the BiGG organisms already catalogued in
   `projects/COT_Fundamental_Generators_Complex/network_catalogue_bigg_organisms.csv`,
   to check whether the top-ranked fusigenic species are recognisable
   cofactors/precursors (consistency check against the E. coli paper) and
   whether any qualitatively different classes of fusigenic inflow turn up.

3. **Real-data endosymbiosis merger(s)**, built from real BiGG genome-scale
   reconstructions already available in this repository:
   - A **mitochondrial-type** merger: a fermentation-restricted host network
     (a subset of a heterotrophic organism's reconstruction with oxidative
     phosphorylation reactions removed, standing in for an anaerobic/
     microaerophilic host lineage) plus a full aerobic respirer as symbiont.
   - A **chloroplast-type** merger: a heterotrophic eukaryote host (e.g.
     *S. cerevisiae*, iND750 — already used in the companion paper) plus
     *Synechocystis* sp. PCC 6803 (iJN678 — also already used there) as a
     photosynthetic symbiont, directly mirroring the real evolutionary
     origin of the plastid.

   Both use the *same* generic merge/interface/hybrid-EPM tooling built for
   the toy model — the toy model is a correctness check on the machinery, not
   a separate implementation.

4. **Report** synthesising the theory, the toy-model validation, the real-data
   results, and an honest discussion of what is a computed result versus a
   curated modelling choice (the interface and any synergy reactions), with
   explicit falsifiable predictions this framework makes that could in
   principle be checked against the real evolutionary/genomic record (e.g.
   gene-retention patterns in reduced organelle/endosymbiont genomes should
   be enriched for reactions that participate in hybrid EPMs, if the
   framework's notion of "load-bearing interface" tracks real selective
   pressure).

## 4. Honesty notes

- This is a *structural, steady-state* framework, exactly like the companion
  paper: it says nothing about the population-genetic or ecological dynamics
  of how an endosymbiotic relationship is established or fixed, only about
  what happens to the elementary-persistent-module structure of the combined
  chemistry once a given interface exists.
- The interface and any synergy reactions are modelling choices, not data.
  Every one used in the real-data section will be stated explicitly with its
  biochemical rationale, and the report will present results conditional on
  those choices rather than as unconditional facts about the named organisms.
- BiGG genome-scale reconstructions are not literal models of the specific
  ancestral lineages involved in the Proterozoic mitochondrial or plastid
  endosymbiosis events; they are modern free-living relatives used as
  structural stand-ins for the metabolic *capabilities* (glycolysis-only vs.
  full oxidative phosphorylation; heterotrophy vs. oxygenic photosynthesis)
  that the real endosymbiosis narrative turns on. The report will say this
  plainly rather than implying a literal historical reconstruction.
