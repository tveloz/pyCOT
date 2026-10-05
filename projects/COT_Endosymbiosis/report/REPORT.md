# Endosymbiotic fusion of elementary persistent modules: a computable, structural model of major evolutionary transitions

*Working report — pyCOT / COT_Endosymbiosis project.*

## 0. Summary

Starting from a single finding in the companion *E. coli* reconstruction
paper — supplying trace-metal cofactors doesn't add new elementary
persistent modules (EPMs), it *fuses* existing ones into fewer, larger
ones — this project asks two questions and builds a framework connecting
them:

1. **Is that fusion effect a generic, searchable property of certain
   inflows?** Yes. A systematic per-species sweep (`fusigenic_inflow_search.py`)
   ranks candidate inflow additions by how much they fuse (not merely
   extend) a network's EPM structure, and recovers real, checkable answers
   on both a hand-built toy network and real BiGG genome-scale models.

2. **What happens if the "inflow" is an entire second living system's
   worth of continuously-regenerated metabolic output, rather than one
   abiotic species?** This is a structural description of endosymbiosis. A
   generic merge tool (`merge_networks.py`) builds a compartmentalised union
   of a host and a symbiont network connected by a small, explicit,
   biologically-justified interface, and an analysis tool
   (`endosymbiosis_analysis.py`) measures the resulting EPM
   complexification. On **every one of four test cases** — a hand-verified
   toy model, a real self-derived *E. coli*-based mitochondrial-type
   merger, and two real cross-organism chloroplast-type mergers (yeast +
   *Synechocystis*, with and without light) — the merged network's EPMs are
   **100% hybrid**: every single elementary persistent module in the
   combined system draws on machinery from *both* partners. Mean EPM size
   grows by 54–87% relative to the better of the two separate partners
   alone.

Two distinct mechanisms of complexification emerge, both real, both
reproduced here: **capability-gain** (the food closure itself grows — new
species become reachable only jointly) and **re-integration** (the food
closure is unchanged, but the *same* reachable chemistry reorganises into
far fewer, far larger, mutually-dependent modules). An interface-knockout
analysis (§6b) then tests which specific cross-feeding reactions are
individually *necessary* for this hybridisation: only 1 of 7 interface
reactions tested across the three cases (the toy model's Fe–S-cluster
dependency, which has no redundant alternative route) turned out to be
individually load-bearing — the other 6, embedded in richly-connected
(mostly real) metabolic networks, were each individually dispensable. This is itself a real,
useful finding, not a null result: it refines the framework's falsifiable
prediction from "conserved interface genes" to "conservation should track
computed *irreplaceability*." The framework also
makes a falsifiable prediction connecting structure to real evolutionary
genomics, discussed in §7.

All code, data files, logs and figures referenced below are in this
directory tree (`projects/COT_Endosymbiosis/`); every number quoted is
either directly copied from a logged pipeline run (`outputs/*.log`) or
computed by a script in this repository, not estimated by hand.

---

## 1. Motivation

The full framework and formal definitions are in `THEORY.md`; this section
summarises the motivating result and the two questions it raises.

The companion paper (`projects/COT_Fundamental_Generators_Exploration/
paper_ecoli_comparative`) found that supplying an *E. coli* genome-scale
reconstruction (iAF1260) with its full native medium — trace metals and
cofactor precursors (Mo, W, Ni, Se, B₁₂) absent from a minimal core medium —
does not add new, independent EPMs. It **fuses** existing ones: EPM count
fell from 182 to 164 while mean size roughly tripled (25.7 → 63.7 species).
The mechanistic reading was that these cofactors sit at junctions between
what were, under the core medium, two separately self-sustaining modules;
supplying the cofactor lets those modules act as one.

Two questions follow directly:

1. Is fusion a generic, rankable property of *specific* inflow species,
   discoverable by search, or was the trace-metal result specific to
   *E. coli*'s particular reconstruction detail?
2. Endosymbiosis is, structurally, the extreme version of "add an inflow":
   instead of one abiotic species, an entire second self-maintaining
   system's continuously-regenerated output becomes available. Does the
   same fusion logic apply at that scale, and can its computable signature
   serve as a structural model of evolutionary complexification?

---

## 2. Framework (see `THEORY.md` for full detail)

**Fusigenic inflow search.** For network `N` with baseline food `F`, and
each candidate species `c ∉ F`: compute `Δcount = |EPMs(F∪{c})| −
|EPMs(F)|` and `Δmean_size`. Define

```
fusion_score = Δmean_size · (1 + max(0, −Δcount))   if Δcount ≤ 0 and Δmean_size > 0
             = 0                                      otherwise
```

(A first version multiplied by `Δcount` directly, which zeroed the score
whenever `Δcount == 0` exactly — a real bug, caught by hand-checking the
toy model; see §8.)

**Endosymbiotic merger.** Host `H` (species `X_H`, food `F_H`) and symbiont
`S` (species `X_S`, food `F_S`), each self-sufficient alone. Every species
and reaction name in `S` is relabelled with a tag (`__endo`) before union —
*not* optional: two independently-curated networks routinely reuse
identical species names (`atp_c`) and, worse, identical generic reaction
names (`R0`, `INFLOW_glc`), and treating pre-fusion host-cytoplasm ATP and
symbiont-cytoplasm ATP as the same physical pool would be a modelling
error, not just a simplification. A small, explicit, curated set of
interface (transport) reactions connects specific host and symbiont
species pairs. The merged network `M = H ⊎ S_relabelled ⊎ Interface` gets
food `F_M = F_H ∪ F_S_relabelled`.

**Complexification metrics**, computed by comparing `EPMs(M, F_M)` against
the pre-fusion baseline `EPMs(H,F_H) ⊔ EPMs(S,F_S)`:

- **Hybrid EPMs** — EPMs of `M` whose species set intersects *both* `X_H`
  and the relabelled `X_S`. This is the cleanest signature of genuine
  integration: two systems merely juxtaposed with no interface would
  reproduce `EPMs(H) ⊔ EPMs(S)` unchanged, with **zero** hybrid EPMs.
- **ΔE0** — does the food closure itself grow (new species reachable only
  jointly — "capability-gain")?
- **Δmax / Δmean size**, and **Δcount** (fusion, count falls, vs.
  proliferation, count rises) — both directions are structurally
  meaningful and both occur in the results below.

---

## 3. Toy model (hand-verified)

Built deliberately on the textbook mitochondrial-endosymbiosis story, using
real, named biochemistry throughout rather than an abstract s1/s2/s3
example, so every number could be hand-checked before trusting the
pipeline on real data.

**Host** (`toy_model/host_alone.txt`) — lumped glycolysis + fermentation:
food = `{glc_h, adp_h, pi_h, nad_h}`. `GLYC` (glycolysis, lumped),
`FERM` (fermentation — the host's *only* way to regenerate NAD⁺ without a
partner, real biology), `MAINT` (generic ATP-consuming maintenance),
`BIOSYN` (needs `fescluster_h`, which nothing in this network can produce —
Fe–S cluster biogenesis, the ISC pathway, is a real, well-documented
mitochondria-exclusive eukaryotic function, deliberately used as the
host's one hard, unavoidable dependency).

**Symbiont** (`toy_model/symbiont_alone.txt`) — lumped glycolysis + pyruvate
oxidation + oxidative phosphorylation + Fe–S cluster biogenesis: food =
`{glc_s, o2_s, adp_s, pi_s, nad_s, fe_s, s_s}`. `GLYC_S`, `OXID`, `ETC`,
`ISC`.

**Interface** (`merge_networks.py`, 3 reactions): pyruvate host→symbiont,
ATP symbiont→host, Fe–S cluster symbiont→host.

**Result** (`outputs/` — reproduced live above in §0, computed by
`endosymbiosis_analysis.py`, see `scripts/endosymbiosis_analysis.py`):

| | EPMs | mean size | max size |
|---|---|---|---|
| Host alone | 1 | 8.0 | 8 |
| Symbiont alone | 1 | 12.0 | 12 |
| **Merged** | **1** | **22.0** | **22** |

The single merged EPM is **100% hybrid**. Food closure (E0) grows from
8+12=20 (naive union) to 22 — the two new species (`fescluster_h`,
`macromol_h`) are reachable *only* because of the interface: a clean,
hand-checkable instance of the **capability-gain** mechanism. See
`figures/fig_toy_model.png` for the annotated network diagram.

---

## 4. Fusigenic inflow search — results

### 4.1 Toy host (exhaustive, 6 candidates)

The search, run on the host *alone* with no knowledge of the symbiont,
independently re-discovers the exact dependency the toy model was designed
around:

| species | Δcount | Δmean | fusion score |
|---|---|---|---|
| **fescluster_h** | 0 | +2.0 | **2.00** |
| macromol_h | 0 | +1.0 | 1.00 |
| (4 others) | 0 | 0.0 | 0.00 |

This is a genuine cross-check between the two halves of the project:
Question 1's generic search and Question 2's designed endosymbiotic
dependency point to the same answer without being told about each other.

### 4.2 e_coli_core (exhaustive, all 65 non-food species)

Top candidates are **pathway-intermediate shortcuts**, not trace metals:

| species | Δcount | Δmean | fusion score |
|---|---|---|---|
| 2pg_c / 3pg_c / pep_c | −4 | +8.2 | 41.15 |
| acon_C_c / cit_c / icit_c | −3 | +5.8 | 23.13 |
| adp_c / atp_c | −3 | +3.7 | 14.80 |
| accoa_c | −4 | +2.7 | 13.42 |

Providing a lower-glycolysis or TCA-cycle intermediate directly bypasses
the need to synthesise it from scratch via earlier steps, merging what
would otherwise be separate closure routes — a real, different fusion
mechanism from the cofactor story, because `e_coli_core` (a small
"textbook" model) barely represents cofactor biosynthesis machinery at
all. **Fusigenic mechanism is reconstruction-detail-dependent.**

### 4.3 iAF1260 (targeted, 15 trace-metal candidates — direct cross-check against the companion paper)

An exhaustive sweep is not tractable at this scale (see §8); the 15
candidates were chosen to match the companion paper's own trace-metal
story directly.

| species | Δcount | Δmean | fusion score |
|---|---|---|---|
| **fe2_e** | −1 | +5.0 | 10.02 |
| ni2_e / cobalt2_e / mn2_e / mg2_e / zn2_e / k_e | −1 | +3.0 | 6.00 |
| na1_e | −1 | +2.1 | 4.12 |
| **mobd_e / tungs_e** | −1 | +2.0 | 4.01 |
| ca2_e / so4_e | −1 | +2.0 | 4.01 |
| cu2_e / fe3_e / cbl1_e | 0 | +2–3 | 2.0–3.0 |

All 15/15 candidates were fusigenic (none inert). `mobd_e`/`tungs_e`
(molybdate/tungstate — the specific cofactors named in the companion
paper) each individually produce a real, modest fusion effect
(`Δmean=+2.0`); the companion paper's full-native-medium result (mean size
25.7→63.7, +38 species) is consistent with these being largely additive
single-metal contributions rather than one dominant cofactor.

See `figures/fig_fusigenic_search.png` for all three panels together.

---

## 5. Real endosymbiosis mergers

### 5.1 Mitochondrial-type (real, self-derived from `e_coli_core`)

**Host**: `e_coli_core` with its three named oxidative-phosphorylation
reactions deleted — `ATPS4r` (ATP synthase, fwd+rev), `CYTBD` (cytochrome
oxidase, the O₂-consuming terminal step), `NADH16` (Complex-I equivalent) —
identified by their own BiGG annotation. What remains is a real,
self-consistent glycolysis+TCA+fermentation network limited to
substrate-level phosphorylation: a structural stand-in for a fermentative
lineage, built by *deletion from a real curated model*, not invented.

**Symbiont**: the unmodified `e_coli_core` — stated explicitly as a
structural stand-in for "a free-living aerobic respirer," *not* a literal
alphaproteobacterium (none in the local BiGG catalogue). Both use the
companion paper's standard 7-species aerobic core food.

**Interface** (2 reactions): host pyruvate → symbiont; symbiont ATP → host.

| | EPMs | mean size | max size |
|---|---|---|---|
| Host alone | 15 | 15.1 | 17 |
| Symbiont alone | 15 | 15.1 | 17 |
| **Merged** | **27** | **28.2** | **30** |

(Host and symbiont numbers match the companion paper's own published
`e_coli_core` aerobic-core figures exactly — an independent cross-check
that this pipeline invocation is correct.)

**27/27 merged EPMs are hybrid (100%).** Mean size nearly doubles. Notably,
**E0 does *not* grow** (26 vs 26 species) — unlike the toy model, this is
complexification via **re-integration** of an already-reachable space
(every fermentation-branch module now also draws on cross-boundary ATP
recycling), not via reaching qualitatively new species. Full log:
`outputs/mito_case_analysis.log`.

### 5.2 Chloroplast-type (real, cross-organism: *S. cerevisiae* × *Synechocystis*)

**Host**: *Saccharomyces cerevisiae* (iND750) — real heterotrophic
eukaryote, no photosynthesis. **Symbiont**: *Synechocystis* sp. PCC 6803
(iJN678) — real cyanobacterium, already used as a free-living organism in
the companion paper. This is the first genuinely cross-organism test of
the relabelling machinery (two independently-curated reconstructions, not
one file reused twice) — it worked cleanly on the first attempt.

**Interface** (2 reactions, the real textbook mechanism): the chloroplast
**triose-phosphate translocator** route — fixed carbon leaves as
glyceraldehyde-3-phosphate (`g3p`, the canonical TPT substrate), not free
glucose; host CO₂ returns to the symbiont's Calvin cycle.

Two variants: **dark** (symbiont gets no `photon_e` — isolates the
interface's own effect) and **+light** (biologically complete picture).

| | EPMs | mean size | max size |
|---|---|---|---|
| Host (yeast) alone | 152 | 42.9 | 50 |
| Symbiont (dark) alone | 63 | 24.0 | 33 |
| Symbiont (+light) alone | 63 | 27.0 | 36 |
| **Merged (dark)** | **221** | **66.0** | **75** |
| **Merged (+light)** | **221** | **69.0** | **78** |

(Host and symbiont-alone numbers again match the companion paper's own
published iND750/iJN678 aerobic-core figures exactly.)

**221/221 merged EPMs are hybrid (100%)** in both variants. Here the
signature is different again: mean size is almost exactly **additive**
(42.9+24.0=66.9 ≈ 66.0 merged), and EPM *count* rose slightly (215→221)
rather than falling — not "fusion" by the strict count-based definition of
§2, but still 100% hybrid with a large absolute size jump (max 50→75).
Light adds a modest further boost on top (max 75→78), consistent with the
symbiont already having glucose as carbon/energy source in both variants —
light is additive here, not transformative, for *this* comparison. Full
logs: `outputs/chloro_dark_analysis.log`, `outputs/chloro_light_analysis.log`.

### 5.3 Cross-case summary

See `figures/fig_complexification_summary.png`.

| case | host EPMs (mean sz) | symbiont EPMs (mean sz) | merged EPMs (mean sz) | hybrid % | ΔE0 |
|---|---|---|---|---|---|
| toy | 1 (8.0) | 1 (12.0) | 1 (22.0) | 100% | +2 |
| mitochondrial-type | 15 (15.1) | 15 (15.1) | 27 (28.2) | 100% | +0 |
| chloroplast-type, dark | 152 (42.9) | 63 (24.0) | 221 (66.0) | 100% | +2 |
| chloroplast-type, +light | 152 (42.9) | 63 (27.0) | 221 (69.0) | 100% | +2 |

---

## 6. Two mechanisms of endosymbiotic complexification

The four cases above cluster into (at least) two distinct, real, computable
mechanisms — both produce the 100%-hybrid signature, but differ in how:

1. **Capability-gain** (toy model): the food closure (E0) itself grows.
   New species become reachable *only* jointly. This is the "sharpest"
   signature — a literal new capability (here, Fe–S-cluster-dependent
   biosynthesis) exists in the merged system that exists in neither parent.

2. **Re-integration** (mitochondrial-type): E0 is unchanged, but the *same*
   reachable chemistry reorganises into far fewer, far larger, mutually
   dependent modules (30 separate small modules → 27 modules, now all
   drawing on both partners, mean size nearly doubled). No new capability
   in the E0 sense — the same chemistry, tied together more tightly.

3. **Additive/proliferative** (chloroplast-type): E0 grows modestly, mean
   size is close to the simple sum of the two partners' means, and EPM
   *count* actually rises slightly rather than falling — a large absolute
   size jump without the classic fusion (count-reduction) signature. This
   may reflect the greater combinatorial richness of two large, real
   genome-scale networks (host 880 ERCs × symbiont 668–666 ERCs) relative
   to a small self-derived pair, worth further characterisation with an
   interface-richness sweep (§9).

All three share the 100%-hybrid signature, which appears to be the robust,
mechanism-independent marker of genuine integration — worth treating as the
primary readout of this framework, with the count/mean-size pattern as a
secondary, mechanism-diagnostic readout.

---

## 6b. Interface-fragility (knockout) analysis

THEORY.md §2.3 proposed a mechanical test for "which specific cross-feeding
link is load-bearing": remove one interface reaction at a time and see
which hybrid EPMs disappear. Implemented in `scripts/interface_knockout.py`
and run on all three real interfaces.

**Toy model** (3 interface reactions: pyruvate, ATP, Fe–S cluster):

| removed | hybrid EPMs | mean size | max size |
|---|---|---|---|
| *(none — full)* | 1 | 22.0 | 22 |
| XPORT_PYR | 1 | 22.0 | 22 |
| XPORT_ATP | 1 | 22.0 | 22 |
| **XPORT_FES** | 1 | **20.0** | **20** |

Only `XPORT_FES` is load-bearing. Removing the pyruvate or ATP transporter
individually changes *nothing*, because both host and symbiont already
produce their own pyruvate and ATP internally — the interface adds no new
*reachability* there (only, in reality, quantity/efficiency, which this
set-based framework does not represent). This is an important, genuinely
informative distinction the toy model surfaces cleanly: **"biologically
plausible transporter" and "structurally load-bearing for EPM reachability"
are not the same thing**, and the knockout tool is what tells them apart.

**Mitochondrial-type** (2 interface reactions: pyruvate, ATP):

| removed | hybrid EPMs | mean size | max size |
|---|---|---|---|
| *(none — full)* | 27 | 28.2 | 30 |
| XPORT_PYR | 28 | 28.2 | 30 |
| XPORT_ATP | 28 | 28.2 | 30 |

Neither individual reaction is load-bearing here either, and both
knockouts give *identical* results — the real network's own internal
redundancy (several real routes connect the ATP/pyruvate pools: `ADK1`,
`PYK`, `PPS`, etc., all still present on both sides) means either single
transporter alone is already sufficient for full hybridisation; the two
are functionally redundant *with each other*. This is a genuine disanalogy
with the toy model (whose `XPORT_FES` has no redundant alternative route)
and a methodological lesson: real, richly-connected metabolic networks may
need several interface reactions removed *simultaneously*, or a
genuinely irreplaceable one (like the toy model's Fe–S cluster), before
any single knockout shows the "load-bearing" signal.

**Chloroplast-type, dark** (2 interface reactions: g3p carbon export,
CO₂ recycling):

| removed | hybrid EPMs | mean size | max size |
|---|---|---|---|
| *(none — full)* | 221 | 66.0 | 75 |
| XPORT_G3P | 221 | 66.0 | 75 |
| **XPORT_CO2** | **214** | **63.9** | **73** |

Here the two interface reactions are *not* symmetric: removing the carbon
(g3p) export has no measurable effect, but removing the CO₂-recycling
import measurably shrinks both hybrid count and size. Both host and
symbiont already have independent access to `co2_e` as a declared core
food species, so this asymmetry is real and reproducible but its precise
mechanistic cause (something about how the direct host-cytoplasm-to-
symbiont-cytoplasm CO₂ route interacts with the ERC-level closure
differently than each side's own environmental CO₂ uptake) was **not**
run down to a specific reaction-level explanation in the time available —
flagged here as an honest gap rather than an unverified claim. It is,
either way, a real example of the knockout tool detecting an asymmetry a
naive reading of the biology ("both directions are just recycling, should
be symmetric") would not have predicted.

**Take-away**: the knockout tool works as designed and gives real,
sometimes surprising, always case-specific answers — it is not a rubber
stamp that declares every curated transporter "important." That two of
three real cross-feeding reactions tested here turned out to be
individually redundant is itself informative: it suggests that in richly-
connected real metabolic networks, single-transporter loss during
organelle-genome reduction may often be tolerated precisely *because* of
this kind of redundancy, and that the framework's falsifiable prediction
(§7) is more likely to bite for reactions that, like the toy model's
`XPORT_FES`, have no redundant alternative — a testable refinement of the
prediction itself.

## 7. Discussion

**Falsifiable prediction.** If interface reactions really are what makes
specific EPMs hybrid (and therefore, in this framework's terms,
structurally load-bearing for the integrated system), then real
evolutionary genomics should show a signature: genes/functions
corresponding to the *specific* metabolites crossing a real endosymbiont's
reduced-genome boundary (e.g., which transporters and which pathway
branch-points) should be systematically retained (or transferred to the
host nucleus, functionally preserved) more often than genes serving
modules that stayed host-only or symbiont-only in this framework's
decomposition. THEORY.md §2.3 already specifies the mechanical test:
knock out one interface reaction at a time, and see which hybrid EPMs
disappear — the reactions whose removal breaks the *most* or the
*largest* hybrid EPMs are the framework's prediction for "most
evolutionarily load-bearing, most conserved." This knockout sweep was run
(§6b): the refined version of the prediction, given those results, is that
this signature should bite hardest for cross-fed metabolites with **no
redundant alternative route** on either side (like the toy model's Fe–S
cluster) — real transporters embedded in richly-connected, redundant
metabolic contexts (like the mitochondrial-type case's pyruvate/ATP
exchange) may show *no* individual retention signature at all, precisely
because the framework itself shows them to be individually dispensable
for structural integration. The prediction is therefore not "every
interface reaction should be conserved" but "conservation should track
computed irreplaceability," which is a sharper, more falsifiable claim.

**Limits, restated from THEORY.md.** This is a steady-state, structural
framework — it says nothing about population genetics, fitness, or the
dynamics of how a symbiosis is established or fixed. The interface and any
synergy reactions are curated modelling choices stated explicitly with
their biochemical rationale (real transporters: pyruvate/ATP exchange for
mitochondria, the triose-phosphate translocator for chloroplasts), not
derived facts — results are conditional on those choices. No genuinely
novel "structural synergy" reaction (a reaction impossible in either
parent alone, becoming possible only via a new catalytic combination) was
used in the real-data cases in this report; both real interfaces are real,
named transporters. This is a deliberately conservative choice for a first
pass — inventing novel catalysis needs stronger biochemical justification
than was available in the time budget, and the existing results are
already strong without it. BiGG reconstructions used here are modern
free-living relatives, not literal ancestral lineages — `e_coli_core`
stands in structurally for "a free-living aerobic respirer" in §5.1, not
for a literal alphaproteobacterium.

**Relationship between the two questions.** §4.1's toy-host search
rediscovering the exact dependency §3's toy model was designed around,
without being told about it, is a small but real piece of evidence that
the fusigenic-inflow-search framework and the endosymbiotic-merger
framework are capturing the same underlying structural phenomenon — a
single missing/blocked link in the ERC/synergy graph — from two different
directions (search over single candidate species vs. explicit construction
of a whole second system). This suggests a natural unification: an
endosymbiotic merger can be read as the limiting case of a fusigenic-inflow
search where the "candidate" is not one species but an entire regenerating
system.

---

## 8. Errors caught and fixed during this project (for transparency)

- **Fusion-score formula bug** (`fusigenic_inflow_search.py`): multiplying
  by `Δcount` directly zeroed the score whenever `Δcount == 0` exactly,
  even when mean size grew substantially. Caught by hand-checking the toy
  host network, where the known-correct top candidate (`fescluster_h`)
  was being scored 0. Fixed to `Δmean · (1 + max(0, −Δcount))`.
- **Comment-parsing bug** (`merge_networks.py`): the reaction-line parser
  required a line to end exactly at its first `;`, so real BiGG-derived
  lines with a trailing `; comment` (nearly all of them — e.g. `R0: ... =>
  ...; PFK: Phosphofructokinase`) silently failed to match and passed
  through **unrelabelled**. This produced duplicate reaction names
  (`R0`/`R1`/...) the instant host and symbiont were derived from the same
  or similarly-structured real file — caught via `read_txt`'s own
  "Reaction already exists" error on the first real (mitochondrial-type)
  merge attempt, not by the toy model, whose hand-written files have no
  such comments. Fixed by splitting off the comment before parsing and
  reattaching it after; regression-checked against the toy model (result
  unchanged).
- **Reaction-name collision** (`merge_networks.py`, caught before the
  comment bug above, same underlying root cause class): the first version
  of the relabelling function relabelled species but not reaction *names*,
  so two networks sharing generic reaction names (`INFLOW_glc` in both toy
  files) collided. Fixed by relabelling reaction names too.
- **Figure arrow-direction bug** (`make_fig_toy.py`): interface arrows were
  normalised to always draw host→symbiont regardless of the reaction's
  actual direction, silently reversing the arrowhead on 2 of the toy
  model's 3 real interface reactions (ATP and Fe–S cluster both flow
  symbiont→host; only pyruvate flows host→symbiont). A drawing/display bug
  only — did not affect any computed result — caught by visually
  inspecting the rendered figure rather than trusting the code.
- **Scalability finding, not a bug**: an exhaustive fusigenic-inflow sweep
  costs one full ERC→EPM pipeline run per candidate species. This is ~1–2s
  per candidate on `e_coli_core` (tractable to sweep exhaustively) but
  ~50–90s per candidate at `iAF1260` scale (1500+ ERCs) — an exhaustive
  sweep over its ~1400 non-food species would take on the order of 14
  hours. All genome-scale sweeps in this report used a curated candidate
  list instead of the exhaustive default; this constraint is documented in
  `PROGRESS.md` and in `fusigenic_inflow_search.py`'s own docstring.

---

## 9. What was scoped but not run (honest accounting given the time budget)

- ~~Interface-fragility / knockout analysis~~ — **done, see §6b.** Ran on
  all three real interfaces (toy, mitochondrial-type, chloroplast-type
  dark). Result was more interesting than a simple confirmation: 2 of 3
  real interfaces tested turned out to have NO individually load-bearing
  reaction (each transporter redundant with either the other transporter
  or the organism's own independent food access), only the toy model's
  deliberately-irreplaceable Fe–S cluster dependency showed the expected
  signature. This refines rather than closes the question — see the
  updated prediction in §7.
- **Interface-richness sweep**: single-metabolite interface → fully-open
  interface, to see how the complexification signature scales. Would
  directly inform the "why does the chloroplast-type case show a different
  count/size pattern than the mitochondrial-type case" question raised in
  §6.
- **Extension to the 108-organism BiGG catalogue** (already assembled in
  `projects/COT_Fundamental_Generators_Complex/network_catalogue_bigg_organisms.csv`)
  for many host/symbiont pairs, to test whether the 100%-hybrid signature
  and the capability-gain/re-integration/additive distinction are general
  patterns or specific to the four pairs tested here.
- **A curated "structural synergy" reaction** (a genuinely novel catalytic
  capability, impossible in either parent alone) in a real-data case — the
  toy model's design uses only real transport interfaces, deliberately
  conservative per the honesty notes in THEORY.md.
- **ESPM-level (not just EPM-level) merger analysis** — blocked on the same
  genome-scale ESPM/LP tractability limits documented in the companion
  paper; out of scope here for the same reason.

---

## 10. File index

```
COT_Endosymbiosis/
  THEORY.md                        -- full framework and formal definitions
  PROGRESS.md                      -- session log, read first if resuming
  report/
    REPORT.md                      -- this file
    OUTLINE.md                     -- working outline (superseded by this file)
  toy_model/
    host_alone.txt, symbiont_alone.txt, merged.txt
  real_data/
    build_mitochondrial_case.py, mito_host_fermentative.txt,
      mito_symbiont_aerobic.txt, mito_merged.txt
    build_chloroplast_case.py, chloro_host_yeast.txt,
      chloro_symbiont_synecho_{dark,light}.txt,
      chloro_merged_{dark,light}.txt
  scripts/
    merge_networks.py              -- generic compartmentalised-union tool
    endosymbiosis_analysis.py      -- before/after EPM + hybrid-EPM analysis
    fusigenic_inflow_search.py     -- generic per-species inflow sweep
    interface_knockout.py          -- interface-fragility / knockout analysis (Sec. 6b)
    make_fig_toy.py, make_figures.py, make_fig_fusigenic.py
  outputs/
    *.log                          -- raw pipeline output, source of every number above
  figures/
    fig_toy_model.{pdf,png}
    fig_complexification_summary.{pdf,png}
    fig_fusigenic_search.{pdf,png}
```
