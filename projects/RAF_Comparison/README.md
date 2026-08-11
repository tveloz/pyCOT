# RAF_Comparison

Code and worked examples supporting **"Decomposing Autocatalytic Sets through
Chemical Organization Theory"** (T. Veloz). Two things live here:

1. **Section 4's two computational illustrations** — a synthetic binary
   heteropolymer network (§4.1) and the core *E. coli* metabolic
   reconstruction analysed with no catalyst assumption at all (§4.2) — with
   the exact scripts that produced every figure and table in that section.
2. **Validation tests** for the worked examples of Sections 2–3 (the toy
   $\{x,y\}$ semi-organization-but-not-organization example, the running
   network $\{f,a,b,g,c,d,h\}$ used for both the relational-bridge theorem
   and the food-dependency DAG), checking the code against the paper's own
   hand-derivable numbers rather than assuming the implementation is correct.

## This folder does not stand alone

`RAF_Comparison` is one project inside the `pyCOT` monorepo, and it is not
self-contained: it implements the **RAF layer** (catalytic reaction systems,
maxRAF, catalyst induction) and reuses two sibling projects for everything
else it needs:

| What | Where | Used for |
|---|---|---|
| Reaction-network I/O (`pyCOT.io`) | `src/pyCOT/` | reading the `.txt` network files |
| ERC / hierarchy / synergy / EPM / ESPM engine (`cot_gen`) | `projects/COT_Fundamental_Generators_Exploration/` | enumerating the full semi-organization lattice (Veloz, "Computing chemical organizations efficiently using a minimal generative structure", submitted 2026) |
| E/F/fragile-circuit decomposition + Hasse-diagram plotting (`decomp`) | `projects/Decomposition_Theorem/` | the paper's decomposition theorem and every Hasse-diagram figure |

Scripts here add all three sibling folders to `sys.path` at import time,
relative to the repository root (see the top of `scripts/run_full_analysis.py`).
**Clone the whole `pyCOT` repository and run these scripts from within it** —
copying just this folder elsewhere will not work.

## Setup

From the repository root:

```bash
poetry install          # or: pip install numpy networkx matplotlib pandas bitarray scipy rustworkx
```

Python ≥3.11 (see the root `pyproject.toml` for exact package versions).
No credentials or network access needed — the one external data file this
project depends on (`bigg_e_coli_core.txt`, the BiGG core *E. coli*
reconstruction) is vendored under `data/bigg/`, not read from the shared
`data/` tree elsewhere in the repo.

## Layout

```
raf/                          RAF engine (this project's own code)
  crs.py                        CRS / Reaction types, gen()
  raf_algo.py                   is_raf, compute_maxRAF, close_raf, irrRAF search
  biomodel_crs.py                three catalyst-induction modes (see below)
  cofactor_pools.py              BiGG cofactor-pool naming-convention catalysis
  heteropolymer.py               binary heteropolymer network generator (§4.1)
  cot_bridge.py / decomp_shim.py / dependency.py
                                 the Omega-outflow COT translation and the
                                 food-dependency DAG, bridging into decomp
  inflow_scenarios.py            rewrite a network's food/outflow set

scripts/
  generate_heteropolymer_network.py   builds data/heteropolymer/heteropolymer_L3_uniform_cleave0.5.txt
  run_full_analysis.py                RAF layer + full semi-organization lattice + figures,
                                       for ONE network selected by NETWORK_CHOICE at the top

tests/                         validate against the paper's own worked examples
  validate_raf_algo.py           the paper's toy networks N/N' (RAF/maxRAF/irrRAF)
  validate_bridge.py              the {f,a,b,g,c,d,h} running example (relational-bridge theorem)
  validate_dependency.py          the same running example's dependency DAG (depth, indecomposability)

data/
  bigg/bigg_e_coli_core.txt      vendored BiGG core E. coli reconstruction (§4.2)
  heteropolymer/*.txt            generated network file for §4.1 (regenerable)

figures/                       exactly the three images the paper's LaTeX includes
```

## Catalyst induction: three modes, why three

RAF theory needs a catalyst assignment $\mathcal C: R \to 2^M$; most real
network sources don't carry one. `raf/biomodel_crs.py` offers three ways to
get one, used in different places in the paper:

- **`net_zero`** — a species is a catalyst of a reaction if it appears with
  equal, cancelling coefficient on both sides of *that one reaction*. This is
  how the heteropolymer network's catalysts (written directly into its `.txt`
  file) are recovered — §4.1 uses this mode.
- **`cofactor_pools`** — matches BiGG naming conventions (`nad_c`/`nadh_c`,
  `atp_c`/`adp_c`, ...) to flag cofactor-pool interconversions as catalysed.
  Not used for any number reported in the final paper, but kept because it's
  a real, documented alternative for BiGG-style networks where `net_zero`
  structurally finds nothing (see the module docstring).
- **`none`** — no induction at all, $\mathcal C(r) = \emptyset$ for every
  $r$. This is what §4.2 uses for *E. coli*: a flux-balance metabolite
  reconstruction genuinely carries no per-reaction catalyst information, so
  taking that at face value (rather than approximating one) is the point of
  that subsection, not a simplification.

## Reproducing the paper's figures and tables

Everything in Section 4 comes from `scripts/run_full_analysis.py`. Edit
`NETWORK_CHOICE` at the top of the file to `"heteropolymer"` or
`"ecoli_none"` and run it:

```bash
python projects/RAF_Comparison/scripts/run_full_analysis.py
```

| Paper reference | `NETWORK_CHOICE` | Output |
|---|---|---|
| §4.1, `\Cref{fig:polymer-hasse}` | `"heteropolymer"` | `figures/heteropolymer_L3_uniform_cleave0.5_organization_hasse.png` |
| §4.1, `\Cref{fig:polymer-raf-per-org}` | `"heteropolymer"` | `figures/heteropolymer_L3_uniform_cleave0.5_raf_vs_organization_per_org.png` |
| §4.1, `\Cref{tab:polymer-raf-per-org}` | `"heteropolymer"` | printed to stdout as "Per-organization RAF vs COT comparison" |
| §4.2, `\Cref{fig:ecoli-clean-breakdown}` | `"ecoli_none"` | `figures/bigg_e_coli_core_organization_hasse.png` |
| §4.2, `\Cref{tab:ecoli-clean-per-org}` | `"ecoli_none"` | printed to stdout as "Per-organization RAF vs COT comparison", including `|E|`, `|F|` and circuit sizes; only the paper table's free-text "biochemical identity" column (e.g. "hexose-phosphate pair (F6P/G6P)") is added by hand, by reading each circuit's species names off `r.circuits[i].species_mask` via `rn_data.bitset_to_names(...)` |

The heteropolymer network file is already committed under
`data/heteropolymer/`, but if you want to regenerate it (e.g. after editing
`raf/heteropolymer.py`), its exact parameters — `max_length=3`,
`catalysis_mode="uniform"`, `p_catalyst=0.2`, `cleavage_prob=0.5`, `seed=3`
— are locked into `scripts/generate_heteropolymer_network.py`'s
configuration block, with the reasoning for `cleavage_prob=0.5` (why full
reversibility trivializes the semi-organization/organization distinction)
documented right above it:

```bash
python projects/RAF_Comparison/scripts/generate_heteropolymer_network.py
```

Both scripts print a `[RAF layer]` block (catalyst counts, global maxRAF)
and a `[COT layer]` block (ERC/semi-organization/organization counts) before
writing figures — the run is not silent, so a wrong number surfaces
immediately rather than only in a plot nobody re-checks.

## Running the validation tests

```bash
python projects/RAF_Comparison/tests/validate_raf_algo.py
python projects/RAF_Comparison/tests/validate_bridge.py
python projects/RAF_Comparison/tests/validate_dependency.py
```

Each prints `[OK]` per check and `ALL CHECKS PASSED` at the end; each checks
the code against a number derivable by hand from the paper's own examples,
not against another run of the same code.
