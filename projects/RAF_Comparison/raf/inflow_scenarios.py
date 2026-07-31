"""
inflow_scenarios.py -- build a variant of a network's own .txt encoding
with a CHOSEN inflow (food) set, so the same network can be re-analyzed
under several inflow scenarios without touching cot_gen internals.

`cot_gen.io_pyCOT.build_rndata` derives the food/E0 set purely from
whichever reactions in the given network have empty support ("=> s");
there is no separate F0 parameter to override. So to vary the inflow
scenario we build a genuinely different network: the SAME real (non-
inflow) reactions, with a different set of native inflow reactions
swapped in. This is done at the TXT-text level (strip existing "=> s"
lines, append new ones for the chosen food species) rather than via the
ReactionNetwork object API, reusing the exact grammar already validated
in tools/sbml_reconstruct.py's round-trip checks.

Outflow is left untouched in every scenario -- this suite varies only
inflow (see scripts/run_suite.py for the reasoning on scope).
"""
from __future__ import annotations


def strip_native_inflows(txt: str) -> list[str]:
    """Return every reaction line of `txt` EXCEPT its native inflow
    reactions (empty left-hand side), unchanged otherwise (including any
    outflow reactions, which this suite does not vary)."""
    kept = []
    for line in txt.splitlines():
        line = line.rstrip()
        if not line.strip():
            continue
        body = line.split(";", 1)[0]
        if ":" not in body or "=>" not in body:
            continue
        lhs = body.split(":", 1)[1].split("=>", 1)[0].strip()
        if lhs == "":
            continue  # drop native inflow
        kept.append(line)
    return kept


def build_scenario_txt(base_txt: str, food_tokens: list[str],
                        outflow_tokens: list[str] | None = None) -> str:
    """`food_tokens`/`outflow_tokens` must be the EXACT species tokens as
    they already appear elsewhere in `base_txt` (e.g. "glc__D_e" for
    BiGG-bare convention, or "C (Cyclin)" for BIOMD-decorated convention)
    -- reusing an existing token guarantees it resolves to the same
    species pyCOT already knows about, rather than accidentally minting a
    new one.

    `outflow_tokens` implements Def:translation's Omega (drained species):
    an outflow reaction s->0 is added for each one, ON TOP OF whatever
    outflow reactions the base network already has (real BiGG exchange
    reactions with a one-way "s =>" form are left untouched, matching the
    "vary inflow/outflow as an explicit, separate dial" reading -- this
    function only ever ADDS outflows, never strips existing ones, unlike
    the inflow side which replaces the network's native food)."""
    lines = strip_native_inflows(base_txt)
    for i, tok in enumerate(food_tokens):
        lines.append(f"INFLOW_{i}:  => 1 {tok};")
    for i, tok in enumerate(outflow_tokens or []):
        lines.append(f"OUTFLOW_{i}: 1 {tok} => ;")
    return "\n".join(lines) + "\n"
