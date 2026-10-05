"""
merge_networks.py -- generic compartmentalised-union endosymbiotic merge tool.

Implements the "Endosymbiotic merger" construction from THEORY.md Section 2.2:

    M = H (unchanged) (+) S_relabelled (+) Interface (+) Synergy

Given a host network file and a symbiont network file (both in the simple
"Name: lhs => rhs;" text format read by pyCOT.io.functions.read_txt), this
module:

  1. Parses each reaction line into (name, [(coeff, species), ...] lhs,
     [(coeff, species), ...] rhs).
  2. Relabels every symbiont species by appending a tag (default "__endo"),
     so host and symbiont internal chemistries stay physically distinct even
     where BiGG-style species tokens coincide (e.g. both use "atp_c") --
     see THEORY.md for why this matters: pre-fusion, host-cytoplasm ATP and
     symbiont-cytoplasm ATP are different physical pools.
  3. Concatenates host reactions (unchanged) + relabelled symbiont reactions
     + explicit interface reactions (transport/exchange across the new
     compartment boundary) + explicit synergy reactions (novel catalytic
     capabilities unlocked only by co-localisation -- curated, not derived).
  4. Writes the merged network to a new .txt file.

Interface and synergy reactions are given as plain reaction-line strings
using the ALREADY-RELABELLED symbiont species names (e.g.
"XPORT_PYR: 1 pyr_h => 1 pyr_s__endo;"), so the caller has full, explicit
control over exactly what crosses the boundary -- this is deliberately not
automated, per THEORY.md's honesty notes: the interface is a modelling
choice, not a derived fact.

Usage
-----
    from merge_networks import merge_networks, relabel_species_in_text

    merge_networks(
        host_path="host_alone.txt",
        symbiont_path="symbiont_alone.txt",
        interface_reactions=[
            "XPORT_PYR: 1 pyr_h => 1 pyr_s__endo;",
            "XPORT_ATP: 1 atp_s__endo => 1 atp_h;",
            "XPORT_FES: 1 fescluster_s__endo => 1 fescluster_h;",
        ],
        out_path="merged.txt",
    )
"""
from __future__ import annotations

import os
import re

_LINE_RE = re.compile(r"^\s*([^:;]+):\s*(.*?)\s*$")
_TERM_RE = re.compile(r"^\s*(?:(\d+(?:\.\d+)?)\s+)?(\S+)\s*$")


def parse_reaction_line(line: str):
    """Parse 'Name: lhs => rhs; optional trailing comment' into
    (name, lhs_terms, rhs_terms, comment). Real BiGG-derived files (e.g.
    e_coli_core) routinely have a human-readable comment after the
    semicolon that closes the reaction (e.g. "R0: ... => ...; PFK:
    Phosphofructokinase") -- an earlier version of this regex required the
    line to end exactly at the first ';', so those lines silently failed to
    match and were passed through UNRELABELLED, which produced duplicate
    "R0"/"R1"/... reaction names (e_coli_core's own generic per-line names)
    the moment host and symbiont were both derived from e_coli_core --
    caught via read_txt's own "Reaction already exists" error on the real
    mitochondrial-type merge, not by the toy model (whose hand-written
    files have no trailing comments, so this path was never exercised).

    Each side is a list of (coeff:str, species:str) tuples. Returns None for
    blank/unparseable lines (passed through unchanged by callers)."""
    if ";" not in line:
        return None
    body_part, comment = line.split(";", 1)
    m = _LINE_RE.match(body_part)
    if not m:
        return None
    name, body = m.group(1).strip(), m.group(2)
    if "=>" not in body:
        return None
    lhs_str, rhs_str = body.split("=>", 1)

    def _terms(side: str):
        side = side.strip()
        if not side:
            return []
        out = []
        for term in side.split("+"):
            tm = _TERM_RE.match(term)
            if not tm:
                continue
            coeff = tm.group(1) or "1"
            species = tm.group(2)
            out.append((coeff, species))
        return out

    return name, _terms(lhs_str), _terms(rhs_str), comment


def format_reaction_line(name: str, lhs, rhs, comment: str | None = None) -> str:
    def _fmt(terms):
        return " + ".join(f"{c} {s}" for c, s in terms)
    line = f"{name}: {_fmt(lhs)} => {_fmt(rhs)};"
    if comment:
        line += comment
    return line


def relabel_species_in_text(text: str, tag: str, exclude: set[str] | None = None) -> str:
    """Append `tag` to every species token AND every reaction name in every
    line of `text` (reaction names also need relabelling: two independently-
    authored networks routinely reuse generic names like "INFLOW_glc" or, at
    genome scale, real BiGG reaction IDs -- read_txt raises on a duplicate
    reaction name, so an unrelabelled collision is a hard parse error, not
    just a modelling nicety).

    `exclude` species are left unrelabelled (use for species that should
    stay shared/identified across host and symbiont namespaces, if any --
    normally empty; prefer explicit interface reactions instead)."""
    exclude = exclude or set()
    out_lines = []
    for line in text.splitlines():
        parsed = parse_reaction_line(line)
        if parsed is None:
            out_lines.append(line)
            continue
        name, lhs, rhs, comment = parsed

        def _relabel(terms):
            return [(c, s if s in exclude else f"{s}{tag}") for c, s in terms]

        out_lines.append(format_reaction_line(f"{name}{tag}", _relabel(lhs), _relabel(rhs), comment))
    return "\n".join(out_lines) + "\n"


def merge_networks(
    host_path: str,
    symbiont_path: str,
    interface_reactions: list[str],
    out_path: str,
    *,
    synergy_reactions: list[str] | None = None,
    symbiont_tag: str = "__endo",
    exclude_from_relabel: set[str] | None = None,
) -> str:
    """Build the compartmentalised-union merged network and write it to
    `out_path`. Returns `out_path`."""
    with open(host_path, "r", encoding="utf-8") as f:
        host_txt = f.read()
    with open(symbiont_path, "r", encoding="utf-8") as f:
        symbiont_txt = f.read()

    symbiont_relabelled = relabel_species_in_text(
        symbiont_txt, symbiont_tag, exclude=exclude_from_relabel
    )

    parts = [
        f"# --- host (from {os.path.basename(host_path)}) ---",
        host_txt.rstrip("\n"),
        f"# --- symbiont (from {os.path.basename(symbiont_path)}, "
        f"relabelled with tag '{symbiont_tag}') ---",
        symbiont_relabelled.rstrip("\n"),
        "# --- interface (curated transport/exchange reactions) ---",
        "\n".join(interface_reactions),
    ]
    if synergy_reactions:
        parts.append("# --- structural synergy (curated, hypothesis-flagged) ---")
        parts.append("\n".join(synergy_reactions))

    merged_txt = "\n".join(parts) + "\n"

    out_dir = os.path.dirname(out_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(merged_txt)
    return out_path


if __name__ == "__main__":
    # Smoke test on the toy model.
    _here = os.path.dirname(os.path.abspath(__file__))
    _toy = os.path.join(_here, "..", "toy_model")
    out = merge_networks(
        host_path=os.path.join(_toy, "host_alone.txt"),
        symbiont_path=os.path.join(_toy, "symbiont_alone.txt"),
        interface_reactions=[
            "XPORT_PYR: 1 pyr_h => 1 pyr_s__endo;",
            "XPORT_ATP: 1 atp_s__endo => 1 atp_h;",
            "XPORT_FES: 1 fescluster_s__endo => 1 fescluster_h;",
        ],
        out_path=os.path.join(_toy, "merged.txt"),
    )
    print(f"wrote {out}")
    with open(out, encoding="utf-8") as f:
        print(f.read())
