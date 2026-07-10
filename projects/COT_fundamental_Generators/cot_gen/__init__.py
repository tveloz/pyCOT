"""
cot_gen — Efficient generative computation of persistent modules (COT).

Implements the theory from Veloz & Bassi (2025) "Synergy and Complementarity:
The Generative Basis of Chemical Organizations".

Build order:
  Stage 0: Preprocessing & E∅  (types, io, metrics)
  Stage 1: Closure + ERCs + MinBas  (closure, erc)
  Stage 2: ERC Hierarchy as query index  (hierarchy)
  Stage 3: Fundamental synergies  (synergy)
  Stage 4: Fundamental complementarities  (complementarity)
  Stage 5: Generators → EPMs → ESPMs  (generators)
"""
__version__ = "0.1.0"
