"""
gen_v2 -- experimental SO-search engine, kept fully separate from the
trusted pyCOT.analysis.organizations.so_search.

Built from the discussion with Tomas (session of 2026-10-03/04) about:
  1. The generator-vs-closure gap: the production search decides "what to
     try next" (for synergy) by asking only the ERCs EXPLICITLY listed in
     the current generator, not everything the generator's closure
     actually covers. Confirmed to be a real asymmetry in the code (traced
     to exact lines); an initial A/B test across 17 synergy-bearing real
     networks (<=100 reactions) showed no observed difference, which in
     hindsight was a corpus-size artifact, not evidence of safety -- see
     the next paragraph.

     gen_v2's fix (ClosureCompleteGraph, engine.py) propagates each ERC's
     registered synergy partnerships up to its hierarchy ancestors. The
     FIRST version of this fix was itself buggy: it didn't check that the
     partner stays incomparable to the ancestor it's copied to (synergy is
     only defined between incomparable ERCs). Caught on a LARGER network
     (BIOMD0000000109, 19 ERCs, 80-800-reaction sweep of 2026-10-04): the
     unfiltered version caused gen_v2 to silently drop 2 genuine
     semi-organizations the production engine correctly finds (confirmed
     via a from-scratch closed+SSM check against the raw reaction data,
     bypassing both engines). Fixed by requiring the inherited partner to
     remain incomparable to its new ancestor. Lesson: the earlier
     "no difference on 17 networks" result was real but not sufficient --
     it was silent about exactly this failure mode because none of those
     17 networks happened to have the right hierarchy shape to trigger it.
     On the same 80-800-reaction sweep, post-fix, the corrected version
     recovered hundreds of genuine previously-missed SOs on several
     networks (e.g. +742 on BIOMD0000000943) with zero soundness failures
     on independent verification.
  2. Irreducible generators: a vertical-lift step always makes its source
     ERC redundant (Lemma 5, proof is immediate from hierarchy
     monotonicity) -- gen_v2 drops it the moment the lift fires, cheaply
     and exactly, rather than leaving dead weight in the generator.
  3. Ordered construction history: the paper's own Definition 29
     ("fundamental generator") is an ORDERED partition of extension steps,
     not a flat set. The production engine computes the right final
     species sets but never retains the order/reasons -- gen_v2 keeps a
     full step-by-step log (seed / complementarity / synergy / lift, and
     what each step made redundant) alongside the existing bookkeeping.

See engine.py for the implementation and validate_small.py for the
correctness check against the SAME oracle-validated small-network corpus
so_search.py is checked against, before this is trusted on anything larger.
compare_old_new.py runs the old (production) and new (this) engines side
by side on real, large (200-1000 reaction) networks, where no brute-force
oracle is feasible -- old vs new is the only available cross-check there.
"""
