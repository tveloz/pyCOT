"""
utils_ercs.py
=============
Shared ERC pkl-cache loader for the Generative_Structure_Orgs scripts.

Public API
----------
load_ercs(txt_path, RN, ERC_class) -> (ercs, from_cache)

    Load ERCs from <txt_path>.pkl if available; otherwise compute via
    ERC_class.ERCs(RN) and write a fresh cache.

    The returned ERC objects have _closure_names pre-populated so that
    get_closure_names(RN) short-circuits without needing the actual RN
    species-object graph.  All synergy / complementarity / hierarchy
    operations that only use closure_names and min_generator species
    names work correctly with these objects.

    Old-format pkl files (list-of-4 records) are silently upgraded to
    the current dict format on first read.

    NOTE: The E_∅ filter (ercs with empty closure) is NOT applied here.
    Each calling script applies it according to its own semantics.
"""
import os
import pickle


# ---------------------------------------------------------------------------
# Minimal species stub — carries only .name, which is all that
# species_list_to_names() and ERC_Hierarchy.build_hierarchy_graph() need.
# ---------------------------------------------------------------------------
class _S:
    __slots__ = ('name',)

    def __init__(self, name):
        self.name = name


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _save_pkl(txt_path, ercs, RN):
    """Write the cache in the current dict format."""
    pkl_path = txt_path.replace('.txt', '.pkl')
    raw = []
    for e in ercs:
        raw.append({
            'label': e.label,
            'closure_names': sorted(e.get_closure_names(RN)),
            'min_generator_names': [
                sorted([sp.name for sp in gen])
                for gen in e.min_generators
            ],
        })
    try:
        with open(pkl_path, 'wb') as f:
            pickle.dump({'ERCs': raw}, f)
    except Exception:
        pass


def _from_raw(raw, ERC_class):
    """
    Reconstruct ERC objects from cached records.

    Handles two formats:
      current  — list of dicts: {'label', 'closure_names', 'min_generator_names'}
      legacy   — list of lists: [min_generators_pickled, closure_names_list, [], label]
    """
    ercs = []
    needs_upgrade = False
    for rec in raw:
        if isinstance(rec, dict):
            label         = rec['label']
            closure_names = set(rec['closure_names'])
            gen_objs      = [[_S(n) for n in gen] for gen in rec['min_generator_names']]
        else:
            # Legacy format written by older versions of the stats scripts
            min_gens, cl_names, _, label = rec
            closure_names = set(cl_names)
            # Pickled species objects still carry .name — reuse them directly
            gen_objs      = min_gens
            needs_upgrade = True

        e = ERC_class(min_generators=gen_objs, label=label)
        e._closure_names = closure_names
        ercs.append(e)

    return ercs, needs_upgrade


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def load_ercs(txt_path, RN, ERC_class):
    """
    Return ``(ercs, from_cache)``.

    Parameters
    ----------
    txt_path  : str  — path to the .txt reaction-network file
    RN        : ReactionNetwork  — loaded network (needed if recomputing)
    ERC_class : the ERC class (passed in to avoid a circular import)

    Returns
    -------
    ercs       : list of ERC objects with _closure_names pre-populated
    from_cache : bool — True if loaded from pkl, False if freshly computed
    """
    pkl_path = txt_path.replace('.txt', '.pkl')

    if os.path.exists(pkl_path):
        try:
            with open(pkl_path, 'rb') as f:
                cache = pickle.load(f)
            raw = cache.get('ERCs', [])
            if raw:
                ercs, needs_upgrade = _from_raw(raw, ERC_class)
                if needs_upgrade:
                    # Silently rewrite in current format
                    _save_pkl(txt_path, ercs, RN)
                return ercs, True
        except Exception:
            pass  # corrupt / incompatible cache — fall through to recompute

    ercs = ERC_class.ERCs(RN)
    _save_pkl(txt_path, ercs, RN)
    return ercs, False
