#!/usr/bin/env python3
"""
bigg_to_pycot.py
================================================================
Downloads models from the BiGG database and converts them to
pyCOT .txt format.

Usage (run from any directory):
    # List all BiGG models with sizes (no download):
    python bigg_to_pycot.py --list

    # Download and convert ALL BiGG models:
    python bigg_to_pycot.py --all

    # Download specific models by BiGG ID:
    python bigg_to_pycot.py iJO1366 iML1515 iMM904

    # Convert a locally-saved BiGG JSON file:
    python bigg_to_pycot.py --local path/to/model.json

    # Download only models in a reaction-count range:
    python bigg_to_pycot.py --all --min-rxns 50 --max-rxns 500

Output directory:
    data/biomodels/biomodels_all_txt/bigg_{model_id}.txt

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
PARSER SAFETY NOTE — digit-starting BiGG metabolite IDs
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
The pyCOT parser uses a compact-notation regex to strip leading
integer coefficients from species tokens.  BiGG metabolite IDs
that start with a digit (e.g. 3pg_c, 2pg_c, 13dpg_c) would be
misinterpreted when written bare:

    BAD:  3pg_c + atp_c => ...   →  parser sees coef=3, species=pg_c
    OK:   1 3pg_c + 1 atp_c => ...  →  parser sees coef=1, species=3pg_c

This converter therefore ALWAYS writes explicit stoichiometric
coefficients in space-separated form, even when the coefficient
is 1.  The pyCOT space-match rule handles all cases correctly.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
REVERSIBLE REACTIONS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
BiGG encodes reversibility via flux bounds.  Any reaction with
lower_bound < 0 is split into two irreversible reactions:
    R{i}_fwd:  lhs => rhs
    R{i}_rev:  rhs => lhs
Exchange reactions (one metabolite, id starts with EX_) that
are bidirectional are split into import and export lines.
"""

import json
import os
import sys
import time
import argparse
import urllib.request
import urllib.error

# ── Paths ────────────────────────────────────────────────────────────────────
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_PYCOT_ROOT = os.path.normpath(os.path.join(_SCRIPT_DIR, '..', '..', '..'))

# Use HTTP instead of HTTPS due to firewall blocking port 443 for Python applications
BIGG_API_BASE = 'http://bigg.ucsd.edu/api/v2'

DEFAULT_OUT_DIR = os.path.join(
    _PYCOT_ROOT, 'data', 'biomodels', 'biomodels_all_txt'
)

# ── HTTP helpers ─────────────────────────────────────────────────────────────
def _get_json(url, timeout=90, retries=3, retry_delay=5):
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                url,
                headers={'User-Agent': 'pyCOT-bigg-converter/1.0'}
            )
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return json.loads(resp.read().decode('utf-8'))
        except urllib.error.URLError as e:
            if attempt < retries - 1:
                print(f'    [retry {attempt+1}] {e}', flush=True)
                time.sleep(retry_delay)
            else:
                raise


def list_bigg_models():
    """Return list of model metadata dicts, sorted by reaction count."""
    data = _get_json(f'{BIGG_API_BASE}/models')
    return sorted(data['results'], key=lambda m: m.get('reaction_count', 0))


def download_bigg_json(model_id):
    """Download and return a single BiGG model as a JSON dict."""
    url = f'{BIGG_API_BASE}/models/{model_id}/download?model_format=JSON'
    return _get_json(url, timeout=180)


# ── Conversion ───────────────────────────────────────────────────────────────
def _coef_str(c):
    """Format a numeric stoichiometry as a clean string."""
    if float(c) == int(float(c)):
        return str(int(float(c)))
    return str(c)


def _fmt_term(coef, met_id):
    """
    Format one stoichiometric term ALWAYS with an explicit coefficient.

    Using '1 met_id' (space-separated) even for coef=1 prevents the pyCOT
    compact-notation regex from incorrectly stripping leading digits from
    BiGG IDs like '3pg_c', '2pg_c', '13dpg_c'.
    """
    return f'{_coef_str(coef)} {met_id}'


def _fmt_side(species_list):
    """
    Format one reaction side from a list of (coef, met_id) tuples.
    Returns '' for an empty side (source/sink reactions).
    """
    if not species_list:
        return ''
    return ' + '.join(_fmt_term(c, m) for c, m in species_list)


def convert_model(model_data):
    """
    Convert a BiGG JSON model dict to pyCOT .txt reaction lines.

    Parameters
    ----------
    model_data : dict
        Parsed BiGG JSON model.

    Returns
    -------
    lines : list[str]
        One pyCOT reaction string per entry, ending with ';'.
    stats : dict
        Summary counts.
    """
    lines = []
    idx = 0          # running reaction index for unique Rn names
    n_split = 0      # number of reversible splits

    for rxn in model_data.get('reactions', []):
        rxn_id   = rxn.get('id', f'R{idx}')
        mets     = rxn.get('metabolites', {})   # met_id → stoich (neg=reactant)
        lb       = float(rxn.get('lower_bound',   0))
        ub       = float(rxn.get('upper_bound', 1000))
        rxn_name = rxn.get('name', rxn_id)

        # Separate into reactants (stoich < 0) and products (stoich > 0).
        # Tuples are (coef, met_id) so _fmt_side unpacks correctly.
        reactants = sorted(
            [(-s, mid) for mid, s in mets.items() if s < 0],
            key=lambda x: x[1]   # sort by metabolite ID for deterministic output
        )
        products = sorted(
            [( s, mid) for mid, s in mets.items() if s > 0],
            key=lambda x: x[1]
        )

        lhs = _fmt_side(reactants)
        rhs = _fmt_side(products)

        if lb < 0 < ub:
            # Bidirectional: split into fwd and rev
            lines.append(
                f'R{idx}_fwd: {lhs} => {rhs}; {rxn_id}: {rxn_name}'
            )
            idx += 1
            lines.append(
                f'R{idx}_rev: {rhs} => {lhs}; {rxn_id} (reverse): {rxn_name}'
            )
            n_split += 1

        elif lb < 0 and ub <= 0:
            # Runs in reverse direction only
            lines.append(
                f'R{idx}: {rhs} => {lhs}; {rxn_id} (rev-only): {rxn_name}'
            )

        else:
            # Forward only (lb >= 0)
            lines.append(
                f'R{idx}: {lhs} => {rhs}; {rxn_id}: {rxn_name}'
            )

        idx += 1

    stats = {
        'n_bigg_reactions': len(model_data.get('reactions', [])),
        'n_pycot_reactions': idx,
        'n_metabolites': len(model_data.get('metabolites', [])),
        'n_reversible_splits': n_split,
    }
    return lines, stats


def convert_and_save(model_id, model_data, out_dir):
    """Convert a BiGG model dict and write the pyCOT .txt file."""
    os.makedirs(out_dir, exist_ok=True)
    lines, stats = convert_model(model_data)

    out_path = os.path.join(out_dir, f'bigg_{model_id}.txt')
    with open(out_path, 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(lines) + '\n')

    return out_path, stats


# ── CLI ───────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(
        description='Download BiGG GEM models and convert to pyCOT .txt format.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument(
        'model_ids', nargs='*',
        help='BiGG model IDs to download (e.g. iJO1366 iML1515 iMM904)',
    )
    ap.add_argument(
        '--list', action='store_true',
        help='Print all available BiGG models sorted by reaction count and exit.',
    )
    ap.add_argument(
        '--all', action='store_true',
        help='Download and convert ALL BiGG models (respects --min-rxns / --max-rxns).',
    )
    ap.add_argument(
        '--local', metavar='JSON_FILE',
        help='Convert a locally-saved BiGG JSON file instead of downloading.',
    )
    ap.add_argument(
        '--out-dir', default=DEFAULT_OUT_DIR,
        help=(
            'Output directory for .txt files '
            f'(default: {DEFAULT_OUT_DIR})'
        ),
    )
    ap.add_argument(
        '--skip-existing', action='store_true', default=True,
        help='Skip models whose .txt file already exists (default: on).',
    )
    ap.add_argument(
        '--no-skip', dest='skip_existing', action='store_false',
        help='Reconvert even if the .txt file already exists.',
    )
    ap.add_argument(
        '--min-rxns', type=int, default=0,
        help='Skip models with fewer reactions than this.',
    )
    ap.add_argument(
        '--max-rxns', type=int, default=999_999,
        help='Skip models with more reactions than this.',
    )
    ap.add_argument(
        '--delay', type=float, default=0.5,
        help='Seconds to wait between consecutive downloads (default: 0.5).',
    )
    args = ap.parse_args()

    out_dir = args.out_dir

    # ── Local JSON conversion ────────────────────────────────────────────────
    if args.local:
        print(f'Converting local file: {args.local}')
        with open(args.local, 'r', encoding='utf-8') as fh:
            model_data = json.load(fh)
        model_id = model_data.get(
            'id',
            os.path.splitext(os.path.basename(args.local))[0]
        )
        out_path, stats = convert_and_save(model_id, model_data, out_dir)
        print(f'  Saved → {out_path}')
        print(f'  {stats}')
        return

    # ── List mode ────────────────────────────────────────────────────────────
    if args.list or (not args.model_ids and not args.all):
        print('Fetching BiGG model list …')
        try:
            models = list_bigg_models()
        except Exception as e:
            sys.exit(f'Cannot reach BiGG API: {e}')

        header = f'{"BiGG ID":<20} {"Organism":<50} {"Rxns":>6} {"Mets":>6}'
        print(header)
        print('-' * len(header))
        for m in models:
            print(
                f"{m['bigg_id']:<20} {m.get('organism','?')[:50]:<50} "
                f"{m.get('reaction_count',0):>6} {m.get('metabolite_count',0):>6}"
            )
        print(f'\nTotal: {len(models)} models')
        return

    # ── Download mode ────────────────────────────────────────────────────────
    if args.all:
        print('Fetching BiGG model list …')
        try:
            all_models = list_bigg_models()
        except Exception as e:
            sys.exit(f'Cannot reach BiGG API: {e}')

        target_ids = [
            m['bigg_id'] for m in all_models
            if args.min_rxns <= m.get('reaction_count', 0) <= args.max_rxns
        ]
        print(
            f'Targeting {len(target_ids)} models '
            f'(reactions {args.min_rxns}–{args.max_rxns})'
        )
    else:
        target_ids = args.model_ids

    if not target_ids:
        ap.print_help()
        return

    ok = skipped = failed = 0
    for model_id in target_ids:
        out_path = os.path.join(out_dir, f'bigg_{model_id}.txt')
        if args.skip_existing and os.path.exists(out_path):
            print(f'  SKIP  {model_id}  (already exists)')
            skipped += 1
            continue

        print(f'  {model_id} … ', end='', flush=True)
        try:
            model_data = download_bigg_json(model_id)
            saved, stats = convert_and_save(model_id, model_data, out_dir)
            print(
                f"OK  {stats['n_pycot_reactions']} rxns "
                f"({stats['n_reversible_splits']} splits, "
                f"{stats['n_metabolites']} metabolites)"
            )
            ok += 1
            time.sleep(args.delay)
        except KeyboardInterrupt:
            print('\nInterrupted.')
            break
        except Exception as e:
            print(f'FAILED — {e}')
            failed += 1

    print(f'\nFinished: {ok} converted, {skipped} skipped, {failed} failed.')
    print(f'Output directory: {out_dir}')


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf-8')
    main()
