"""
run_all_networks.py — Batch COT pipeline over all available networks.

Runs the full analysis (ERC → Hierarchy → Synergy → Complementarity →
Generators → EPM → ESPM) on every reaction network found under the data
directory, subject to reaction-count and ERC-count filters.

Results are written to a single self-extending CSV (RESULTS_CSV).  Dynamic
columns (epms_lenN, espm_oK, leaves_lenN) grow automatically as larger
networks are processed.  An Excel copy is also produced if openpyxl is
installed.

CONFIGURATION
-------------
Edit the block below, then run:
    python projects/COT_Fundamental_Generators_Complex/scripts/run_all_networks.py
"""

# ── Path setup ────────────────────────────────────────────────────────────────
from __future__ import annotations
import os, sys, time, traceback

_here = os.path.dirname(os.path.abspath(__file__))
_proj = os.path.normpath(os.path.join(_here, ".."))
_repo = os.path.normpath(os.path.join(_here, "..", "..", ".."))
for _p in (_proj, os.path.join(_repo, "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  CONFIGURATION                                                           ║
# ╚══════════════════════════════════════════════════════════════════════════╝

MIN_REACTIONS  = 100        # skip networks with fewer reactions (0 = no limit)
MAX_REACTIONS  = 1000      # skip networks with more reactions (0 = no limit)
MAX_ERCS       = 1000      # skip synergy if more ERCs than this (0 = no limit)
COMPUTE_ESPM   = True
ESPM_MAX_ORDER = 3000

# Primary output — one growing CSV shared with run_network.py
RESULTS_CSV = os.path.join(_here, "..", "outputs", "cot_results.csv")

# Optional Excel copy (set to None to skip)
OUTPUT_XLS  = os.path.join(_here, "..", "outputs", "cot_results.xlsx")

# Optional: restrict to one subfolder of data/biomodels/
# e.g. "biomodels_interesting" or "" for all
DATA_SUBFOLDER = "BiGG"

# ╔══════════════════════════════════════════════════════════════════════════╝

# ── Imports ───────────────────────────────────────────────────────────────────
from pyCOT.io.functions      import read_txt
from pyCOT.analysis.organizations.io_pyCOT        import build_rndata
from pyCOT.analysis.organizations.erc              import compute_ercs
from pyCOT.analysis.organizations.hierarchy        import build_hierarchy
from pyCOT.analysis.organizations.synergy          import compute_synergies_basis_first
from pyCOT.analysis.organizations.complementarity  import compute_complementarities
from pyCOT.analysis.organizations.generators       import compute_generators
from pyCOT.analysis.organizations.epm              import compute_epms, compute_espm
from cot_gen.results_io      import make_row, update_results_csv

try:
    import openpyxl
    from openpyxl.styles import Font, PatternFill, Alignment
    from openpyxl.utils  import get_column_letter
    HAS_OPENPYXL = True
except ImportError:
    HAS_OPENPYXL = False

# ── Network discovery ─────────────────────────────────────────────────────────
_DATA_ROOT   = os.path.join(_repo, "data", "biomodels")
_SEARCH_ROOT = os.path.join(_DATA_ROOT, DATA_SUBFOLDER) if DATA_SUBFOLDER else _DATA_ROOT


def _discover_networks(root: str) -> list[tuple[str, str]]:
    found: dict[str, tuple[str, int]] = {}
    for dirpath, _dirs, files in os.walk(root):
        for fname in files:
            if not fname.endswith(".txt"):
                continue
            full = os.path.join(dirpath, fname)
            name = os.path.splitext(fname)[0]
            if name.startswith("bigg_"):
                name = name[5:]
            if name not in found or len(full) < found[name][1]:
                found[name] = (full, len(full))
    return sorted((name, path) for name, (path, _) in found.items())


def _fast_reaction_count(path: str) -> int:
    count = 0
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if "=>" in line:
                count += 1
    return count


# ── Single network pipeline ───────────────────────────────────────────────────

def run_one(name: str, path: str) -> dict:
    """
    Run the full pipeline on one network.
    Writes the result row to RESULTS_CSV immediately (so partial runs
    recover) UNLESS the network exceeds MIN/MAX_REACTIONS or MAX_ERCS, in
    which case nothing is computed and no row is written -- only networks
    actually computed end up in the CSV.
    Returns the row dict (used for console summary and XLS).
    """
    timing: dict[str, float] = {}
    status = "ok"
    rn = ercs_list = hier = syn_fund = comp = epm = espm = None

    try:
        # ── Fast reaction pre-check ──────────────────────────────────────
        # Exceeding MIN/MAX_REACTIONS means this network falls outside the
        # scope of this run entirely: don't compute anything, don't write a
        # row.  Only the console summary (built from this returned dict)
        # ever sees it.
        if MIN_REACTIONS > 0 or MAX_REACTIONS > 0:
            n_rxn = _fast_reaction_count(path)
            if (MIN_REACTIONS > 0 and n_rxn < MIN_REACTIONS) or \
               (MAX_REACTIONS > 0 and n_rxn > MAX_REACTIONS):
                return {"name": name, "reactions": n_rxn, "status": "skipped_rxn",
                        "_written": False}

        # ── Load ──────────────────────────────────────────────────────────
        t0 = time.perf_counter()
        rn_pycot = read_txt(path)
        rn       = build_rndata(rn_pycot, network_id=name)
        timing["load"] = (time.perf_counter() - t0) * 1000

        # ── ERCs ──────────────────────────────────────────────────────────
        t0 = time.perf_counter()
        ercs_list = compute_ercs(rn)
        timing["erc"] = (time.perf_counter() - t0) * 1000

        if not ercs_list:
            status = "no_ercs"
            row = make_row(name, rn, ercs_list, None, None, None, None,
                           timing, status=status)
            update_results_csv(RESULTS_CSV, row)
            row["_written"] = True
            return row

        # Exceeding MAX_ERCS also puts the whole network outside scope:
        # skip the rest of the pipeline and don't write a (degraded) row.
        if MAX_ERCS > 0 and len(ercs_list) > MAX_ERCS:
            return {"name": name, "reactions": rn.n_reactions, "ercs": len(ercs_list),
                    "status": "skipped_ercs", "_written": False}

        # ── Hierarchy ─────────────────────────────────────────────────────
        t0 = time.perf_counter()
        hier = build_hierarchy(ercs_list)
        timing["hier"] = (time.perf_counter() - t0) * 1000

        # ── Synergy ───────────────────────────────────────────────────────
        t0 = time.perf_counter()
        syn_fund = compute_synergies_basis_first(ercs_list, hier)
        timing["syn"] = (time.perf_counter() - t0) * 1000

        # ── Complementarity ───────────────────────────────────────────────
        t0 = time.perf_counter()
        comp = compute_complementarities(ercs_list, hier, syn_result=syn_fund)
        timing["comp"] = (time.perf_counter() - t0) * 1000

        # ── Generators ────────────────────────────────────────────────────
        if syn_fund is not None:
            t0 = time.perf_counter()
            try:
                compute_generators(ercs_list, hier, syn_fund)  # not stored in CSV yet
            except Exception:
                pass
            timing["gen"] = (time.perf_counter() - t0) * 1000
        else:
            timing["gen"] = 0.0

        # ── EPMs ──────────────────────────────────────────────────────────
        t0 = time.perf_counter()
        epm = compute_epms(rn, ercs_list, hier, syn_result=syn_fund,
                           comp_result=comp)
        timing["epm"] = (time.perf_counter() - t0) * 1000

        # ── ESPMs ─────────────────────────────────────────────────────────
        if COMPUTE_ESPM and syn_fund is not None:
            t0 = time.perf_counter()
            try:
                espm = compute_espm(rn, ercs_list, hier, syn_fund, comp, epm,
                                    max_order=ESPM_MAX_ORDER)
            except Exception:
                status = "espm_error"
            timing["espm"] = (time.perf_counter() - t0) * 1000
        else:
            timing["espm"] = 0.0

    except Exception as exc:
        status = f"error:{type(exc).__name__}"
        traceback.print_exc()

    row = make_row(name, rn, ercs_list, syn_fund, comp, epm, espm,
                   timing, status=status)
    update_results_csv(RESULTS_CSV, row)
    row["_written"] = True
    return row


# ── Excel export (regenerated from CSV at the end) ────────────────────────────

def _write_xls_from_csv(csv_path: str, xls_path: str) -> None:
    if not HAS_OPENPYXL:
        return
    import csv as csv_mod

    with open(csv_path, newline="", encoding="utf-8") as fh:
        reader = csv_mod.DictReader(fh)
        headers = reader.fieldnames or []
        rows    = list(reader)

    HEADER_FILL = PatternFill("solid", fgColor="1F4E79")
    HEADER_FONT = Font(bold=True, color="FFFFFF")
    ALT_FILL    = PatternFill("solid", fgColor="D6E4F0")

    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "COT Results"

    for col_i, h in enumerate(headers, 1):
        cell = ws.cell(row=1, column=col_i, value=h)
        cell.font      = HEADER_FONT
        cell.fill      = HEADER_FILL
        cell.alignment = Alignment(horizontal="center", wrap_text=True)
    ws.row_dimensions[1].height = 30

    for row_i, data in enumerate(rows, 2):
        fill = ALT_FILL if row_i % 2 == 0 else None
        for col_i, h in enumerate(headers, 1):
            val = data.get(h, "")
            # coerce numeric strings back to numbers for Excel
            if val not in ("", None):
                try:
                    val = int(val) if "." not in str(val) else float(val)
                except (ValueError, TypeError):
                    pass
            cell = ws.cell(row=row_i, column=col_i, value=val)
            cell.alignment = Alignment(horizontal="center")
            if fill:
                cell.fill = fill

    ws.freeze_panes = "B2"
    for col_cells in ws.columns:
        max_len = max((len(str(c.value or "")) for c in col_cells), default=8)
        ws.column_dimensions[
            get_column_letter(col_cells[0].column)
        ].width = min(max_len + 2, 28)

    os.makedirs(os.path.dirname(os.path.abspath(xls_path)), exist_ok=True)
    wb.save(xls_path)


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    networks = _discover_networks(_SEARCH_ROOT)
    if not networks:
        print(f"No networks found under: {_SEARCH_ROOT}")
        return

    print(f"Found {len(networks)} networks under {_SEARCH_ROOT}")
    print(f"MIN_REACTIONS={MIN_REACTIONS}  MAX_REACTIONS={MAX_REACTIONS}  "
          f"MAX_ERCS={MAX_ERCS}  COMPUTE_ESPM={COMPUTE_ESPM}  "
          f"ESPM_MAX_ORDER={ESPM_MAX_ORDER}")
    print(f"CSV output: {RESULTS_CSV}\n")
    print(f"{'#':>4}  {'Network':<40}  {'Rxn':>5}  {'ERCs':>5}  "
          f"{'EPMs':>7}  {'ESPMs':>7}  {'t_total':>9}  Status")
    print("-" * 90)

    all_rows: list[dict] = []

    for idx, (name, path) in enumerate(networks, 1):
        sys.stdout.write(f"{idx:>4}  {name:<40}")
        sys.stdout.flush()

        row       = run_one(name, path)
        all_rows.append(row)

        rxn    = row.get("reactions",  "?")
        ercs   = row.get("ercs",       "-")
        epms   = row.get("n_epms",     "-")
        espms  = row.get("n_espms",    "-")
        t_tot  = float(row.get("t_total_ms", 0))
        status = row.get("status", "?")
        note   = "" if row.get("_written") else "  (not written to CSV)"

        print(f"  {str(rxn):>5}  {str(ercs):>5}  "
              f"{str(epms):>7}  {str(espms):>7}  "
              f"{t_tot/1000:>8.2f}s  [{status}]{note}")

    # ── Excel export ──────────────────────────────────────────────────────
    if OUTPUT_XLS and HAS_OPENPYXL and os.path.exists(RESULTS_CSV):
        _write_xls_from_csv(RESULTS_CSV, OUTPUT_XLS)
        print(f"\nExcel copy written to: {OUTPUT_XLS}")
    elif OUTPUT_XLS and not HAS_OPENPYXL:
        print("\nSkipping Excel output (openpyxl not installed).")

    # ── Summary ───────────────────────────────────────────────────────────
    ok        = sum(1 for r in all_rows if r.get("status") == "ok")
    skipped   = sum(1 for r in all_rows if "skip" in str(r.get("status", "")))
    errors    = len(all_rows) - ok - skipped
    written   = sum(1 for r in all_rows if r.get("_written"))
    not_written = len(all_rows) - written
    print(f"\nDone: {ok} ok  {skipped} skipped  {errors} errors"
          f"  (of {len(all_rows)} total)")
    print(f"CSV rows written: {written}  "
          f"(outside MIN/MAX_REACTIONS or MAX_ERCS, not written: {not_written})")
    print(f"CSV: {RESULTS_CSV}")


if __name__ == "__main__":
    main()
