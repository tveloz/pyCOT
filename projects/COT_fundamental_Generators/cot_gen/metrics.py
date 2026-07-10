"""
metrics.py — Instrumentation for cot_gen stages.

Every stage wraps its work in a StageContext (used as a context manager or
via the @stage decorator).  On exit it records:

  • wall_s        — wall-clock seconds
  • peak_mb       — peak resident memory increase (tracemalloc)
  • counters      — user-defined named counters (see Counters below)
  • input fields  — network_id, |M|, |R|, |E|, etc. as passed by the caller

All rows are appended to the module-level LOG list.  Call export_csv() or
export_json() to write them out.  Every run stamps the git SHA, hostname,
and a timestamp.

Usage:
    from cot_gen.metrics import StageContext, Counters

    c = Counters()
    with StageContext("closure", network_id="B237", n_species=26, n_reactions=31,
                      counters=c) as ctx:
        result = my_algo(...)
        c.inc("closures_computed")
    # ctx.row now has the metrics
"""
from __future__ import annotations

import os
import json
import socket
import subprocess
import time
import tracemalloc
from contextlib import contextmanager
from dataclasses import dataclass, field, asdict
from datetime import datetime, timezone
from typing import Any


# ---------------------------------------------------------------------------
# Counters
# ---------------------------------------------------------------------------

class Counters:
    """
    Thread-local (per-stage) accumulator for named integer and float counters.

    Keys are plain strings; use '.' for ad-hoc namespacing, e.g.
    "closure.firings", "erc.dedup_hits".
    """
    def __init__(self):
        self._data: dict[str, int | float] = {}

    def inc(self, key: str, n: int | float = 1) -> None:
        self._data[key] = self._data.get(key, 0) + n

    def set(self, key: str, val: int | float) -> None:
        self._data[key] = val

    def get(self, key: str, default: int | float = 0) -> int | float:
        return self._data.get(key, default)

    def as_dict(self) -> dict[str, int | float]:
        return dict(self._data)

    def __repr__(self) -> str:
        return f"Counters({self._data})"


# ---------------------------------------------------------------------------
# MetricsRow
# ---------------------------------------------------------------------------

@dataclass
class MetricsRow:
    """One row in the global log (one stage run on one network)."""
    timestamp: str
    git_sha: str
    hostname: str
    network_id: str
    stage: str
    wall_s: float
    peak_mb: float
    n_species: int
    n_reactions: int
    n_ercs: int                      # -1 if not known yet at this stage
    counters: dict[str, int | float]
    extra: dict[str, Any]            # stage-specific structured fields


# ---------------------------------------------------------------------------
# Global log
# ---------------------------------------------------------------------------

LOG: list[MetricsRow] = []


def _git_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=os.path.dirname(__file__),
            stderr=subprocess.DEVNULL,
            text=True,
        ).strip()
    except Exception:
        return "unknown"


def _hostname() -> str:
    try:
        return socket.gethostname()
    except Exception:
        return "unknown"


# cached once per process
_GIT_SHA = _git_sha()
_HOSTNAME = _hostname()


# ---------------------------------------------------------------------------
# StageContext
# ---------------------------------------------------------------------------

class StageContext:
    """
    Context manager that records wall time and peak memory for a stage.

    Parameters
    ----------
    stage : str
        Name of the stage (e.g. "closure", "erc", "hierarchy").
    network_id : str
        Identifier of the network being processed (e.g. "BIOMD0000000237").
    n_species : int
        Number of species in the network (before E0 quotient).
    n_reactions : int
        Number of reactions.
    n_ercs : int
        Number of ERCs (-1 if not yet known).
    counters : Counters
        User-defined counter object (shared with the caller so they can
        increment inside the `with` block).
    extra : dict
        Any additional structured fields to include in the row.
    append_log : bool
        If True (default), append the finished row to the global LOG.

    Usage
    -----
    c = Counters()
    with StageContext("erc", "B237", 26, 31, counters=c) as ctx:
        ...do work...
        c.inc("closures_computed", 31)
    row = ctx.row  # MetricsRow populated after __exit__
    """

    def __init__(
        self,
        stage: str,
        network_id: str = "",
        n_species: int = -1,
        n_reactions: int = -1,
        n_ercs: int = -1,
        counters: Counters | None = None,
        extra: dict | None = None,
        append_log: bool = True,
    ):
        self.stage = stage
        self.network_id = network_id
        self.n_species = n_species
        self.n_reactions = n_reactions
        self.n_ercs = n_ercs
        self.counters = counters if counters is not None else Counters()
        self.extra = extra or {}
        self.append_log = append_log
        self.row: MetricsRow | None = None

    def __enter__(self) -> 'StageContext':
        tracemalloc.start()
        self._t0 = time.perf_counter()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> bool:
        wall_s = time.perf_counter() - self._t0
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        peak_mb = peak / (1024 * 1024)

        self.row = MetricsRow(
            timestamp=datetime.now(timezone.utc).isoformat(),
            git_sha=_GIT_SHA,
            hostname=_HOSTNAME,
            network_id=self.network_id,
            stage=self.stage,
            wall_s=wall_s,
            peak_mb=peak_mb,
            n_species=self.n_species,
            n_reactions=self.n_reactions,
            n_ercs=self.n_ercs,
            counters=self.counters.as_dict(),
            extra=dict(self.extra),
        )
        if self.append_log:
            LOG.append(self.row)
        return False  # do not suppress exceptions

    # ---- output helpers ----------------------------------------------------

    def summary(self) -> str:
        if self.row is None:
            return f"StageContext({self.stage!r}) — not yet completed"
        r = self.row
        ctr_str = ", ".join(f"{k}={v}" for k, v in sorted(r.counters.items()))
        return (
            f"[{r.stage}] {r.network_id}  "
            f"|M|={r.n_species} |R|={r.n_reactions} |E|={r.n_ercs}  "
            f"wall={r.wall_s:.4f}s  peak={r.peak_mb:.2f}MB  "
            f"counters: {ctr_str}"
        )


# ---------------------------------------------------------------------------
# Export helpers
# ---------------------------------------------------------------------------

def export_csv(path: str, rows: list[MetricsRow] | None = None) -> None:
    """Write LOG (or a subset) to a tidy CSV."""
    import csv
    rows = rows if rows is not None else LOG
    if not rows:
        return

    # Flatten each row: counters and extra become individual columns
    all_counter_keys: set[str] = set()
    all_extra_keys: set[str] = set()
    for r in rows:
        all_counter_keys.update(r.counters.keys())
        all_extra_keys.update(r.extra.keys())

    base_fields = [
        "timestamp", "git_sha", "hostname", "network_id", "stage",
        "wall_s", "peak_mb", "n_species", "n_reactions", "n_ercs",
    ]
    counter_cols = sorted(all_counter_keys)
    extra_cols = sorted(all_extra_keys)
    fieldnames = base_fields + [f"ctr.{k}" for k in counter_cols] + extra_cols

    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            row_dict = {
                "timestamp": r.timestamp, "git_sha": r.git_sha,
                "hostname": r.hostname, "network_id": r.network_id,
                "stage": r.stage, "wall_s": r.wall_s, "peak_mb": r.peak_mb,
                "n_species": r.n_species, "n_reactions": r.n_reactions,
                "n_ercs": r.n_ercs,
            }
            for k in counter_cols:
                row_dict[f"ctr.{k}"] = r.counters.get(k, "")
            for k in extra_cols:
                row_dict[k] = r.extra.get(k, "")
            writer.writerow(row_dict)


def export_json(path: str, rows: list[MetricsRow] | None = None) -> None:
    """Write LOG (or a subset) to JSON."""
    rows = rows if rows is not None else LOG
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        json.dump([asdict(r) for r in rows], f, indent=2)


def print_summary(rows: list[MetricsRow] | None = None) -> None:
    """Print a human-readable summary of all logged stages."""
    rows = rows if rows is not None else LOG
    for r in rows:
        ctr_str = "  ".join(f"{k}={v}" for k, v in sorted(r.counters.items()))
        print(
            f"[{r.stage:<16}] {r.network_id:<20}  "
            f"|M|={r.n_species:>4} |R|={r.n_reactions:>4} |E|={r.n_ercs:>4}  "
            f"wall={r.wall_s:>8.4f}s  peak={r.peak_mb:>7.2f}MB  {ctr_str}"
        )
