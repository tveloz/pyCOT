"""
conftest.py — Shared pytest fixtures for cot_gen tests.

Fixtures
--------
biomd237_path   : path to BIOMD0000000237.txt
biomd237_rn     : pyCOT ReactionNetwork for Biomodel 237
biomd237_rndata : cot_gen RNData for Biomodel 237
"""
from __future__ import annotations

import os
import sys

import pytest

# Layout: pyCOT/projects/COT_fundamental_Generators/tests/
# _proj = .../COT_fundamental_Generators/  (cot_gen + oracles importable from here)
# _repo = .../pyCOT/                       (pyCOT importable from _repo/src, data from _repo/data)
_here = os.path.normpath(os.path.dirname(os.path.abspath(__file__)))
_proj = os.path.normpath(os.path.join(_here, ".."))
_repo = os.path.normpath(os.path.join(_here, "..", "..", ".."))
for _p in [_proj, os.path.join(_repo, "src")]:
    if _p not in sys.path:
        sys.path.insert(0, _p)

from pyCOT.io.functions import read_txt
from pyCOT.analysis.organizations.io_pyCOT import build_rndata

_DATA_DIR = os.path.join(_repo, "data", "biomodels", "BioMD_other")
_BIOMD237 = os.path.join(_DATA_DIR, "BIOMD0000000237.txt")


@pytest.fixture(scope="session")
def biomd237_path():
    if not os.path.exists(_BIOMD237):
        pytest.skip(f"Biomodel 237 not found at {_BIOMD237}")
    return _BIOMD237


@pytest.fixture(scope="session")
def biomd237_rn(biomd237_path):
    return read_txt(biomd237_path)


@pytest.fixture(scope="session")
def biomd237_rndata(biomd237_rn):
    return build_rndata(biomd237_rn, network_id="BIOMD0000000237")
