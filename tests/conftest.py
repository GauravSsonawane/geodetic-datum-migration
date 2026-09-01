"""Fixtures shared by the host-side suite."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data"
GROUND_TRUTH_PATH = REPO_ROOT / "acceptance" / "ground_truth.json"


@pytest.fixture(scope="session")
def data_dir() -> Path:
    """The four shipped input files."""
    assert DATA_DIR.is_dir(), f"the shipped dataset is missing from {DATA_DIR}"
    return DATA_DIR


@pytest.fixture(scope="session")
def ground_truth() -> dict[str, Any]:
    """The generating parameters and true positions, as written by tools/."""
    with GROUND_TRUTH_PATH.open(encoding="utf-8") as handle:
        truth: dict[str, Any] = json.load(handle)
    return truth
