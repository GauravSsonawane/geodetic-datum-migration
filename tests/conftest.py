"""Fixtures shared by the host-side suite."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from geodetic_archive.adjustment import Station
from geodetic_archive.archive import read_tie_points
from geodetic_archive.ellipsoid import KID62, WGS84, geodetic_to_ecef

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


@pytest.fixture(scope="session")
def connection_stations(data_dir: Path) -> list[Station]:
    """The shipped connection stations, both datums, in geocentric metres."""
    stations = []
    for tie in read_tie_points(data_dir / "tiepoints.csv"):
        latitude, longitude = tie.legacy_position()
        stations.append(
            Station(
                station_id=tie.station_id,
                legacy_xyz=geodetic_to_ecef(
                    latitude, longitude, tie.legacy_ellipsoidal_height_m, KID62
                ),
                wgs84_xyz=geodetic_to_ecef(
                    tie.wgs84_lat_deg, tie.wgs84_lon_deg, tie.wgs84_ellipsoidal_height_m, WGS84
                ),
                wgs84_lat_deg=tie.wgs84_lat_deg,
                wgs84_lon_deg=tie.wgs84_lon_deg,
            )
        )
    return stations
