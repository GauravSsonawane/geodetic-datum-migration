"""The Kaveri Grid 1962 Transverse Mercator projection."""

from __future__ import annotations

import pytest

from geodetic_archive.ellipsoid import WGS84
from geodetic_archive.projection import KAVERI_GRID_1962, TransverseMercator

# The archive's extent, from tools/generate_dataset.py.
EXTENT = [(7.62, 72.44), (8.0, 73.5), (8.5, 73.9), (9.88, 74.86)]


def test_grid_definition_matches_the_format_spec() -> None:
    """data/FORMAT.txt section 1.2 is normative for the grid."""
    grid = KAVERI_GRID_1962
    assert grid.lat_origin_deg == pytest.approx(8.0)
    assert grid.lon_origin_deg == pytest.approx(73.5)
    assert grid.scale_factor == pytest.approx(0.99965)
    assert grid.false_easting_m == pytest.approx(250000.0)
    assert grid.false_northing_m == pytest.approx(150000.0)


def test_the_grid_is_computed_on_the_legacy_ellipsoid() -> None:
    """A legacy grid belongs to its own datum; WGS84 here shifts every plane point."""
    assert KAVERI_GRID_1962.ellipsoid.a == pytest.approx(6377295.664)


def test_grid_origin_lands_on_the_false_origin() -> None:
    easting, northing = KAVERI_GRID_1962.from_geodetic(8.0, 73.5)
    assert easting == pytest.approx(250000.0, abs=1e-6)
    assert northing == pytest.approx(150000.0, abs=1e-6)


@pytest.mark.parametrize(("lat", "lon"), EXTENT)
def test_projection_round_trip(lat: float, lon: float) -> None:
    """Sub-millimetre over the archive extent, so the series adds no real error."""
    back_lat, back_lon = KAVERI_GRID_1962.to_geodetic(*KAVERI_GRID_1962.from_geodetic(lat, lon))
    assert back_lat * 111_000 == pytest.approx(lat * 111_000, abs=1e-3)
    assert back_lon * 111_000 == pytest.approx(lon * 111_000, abs=1e-3)


def test_northing_increases_northwards_and_easting_eastwards() -> None:
    south = KAVERI_GRID_1962.from_geodetic(8.0, 73.5)
    north = KAVERI_GRID_1962.from_geodetic(9.0, 73.5)
    east = KAVERI_GRID_1962.from_geodetic(8.0, 74.5)
    assert north[1] > south[1]
    assert east[0] > south[0]


def test_scale_factor_shortens_the_central_meridian() -> None:
    """One degree of latitude on the CM comes out short by the scale factor."""
    _, n0 = KAVERI_GRID_1962.from_geodetic(8.0, 73.5)
    _, n1 = KAVERI_GRID_1962.from_geodetic(9.0, 73.5)
    unscaled = TransverseMercator(
        ellipsoid=KAVERI_GRID_1962.ellipsoid,
        lat_origin_deg=8.0,
        lon_origin_deg=73.5,
        scale_factor=1.0,
        false_easting_m=250000.0,
        false_northing_m=150000.0,
    )
    _, u1 = unscaled.from_geodetic(9.0, 73.5)
    assert (n1 - n0) / (u1 - 150000.0) == pytest.approx(0.99965, rel=1e-9)


def test_a_different_ellipsoid_moves_the_plane_coordinates() -> None:
    """The trap this guards: projecting the legacy grid on the wrong ellipsoid."""
    on_wgs84 = TransverseMercator(
        ellipsoid=WGS84,
        lat_origin_deg=8.0,
        lon_origin_deg=73.5,
        scale_factor=0.99965,
        false_easting_m=250000.0,
        false_northing_m=150000.0,
    )
    legacy = KAVERI_GRID_1962.from_geodetic(8.5, 73.9)
    wrong = on_wgs84.from_geodetic(8.5, 73.9)
    assert abs(legacy[1] - wrong[1]) > 1.0
