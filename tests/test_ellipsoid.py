"""The ellipsoid definitions and the geodetic/geocentric pair."""

from __future__ import annotations

import math

import numpy as np
import pytest

from geodetic_archive.ellipsoid import (
    KID62,
    WGS84,
    Ellipsoid,
    ecef_to_geodetic,
    geodetic_to_ecef,
)

# Somewhere on Kaveri Island, in the middle of the archive's extent.
SAMPLE = (8.5, 73.5, 120.0)


def test_shipped_ellipsoids_match_the_format_spec() -> None:
    """data/FORMAT.txt section 1.1 is normative for the legacy ellipsoid."""
    assert KID62.a == pytest.approx(6377295.664)
    assert KID62.inv_f == pytest.approx(300.8017)
    assert WGS84.a == pytest.approx(6378137.0)
    assert WGS84.inv_f == pytest.approx(298.257223563)


def test_legacy_ellipsoid_differs_enough_to_matter() -> None:
    """The 841 m axis difference is the whole reason the legacy side is separate."""
    assert WGS84.a - KID62.a == pytest.approx(841.336, abs=1e-3)


def test_derived_ellipsoid_quantities() -> None:
    assert KID62.f == pytest.approx(1.0 / 300.8017)
    assert KID62.b == pytest.approx(KID62.a * (1.0 - KID62.f))
    assert KID62.e2 == pytest.approx(2.0 * KID62.f - KID62.f**2)
    assert 0.0 < KID62.e2 < 0.01


def test_prime_vertical_radius_grows_towards_the_pole() -> None:
    assert KID62.prime_vertical_radius(0.0) == pytest.approx(KID62.a)
    assert KID62.prime_vertical_radius(60.0) > KID62.prime_vertical_radius(0.0)


@pytest.mark.parametrize(
    ("a", "inv_f"),
    [(0.0, 300.0), (-1.0, 300.0), (6378137.0, 0.0), (6378137.0, -298.0)],
)
def test_degenerate_ellipsoids_are_rejected(a: float, inv_f: float) -> None:
    with pytest.raises(ValueError, match="must be positive"):
        Ellipsoid(name="broken", a=a, inv_f=inv_f)


@pytest.mark.parametrize("ell", [KID62, WGS84])
@pytest.mark.parametrize(
    "position",
    [SAMPLE, (7.62, 72.44, -5.0), (9.88, 74.86, 2180.0), (-8.5, -73.5, 0.0)],
)
def test_geodetic_ecef_round_trip(ell: Ellipsoid, position: tuple[float, float, float]) -> None:
    """Sub-micrometre agreement, so the pair adds no error of its own."""
    lat, lon, height = ecef_to_geodetic(geodetic_to_ecef(*position, ell), ell)
    assert lat == pytest.approx(position[0], abs=1e-11)
    assert lon == pytest.approx(position[1], abs=1e-11)
    assert height == pytest.approx(position[2], abs=1e-6)


def test_equatorial_point_sits_on_the_semi_major_axis() -> None:
    xyz = geodetic_to_ecef(0.0, 0.0, 0.0, WGS84)
    assert xyz[0] == pytest.approx(WGS84.a)
    assert xyz[1] == pytest.approx(0.0, abs=1e-9)
    assert xyz[2] == pytest.approx(0.0, abs=1e-9)


def test_height_moves_the_point_along_the_normal() -> None:
    base = geodetic_to_ecef(*SAMPLE[:2], 0.0, KID62)
    raised = geodetic_to_ecef(*SAMPLE[:2], 1000.0, KID62)
    assert float(np.linalg.norm(raised - base)) == pytest.approx(1000.0, abs=1e-6)


def test_the_two_ellipsoids_disagree_by_hundreds_of_metres() -> None:
    """Running the legacy side through WGS84 is the failure this guards against."""
    correct = geodetic_to_ecef(*SAMPLE, KID62)
    wrong = geodetic_to_ecef(*SAMPLE, WGS84)
    assert float(np.linalg.norm(correct - wrong)) > 100.0


@pytest.mark.parametrize(
    "bad",
    [
        [1.0, 2.0],
        [1.0, 2.0, 3.0, 4.0],
        [float("nan"), 0.0, 0.0],
        [float("inf"), 0.0, 0.0],
    ],
)
def test_ecef_to_geodetic_rejects_malformed_input(bad: list[float]) -> None:
    with pytest.raises(ValueError, match="three finite coordinates"):
        ecef_to_geodetic(np.asarray(bad, dtype=float), WGS84)


def test_ecef_to_geodetic_rejects_the_polar_axis() -> None:
    with pytest.raises(ValueError, match="polar axis"):
        ecef_to_geodetic(np.array([0.0, 0.0, KID62.b]), KID62)


def test_longitude_sign_convention_is_east_positive() -> None:
    east = geodetic_to_ecef(0.0, 90.0, 0.0, WGS84)
    assert east[1] == pytest.approx(WGS84.a)
    assert math.copysign(1.0, east[1]) > 0.0
