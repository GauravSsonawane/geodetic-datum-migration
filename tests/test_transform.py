"""The seven-parameter transformation, its conventions and the topocentric frame."""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from geodetic_archive.ellipsoid import KID62, geodetic_to_ecef
from geodetic_archive.transform import (
    PARAMETER_KEYS,
    ConventionError,
    HelmertParameters,
    apply_helmert,
    enu_basis,
    to_enu,
)

SAMPLE = HelmertParameters(
    dx_m=317.812,
    dy_m=708.176,
    dz_m=234.291,
    rx_arcsec=-1.8734,
    ry_arcsec=2.4171,
    rz_arcsec=-0.9265,
    scale_ppm=2.8137,
)
STATION = geodetic_to_ecef(8.5, 73.5, 120.0, KID62)


def test_parameter_keys_are_the_documented_seven() -> None:
    assert PARAMETER_KEYS == (
        "dx_m",
        "dy_m",
        "dz_m",
        "rx_arcsec",
        "ry_arcsec",
        "rz_arcsec",
        "scale_ppm",
    )


def test_vector_round_trip() -> None:
    assert HelmertParameters.from_vector(SAMPLE.as_vector()) == SAMPLE


def test_from_vector_rejects_the_wrong_length() -> None:
    with pytest.raises(ValueError, match="expected seven parameters"):
        HelmertParameters.from_vector(np.zeros(6))


def test_dict_round_trip() -> None:
    assert HelmertParameters.from_dict(SAMPLE.to_dict()) == SAMPLE


def test_to_dict_rounds_for_serialisation() -> None:
    assert SAMPLE.to_dict(decimals=2)["rx_arcsec"] == pytest.approx(-1.87)


def test_from_dict_names_every_missing_key() -> None:
    payload: dict[str, Any] = SAMPLE.to_dict()
    del payload["rz_arcsec"]
    del payload["scale_ppm"]
    with pytest.raises(KeyError, match="rz_arcsec, scale_ppm"):
        HelmertParameters.from_dict(payload)


def test_translation_only_parameters_shift_by_the_translation() -> None:
    shifted = apply_helmert(STATION, HelmertParameters(1.0, 2.0, 3.0, 0.0, 0.0, 0.0, 0.0))
    assert shifted == pytest.approx(STATION + np.array([1.0, 2.0, 3.0]))


def test_identity_parameters_leave_the_position_alone() -> None:
    identity = HelmertParameters(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    assert apply_helmert(STATION, identity) == pytest.approx(STATION)


def test_scale_stretches_from_the_geocentre() -> None:
    scaled = apply_helmert(STATION, HelmertParameters(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0))
    assert float(np.linalg.norm(scaled)) == pytest.approx(
        float(np.linalg.norm(STATION)) * (1.0 + 1e-6)
    )


def test_the_two_rotation_conventions_differ_by_the_rotation_signs() -> None:
    """EPSG 9607 is EPSG 9606 with all three rotations negated, nothing else."""
    other = SAMPLE.in_convention("coordinate_frame")
    assert other.rx_arcsec == pytest.approx(-SAMPLE.rx_arcsec)
    assert other.ry_arcsec == pytest.approx(-SAMPLE.ry_arcsec)
    assert other.rz_arcsec == pytest.approx(-SAMPLE.rz_arcsec)
    assert (other.dx_m, other.dy_m, other.dz_m) == (SAMPLE.dx_m, SAMPLE.dy_m, SAMPLE.dz_m)
    assert other.scale_ppm == pytest.approx(SAMPLE.scale_ppm)


def test_position_vector_convention_is_the_identity() -> None:
    assert SAMPLE.in_convention("Position_Vector") is SAMPLE


def test_reporting_the_wrong_convention_moves_the_station_metres() -> None:
    """Silent on the horizontal fit, so worth pinning explicitly."""
    correct = apply_helmert(STATION, SAMPLE)
    flipped = apply_helmert(STATION, SAMPLE.in_convention("coordinate_frame"))
    assert float(np.linalg.norm(correct - flipped)) > 50.0


def test_unknown_convention_is_named_in_the_error() -> None:
    with pytest.raises(ConventionError, match="EPSG 9606"):
        SAMPLE.in_convention("bursa_wolf_ish")


def test_enu_basis_is_orthonormal_and_right_handed() -> None:
    basis = enu_basis(8.5, 73.5)
    assert basis @ basis.T == pytest.approx(np.eye(3), abs=1e-12)
    assert float(np.linalg.det(basis)) == pytest.approx(1.0)


def test_up_at_the_equator_points_out_along_the_meridian() -> None:
    east, north, up = enu_basis(0.0, 0.0)
    assert up == pytest.approx(np.array([1.0, 0.0, 0.0]))
    assert east == pytest.approx(np.array([0.0, 1.0, 0.0]))
    assert north == pytest.approx(np.array([0.0, 0.0, 1.0]))


def test_to_enu_preserves_length() -> None:
    delta = np.array([1.5, -2.5, 0.75])
    assert float(np.linalg.norm(to_enu(delta, 8.5, 73.5))) == pytest.approx(
        float(np.linalg.norm(delta))
    )


def test_to_enu_splits_a_purely_vertical_offset_onto_up() -> None:
    up_vector = enu_basis(8.5, 73.5)[2]
    east, north, up = to_enu(up_vector * 3.0, 8.5, 73.5)
    assert up == pytest.approx(3.0)
    assert east == pytest.approx(0.0, abs=1e-12)
    assert north == pytest.approx(0.0, abs=1e-12)
