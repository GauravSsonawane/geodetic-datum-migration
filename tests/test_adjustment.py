"""The two-stage blunder screen and the seven-parameter adjustment.

The dataset's whole point lives here, so these tests assert the two properties the
migration turns on: that one screening pass is not enough, and that the post-fit
cut has to be read off the residuals rather than remembered.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from geodetic_archive.adjustment import (
    Adjustment,
    AdjustmentError,
    Station,
    adjust,
    fit_helmert,
    screen_gross_errors,
)
from geodetic_archive.ellipsoid import KID62, WGS84, ecef_to_geodetic, geodetic_to_ecef
from geodetic_archive.transform import HelmertParameters, apply_helmert

TRUE_PARAMETERS = HelmertParameters(
    dx_m=317.812,
    dy_m=708.176,
    dz_m=234.291,
    rx_arcsec=-1.8734,
    ry_arcsec=2.4171,
    rz_arcsec=-0.9265,
    scale_ppm=2.8137,
)


def _synthetic(count: int = 20, noise_m: float = 0.0, seed: int = 7) -> list[Station]:
    """Clean control generated from TRUE_PARAMETERS, optionally with GNSS noise."""
    rng = np.random.default_rng(seed)
    stations = []
    for index in range(count):
        latitude = 7.7 + 2.0 * rng.random()
        longitude = 72.5 + 2.2 * rng.random()
        height = 2000.0 * rng.random()
        legacy = geodetic_to_ecef(latitude, longitude, height, KID62)
        wgs84 = apply_helmert(legacy, TRUE_PARAMETERS) + rng.normal(0.0, noise_m, 3)
        stations.append(
            Station(
                station_id=f"SY{index:03d}",
                legacy_xyz=legacy,
                wgs84_xyz=wgs84,
                wgs84_lat_deg=latitude,
                wgs84_lon_deg=longitude,
            )
        )
    return stations


def _displaced(station: Station, metres: float) -> Station:
    """The same station with its legacy coordinate miskeyed by ``metres``."""
    offset = np.array([metres, 0.0, 0.0])
    return Station(
        station_id=station.station_id,
        legacy_xyz=station.legacy_xyz + offset,
        wgs84_xyz=station.wgs84_xyz,
        wgs84_lat_deg=station.wgs84_lat_deg,
        wgs84_lon_deg=station.wgs84_lon_deg,
    )


# --- the fit -----------------------------------------------------------------
def test_clean_control_recovers_the_generating_parameters_exactly() -> None:
    fitted = fit_helmert(_synthetic())
    assert fitted.dx_m == pytest.approx(TRUE_PARAMETERS.dx_m, abs=1e-4)
    assert fitted.rx_arcsec == pytest.approx(TRUE_PARAMETERS.rx_arcsec, abs=1e-6)
    assert fitted.scale_ppm == pytest.approx(TRUE_PARAMETERS.scale_ppm, abs=1e-6)


def test_noisy_control_recovers_the_parameters_within_the_formal_sigma() -> None:
    fitted = fit_helmert(_synthetic(noise_m=0.15))
    assert fitted.dx_m == pytest.approx(TRUE_PARAMETERS.dx_m, abs=10.0)
    assert fitted.rx_arcsec == pytest.approx(TRUE_PARAMETERS.rx_arcsec, abs=0.35)
    assert fitted.scale_ppm == pytest.approx(TRUE_PARAMETERS.scale_ppm, abs=1.0)


def test_a_fit_with_too_little_control_is_refused() -> None:
    with pytest.raises(AdjustmentError, match="needs at least 3 stations"):
        fit_helmert(_synthetic(count=2))


def test_one_gross_error_wrecks_an_unscreened_fit() -> None:
    """Why the screen exists: least squares has no defence against this."""
    stations = _synthetic()
    stations[3] = _displaced(stations[3], 16_000_000.0)
    fitted = fit_helmert(stations)
    assert abs(fitted.dx_m - TRUE_PARAMETERS.dx_m) > 1000.0


# --- the pre-screen ----------------------------------------------------------
def test_the_screen_cuts_gross_errors_and_keeps_the_rest() -> None:
    stations = _synthetic(noise_m=0.15)
    stations[2] = _displaced(stations[2], 73_000.0)
    stations[9] = _displaced(stations[9], 998.0)
    kept, rejected, cut = screen_gross_errors(stations)
    assert {station.station_id for station in rejected} == {"SY002", "SY009"}
    assert len(kept) == len(stations) - 2
    assert cut >= 5.0


def test_the_screen_never_cuts_below_its_floor() -> None:
    """Honest difference-vector scatter reaches metres; the floor protects it."""
    _, rejected, cut = screen_gross_errors(_synthetic(noise_m=0.15))
    assert rejected == []
    assert cut == pytest.approx(5.0)


def test_the_screen_cannot_see_the_quiet_blunders() -> None:
    """The reason a second, post-fit pass is load-bearing rather than belt-and-braces."""
    stations = _synthetic(noise_m=0.15)
    stations[5] = _displaced(stations[5], 1.5)
    stations[6] = _displaced(stations[6], 2.5)
    _, rejected, _ = screen_gross_errors(stations)
    assert rejected == [], "a metre-scale error is inside the honest scatter"


# --- the two stages together -------------------------------------------------
def test_the_post_fit_pass_catches_what_the_screen_left() -> None:
    stations = _synthetic(noise_m=0.05, count=24)
    stations[1] = _displaced(stations[1], 16_000_000.0)
    stations[5] = _displaced(stations[5], 1.6)
    stations[6] = _displaced(stations[6], 3.1)
    result = adjust(stations)
    assert set(result.rejected_station_ids) == {"SY001", "SY005", "SY006"}
    assert result.snoop_passes == 2


def test_control_too_contaminated_to_adjust_is_refused() -> None:
    stations = _synthetic(count=12)
    for index in range(5):
        stations[index] = _displaced(stations[index], 50_000.0 * (index + 1))
    with pytest.raises(AdjustmentError, match="too contaminated to adjust"):
        adjust(stations)


def test_the_snoop_loop_preserves_the_minimum_redundancy() -> None:
    """It stops at the redundancy floor rather than fitting seven parameters to noise."""
    result = adjust(_synthetic(count=10, noise_m=2.0), snoop_floor_m=0.01)
    assert len(result.used_station_ids) == 10
    assert result.snoop_passes == 0


# --- the shipped control -----------------------------------------------------
def test_the_shipped_control_rejects_exactly_the_miskeyed_stations(
    connection_stations: list[Station], ground_truth: dict[str, Any]
) -> None:
    result = adjust(connection_stations)
    assert list(result.rejected_station_ids) == ground_truth["blunder_tie_points"]
    assert len(result.used_station_ids) == 24


def test_the_shipped_adjustment_recovers_the_true_parameters(
    connection_stations: list[Station], ground_truth: dict[str, Any]
) -> None:
    """Bands are the formal sigma of the fit, roughly four sigma; see docs/design.md."""
    recovered = adjust(connection_stations).parameters
    truth = ground_truth["true_params"]
    for key, true_key, tolerance in [
        ("dx_m", "dx", 10.0),
        ("dy_m", "dy", 10.0),
        ("dz_m", "dz", 10.0),
        ("rx_arcsec", "rx", 0.35),
        ("ry_arcsec", "ry", 0.35),
        ("rz_arcsec", "rz", 0.35),
        ("scale_ppm", "s", 1.0),
    ]:
        assert getattr(recovered, key) == pytest.approx(truth[true_key], abs=tolerance)


def test_the_shipped_adjustment_fits_to_gnss_quality(
    connection_stations: list[Station],
) -> None:
    """Decimetre control cannot fit much better or much worse than this."""
    result = adjust(connection_stations)
    assert 0.02 <= result.rms_residual_m <= 1.00
    assert result.rms_residual_m == pytest.approx(0.19, abs=0.02)


def test_a_habitual_two_metre_cut_returns_nine_stations_instead_of_ten(
    connection_stations: list[Station],
) -> None:
    """The crux. A round number that sounds right leaves the 1.5 m station in."""
    careless = adjust(connection_stations, snoop_floor_m=2.0, snoop_sigma=2.0)
    assert len(careless.rejected_station_ids) == 9
    assert "TP004" not in careless.rejected_station_ids


@pytest.mark.parametrize("floor_m", [0.5, 1.0])
@pytest.mark.parametrize("sigma", [3.0, 4.0, 5.0, 6.0])
def test_any_cut_inside_the_window_returns_the_same_set(
    connection_stations: list[Station],
    ground_truth: dict[str, Any],
    floor_m: float,
    sigma: float,
) -> None:
    """Honest residuals top out at 0.43 m and the smallest blunder sits at 1.61 m."""
    result = adjust(connection_stations, snoop_sigma=sigma, snoop_floor_m=floor_m)
    assert list(result.rejected_station_ids) == ground_truth["blunder_tie_points"]


def test_residuals_are_reported_for_every_used_station(
    connection_stations: list[Station],
) -> None:
    result = adjust(connection_stations)
    assert set(result.residuals_enu_m) == set(result.used_station_ids)
    assert all(len(triple) == 3 for triple in result.residuals_enu_m.values())


def test_the_reported_rms_matches_the_reported_residuals(
    connection_stations: list[Station],
) -> None:
    result = adjust(connection_stations)
    triples = np.array(list(result.residuals_enu_m.values()))
    assert float(np.sqrt(np.mean(np.sum(triples**2, axis=1)))) == pytest.approx(
        result.rms_residual_m, rel=1e-9
    )


def test_the_adjustment_records_both_cuts(connection_stations: list[Station]) -> None:
    """Both numbers are evidence for the result, so they are on the return value."""
    result: Adjustment = adjust(connection_stations)
    assert result.screen_cut_m >= 5.0
    assert 0.43 < result.snoop_cut_m < 1.61
    assert result.snoop_passes == 4


def test_a_wgs84_legacy_side_biases_every_parameter(
    connection_stations: list[Station],
) -> None:
    """Carrying the legacy coordinates through the wrong ellipsoid, the loudest trap."""
    wrong = [
        Station(
            station_id=station.station_id,
            legacy_xyz=geodetic_to_ecef(*ecef_to_geodetic(station.legacy_xyz, KID62), WGS84),
            wgs84_xyz=station.wgs84_xyz,
            wgs84_lat_deg=station.wgs84_lat_deg,
            wgs84_lon_deg=station.wgs84_lon_deg,
        )
        for station in connection_stations
    ]
    biased = adjust(wrong).parameters
    assert abs(biased.dy_m - TRUE_PARAMETERS.dy_m) > 100.0
    assert abs(biased.scale_ppm - TRUE_PARAMETERS.scale_ppm) > 100.0
