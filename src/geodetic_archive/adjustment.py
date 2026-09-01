"""Recovering the seven parameters from blunder-contaminated connection stations.

Least squares has no defence against a gross error. One station kilometres out of
place spreads itself across every residual, so fitting once and screening the
residuals fingers the *wrong* stations. The errors here run down a ladder --
16 000 km, 239 km, 73 km, 998 m, 2.7 km, then 27, 8, 3.1, 2.5 and 1.5 m -- and no
single pass separates them:

* A robust screen of the raw difference vectors never sees a least squares
  solution, so it survives the gross errors. But it has to tolerate the honest
  scatter of the difference vectors, which reaches about 2.5 m, so the four
  quietest blunders pass straight through it.
* Those four are only separable *after* a clean fit, where honest residuals top
  out at 0.43 m and the smallest blunder stands at 1.61 m.

So the post-fit cut has to be read off the residual spread in hand. Scaled to the
robust spread of the residuals, or floored near the decimetre quality the GNSS side
is quoted at, it lands inside the (0.43, 1.61) m window comfortably. A habitual
round "2 m" does not: it leaves the 1.5 m station in the fit and returns nine
rejected stations instead of ten.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from geodetic_archive.transform import (
    ARCSEC,
    HelmertParameters,
    apply_helmert,
    to_enu,
)

logger = logging.getLogger(__name__)

Vector3 = NDArray[np.float64]

#: Gauss-Newton iterations. The observation equations are near-linear about the
#: identity for a regional shift, so this converges in one; four settles it.
_FIT_ITERATIONS = 4

#: Scale factor from a median absolute deviation to a standard deviation, for a
#: normally distributed population.
_MAD_TO_SIGMA = 1.4826

#: Multiple of the robust spread at which a raw difference vector is a gross error.
SCREEN_SIGMA = 8.0

#: Floor under the pre-screen cut, at the level honest difference-vector scatter
#: reaches. Below this the screen starts rejecting sound control.
SCREEN_FLOOR_M = 5.0

#: Multiple of the robust residual spread at which a post-fit residual is snooped.
SNOOP_SIGMA = 4.0

#: Floor under the post-fit cut, near the decimetre quality quoted for the GNSS
#: side. This is the number that has to land inside (0.43, 1.61) m; a remembered
#: "2 m" does not.
SNOOP_FLOOR_M = 0.5

#: Redundancy a seven-parameter fit needs before its residuals mean anything.
MIN_ADJUSTED_STATIONS = 10


class AdjustmentError(RuntimeError):
    """Raised when the control cannot support a seven-parameter adjustment."""


@dataclass(frozen=True)
class Station:
    """One connection station, both datums, in geocentric Cartesian metres."""

    station_id: str
    legacy_xyz: Vector3
    wgs84_xyz: Vector3
    wgs84_lat_deg: float
    wgs84_lon_deg: float


@dataclass(frozen=True)
class Adjustment:
    """The recovered transformation and the evidence for it.

    Args:
        parameters: The seven parameters, position-vector convention.
        used_station_ids: Stations carried into the final fit, in input order.
        rejected_station_ids: Stations judged miskeyed, sorted.
        residuals_enu_m: East, north, up residual per used station, metres,
            as GNSS position minus transformed legacy position.
        rms_residual_m: Three-dimensional RMS residual over the used stations.
        screen_cut_m: The pre-screen cut actually applied.
        snoop_cut_m: The post-fit cut at the point the snooping loop stopped.
        snoop_passes: How many stations the post-fit loop removed.
    """

    parameters: HelmertParameters
    used_station_ids: tuple[str, ...]
    rejected_station_ids: tuple[str, ...]
    residuals_enu_m: dict[str, tuple[float, float, float]]
    rms_residual_m: float
    screen_cut_m: float
    snoop_cut_m: float
    snoop_passes: int


def _robust_scale(values: NDArray[np.float64]) -> tuple[float, float]:
    """``(median, sigma)`` of ``values``, with sigma from the median deviation.

    The median absolute deviation is used rather than a standard deviation because
    the blunders being looked for are exactly what would inflate a standard
    deviation past any cut derived from it.
    """
    median = float(np.median(values))
    return median, _MAD_TO_SIGMA * float(np.median(np.abs(values - median)))


def fit_helmert(stations: list[Station]) -> HelmertParameters:
    """Least-squares fit of the seven parameters over ``stations``.

    Linearised about the identity, the observation equations are linear in the
    three translations, the three small rotations and the scale, so a single solve
    would do; the iteration costs nothing and removes the question. Rotations are
    carried in arc-seconds and scale in parts per million inside the design matrix
    so the normal matrix stays well conditioned.

    Args:
        stations: The control to fit. No screening is done here.

    Returns:
        The fitted parameters.

    Raises:
        AdjustmentError: If fewer than seven parameters' worth of control is given.
    """
    if len(stations) < 3:
        raise AdjustmentError(
            f"a seven-parameter fit needs at least 3 stations, got {len(stations)}"
        )

    source = [station.legacy_xyz for station in stations]
    target = [station.wgs84_xyz for station in stations]
    parameters = np.zeros(7, dtype=np.float64)

    for _ in range(_FIT_ITERATIONS):
        design = np.zeros((3 * len(source), 7), dtype=np.float64)
        misclosure = np.zeros(3 * len(source), dtype=np.float64)
        current = HelmertParameters.from_vector(parameters)
        for index, (src, dst) in enumerate(zip(source, target, strict=True)):
            x, y, z = (float(v) for v in src)
            row = 3 * index
            design[row + 0, 0] = design[row + 1, 1] = design[row + 2, 2] = 1.0
            design[row + 0, 4], design[row + 0, 5] = z * ARCSEC, -y * ARCSEC
            design[row + 1, 3], design[row + 1, 5] = -z * ARCSEC, x * ARCSEC
            design[row + 2, 3], design[row + 2, 4] = y * ARCSEC, -x * ARCSEC
            design[row + 0, 6] = x * 1e-6
            design[row + 1, 6] = y * 1e-6
            design[row + 2, 6] = z * 1e-6
            misclosure[row : row + 3] = dst - apply_helmert(src, current)
        parameters = parameters + np.linalg.lstsq(design, misclosure, rcond=None)[0]

    return HelmertParameters.from_vector(parameters)


def screen_gross_errors(
    stations: list[Station],
    sigma: float = SCREEN_SIGMA,
    floor_m: float = SCREEN_FLOOR_M,
) -> tuple[list[Station], list[Station], float]:
    """Cut stations whose raw geocentric difference vector leaves the cluster.

    Over a region this small the rotation and scale terms only move a station a few
    metres relative to its neighbours, so the difference vectors of the honest
    stations sit tightly around their component-wise median. A miskeyed coordinate
    throws its difference vector out by whatever the keying error was worth, which
    is orders of magnitude further.

    The statistic never sees a least squares solution, which is why it survives the
    blunders it is looking for. The cut comes from the spread of the difference
    vectors themselves, on a median absolute deviation the blunders cannot inflate,
    floored at the level honest scatter reaches.

    Args:
        stations: All available control.
        sigma: Multiple of the robust spread at which to cut.
        floor_m: Lower bound on the cut, in metres.

    Returns:
        ``(kept, rejected, cut_m)``.
    """
    differences = np.array(
        [station.wgs84_xyz - station.legacy_xyz for station in stations], dtype=np.float64
    )
    distances = np.linalg.norm(differences - np.median(differences, axis=0), axis=1)
    median, scale = _robust_scale(distances)
    cut = max(median + sigma * scale, floor_m)

    kept = [station for station, value in zip(stations, distances, strict=True) if value <= cut]
    rejected = [station for station, value in zip(stations, distances, strict=True) if value > cut]
    logger.debug(
        "pre-screen cut %.3f m (median %.3f, sigma %.3f) rejected %d of %d stations",
        cut,
        median,
        scale,
        len(rejected),
        len(stations),
    )
    return kept, rejected, cut


def adjust(
    stations: list[Station],
    snoop_sigma: float = SNOOP_SIGMA,
    snoop_floor_m: float = SNOOP_FLOOR_M,
    min_stations: int = MIN_ADJUSTED_STATIONS,
) -> Adjustment:
    """Screen gross errors, fit, then snoop the residuals for what the screen left.

    Args:
        stations: All available control.
        snoop_sigma: Multiple of the robust residual spread at which to snoop.
        snoop_floor_m: Floor under the post-fit cut, metres.
        min_stations: Redundancy to preserve; the loop stops rather than going below.

    Returns:
        The adjustment, its residuals and both cuts.

    Raises:
        AdjustmentError: If the pre-screen leaves fewer than ``min_stations``, which
            means the control is too contaminated to support the fit at all.
    """
    kept, rejected, screen_cut = screen_gross_errors(stations)
    if len(kept) < min_stations:
        raise AdjustmentError(
            f"the gross-error screen left {len(kept)} of {len(stations)} stations, "
            f"below the {min_stations} a seven-parameter fit needs; the control is "
            "too contaminated to adjust"
        )

    snoop_cut = float("inf")
    passes = 0
    while True:
        parameters = fit_helmert(kept)
        residuals = np.array(
            [
                float(np.linalg.norm(station.wgs84_xyz - apply_helmert(station.legacy_xyz, parameters)))
                for station in kept
            ]
        )
        median, scale = _robust_scale(residuals)
        snoop_cut = max(median + snoop_sigma * scale, snoop_floor_m)
        worst = int(np.argmax(residuals))

        if len(kept) <= min_stations or residuals[worst] <= snoop_cut:
            break

        logger.debug(
            "snoop pass %d: cut %.3f m, dropping %s at %.3f m",
            passes + 1,
            snoop_cut,
            kept[worst].station_id,
            residuals[worst],
        )
        rejected.append(kept.pop(worst))
        passes += 1

    residuals_enu: dict[str, tuple[float, float, float]] = {}
    sum_of_squares = 0.0
    for station in kept:
        delta = station.wgs84_xyz - apply_helmert(station.legacy_xyz, parameters)
        east, north, up = to_enu(delta, station.wgs84_lat_deg, station.wgs84_lon_deg)
        residuals_enu[station.station_id] = (float(east), float(north), float(up))
        sum_of_squares += float(delta @ delta)
    rms = float(np.sqrt(sum_of_squares / len(kept)))

    logger.info(
        "adjusted %d stations, rejected %d, 3D RMS %.4f m (screen %.2f m, snoop %.2f m)",
        len(kept),
        len(rejected),
        rms,
        screen_cut,
        snoop_cut,
    )
    return Adjustment(
        parameters=parameters,
        used_station_ids=tuple(station.station_id for station in kept),
        rejected_station_ids=tuple(sorted(station.station_id for station in rejected)),
        residuals_enu_m=residuals_enu,
        rms_residual_m=rms,
        screen_cut_m=screen_cut,
        snoop_cut_m=snoop_cut,
        snoop_passes=passes,
    )
