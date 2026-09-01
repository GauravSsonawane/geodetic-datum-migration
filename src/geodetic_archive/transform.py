"""The seven-parameter geocentric similarity transformation.

Reported in the position-vector convention (Bursa-Wolf, EPSG 9606), the one
docs/brief.md names. The coordinate-frame convention (EPSG 9607) is the same
transformation with all three rotations negated: mixing them up leaves the
horizontal fit looking healthy while every rotation is reported with the wrong
sign, which is why :class:`HelmertParameters` records the convention explicitly.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray

Vector3 = NDArray[np.float64]
Matrix3 = NDArray[np.float64]

#: Radians per arc-second. Rotations are carried in arc-seconds and scale in parts
#: per million so the normal matrix of the adjustment stays well conditioned.
ARCSEC = math.pi / (180.0 * 3600.0)

#: The keys docs/brief.md fixes for datum_params.json, in order.
PARAMETER_KEYS = (
    "dx_m",
    "dy_m",
    "dz_m",
    "rx_arcsec",
    "ry_arcsec",
    "rz_arcsec",
    "scale_ppm",
)


class ConventionError(ValueError):
    """Raised when a rotation convention is named that this module does not know."""


@dataclass(frozen=True)
class HelmertParameters:
    """Seven parameters taking one geocentric frame onto another.

    Args:
        dx_m: X translation, metres.
        dy_m: Y translation, metres.
        dz_m: Z translation, metres.
        rx_arcsec: Rotation about X, arc-seconds, position-vector sense.
        ry_arcsec: Rotation about Y, arc-seconds, position-vector sense.
        rz_arcsec: Rotation about Z, arc-seconds, position-vector sense.
        scale_ppm: Scale correction, parts per million.
    """

    dx_m: float
    dy_m: float
    dz_m: float
    rx_arcsec: float
    ry_arcsec: float
    rz_arcsec: float
    scale_ppm: float

    @classmethod
    def from_vector(cls, values: NDArray[np.float64]) -> HelmertParameters:
        """Build from the seven-element vector the adjustment solves for."""
        flat = np.asarray(values, dtype=np.float64).reshape(-1)
        if flat.size != 7:
            raise ValueError(f"expected seven parameters, got {flat.size}")
        return cls(*(float(v) for v in flat))

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> HelmertParameters:
        """Build from a datum_params.json payload.

        Raises:
            KeyError: If any of the seven documented keys is absent.
        """
        missing = [key for key in PARAMETER_KEYS if key not in payload]
        if missing:
            raise KeyError(f"datum parameters are missing {', '.join(missing)}")
        return cls(*(float(payload[key]) for key in PARAMETER_KEYS))

    def as_vector(self) -> NDArray[np.float64]:
        """The seven parameters in :data:`PARAMETER_KEYS` order."""
        return np.array(
            [
                self.dx_m,
                self.dy_m,
                self.dz_m,
                self.rx_arcsec,
                self.ry_arcsec,
                self.rz_arcsec,
                self.scale_ppm,
            ],
            dtype=np.float64,
        )

    def to_dict(self, decimals: int | None = None) -> dict[str, float]:
        """As a datum_params.json payload, optionally rounded for serialisation."""
        values = self.as_vector()
        if decimals is not None:
            values = np.round(values, decimals)
        return dict(zip(PARAMETER_KEYS, (float(v) for v in values), strict=True))

    def in_convention(self, convention: str) -> HelmertParameters:
        """Return the same transformation under a named rotation convention.

        Args:
            convention: ``"position_vector"`` (EPSG 9606) or ``"coordinate_frame"``
                (EPSG 9607). Case-insensitive.

        Raises:
            ConventionError: If the name is neither of those.
        """
        name = convention.strip().lower()
        if name == "position_vector":
            return self
        if name == "coordinate_frame":
            return HelmertParameters(
                self.dx_m,
                self.dy_m,
                self.dz_m,
                -self.rx_arcsec,
                -self.ry_arcsec,
                -self.rz_arcsec,
                self.scale_ppm,
            )
        raise ConventionError(
            f"unknown rotation convention {convention!r}; "
            "expected 'position_vector' (EPSG 9606) or 'coordinate_frame' (EPSG 9607)"
        )


def apply_helmert(xyz: Vector3, params: HelmertParameters) -> Vector3:
    """Apply the transformation to one geocentric position.

    Position-vector convention: ``X2 = T + (1 + s) R X1``, with ``R`` linearised
    for the small rotations a regional datum shift involves.

    Args:
        xyz: Source geocentric ``(X, Y, Z)``, metres.
        params: The seven parameters.

    Returns:
        The transformed ``(X, Y, Z)``, metres.
    """
    rx = params.rx_arcsec * ARCSEC
    ry = params.ry_arcsec * ARCSEC
    rz = params.rz_arcsec * ARCSEC
    k = 1.0 + params.scale_ppm * 1e-6
    x, y, z = (float(v) for v in np.asarray(xyz, dtype=np.float64).reshape(-1))
    return np.array(
        [
            params.dx_m + k * (x - rz * y + ry * z),
            params.dy_m + k * (rz * x + y - rx * z),
            params.dz_m + k * (-ry * x + rx * y + z),
        ],
        dtype=np.float64,
    )


def enu_basis(lat_deg: float, lon_deg: float) -> Matrix3:
    """Rows are the east, north and up unit vectors at ``(lat_deg, lon_deg)``.

    Left-multiplying a geocentric difference vector by this puts it in the local
    topocentric frame, which is how residuals are reported.
    """
    lat, lon = math.radians(lat_deg), math.radians(lon_deg)
    sin_lat, cos_lat = math.sin(lat), math.cos(lat)
    sin_lon, cos_lon = math.sin(lon), math.cos(lon)
    return np.array(
        [
            [-sin_lon, cos_lon, 0.0],
            [-sin_lat * cos_lon, -sin_lat * sin_lon, cos_lat],
            [cos_lat * cos_lon, cos_lat * sin_lon, sin_lat],
        ],
        dtype=np.float64,
    )


def to_enu(delta_xyz: Vector3, lat_deg: float, lon_deg: float) -> Vector3:
    """Express a geocentric difference vector as east, north, up metres."""
    basis: Vector3 = enu_basis(lat_deg, lon_deg) @ np.asarray(delta_xyz, dtype=np.float64)
    return basis
