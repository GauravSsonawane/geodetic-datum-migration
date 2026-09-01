"""Reference ellipsoids and the geodetic/geocentric conversions between them.

A datum transformation of the Helmert kind lives in geocentric Cartesian space, so
every legacy coordinate has to be carried onto its *own* ellipsoid on the way in.
KID62's semi-major axis is 841 m short of WGS84's; substituting a familiar default
here biases every recovered parameter by more than a hundred metres while leaving
residuals that still look plausible.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

Vector3 = NDArray[np.float64]

#: Iterations of Bowring's method. The series converges to well below a micrometre
#: within four passes at these eccentricities; twelve is free and leaves no doubt.
_GEODETIC_ITERATIONS = 12


@dataclass(frozen=True)
class Ellipsoid:
    """A biaxial reference ellipsoid, defined the way survey documents define one.

    Args:
        name: Human-readable identifier, used in error messages.
        a: Semi-major axis, metres.
        inv_f: Inverse flattening, ``1/f``.
    """

    name: str
    a: float
    inv_f: float

    def __post_init__(self) -> None:
        if self.a <= 0.0:
            raise ValueError(f"{self.name}: semi-major axis must be positive, got {self.a}")
        if self.inv_f <= 0.0:
            raise ValueError(f"{self.name}: inverse flattening must be positive, got {self.inv_f}")

    @property
    def f(self) -> float:
        """Flattening."""
        return 1.0 / self.inv_f

    @property
    def b(self) -> float:
        """Semi-minor axis, metres."""
        return self.a * (1.0 - self.f)

    @property
    def e2(self) -> float:
        """First eccentricity squared."""
        f = self.f
        return 2.0 * f - f * f

    def prime_vertical_radius(self, lat_deg: float) -> float:
        """Radius of curvature in the prime vertical at ``lat_deg``, metres."""
        sin_lat = math.sin(math.radians(lat_deg))
        return self.a / math.sqrt(1.0 - self.e2 * sin_lat * sin_lat)


#: Kaveri Island Datum 1962, from data/FORMAT.txt section 1.1.
KID62 = Ellipsoid(name="Kaveri Island Datum 1962", a=6377295.664, inv_f=300.8017)

#: WGS 84, the target datum.
WGS84 = Ellipsoid(name="WGS 84", a=6378137.0, inv_f=298.257223563)


def geodetic_to_ecef(lat_deg: float, lon_deg: float, height_m: float, ell: Ellipsoid) -> Vector3:
    """Convert geodetic coordinates to earth-centred, earth-fixed Cartesian metres.

    Args:
        lat_deg: Geodetic latitude, degrees, positive north.
        lon_deg: Geodetic longitude, degrees, positive east.
        height_m: Height above the ellipsoid, metres. Not orthometric height.
        ell: The ellipsoid ``lat_deg``, ``lon_deg`` and ``height_m`` are referred to.

    Returns:
        The ``(X, Y, Z)`` position as a length-3 array of metres.
    """
    lat, lon = math.radians(lat_deg), math.radians(lon_deg)
    sin_lat, cos_lat = math.sin(lat), math.cos(lat)
    n = ell.a / math.sqrt(1.0 - ell.e2 * sin_lat * sin_lat)
    return np.array(
        [
            (n + height_m) * cos_lat * math.cos(lon),
            (n + height_m) * cos_lat * math.sin(lon),
            (n * (1.0 - ell.e2) + height_m) * sin_lat,
        ],
        dtype=np.float64,
    )


def ecef_to_geodetic(xyz: Vector3, ell: Ellipsoid) -> tuple[float, float, float]:
    """Invert :func:`geodetic_to_ecef` by Bowring's iteration.

    Args:
        xyz: Geocentric ``(X, Y, Z)`` in metres.
        ell: The ellipsoid to express the result on.

    Returns:
        ``(latitude_deg, longitude_deg, ellipsoidal_height_m)``.

    Raises:
        ValueError: If ``xyz`` is not three finite numbers, or sits on the polar
            axis, where longitude is undefined and the iteration cannot start.
    """
    values = np.asarray(xyz, dtype=np.float64).reshape(-1)
    if values.size != 3 or not np.all(np.isfinite(values)):
        raise ValueError(f"expected three finite coordinates, got {xyz!r}")

    x, y, z = (float(v) for v in values)
    p = math.hypot(x, y)
    if p == 0.0:
        raise ValueError("longitude is undefined on the polar axis (X = Y = 0)")

    lon = math.atan2(y, x)
    lat = math.atan2(z, p * (1.0 - ell.e2))
    height = 0.0
    for _ in range(_GEODETIC_ITERATIONS):
        sin_lat = math.sin(lat)
        n = ell.a / math.sqrt(1.0 - ell.e2 * sin_lat * sin_lat)
        height = p / math.cos(lat) - n
        lat = math.atan2(z, p * (1.0 - ell.e2 * n / (n + height)))
    return math.degrees(lat), math.degrees(lon), height
