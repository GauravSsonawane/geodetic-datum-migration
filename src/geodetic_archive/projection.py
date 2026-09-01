"""The Kaveri Grid 1962 Transverse Mercator projection.

One of the three coordinate encodings in the register is plane coordinates on a
local Transverse Mercator grid, which have to come back through an inverse
projection *on the legacy ellipsoid* before they are any use. The Redfearn series
is the formulation the 1962 grid was defined with, and is what data/FORMAT.txt
section 1.2 describes.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from geodetic_archive.ellipsoid import KID62, Ellipsoid

#: Footpoint latitude iterations. The correction is metres on the first pass and
#: sub-nanometre by the fourth; the loop exits on tolerance well before this.
_FOOTPOINT_ITERATIONS = 30

#: Metres of northing error at which the footpoint latitude is settled.
_FOOTPOINT_TOLERANCE_M = 1e-9


@dataclass(frozen=True)
class TransverseMercator:
    """A Transverse Mercator grid, defined as a survey document defines one.

    Args:
        ellipsoid: The ellipsoid the grid is computed on. For a legacy grid this
            is the legacy ellipsoid, never WGS84.
        lat_origin_deg: Latitude of the grid origin, degrees.
        lon_origin_deg: Central meridian, degrees.
        scale_factor: Scale factor on the central meridian.
        false_easting_m: Easting added to the central meridian.
        false_northing_m: Northing added to the latitude of origin.
    """

    ellipsoid: Ellipsoid
    lat_origin_deg: float
    lon_origin_deg: float
    scale_factor: float
    false_easting_m: float
    false_northing_m: float

    def _curvature(self, lat: float) -> tuple[float, float, float]:
        """``(nu, rho, eta2)`` at latitude ``lat`` in radians, scaled by ``k0``."""
        a, e2, k0 = self.ellipsoid.a, self.ellipsoid.e2, self.scale_factor
        sin_lat = math.sin(lat)
        w = 1.0 - e2 * sin_lat * sin_lat
        nu = a * k0 / math.sqrt(w)
        rho = a * k0 * (1.0 - e2) / w**1.5
        return nu, rho, nu / rho - 1.0

    def _meridional_arc(self, lat: float) -> float:
        """Meridian distance from the latitude of origin to ``lat``, in metres."""
        a, b = self.ellipsoid.a, self.ellipsoid.b
        n = (a - b) / (a + b)
        lat0 = math.radians(self.lat_origin_deg)
        diff, total = lat - lat0, lat + lat0
        return (
            b
            * self.scale_factor
            * (
                (1 + n + 1.25 * n**2 + 1.25 * n**3) * diff
                - (3 * n + 3 * n**2 + 2.625 * n**3) * math.sin(diff) * math.cos(total)
                + (1.875 * n**2 + 1.875 * n**3) * math.sin(2 * diff) * math.cos(2 * total)
                - (35.0 / 24.0) * n**3 * math.sin(3 * diff) * math.cos(3 * total)
            )
        )

    def from_geodetic(self, lat_deg: float, lon_deg: float) -> tuple[float, float]:
        """Project geodetic coordinates onto the grid.

        Args:
            lat_deg: Latitude on :attr:`ellipsoid`, degrees.
            lon_deg: Longitude on :attr:`ellipsoid`, degrees.

        Returns:
            ``(easting_m, northing_m)``.
        """
        lat, lon = math.radians(lat_deg), math.radians(lon_deg)
        nu, rho, eta2 = self._curvature(lat)
        sin_lat, cos_lat, t = math.sin(lat), math.cos(lat), math.tan(lat)
        d_lon = lon - math.radians(self.lon_origin_deg)

        i = self._meridional_arc(lat) + self.false_northing_m
        ii = nu / 2 * sin_lat * cos_lat
        iii = nu / 24 * sin_lat * cos_lat**3 * (5 - t**2 + 9 * eta2)
        iiia = nu / 720 * sin_lat * cos_lat**5 * (61 - 58 * t**2 + t**4)
        iv = nu * cos_lat
        v = nu / 6 * cos_lat**3 * (nu / rho - t**2)
        vi = nu / 120 * cos_lat**5 * (5 - 18 * t**2 + t**4 + 14 * eta2 - 58 * t**2 * eta2)

        northing = i + ii * d_lon**2 + iii * d_lon**4 + iiia * d_lon**6
        easting = self.false_easting_m + iv * d_lon + v * d_lon**3 + vi * d_lon**5
        return easting, northing

    def to_geodetic(self, easting_m: float, northing_m: float) -> tuple[float, float]:
        """Invert :meth:`from_geodetic`.

        Args:
            easting_m: Grid easting, metres.
            northing_m: Grid northing, metres.

        Returns:
            ``(latitude_deg, longitude_deg)`` on :attr:`ellipsoid`.
        """
        lat = math.radians(self.lat_origin_deg)
        a_k0 = self.ellipsoid.a * self.scale_factor
        for _ in range(_FOOTPOINT_ITERATIONS):
            delta_n = northing_m - self.false_northing_m - self._meridional_arc(lat)
            if abs(delta_n) < _FOOTPOINT_TOLERANCE_M:
                break
            lat += delta_n / a_k0

        nu, rho, eta2 = self._curvature(lat)
        t, sec = math.tan(lat), 1.0 / math.cos(lat)
        vii = t / (2 * rho * nu)
        viii = t / (24 * rho * nu**3) * (5 + 3 * t**2 + eta2 - 9 * t**2 * eta2)
        ix = t / (720 * rho * nu**5) * (61 + 90 * t**2 + 45 * t**4)
        x = sec / nu
        xi = sec / (6 * nu**3) * (nu / rho + 2 * t**2)
        xii = sec / (120 * nu**5) * (5 + 28 * t**2 + 24 * t**4)
        xiia = sec / (5040 * nu**7) * (61 + 662 * t**2 + 1320 * t**4 + 720 * t**6)

        d_e = easting_m - self.false_easting_m
        latitude = lat - vii * d_e**2 + viii * d_e**4 - ix * d_e**6
        longitude = (
            math.radians(self.lon_origin_deg)
            + x * d_e
            - xi * d_e**3
            + xii * d_e**5
            - xiia * d_e**7
        )
        return math.degrees(latitude), math.degrees(longitude)


#: Kaveri Grid 1962, from data/FORMAT.txt section 1.2.
KAVERI_GRID_1962 = TransverseMercator(
    ellipsoid=KID62,
    lat_origin_deg=8.0,
    lon_origin_deg=73.5,
    scale_factor=0.99965,
    false_easting_m=250000.0,
    false_northing_m=150000.0,
)
