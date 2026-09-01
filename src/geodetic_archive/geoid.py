"""The geoid-ellipsoid separation model, and the vertical datum it resolves.

Heights in the register are orthometric over local mean sea level, while both the
transformation and the requested output are ellipsoidal. Skipping this step leaves
the horizontal work checking out perfectly and every height wrong by about 38 m,
which is the quietest of the failure modes in the archive.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass
from pathlib import Path

logger = logging.getLogger(__name__)

#: Header keys data/FORMAT.txt section 3 requires of a separation model.
_REQUIRED_HEADER_KEYS = ("lat_min", "lon_min", "d_lat", "d_lon", "n_lat", "n_lon")


class GeoidGridError(RuntimeError):
    """Raised when a separation model cannot be read or does not describe a grid."""


@dataclass(frozen=True)
class GeoidGrid:
    """A regular latitude/longitude grid of separation values.

    Args:
        lat_min: Latitude of the first row, degrees.
        lon_min: Longitude of the first column, degrees.
        d_lat: Row spacing, degrees.
        d_lon: Column spacing, degrees.
        values: ``n_lat`` rows of ``n_lon`` separations ``N = h - H``, metres.
        source: Where the grid was read from, carried for error messages.
    """

    lat_min: float
    lon_min: float
    d_lat: float
    d_lon: float
    values: tuple[tuple[float, ...], ...]
    source: str = "<memory>"

    @property
    def n_lat(self) -> int:
        """Number of rows."""
        return len(self.values)

    @property
    def n_lon(self) -> int:
        """Number of columns."""
        return len(self.values[0])

    @classmethod
    def from_file(cls, path: str | Path) -> GeoidGrid:
        """Read a separation model in the KA-71 grid format.

        Args:
            path: The ``.grd`` file to read.

        Returns:
            The parsed grid.

        Raises:
            GeoidGridError: If the file is missing or unreadable, a header key is
                absent, a value is not a number, or the body does not match the
                ``n_lat`` by ``n_lon`` shape the header declares.
        """
        grid_path = Path(path)
        try:
            text = grid_path.read_text(encoding="utf-8")
        except OSError as exc:
            raise GeoidGridError(f"cannot read separation model {grid_path}: {exc}") from exc

        header: dict[str, float] = {}
        rows: list[tuple[float, ...]] = []
        for number, line in enumerate(text.splitlines(), start=1):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.split()
            try:
                if len(fields) == 2 and fields[0][0].isalpha():
                    header[fields[0]] = float(fields[1])
                else:
                    rows.append(tuple(float(value) for value in fields))
            except ValueError as exc:
                raise GeoidGridError(
                    f"{grid_path} line {number}: expected numbers, got {line.strip()!r}"
                ) from exc

        missing = [key for key in _REQUIRED_HEADER_KEYS if key not in header]
        if missing:
            raise GeoidGridError(f"{grid_path}: header is missing {', '.join(missing)}")

        n_lat, n_lon = int(header["n_lat"]), int(header["n_lon"])
        if n_lat < 2 or n_lon < 2:
            raise GeoidGridError(
                f"{grid_path}: bilinear interpolation needs at least 2x2, "
                f"header declares {n_lat}x{n_lon}"
            )
        if len(rows) != n_lat or any(len(row) != n_lon for row in rows):
            shape = f"{len(rows)} rows of {sorted({len(row) for row in rows})}"
            raise GeoidGridError(
                f"{grid_path}: header declares {n_lat}x{n_lon} but the body holds {shape}"
            )

        logger.debug(
            "read separation model %s: %dx%d from (%.3f, %.3f) at %.3f deg",
            grid_path,
            n_lat,
            n_lon,
            header["lat_min"],
            header["lon_min"],
            header["d_lat"],
        )
        return cls(
            lat_min=header["lat_min"],
            lon_min=header["lon_min"],
            d_lat=header["d_lat"],
            d_lon=header["d_lon"],
            values=tuple(rows),
            source=str(grid_path),
        )

    def separation(self, lat_deg: float, lon_deg: float) -> float:
        """Bilinearly interpolated separation ``N = h - H`` at a position, metres.

        The cell index is clamped to the grid, so a position outside it extrapolates
        linearly from the nearest cell rather than raising. The archive's extent sits
        well inside the model, so this never arises on the shipped data; it keeps the
        caller from having to guard every one of 1200 records.
        """
        fi = (lat_deg - self.lat_min) / self.d_lat
        fj = (lon_deg - self.lon_min) / self.d_lon
        i = min(max(math.floor(fi), 0), self.n_lat - 2)
        j = min(max(math.floor(fj), 0), self.n_lon - 2)
        di, dj = fi - i, fj - j
        return (
            (1 - di) * (1 - dj) * self.values[i][j]
            + (1 - di) * dj * self.values[i][j + 1]
            + di * (1 - dj) * self.values[i + 1][j]
            + di * dj * self.values[i + 1][j + 1]
        )

    def to_ellipsoidal(self, orthometric_m: float, lat_deg: float, lon_deg: float) -> float:
        """Convert an orthometric height to an ellipsoidal one: ``h = H + N``."""
        return orthometric_m + self.separation(lat_deg, lon_deg)
