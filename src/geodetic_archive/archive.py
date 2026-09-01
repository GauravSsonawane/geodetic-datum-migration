"""Readers for the KA-71 station register and the datum connection stations.

The register is fixed-width, carries three different coordinate encodings and two
different height units, spreads long names across continuation lines, and marks
never-observed quantities with a field of nines. All four of those feed the record
counts a migration is judged on. Column positions here are the 1-based, inclusive
ones data/FORMAT.txt section 2 declares, translated to Python slices.
"""

from __future__ import annotations

import csv
import logging
from dataclasses import dataclass
from pathlib import Path

from geodetic_archive.projection import KAVERI_GRID_1962, TransverseMercator

logger = logging.getLogger(__name__)

#: Record types that carry a station, as opposed to a continuation line.
RECORD_TYPES: tuple[str, ...] = ("G", "D", "P")
CONTINUATION_TYPE = "+"

#: Slices of a register line, from FORMAT.txt section 2.
_TYPE = slice(0, 1)
_IDENTIFIER = slice(1, 8)
_NAME = slice(8, 32)
_COORDINATE = slice(32, 62)
_HEIGHT = slice(62, 71)

#: Width of the coordinate field, so a trailing-blank-trimmed field can be restored.
COORDINATE_FIELD_WIDTH = _COORDINATE.stop - _COORDINATE.start

#: Height field units per record type, in metres per recorded integer unit.
_HEIGHT_UNIT_M = {"G": 0.1, "P": 0.1, "D": 0.01}

_TIE_POINT_COLUMNS = (
    "station_id",
    "legacy_encoding",
    "legacy_coord_field",
    "legacy_ellip_height_mm",
    "wgs84_lat_deg",
    "wgs84_lon_deg",
    "wgs84_ellip_height_m",
)


class ArchiveFormatError(ValueError):
    """Raised when a file does not conform to the KA-71 interchange format."""


def is_void(field: str) -> bool:
    """Whether a field is a void marker: non-empty and all the digit 9.

    A void marker means the quantity was never observed or has been lost, so the
    record cannot be migrated. See FORMAT.txt section 2.4.
    """
    stripped = field.strip()
    return bool(stripped) and set(stripped) == {"9"}


def decode_coordinate(
    field: str,
    encoding: str,
    grid: TransverseMercator = KAVERI_GRID_1962,
) -> tuple[float, float]:
    """Decode one coordinate field to legacy latitude and longitude in degrees.

    Args:
        field: The coordinate field. Trailing blanks may be absent; it is padded
            back to :data:`COORDINATE_FIELD_WIDTH` before slicing.
        encoding: ``"G"`` packed sexagesimal, ``"D"`` decimal degrees, or ``"P"``
            plane coordinates on ``grid``.
        grid: The projection type ``P`` fields are referred to. Defaults to the
            Kaveri Grid 1962, which is defined on the *legacy* ellipsoid.

    Returns:
        ``(latitude_deg, longitude_deg)`` on the legacy datum.

    Raises:
        ArchiveFormatError: If ``encoding`` is not one of the three, or the field
            does not hold what the layout for that encoding declares.
    """
    if encoding not in RECORD_TYPES:
        raise ArchiveFormatError(
            f"unknown coordinate encoding {encoding!r}; expected one of {', '.join(RECORD_TYPES)}"
        )

    padded = field.ljust(COORDINATE_FIELD_WIDTH)
    try:
        if encoding == "G":
            lat_field, lon_field = padded[0:10], padded[11:22]
            latitude = (
                int(lat_field[0:2]) + int(lat_field[2:4]) / 60.0 + int(lat_field[4:9]) / 3.6e6
            )
            longitude = (
                int(lon_field[0:3]) + int(lon_field[3:5]) / 60.0 + int(lon_field[5:10]) / 3.6e6
            )
            return (
                -latitude if lat_field[9] == "S" else latitude,
                -longitude if lon_field[10] == "W" else longitude,
            )
        if encoding == "D":
            return float(padded[0:11]), float(padded[12:24])
        easting_mm, northing_mm = int(padded[0:12]), int(padded[13:25])
        return grid.to_geodetic(easting_mm / 1000.0, northing_mm / 1000.0)
    except ValueError as exc:
        raise ArchiveFormatError(
            f"type {encoding} coordinate field {field.strip()!r} does not match its layout: {exc}"
        ) from exc


@dataclass(frozen=True)
class ArchiveRecord:
    """One assembled station record from the register.

    Args:
        record_id: The identifier from columns 2-8.
        encoding: The record type, which fixes both the coordinate layout and the
            height unit.
        name: The full name, assembled across any continuation lines.
        coordinate_field: The raw coordinate field, as read.
        height_field: The raw height field, as read.
        line_number: Where the record began, for error messages.
    """

    record_id: str
    encoding: str
    name: str
    coordinate_field: str
    height_field: str
    line_number: int

    @property
    def is_complete(self) -> bool:
        """Whether both the position and the height were observed."""
        return not is_void(self.coordinate_field) and not is_void(self.height_field)

    def position(self, grid: TransverseMercator = KAVERI_GRID_1962) -> tuple[float, float]:
        """Legacy latitude and longitude in degrees.

        Raises:
            ArchiveFormatError: If the coordinate field is void or malformed.
        """
        if is_void(self.coordinate_field):
            raise ArchiveFormatError(
                f"record {self.record_id} (line {self.line_number}) has a void coordinate"
            )
        return decode_coordinate(self.coordinate_field, self.encoding, grid)

    def orthometric_height_m(self) -> float:
        """Archived height above local mean sea level, metres.

        The recorded integer is decimetres on type G and P records and centimetres
        on type D records, per FORMAT.txt section 2.2.

        Raises:
            ArchiveFormatError: If the height field is void or not an integer.
        """
        if is_void(self.height_field):
            raise ArchiveFormatError(
                f"record {self.record_id} (line {self.line_number}) has a void height"
            )
        try:
            return int(self.height_field) * _HEIGHT_UNIT_M[self.encoding]
        except ValueError as exc:
            raise ArchiveFormatError(
                f"record {self.record_id} (line {self.line_number}): "
                f"height field {self.height_field!r} is not a signed integer"
            ) from exc


def read_archive(path: str | Path) -> list[ArchiveRecord]:
    """Read the station register, assembling names across continuation lines.

    Args:
        path: The ``archive.dat`` file to read.

    Returns:
        The records in file order. Continuation lines are folded into the record
        they belong to and are not records of their own.

    Raises:
        ArchiveFormatError: If the file cannot be read, a line carries an unknown
            record type, an identifier repeats, or a continuation line refers to a
            record that has not been seen.
    """
    archive_path = Path(path)
    try:
        text = archive_path.read_text(encoding="utf-8")
    except OSError as exc:
        raise ArchiveFormatError(f"cannot read station register {archive_path}: {exc}") from exc

    name_segments: dict[str, list[str]] = {}
    started: dict[str, ArchiveRecord] = {}
    order: list[str] = []

    for number, line in enumerate(text.splitlines(), start=1):
        if not line.strip():
            continue
        record_type = line[_TYPE]
        record_id = line[_IDENTIFIER]

        if record_type == CONTINUATION_TYPE:
            if record_id not in name_segments:
                raise ArchiveFormatError(
                    f"{archive_path} line {number}: continuation for unknown record {record_id!r}"
                )
            name_segments[record_id].append(line[_NAME])
            continue

        if record_type not in RECORD_TYPES:
            raise ArchiveFormatError(
                f"{archive_path} line {number}: unknown record type {record_type!r}; "
                f"expected one of {', '.join((*RECORD_TYPES, CONTINUATION_TYPE))}"
            )
        if record_id in name_segments:
            raise ArchiveFormatError(
                f"{archive_path} line {number}: record {record_id!r} appears more than once"
            )

        name_segments[record_id] = [line[_NAME]]
        started[record_id] = ArchiveRecord(
            record_id=record_id,
            encoding=record_type,
            name="",
            coordinate_field=line[_COORDINATE],
            height_field=line[_HEIGHT],
            line_number=number,
        )
        order.append(record_id)

    records = [
        ArchiveRecord(
            record_id=record_id,
            encoding=started[record_id].encoding,
            name="".join(name_segments[record_id]).rstrip(),
            coordinate_field=started[record_id].coordinate_field,
            height_field=started[record_id].height_field,
            line_number=started[record_id].line_number,
        )
        for record_id in order
    ]
    logger.debug("read %d records from %s", len(records), archive_path)
    return records


@dataclass(frozen=True)
class TiePoint:
    """A datum connection station: the same mark, coordinated on both datums.

    The legacy side is a hand transcription from the 1962-1971 field books and may
    be grossly wrong. The WGS84 side is a modern GNSS result, good to about a
    decimetre. Both heights are already ellipsoidal, so no separation model applies
    to either column.
    """

    station_id: str
    encoding: str
    legacy_coordinate_field: str
    legacy_ellipsoidal_height_m: float
    wgs84_lat_deg: float
    wgs84_lon_deg: float
    wgs84_ellipsoidal_height_m: float

    def legacy_position(self, grid: TransverseMercator = KAVERI_GRID_1962) -> tuple[float, float]:
        """Legacy latitude and longitude in degrees, as transcribed."""
        return decode_coordinate(self.legacy_coordinate_field, self.encoding, grid)


def read_tie_points(path: str | Path) -> list[TiePoint]:
    """Read the datum connection stations.

    Args:
        path: The ``tiepoints.csv`` file to read.

    Returns:
        The stations in file order.

    Raises:
        ArchiveFormatError: If the file cannot be read, a documented column is
            absent, or a numeric column does not hold a number.
    """
    tie_path = Path(path)
    try:
        text = tie_path.read_text(encoding="utf-8")
    except OSError as exc:
        raise ArchiveFormatError(f"cannot read connection stations {tie_path}: {exc}") from exc

    reader = csv.DictReader(text.splitlines())
    present = reader.fieldnames or ()
    missing = [column for column in _TIE_POINT_COLUMNS if column not in present]
    if missing:
        raise ArchiveFormatError(f"{tie_path}: missing column(s) {', '.join(missing)}")

    stations: list[TiePoint] = []
    for number, row in enumerate(reader, start=2):
        try:
            stations.append(
                TiePoint(
                    station_id=row["station_id"],
                    encoding=row["legacy_encoding"],
                    legacy_coordinate_field=row["legacy_coord_field"],
                    legacy_ellipsoidal_height_m=int(row["legacy_ellip_height_mm"]) / 1000.0,
                    wgs84_lat_deg=float(row["wgs84_lat_deg"]),
                    wgs84_lon_deg=float(row["wgs84_lon_deg"]),
                    wgs84_ellipsoidal_height_m=float(row["wgs84_ellip_height_m"]),
                )
            )
        except (TypeError, ValueError) as exc:
            raise ArchiveFormatError(f"{tie_path} line {number}: {exc}") from exc

    logger.debug("read %d connection stations from %s", len(stations), tie_path)
    return stations
