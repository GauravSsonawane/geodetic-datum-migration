"""End-to-end migration of the register from the legacy datum onto WGS84.

Order matters. The transformation has to be recovered from the connection stations
before a single record can be moved, and each record has to be taken from its
archived orthometric height to an ellipsoidal one before it is transformed, because
the transformation lives in geocentric space and knows nothing of mean sea level.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from geodetic_archive.adjustment import Adjustment, Station, adjust
from geodetic_archive.archive import ArchiveFormatError, read_archive, read_tie_points
from geodetic_archive.ellipsoid import KID62, WGS84, ecef_to_geodetic, geodetic_to_ecef
from geodetic_archive.geoid import GeoidGrid
from geodetic_archive.transform import apply_helmert

logger = logging.getLogger(__name__)

#: Input file names, fixed by the KA-71 interchange format.
ARCHIVE_FILE = "archive.dat"
TIE_POINTS_FILE = "tiepoints.csv"
GEOID_FILE = "geoid_kid62.grd"

#: Output file names, fixed by docs/brief.md.
PARAMS_FILE = "datum_params.json"
GEOJSON_FILE = "archive_wgs84.geojson"
QA_FILE = "qa_report.json"

#: Rounding applied on serialisation. Positions keep nine decimal places, which is
#: well below a millimetre of latitude and far finer than the data warrants; the
#: point is not to lose anything to the file format.
_TRANSLATION_DECIMALS = 4
_ANGLE_DECIMALS = 6
_DEGREE_DECIMALS = 9
_HEIGHT_DECIMALS = 4
_ORTHOMETRIC_DECIMALS = 3
_RESIDUAL_DECIMALS = 4


@dataclass(frozen=True)
class MigrationResult:
    """What a run produced, for a caller that wants it without re-reading the files.

    Args:
        adjustment: The recovered transformation and the evidence for it.
        migrated_record_count: Records written to the FeatureCollection.
        excluded_record_ids: Records left out for carrying a void field, sorted.
        outputs: The three files written, by their file name.
    """

    adjustment: Adjustment
    migrated_record_count: int
    excluded_record_ids: tuple[str, ...]
    outputs: dict[str, Path]


def load_connection_stations(data_dir: Path) -> list[Station]:
    """Read the connection stations and put both sides into geocentric metres.

    The legacy side goes through the legacy ellipsoid and the GNSS side through
    WGS84. Both heights in ``tiepoints.csv`` are already ellipsoidal, so the
    separation model does not apply here.
    """
    stations: list[Station] = []
    for tie in read_tie_points(data_dir / TIE_POINTS_FILE):
        latitude, longitude = tie.legacy_position()
        stations.append(
            Station(
                station_id=tie.station_id,
                legacy_xyz=geodetic_to_ecef(
                    latitude, longitude, tie.legacy_ellipsoidal_height_m, KID62
                ),
                wgs84_xyz=geodetic_to_ecef(
                    tie.wgs84_lat_deg,
                    tie.wgs84_lon_deg,
                    tie.wgs84_ellipsoidal_height_m,
                    WGS84,
                ),
                wgs84_lat_deg=tie.wgs84_lat_deg,
                wgs84_lon_deg=tie.wgs84_lon_deg,
            )
        )
    return stations


def build_feature_collection(
    data_dir: Path, adjustment: Adjustment
) -> tuple[dict[str, Any], list[str]]:
    """Transform every complete record and collect the void ones.

    Args:
        data_dir: Directory holding the register and the separation model.
        adjustment: The recovered transformation to apply.

    Returns:
        ``(feature_collection, excluded_record_ids)``.
    """
    geoid = GeoidGrid.from_file(data_dir / GEOID_FILE)
    features: list[dict[str, Any]] = []
    excluded: list[str] = []

    for record in read_archive(data_dir / ARCHIVE_FILE):
        if not record.is_complete:
            excluded.append(record.record_id)
            continue

        latitude, longitude = record.position()
        orthometric = record.orthometric_height_m()
        ellipsoidal = geoid.to_ellipsoidal(orthometric, latitude, longitude)
        legacy_xyz = geodetic_to_ecef(latitude, longitude, ellipsoidal, KID62)
        new_lat, new_lon, new_height = ecef_to_geodetic(
            apply_helmert(legacy_xyz, adjustment.parameters), WGS84
        )
        features.append(
            {
                "type": "Feature",
                "geometry": {
                    "type": "Point",
                    "coordinates": [
                        round(new_lon, _DEGREE_DECIMALS),
                        round(new_lat, _DEGREE_DECIMALS),
                        round(new_height, _HEIGHT_DECIMALS),
                    ],
                },
                "properties": {
                    "record_id": record.record_id,
                    "name": record.name,
                    "height_orthometric_m": round(orthometric, _ORTHOMETRIC_DECIMALS),
                },
            }
        )

    logger.info("migrated %d records, excluded %d for void fields", len(features), len(excluded))
    return {"type": "FeatureCollection", "features": features}, excluded


def migrate(data_dir: str | Path, out_dir: str | Path) -> MigrationResult:
    """Recover the transformation, apply it to the register, write the three outputs.

    Args:
        data_dir: Directory holding archive.dat, tiepoints.csv and geoid_kid62.grd.
        out_dir: Directory to write datum_params.json, archive_wgs84.geojson and
            qa_report.json into. Created if it does not exist.

    Returns:
        The adjustment and the record counts, so a caller need not re-read the JSON.

    Raises:
        ArchiveFormatError: If ``data_dir`` is not a directory, or an input file is
            missing or malformed.
        AdjustmentError: If the control cannot support the fit.
        GeoidGridError: If the separation model cannot be read.
    """
    source = Path(data_dir)
    destination = Path(out_dir)
    if not source.is_dir():
        raise ArchiveFormatError(f"input directory {source} does not exist")

    stations = load_connection_stations(source)
    adjustment = adjust(stations)
    collection, excluded = build_feature_collection(source, adjustment)

    destination.mkdir(parents=True, exist_ok=True)
    parameters = adjustment.parameters.to_dict()
    payload = {
        "dx_m": round(parameters["dx_m"], _TRANSLATION_DECIMALS),
        "dy_m": round(parameters["dy_m"], _TRANSLATION_DECIMALS),
        "dz_m": round(parameters["dz_m"], _TRANSLATION_DECIMALS),
        "rx_arcsec": round(parameters["rx_arcsec"], _ANGLE_DECIMALS),
        "ry_arcsec": round(parameters["ry_arcsec"], _ANGLE_DECIMALS),
        "rz_arcsec": round(parameters["rz_arcsec"], _ANGLE_DECIMALS),
        "scale_ppm": round(parameters["scale_ppm"], _ANGLE_DECIMALS),
    }
    qa_report = {
        "n_records_migrated": len(collection["features"]),
        "n_records_excluded_void": len(excluded),
        "excluded_record_ids": sorted(excluded),
        "rejected_tie_points": list(adjustment.rejected_station_ids),
        "tie_point_residuals": {
            station_id: [round(value, _RESIDUAL_DECIMALS) for value in triple]
            for station_id, triple in adjustment.residuals_enu_m.items()
        },
        "rms_residual_m": round(adjustment.rms_residual_m, _RESIDUAL_DECIMALS),
    }

    outputs = {
        PARAMS_FILE: destination / PARAMS_FILE,
        GEOJSON_FILE: destination / GEOJSON_FILE,
        QA_FILE: destination / QA_FILE,
    }
    _write_json(outputs[PARAMS_FILE], payload, indent=1)
    _write_json(outputs[GEOJSON_FILE], collection, indent=None)
    _write_json(outputs[QA_FILE], qa_report, indent=1)

    logger.info("wrote %s", ", ".join(str(path) for path in outputs.values()))
    return MigrationResult(
        adjustment=adjustment,
        migrated_record_count=len(collection["features"]),
        excluded_record_ids=tuple(sorted(excluded)),
        outputs=outputs,
    )


def _write_json(path: Path, payload: object, indent: int | None) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(payload, handle, indent=indent)
