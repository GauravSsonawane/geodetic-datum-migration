"""Command line entry point: ``geodetic-archive``."""

from __future__ import annotations

import argparse
import logging
import sys
from collections.abc import Sequence
from pathlib import Path

from geodetic_archive import __version__
from geodetic_archive.adjustment import AdjustmentError
from geodetic_archive.archive import ArchiveFormatError
from geodetic_archive.geoid import GeoidGridError
from geodetic_archive.pipeline import migrate

#: Errors that mean the inputs or the control are unusable. They carry a message
#: written for a person, so the CLI prints that rather than a traceback.
_EXPECTED_ERRORS = (ArchiveFormatError, AdjustmentError, GeoidGridError)

DEFAULT_DATA_DIR = "data"
DEFAULT_OUT_DIR = "out"

EXIT_OK = 0
EXIT_FAILED = 1


def build_parser() -> argparse.ArgumentParser:
    """The argument parser, exposed so the tests can exercise it directly."""
    parser = argparse.ArgumentParser(
        prog="geodetic-archive",
        description=(
            "Recover the KID62 to WGS84 transformation from the archive's own "
            "connection stations and migrate the station register onto WGS84."
        ),
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    migrate_parser = subparsers.add_parser(
        "migrate",
        help="recover the transformation and migrate the register",
        description=(
            "Writes datum_params.json, archive_wgs84.geojson and qa_report.json "
            "into the output directory."
        ),
    )
    migrate_parser.add_argument(
        "--data",
        type=Path,
        default=Path(DEFAULT_DATA_DIR),
        metavar="DIR",
        help="directory holding archive.dat, tiepoints.csv and geoid_kid62.grd "
        f"(default: {DEFAULT_DATA_DIR})",
    )
    migrate_parser.add_argument(
        "--out",
        type=Path,
        default=Path(DEFAULT_OUT_DIR),
        metavar="DIR",
        help=f"directory to write the three outputs into (default: {DEFAULT_OUT_DIR})",
    )
    migrate_parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="verbosity on stderr (default: INFO)",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run the CLI.

    Args:
        argv: Arguments to parse. Defaults to ``sys.argv[1:]``.

    Returns:
        ``0`` on success, ``1`` if the inputs or the control turned out unusable.
    """
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(levelname)-7s %(name)s: %(message)s",
        stream=sys.stderr,
    )

    try:
        result = migrate(args.data, args.out)
    except _EXPECTED_ERRORS as exc:
        print(f"geodetic-archive: {exc}", file=sys.stderr)
        return EXIT_FAILED

    adjustment = result.adjustment
    print(f"parameters : {adjustment.parameters.to_dict(4)}")
    print(f"adjusted   : {len(adjustment.used_station_ids)} stations")
    print(f"rejected   : {', '.join(adjustment.rejected_station_ids) or 'none'}")
    print(f"rms        : {adjustment.rms_residual_m:.4f} m")
    print(
        f"migrated   : {result.migrated_record_count} records, "
        f"{len(result.excluded_record_ids)} excluded for void fields"
    )
    print(f"written    : {args.out}")
    return EXIT_OK


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
