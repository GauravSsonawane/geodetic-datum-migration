"""The end-to-end migration and the command line that drives it."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import pytest

from geodetic_archive.archive import ArchiveFormatError
from geodetic_archive.cli import EXIT_FAILED, EXIT_OK, build_parser, main
from geodetic_archive.geoid import GeoidGridError
from geodetic_archive.pipeline import (
    GEOJSON_FILE,
    PARAMS_FILE,
    QA_FILE,
    MigrationResult,
    migrate,
)


@pytest.fixture(scope="session")
def migration(data_dir: Path, tmp_path_factory: pytest.TempPathFactory) -> MigrationResult:
    """One full run over the shipped register, shared by the assertions below."""
    return migrate(data_dir, tmp_path_factory.mktemp("out"))


@pytest.fixture(scope="session")
def outputs(migration: MigrationResult) -> dict[str, Any]:
    """The three written files, parsed."""
    return {
        name: json.loads(path.read_text(encoding="utf-8"))
        for name, path in migration.outputs.items()
    }


def _horizontal_m(lat: float, lon: float, true_lat: float, true_lon: float) -> float:
    return math.hypot(
        (lat - true_lat) * 110_900.0,
        (lon - true_lon) * 111_320.0 * math.cos(math.radians(true_lat)),
    )


# --- what a run produces -----------------------------------------------------
def test_all_three_outputs_are_written(migration: MigrationResult) -> None:
    assert set(migration.outputs) == {PARAMS_FILE, GEOJSON_FILE, QA_FILE}
    assert all(path.is_file() for path in migration.outputs.values())


def test_the_output_directory_is_created_if_absent(data_dir: Path, tmp_path: Path) -> None:
    destination = tmp_path / "nested" / "out"
    result = migrate(data_dir, destination)
    assert destination.is_dir()
    assert result.migrated_record_count == 1170


def test_record_counts_split_the_register(migration: MigrationResult) -> None:
    assert migration.migrated_record_count == 1170
    assert len(migration.excluded_record_ids) == 30
    assert migration.migrated_record_count + len(migration.excluded_record_ids) == 1200


def test_parameters_carry_exactly_the_documented_keys(outputs: dict[str, Any]) -> None:
    assert set(outputs[PARAMS_FILE]) == {
        "dx_m",
        "dy_m",
        "dz_m",
        "rx_arcsec",
        "ry_arcsec",
        "rz_arcsec",
        "scale_ppm",
    }
    assert all(math.isfinite(value) for value in outputs[PARAMS_FILE].values())


def test_the_geojson_is_a_point_feature_collection(outputs: dict[str, Any]) -> None:
    collection = outputs[GEOJSON_FILE]
    assert collection["type"] == "FeatureCollection"
    assert len(collection["features"]) == 1170
    for feature in collection["features"][:20]:
        assert feature["geometry"]["type"] == "Point"
        assert len(feature["geometry"]["coordinates"]) == 3
        assert set(feature["properties"]) == {"record_id", "name", "height_orthometric_m"}


def test_coordinates_are_longitude_latitude_height_in_that_order(
    outputs: dict[str, Any],
) -> None:
    """GeoJSON is x, y, z; getting this backwards is silent and total."""
    for feature in outputs[GEOJSON_FILE]["features"][:50]:
        longitude, latitude, _ = feature["geometry"]["coordinates"]
        assert 72.0 < longitude < 75.5
        assert 7.0 < latitude < 10.5


def test_qa_report_carries_the_documented_keys(outputs: dict[str, Any]) -> None:
    assert set(outputs[QA_FILE]) == {
        "n_records_migrated",
        "n_records_excluded_void",
        "excluded_record_ids",
        "rejected_tie_points",
        "tie_point_residuals",
        "rms_residual_m",
    }


def test_qa_counts_agree_with_the_geojson(outputs: dict[str, Any]) -> None:
    report = outputs[QA_FILE]
    assert report["n_records_migrated"] == len(outputs[GEOJSON_FILE]["features"])
    assert report["n_records_excluded_void"] == len(report["excluded_record_ids"])


# --- accuracy against ground truth -------------------------------------------
def test_the_rejected_stations_are_the_miskeyed_ones(
    outputs: dict[str, Any], ground_truth: dict[str, Any]
) -> None:
    assert outputs[QA_FILE]["rejected_tie_points"] == ground_truth["blunder_tie_points"]


def test_the_excluded_records_are_the_void_ones(
    outputs: dict[str, Any], ground_truth: dict[str, Any]
) -> None:
    assert outputs[QA_FILE]["excluded_record_ids"] == ground_truth["void_record_ids"]


@pytest.mark.slow
def test_migrated_positions_match_ground_truth(
    outputs: dict[str, Any], ground_truth: dict[str, Any]
) -> None:
    """0.30 m RMS and 0.75 m worst case are the acceptance bands; see docs/design.md."""
    truth = ground_truth["records"]
    features = {
        feature["properties"]["record_id"]: feature for feature in outputs[GEOJSON_FILE]["features"]
    }
    assert set(features) == set(truth)

    errors = []
    for record_id, feature in features.items():
        true_lon, true_lat, _, _ = truth[record_id]
        longitude, latitude, _ = feature["geometry"]["coordinates"]
        errors.append(_horizontal_m(latitude, longitude, true_lat, true_lon))

    assert math.sqrt(sum(e**2 for e in errors) / len(errors)) < 0.30
    assert max(errors) < 0.75


@pytest.mark.slow
def test_migrated_heights_match_ground_truth(
    outputs: dict[str, Any], ground_truth: dict[str, Any]
) -> None:
    """Skipping the separation model puts every height 38 m out; the band is 1.50 m."""
    truth = ground_truth["records"]
    errors = []
    orthometric_errors = []
    for feature in outputs[GEOJSON_FILE]["features"]:
        record_id = feature["properties"]["record_id"]
        _, _, true_height, true_orthometric = truth[record_id]
        errors.append(feature["geometry"]["coordinates"][2] - true_height)
        orthometric_errors.append(feature["properties"]["height_orthometric_m"] - true_orthometric)

    assert math.sqrt(sum(e**2 for e in errors) / len(errors)) < 1.50
    assert max(abs(e) for e in errors) < 4.00
    assert max(abs(e) for e in orthometric_errors) < 0.05


@pytest.mark.slow
def test_station_names_survive_the_continuation_lines(
    outputs: dict[str, Any], ground_truth: dict[str, Any]
) -> None:
    names = ground_truth["names"]
    for feature in outputs[GEOJSON_FILE]["features"]:
        record_id = feature["properties"]["record_id"]
        assert " ".join(feature["properties"]["name"].split()) == " ".join(names[record_id].split())


# --- failure paths -----------------------------------------------------------
def test_a_missing_input_directory_is_reported(tmp_path: Path) -> None:
    with pytest.raises(ArchiveFormatError, match="does not exist"):
        migrate(tmp_path / "absent", tmp_path / "out")


def test_a_directory_without_the_connection_stations_is_reported(tmp_path: Path) -> None:
    empty = tmp_path / "data"
    empty.mkdir()
    with pytest.raises(ArchiveFormatError, match="cannot read connection stations"):
        migrate(empty, tmp_path / "out")


def test_a_directory_without_the_separation_model_is_reported(
    tmp_path: Path, data_dir: Path
) -> None:
    partial = tmp_path / "data"
    partial.mkdir()
    for name in ("tiepoints.csv", "archive.dat"):
        (partial / name).write_bytes((data_dir / name).read_bytes())
    with pytest.raises(GeoidGridError, match="cannot read separation model"):
        migrate(partial, tmp_path / "out")


# --- the command line --------------------------------------------------------
def test_the_parser_defaults_to_data_and_out() -> None:
    args = build_parser().parse_args(["migrate"])
    assert (str(args.data), str(args.out)) == ("data", "out")
    assert args.log_level == "INFO"


def test_a_command_is_required() -> None:
    with pytest.raises(SystemExit) as exit_info:
        build_parser().parse_args([])
    assert exit_info.value.code == 2


def test_the_cli_runs_the_migration(
    data_dir: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    destination = tmp_path / "out"
    code = main(
        ["migrate", "--data", str(data_dir), "--out", str(destination), "--log-level", "WARNING"]
    )
    assert code == EXIT_OK

    printed = capsys.readouterr().out
    assert "adjusted   : 24 stations" in printed
    assert "migrated   : 1170 records" in printed
    assert (destination / GEOJSON_FILE).is_file()


def test_the_cli_reports_a_bad_input_directory_without_a_traceback(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    code = main(["migrate", "--data", str(tmp_path / "absent"), "--out", str(tmp_path / "out")])
    assert code == EXIT_FAILED
    assert "geodetic-archive:" in capsys.readouterr().err


def test_the_cli_reports_its_version(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as exit_info:
        main(["--version"])
    assert exit_info.value.code == 0
    assert "geodetic-archive" in capsys.readouterr().out
