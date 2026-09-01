"""The KA-71 register reader, its three encodings and its void markers."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

from geodetic_archive.archive import (
    ArchiveFormatError,
    ArchiveRecord,
    decode_coordinate,
    is_void,
    read_archive,
    read_tie_points,
)

# Three encodings of one position, laid out as FORMAT.txt section 2.1 declares.
SEXAGESIMAL = "091207020N 0744459239E"
DECIMAL = "+08.0697278 +073.3055309"
PLANE = "000361739023 000257797238"


def _record(**overrides: Any) -> ArchiveRecord:
    fields: dict[str, Any] = {
        "record_id": "KA00001",
        "encoding": "G",
        "name": "KAVERI TRIG PILLAR",
        "coordinate_field": SEXAGESIMAL.ljust(30),
        "height_field": "+00010392",
        "line_number": 1,
    }
    fields.update(overrides)
    return ArchiveRecord(**fields)


# --- void markers ------------------------------------------------------------
@pytest.mark.parametrize("field", ["9" * 30, "9" * 9, "  999  "])
def test_void_markers_are_recognised(field: str) -> None:
    assert is_void(field)


@pytest.mark.parametrize("field", ["", "   ", "+00099999", "9999999X", SEXAGESIMAL])
def test_ordinary_fields_are_not_void(field: str) -> None:
    """A height that merely contains nines is an observation, not a void marker."""
    assert not is_void(field)


def test_a_record_is_incomplete_if_either_field_is_void() -> None:
    assert _record().is_complete
    assert not _record(coordinate_field="9" * 30).is_complete
    assert not _record(height_field="9" * 9).is_complete


# --- coordinate decoding -----------------------------------------------------
def test_sexagesimal_decoding_unpacks_degrees_minutes_and_seconds() -> None:
    latitude, longitude = decode_coordinate(SEXAGESIMAL, "G")
    assert latitude == pytest.approx(9 + 12 / 60 + 7.020 / 3600)
    assert longitude == pytest.approx(74 + 44 / 60 + 59.239 / 3600)


def test_hemisphere_letters_set_the_sign() -> None:
    southern = SEXAGESIMAL.replace("N", "S")
    western = SEXAGESIMAL.replace("E", "W")
    assert decode_coordinate(southern, "G")[0] < 0
    assert decode_coordinate(western, "G")[1] < 0


def test_decimal_degree_decoding() -> None:
    latitude, longitude = decode_coordinate(DECIMAL, "D")
    assert latitude == pytest.approx(8.0697278)
    assert longitude == pytest.approx(73.3055309)


def test_plane_decoding_runs_the_inverse_projection() -> None:
    """Millimetre integers on the Kaveri Grid, back through the legacy ellipsoid."""
    latitude, longitude = decode_coordinate(PLANE, "P")
    assert 7.5 < latitude < 10.0
    assert 72.0 < longitude < 75.0


def test_a_trimmed_field_decodes_the_same_as_a_padded_one() -> None:
    assert decode_coordinate(SEXAGESIMAL, "G") == decode_coordinate(SEXAGESIMAL.ljust(30), "G")


def test_an_unknown_encoding_is_named_in_the_error() -> None:
    with pytest.raises(ArchiveFormatError, match="unknown coordinate encoding 'X'"):
        decode_coordinate(SEXAGESIMAL, "X")


@pytest.mark.parametrize(
    ("field", "encoding"),
    [("not a coordinate", "G"), ("abc.def ghi.jkl", "D"), ("00036173xxxx 000257797238", "P")],
)
def test_a_field_that_defies_its_layout_is_rejected(field: str, encoding: str) -> None:
    with pytest.raises(ArchiveFormatError, match="does not match its layout"):
        decode_coordinate(field, encoding)


# --- height units ------------------------------------------------------------
def test_height_is_decimetres_on_type_g_and_p() -> None:
    assert _record(encoding="G", height_field="+00010392").orthometric_height_m() == pytest.approx(
        1039.2
    )
    assert _record(
        encoding="P", coordinate_field=PLANE.ljust(30), height_field="+00007723"
    ).orthometric_height_m() == pytest.approx(772.3)


def test_height_is_centimetres_on_type_d() -> None:
    """Reading every record in decimetres is a factor of ten on a third of them."""
    record = _record(encoding="D", coordinate_field=DECIMAL.ljust(30), height_field="+00193667")
    assert record.orthometric_height_m() == pytest.approx(1936.67)


def test_negative_heights_survive_the_sign() -> None:
    assert _record(height_field="-00000031").orthometric_height_m() == pytest.approx(-3.1)


def test_a_void_height_is_reported_with_its_record_and_line() -> None:
    with pytest.raises(ArchiveFormatError, match="KA00001 .line 1. has a void height"):
        _record(height_field="9" * 9).orthometric_height_m()


def test_a_void_coordinate_is_reported_with_its_record_and_line() -> None:
    with pytest.raises(ArchiveFormatError, match="KA00001 .line 1. has a void coordinate"):
        _record(coordinate_field="9" * 30).position()


def test_a_non_integer_height_is_rejected() -> None:
    with pytest.raises(ArchiveFormatError, match="is not a signed integer"):
        _record(height_field="+0001.39").orthometric_height_m()


# --- the shipped register ----------------------------------------------------
def test_shipped_register_holds_the_expected_populations(data_dir: Path) -> None:
    records = read_archive(data_dir / "archive.dat")
    assert len(records) == 1200
    assert sum(record.is_complete for record in records) == 1170
    assert sum(not record.is_complete for record in records) == 30


def test_record_identifiers_are_unique_and_in_file_order(data_dir: Path) -> None:
    ids = [record.record_id for record in read_archive(data_dir / "archive.dat")]
    assert len(set(ids)) == len(ids)
    assert ids == sorted(ids)


def test_all_three_encodings_appear(data_dir: Path) -> None:
    encodings = {record.encoding for record in read_archive(data_dir / "archive.dat")}
    assert encodings == {"G", "D", "P"}


def test_long_names_are_assembled_across_continuation_lines(data_dir: Path) -> None:
    """A continuation line is not a record; its 24 characters extend the name."""
    records = read_archive(data_dir / "archive.dat")
    long_names = [record.name for record in records if len(record.name) > 24]
    assert long_names, "the register is expected to exercise continuation lines"
    assert "THEKKADY ECCENTRIC STATION" in long_names
    assert not any(name != name.rstrip() for name in long_names)


def test_every_complete_record_decodes(data_dir: Path) -> None:
    for record in read_archive(data_dir / "archive.dat"):
        if not record.is_complete:
            continue
        latitude, longitude = record.position()
        assert 7.0 < latitude < 10.5
        assert 72.0 < longitude < 75.5
        assert -10.0 < record.orthometric_height_m() < 2500.0


# --- malformed registers -----------------------------------------------------
def test_a_missing_register_names_the_path(tmp_path: Path) -> None:
    with pytest.raises(ArchiveFormatError, match="cannot read station register"):
        read_archive(tmp_path / "absent.dat")


def test_an_unknown_record_type_names_the_line(tmp_path: Path) -> None:
    path = tmp_path / "archive.dat"
    path.write_text("XKA00001NAME\n", encoding="utf-8")
    with pytest.raises(ArchiveFormatError, match="line 1: unknown record type 'X'"):
        read_archive(path)


def test_an_orphan_continuation_line_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "archive.dat"
    path.write_text("+KA00099TRAILING\n", encoding="utf-8")
    with pytest.raises(ArchiveFormatError, match="continuation for unknown record 'KA00099'"):
        read_archive(path)


def test_a_repeated_identifier_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "archive.dat"
    line = "GKA00001NAME\n"
    path.write_text(line + line, encoding="utf-8")
    with pytest.raises(ArchiveFormatError, match="'KA00001' appears more than once"):
        read_archive(path)


def test_blank_lines_are_ignored(tmp_path: Path) -> None:
    path = tmp_path / "archive.dat"
    path.write_text("\nGKA00001NAME\n\n", encoding="utf-8")
    assert len(read_archive(path)) == 1


# --- connection stations -----------------------------------------------------
def test_shipped_connection_stations(data_dir: Path) -> None:
    stations = read_tie_points(data_dir / "tiepoints.csv")
    assert len(stations) == 34
    assert stations[0].station_id == "TP001"
    assert {station.encoding for station in stations} == {"G", "P"}


def test_legacy_tie_heights_are_millimetres_and_already_ellipsoidal(data_dir: Path) -> None:
    """No separation model applies to either side of a connection station."""
    stations = read_tie_points(data_dir / "tiepoints.csv")
    assert stations[0].legacy_ellipsoidal_height_m == pytest.approx(550.493)
    assert all(-10.0 < station.legacy_ellipsoidal_height_m < 2300.0 for station in stations)


def test_every_legacy_tie_coordinate_decodes(data_dir: Path) -> None:
    """Including the miskeyed ones: a blunder is a plausible field, just wrong."""
    for station in read_tie_points(data_dir / "tiepoints.csv"):
        latitude, longitude = station.legacy_position()
        assert -180.0 <= longitude <= 180.0
        assert -90.0 <= latitude <= 90.0


def test_a_missing_tie_file_names_the_path(tmp_path: Path) -> None:
    with pytest.raises(ArchiveFormatError, match="cannot read connection stations"):
        read_tie_points(tmp_path / "absent.csv")


def test_missing_tie_columns_are_all_named(tmp_path: Path) -> None:
    path = tmp_path / "tiepoints.csv"
    path.write_text("station_id,legacy_encoding\nTP001,G\n", encoding="utf-8")
    with pytest.raises(ArchiveFormatError, match="missing column"):
        read_tie_points(path)


def test_a_non_numeric_tie_column_names_the_line(tmp_path: Path, data_dir: Path) -> None:
    source = (data_dir / "tiepoints.csv").read_text(encoding="utf-8").splitlines()
    broken = source[1].replace("+000550493", "not-a-height")
    path = tmp_path / "tiepoints.csv"
    path.write_text("\n".join([source[0], broken]) + "\n", encoding="utf-8")
    with pytest.raises(ArchiveFormatError, match="line 2"):
        read_tie_points(path)
