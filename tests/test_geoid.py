"""The separation model reader and the vertical datum conversion."""

from __future__ import annotations

from pathlib import Path

import pytest

from geodetic_archive.geoid import GeoidGrid, GeoidGridError

GRID_TEXT = """\
# a two-by-two test grid
lat_min 8.0
lon_min 73.0
d_lat 1.0
d_lon 1.0
n_lat 2
n_lon 2
  0.0   10.0
 20.0   30.0
"""


@pytest.fixture
def tiny_grid(tmp_path: Path) -> GeoidGrid:
    path = tmp_path / "tiny.grd"
    path.write_text(GRID_TEXT, encoding="utf-8")
    return GeoidGrid.from_file(path)


def _write(tmp_path: Path, text: str) -> Path:
    path = tmp_path / "broken.grd"
    path.write_text(text, encoding="utf-8")
    return path


# --- the shipped model -------------------------------------------------------
def test_shipped_grid_matches_its_declared_shape(data_dir: Path) -> None:
    grid = GeoidGrid.from_file(data_dir / "geoid_kid62.grd")
    assert (grid.n_lat, grid.n_lon) == (7, 7)
    assert (grid.lat_min, grid.lon_min) == (7.0, 72.0)
    assert (grid.d_lat, grid.d_lon) == (0.5, 0.5)
    assert grid.source.endswith("geoid_kid62.grd")


def test_shipped_separations_are_the_documented_tens_of_metres(data_dir: Path) -> None:
    """Roughly -42 to -34 m across the island; this is the 38 m the vertical turns on."""
    grid = GeoidGrid.from_file(data_dir / "geoid_kid62.grd")
    flat = [value for row in grid.values for value in row]
    assert min(flat) == pytest.approx(-42.376)
    assert max(flat) == pytest.approx(-33.610)
    assert grid.separation(8.5, 73.5) < -30.0


def test_orthometric_to_ellipsoidal_applies_the_separation(data_dir: Path) -> None:
    grid = GeoidGrid.from_file(data_dir / "geoid_kid62.grd")
    assert grid.to_ellipsoidal(100.0, 8.5, 73.5) == pytest.approx(
        100.0 + grid.separation(8.5, 73.5)
    )


# --- interpolation -----------------------------------------------------------
def test_grid_nodes_interpolate_to_themselves(tiny_grid: GeoidGrid) -> None:
    assert tiny_grid.separation(8.0, 73.0) == pytest.approx(0.0)
    assert tiny_grid.separation(8.0, 74.0) == pytest.approx(10.0)
    assert tiny_grid.separation(9.0, 73.0) == pytest.approx(20.0)
    assert tiny_grid.separation(9.0, 74.0) == pytest.approx(30.0)


def test_cell_centre_is_the_mean_of_its_corners(tiny_grid: GeoidGrid) -> None:
    assert tiny_grid.separation(8.5, 73.5) == pytest.approx(15.0)


def test_interpolation_is_bilinear_not_nearest_neighbour(tiny_grid: GeoidGrid) -> None:
    assert tiny_grid.separation(8.25, 73.0) == pytest.approx(5.0)
    assert tiny_grid.separation(8.0, 73.75) == pytest.approx(7.5)


def test_positions_outside_the_grid_extrapolate_from_the_nearest_cell(
    tiny_grid: GeoidGrid,
) -> None:
    """The index is clamped, the fraction is not, so the edge gradient continues."""
    assert tiny_grid.separation(7.0, 73.0) == pytest.approx(-20.0)
    assert tiny_grid.separation(10.0, 75.0) == pytest.approx(60.0)


# --- malformed input ---------------------------------------------------------
def test_missing_file_names_the_path(tmp_path: Path) -> None:
    with pytest.raises(GeoidGridError, match="cannot read separation model"):
        GeoidGrid.from_file(tmp_path / "absent.grd")


def test_missing_header_keys_are_all_named(tmp_path: Path) -> None:
    path = _write(tmp_path, "lat_min 8.0\nlon_min 73.0\nn_lat 2\nn_lon 2\n1 2\n3 4\n")
    with pytest.raises(GeoidGridError, match="missing d_lat, d_lon"):
        GeoidGrid.from_file(path)


def test_a_non_numeric_body_names_the_line(tmp_path: Path) -> None:
    path = _write(tmp_path, GRID_TEXT.replace(" 20.0   30.0", " 20.0   thirty"))
    with pytest.raises(GeoidGridError, match="line 9: expected numbers"):
        GeoidGrid.from_file(path)


def test_a_body_that_contradicts_the_header_is_rejected(tmp_path: Path) -> None:
    path = _write(tmp_path, GRID_TEXT.replace("n_lat 2", "n_lat 3"))
    with pytest.raises(GeoidGridError, match="declares 3x2 but the body holds 2 rows"):
        GeoidGrid.from_file(path)


def test_ragged_rows_are_rejected(tmp_path: Path) -> None:
    path = _write(tmp_path, GRID_TEXT.replace(" 20.0   30.0", " 20.0   30.0   40.0"))
    with pytest.raises(GeoidGridError, match="but the body holds"):
        GeoidGrid.from_file(path)


def test_a_grid_too_small_to_interpolate_is_rejected(tmp_path: Path) -> None:
    path = _write(
        tmp_path,
        "lat_min 8.0\nlon_min 73.0\nd_lat 1.0\nd_lon 1.0\nn_lat 1\nn_lon 1\n5.0\n",
    )
    with pytest.raises(GeoidGridError, match="needs at least 2x2"):
        GeoidGrid.from_file(path)
