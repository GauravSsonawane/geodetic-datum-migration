"""Verifier for dynamo/migrate-geodetic-archive.

Grades the three files instruction.md asks for against ground truth held only in
/tests. The bands are justified in task.toml under verification_explanation. Each
was checked to admit other defensible adjustments (different tie point subsets,
different interpolation of the separation grid) and to reject the methodological
failures the task is built around: wrong ellipsoid on the legacy side, translation
only shift, separation grid ignored, blunders left in the fit, and rotations
reported under the other sign convention.
"""
import json
import math
import os

import numpy as np
import pytest

OUT = "/app/out"
DATA = "/app/data"
GROUND_TRUTH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "ground_truth.json")

ARCSEC = math.pi / (180.0 * 3600.0)
KID62_A, KID62_INVF = 6377295.664, 300.8017
WGS84_A, WGS84_INVF = 6378137.0, 298.257223563
GRID = dict(lat0=8.0, lon0=73.5, k0=0.99965, e0=250000.0, n0=150000.0)

# --- tolerances --------------------------------------------------------------
TOL_TRANSLATION_M = 10.0     # formal sigma of the adjustment is ~2.2-2.7 m
TOL_ROTATION_ARCSEC = 0.35   # formal sigma ~0.06-0.08"
TOL_SCALE_PPM = 1.0          # formal sigma ~0.27 ppm
TOL_HORIZ_RMS_M = 0.30
TOL_HORIZ_MAX_M = 0.75
TOL_VERT_RMS_M = 1.50
TOL_VERT_MAX_M = 4.00
TOL_ORTHO_M = 0.05           # the archived height, only decoded and carried through
TOL_RESIDUAL_M = 0.05        # self-consistency of the reported residuals
RMS_BAND_M = (0.02, 1.00)    # a decimetre-quality GNSS side cannot fit better/worse

PARAM_KEYS = ["dx_m", "dy_m", "dz_m", "rx_arcsec", "ry_arcsec", "rz_arcsec", "scale_ppm"]
TRUTH_KEYS = ["dx", "dy", "dz", "rx", "ry", "rz", "s"]


# --- fixtures ----------------------------------------------------------------
def _load_json(path):
    if not os.path.isfile(path):
        pytest.fail("required output file %s was not produced" % path)
    with open(path) as fh:
        try:
            return json.load(fh)
        except json.JSONDecodeError as exc:
            pytest.fail("%s is not valid JSON: %s" % (path, exc))


@pytest.fixture(scope="module")
def truth():
    with open(GROUND_TRUTH) as fh:
        return json.load(fh)


@pytest.fixture(scope="module")
def params():
    return _load_json(os.path.join(OUT, "datum_params.json"))


@pytest.fixture(scope="module")
def qa():
    return _load_json(os.path.join(OUT, "qa_report.json"))


@pytest.fixture(scope="module")
def features():
    fc = _load_json(os.path.join(OUT, "archive_wgs84.geojson"))
    assert isinstance(fc, dict) and fc.get("type") == "FeatureCollection", \
        "archive_wgs84.geojson must be a GeoJSON FeatureCollection"
    feats = fc.get("features")
    assert isinstance(feats, list) and feats, "the FeatureCollection carries no features"
    return feats


@pytest.fixture(scope="module")
def by_id(features):
    out = {}
    for f in features:
        rid = f.get("properties", {}).get("record_id")
        assert rid is not None, "a feature is missing its record_id property"
        assert rid not in out, "record %s appears more than once" % rid
        out[rid] = f
    return out


# --- geodesy needed to re-derive the agent's own numbers ---------------------
def _e2(a, invf):
    f = 1.0 / invf
    return 2.0 * f - f * f


KID62_E2, WGS84_E2 = _e2(KID62_A, KID62_INVF), _e2(WGS84_A, WGS84_INVF)


def _to_ecef(lat_deg, lon_deg, h, a, e2):
    lat, lon = math.radians(lat_deg), math.radians(lon_deg)
    sl, cl = math.sin(lat), math.cos(lat)
    n = a / math.sqrt(1.0 - e2 * sl * sl)
    return np.array([(n + h) * cl * math.cos(lon), (n + h) * cl * math.sin(lon),
                     (n * (1.0 - e2) + h) * sl])


def _helmert(xyz, p):
    rx, ry, rz = (p["rx_arcsec"] * ARCSEC, p["ry_arcsec"] * ARCSEC, p["rz_arcsec"] * ARCSEC)
    k = 1.0 + p["scale_ppm"] * 1e-6
    x, y, z = xyz
    return np.array([p["dx_m"] + k * (x - rz * y + ry * z),
                     p["dy_m"] + k * (rz * x + y - rx * z),
                     p["dz_m"] + k * (-ry * x + rx * y + z)])


def _mer_arc(lat):
    b = KID62_A * (1.0 - 1.0 / KID62_INVF)
    n = (KID62_A - b) / (KID62_A + b)
    lat0 = math.radians(GRID["lat0"])
    d, s = lat - lat0, lat + lat0
    return b * GRID["k0"] * (
        (1 + n + 1.25 * n**2 + 1.25 * n**3) * d
        - (3 * n + 3 * n**2 + 2.625 * n**3) * math.sin(d) * math.cos(s)
        + (1.875 * n**2 + 1.875 * n**3) * math.sin(2 * d) * math.cos(2 * s)
        - (35.0 / 24.0) * n**3 * math.sin(3 * d) * math.cos(3 * s))


def _grid_to_geodetic(east, north):
    lat = math.radians(GRID["lat0"])
    for _ in range(30):
        dn = north - GRID["n0"] - _mer_arc(lat)
        if abs(dn) < 1e-9:
            break
        lat += dn / (KID62_A * GRID["k0"])
    sl = math.sin(lat)
    w = 1.0 - KID62_E2 * sl * sl
    nu = KID62_A * GRID["k0"] / math.sqrt(w)
    rho = KID62_A * GRID["k0"] * (1.0 - KID62_E2) / w ** 1.5
    eta2, t, sec = nu / rho - 1.0, math.tan(lat), 1.0 / math.cos(lat)
    de = east - GRID["e0"]
    lat_f = (lat
             - t / (2 * rho * nu) * de**2
             + t / (24 * rho * nu**3) * (5 + 3 * t**2 + eta2 - 9 * t**2 * eta2) * de**4
             - t / (720 * rho * nu**5) * (61 + 90 * t**2 + 45 * t**4) * de**6)
    lon_f = (math.radians(GRID["lon0"])
             + sec / nu * de
             - sec / (6 * nu**3) * (nu / rho + 2 * t**2) * de**3
             + sec / (120 * nu**5) * (5 + 28 * t**2 + 24 * t**4) * de**5
             - sec / (5040 * nu**7) * (61 + 662 * t**2 + 1320 * t**4 + 720 * t**6) * de**7)
    return math.degrees(lat_f), math.degrees(lon_f)


def _decode_legacy(field, kind):
    field = field.ljust(30)
    if kind == "G":
        a, b = field[0:10], field[11:22]
        lat = int(a[0:2]) + int(a[2:4]) / 60.0 + int(a[4:9]) / 3.6e6
        lon = int(b[0:3]) + int(b[3:5]) / 60.0 + int(b[5:10]) / 3.6e6
        return (-lat if a[9] == "S" else lat), (-lon if b[10] == "W" else lon)
    return _grid_to_geodetic(int(field[0:12]) / 1000.0, int(field[13:25]) / 1000.0)


@pytest.fixture(scope="module")
def tie_stations():
    """The connection stations as shipped, indexed by identifier."""
    import csv
    out = {}
    with open(os.path.join(DATA, "tiepoints.csv")) as fh:
        for row in csv.DictReader(fh):
            lat, lon = _decode_legacy(row["legacy_coord_field"], row["legacy_encoding"])
            out[row["station_id"]] = dict(
                src=_to_ecef(lat, lon, int(row["legacy_ellip_height_mm"]) / 1000.0,
                             KID62_A, KID62_E2),
                dst=_to_ecef(float(row["wgs84_lat_deg"]), float(row["wgs84_lon_deg"]),
                             float(row["wgs84_ellip_height_m"]), WGS84_A, WGS84_E2))
    return out


def _horizontal_m(lat, lon, tlat, tlon):
    return math.hypot((lat - tlat) * 110900.0,
                      (lon - tlon) * 111320.0 * math.cos(math.radians(tlat)))


# --- datum_params.json -------------------------------------------------------
def test_params_schema(params):
    """datum_params.json carries exactly the seven documented numeric keys."""
    assert isinstance(params, dict), "datum_params.json must be a JSON object"
    missing = [k for k in PARAM_KEYS if k not in params]
    assert not missing, "datum_params.json is missing keys: %s" % missing
    for k in PARAM_KEYS:
        assert isinstance(params[k], (int, float)) and not isinstance(params[k], bool), \
            "%s must be a number, got %r" % (k, params[k])


def test_translation_parameters(params, truth):
    """The recovered geocentric translation matches the true KID62 to WGS84 shift."""
    for pk, tk in zip(PARAM_KEYS[:3], TRUTH_KEYS[:3]):
        err = abs(params[pk] - truth["true_params"][tk])
        assert err <= TOL_TRANSLATION_M, \
            "%s off by %.3f m (tolerance %.1f m)" % (pk, err, TOL_TRANSLATION_M)


def test_rotation_parameters(params, truth):
    """Rotations match the truth in the position-vector convention the prompt names."""
    for pk, tk in zip(PARAM_KEYS[3:6], TRUTH_KEYS[3:6]):
        err = abs(params[pk] - truth["true_params"][tk])
        assert err <= TOL_ROTATION_ARCSEC, \
            "%s off by %.4f arcsec (tolerance %.2f); note the sign convention" \
            % (pk, err, TOL_ROTATION_ARCSEC)


def test_scale_parameter(params, truth):
    """The recovered scale difference matches the true value."""
    err = abs(params["scale_ppm"] - truth["true_params"]["s"])
    assert err <= TOL_SCALE_PPM, \
        "scale_ppm off by %.4f ppm (tolerance %.1f)" % (err, TOL_SCALE_PPM)


# --- qa_report.json ----------------------------------------------------------
def test_blunder_stations_rejected(qa, truth):
    """Exactly the mis-transcribed connection stations are excluded from the fit."""
    assert "rejected_tie_points" in qa, "qa_report.json is missing rejected_tie_points"
    got = set(qa["rejected_tie_points"])
    expected = set(truth["blunder_tie_points"])
    assert got == expected, \
        "rejected %s; missed %s and wrongly rejected %s" \
        % (sorted(got), sorted(expected - got), sorted(got - expected))


def test_residuals_reported_for_kept_stations(qa, truth, tie_stations):
    """Residuals are reported for every station kept, and for none that was rejected."""
    res = qa.get("tie_point_residuals")
    assert isinstance(res, dict) and res, "tie_point_residuals must be a non-empty object"
    kept = set(tie_stations) - set(truth["blunder_tie_points"])
    assert set(res) == kept, \
        "residuals cover %s; expected the %d kept stations (missing %s, extra %s)" \
        % (len(res), len(kept), sorted(kept - set(res)), sorted(set(res) - kept))
    for sid, triple in res.items():
        assert isinstance(triple, (list, tuple)) and len(triple) == 3, \
            "residual for %s must be an east, north, up triple" % sid


def test_residuals_consistent_with_parameters(qa, params, tie_stations):
    """Reported residuals are the ones the reported parameters actually produce.

    Compared as vector magnitudes, so either residual sign convention is accepted.
    """
    for sid, triple in qa["tie_point_residuals"].items():
        st = tie_stations[sid]
        expect = float(np.linalg.norm(st["dst"] - _helmert(st["src"], params)))
        got = float(np.linalg.norm(np.asarray(triple, dtype=float)))
        assert abs(got - expect) <= TOL_RESIDUAL_M, \
            "residual for %s is %.4f m but the reported parameters give %.4f m" \
            % (sid, got, expect)


def test_rms_residual(qa):
    """rms_residual_m is the 3D RMS of the reported residuals, at GNSS accuracy."""
    triples = [np.asarray(v, dtype=float) for v in qa["tie_point_residuals"].values()]
    expect = math.sqrt(sum(float(v @ v) for v in triples) / len(triples))
    got = qa.get("rms_residual_m")
    assert isinstance(got, (int, float)), "rms_residual_m must be a number"
    assert abs(got - expect) <= 1e-3 + 1e-2 * expect, \
        "rms_residual_m is %.4f but the reported residuals give %.4f" % (got, expect)
    assert RMS_BAND_M[0] <= got <= RMS_BAND_M[1], \
        "an adjustment of decimetre-quality GNSS ties cannot have a 3D RMS of %.4f m" % got


def test_void_records_excluded(qa, truth, by_id):
    """Records carrying a void field are reported as excluded and are not migrated."""
    expected = set(truth["void_record_ids"])
    got = set(qa.get("excluded_record_ids", []))
    assert got == expected, \
        "excluded %d records; missed %s and wrongly excluded %s" \
        % (len(got), sorted(expected - got), sorted(got - expected))
    assert qa.get("n_records_excluded_void") == len(expected), \
        "n_records_excluded_void is %r, expected %d" \
        % (qa.get("n_records_excluded_void"), len(expected))
    leaked = expected & set(by_id)
    assert not leaked, "void records were migrated anyway: %s" % sorted(leaked)


def test_migrated_record_count(qa, truth, by_id):
    """Every complete record appears once, and the reported count agrees."""
    expected = set(truth["records"])
    assert set(by_id) == expected, \
        "migrated %d records; missing %s, unexpected %s" \
        % (len(by_id), sorted(expected - set(by_id))[:5], sorted(set(by_id) - expected)[:5])
    assert qa.get("n_records_migrated") == len(expected), \
        "n_records_migrated is %r but %d records were expected" \
        % (qa.get("n_records_migrated"), len(expected))


# --- archive_wgs84.geojson ---------------------------------------------------
def test_feature_structure(features):
    """Each feature is a Point with three ordinates and the documented properties."""
    for f in features:
        assert f.get("type") == "Feature", "every entry must be a GeoJSON Feature"
        geom = f.get("geometry") or {}
        assert geom.get("type") == "Point", "every geometry must be a Point"
        coords = geom.get("coordinates")
        assert isinstance(coords, (list, tuple)) and len(coords) == 3, \
            "coordinates must be [longitude, latitude, ellipsoidal_height]"
        for v in coords:
            assert isinstance(v, (int, float)) and not isinstance(v, bool), \
                "coordinates must be numbers, got %r" % (coords,)
        props = f.get("properties") or {}
        for key in ("record_id", "name", "height_orthometric_m"):
            assert key in props, "feature %r is missing property %s" % (props, key)


def test_horizontal_positions(by_id, truth):
    """Migrated longitudes and latitudes match the true WGS84 positions."""
    errs = []
    for rid, f in by_id.items():
        lon, lat, _ = f["geometry"]["coordinates"][:3]
        tlon, tlat = truth["records"][rid][0], truth["records"][rid][1]
        errs.append(_horizontal_m(lat, lon, tlat, tlon))
    rms = math.sqrt(sum(e * e for e in errs) / len(errs))
    assert rms <= TOL_HORIZ_RMS_M, \
        "horizontal RMS error %.3f m exceeds %.2f m" % (rms, TOL_HORIZ_RMS_M)
    assert max(errs) <= TOL_HORIZ_MAX_M, \
        "worst horizontal error %.3f m exceeds %.2f m" % (max(errs), TOL_HORIZ_MAX_M)


def test_ellipsoidal_heights(by_id, truth):
    """Third ordinates are WGS84 ellipsoidal heights, i.e. the separation model was
    applied to the archived orthometric heights."""
    errs = [abs(f["geometry"]["coordinates"][2] - truth["records"][rid][2])
            for rid, f in by_id.items()]
    rms = math.sqrt(sum(e * e for e in errs) / len(errs))
    assert rms <= TOL_VERT_RMS_M, \
        "ellipsoidal height RMS error %.3f m exceeds %.2f m" % (rms, TOL_VERT_RMS_M)
    assert max(errs) <= TOL_VERT_MAX_M, \
        "worst ellipsoidal height error %.3f m exceeds %.2f m" % (max(errs), TOL_VERT_MAX_M)


def test_orthometric_heights_carried_through(by_id, truth):
    """height_orthometric_m reproduces the archived height, decoded in the units the
    record type declares."""
    bad = []
    for rid, f in by_id.items():
        got = f["properties"]["height_orthometric_m"]
        if not isinstance(got, (int, float)) or abs(got - truth["records"][rid][3]) > TOL_ORTHO_M:
            bad.append((rid, got, truth["records"][rid][3]))
    assert not bad, "%d records have a wrong archived height, e.g. %s" % (len(bad), bad[:4])


def test_station_names(by_id, truth):
    """Station names are complete, including those carried onto continuation lines."""
    norm = lambda s: " ".join(str(s).split())
    bad = [(rid, f["properties"]["name"], truth["names"][rid])
           for rid, f in by_id.items()
           if norm(f["properties"]["name"]) != norm(truth["names"][rid])]
    assert not bad, "%d records have a wrong station name, e.g. %s" % (len(bad), bad[:3])
