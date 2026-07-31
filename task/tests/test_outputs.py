"""Verifier for dynamo/migrate-geodetic-archive.

Grades the three files instruction.md asks for against ground truth held only in
/tests. The bands are justified in task.toml under verification_explanation. Each
was checked to admit other defensible adjustments (a different subset of the sound
stations carried into the fit, different interpolation of the separation grid) and
to reject the methodological failures the task is built around: wrong ellipsoid on
the legacy side, translation only shift, separation grid ignored, blunders left in
the fit, and rotations reported under the other sign convention.
"""
import json
import math
import os

import numpy as np
import pytest

OUT = "/app/out"
GROUND_TRUTH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "ground_truth.json")

ARCSEC = math.pi / (180.0 * 3600.0)
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
MIN_ADJUSTED_STATIONS = 10   # seven parameters need real redundancy behind them

PARAM_KEYS = ["dx_m", "dy_m", "dz_m", "rx_arcsec", "ry_arcsec", "rz_arcsec", "scale_ppm"]
TRUTH_KEYS = ["dx", "dy", "dz", "rx", "ry", "rz", "s"]


# --- fixtures ----------------------------------------------------------------
def _finite(v):
    """A real, finite number. bool is not a number here, and NaN is not a value.

    NaN matters: every comparison against it is False, so a check written as
    "flag it when the error exceeds the tolerance" silently accepts NaN. Tests
    below either use this helper or assert the passing condition directly.
    """
    return (isinstance(v, (int, float)) and not isinstance(v, bool)
            and math.isfinite(v))


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
def _helmert(xyz, p):
    rx, ry, rz = (p["rx_arcsec"] * ARCSEC, p["ry_arcsec"] * ARCSEC, p["rz_arcsec"] * ARCSEC)
    k = 1.0 + p["scale_ppm"] * 1e-6
    x, y, z = xyz
    return np.array([p["dx_m"] + k * (x - rz * y + ry * z),
                     p["dy_m"] + k * (rz * x + y - rx * z),
                     p["dz_m"] + k * (-ry * x + rx * y + z)])


@pytest.fixture(scope="module")
def tie_stations(truth):
    """The connection stations as recorded when the archive was built.

    Taken from ground truth rather than re-read from /app/data/tiepoints.csv. That
    file sits on the agent's own writable filesystem, so certifying the reported
    residuals against it would certify them against a file the agent could have
    rewritten. Nothing the verifier grades is read from a path the agent can touch,
    apart from the outputs it is asked to produce.
    """
    return {sid: dict(src=np.asarray(v["src"], dtype=float),
                      dst=np.asarray(v["dst"], dtype=float),
                      wlat=v["wlat"], wlon=v["wlon"])
            for sid, v in truth["tie_stations"].items()}


def _enu_basis(lat_deg, lon_deg):
    """Rows are the east, north and up unit vectors at the station."""
    lat, lon = math.radians(lat_deg), math.radians(lon_deg)
    sp, cp, sl, cl = math.sin(lat), math.cos(lat), math.sin(lon), math.cos(lon)
    return np.array([[-sl, cl, 0.0],
                     [-sp * cl, -sp * sl, cp],
                     [cp * cl, cp * sl, sp]])


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
        assert _finite(params[k]), "%s must be a finite number, got %r" % (k, params[k])


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
    """Exactly the mis-transcribed connection stations are identified as such.

    This is an objective property of the archive, not a judgement call: a station's
    legacy coordinate either was miskeyed or was not. Which sound stations to carry
    into the adjustment is a separate, freer choice, checked by the residual test.
    """
    assert "rejected_tie_points" in qa, "qa_report.json is missing rejected_tie_points"
    got = set(qa["rejected_tie_points"])
    expected = set(truth["blunder_tie_points"])
    assert got == expected, \
        "rejected %s; missed %s and wrongly rejected %s" \
        % (sorted(got), sorted(expected - got), sorted(got - expected))


def test_residuals_reported_for_adjusted_stations(qa, truth, tie_stations):
    """Residuals cover the stations the adjustment was computed from.

    Which sound stations to carry into the adjustment is the solver's call, so this
    does not demand all of them: dropping a few extra as a precaution is ordinary
    practice. What it does demand is that no mis-transcribed station was used, and
    that enough were kept for a seven-parameter fit to mean anything.
    """
    res = qa.get("tie_point_residuals")
    assert isinstance(res, dict) and res, "tie_point_residuals must be a non-empty object"
    sound = set(tie_stations) - set(truth["blunder_tie_points"])
    used_bad = set(res) & set(truth["blunder_tie_points"])
    assert not used_bad, \
        "residuals reported for mis-transcribed stations, which cannot have been " \
        "used in the adjustment: %s" % sorted(used_bad)
    unknown = set(res) - set(tie_stations)
    assert not unknown, "residuals reported for unknown stations: %s" % sorted(unknown)
    assert len(res) >= MIN_ADJUSTED_STATIONS, \
        "only %d stations carried into the adjustment; a seven-parameter fit needs " \
        "at least %d to be meaningful" % (len(res), MIN_ADJUSTED_STATIONS)
    assert set(res) <= sound  # implied by the two checks above; stated for clarity
    for sid, triple in res.items():
        assert isinstance(triple, (list, tuple)) and len(triple) == 3, \
            "residual for %s must be an east, north, up triple" % sid


def test_residuals_consistent_with_parameters(qa, params, tie_stations):
    """Reported east, north, up residuals are the ones the reported parameters
    actually produce at that station.

    Compared as whole vectors, allowing an overall sign flip so either direction
    convention is accepted. Comparing magnitudes alone would not do: any triple of
    the right length would pass, including a fabricated one that puts the whole
    residual on a single axis.
    """
    for sid, triple in qa["tie_point_residuals"].items():
        st = tie_stations[sid]
        expect = _enu_basis(st["wlat"], st["wlon"]) @ (st["dst"] - _helmert(st["src"], params))
        got = np.asarray(triple, dtype=float)
        err = min(float(np.linalg.norm(got - expect)), float(np.linalg.norm(got + expect)))
        assert err <= TOL_RESIDUAL_M, \
            "residual for %s is %s but the reported parameters give %s (off by %.4f m)" \
            % (sid, np.round(got, 4).tolist(), np.round(expect, 4).tolist(), err)


def test_rms_residual(qa):
    """rms_residual_m is the 3D RMS of the reported residuals, at GNSS accuracy."""
    triples = [np.asarray(v, dtype=float) for v in qa["tie_point_residuals"].values()]
    expect = math.sqrt(sum(float(v @ v) for v in triples) / len(triples))
    got = qa.get("rms_residual_m")
    assert _finite(got), "rms_residual_m must be a finite number, got %r" % (got,)
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
            assert _finite(v), "coordinates must be finite numbers, got %r" % (coords,)
        props = f.get("properties") or {}
        for key in ("record_id", "name", "height_orthometric_m"):
            assert key in props, "feature %r is missing property %s" % (props, key)
        assert _finite(props["height_orthometric_m"]), \
            "height_orthometric_m must be a finite number, got %r" \
            % (props["height_orthometric_m"],)


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
        want = truth["records"][rid][3]
        if not (_finite(got) and abs(got - want) <= TOL_ORTHO_M):
            bad.append((rid, got, want))
    assert not bad, "%d records have a wrong archived height, e.g. %s" % (len(bad), bad[:4])


def test_station_names(by_id, truth):
    """Station names are complete, including those carried onto continuation lines."""
    norm = lambda s: " ".join(str(s).split())
    bad = [(rid, f["properties"]["name"], truth["names"][rid])
           for rid, f in by_id.items()
           if norm(f["properties"]["name"]) != norm(truth["names"][rid])]
    assert not bad, "%d records have a wrong station name, e.g. %s" % (len(bad), bad[:3])
