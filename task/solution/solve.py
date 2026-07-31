"""Migrate the KID62 geodetic archive to WGS84.

The seven-parameter transformation is unknown, so it is estimated from the
connection stations by a linearised geocentric similarity adjustment, after the
miskeyed legacy coordinates have been screened out. The recovered parameters are
then applied to the register.
"""
import csv
import json
import math
import os

import numpy as np

DATA = "/app/data"
OUT = "/app/out"

ARCSEC = math.pi / (180.0 * 3600.0)

# Reference systems, from /app/data/FORMAT.txt sections 1.1 and 1.2.
KID62_A, KID62_INVF = 6377295.664, 300.8017
WGS84_A, WGS84_INVF = 6378137.0, 298.257223563
GRID = dict(lat0=8.0, lon0=73.5, k0=0.99965, e0=250000.0, n0=150000.0)


def e2_of(a, invf):
    f = 1.0 / invf
    return 2.0 * f - f * f


KID62_E2 = e2_of(KID62_A, KID62_INVF)
WGS84_E2 = e2_of(WGS84_A, WGS84_INVF)


# --------------------------------------------------------------- geodesy core
def to_ecef(lat_deg, lon_deg, h, a, e2):
    lat, lon = math.radians(lat_deg), math.radians(lon_deg)
    sl, cl = math.sin(lat), math.cos(lat)
    n = a / math.sqrt(1.0 - e2 * sl * sl)
    return np.array([(n + h) * cl * math.cos(lon),
                     (n + h) * cl * math.sin(lon),
                     (n * (1.0 - e2) + h) * sl])


def to_geodetic(x, y, z, a, e2):
    p = math.hypot(x, y)
    lon = math.atan2(y, x)
    lat = math.atan2(z, p * (1.0 - e2))
    h = 0.0
    for _ in range(12):
        sl = math.sin(lat)
        n = a / math.sqrt(1.0 - e2 * sl * sl)
        h = p / math.cos(lat) - n
        lat = math.atan2(z, p * (1.0 - e2 * n / (n + h)))
    return math.degrees(lat), math.degrees(lon), h


def helmert(xyz, p):
    """Position-vector convention: X2 = T + (1+s) R X1."""
    rx, ry, rz = p[3] * ARCSEC, p[4] * ARCSEC, p[5] * ARCSEC
    k = 1.0 + p[6] * 1e-6
    x, y, z = xyz
    return np.array([p[0] + k * (x - rz * y + ry * z),
                     p[1] + k * (rz * x + y - rx * z),
                     p[2] + k * (-ry * x + rx * y + z)])


def enu_basis(lat_deg, lon_deg):
    lat, lon = math.radians(lat_deg), math.radians(lon_deg)
    sp, cp, sl, cl = math.sin(lat), math.cos(lat), math.sin(lon), math.cos(lon)
    return np.array([[-sl, cl, 0.0],
                     [-sp * cl, -sp * sl, cp],
                     [cp * cl, cp * sl, sp]])


# ------------------------------------------- Transverse Mercator on KID62
def _tm_consts(lat):
    sl = math.sin(lat)
    w = 1.0 - KID62_E2 * sl * sl
    nu = KID62_A * GRID["k0"] / math.sqrt(w)
    rho = KID62_A * GRID["k0"] * (1.0 - KID62_E2) / w ** 1.5
    return nu, rho, nu / rho - 1.0


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


def grid_to_geodetic(east, north):
    lat = math.radians(GRID["lat0"])
    ak0 = KID62_A * GRID["k0"]
    for _ in range(30):
        dn = north - GRID["n0"] - _mer_arc(lat)
        if abs(dn) < 1e-9:
            break
        lat += dn / ak0
    nu, rho, eta2 = _tm_consts(lat)
    t, sec = math.tan(lat), 1.0 / math.cos(lat)
    vii = t / (2 * rho * nu)
    viii = t / (24 * rho * nu**3) * (5 + 3 * t**2 + eta2 - 9 * t**2 * eta2)
    ix = t / (720 * rho * nu**5) * (61 + 90 * t**2 + 45 * t**4)
    x1 = sec / nu
    xi = sec / (6 * nu**3) * (nu / rho + 2 * t**2)
    xii = sec / (120 * nu**5) * (5 + 28 * t**2 + 24 * t**4)
    xiia = sec / (5040 * nu**7) * (61 + 662 * t**2 + 1320 * t**4 + 720 * t**6)
    de = east - GRID["e0"]
    return (math.degrees(lat - vii * de**2 + viii * de**4 - ix * de**6),
            math.degrees(math.radians(GRID["lon0"])
                         + x1 * de - xi * de**3 + xii * de**5 - xiia * de**7))


# ------------------------------------------------------- separation model
class Geoid:
    def __init__(self, path):
        hdr, rows = {}, []
        with open(path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split()
                if len(parts) == 2 and parts[0][0].isalpha():
                    hdr[parts[0]] = float(parts[1])
                else:
                    rows.append([float(v) for v in parts])
        self.lat0, self.lon0 = hdr["lat_min"], hdr["lon_min"]
        self.dlat, self.dlon = hdr["d_lat"], hdr["d_lon"]
        self.nlat, self.nlon = int(hdr["n_lat"]), int(hdr["n_lon"])
        self.v = rows
        assert len(rows) == self.nlat and all(len(r) == self.nlon for r in rows)

    def sep(self, lat, lon):
        fi = (lat - self.lat0) / self.dlat
        fj = (lon - self.lon0) / self.dlon
        i = min(max(int(math.floor(fi)), 0), self.nlat - 2)
        j = min(max(int(math.floor(fj)), 0), self.nlon - 2)
        di, dj = fi - i, fj - j
        return ((1 - di) * (1 - dj) * self.v[i][j] + (1 - di) * dj * self.v[i][j + 1]
                + di * (1 - dj) * self.v[i + 1][j] + di * dj * self.v[i + 1][j + 1])


# ------------------------------------------------------------ legacy decoding
def decode_coord(field, kind):
    """Decode a KA-71 coordinate field (FORMAT.txt section 2.1) -> KID62 lat, lon."""
    if kind == "G":
        a, b = field[0:10], field[11:22]
        lat = int(a[0:2]) + int(a[2:4]) / 60.0 + int(a[4:9]) / 3.6e6
        lon = int(b[0:3]) + int(b[3:5]) / 60.0 + int(b[5:10]) / 3.6e6
        return (-lat if a[9] == "S" else lat), (-lon if b[10] == "W" else lon)
    if kind == "D":
        return float(field[0:11]), float(field[12:24])
    if kind == "P":
        return grid_to_geodetic(int(field[0:12]) / 1000.0, int(field[13:25]) / 1000.0)
    raise ValueError(kind)


def is_void(field):
    s = field.strip()
    return bool(s) and set(s) == {"9"}


def read_archive(path):
    """Assemble records across continuation lines (FORMAT.txt section 2.3)."""
    records, order = {}, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            kind, rid = line[0], line[1:8]
            if kind == "+":
                records[rid]["name"] += line[8:32]
                continue
            records[rid] = dict(kind=kind, rid=rid, name=line[8:32],
                                coord=line[32:62], height=line[62:71])
            order.append(rid)
    for rid in order:
        records[rid]["name"] = records[rid]["name"].rstrip()
    return [records[r] for r in order]


# ------------------------------------------------------------- the adjustment
def fit_helmert(src, dst):
    """Gauss-Newton least squares for the 7 position-vector parameters."""
    p = np.zeros(7)
    for _ in range(4):
        a = np.zeros((3 * len(src), 7))
        b = np.zeros(3 * len(src))
        for i, (s, d) in enumerate(zip(src, dst)):
            x, y, z = s
            r = 3 * i
            a[r + 0, 0] = a[r + 1, 1] = a[r + 2, 2] = 1.0
            a[r + 0, 4], a[r + 0, 5] = z * ARCSEC, -y * ARCSEC
            a[r + 1, 3], a[r + 1, 5] = -z * ARCSEC, x * ARCSEC
            a[r + 2, 3], a[r + 2, 4] = y * ARCSEC, -x * ARCSEC
            a[r + 0, 6], a[r + 1, 6], a[r + 2, 6] = x * 1e-6, y * 1e-6, z * 1e-6
            b[r:r + 3] = d - helmert(s, p)
        p = p + np.linalg.lstsq(a, b, rcond=None)[0]
    return p


# Least squares has no defence against gross errors. One station kilometres out of
# place spreads itself across every residual, so screening the residuals of a fit
# that still holds the blunders fingers the wrong stations. Cut them beforehand
# instead, on a statistic that never sees a least squares solution.
SCREEN_FLOOR_M = 5.0
# The post-fit cut has to come from the residuals in hand, not from a remembered
# number. The smallest transcription errors here leave residuals of a metre or two,
# which a rule-of-thumb cut in the low metres sails straight past, while the honest
# residuals sit near a decimetre. Scaling to the robust spread of the residuals
# themselves puts the cut just above the honest floor wherever that floor happens to
# be, and the blunders stand at eight sigma and beyond.
SNOOP_SIGMA = 4.0
SNOOP_FLOOR_M = 0.5


def screen_blunders(stations):
    """Cut stations whose raw geocentric difference vector leaves the cluster.

    Over a region this small the rotation and scale terms only move a station a few
    metres relative to its neighbours, so the difference vectors of the honest
    stations sit tightly around their component-wise median. A miskeyed coordinate
    throws its difference vector out by whatever the keying error was worth, which
    is orders of magnitude further. The cut comes from the spread of the difference
    vectors themselves, via a median absolute deviation the blunders cannot shift,
    with a floor at the level honest scatter can reach.
    """
    d = np.array([s["dst"] - s["src"] for s in stations])
    norms = np.linalg.norm(d - np.median(d, axis=0), axis=1)
    mad = 1.4826 * float(np.median(np.abs(norms - np.median(norms))))
    cut = max(float(np.median(norms)) + 8.0 * mad, SCREEN_FLOOR_M)
    return [s for s, v in zip(stations, norms) if v <= cut], cut


def robust_fit(stations):
    """Screen gross errors, adjust, then snoop the residuals for anything left."""
    keep, cut = screen_blunders(stations)
    while True:
        src = [s["src"] for s in keep]
        dst = [s["dst"] for s in keep]
        p = fit_helmert(src, dst)
        norms = np.array([np.linalg.norm(d - helmert(s, p)) for s, d in zip(src, dst)])
        med = float(np.median(norms))
        sigma = 1.4826 * float(np.median(np.abs(norms - med)))
        snoop = max(med + SNOOP_SIGMA * sigma, SNOOP_FLOOR_M)
        worst = int(np.argmax(norms))
        if len(keep) > 10 and norms[worst] > snoop:
            keep.pop(worst)
            continue
        return p, keep, cut


def main():
    geoid = Geoid(os.path.join(DATA, "geoid_kid62.grd"))

    stations = []
    with open(os.path.join(DATA, "tiepoints.csv")) as fh:
        for row in csv.DictReader(fh):
            lat, lon = decode_coord(row["legacy_coord_field"].ljust(30),
                                    row["legacy_encoding"])
            h = int(row["legacy_ellip_height_mm"]) / 1000.0
            stations.append(dict(
                sid=row["station_id"],
                src=to_ecef(lat, lon, h, KID62_A, KID62_E2),
                dst=to_ecef(float(row["wgs84_lat_deg"]), float(row["wgs84_lon_deg"]),
                            float(row["wgs84_ellip_height_m"]), WGS84_A, WGS84_E2),
                wlat=float(row["wgs84_lat_deg"]), wlon=float(row["wgs84_lon_deg"])))

    params, kept, screen_cut = robust_fit(stations)
    kept_ids = {s["sid"] for s in kept}
    rejected = sorted(s["sid"] for s in stations if s["sid"] not in kept_ids)

    residuals, sq = {}, 0.0
    for s in kept:
        d = s["dst"] - helmert(s["src"], params)
        enu = enu_basis(s["wlat"], s["wlon"]) @ d
        residuals[s["sid"]] = [round(float(v), 4) for v in enu]
        sq += float(d @ d)
    rms = math.sqrt(sq / len(kept))

    features, excluded = [], []
    for rec in read_archive(os.path.join(DATA, "archive.dat")):
        if is_void(rec["coord"]) or is_void(rec["height"]):
            excluded.append(rec["rid"])
            continue
        lat, lon = decode_coord(rec["coord"], rec["kind"])
        h_orth = int(rec["height"]) / (100.0 if rec["kind"] == "D" else 10.0)
        h_ell = h_orth + geoid.sep(lat, lon)
        wlat, wlon, wh = to_geodetic(*helmert(to_ecef(lat, lon, h_ell, KID62_A, KID62_E2),
                                              params), WGS84_A, WGS84_E2)
        features.append(dict(
            type="Feature",
            geometry=dict(type="Point",
                          coordinates=[round(wlon, 9), round(wlat, 9), round(wh, 4)]),
            properties=dict(record_id=rec["rid"], name=rec["name"],
                            height_orthometric_m=round(h_orth, 3))))

    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "datum_params.json"), "w") as fh:
        json.dump(dict(dx_m=round(params[0], 4), dy_m=round(params[1], 4),
                       dz_m=round(params[2], 4), rx_arcsec=round(params[3], 6),
                       ry_arcsec=round(params[4], 6), rz_arcsec=round(params[5], 6),
                       scale_ppm=round(params[6], 6)), fh, indent=1)
    with open(os.path.join(OUT, "archive_wgs84.geojson"), "w") as fh:
        json.dump(dict(type="FeatureCollection", features=features), fh)
    with open(os.path.join(OUT, "qa_report.json"), "w") as fh:
        json.dump(dict(n_records_migrated=len(features),
                       n_records_excluded_void=len(excluded),
                       excluded_record_ids=sorted(excluded),
                       rejected_tie_points=rejected,
                       tie_point_residuals=residuals,
                       rms_residual_m=round(rms, 4)), fh, indent=1)

    print("parameters :", np.round(params, 4).tolist())
    print("screen cut : %.2f m" % screen_cut)
    print("rejected   :", rejected)
    print("rms        : %.4f m" % rms)
    print("migrated   :", len(features), "excluded:", len(excluded))


if __name__ == "__main__":
    main()
