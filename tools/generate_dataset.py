"""Dataset generator for dynamo/migrate-geodetic-archive.

Deliberately outside task/: this is authoring tooling, it is never built into the
image, and the task does not reference it. It is kept in the repo so a reviewer can
see where tests/ground_truth.json came from rather than take 1170 hardcoded records
on trust.

Writes the agent-visible inputs into task/environment/data/ and the verifier's
ground truth into task/tests/ground_truth.json. Fixed seed, so regenerating
reproduces both byte for byte:

    python tools/generate_dataset.py .

The true transformation parameters live in tools/geodesy.py. Regenerating the
inputs regenerates the ground truth alongside them, so the two cannot drift apart.
"""
import json
import math
import os
import random
import sys

import geodesy as g

SEED = 20260730
N_TIE = 34
N_REC = 1200
N_VOID_COORD = 20
N_VOID_HEIGHT = 10
P_LONG_NAME = 0.10

LAT_LO, LAT_HI = 7.62, 9.88
LON_LO, LON_HI = 72.44, 74.86

rnd = random.Random(SEED)

# ---------------------------------------------------------------- geoid model
GRID_LAT0, GRID_LON0, GRID_D, GRID_N = 7.0, 72.0, 0.5, 7


def _geoid_analytic(lat, lon):
    return (-36.0
            + 2.8 * math.sin(0.9 * (lat - 8.5))
            - 3.4 * math.cos(0.7 * (lon - 73.5))
            + 0.6 * (lat - 8.5) * (lon - 73.5))


GEOID = [[round(_geoid_analytic(GRID_LAT0 + i * GRID_D, GRID_LON0 + j * GRID_D), 3)
          for j in range(GRID_N)] for i in range(GRID_N)]


def geoid_sep(lat, lon):
    """Bilinear interpolation of the shipped separation grid."""
    fi = (lat - GRID_LAT0) / GRID_D
    fj = (lon - GRID_LON0) / GRID_D
    i = min(max(int(math.floor(fi)), 0), GRID_N - 2)
    j = min(max(int(math.floor(fj)), 0), GRID_N - 2)
    di, dj = fi - i, fj - j
    return ((1 - di) * (1 - dj) * GEOID[i][j] + (1 - di) * dj * GEOID[i][j + 1]
            + di * (1 - dj) * GEOID[i + 1][j] + di * dj * GEOID[i + 1][j + 1])


# ------------------------------------------------------------------- encoders
def _dms_parts(v):
    d = int(v)
    m = int((v - d) * 60)
    si = int(round((v - d - m / 60.0) * 3600.0 * 1000))
    if si >= 60000:
        si -= 60000
        m += 1
    if m >= 60:
        m -= 60
        d += 1
    return d, m, si


def enc_g(lat, lon):
    d, m, s = _dms_parts(abs(lat))
    latf = "%02d%02d%05d%s" % (d, m, s, "N" if lat >= 0 else "S")
    d, m, s = _dms_parts(abs(lon))
    lonf = "%03d%02d%05d%s" % (d, m, s, "E" if lon >= 0 else "W")
    return (latf + " " + lonf).ljust(30)


def dec_g(field):
    a, b = field[0:10], field[11:22]
    lat = int(a[0:2]) + int(a[2:4]) / 60.0 + int(a[4:9]) / 3.6e6
    lon = int(b[0:3]) + int(b[3:5]) / 60.0 + int(b[5:10]) / 3.6e6
    return (-lat if a[9] == "S" else lat), (-lon if b[10] == "W" else lon)


def enc_d(lat, lon):
    return ("%+011.7f %+012.7f" % (lat, lon)).ljust(30)


def dec_d(field):
    return float(field[0:11]), float(field[12:24])


def enc_p(lat, lon):
    e, n = g.tm_forward(lat, lon)
    return ("%012d %012d" % (round(e * 1000), round(n * 1000))).ljust(30)


def dec_p(field):
    return g.tm_inverse(int(field[0:12]) / 1000.0, int(field[13:25]) / 1000.0)


DECODE = {"G": dec_g, "D": dec_d, "P": dec_p}


# -------------------------------------------------------------------- names
PLACES = ["KAVERI", "PONMUDI", "ALANGAD", "THEKKADY", "VARKALA", "MANJERI",
          "KUTTANAD", "ANJENGO", "KOLLAM", "NEDUMBASSERY", "PALLIPURAM",
          "CHERTHALA", "MUNNAR", "IDUKKI", "ARUVIKKARA", "VAIKOM",
          "KARUNAGAPALLY", "MAVELIKKARA", "THIRUVALLA", "CHENGANNUR",
          "PANDALAM", "ADOOR", "KOTTARAKKARA", "PUNALUR", "ANCHAL",
          "NEDUMANGAD", "ATTINGAL", "PARAVUR", "CHIRAYINKEEZHU", "VELLANAD"]
KINDS = ["TRIG PILLAR", "BENCH MARK", "AZIMUTH MARK", "TIDAL STATION",
         "REVENUE STONE", "SURVEY PILLAR", "TOWER STATION",
         "ECCENTRIC STATION", "BASE TERMINAL"]


def make_name(long_wanted):
    """Compose a station name; long names never split on a space at column 24."""
    while True:
        nm = "%s %s" % (rnd.choice(PLACES), rnd.choice(KINDS))
        if rnd.random() < 0.45:
            nm += " NO %d" % rnd.randint(2, 97)
        if (len(nm) > 24) != long_wanted:
            continue
        if long_wanted and (nm[23] == " " or nm[24] == " "):
            continue
        return nm


# ----------------------------------------------------------------- tie points
# Gross keying errors (km scale) and two quiet ones. The gross errors are what
# wreck a contaminated least squares fit, so they have to stay: they are the reason
# screening post-fit residuals identifies the wrong stations. The quiet pair sits
# below a typical robust pre-screen and only shows up in the residuals of a fit
# that is already clean, so a solver who screens once and stops carries one into
# the adjustment and lands an RMS well above the decimetre the GNSS side promises.
BLUNDERS = ["TP007", "TP011", "TP016", "TP022", "TP024", "TP027", "TP030"]


def tie_encoding(k):
    return "P" if k % 3 == 0 else "G"


def gen_tie_points():
    rows = []
    for k in range(1, N_TIE + 1):
        sid = "TP%03d" % k
        lat = rnd.uniform(LAT_LO, LAT_HI)
        lon = rnd.uniform(LON_LO, LON_HI)
        h_kid = rnd.uniform(-5.0, 2180.0)
        enc = tie_encoding(k)
        field = enc_g(lat, lon) if enc == "G" else enc_p(lat, lon)
        rlat, rlon = DECODE[enc](field)
        h_mm = round(h_kid * 1000)
        wlat, wlon, wh = g.kid62_to_wgs84(rlat, rlon, h_mm / 1000.0)
        wlat += rnd.gauss(0, 0.08) / 110900.0
        wlon += rnd.gauss(0, 0.08) / (111320.0 * math.cos(math.radians(rlat)))
        wh += rnd.gauss(0, 0.15)
        rows.append([sid, enc, field.rstrip(), "%+010d" % h_mm,
                     "%+.8f" % wlat, "%+.8f" % wlon, "%+.3f" % wh])
    return rows


def apply_blunders(rows):
    """Inject realistic transcription errors into the legacy side of BLUNDERS."""
    by_id = {r[0]: r for r in rows}
    report = {}

    for sid in BLUNDERS:
        r = by_id[sid]
        enc, f = r[1], r[2]
        before = DECODE[enc](f.ljust(30))
        if sid == "TP007":                        # seconds digits transposed
            f = f[:4] + f[6] + f[5] + f[4] + f[7:]
        elif sid == "TP011":                      # milliarcsecond digits transposed
            f = f[:6] + f[8] + f[7] + f[6] + f[9:]
        elif sid == "TP027":                      # single-digit slip, metres column
            f = "%012d" % (int(f[0:12]) + 8000) + f[12:]
        elif sid == "TP016":                      # longitude hemisphere flipped
            f = f[:21] + "W" + f[22:]
        elif sid == "TP022":                      # decimal degrees keyed into DMS slot
            f = ("%09d" % round(before[0] * 1e7))[:9] + f[9:]
        elif sid == "TP024":                      # easting/northing swapped
            f = f[13:25] + " " + f[0:12]
        elif sid == "TP030":                      # +1 km typo in the easting
            f = "%012d" % (int(f[0:12]) + 1_000_000) + f[12:]
        r[2] = f
        after = DECODE[enc](f.ljust(30))
        d = math.hypot((after[0] - before[0]) * 110900.0,
                       (after[1] - before[1]) * 110000.0)
        report[sid] = d
        assert enc == tie_encoding(int(sid[2:])), sid
    return report


# -------------------------------------------------------------------- archive
def gen_archive():
    lines, truth, names, void_ids = [], {}, {}, []
    void_coord = set(rnd.sample(range(1, N_REC + 1), N_VOID_COORD))
    void_height = set(rnd.sample([i for i in range(1, N_REC + 1)
                                  if i not in void_coord], N_VOID_HEIGHT))

    for k in range(1, N_REC + 1):
        rid = "KA%05d" % k
        name = make_name(rnd.random() < P_LONG_NAME)
        lat = rnd.uniform(LAT_LO, LAT_HI)
        lon = rnd.uniform(LON_LO, LON_HI)
        H = rnd.uniform(-3.0, 2210.0)
        rtype = "GDP"[k % 3]

        if k in void_coord:
            field = "9" * 30
        else:
            field = {"G": enc_g, "D": enc_d, "P": enc_p}[rtype](lat, lon)

        if k in void_height:
            hfield = "9" * 9
        else:
            hfield = "%+09d" % round(H * (100 if rtype == "D" else 10))

        lines.append(rtype + rid + name[:24].ljust(24) + field + hfield)
        if len(name) > 24:
            lines.append("+" + rid + name[24:].ljust(24) + " " * 39)

        if k in void_coord or k in void_height:
            void_ids.append(rid)
            continue

        rlat, rlon = DECODE[rtype](field)
        rH = int(hfield) / (100.0 if rtype == "D" else 10.0)
        wlat, wlon, wh = g.kid62_to_wgs84(rlat, rlon, rH + geoid_sep(rlat, rlon))
        truth[rid] = [round(wlon, 9), round(wlat, 9), round(wh, 4), round(rH, 3)]
        names[rid] = name
    return lines, truth, names, void_ids


# ----------------------------------------------------------------------- main
def main(root):
    data_dir = os.path.join(root, "task", "environment", "data")
    tests_dir = os.path.join(root, "task", "tests")
    os.makedirs(data_dir, exist_ok=True)

    tie_rows = gen_tie_points()
    blunder_report = apply_blunders(tie_rows)
    with open(os.path.join(data_dir, "tiepoints.csv"), "w", newline="\n") as fh:
        fh.write("station_id,legacy_encoding,legacy_coord_field,legacy_ellip_height_mm,"
                 "wgs84_lat_deg,wgs84_lon_deg,wgs84_ellip_height_m\n")
        for r in tie_rows:
            fh.write(",".join(r) + "\n")

    lines, truth, names, void_ids = gen_archive()
    with open(os.path.join(data_dir, "archive.dat"), "w", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")

    with open(os.path.join(data_dir, "geoid_kid62.grd"), "w", newline="\n") as fh:
        fh.write("# Kaveri Island geoid-ellipsoid separation model KIG61\n")
        fh.write("# N = h(ellipsoidal, KID62) - H(orthometric, Kaveri MSL 1961), metres\n")
        fh.write("lat_min %.1f\nlat_max %.1f\nlon_min %.1f\nlon_max %.1f\n"
                 % (GRID_LAT0, GRID_LAT0 + (GRID_N - 1) * GRID_D,
                    GRID_LON0, GRID_LON0 + (GRID_N - 1) * GRID_D))
        fh.write("d_lat %.1f\nd_lon %.1f\nn_lat %d\nn_lon %d\n"
                 % (GRID_D, GRID_D, GRID_N, GRID_N))
        fh.write("# one row per latitude, ascending from lat_min;\n")
        fh.write("# columns ascending in longitude from lon_min\n")
        for row in GEOID:
            fh.write(" ".join("%9.3f" % v for v in row) + "\n")

    with open(os.path.join(tests_dir, "ground_truth.json"), "w", newline="\n") as fh:
        json.dump(dict(true_params=g.TRUE_PARAMS,
                       blunder_tie_points=sorted(BLUNDERS),
                       void_record_ids=sorted(void_ids),
                       expected_record_count=len(truth),
                       records=truth,
                       names=names), fh, indent=1, sort_keys=True)

    print("archive lines       :", len(lines))
    print("records migratable  :", len(truth))
    print("void records        :", len(void_ids))
    print("continuation lines  :", len(lines) - N_REC)
    print("geoid range         : %.3f .. %.3f"
          % (min(min(r) for r in GEOID), max(max(r) for r in GEOID)))
    for sid, d in sorted(blunder_report.items()):
        print("blunder %s displaced %11.1f m" % (sid, d))


if __name__ == "__main__":
    main(sys.argv[1])
