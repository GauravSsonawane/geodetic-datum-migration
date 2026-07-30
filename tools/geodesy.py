"""Geodesy primitives and the true datum parameters used to build the dataset.

Authoring tooling, outside task/ and never built into the image. See
tools/generate_dataset.py."""
import math

KID62 = dict(a=6377295.664, invf=300.8017)
WGS84 = dict(a=6378137.0, invf=298.257223563)

TM = dict(lat0=8.0, lon0=73.5, k0=0.99965, E0=250000.0, N0=150000.0)

# True KID62 -> WGS84 Helmert, position-vector (Bursa-Wolf / EPSG 9606) convention.
TRUE_PARAMS = dict(
    dx=317.812, dy=708.176, dz=234.291,     # metres
    rx=-1.8734, ry=2.4171, rz=-0.9265,      # arc-seconds
    s=2.8137,                               # ppm
)

ARCSEC = math.pi / (180.0 * 3600.0)


def e2_of(ell):
    f = 1.0 / ell["invf"]
    return 2.0 * f - f * f


def geodetic_to_ecef(lat_deg, lon_deg, h, ell):
    a, e2 = ell["a"], e2_of(ell)
    lat, lon = math.radians(lat_deg), math.radians(lon_deg)
    sl, cl = math.sin(lat), math.cos(lat)
    n = a / math.sqrt(1.0 - e2 * sl * sl)
    return (
        (n + h) * cl * math.cos(lon),
        (n + h) * cl * math.sin(lon),
        (n * (1.0 - e2) + h) * sl,
    )


def ecef_to_geodetic(x, y, z, ell):
    a, e2 = ell["a"], e2_of(ell)
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


def helmert_pv(x, y, z, p):
    """Position-vector convention: X2 = T + (1+s)*R*X1."""
    rx, ry, rz = p["rx"] * ARCSEC, p["ry"] * ARCSEC, p["rz"] * ARCSEC
    k = 1.0 + p["s"] * 1e-6
    return (
        p["dx"] + k * (x - rz * y + ry * z),
        p["dy"] + k * (rz * x + y - rx * z),
        p["dz"] + k * (-ry * x + rx * y + z),
    )


def kid62_to_wgs84(lat, lon, h):
    x, y, z = geodetic_to_ecef(lat, lon, h, KID62)
    return ecef_to_geodetic(*helmert_pv(x, y, z, TRUE_PARAMS), WGS84)


# --- Transverse Mercator (Redfearn series) on the KID62 ellipsoid -------------

def _tm_consts(lat):
    a, e2 = KID62["a"], e2_of(KID62)
    k0 = TM["k0"]
    sl = math.sin(lat)
    nu = a * k0 / math.sqrt(1.0 - e2 * sl * sl)
    rho = a * k0 * (1.0 - e2) / (1.0 - e2 * sl * sl) ** 1.5
    eta2 = nu / rho - 1.0
    return nu, rho, eta2


def _meridional_arc(lat):
    a = KID62["a"]
    b = a * (1.0 - 1.0 / KID62["invf"])
    n = (a - b) / (a + b)
    lat0 = math.radians(TM["lat0"])
    dl, sl = lat - lat0, lat + lat0
    return b * TM["k0"] * (
        (1 + n + 1.25 * n**2 + 1.25 * n**3) * dl
        - (3 * n + 3 * n**2 + 2.625 * n**3) * math.sin(dl) * math.cos(sl)
        + (1.875 * n**2 + 1.875 * n**3) * math.sin(2 * dl) * math.cos(2 * sl)
        - (35.0 / 24.0) * n**3 * math.sin(3 * dl) * math.cos(3 * sl)
    )


def tm_forward(lat_deg, lon_deg):
    lat, lon = math.radians(lat_deg), math.radians(lon_deg)
    nu, rho, eta2 = _tm_consts(lat)
    m = _meridional_arc(lat)
    sl, cl, t = math.sin(lat), math.cos(lat), math.tan(lat)
    dl = lon - math.radians(TM["lon0"])
    i = m + TM["N0"]
    ii = nu / 2 * sl * cl
    iii = nu / 24 * sl * cl**3 * (5 - t**2 + 9 * eta2)
    iiia = nu / 720 * sl * cl**5 * (61 - 58 * t**2 + t**4)
    iv = nu * cl
    v = nu / 6 * cl**3 * (nu / rho - t**2)
    vi = nu / 120 * cl**5 * (5 - 18 * t**2 + t**4 + 14 * eta2 - 58 * t**2 * eta2)
    north = i + ii * dl**2 + iii * dl**4 + iiia * dl**6
    east = TM["E0"] + iv * dl + v * dl**3 + vi * dl**5
    return east, north


def tm_inverse(east, north):
    lat = math.radians(TM["lat0"])
    a_k0 = KID62["a"] * TM["k0"]
    for _ in range(30):
        m = _meridional_arc(lat)
        dn = north - TM["N0"] - m
        if abs(dn) < 1e-8:
            break
        lat += dn / a_k0
    nu, rho, eta2 = _tm_consts(lat)
    t, sec = math.tan(lat), 1.0 / math.cos(lat)
    vii = t / (2 * rho * nu)
    viii = t / (24 * rho * nu**3) * (5 + 3 * t**2 + eta2 - 9 * t**2 * eta2)
    ix = t / (720 * rho * nu**5) * (61 + 90 * t**2 + 45 * t**4)
    x10 = sec / nu
    xi = sec / (6 * nu**3) * (nu / rho + 2 * t**2)
    xii = sec / (120 * nu**5) * (5 + 28 * t**2 + 24 * t**4)
    xiia = sec / (5040 * nu**7) * (61 + 662 * t**2 + 1320 * t**4 + 720 * t**6)
    de = east - TM["E0"]
    lat_f = lat - vii * de**2 + viii * de**4 - ix * de**6
    lon_f = math.radians(TM["lon0"]) + x10 * de - xi * de**3 + xii * de**5 - xiia * de**7
    return math.degrees(lat_f), math.degrees(lon_f)
