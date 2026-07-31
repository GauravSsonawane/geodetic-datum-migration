The Kaveri Island Survey Department is retiring its 1962 local datum and needs its
geodetic station register republished on WGS84. The archive is in `/app/data`: the
station register `/app/data/archive.dat`, the datum connection stations
`/app/data/tiepoints.csv`, the separation model `/app/data/geoid_kid62.grd`, and
`/app/data/FORMAT.txt`, which is the normative specification of all three.

The transformation from KID62 to WGS84 was never published and the department's
working copy is lost, so recover it from the connection stations. Their WGS84
coordinates come from a modern GNSS campaign and are good to about a decimetre.
Their legacy coordinates were transcribed by hand from the 1962-1971 field books
and a number of them are grossly wrong; those stations must not influence the
result.

Write three files. The schemas below are normative. Use the keys exactly as written.

`/app/out/datum_params.json`, the recovered KID62 to WGS84 transformation as a
seven-parameter geocentric similarity transformation, in the position-vector
(Bursa-Wolf, EPSG 9606) convention:

```json
{
  "dx_m": 0.0, "dy_m": 0.0, "dz_m": 0.0,
  "rx_arcsec": 0.0, "ry_arcsec": 0.0, "rz_arcsec": 0.0,
  "scale_ppm": 0.0
}
```

`/app/out/archive_wgs84.geojson`, a GeoJSON `FeatureCollection` with one `Point`
feature per migrated record. Every archive record that carries both a complete
position and a complete height must appear exactly once; records with any void
field are excluded. Feature coordinates are `[longitude, latitude,
ellipsoidal_height]` in degrees and metres on WGS84, and each feature carries:

```json
{
  "record_id": "KA00001",
  "name": "full station name",
  "height_orthometric_m": 0.0
}
```

`height_orthometric_m` is the record's height as archived, in metres.

`/app/out/qa_report.json`:

```json
{
  "n_records_migrated": 0,
  "n_records_excluded_void": 0,
  "excluded_record_ids": ["KA00001"],
  "rejected_tie_points": ["TP001"],
  "tie_point_residuals": {"TP001": [0.0, 0.0, 0.0]},
  "rms_residual_m": 0.0
}
```

`rejected_tie_points` lists the stations whose legacy coordinates are
mis-transcribed. `tie_point_residuals` covers the stations your final adjustment
was computed from, each as an east, north, up triple in metres in a local
topocentric frame at that station, taken as the GNSS position minus the
transformed legacy position. `rms_residual_m` is the three-dimensional RMS
residual over those same stations.
