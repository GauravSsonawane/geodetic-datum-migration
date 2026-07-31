# dynamo/migrate-geodetic-archive

Development notes for reviewers. The task itself is entirely under [`task/`](task/).
This file is not passed to the agent.

## Overview

A fictional national mapping agency is retiring the Kaveri Island Datum 1962
(KID62) and needs its geodetic station register republished on WGS84. The
transformation between the two datums was never published, so it has to be
recovered from the archive's own tie stations before anything can be migrated, and
then applied to 1200 legacy records.

The agent gets four files in `/app/data`: the fixed-width station register
`archive.dat`, the tie stations `tiepoints.csv`, the geoid separation grid
`geoid_kid62.grd`, and `FORMAT.txt`, which is the normative spec for all three. It
writes the recovered parameters, a migrated GeoJSON `FeatureCollection`, and a QA
report.

## Where it bites

Five separate ways to produce tidy output that is wrong:

1. **Legacy ellipsoid.** Both the geocentric conversion and the inverse projection
   have to use the KID62 ellipsoid from `FORMAT.txt`, whose semi-major axis is
   841 m off WGS84. A familiar default biases every parameter by over 100 m and
   still leaves plausible-looking residuals.
2. **Three coordinate encodings.** Packed sexagesimal, decimal degrees, and plane
   coordinates on a local Transverse Mercator grid. The last needs an inverse
   projection before it is any use.
3. **Errors in the tie stations, loud and quiet.** Seven legacy coordinates were
   miskeyed: transposed seconds, a flipped hemisphere letter, decimal degrees typed
   into a sexagesimal field, easting and northing the wrong way round, a kilometre
   slip, and two smaller digit slips. Least squares has no defence against the
   large ones — fitting once and screening the residuals fingers the *wrong*
   stations, and can settle on a parameter set that fits the bad control and
   discards the good. But a robust pre-screen alone is not enough either: the
   smallest error is worth 8 m and hides inside the scatter the screen must
   tolerate, surfacing only in the residuals of an already-clean fit. Screening
   once and stopping carries it into the parameters, and is rejected three ways —
   wrong rejected set, RMS 1.42 m outside the [0.02, 1.00] band, and translations
   off by 20.9 m against a 10.0 m tolerance.
4. **Vertical datum.** Archived heights are orthometric over local MSL, while the
   transformation and the requested output are ellipsoidal. Skip the separation
   grid and every height is out by about 38 m with the horizontal work still
   checking out.
5. **Rotation sign convention.** The prompt names position vector (EPSG 9606). The
   coordinate frame convention differs by the sign of all three rotations.

Record assembly across continuation lines and void marker handling sit on top of
that, and both feed the record counts the verifier checks.

## Environment

Single image on the pre-approved `ubuntu:24.04` digest, with `numpy`, `pytest` and
`pytest-json-ctrf` pinned and baked so the verifier installs nothing. `COPY data
/app/data` ships the four inputs. `solution/` and `tests/` are never copied in, so
ground truth is out of reach during the agent's run; this was confirmed by
inspecting the built image.

## Verification

`tests/ground_truth.json`, overlaid at `/tests` only at verify time, holds the
generating parameters, the blunder and void identifiers, and the true position,
ellipsoidal height, archived orthometric height and full name of every complete
record. Fifteen atomic tests grade the three requested output files and nothing
else.

Every row below was written out as a full set of output files and put through the
real verifier, rather than argued about on paper:

| approach | horiz RMS | vert RMS | worst parameter error | reward |
|---|---|---|---|---|
| reference solution | 0.018 m | 0.033 m | 2.45 m, 0.089″, 0.08 ppm | **1** |
| nearest neighbour geoid interpolation | 0.018 m | 0.378 m | same | **1** |
| only 12 tie stations used | 0.052 m | 0.075 m | 4.20 m, 0.155″, 0.226 ppm | **1** |
| translation only (3 parameter) shift | 0.894 m | 0.994 m | 73 m, 2.42″ | 0 |
| WGS84 ellipsoid on the legacy side | 3.011 m | 0.084 m | 139 m, 110 ppm | 0 |
| separation grid ignored | 0.017 m | 38.430 m | same as reference | 0 |
| coordinate frame rotation convention | 0.018 m | 0.033 m | 4.84″ | 0 |
| blunders left in the adjustment | 138 km | 698 m | 6026 km | 0 |

The parameter bands are deliberately loose and the coordinate bands tight. Formal
sigmas on clean control are 2.2 to 2.7 m on the translations and 0.06 to 0.08″ on
the rotations, because the seven parameters of a small region similarity
transformation are correlated. Grading the parameters at half a metre would fail
correct work.

The rejected station set is graded on exact equality, which is safe because the
populations do not overlap. Six blunders move their stations by 27 m to 16 000 km
against a worst honest screening statistic of 2.52 m. The seventh is the 8 m one,
which a single screening pass may or may not catch depending on where its cut
falls — either way it ends up rejected, because after a clean fit it sits ~20×
above the worst honest post-fit residual of 0.41 m. The same set comes out
regardless of which stage removes each station.

The verifier reads no agent-writable path except the three requested outputs. The
tie stations it certifies residuals against are recorded in `ground_truth.json` when
the archive is generated, not re-read from `/app/data/tiepoints.csv` — that file is
on the agent's own filesystem, so grading residuals against it would let a rewritten
copy certify fabricated ones. Verified: destroying every legacy coordinate in that
file leaves the verdict unchanged.

Reported residuals are recomputed from the agent's *own* parameters and compared as
whole east/north/up vectors, allowing an overall sign flip so either direction
convention passes. Magnitude-only comparison was tried first and is not enough: any
triple of the right length satisfies it, including `[‖residual‖, 0, 0]`, which
fabricates a zero north and zero vertical component.

Every numeric field must be finite. NaN compares `False` against any bound, so a
check phrased as "flag it when the error exceeds the tolerance" accepts NaN
silently — the tests assert the passing condition instead.

Under Harbor 0.20.0:

```
harbor run -p task --agent oracle   ->  1 trial, 0 exceptions, mean 1.000
harbor run -p task --agent nop      ->  1 trial, 0 exceptions, mean 0.000
```

## Data provenance

Synthetic, from an invented datum, so no parameter set or answer is retrievable
online. Extents, ellipsoid, grid definition, geoid gradient (−42.4 m to −33.6 m
across the island), decimetre GNSS noise and the list of transcription errors are
modelled on real national archives of the period.

[`tools/generate_dataset.py`](tools/generate_dataset.py) builds the inputs and the
ground truth together from a fixed seed, with the true parameters in
[`tools/geodesy.py`](tools/geodesy.py):

```bash
python tools/generate_dataset.py .
```

It sits outside `task/` on purpose. It is never built into the image and the task
does not reference it, but keeping it in the repo means the 1170 records in
`tests/ground_truth.json` can be traced back to what produced them instead of being
taken on trust, and the inputs and the ground truth cannot drift apart. Verified:
regenerating reproduces `archive.dat`, `tiepoints.csv`, `geoid_kid62.grd` and
`ground_truth.json` byte for byte. `FORMAT.txt` is hand-written and static.
