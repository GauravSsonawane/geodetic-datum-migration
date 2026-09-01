# Design notes

Why the migration is built the way it is, and why the numbers that grade it sit
where they do. The problem statement is in [brief.md](brief.md).

## The shape of the job

Every country that has moved off a local datum has paid someone to do this job. The
awkward version, and a common one in older archives, is that the transformation to
WGS84 was never written down, so it has to be pulled back out of the register
before a single record can be moved. The three outputs are what a geodesist hands
over at the end of it: the parameters, the migrated register, and a QA sheet
showing which control was used and how well it fitted.

The difficulty is geodetic rather than programmatic. There is no clever algorithm
here and nothing is computationally heavy; the whole run takes under two seconds.
What it takes is getting five separate things right, each of which produces tidy,
plausible, wrong output when got wrong.

## Five ways to produce output that is wrong

### 1. The legacy ellipsoid

A Helmert transformation lives in geocentric Cartesian space, so the legacy
positions have to go through the *legacy* ellipsoid on the way in. KID62's
semi-major axis is 841 m off WGS84's. Reach for a familiar default instead and
every recovered parameter is out by more than a hundred metres while the residuals
still look reasonable — the sort of error that survives review.

Guarded by `Ellipsoid` being a value the geodesy functions take rather than a
module constant they close over, and by
`test_the_two_ellipsoids_disagree_by_hundreds_of_metres`.

### 2. Three coordinate encodings

Legacy positions arrive as packed sexagesimal, decimal degrees, and plane
coordinates on a local Transverse Mercator grid. The plane subset has to come back
through an inverse projection on that same legacy ellipsoid before it is any use.

### 3. Errors in the connection stations, on a ladder

Ten legacy coordinates were miskeyed, running 16 000 km, 239 km, 73 km, 998 m,
2.7 km, then 27, 8, 3.1, 2.5 and 1.5 m.

Least squares has no defence against the large ones: one station kilometres out of
place spreads itself across every residual, so fitting once and screening the
residuals fingers the *wrong* stations. But a robust pre-screen is not enough
either. Its honest scatter reaches 2.5 m, so the four quiet errors pass straight
through into the adjustment, and are only separable afterwards, where honest
residuals top out at 0.43 m and the smallest blunder sits at 1.61 m.

**The cut for that second pass is the crux.** It has to be read off the residual
spread and land in the window (0.43, 1.61) m. A cut scaled to the robust spread, or
floored near the decimetre quoted for the GNSS side, lands inside it. A round "2 m"
— a natural number to reach for — does not: it leaves the 1.5 m station in the fit
and returns nine stations instead of ten.

Fifteen cut rules a practitioner might plausibly use were run against the data.
Eleven return the correct set — every MAD-based rule from 3σ to 6σ, 3× RMS,
5× median, mean+3σ, and fixed cuts at 0.5 and 1.0 m. The four that fail are all
methodologically unsound: 4× RMS and mean+4σ use a **non-robust scale inflated by
the very blunders they are meant to find**, a fixed 0.3 m sits below the honest
floor, and a fixed 2 m is the habitual constant this data is built to punish.

Separately, *no* rule works without the raw-difference pre-screen — robust rules
then reject nothing, because the 16 000 km error inflates the MAD past every cut,
and fixed cuts reject 24 of 34. The two-stage structure is load-bearing, which is
why `test_the_screen_cannot_see_the_quiet_blunders` asserts it directly rather than
leaving it as a comment.

### 4. The vertical datum

Archived heights are orthometric over local MSL; the transformation and the
requested output are both ellipsoidal. Skip the separation grid and every height is
out by about 38 m *with the horizontal work still checking out*. This is the
quietest failure in the archive.

### 5. The rotation sign convention

The brief names position vector (EPSG 9606). The coordinate frame convention (EPSG
9607) differs by the sign of all three rotations, and the difference is invisible in
the horizontal fit. `HelmertParameters.in_convention()` makes the choice explicit;
`test_reporting_the_wrong_convention_moves_the_station_metres` pins what it costs.

Record assembly across continuation lines and void marker handling sit on top of
all that, and both feed the record counts.

## Where the tolerances come from

The acceptance suite in [`acceptance/`](../acceptance) grades the three output
files and nothing else. Every band below was checked in both directions: that it
admits other defensible adjustments, and that it rejects the methodological
failures above. Each variant was written out as a full set of output files and put
through the real suite rather than argued about on paper.

| approach | horiz RMS | vert RMS | worst parameter error | verdict |
|---|---|---|---|---|
| this implementation | 0.018 m | 0.033 m | 2.45 m, 0.089″, 0.08 ppm | **pass** |
| nearest neighbour geoid interpolation | 0.018 m | 0.378 m | same | **pass** |
| only 12 connection stations used | 0.052 m | 0.075 m | 4.20 m, 0.155″, 0.226 ppm | **pass** |
| translation only (3 parameter) shift | 0.894 m | 0.994 m | 73 m, 2.42″ | fail |
| WGS84 ellipsoid on the legacy side | 3.011 m | 0.084 m | 139 m, 110 ppm | fail |
| separation grid ignored | 0.017 m | 38.430 m | same as reference | fail |
| coordinate frame rotation convention | 0.018 m | 0.033 m | 4.84″ | fail |
| blunders left in the adjustment | 138 km | 698 m | 6026 km | fail |

**Parameters** are graded at 10.0 m, 0.35″ and 1.0 ppm — roughly four sigma. Formal
sigmas on clean control are 2.2 to 2.7 m on the translations, 0.06 to 0.08″ on the
rotations and 0.27 ppm on the scale, because the seven parameters of a small-region
similarity transformation are correlated. The bands have to bracket *parameter*
uncertainty, not the much tighter coordinate uncertainty; grading them at half a
metre would fail correct work.

**Positions** are graded at 0.30 m RMS and 0.75 m worst case. This implementation
gives 0.018 m, the loosest defensible variant 0.052 m, and the nearest failure
0.894 m RMS with 1.54 m worst case — so the band sits about six times above the
sound spread and three times below the nearest failure.

**Heights** are graded at 1.50 m RMS and 4.00 m worst case. Nearest-neighbour
interpolation of the separation grid gives 0.378 m RMS and 1.26 m worst case, which
is a defensible reading even though the model documents bilinear, so the band lets
any reasonable interpolation through and still rejects skipping the step by a
factor of twenty five.

**The rejected station set** is graded on exact equality, which is safe because the
populations do not overlap: against a clean fit, honest residuals run 0.03–0.43 m
and all ten blunders sit at 1.61 m or above. Any cut inside that window returns the
same set, whatever route the solver took to it. Note this failure is caught by
`test_blunder_stations_rejected` *alone* — such a fit still has in-tolerance
parameters and an in-band RMS, which is why the set is compared exactly rather than
by count.

**Objective vs. judgement.** `rejected_tie_points` asks whether a coordinate was
miskeyed, which has one right answer, so it is graded exactly. Which *sound*
stations to carry into the adjustment is the solver's call — dropping a few extra
as a precaution is ordinary practice — so the residual test accepts any subset of
the sound stations, provided no miskeyed one was used and at least 10 remain. A
conservative adjustment using 20 of the 24 sound stations passes 15/15.

## What the acceptance suite refuses to read

The suite reads no input file. The connection stations it certifies residuals
against are recorded in `ground_truth.json` when the archive is generated, not
re-read from `data/tiepoints.csv` — that file sits on a filesystem the
implementation can write to, so grading residuals against it would let a rewritten
copy certify fabricated ones. Verified: destroying every legacy coordinate in that
file leaves the verdict unchanged.

Reported residuals are recomputed from the implementation's *own* parameters and
compared as whole east/north/up vectors, allowing an overall sign flip so either
direction convention passes. Magnitude-only comparison was tried first and is not
enough: any triple of the right length satisfies it, including `[‖residual‖, 0, 0]`,
which fabricates a zero north and zero vertical component.

Every numeric field must be finite. NaN compares `False` against any bound, so a
check phrased as "flag it when the error exceeds the tolerance" accepts NaN
silently — the tests assert the passing condition instead.

**Known and deliberate:** the JSON schema checks require the documented keys but do
not reject *extra* ones. Enforcing exact key sets would fail an otherwise-correct
submission that adds a harmless diagnostic field, turning a cosmetic looseness into
a formatting-based failure mode. Trivial to tighten if you would rather it were.

## Data provenance

Synthetic, from an invented datum, so no parameter set or answer is retrievable
online. Extents, ellipsoid, grid definition, geoid gradient (−42.4 m to −33.6 m
across the island), decimetre GNSS noise and the list of transcription errors are
modelled on real national archives of the period.

[`tools/generate_dataset.py`](../tools/generate_dataset.py) builds the inputs and
the ground truth together from a fixed seed, with the true parameters in
[`tools/geodesy.py`](../tools/geodesy.py):

```bash
make dataset          # python tools/generate_dataset.py .
```

It sits outside the installable package on purpose, and the package never imports
it: the true parameters live in `tools/` and nowhere else, so an implementation
cannot read the answer out of the library it is using. Keeping it in the repository
means the 1170 records in `acceptance/ground_truth.json` can be traced back to what
produced them instead of being taken on trust, and the inputs and the ground truth
cannot drift apart. CI regenerates on every push and fails if a single byte moves.
`FORMAT.txt` is hand-written and static.
