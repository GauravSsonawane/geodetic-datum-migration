# geodetic-datum-migration

[![CI](https://github.com/GauravSsonawane/geodetic-datum-migration/actions/workflows/ci.yml/badge.svg)](https://github.com/GauravSsonawane/geodetic-datum-migration/actions/workflows/ci.yml)
[![Python 3.11+](https://img.shields.io/badge/python-3.11%20%7C%203.12%20%7C%203.13-blue)](pyproject.toml)
[![License: MIT](https://img.shields.io/badge/license-MIT-green)](LICENSE)

Recover an unpublished geodetic datum transformation from blunder-contaminated
control, and migrate a legacy survey archive onto WGS84.

A fictional national mapping agency is retiring the Kaveri Island Datum 1962
(KID62) and needs its geodetic station register republished on WGS84. The
transformation between the two datums was never published, so it has to be
recovered from the archive's own connection stations before anything can be
migrated — and ten of those stations were miskeyed by hand in the 1960s. The
package fits the seven parameters, screens the miskeyed control out in two stages,
and applies the result to 1200 legacy records.

**Project type:** a Python library and CLI, plus a synthetic dataset and an
acceptance suite that grades a run against ground truth. It is not a deployed
service: there is no server, no database and no network I/O, so there is nothing
here to health-check or scrape metrics from. Everything runs offline in about two
seconds.

---

## Quick start

```bash
git clone https://github.com/GauravSsonawane/geodetic-datum-migration.git
cd geodetic-datum-migration

make install          # creates ./.venv and installs the package with dev extras
make check            # ruff, mypy, and the test suite with coverage gated at 85%
make migrate          # writes the three outputs into ./out
make verify           # grades ./out against ground truth
```

`make help` lists every target. If you would rather not use make:

```bash
python -m venv .venv && .venv/bin/python -m pip install -e ".[dev]"
.venv/bin/geodetic-archive migrate --data data --out out
GEODETIC_OUT_DIR=out .venv/bin/python -m pytest acceptance -rA --no-cov
```

For a reproducible install, use the hash-pinned lockfile instead of the extras:

```bash
python -m pip install --require-hashes --no-deps -r requirements-dev.txt
python -m pip install --no-deps -e .
```

### With Docker

No local Python needed:

```bash
docker compose run --rm migrate    # writes ./out on the host
docker compose run --rm verify     # grades it, 15 checks
```

The image is built from a digest-pinned `python:3.13-slim` base and installs from
the same hash-pinned lockfile, so it and a host install resolve to the same
dependency tree.

## The command line

```
geodetic-archive migrate [--data DIR] [--out DIR] [--log-level LEVEL]
```

Reads `archive.dat`, `tiepoints.csv` and `geoid_kid62.grd` from `--data`, and
writes three files into `--out`:

| file | contents |
|---|---|
| `datum_params.json` | the seven recovered parameters, position-vector convention (EPSG 9606) |
| `archive_wgs84.geojson` | a `FeatureCollection` of `Point` features, one per migrated record |
| `qa_report.json` | rejected stations, per-station residuals, RMS, and the record counts |

A run over the shipped archive prints:

```
adjusted   : 24 stations
rejected   : TP004, TP007, TP009, TP011, TP014, TP016, TP022, TP024, TP027, TP030
rms        : 0.1894 m
migrated   : 1170 records, 30 excluded for void fields
```

## How it works

```
data/                     read
  tiepoints.csv  ──▶ archive.read_tie_points ──▶ ellipsoid.geodetic_to_ecef (KID62 | WGS84)
                                                          │
                                                          ▼
                                            adjustment.screen_gross_errors     stage 1: raw
                                                          │                    difference vectors
                                                          ▼
                                            adjustment.adjust ⟳ fit + snoop    stage 2: post-fit
                                                          │                    residuals
                                                          ▼
                                              HelmertParameters (7)
                                                          │
  archive.dat    ──▶ archive.read_archive ────────────────┤
  geoid_kid62.grd ─▶ geoid.GeoidGrid ─── H → h ───────────┤
                                                          ▼
                                              pipeline.migrate ──▶ out/
```

| module | responsibility |
|---|---|
| [`ellipsoid.py`](src/geodetic_archive/ellipsoid.py) | reference ellipsoids; geodetic ⇄ geocentric |
| [`projection.py`](src/geodetic_archive/projection.py) | the Kaveri Grid 1962 Transverse Mercator, forward and inverse |
| [`transform.py`](src/geodetic_archive/transform.py) | the seven-parameter transformation, its conventions, the topocentric frame |
| [`geoid.py`](src/geodetic_archive/geoid.py) | the separation model; orthometric → ellipsoidal heights |
| [`archive.py`](src/geodetic_archive/archive.py) | the fixed-width register, its three encodings and its void markers |
| [`adjustment.py`](src/geodetic_archive/adjustment.py) | the two-stage blunder screen and the least-squares fit |
| [`pipeline.py`](src/geodetic_archive/pipeline.py) | the end-to-end run and the three output files |
| [`cli.py`](src/geodetic_archive/cli.py) | argument parsing, logging, exit codes |

The interesting part is `adjustment.py`. Least squares has no defence against a
gross error, so fitting once and screening the residuals fingers the *wrong*
stations. But a single robust pre-screen is not enough either: its honest scatter
reaches 2.5 m, so the four quietest blunders survive into the adjustment and are
only separable afterwards. The post-fit cut has to land between 0.43 m and 1.61 m,
and a habitual round "2 m" does not.
[docs/design.md](docs/design.md) works through all of that, and through where every
grading tolerance comes from.

## Layout

```
src/geodetic_archive/   the library and its CLI
tests/                  host-side suite: unit tests plus an end-to-end run
acceptance/             grades a run's three output files against ground truth
data/                   the archive, the connection stations, the separation model,
                        and FORMAT.txt, which is normative for all three
tools/                  dataset generator; builds data/ and ground truth together
docs/                   the problem statement and the design notes
docker/                 runtime and verification images
```

Two suites, deliberately separate. `tests/` exercises the package: 152 tests, 99%
line and branch coverage, gated at 85%. `acceptance/` grades only the three written
files — it imports nothing from the package and reads no input file, so a run that
rewrote its own inputs could not certify itself. See
[docs/design.md](docs/design.md#what-the-acceptance-suite-refuses-to-read).

## Development

| command | what it does |
|---|---|
| `make check` | lint, type check and the test suite — what CI gates on |
| `make lint` / `make format` | `ruff check` and `ruff format` |
| `make typecheck` | `mypy` over the package, its tests and `tools/` |
| `make test` | pytest with coverage, gated at 85% |
| `make audit` | `pip-audit` over both requirement files |
| `make lock` | recompile the hash-pinned lockfiles from `pyproject.toml` |
| `make dataset` | regenerate `data/` and ground truth from the fixed seed |
| `make docker-verify` | migrate and grade inside containers |

CI runs all of it on every push and pull request, across Python 3.11, 3.12 and
3.13, and additionally checks that the lockfiles still match `pyproject.toml` and
that regenerating the dataset changes no committed byte.

Contributions: see [CONTRIBUTING.md](CONTRIBUTING.md). Changes:
[CHANGELOG.md](CHANGELOG.md).

## Data

Synthetic, from an invented island datum, so no parameter set or answer is
retrievable online. The extents, ellipsoid, grid definition, geoid gradient,
decimetre GNSS noise and the list of transcription errors are modelled on real
national archives of the period. `tools/generate_dataset.py` builds the inputs and
the ground truth together from a fixed seed; CI fails if regenerating changes a
byte.

## Licence

MIT — see [LICENSE](LICENSE).
