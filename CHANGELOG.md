# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project uses
[semantic versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

Nothing yet.

## [0.1.0] — 2026-09-01

First release as a standalone toolkit. The geodesy and the dataset are unchanged;
everything around them is new.

### Added

- `geodetic_archive`, an installable package with a `geodetic-archive` CLI,
  extracted from what had been a single 329-line script. Seven modules, each with a
  typed public surface: `ellipsoid`, `projection`, `transform`, `geoid`, `archive`,
  `adjustment` and `pipeline`.
- Named exceptions carrying the file, record and line at fault —
  `ArchiveFormatError`, `GeoidGridError`, `AdjustmentError`, `ConventionError` —
  in place of raw `ValueError`, `IndexError` and `KeyError` from library internals.
- A host-side test suite: 152 tests at 99% line and branch coverage, gated at 85%.
  It pins the properties the migration turns on, including that the pre-screen
  alone cannot see the four quietest miskeyed stations, and that a habitual 2 m
  post-fit cut returns nine rejected stations instead of ten.
- `logging` through the adjustment and the pipeline, with `--log-level` on the CLI.
- `pyproject.toml` configuring ruff, mypy, pytest and coverage, with three
  committed lockfiles: `uv.lock`, and hash-pinned `requirements.txt` and
  `requirements-dev.txt`.
- A `Makefile` covering the whole fresh-clone path, and `docker-compose.yml` with
  `migrate` and `verify` services.
- `.github/workflows/ci.yml`: lint, type check, tests across Python 3.11–3.13, a
  lockfile currency check, a dependency audit, a byte-for-byte dataset
  regeneration check, and a containerised end-to-end run. Plus Dependabot.
- `docs/brief.md` and `docs/design.md`, `CONTRIBUTING.md`, this changelog, and an
  MIT `LICENSE`.

### Changed

- Repository laid out as a Python project: `src/`, `tests/`, `acceptance/`,
  `data/`, `docs/`, `docker/`, `tools/`.
- The acceptance suite takes its output directory from `GEODETIC_OUT_DIR` and
  defaults to `./out`, so it runs on the host as readily as in a container.
- The container is built from a digest-pinned `python:3.13-slim` base and installs
  from the lockfile with `--require-hashes`, replacing an Ubuntu image with
  dependencies pinned inline and `--break-system-packages`.
- Development dependency on pytest moved to `>=9.0.3` to clear PYSEC-2026-1845.

### Removed

- Four workflows that pointed into an external review pipeline and ran none of this
  project's own checks, along with that pipeline's rubric, base-image policy and
  label taxonomy under `references/`.
- The submission-bundle scaffolding: `task.toml`, `solve.sh` and `test.sh`.
- `tools/reference_solve.py`, superseded by the package extracted from it. Keeping
  both would have left two implementations of the same geodesy to drift apart.

[Unreleased]: https://github.com/GauravSsonawane/geodetic-datum-migration/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/GauravSsonawane/geodetic-datum-migration/releases/tag/v0.1.0
