# Contributing

## Getting set up

```bash
make install    # creates ./.venv and installs the package with its dev extras
make check      # lint, type check and the test suite — what CI gates on
```

If `make check` is green on a fresh clone and red after your change, the change is
the problem. If it is red on a fresh clone, that is a bug and worth an issue.

## What CI will hold you to

Seven jobs, each of which is a Makefile target you can run yourself:

| job | command | gate |
|---|---|---|
| lint | `make lint` | `ruff check` and `ruff format --check` clean across the tree |
| typecheck | `make typecheck` | `mypy` clean; new code in `src/` needs annotations |
| test | `make test` | the suite passes on 3.11, 3.12 and 3.13, coverage ≥ 85% |
| lockfile | `make lock-check` | the committed hashes install, and pin what `pyproject.toml` declares |
| audit | `make audit` | `pip-audit` reports nothing against either requirement file |
| dataset | `make dataset` | regenerating changes no committed byte |
| docker | `make docker-verify` | the container migrates and the acceptance suite passes |

## Shape of a change

**One change per commit, with the tests that pin it in the same commit.** A feature
that arrives without its test is not finished, and a commit that mixes a refactor,
a formatting pass and a behaviour change is hard to review and harder to revert.
Formatting-only changes get their own commit.

Write the commit subject as what the change does, in the imperative, under about
70 characters. The body is for *why* — what was wrong before, what you considered,
what you verified. If you ran something to check the change, say what it printed.

## Touching the geodesy

Two rules, both documented at length in [docs/design.md](docs/design.md):

1. **Ellipsoids are arguments, never module constants.** Legacy coordinates have to
   go through the legacy ellipsoid; substituting WGS84 biases every recovered
   parameter by over 100 m while leaving residuals that look fine.
2. **The blunder screen is two stages, and both are load-bearing.** The pre-screen
   cannot see the four quietest miskeyed stations, and the post-fit cut has to land
   between 0.43 m and 1.61 m. If you change either, run
   `tests/test_adjustment.py` and check that
   `test_a_habitual_two_metre_cut_returns_nine_stations_instead_of_ten` still means
   what its name says.

## Changing the dataset

`data/` and `acceptance/ground_truth.json` are generated together by
`tools/generate_dataset.py` from a fixed seed, and CI fails if regenerating changes
a byte. So: never hand-edit either. Change the generator, run `make dataset`, and
commit the generator and everything it regenerated in the same commit.

The true transformation parameters live in `tools/geodesy.py`, and the package
never imports it. Keep it that way — the moment the library can reach the answer,
the acceptance suite stops meaning anything.

## Dependencies

Add them to `pyproject.toml`, then run `make lock` and commit all three lockfiles
with the change. Never edit `requirements*.txt` or `uv.lock` by hand; the lockfile
job recompiles them and diffs.

`make lock-check` (part of `make check`, and its own CI job) verifies that every
dependency `pyproject.toml` declares is pinned in both lockfiles at a version its
specifier admits, and that the dev lock covers the runtime one. It deliberately
does *not* recompile and diff: `uv pip compile` resolves differently depending on
the host operating system, so that check failed on whoever ran `make lock` rather
than on a stale lockfile. `UV_VERSION` in the Makefile records the uv the committed
lockfiles were produced with, so a regeneration can be reproduced if it ever
matters.

## Releases

Update [CHANGELOG.md](CHANGELOG.md) under `Unreleased` as you go. A release moves
that section under a version heading, bumps `version` in `pyproject.toml` and
`__version__` in `src/geodetic_archive/__init__.py`, and gets a tag once CI is
green on `main`.
