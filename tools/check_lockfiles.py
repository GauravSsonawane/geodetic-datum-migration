"""Check that the committed lockfiles still describe what pyproject.toml declares.

The obvious check -- recompile and diff -- turns out not to be reproducible: `uv pip
compile` resolves differently depending on the host operating system, so the same
uv version on Windows and on Linux disagree about numpy even with `--universal`.
That made the check fail on whoever ran `make lock` rather than on a stale
lockfile, which is the opposite of useful.

This checks the property that actually matters instead, and it is deterministic
everywhere: every dependency pyproject.toml declares is pinned in the lockfiles at
a version its specifier admits, and the dev lockfile covers the runtime one. That
catches the mistake the check exists for -- adding or changing a dependency and
forgetting to re-lock -- without depending on resolver behaviour.

    python tools/check_lockfiles.py
"""

from __future__ import annotations

import re
import sys
import tomllib
from pathlib import Path

from packaging.requirements import Requirement
from packaging.utils import canonicalize_name
from packaging.version import Version

REPO_ROOT = Path(__file__).resolve().parent.parent
PYPROJECT = REPO_ROOT / "pyproject.toml"
RUNTIME_LOCK = REPO_ROOT / "requirements.txt"
DEV_LOCK = REPO_ROOT / "requirements-dev.txt"

#: A pinned line in a compiled lockfile: `name==version` with an optional marker and
#: an optional trailing backslash before the hashes.
_PIN = re.compile(r"^(?P<name>[A-Za-z0-9._-]+)==(?P<version>[^\s;\\]+)")


def read_lock(path: Path) -> dict[str, Version]:
    """The pinned versions in a compiled requirements file, by canonical name."""
    pins: dict[str, Version] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith((" ", "\t", "#")) or not line.strip():
            continue
        match = _PIN.match(line.strip())
        if match:
            pins[canonicalize_name(match.group("name"))] = Version(match.group("version"))
    if not pins:
        raise SystemExit(f"{path.name}: no pinned requirements found; is it a compiled lockfile?")
    return pins


def declared_dependencies() -> tuple[list[Requirement], list[Requirement]]:
    """``(runtime, dev)`` requirements as pyproject.toml declares them."""
    with PYPROJECT.open("rb") as handle:
        config = tomllib.load(handle)
    project = config["project"]
    runtime = [Requirement(entry) for entry in project.get("dependencies", [])]
    dev = [Requirement(entry) for entry in project.get("optional-dependencies", {}).get("dev", [])]
    return runtime, dev


def check(requirements: list[Requirement], pins: dict[str, Version], lock_name: str) -> list[str]:
    """Problems with ``requirements`` against ``pins``, as messages."""
    problems: list[str] = []
    for requirement in requirements:
        name = canonicalize_name(requirement.name)
        pinned = pins.get(name)
        if pinned is None:
            problems.append(f"{lock_name}: {requirement.name} is declared but not pinned")
        elif not requirement.specifier.contains(pinned, prereleases=True):
            problems.append(
                f"{lock_name}: {requirement.name} is pinned at {pinned}, "
                f"which {requirement.specifier} does not admit"
            )
    return problems


def main() -> int:
    runtime_requirements, dev_requirements = declared_dependencies()
    runtime_pins = read_lock(RUNTIME_LOCK)
    dev_pins = read_lock(DEV_LOCK)

    problems = check(runtime_requirements, runtime_pins, RUNTIME_LOCK.name)
    problems += check(runtime_requirements + dev_requirements, dev_pins, DEV_LOCK.name)

    # The dev lockfile is compiled with the dev extra on top of the runtime one, so
    # anything the runtime needs has to be in it, at the same version.
    for name, version in runtime_pins.items():
        if dev_pins.get(name) != version:
            problems.append(
                f"{DEV_LOCK.name}: {name} is {dev_pins.get(name, 'absent')} "
                f"but {RUNTIME_LOCK.name} pins {version}"
            )

    if problems:
        print("The lockfiles no longer match pyproject.toml:", file=sys.stderr)
        for problem in problems:
            print(f"  - {problem}", file=sys.stderr)
        print("\nRun `make lock` and commit the result.", file=sys.stderr)
        return 1

    print(
        f"lockfiles current: {len(runtime_requirements)} runtime and "
        f"{len(dev_requirements)} dev dependencies pinned within their specifiers "
        f"({len(runtime_pins)} and {len(dev_pins)} packages resolved)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
