# Runtime image for the migration. Build from the repository root:
#
#     docker build -f docker/Dockerfile -t geodetic-archive .
#
# The dataset is baked in so the image runs with no bind mount, and the outputs go
# to /app/out, which is where a mount belongs if you want them on the host.
#
# The `verify` target adds the test dependencies and the acceptance suite, so the
# grading runs against the same image that produced the outputs:
#
#     docker build -f docker/Dockerfile --target verify -t geodetic-archive:verify .
FROM python:3.15-slim-bookworm@sha256:2a41267d9de9d358f77a1de356e0f16ecba8d39827855be0ce9ef4792edce402 AS runtime

# Fail fast, and keep pip from writing anything it will not need again.
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

WORKDIR /app

# Dependencies first, from the hash-pinned lockfile, so a source change does not
# invalidate this layer. --require-hashes makes an unpinned transitive dependency
# a build failure rather than a silent substitution.
COPY requirements.txt ./
RUN pip install --require-hashes --no-deps -r requirements.txt

# Then the package itself.
COPY pyproject.toml README.md ./
COPY src ./src
RUN pip install --no-deps .

# The legacy archive to migrate, and the directory the outputs go to.
COPY data ./data
RUN mkdir -p /app/out

ENTRYPOINT ["geodetic-archive"]
CMD ["migrate", "--data", "/app/data", "--out", "/app/out"]


# --- verify -------------------------------------------------------------------
# Grades the three output files against ground truth. It reads nothing else: the
# inputs are on a filesystem a solver can write to, so certifying results against
# them would certify nothing.
FROM runtime AS verify

COPY requirements-dev.txt ./
RUN pip install --require-hashes --no-deps -r requirements-dev.txt

COPY acceptance ./acceptance
ENV GEODETIC_OUT_DIR=/app/out

ENTRYPOINT ["python", "-m", "pytest"]
# --no-cov: pyproject measures coverage of the package, and the acceptance suite
# deliberately exercises none of it -- it reads the written files and nothing else.
CMD ["/app/acceptance", "-rA", "--no-cov", "-p", "no:cacheprovider"]
