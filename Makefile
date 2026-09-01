# Every target here is what CI runs, so a green `make check` locally means a green
# pipeline. Run `make help` for the list.
#
# Tools are taken from ./.venv when it exists and from PATH otherwise, so the same
# targets work for a developer with a virtualenv and for CI, which installs into
# the runner's own Python.

PYTHON ?= python3
VENV ?= .venv

ifeq ($(OS),Windows_NT)
VENV_BIN := $(VENV)/Scripts
else
VENV_BIN := $(VENV)/bin
endif

ifneq ($(wildcard $(VENV_BIN)/.),)
PY := $(VENV_BIN)/python
TOOL := $(VENV_BIN)/
else
PY := $(PYTHON)
TOOL :=
endif

DATA_DIR ?= data
OUT_DIR ?= out
IMAGE ?= geodetic-archive

.DEFAULT_GOAL := help
.PHONY: help install lock lint format typecheck test audit check migrate verify \
        dataset docker-build docker-verify clean

help: ## Show this help
	@grep -hE '^[a-zA-Z_-]+:.*?## ' $(MAKEFILE_LIST) \
		| awk 'BEGIN {FS = ":.*?## "}; {printf "  %-14s %s\n", $$1, $$2}'

install: ## Create ./.venv and install the package with its dev extras
	$(PYTHON) -m venv $(VENV)
	$(VENV_BIN)/python -m pip install --upgrade pip
	$(VENV_BIN)/python -m pip install -e ".[dev]"

lock: ## Recompile the hash-pinned lockfiles from pyproject.toml (needs uv)
	# --no-header: the header records the -o path, so without it the CI currency
	# check would diff on where it recompiled rather than on what it resolved.
	uv pip compile pyproject.toml --python-version 3.11 --generate-hashes \
		--no-header -o requirements.txt
	uv pip compile pyproject.toml --extra dev --python-version 3.11 --generate-hashes \
		--no-header -o requirements-dev.txt
	uv lock

lint: ## Check lint rules and formatting
	$(TOOL)ruff check .
	$(TOOL)ruff format --check .

format: ## Apply the formatter and the safe lint fixes
	$(TOOL)ruff check . --fix
	$(TOOL)ruff format .

typecheck: ## Run mypy over the package, its tests and the authoring tools
	$(TOOL)mypy

test: ## Run the host-side suite with coverage, gated at 85%
	$(PY) -m pytest

audit: ## Check the pinned dependencies for known vulnerabilities
	$(TOOL)pip-audit --requirement requirements.txt --requirement requirements-dev.txt

check: lint typecheck test ## Everything CI gates on, in one go

migrate: ## Run the migration into $(OUT_DIR)
	$(PY) -m geodetic_archive.cli migrate --data $(DATA_DIR) --out $(OUT_DIR)

verify: ## Grade the files in $(OUT_DIR) against ground truth
	GEODETIC_OUT_DIR=$(OUT_DIR) $(PY) -m pytest acceptance -rA --no-cov \
		-p no:cacheprovider

dataset: ## Regenerate the inputs and ground truth from the fixed seed
	$(PY) tools/generate_dataset.py .

docker-build: ## Build the runtime image
	docker build -f docker/Dockerfile --target runtime -t $(IMAGE):local .

docker-verify: ## Migrate and grade inside containers
	docker compose run --rm migrate
	docker compose run --rm verify

clean: ## Remove build artefacts, caches and the run output
	rm -rf $(OUT_DIR) build dist htmlcov .coverage coverage.xml
	rm -rf .pytest_cache .ruff_cache .mypy_cache src/*.egg-info
	find . -name __pycache__ -type d -prune -exec rm -rf {} +
