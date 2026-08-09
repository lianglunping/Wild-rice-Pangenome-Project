#!/usr/bin/env bash
set -euo pipefail

ENGINE_ENV="${PANFAMFLOW_ENGINE_ENV:-panfamflow-engine}"

for executable in mamba uv; do
  if ! command -v "$executable" >/dev/null 2>&1; then
    echo "$executable is required. Install Miniforge (mamba) and uv first." >&2
    exit 2
  fi
done

if mamba env list | awk '{print $1}' | grep -Fxq "$ENGINE_ENV"; then
  mamba env update -n "$ENGINE_ENV" -f environment.yaml --prune
else
  mamba env create -n "$ENGINE_ENV" -f environment.yaml
fi

uv sync --all-extras --dev

cat <<MSG
Bootstrap completed.
Python CLI:
  uv run panfamflow --help
Workflow engine:
  mamba run -n "$ENGINE_ENV" snakemake --version
Validation:
  uv run panfamflow validate -c examples/toy/config.yaml
Resume after a corrected failure:
  uv run panfamflow resume -c config.yaml
MSG
