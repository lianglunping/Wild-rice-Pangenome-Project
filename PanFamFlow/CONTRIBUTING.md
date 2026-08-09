# Contributing

Create a feature branch and keep each change scoped. Do not commit raw biological data, private accessions, credentials, Conda caches, `.snakemake/`, `work/`, or generated `results/`.

Before opening a pull request:

```bash
uv lock --check
uv sync --locked --dev
uv run ruff check .
uv run ruff format --check .
uv run mypy src/panfamflow
uv run pytest -q
uv build
```

Workflow changes should include a unit test or a toy dry-run target. Any change to a biological threshold must update the config template, documentation, changelog and source-traceability rationale.
