# Remote validation contract

This file defines the minimum evidence required before PanFamFlow is described as remotely published or validated.

## Scope

PanFamFlow is a **target pan-gene-family analysis workflow**. It consumes already assembled and annotated genomes or materials. It does not assemble genomes, build graph pangenomes, call genome-wide SV/PAV, or classify every genome-wide HOG as a general pangenome product.

## Required repository checks

- `PanFamFlow/` is directly browsable on the feature branch.
- The canonical module and output names are `pan_family` and `06_pan_family`.
- The canonical configuration rejects `whole_genome` scope.
- `_panfamflow_publish/` and `_panfamflow_staging/` are absent.
- A durable root-level `panfamflow-ci.yml` workflow exists.

## Required CI checks

- locked `uv` environment sync;
- Ruff lint and format checks;
- mypy;
- pytest;
- wheel and source-distribution build;
- CLI scope, validation and planning checks;
- Snakemake toy DAG dry-run;
- an executed smart-resume smoke test proving that a second identical run skips valid outputs and that deletion of one terminal marker only regenerates that marker and its required downstream state.

## Biological validation boundary

Software CI and toy recovery tests do not establish biological validity. A real-data benchmark on a manually reviewed target rice gene family remains required before production or publication-grade biological claims are made.
