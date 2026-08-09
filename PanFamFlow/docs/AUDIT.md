# Independent-style audit and corrective actions

This document records issues that a third-party reviewer would identify, including problems not visible from the analysis template screenshots.

## Release and repository integrity

### Finding A1 — source branch could display a 404

A release claim is invalid if the branch does not contain the documented `PanFamFlow/` path. Temporary Base64 chunks, one-time bootstrap workflows or external archives are not a substitute for committed source files.

Corrective action:

- commit the complete source tree directly under `PanFamFlow/`;
- remove `_panfamflow_staging/` and one-time bootstrap/finalize/recover workflows;
- keep one root-level long-lived CI workflow;
- verify the exact branch URL through GitHub's contents API before reporting publication.

Acceptance:

```text
GET /contents/PanFamFlow?ref=<branch> -> directory listing
GET /contents/PanFamFlow/pyproject.toml?ref=<branch> -> 200
GET /contents/PanFamFlow/src/panfamflow/cli.py?ref=<branch> -> 200
```

### Finding A2 — nested workflow files are inert in a monorepo

A file at `PanFamFlow/.github/workflows/*.yml` is not executed by GitHub Actions because workflows are only discovered at the repository root `.github/workflows/`.

Corrective action: install `/.github/workflows/panfamflow-ci.yml` and remove misleading nested workflow copies.

## Scientific scope and terminology

### Finding B1 — project wording drifted toward general pangenome analysis

The original project objective is target pan-gene-family analysis, not pangenome assembly or whole-genome HOG classification.

Corrective action:

- title and package description changed to target pan-gene-family analysis;
- canonical module renamed from `pangenome` to `pan_family`;
- `whole_genome` scope removed and rejected;
- output directory renamed to `06_pan_family/`;
- all HOG rows are intersected with target-family stable IDs.

### Finding B2 — OGG, HOG, clade and pan-locus were at risk of conflation

A family-tree clade is not automatically an orthogroup; a HOG is not automatically a material-level syntenic pan-locus.

Corrective action:

- canonical output uses `HOG_ID`, not `OGG_ID`;
- family `subfamily` remains a separate field;
- docs state that pan-locus inference is not implemented in v0.1.1;
- no gene-loss conclusion is generated automatically.

### Finding B3 — family genes absent from the selected HOG node could disappear silently

A stable-ID mismatch or inappropriate HOG node could cause target-family members to be omitted without an explicit audit table.

Corrective action: write `06_pan_family/unassigned_family_members.tsv` with a reason and selected node.

### Finding B4 — annotation absence was too easy to overinterpret

Binary HOG occupancy is influenced by assembly and annotation quality.

Corrective action: each classification row now records that occupancy is annotation/HOG based and not genome-rescued. The report cannot label zeroes as confirmed gene losses.

## Workflow engineering

### Finding C1 — retry/resume claims must be tested against actual outputs

Command-construction unit tests alone do not prove skip/resume behavior.

Required CI gates:

1. unit tests and type/lint checks;
2. Snakemake lint and dry-run;
3. a real toy `qc` execution;
4. a second identical execution with unchanged output hashes/mtimes;
5. deletion of one downstream marker followed by local reconstruction only.

### Finding C2 — list-valued Snakemake inputs were handled as scalar paths in Ka/Ks

Named list inputs can be list-like, so `Path(snakemake.input.membership)` is not universally safe.

Corrective action: normalize scalar/list input values before creating `Path` objects.

### Finding C3 — cached Ka/Ks results did not report cache reuse

Corrective action: set `resumed_from_cache = true` when a validated pair cache is returned.

### Finding C4 — provenance records only described launch intent

Corrective action: run records now start as `RUNNING` and are atomically finalized with `COMPLETED`/`FAILED`, exit code and finish time. A failed launcher invocation is recorded as exit code 127. Per-rule failure summaries remain a future enhancement.

## Biological limitations still open

The following are not represented as completed features:

- genome-level rescue of apparent absent family genes;
- syntenic pan-locus construction for intraspecific accessions;
- DESeq2/edgeR inferential expression analysis;
- codeml/HyPhy positive-selection models;
- systematic matched-background motif enrichment;
- full-scale OrthoFinder/IQ-TREE/KaKs interruption tests on real rice data.

These remain explicit roadmap items rather than implied capabilities.
