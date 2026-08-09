# Changelog

## 0.1.1 - 2026-08-09
- Added immutable `project.analysis_scope: target_pan_gene_family` validation and provenance fingerprinting.

- Corrected project scope to target pan-gene-family analysis; explicitly excluded pangenome assembly and whole-genome HOG classification.
- Renamed the canonical `pangenome` module to `pan_family` and outputs to `06_pan_family/`; retained a narrow pre-release alias for `scope: target_family` only.
- Rejected the former `whole_genome` scope at configuration validation.
- Replaced ambiguous `OGG_ID`/`pangenome_class` output fields with `HOG_ID`/`pan_family_class`.
- Added `unassigned_family_members.tsv` so stable-ID or HOG-node mismatches cannot silently drop target-family genes.
- Added explicit occupancy/absence interpretation flags to prevent annotation absence from being reported as validated gene loss.
- Fixed list-valued Ka/Ks Snakemake input handling and cache-reuse reporting.
- Added smart resume defaults with rerun-incomplete, keep-going, retries and explicit rerun triggers.
- Added `status`, `resume` and `retry` CLI commands.
- Added analysis/execution configuration fingerprints and atomic provenance snapshots.
- Added atomic file publication to prevent partial outputs from being mistaken for completed results.
- Added IQ-TREE checkpoint reuse without unconditional `--redo`.
- Added signature-based OrthoFinder work directories and WorkingDirectory reuse.
- Added stable per-pair Ka/Ks caches and non-destructive BUSCO run directories.
- Added failure-injection and resume regression tests.

## 0.1.0 - 2026-08-08

- Initial public alpha implementation.
- Added one-file YAML configuration validated with Pydantic.
- Added a Python CLI for initialization, validation, planning, selective module execution, diagnostics and schema export.
- Added a Snakemake 9 workflow with rule-specific Conda environments.
- Added input audit, optional BUSCO, canonical transcript generation, family discovery, phylogeny, gene structure, OrthoFinder/HOG target pan-family classification, chromosome summaries, duplication, Ka/Ks, promoter/FIMO, expression and report modules.
- Added toy data, unit tests, CI and local/SLURM profiles.
