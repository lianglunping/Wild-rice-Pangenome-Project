# Independent third-party audit

## Audit perspective

This audit evaluates the repository as a new external user would encounter it. It is not limited to the previously reported 404 screenshot.

## Resolved blocking defects

1. `PanFamFlow/` is now a directly browsable source directory rather than an unpublished archive or staging placeholder.
2. Temporary publication directories and one-time bootstrap/recovery workflows are not part of the usable source tree.
3. The canonical title is **PanFamFlow — 泛基因家族分析流程**.
4. The canonical module is `pan_family`; the canonical output directory is `06_pan_family`.
5. The canonical configuration uses `pan_family:` and rejects a whole-genome analysis scope.
6. The workflow consumes already assembled and annotated genomes or materials. Genome assembly, graph-pangenome construction, genome-wide SV/PAV calling and general whole-genome pangenome catalog construction are outside scope.
7. OrthoFinder may use complete proteomes to infer orthology, but downstream occupancy and Core/Soft-core/Shell/Cloud classification are restricted to the configured target-family members.
8. Target-family members missing from the selected HOG table are retained in an explicit unassigned-members audit rather than being silently discarded.
9. Annotation/HOG absence is not labelled as validated gene loss without genome-level rescue evidence.
10. Smart resume, incomplete-job rerun, valid-output skipping, provenance fingerprints and atomic outputs are part of the software contract.

## Remaining biological validation boundary

Repository checks, unit tests, toy DAG execution and resume smoke tests do not establish biological validity. Before production use, one manually reviewed rice gene family should be benchmarked across 5–10 high-quality assembled and annotated genomes. The benchmark must review family-member precision/recall, selected HOG node, unassigned members, genome-level absence rescue, duplication assignments, constrained Ka/Ks pairs and expression identifiers.

## Hosting observation

The current feature branch is hosted inside the pre-existing `Wild-rice-Pangenome-Project` repository. This preserves development history but can confuse users because the host repository contains unrelated whole-project pangenome code. A standalone `PanFamFlow` repository is the preferred release target once the real-data benchmark is passed. Until then, the root branch README must keep a prominent scope statement and link directly to `PanFamFlow/`.
