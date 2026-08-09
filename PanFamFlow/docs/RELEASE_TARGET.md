# Release target decision

PanFamFlow is developed on this feature branch for traceability, but the preferred public release target is a standalone repository named `PanFamFlow`.

Rationale:

- the software analyzes a configured target gene family across assembled and annotated genomes;
- it is not a genome-assembly or graph-pangenome workflow;
- a standalone repository prevents the host repository title and unrelated directories from implying a different analytical scope;
- releases, issues, citations and CI should belong to the software rather than to the source pangenome project.

The current branch remains the auditable development source until the standalone repository has been created and verified. The real-data biological benchmark remains a separate release gate.