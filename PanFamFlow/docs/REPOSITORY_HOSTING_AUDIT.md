# Repository hosting audit

A third-party user arriving at `Wild-rice-Pangenome-Project` will reasonably expect a whole-project pangenome repository. That host name and its unrelated top-level analysis directories can obscure the narrower PanFamFlow scope.

The current feature branch is acceptable as a development and review location because it preserves the implementation history. It is not the preferred long-term distribution surface.

The preferred release topology is:

```text
lianglunping/PanFamFlow
├── README.md
├── pyproject.toml
├── uv.lock
├── src/panfamflow/
├── tests/
├── examples/
├── docs/
├── profiles/
└── .github/workflows/
```

No genome-assembly, graph-pangenome, whole-genome variant-calling or unrelated domestication-analysis directories should be included in that standalone repository.
