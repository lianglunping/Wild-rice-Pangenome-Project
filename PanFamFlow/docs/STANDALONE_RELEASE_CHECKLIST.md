# Standalone release checklist

- [ ] Create `lianglunping/PanFamFlow` without unrelated pangenome-project directories.
- [ ] Preserve full software source, tests, examples, documentation and durable CI.
- [ ] Set repository description to target pan-gene-family analysis on assembled and annotated genomes.
- [ ] Keep `pan_family` and `06_pan_family` as canonical names.
- [ ] Verify `config.yaml` rejects whole-genome scope.
- [ ] Run locked environment, Ruff, mypy, pytest, build, CLI/DAG and smart-resume integration tests.
- [ ] Add a CI attestation for the exact validated commit.
- [ ] Keep the release pre-1.0 until a real rice-family benchmark is complete.
