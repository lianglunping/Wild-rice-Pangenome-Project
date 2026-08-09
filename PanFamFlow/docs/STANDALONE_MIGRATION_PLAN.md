# Standalone PanFamFlow migration plan

The feature branch is now directly usable and auditable. The long-term public distribution target should be a standalone `lianglunping/PanFamFlow` repository so that the repository name and top-level content cannot be mistaken for a genome-assembly or graph-pangenome project.

The standalone repository should contain only:

```text
README.md
LICENSE
CITATION.cff
pyproject.toml
uv.lock
src/
tests/
examples/
docs/
profiles/
scripts/
.github/workflows/
```

Migration must preserve the exact validated source commit, durable CI, scope guards and smart-resume integration test. It must not include the host repository's unrelated genome assembly, TE, variant calling, domestication or whole-project pangenome directories.

Creation of the standalone repository is an external repository-management action and is kept separate from source correction and biological validation.