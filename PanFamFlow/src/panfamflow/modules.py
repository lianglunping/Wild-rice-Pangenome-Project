"""Workflow module registry and dependency resolution.

PanFamFlow analyses one configured target gene family across multiple assembled
and annotated genomes.  It does not construct a graph pangenome or assemble a
pangenome reference.  The ``pan_family`` module classifies target-family HOGs
by sample/species occupancy.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from panfamflow.config import WorkflowConfig


@dataclass(frozen=True, slots=True)
class ModuleSpec:
    """Metadata for one user-selectable workflow module."""

    name: str
    target: str
    description: str
    direct_dependencies: tuple[str, ...]
    primary_tools: tuple[str, ...]


MODULE_ALIASES: dict[str, str] = {
    # v0.1.0/v0.1.1 pre-release compatibility.  The old name was ambiguous
    # because it sounded like whole-pangenome construction.
    "pangenome": "pan_family",
    "pan_gene_family": "pan_family",
    "panfamily": "pan_family",
}


MODULES: dict[str, ModuleSpec] = {
    "qc": ModuleSpec(
        "qc",
        "00_qc/qc.done",
        "Audit configured assembled genomes, annotations, SHA256 and optional BUSCO results.",
        (),
        ("Python", "BUSCO (optional)"),
    ),
    "normalize": ModuleSpec(
        "normalize",
        "01_normalized/normalized.done",
        "Select one longest-CDS transcript per gene and regenerate canonical sequences.",
        ("qc",),
        ("AGAT", "gffread", "Python"),
    ),
    "family": ModuleSpec(
        "family",
        "02_family/family_members.tsv",
        "Identify target-family members by HMMER and/or BLASTP with evidence auditing.",
        ("normalize",),
        ("HMMER", "BLAST+", "Biopython"),
    ),
    "phylogeny": ModuleSpec(
        "phylogeny",
        "03_phylogeny/family.treefile",
        "Infer a maximum-likelihood tree for the configured target family.",
        ("family",),
        ("MAFFT", "ClipKIT", "IQ-TREE"),
    ),
    "gene_structure": ModuleSpec(
        "gene_structure",
        "04_gene_structure/gene_structure_metrics.tsv",
        "Extract gene, CDS, exon, intron and UTR metrics for target-family members.",
        ("family",),
        ("Python",),
    ),
    "orthology": ModuleSpec(
        "orthology",
        "05_orthology/orthofinder.done",
        "Infer HOG context from canonical proteomes; downstream results are filtered to the target family.",
        ("normalize",),
        ("OrthoFinder",),
    ),
    "pan_family": ModuleSpec(
        "pan_family",
        "06_pan_family/pan_family_classification.tsv",
        "Classify target-family HOG occupancy and compute target-family rarefaction curves.",
        ("orthology", "family"),
        ("Python", "OrthoFinder HOGs"),
    ),
    "chromosome": ModuleSpec(
        "chromosome",
        "07_chromosome/chromosome_distribution.tsv",
        "Summarize target-family chromosome coordinates, counts and density.",
        ("family",),
        ("Python",),
    ),
    "duplication": ModuleSpec(
        "duplication",
        "08_duplication/duplication_mode.tsv",
        "Classify duplication modes for target-family genes or import audited results.",
        ("normalize", "family"),
        ("DIAMOND", "DupGen_finder-unique"),
    ),
    "kaks": ModuleSpec(
        "kaks",
        "09_kaks/kaks_pairs.tsv",
        "Estimate pairwise Ka, Ks and Ka/Ks for constrained target-family pairs.",
        ("family",),
        ("MAFFT", "PAL2NAL", "KaKs_Calculator"),
    ),
    "promoter": ModuleSpec(
        "promoter",
        "10_promoter/promoter_elements.tsv",
        "Extract promoters of target-family genes and scan a versioned motif database.",
        ("family",),
        ("MEME Suite/FIMO", "Python"),
    ),
    "expression": ModuleSpec(
        "expression",
        "11_expression/expression_matrix.tsv",
        "Import or quantify expression for target-family genes across configured samples.",
        ("family", "normalize"),
        ("fastp", "HISAT2", "StringTie", "Python"),
    ),
    "report": ModuleSpec(
        "report",
        "report/index.html",
        "Build the target-family master table, result manifest and static HTML report.",
        (),
        ("Python",),
    ),
}

DEFAULT_MODULES: tuple[str, ...] = (
    "qc",
    "normalize",
    "family",
    "phylogeny",
    "gene_structure",
    "orthology",
    "pan_family",
    "chromosome",
    "report",
)


def module_names() -> tuple[str, ...]:
    """Return canonical module names in stable execution order."""

    return tuple(MODULES)


def canonical_module_name(name: str) -> str:
    """Return the canonical module name, accepting documented legacy aliases."""

    normalized = name.strip().lower().replace("-", "_")
    return MODULE_ALIASES.get(normalized, normalized)


def _dynamic_dependencies(name: str, config: WorkflowConfig) -> tuple[str, ...]:
    dependencies = list(MODULES[name].direct_dependencies)
    if name == "kaks":
        source = config.kaks.pair_source
        if source in {"orthology", "both"}:
            dependencies.append("pan_family")
        if source in {"duplication", "both"}:
            dependencies.append("duplication")
    return tuple(dict.fromkeys(dependencies))


def normalize_requested_modules(modules: Iterable[str]) -> tuple[str, ...]:
    """Normalize names, resolve aliases, reject unknown modules and preserve order."""

    requested = [canonical_module_name(item) for item in modules if item.strip()]
    if not requested:
        raise ValueError("At least one module must be selected.")
    if "all" in requested:
        requested = list(MODULES)
    unknown = sorted(set(requested).difference(MODULES))
    if unknown:
        valid = ", ".join(MODULES)
        aliases = ", ".join(f"{old}->{new}" for old, new in sorted(MODULE_ALIASES.items()))
        raise ValueError(
            f"Unknown module(s): {', '.join(unknown)}. Valid modules: {valid}. "
            f"Accepted aliases: {aliases}."
        )
    selected = set(requested)
    return tuple(name for name in MODULES if name in selected)


def resolve_modules(modules: Iterable[str], config: WorkflowConfig) -> tuple[str, ...]:
    """Expand transitive dependencies and return a stable topological order."""

    requested = normalize_requested_modules(modules)
    resolved: set[str] = set()
    active: set[str] = set()

    def visit(name: str) -> None:
        if name in resolved:
            return
        if name in active:
            raise RuntimeError(f"Circular module dependency detected at {name!r}.")
        active.add(name)
        for dependency in _dynamic_dependencies(name, config):
            visit(dependency)
        active.remove(name)
        resolved.add(name)

    for module in requested:
        visit(module)
    return tuple(name for name in MODULES if name in resolved)


def targets_for_modules(modules: Iterable[str], results_dir: str) -> tuple[str, ...]:
    """Return Snakemake target paths for selected canonical modules."""

    root = results_dir.rstrip("/")
    return tuple(f"{root}/{MODULES[canonical_module_name(name)].target}" for name in modules)
