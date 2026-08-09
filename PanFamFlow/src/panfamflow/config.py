"""Strict configuration model and module-aware input validation."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Annotated, Any, Literal

import yaml
from pydantic import BaseModel, ConfigDict, Field, StringConstraints, model_validator

from panfamflow.modules import DEFAULT_MODULES, resolve_modules

StrictId = Annotated[
    str,
    StringConstraints(pattern=r"^[A-Za-z][A-Za-z0-9_.-]*$", min_length=1, max_length=80),
]

ResumeMode = Literal["smart", "mtime_only", "off"]
RerunTrigger = Literal["mtime", "input", "params", "code", "software-env"]


def _default_rerun_triggers() -> list[Literal["mtime", "input", "params", "code", "software-env"]]:
    return ["mtime", "input", "params", "code", "software-env"]


class StrictModel(BaseModel):
    model_config = ConfigDict(extra="forbid", validate_assignment=True)


class ProjectSettings(StrictModel):
    # Hard scope guard: PanFamFlow analyzes one configured gene family across
    # already assembled and annotated genomes. It is not a genome-assembly or
    # whole-genome pangenome-construction workflow.
    analysis_scope: Literal["target_pan_gene_family"] = "target_pan_gene_family"
    name: StrictId = "gene_family_project"
    root: Path = Path(".")
    seed: int = 20260807
    results_dir: Path = Path("results")
    work_dir: Path = Path("work")
    logs_dir: Path = Path("logs")


class RunSettings(StrictModel):
    modules: list[str] = Field(default_factory=lambda: list(DEFAULT_MODULES))
    cores: int = Field(default=16, ge=1)
    jobs: int = Field(default=16, ge=1)
    engine_runner: Literal["auto", "mamba", "conda", "current"] = "auto"
    engine_env: str | None = "panfamflow-engine"
    use_conda: bool = True
    # Long-running genomics analyses should resume safely by default.  These
    # settings only affect execution and are excluded from the biological
    # analysis fingerprint written before each run.
    resume_mode: ResumeMode = "smart"
    keep_going: bool = True
    rerun_incomplete: bool = True
    latency_wait: int = Field(default=120, ge=0)
    retries: int = Field(default=1, ge=0)
    rerun_triggers: list[RerunTrigger] = Field(default_factory=_default_rerun_triggers)
    printshellcmds: bool = True
    show_failed_logs: bool = True
    profile: Path | None = None
    extra_snakemake_args: list[str] = Field(default_factory=list)

    @model_validator(mode="after")
    def validate_resume_settings(self) -> RunSettings:
        if len(self.rerun_triggers) != len(set(self.rerun_triggers)):
            raise ValueError("run.rerun_triggers contains duplicates")
        if not self.rerun_triggers:
            raise ValueError("run.rerun_triggers cannot be empty")
        return self


class SpeciesInput(StrictModel):
    id: StrictId
    name: str
    genome: Path
    gff3: Path
    protein: Path | None = None
    cds: Path | None = None
    group: str | None = None
    subfamily: str | None = None
    representative: bool = False
    outgroup: StrictId | None = None
    busco_lineage: str | None = None

    @model_validator(mode="after")
    def reject_internal_delimiter(self) -> SpeciesInput:
        if "__" in self.id:
            raise ValueError("Species IDs must not contain '__'; it is reserved for stable IDs.")
        return self


class RNASeqSample(StrictModel):
    id: StrictId
    species_id: StrictId
    condition: str
    tissue: str | None = None
    stress_type: Literal["Control", "Abiotic", "Biotic", "Other"] = "Other"
    timepoint: str | None = None
    replicate: str | int | None = None
    batch: str | None = None
    strandedness: Literal["unstranded", "forward", "reverse"] = "unstranded"
    r1: Path
    r2: Path | None = None


class InputsSettings(StrictModel):
    species: list[SpeciesInput] = Field(min_length=1)
    rnaseq_samples: list[RNASeqSample] = Field(default_factory=list)
    expression_matrix: Path | None = None
    sample_metadata: Path | None = None

    @model_validator(mode="after")
    def validate_references(self) -> InputsSettings:
        species_ids = [species.id for species in self.species]
        duplicates = sorted({item for item in species_ids if species_ids.count(item) > 1})
        if duplicates:
            raise ValueError(f"Duplicate species IDs: {', '.join(duplicates)}")
        known = set(species_ids)
        for species in self.species:
            if species.outgroup is not None and species.outgroup not in known:
                raise ValueError(
                    f"Species {species.id!r} references unknown outgroup {species.outgroup!r}."
                )
            if species.outgroup == species.id:
                raise ValueError(f"Species {species.id!r} cannot be its own outgroup.")
        sample_ids = [sample.id for sample in self.rnaseq_samples]
        if len(sample_ids) != len(set(sample_ids)):
            raise ValueError("RNA-seq sample IDs must be unique.")
        for sample in self.rnaseq_samples:
            if sample.species_id not in known:
                raise ValueError(
                    f"RNA-seq sample {sample.id!r} references unknown species "
                    f"{sample.species_id!r}."
                )
        return self


class BuscoSettings(StrictModel):
    enabled: bool = False
    mode: Literal["genome"] = "genome"
    threads: int = Field(default=8, ge=1)
    offline: bool = False
    download_path: Path | None = None
    extra_args: list[str] = Field(default_factory=list)


class QcSettings(StrictModel):
    calculate_sha256: bool = True
    busco: BuscoSettings = Field(default_factory=BuscoSettings)


class CanonicalTranscriptSettings(StrictModel):
    method: Literal["longest_cds"] = "longest_cds"
    sequence_source: Literal["gffread"] = "gffread"
    stable_id_separator: str = "__"

    @model_validator(mode="after")
    def validate_separator(self) -> CanonicalTranscriptSettings:
        if not self.stable_id_separator or any(char.isspace() for char in self.stable_id_separator):
            raise ValueError("canonical_transcript.stable_id_separator must be non-empty/no-space.")
        return self


class HMMSettings(StrictModel):
    enabled: bool = False
    hmm: Path | None = None
    evalue: float = Field(default=1.0e-5, gt=0)
    domain_evalue: float = Field(default=1.0e-3, gt=0)
    cut_ga: bool = False


class BlastSettings(StrictModel):
    enabled: bool = False
    reference_proteins: Path | None = None
    evalue: float = Field(default=1.0e-5, gt=0)
    min_identity: float = Field(default=30.0, ge=0, le=100)
    min_query_coverage: float = Field(default=50.0, ge=0, le=100)
    max_target_seqs: int = Field(default=100, ge=1)


class FamilySettings(StrictModel):
    name: StrictId = "TARGET_FAMILY"
    combine_evidence: Literal["union", "intersection", "hmm_only", "blast_only"] = "union"
    hmm: HMMSettings = Field(default_factory=HMMSettings)
    blast: BlastSettings = Field(default_factory=BlastSettings)
    calculate_protein_properties: bool = True
    subfamily_assignments: Path | None = None
    domain_validation_table: Path | None = None
    subcellular_localization_table: Path | None = None
    precomputed_members: Path | None = None

    @model_validator(mode="after")
    def validate_evidence(self) -> FamilySettings:
        if self.precomputed_members is None and not (self.hmm.enabled or self.blast.enabled):
            raise ValueError(
                "Family identification requires HMM, BLAST, or family.precomputed_members."
            )
        if self.hmm.enabled and self.hmm.hmm is None:
            raise ValueError("family.hmm.hmm is required when HMM search is enabled.")
        if self.blast.enabled and self.blast.reference_proteins is None:
            raise ValueError(
                "family.blast.reference_proteins is required when BLAST search is enabled."
            )
        if self.combine_evidence == "intersection" and not (
            self.hmm.enabled and self.blast.enabled
        ):
            raise ValueError("intersection evidence requires both HMM and BLAST to be enabled.")
        if self.combine_evidence == "hmm_only" and not self.hmm.enabled:
            raise ValueError("hmm_only evidence requires HMM to be enabled.")
        if self.combine_evidence == "blast_only" and not self.blast.enabled:
            raise ValueError("blast_only evidence requires BLAST to be enabled.")
        return self


class PhylogenySettings(StrictModel):
    mafft_mode: Literal["auto", "linsi", "ginsi", "einsi"] = "auto"
    trim_mode: str = "smart-gap"
    model: str = "MFP"
    ultrafast_bootstrap: int = Field(default=1000, ge=0)
    sh_alrt: int = Field(default=1000, ge=0)
    min_sequences: int = Field(default=4, ge=3)


class OrthoFinderSettings(StrictModel):
    hog_node: str = "auto"
    search_threads: int = Field(default=32, ge=1)
    analysis_threads: int = Field(default=8, ge=1)
    extra_args: list[str] = Field(default_factory=list)


class PanFamilySettings(StrictModel):
    """Occupancy classification for HOGs containing target-family members only."""

    core_min: float = Field(default=0.99, gt=0, le=1)
    soft_core_min: float = Field(default=0.90, gt=0, le=1)
    shell_min: float = Field(default=0.10, gt=0, le=1)
    rarefaction_iterations: int = Field(default=1000, ge=1)
    max_exact_combinations: int = Field(default=5000, ge=1)

    @model_validator(mode="after")
    def validate_threshold_order(self) -> PanFamilySettings:
        if not self.core_min >= self.soft_core_min >= self.shell_min:
            raise ValueError(
                "Pan-family thresholds must satisfy core_min >= soft_core_min >= shell_min."
            )
        return self


class ChromosomeSettings(StrictModel):
    representative_only: bool = False
    density_window_bp: int = Field(default=1_000_000, ge=1)


class DuplicationSettings(StrictModel):
    backend: Literal["dupgen_finder_unique", "precomputed"] = "dupgen_finder_unique"
    targets: list[StrictId] | None = None
    precomputed_table: Path | None = None
    dupgen_executable: str = "DupGen_finder-unique.pl"
    diamond_evalue: float = Field(default=1.0e-10, gt=0)
    max_target_seqs: int = Field(default=5, ge=1)
    proximal_max_gene_distance: int = Field(default=10, ge=1)
    extra_args: list[str] = Field(default_factory=list)


class KaksSettings(StrictModel):
    pair_source: Literal["orthology", "duplication", "both"] = "both"
    reference_species: StrictId | None = None
    method: str = "MA"
    max_pairs_per_group: int | None = Field(default=None, ge=1)
    saturation_ks: float = Field(default=2.0, gt=0)
    workers: int = Field(default=4, ge=1)


class PromoterSettings(StrictModel):
    backend: Literal["fimo", "precomputed_plantcare"] = "fimo"
    upstream_bp: int = Field(default=2000, ge=1)
    downstream_bp: int = Field(default=0, ge=0)
    motif_database: Path | None = None
    category_map: Path | None = None
    precomputed_table: Path | None = None
    fimo_threshold: float = Field(default=1.0e-4, gt=0, le=1)
    top_n_elements: int = Field(default=20, ge=1)


class ExpressionSettings(StrictModel):
    mode: Literal["imported_matrix", "fastq_stringtie"] = "imported_matrix"
    min_tpm_detected: float = Field(default=1.0, ge=0)
    heatmap_transform: Literal["log2_tpm1_zscore", "log2_tpm1"] = "log2_tpm1_zscore"
    fastp_extra_args: list[str] = Field(default_factory=list)
    hisat2_extra_args: list[str] = Field(default_factory=list)
    stringtie_extra_args: list[str] = Field(default_factory=list)


class PlotSettings(StrictModel):
    pdf: bool = True
    png: bool = True
    png_dpi: int = Field(default=600, ge=72)
    language: Literal["English"] = "English"


class ReportSettings(StrictModel):
    title: str | None = None
    include_existing_results: bool = True


class WorkflowConfig(StrictModel):
    @model_validator(mode="before")
    @classmethod
    def migrate_pre_release_names(cls, value: Any) -> Any:
        """Migrate unambiguous pre-release names while rejecting scope drift.

        PanFamFlow is a target pan-gene-family workflow.  It intentionally
        rejects the former ``pangenome.scope: whole_genome`` option because
        that belongs to a whole-pangenome analysis/assembly project, not this
        pipeline.
        """

        if not isinstance(value, dict):
            return value
        data = dict(value)
        if "pangenome" in data:
            if "pan_family" in data:
                raise ValueError(
                    "Configure only pan_family; do not provide both pan_family and pangenome."
                )
            legacy = dict(data.pop("pangenome") or {})
            scope = legacy.pop("scope", "target_family")
            if scope != "target_family":
                raise ValueError(
                    "PanFamFlow does not support whole-genome pangenome analysis or pangenome "
                    "assembly. Use pan_family for the configured target gene family only."
                )
            data["pan_family"] = legacy
        run = data.get("run")
        if isinstance(run, dict) and isinstance(run.get("modules"), list):
            migrated_run = dict(run)
            migrated_run["modules"] = [
                "pan_family" if str(item).strip().lower().replace("-", "_") == "pangenome" else item
                for item in run["modules"]
            ]
            data["run"] = migrated_run
        return data

    schema_version: Literal["1.0"] = "1.0"
    project: ProjectSettings = Field(default_factory=ProjectSettings)
    run: RunSettings = Field(default_factory=RunSettings)
    inputs: InputsSettings
    qc: QcSettings = Field(default_factory=QcSettings)
    canonical_transcript: CanonicalTranscriptSettings = Field(
        default_factory=CanonicalTranscriptSettings
    )
    family: FamilySettings
    phylogeny: PhylogenySettings = Field(default_factory=PhylogenySettings)
    orthofinder: OrthoFinderSettings = Field(default_factory=OrthoFinderSettings)
    pan_family: PanFamilySettings = Field(default_factory=PanFamilySettings)
    chromosome: ChromosomeSettings = Field(default_factory=ChromosomeSettings)
    duplication: DuplicationSettings = Field(default_factory=DuplicationSettings)
    kaks: KaksSettings = Field(default_factory=KaksSettings)
    promoter: PromoterSettings = Field(default_factory=PromoterSettings)
    expression: ExpressionSettings = Field(default_factory=ExpressionSettings)
    plot: PlotSettings = Field(default_factory=PlotSettings)
    report: ReportSettings = Field(default_factory=ReportSettings)

    @model_validator(mode="after")
    def validate_cross_references(self) -> WorkflowConfig:
        known = {species.id for species in self.inputs.species}
        if self.duplication.targets is not None:
            unknown = sorted(set(self.duplication.targets).difference(known))
            if unknown:
                raise ValueError(f"Unknown duplication.targets: {', '.join(unknown)}")
        if self.kaks.reference_species is not None and self.kaks.reference_species not in known:
            raise ValueError(
                f"kaks.reference_species {self.kaks.reference_species!r} is not in inputs.species."
            )
        return self


def _fingerprint_payload(payload: Any) -> str:
    """Return a deterministic SHA256 fingerprint for JSON-compatible data."""

    encoded = json.dumps(
        payload,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def analysis_config_payload(config: WorkflowConfig) -> dict[str, Any]:
    """Return fields that can change biological or statistical results.

    Resource allocation, engine selection, output locations and rendering-only
    settings are deliberately excluded.  Input file content is audited
    separately by the workflow, while configured paths remain part of this
    payload for provenance and change detection.
    """

    payload = config.model_dump(mode="json", exclude_none=False)
    project = payload.get("project", {})
    payload["project"] = {
        "analysis_scope": project.get("analysis_scope"),
        "seed": project.get("seed"),
    }
    payload.pop("run", None)
    payload.pop("plot", None)
    payload.pop("report", None)
    return payload


def execution_config_payload(config: WorkflowConfig) -> dict[str, Any]:
    """Return execution and presentation settings excluded from biology hash."""

    project = config.project.model_dump(mode="json", exclude_none=False)
    return {
        "project": {
            "analysis_scope": project["analysis_scope"],
            "name": project["name"],
            "root": project["root"],
            "results_dir": project["results_dir"],
            "work_dir": project["work_dir"],
            "logs_dir": project["logs_dir"],
        },
        "run": config.run.model_dump(mode="json", exclude_none=False),
        "plot": config.plot.model_dump(mode="json", exclude_none=False),
        "report": config.report.model_dump(mode="json", exclude_none=False),
    }


def analysis_config_hash(config: WorkflowConfig) -> str:
    """Return the biological/statistical configuration fingerprint."""

    return _fingerprint_payload(analysis_config_payload(config))


def execution_config_hash(config: WorkflowConfig) -> str:
    """Return the execution/presentation configuration fingerprint."""

    return _fingerprint_payload(execution_config_payload(config))


@dataclass(frozen=True, slots=True)
class ValidationIssue:
    severity: Literal["ERROR", "WARNING"]
    field: str
    message: str


def load_config(path: Path) -> WorkflowConfig:
    """Load and strictly validate a YAML configuration file."""

    config_path = path.expanduser().resolve()
    if not config_path.is_file():
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    with config_path.open("r", encoding="utf-8") as handle:
        raw = yaml.safe_load(handle)
    if not isinstance(raw, dict):
        raise ValueError(f"Configuration root must be a YAML mapping: {config_path}")
    return WorkflowConfig.model_validate(raw)


def project_root(config: WorkflowConfig, config_path: Path) -> Path:
    """Resolve the configured project root relative to the config file directory."""

    root = config.project.root.expanduser()
    if not root.is_absolute():
        root = config_path.expanduser().resolve().parent / root
    return root.resolve()


def resolve_project_path(path: Path | None, root: Path) -> Path | None:
    """Resolve an optional project-relative path."""

    if path is None:
        return None
    expanded = path.expanduser()
    return expanded.resolve() if expanded.is_absolute() else (root / expanded).resolve()


def validate_input_paths(
    config: WorkflowConfig,
    config_path: Path,
    requested_modules: list[str] | tuple[str, ...] | None = None,
) -> tuple[ValidationIssue, ...]:
    """Validate only inputs needed by the requested module closure."""

    modules = resolve_modules(requested_modules or config.run.modules, config)
    selected = set(modules)
    root = project_root(config, config_path)
    issues: list[ValidationIssue] = []

    def require(path: Path | None, field: str) -> None:
        if path is None:
            issues.append(ValidationIssue("ERROR", field, "Required path is not configured."))
            return
        resolved = resolve_project_path(path, root)
        assert resolved is not None
        if not resolved.is_file():
            issues.append(ValidationIssue("ERROR", field, f"File not found: {resolved}"))
        elif resolved.stat().st_size == 0:
            issues.append(ValidationIssue("ERROR", field, f"File is empty: {resolved}"))

    data_modules = selected.difference({"report"})
    if data_modules:
        for index, species in enumerate(config.inputs.species):
            require(species.genome, f"inputs.species[{index}].genome")
            require(species.gff3, f"inputs.species[{index}].gff3")

    if "qc" in selected and config.qc.busco.enabled:
        for index, species in enumerate(config.inputs.species):
            if not species.busco_lineage:
                issues.append(
                    ValidationIssue(
                        "ERROR",
                        f"inputs.species[{index}].busco_lineage",
                        "BUSCO is enabled; set an appropriate lineage for each species.",
                    )
                )
        if config.qc.busco.offline:
            if config.qc.busco.download_path is None:
                issues.append(
                    ValidationIssue(
                        "ERROR",
                        "qc.busco.download_path",
                        "BUSCO offline mode requires a configured download_path.",
                    )
                )
            else:
                resolved = resolve_project_path(config.qc.busco.download_path, root)
                assert resolved is not None
                if not resolved.is_dir():
                    issues.append(
                        ValidationIssue(
                            "ERROR",
                            "qc.busco.download_path",
                            f"Directory not found: {resolved}",
                        )
                    )

    if "family" in selected and config.family.precomputed_members is not None:
        require(config.family.precomputed_members, "family.precomputed_members")
    elif "family" in selected:
        if config.family.hmm.enabled:
            require(config.family.hmm.hmm, "family.hmm.hmm")
        if config.family.blast.enabled:
            require(config.family.blast.reference_proteins, "family.blast.reference_proteins")
    for field, path in (
        ("family.subfamily_assignments", config.family.subfamily_assignments),
        ("family.domain_validation_table", config.family.domain_validation_table),
        ("family.subcellular_localization_table", config.family.subcellular_localization_table),
    ):
        if "family" in selected and path is not None:
            require(path, field)

    if "orthology" in selected and len(config.inputs.species) < 2:
        issues.append(
            ValidationIssue("ERROR", "inputs.species", "OrthoFinder requires at least two species.")
        )

    if "pan_family" in selected and config.orthofinder.hog_node.lower() == "auto":
        issues.append(
            ValidationIssue(
                "WARNING",
                "orthofinder.hog_node",
                "auto is acceptable for discovery, but a final analysis should record the "
                "target-clade N* node from SpeciesTree_rooted_node_labels.txt.",
            )
        )

    if "duplication" in selected:
        if config.duplication.backend == "precomputed":
            require(config.duplication.precomputed_table, "duplication.precomputed_table")
        else:
            targets = config.duplication.targets or [
                species.id for species in config.inputs.species if species.outgroup is not None
            ]
            if not targets:
                issues.append(
                    ValidationIssue(
                        "ERROR",
                        "duplication.targets",
                        "DupGen_finder-unique requires at least one target with an outgroup.",
                    )
                )
            species_by_id = {species.id: species for species in config.inputs.species}
            for target in targets:
                if species_by_id[target].outgroup is None:
                    issues.append(
                        ValidationIssue(
                            "ERROR",
                            f"inputs.species[{target}].outgroup",
                            "Each DupGen_finder target requires a configured outgroup.",
                        )
                    )

    if "promoter" in selected:
        if config.promoter.backend == "fimo":
            require(config.promoter.motif_database, "promoter.motif_database")
            if config.promoter.category_map is None:
                issues.append(
                    ValidationIssue(
                        "WARNING",
                        "promoter.category_map",
                        "No curated motif-to-category map is configured; major-class summaries "
                        "will remain Unclassified.",
                    )
                )
            else:
                require(config.promoter.category_map, "promoter.category_map")
        else:
            require(config.promoter.precomputed_table, "promoter.precomputed_table")

    if "expression" in selected:
        if config.expression.mode == "imported_matrix":
            require(config.inputs.expression_matrix, "inputs.expression_matrix")
            if config.inputs.sample_metadata is not None:
                require(config.inputs.sample_metadata, "inputs.sample_metadata")
        else:
            if not config.inputs.rnaseq_samples:
                issues.append(
                    ValidationIssue(
                        "ERROR",
                        "inputs.rnaseq_samples",
                        "fastq_stringtie mode requires at least one RNA-seq sample.",
                    )
                )
            for index, sample in enumerate(config.inputs.rnaseq_samples):
                require(sample.r1, f"inputs.rnaseq_samples[{index}].r1")
                if sample.r2 is not None:
                    require(sample.r2, f"inputs.rnaseq_samples[{index}].r2")
            biological_groups: dict[tuple[str, str, str | None], int] = {}
            for sample in config.inputs.rnaseq_samples:
                key = (sample.species_id, sample.condition, sample.tissue)
                biological_groups[key] = biological_groups.get(key, 0) + 1
            if biological_groups and min(biological_groups.values()) < 2:
                issues.append(
                    ValidationIssue(
                        "WARNING",
                        "inputs.rnaseq_samples",
                        "At least one species/condition/tissue group has fewer than two samples. "
                        "TPM quantification remains possible, but inferential differential "
                        "expression is not implemented by this v0.1 module.",
                    )
                )

    if config.run.profile is not None:
        profile = resolve_project_path(config.run.profile, root)
        assert profile is not None
        if not profile.exists():
            issues.append(ValidationIssue("ERROR", "run.profile", f"Profile not found: {profile}"))

    return tuple(issues)
