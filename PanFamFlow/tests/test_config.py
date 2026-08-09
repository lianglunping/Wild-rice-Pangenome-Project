from __future__ import annotations

from pathlib import Path

import yaml
from pydantic import ValidationError

from panfamflow.config import WorkflowConfig, load_config, validate_input_paths
from panfamflow.modules import resolve_modules

TOY_CONFIG = Path(__file__).parents[1] / "examples" / "toy" / "config.yaml"


def test_toy_config_loads() -> None:
    config = load_config(TOY_CONFIG)
    assert config.project.name == "toy_gene_family"
    assert [species.id for species in config.inputs.species] == ["SpA", "SpB"]


def test_module_dependency_expansion() -> None:
    config = load_config(TOY_CONFIG)
    assert resolve_modules(["phylogeny"], config) == (
        "qc",
        "normalize",
        "family",
        "phylogeny",
    )


def test_kaks_dynamic_dependency_expansion() -> None:
    raw = yaml.safe_load(TOY_CONFIG.read_text(encoding="utf-8"))
    raw["kaks"]["pair_source"] = "both"
    config = WorkflowConfig.model_validate(raw)
    resolved = resolve_modules(["kaks"], config)
    assert "pan_family" in resolved
    assert "duplication" in resolved
    assert resolved[-1] == "kaks"


def test_qc_path_validation_passes() -> None:
    config = load_config(TOY_CONFIG)
    issues = validate_input_paths(config, TOY_CONFIG, ["qc"])
    assert not [issue for issue in issues if issue.severity == "ERROR"]


def test_duplicate_species_id_is_rejected() -> None:
    raw = yaml.safe_load(TOY_CONFIG.read_text(encoding="utf-8"))
    raw["inputs"]["species"][1]["id"] = "SpA"
    try:
        WorkflowConfig.model_validate(raw)
    except ValidationError as error:
        assert "Duplicate species IDs" in str(error)
    else:
        raise AssertionError("Duplicate species IDs were accepted")


def test_analysis_and_execution_fingerprints_are_separated() -> None:
    from panfamflow.config import analysis_config_hash, execution_config_hash

    config = load_config(TOY_CONFIG)
    more_cores = config.model_copy(
        update={"run": config.run.model_copy(update={"cores": config.run.cores + 4})}
    )
    assert analysis_config_hash(more_cores) == analysis_config_hash(config)
    assert execution_config_hash(more_cores) != execution_config_hash(config)

    changed_family = config.model_copy(
        update={
            "family": config.family.model_copy(
                update={
                    "hmm": config.family.hmm.model_copy(
                        update={"evalue": config.family.hmm.evalue / 10}
                    )
                }
            )
        }
    )
    assert analysis_config_hash(changed_family) != analysis_config_hash(config)


def test_legacy_pangenome_name_migrates_to_pan_family() -> None:
    raw = yaml.safe_load(TOY_CONFIG.read_text(encoding="utf-8"))
    raw["pangenome"] = {**raw.pop("pan_family"), "scope": "target_family"}
    raw["run"]["modules"] = ["pangenome"]
    config = WorkflowConfig.model_validate(raw)
    assert config.pan_family.core_min == 0.99
    assert resolve_modules(config.run.modules, config)[-1] == "pan_family"


def test_whole_genome_scope_is_rejected() -> None:
    raw = yaml.safe_load(TOY_CONFIG.read_text(encoding="utf-8"))
    raw["pangenome"] = {**raw.pop("pan_family"), "scope": "whole_genome"}
    try:
        WorkflowConfig.model_validate(raw)
    except ValidationError as error:
        assert "does not support whole-genome pangenome analysis" in str(error)
    else:
        raise AssertionError("whole_genome scope was accepted")


def test_analysis_scope_is_hard_guard() -> None:
    raw = yaml.safe_load(TOY_CONFIG.read_text(encoding="utf-8"))
    assert raw["project"]["analysis_scope"] == "target_pan_gene_family"
    raw["project"]["analysis_scope"] = "graph_pangenome"
    try:
        WorkflowConfig.model_validate(raw)
    except ValidationError as error:
        assert "analysis_scope" in str(error)
    else:
        raise AssertionError("Non-family analysis scope was accepted")
