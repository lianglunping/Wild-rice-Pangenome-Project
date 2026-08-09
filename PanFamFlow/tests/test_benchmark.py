from __future__ import annotations

import hashlib
from pathlib import Path

import pytest
from typer.testing import CliRunner

from panfamflow.benchmark import (
    audit_benchmark,
    initialize_benchmark,
    write_benchmark_audit,
)
from panfamflow.cli import app

runner = CliRunner()


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_ready_benchmark(tmp_path: Path) -> Path:
    root = tmp_path / "benchmark"
    root.mkdir()
    (root / "references").mkdir()
    (root / "manual_review").mkdir()
    (root / "references/family.hmm").write_text(
        "HMMER3/f\nNAME  TEST\nLENG  1\n//\n",
        encoding="utf-8",
    )
    (root / "references/family_reference.pep.fa").write_text(
        ">RefGene\nMKM\n",
        encoding="utf-8",
    )

    rows: list[dict[str, str]] = []
    for index, group in ((1, "A"), (2, "B")):
        species_id = f"Rice{index:02d}"
        species_dir = root / "inputs" / species_id
        species_dir.mkdir(parents=True)
        genome = species_dir / "genome.fa"
        gff3 = species_dir / "annotation.gff3"
        protein = species_dir / "protein.fa"
        cds = species_dir / "cds.fa"
        genome.write_text(">Chr1\nATGAAGATG\n", encoding="utf-8")
        gff3.write_text(
            "##gff-version 3\nChr1\ttest\tgene\t1\t9\t.\t+\t.\tID=Gene1\n",
            encoding="utf-8",
        )
        protein.write_text(">Gene1\nMKM\n", encoding="utf-8")
        cds.write_text(">Gene1\nATGAAGATG\n", encoding="utf-8")
        rows.append(
            {
                "species_id": species_id,
                "species_name": f"Oryza sativa {index}",
                "group": group,
                "include": "true",
                "representative": "true" if index == 1 else "false",
                "input_kind": "assembled_genome",
                "assembly_accession": f"TEST{index}.1",
                "assembly_level": "chromosome",
                "annotation_version": "v1",
                "coordinate_system": f"TEST{index}.1",
                "genome": str(genome.relative_to(root)),
                "genome_sha256": _sha256(genome),
                "gff3": str(gff3.relative_to(root)),
                "gff3_sha256": _sha256(gff3),
                "protein": str(protein.relative_to(root)),
                "protein_sha256": _sha256(protein),
                "cds": str(cds.relative_to(root)),
                "cds_sha256": _sha256(cds),
                "source_uri": f"test://{species_id}",
                "outgroup_species_id": "",
                "notes": "fixture",
            }
        )

    headers = list(rows[0])
    with (root / "species.tsv").open("w", encoding="utf-8") as handle:
        handle.write("\t".join(headers) + "\n")
        for row in rows:
            handle.write("\t".join(row[header] for header in headers) + "\n")

    (root / "manual_review/manual_truth_set.tsv").write_text(
        "species_id\tgene_id\texpected_status\tevidence\treviewer\treview_date\tnotes\n"
        "Rice01\tGene1\tPOSITIVE\tcurated domain\tReviewerA\t2026-08-09\t\n"
        "Rice02\tOther1\tNEGATIVE\tcurated non-member\tReviewerA\t2026-08-09\t\n",
        encoding="utf-8",
    )
    manifest = root / "benchmark.yaml"
    manifest.write_text(
        """schema_version: \"1.0\"
project:
  name: test_rice_pilot
  analysis_scope: target_pan_gene_family
  seed: 20260807
family:
  name: TEST_FAMILY
  approval_state: approved
  pfam_ids: [PF00001]
  interpro_ids: []
  hmm: references/family.hmm
  reference_proteins: references/family_reference.pep.fa
  manual_truth_set: manual_review/manual_truth_set.tsv
panel:
  species_table: species.tsv
  min_genomes: 2
  max_genomes: 2
acceptance:
  approval_state: approved
  require_family_approved: true
  require_sha256: true
  require_manual_truth_set: true
  min_positive_controls: 1
  min_negative_controls: 1
  require_representative: true
  require_annotation_version: true
  chromosome_level_policy: block
""",
        encoding="utf-8",
    )
    return manifest


def test_ready_benchmark_audit_and_outputs(tmp_path: Path) -> None:
    manifest = _write_ready_benchmark(tmp_path)
    audit = audit_benchmark(manifest)
    assert audit.overall_status == "READY"
    assert audit.blocking_failures == 0
    assert audit.warnings == 0
    outputs = write_benchmark_audit(audit, tmp_path / "audit")
    for path in outputs.values():
        assert path.is_file()
        assert path.stat().st_size > 0
    html = outputs["html"].read_text(encoding="utf-8")
    assert "生物学验收门" in html
    assert "READY" in html


def test_missing_required_input_blocks_benchmark(tmp_path: Path) -> None:
    manifest = _write_ready_benchmark(tmp_path)
    (manifest.parent / "inputs/Rice02/protein.fa").unlink()
    audit = audit_benchmark(manifest)
    assert audit.overall_status == "BLOCKED"
    assert any(
        check.check_id == "BMG110"
        and check.item_id == "Rice02"
        and check.scope == "species.protein"
        and check.status == "FAIL"
        for check in audit.checks
    )


def test_reference_aligned_sample_cannot_substitute_for_genome(tmp_path: Path) -> None:
    manifest = _write_ready_benchmark(tmp_path)
    species = manifest.parent / "species.tsv"
    text = species.read_text(encoding="utf-8")
    species.write_text(
        text.replace("assembled_genome", "reference_aligned_sample", 1), encoding="utf-8"
    )
    audit = audit_benchmark(manifest)
    assert audit.overall_status == "BLOCKED"
    assert any(check.check_id == "BMG100" and check.status == "FAIL" for check in audit.checks)


def test_checksum_mismatch_blocks_benchmark(tmp_path: Path) -> None:
    manifest = _write_ready_benchmark(tmp_path)
    genome = manifest.parent / "inputs/Rice01/genome.fa"
    genome.write_text(">Chr1\nAAAAAAAAA\n", encoding="utf-8")
    audit = audit_benchmark(manifest)
    assert audit.overall_status == "BLOCKED"
    assert any(check.check_id == "BMG112" and check.status == "FAIL" for check in audit.checks)


def test_benchmark_init_is_non_destructive(tmp_path: Path) -> None:
    target = tmp_path / "pilot"
    written = initialize_benchmark(target)
    assert written
    assert (target / "benchmark.yaml").is_file()
    assert (target / "README.zh-CN.md").is_file()
    with pytest.raises(FileExistsError):
        initialize_benchmark(target)


def test_benchmark_cli_commands_are_registered(tmp_path: Path) -> None:
    target = tmp_path / "pilot"
    result = runner.invoke(app, ["benchmark", "init", str(target)])
    assert result.exit_code == 0
    audit_result = runner.invoke(
        app,
        [
            "benchmark",
            "audit",
            "--manifest",
            str(target / "benchmark.yaml"),
            "--output",
            str(target / "audits/intake_001"),
            "--allow-blocked",
        ],
    )
    assert audit_result.exit_code == 0
    assert "BLOCKED" in audit_result.stdout
    assert (target / "audits/intake_001/benchmark_readiness.html").is_file()
