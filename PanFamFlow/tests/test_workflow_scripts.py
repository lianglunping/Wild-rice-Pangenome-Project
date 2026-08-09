from __future__ import annotations

import runpy
from pathlib import Path
from types import SimpleNamespace

import pandas as pd

SCRIPT_DIR = Path(__file__).parents[1] / "src" / "panfamflow" / "workflow" / "scripts"
TOY_DIR = Path(__file__).parents[1] / "examples" / "toy"


def test_input_audit_script(tmp_path: Path) -> None:
    output = SimpleNamespace(
        tsv=str(tmp_path / "input_audit.tsv"),
        xlsx=str(tmp_path / "input_audit.xlsx"),
        manifest=str(tmp_path / "manifest.json"),
    )
    fake = SimpleNamespace(
        scriptdir=str(SCRIPT_DIR),
        params=SimpleNamespace(
            records=[
                {
                    "species_id": "SpA",
                    "role": "genome",
                    "path": str(TOY_DIR / "data" / "SpA" / "genome.fa"),
                },
                {
                    "species_id": "SpA",
                    "role": "gff3",
                    "path": str(TOY_DIR / "data" / "SpA" / "annotation.gff3"),
                },
            ],
            calculate_sha256=True,
        ),
        config={"project": {"name": "toy"}, "panfamflow_selected_modules": "qc"},
        output=output,
    )
    runpy.run_path(str(SCRIPT_DIR / "input_audit.py"), init_globals={"snakemake": fake})
    table = pd.read_csv(output.tsv, sep="\t")
    assert table.shape[0] == 2
    assert set(table["status"]) == {"PASS"}
    assert table.loc[table["role"] == "genome", "record_count"].iloc[0] == 1


def test_pan_family_parser(tmp_path: Path) -> None:
    result_dir = tmp_path / "orthofinder"
    hog_dir = result_dir / "Phylogenetic_Hierarchical_Orthogroups"
    hog_dir.mkdir(parents=True)
    (hog_dir / "N0.tsv").write_text(
        "HOG\tOG\tGene Tree Parent Clade\tSpA\tSpB\n"
        "N0.HOG0000001\tOG0000001\tN0\tSpA__GeneA1\tSpB__GeneB1\n"
        "N0.HOG0000002\tOG0000002\tN0\tSpA__Other\t\n",
        encoding="utf-8",
    )
    marker = tmp_path / "result_dir.txt"
    marker.write_text(str(result_dir), encoding="utf-8")
    members = tmp_path / "family_members.tsv"
    members.write_text(
        "species_id\tgene_id\tstable_id\nSpA\tGeneA1\tSpA__GeneA1\nSpB\tGeneB1\tSpB__GeneB1\n",
        encoding="utf-8",
    )
    output = SimpleNamespace(
        classification=str(tmp_path / "classification.tsv"),
        membership=str(tmp_path / "membership.tsv"),
        presence=str(tmp_path / "presence.tsv"),
        unassigned_members=str(tmp_path / "unassigned.tsv"),
        rarefaction=str(tmp_path / "rarefaction.tsv"),
        rarefaction_summary=str(tmp_path / "rarefaction_summary.tsv"),
        xlsx=str(tmp_path / "pan_family.xlsx"),
        class_plot_pdf=str(tmp_path / "classes.pdf"),
        class_plot_png=str(tmp_path / "classes.png"),
        rarefaction_plot_pdf=str(tmp_path / "rarefaction.pdf"),
        rarefaction_plot_png=str(tmp_path / "rarefaction.png"),
    )
    fake = SimpleNamespace(
        scriptdir=str(SCRIPT_DIR),
        input=SimpleNamespace(result_dir=str(marker), members=str(members)),
        output=output,
        params=SimpleNamespace(
            species_ids=["SpA", "SpB"],
            separator="__",
            hog_node="N0",
            core_min=0.99,
            soft_core_min=0.90,
            shell_min=0.10,
            rarefaction_iterations=10,
            max_exact_combinations=100,
            seed=20260807,
            png_dpi=120,
        ),
    )
    runpy.run_path(str(SCRIPT_DIR / "parse_pan_family.py"), init_globals={"snakemake": fake})
    classification = pd.read_csv(output.classification, sep="\t")
    assert classification.shape[0] == 1
    assert classification.loc[0, "pan_family_class"] == "Core"
    assert classification.loc[0, "analysis_scope"] == "TARGET_GENE_FAMILY_ONLY"
    assert classification.loc[0, "HOG_ID"] == "N0.HOG0000001"
    unassigned = pd.read_csv(output.unassigned_members, sep="\t")
    assert unassigned.empty
    assert Path(output.class_plot_pdf).is_file()


def test_stringtie_combiner_uses_species_scoped_gene_ids(tmp_path: Path) -> None:
    members = tmp_path / "family_members.tsv"
    members.write_text(
        "stable_id\tspecies_id\tgene_id\tsubfamily\n"
        "SpA__Gene1\tSpA\tGene1\tA\n"
        "SpB__Gene1\tSpB\tGene1\tB\n",
        encoding="utf-8",
    )
    map_a = tmp_path / "SpA.map.tsv"
    map_b = tmp_path / "SpB.map.tsv"
    map_a.write_text(
        "species_id\tgene_id\tstable_id\nSpA\tGene1\tSpA__Gene1\n",
        encoding="utf-8",
    )
    map_b.write_text(
        "species_id\tgene_id\tstable_id\nSpB\tGene1\tSpB__Gene1\n",
        encoding="utf-8",
    )
    abundance_a = tmp_path / "sample_a.tsv"
    abundance_b = tmp_path / "sample_b.tsv"
    abundance_a.write_text("Gene ID\tTPM\nGene1\t3.5\n", encoding="utf-8")
    abundance_b.write_text("Gene ID\tTPM\nGene1\t8.0\n", encoding="utf-8")
    output = SimpleNamespace(
        matrix=str(tmp_path / "expression_matrix.tsv"),
        long=str(tmp_path / "expression_long.tsv"),
        summary=str(tmp_path / "expression_summary.tsv"),
        xlsx=str(tmp_path / "expression.xlsx"),
        plot_pdf=str(tmp_path / "expression.pdf"),
        plot_png=str(tmp_path / "expression.png"),
    )
    fake = SimpleNamespace(
        scriptdir=str(SCRIPT_DIR),
        input=SimpleNamespace(
            members=str(members),
            maps=[str(map_a), str(map_b)],
            abundance=[str(abundance_a), str(abundance_b)],
        ),
        output=output,
        params=SimpleNamespace(
            sample_ids=["SampleA", "SampleB"],
            sample_species_ids=["SpA", "SpB"],
            min_tpm_detected=1.0,
            heatmap_transform="log2_tpm1_zscore",
            png_dpi=120,
        ),
    )
    runpy.run_path(str(SCRIPT_DIR / "combine_stringtie.py"), init_globals={"snakemake": fake})
    matrix = pd.read_csv(output.matrix, sep="\t").set_index("stable_id")
    assert matrix.loc["SpA__Gene1", "SampleA"] == 3.5
    assert matrix.loc["SpA__Gene1", "SampleB"] == 0.0
    assert matrix.loc["SpB__Gene1", "SampleA"] == 0.0
    assert matrix.loc["SpB__Gene1", "SampleB"] == 8.0


def test_orthofinder_preserves_prefixed_stable_ids() -> None:
    source = (SCRIPT_DIR / "run_orthofinder.py").read_text(encoding="utf-8")
    assert '"-X"' in source
