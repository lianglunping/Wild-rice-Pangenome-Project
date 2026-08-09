from __future__ import annotations

import sys
from pathlib import Path as _ScriptPath

sys.path.insert(0, str(_ScriptPath(snakemake.scriptdir)))

import html
import importlib.metadata
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pandas as pd
from workflow_utils import save_table, save_workbook, sha256_file, write_json, write_text_atomic

results_dir = Path(snakemake.params.results_dir)
results_dir.mkdir(parents=True, exist_ok=True)
family_path = results_dir / "02_family" / "family_members.tsv"
if family_path.is_file():
    master = pd.read_csv(family_path, sep="\t")
else:
    master = pd.DataFrame(columns=["stable_id", "species_id", "gene_id"])


def merge_one(path: Path, columns: list[str] | None = None, suffix: str = "") -> None:
    global master
    if not path.is_file() or master.empty:
        return
    table = pd.read_csv(path, sep="\t")
    if "stable_id" not in table.columns:
        return
    if columns is not None:
        table = table[[column for column in columns if column in table.columns]]
    table = table.drop_duplicates("stable_id")
    overlapping = set(master.columns).intersection(table.columns).difference({"stable_id"})
    if overlapping:
        table = table.rename(columns={column: f"{column}{suffix}" for column in overlapping})
    master = master.merge(table, on="stable_id", how="left", validate="one_to_one")


merge_one(results_dir / "04_gene_structure" / "gene_structure_metrics.tsv", suffix="_structure")
merge_one(results_dir / "07_chromosome" / "chromosome_distribution.tsv", suffix="_chromosome")
merge_one(results_dir / "08_duplication" / "duplication_mode.tsv", suffix="_duplication")

membership_path = results_dir / "06_pan_family" / "family_hog_membership.tsv"
classification_path = results_dir / "06_pan_family" / "pan_family_classification.tsv"
if membership_path.is_file() and classification_path.is_file() and not master.empty:
    membership = pd.read_csv(membership_path, sep="\t")
    classification = pd.read_csv(classification_path, sep="\t")
    membership = membership.merge(classification, on="HOG_ID", how="left", validate="many_to_one")
    membership = membership.drop_duplicates("stable_id")
    master = master.merge(membership, on="stable_id", how="left", suffixes=("", "_pan_family"))

promoter_path = results_dir / "10_promoter" / "promoter_elements_per_gene.tsv"
if promoter_path.is_file() and not master.empty:
    promoter = pd.read_csv(promoter_path, sep="\t")
    promoter_total = (
        promoter.groupby("stable_id", as_index=False)
        .agg(
            cis_element_total=("element_count", "sum"),
            cis_element_major_classes=(
                "major_class",
                lambda values: ";".join(sorted(set(values.dropna().astype(str)))),
            ),
        )
        .drop_duplicates("stable_id")
    )
    master = master.merge(promoter_total, on="stable_id", how="left", validate="one_to_one")

expression_path = results_dir / "11_expression" / "expression_summary.tsv"
merge_one(expression_path, suffix="_expression")

kaks_path = results_dir / "09_kaks" / "kaks_pairs.tsv"
if kaks_path.is_file() and not master.empty:
    kaks = pd.read_csv(kaks_path, sep="\t")
    required = {"stable_id_1", "stable_id_2", "Ka", "Ks", "Ka_Ks"}
    if required.issubset(kaks.columns):
        first = kaks[["stable_id_1", "Ka", "Ks", "Ka_Ks"]].rename(
            columns={"stable_id_1": "stable_id"}
        )
        second = kaks[["stable_id_2", "Ka", "Ks", "Ka_Ks"]].rename(
            columns={"stable_id_2": "stable_id"}
        )
        per_gene_kaks = (
            pd.concat([first, second], ignore_index=True)
            .groupby("stable_id", as_index=False)
            .agg(
                kaks_pair_count=("Ka_Ks", "count"),
                median_Ka=("Ka", "median"),
                median_Ks=("Ks", "median"),
                median_Ka_Ks=("Ka_Ks", "median"),
            )
        )
        master = master.merge(per_gene_kaks, on="stable_id", how="left", validate="one_to_one")

master_tsv = Path(snakemake.output.master_tsv)
master_xlsx = Path(snakemake.output.master_xlsx)
save_table(master, master_tsv)
save_workbook({"master_gene_table": master}, master_xlsx)

python_packages = []
for distribution in ("panfamflow", "pandas", "pydantic", "pyyaml", "openpyxl"):
    try:
        version = importlib.metadata.version(distribution)
    except importlib.metadata.PackageNotFoundError:
        version = "not-installed-in-report-environment"
    python_packages.append({"software": distribution, "version": version, "source": "Python"})
software = pd.DataFrame(python_packages)
save_table(software, snakemake.output.software_versions)

excluded_manifest_paths = {
    Path(snakemake.output.index),
    Path(snakemake.output.manifest),
    Path(snakemake.output.run_info),
}
existing_before_run_info = [
    path
    for path in sorted(results_dir.rglob("*"))
    if path.is_file() and path not in excluded_manifest_paths
]
run_info = {
    "project": snakemake.config.get("project", {}),
    "selected_modules": str(snakemake.params.selected_modules).split(","),
    "generated_at_utc": datetime.now(UTC).isoformat(),
    "config_path": snakemake.config.get("panfamflow_config_path"),
    "master_rows": int(master.shape[0]),
    # run_info.json itself is included in the manifest created immediately below.
    "manifest_files": len(existing_before_run_info) + 1,
}
write_json(run_info, snakemake.output.run_info)

manifest_rows: list[dict[str, Any]] = []
for path in sorted(results_dir.rglob("*")):
    if not path.is_file() or path in {
        Path(snakemake.output.index),
        Path(snakemake.output.manifest),
    }:
        continue
    manifest_rows.append(
        {
            "relative_path": str(path.relative_to(results_dir)),
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
            "suffix": path.suffix.lower(),
        }
    )
manifest = pd.DataFrame(manifest_rows)
if int(manifest.shape[0]) != int(run_info["manifest_files"]):
    raise RuntimeError("Result manifest file count changed unexpectedly during report generation.")
save_table(manifest, snakemake.output.manifest)

plot_files = [path for path in sorted(results_dir.rglob("*.png")) if "report" not in path.parts]
table_files = [
    path for path in sorted(results_dir.rglob("*.tsv")) if path != Path(snakemake.output.manifest)
]
report_title = snakemake.params.title or snakemake.config.get("project", {}).get(
    "name", "PanFamFlow pan-gene-family report"
)
report_path = Path(snakemake.output.index)
report_path.parent.mkdir(parents=True, exist_ok=True)

cards = [
    ("Family genes", str(master.shape[0])),
    ("Result files", str(manifest.shape[0])),
    ("Selected modules", html.escape(str(snakemake.params.selected_modules))),
]
card_html = "".join(
    f'<div class="card"><div class="value">{value}</div><div>{label}</div></div>'
    for label, value in cards
)
plot_html = "".join(
    f'<figure><img src="../{html.escape(str(path.relative_to(results_dir)))}" '
    f'alt="{html.escape(path.stem)}"><figcaption>{html.escape(path.stem)}</figcaption></figure>'
    for path in plot_files
)
table_html = "".join(
    f'<li><a href="../{html.escape(str(path.relative_to(results_dir)))}">'
    f"{html.escape(str(path.relative_to(results_dir)))}</a></li>"
    for path in table_files
)
write_text_atomic(
    f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(str(report_title))}</title>
<style>
body {{ font-family: Arial, sans-serif; margin: 2rem auto; max-width: 1200px; color: #222; line-height: 1.5; }}
h1, h2 {{ border-bottom: 1px solid #ddd; padding-bottom: .35rem; }}
.cards {{ display: flex; flex-wrap: wrap; gap: 1rem; }}
.card {{ border: 1px solid #ddd; padding: 1rem; min-width: 180px; border-radius: 4px; }}
.value {{ font-size: 1.6rem; font-weight: bold; }}
.gallery {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(340px, 1fr)); gap: 1rem; }}
figure {{ margin: 0; border: 1px solid #ddd; padding: .7rem; }}
img {{ max-width: 100%; height: auto; }}
code {{ background: #f4f4f4; padding: .1rem .25rem; }}
</style>
</head>
<body>
<h1>{html.escape(str(report_title))}</h1>
<p>Generated by PanFamFlow v0.1.1 (target pan-gene-family analysis) at {html.escape(run_info["generated_at_utc"])}. This report links computational outputs; biological conclusions still require project-specific QC and interpretation.</p>
<div class="cards">{card_html}</div>
<h2>Integrated master table</h2>
<p><a href="../{html.escape(str(master_tsv.relative_to(results_dir)))}">TSV</a> · <a href="../{html.escape(str(master_xlsx.relative_to(results_dir)))}">XLSX</a></p>
<h2>Figures</h2>
<div class="gallery">{plot_html or "<p>No PNG figures are currently available.</p>"}</div>
<h2>Tables</h2>
<ul>{table_html or "<li>No TSV tables are currently available.</li>"}</ul>
<h2>Provenance</h2>
<p><a href="result_manifest.tsv">SHA256 result manifest</a> · <a href="software_versions.tsv">Software versions</a> · <a href="run_info.json">Run information</a></p>
</body>
</html>
""",
    report_path,
)
