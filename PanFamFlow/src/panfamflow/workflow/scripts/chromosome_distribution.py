import sys
from pathlib import Path as _ScriptPath

sys.path.insert(0, str(_ScriptPath(snakemake.scriptdir)))

from typing import Any

import matplotlib.pyplot as plt
import pandas as pd
from workflow_utils import fasta_lengths, save_table, save_workbook

members = pd.read_csv(snakemake.input.members, sep="\t")
representative_only = bool(snakemake.params.representative_only)
representatives = set(str(item) for item in snakemake.params.representatives)
if representative_only:
    members = members.loc[members["species_id"].astype(str).isin(representatives)].copy()

lengths_by_species: dict[str, dict[str, int]] = {}
for species, genome_path in zip(snakemake.params.species_ids, snakemake.input.genomes, strict=True):
    lengths_by_species[str(species)] = dict(fasta_lengths(genome_path))

rows: list[dict[str, Any]] = []
for row in members.to_dict(orient="records"):
    species = str(row["species_id"])
    chromosome = str(row.get("chromosome", ""))
    start = pd.to_numeric(row.get("gene_start"), errors="coerce")
    end = pd.to_numeric(row.get("gene_end"), errors="coerce")
    chromosome_length = lengths_by_species.get(species, {}).get(chromosome)
    if pd.isna(start) or pd.isna(end) or chromosome_length is None:
        status = "MISSING_COORDINATE_OR_CHROMOSOME"
        midpoint = pd.NA
        relative_position = pd.NA
    else:
        midpoint = (float(start) + float(end)) / 2.0
        relative_position = midpoint / chromosome_length
        status = "PASS"
    rows.append(
        {
            "species_id": species,
            "gene_id": row["gene_id"],
            "stable_id": row["stable_id"],
            "subfamily": row.get("subfamily", pd.NA),
            "group": row.get("group", pd.NA),
            "chromosome": chromosome,
            "start": start,
            "end": end,
            "midpoint": midpoint,
            "chromosome_length": chromosome_length if chromosome_length is not None else pd.NA,
            "relative_position": relative_position,
            "coordinate_qc": status,
        }
    )
distribution = pd.DataFrame(rows)
summary = (
    distribution.loc[distribution["coordinate_qc"] == "PASS"]
    .groupby(["species_id", "chromosome", "chromosome_length"], as_index=False)
    .agg(gene_count=("stable_id", "nunique"))
)
summary["gene_density_per_mb"] = summary["gene_count"] / (summary["chromosome_length"] / 1_000_000)
save_table(distribution, snakemake.output.distribution)
save_table(summary, snakemake.output.summary)
save_workbook(
    {"gene_coordinates": distribution, "chromosome_summary": summary}, snakemake.output.xlsx
)

valid = distribution.loc[distribution["coordinate_qc"] == "PASS"].copy()
if valid.empty:
    raise RuntimeError("No family member has a valid chromosome coordinate.")
species_order = list(dict.fromkeys(valid["species_id"].astype(str)))
fig_height = max(3.5, 2.2 * len(species_order))
fig, axes = plt.subplots(len(species_order), 1, figsize=(10, fig_height), squeeze=False)
for axis, species in zip(axes[:, 0], species_order, strict=True):
    subset = valid.loc[valid["species_id"].astype(str) == species].copy()
    chromosome_order = sorted(
        subset["chromosome"].unique(),
        key=lambda value: (len(str(value)), str(value)),
    )
    y_lookup = {chromosome: index for index, chromosome in enumerate(chromosome_order)}
    for chromosome in chromosome_order:
        chromosome_length = float(
            subset.loc[subset["chromosome"] == chromosome, "chromosome_length"].iloc[0]
        )
        axis.hlines(y_lookup[chromosome], 0, chromosome_length / 1_000_000, linewidth=5, alpha=0.25)
    axis.scatter(
        subset["midpoint"].astype(float) / 1_000_000,
        subset["chromosome"].map(y_lookup),
        s=18,
        alpha=0.85,
    )
    axis.set_yticks(range(len(chromosome_order)), chromosome_order)
    axis.set_title(species, loc="left", fontweight="bold")
    axis.set_xlabel("Position (Mb)")
    axis.set_ylabel("Chromosome")
    axis.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(snakemake.output.plot_pdf)
fig.savefig(snakemake.output.plot_png, dpi=int(snakemake.params.png_dpi))
plt.close(fig)
