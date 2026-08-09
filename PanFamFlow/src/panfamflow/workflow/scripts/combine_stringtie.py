import sys
from pathlib import Path as _ScriptPath

sys.path.insert(0, str(_ScriptPath(snakemake.scriptdir)))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from workflow_utils import resolve_column, save_table, save_workbook

members = pd.read_csv(snakemake.input.members, sep="\t")
family_ids = set(members["stable_id"].astype(str))
map_frames = [pd.read_csv(path, sep="\t") for path in snakemake.input.maps]
mapping = pd.concat(map_frames, ignore_index=True)
sample_ids = [str(item) for item in snakemake.params.sample_ids]
sample_species_ids = [str(item) for item in snakemake.params.sample_species_ids]
if not (len(sample_ids) == len(sample_species_ids) == len(snakemake.input.abundance)):
    raise ValueError(
        "StringTie abundance files, sample IDs and sample species IDs have different lengths"
    )

stable_lookup_by_species: dict[str, dict[str, str]] = {}
for species_id, species_mapping in mapping.groupby("species_id", sort=False):
    gene_ids = species_mapping["gene_id"].astype(str)
    duplicated = sorted(gene_ids[gene_ids.duplicated(keep=False)].unique())
    if duplicated:
        raise ValueError(
            f"Canonical mapping for {species_id} contains duplicate gene IDs: "
            + ", ".join(duplicated[:10])
        )
    stable_lookup_by_species[str(species_id)] = dict(
        zip(gene_ids, species_mapping["stable_id"].astype(str), strict=True)
    )

series: list[pd.Series] = []
for sample_id, species_id, path in zip(
    sample_ids, sample_species_ids, snakemake.input.abundance, strict=True
):
    if species_id not in stable_lookup_by_species:
        raise ValueError(f"No canonical ID mapping was found for sample species {species_id!r}")
    table = pd.read_csv(path, sep="\t")
    gene_column = resolve_column(table, ["Gene ID", "gene_id", "gene"])
    tpm_column = resolve_column(table, ["TPM", "tpm"])
    table["stable_id"] = table[gene_column].astype(str).map(stable_lookup_by_species[species_id])
    table = table.dropna(subset=["stable_id"])
    values = table.groupby("stable_id")[tpm_column].sum()
    values.name = sample_id
    series.append(values)
wide = pd.concat(series, axis=1).fillna(0.0)
wide = wide.loc[wide.index.astype(str).isin(family_ids)].reset_index()
wide = members[["stable_id", "species_id", "gene_id", "subfamily"]].merge(
    wide, on="stable_id", how="left", validate="one_to_one"
)
wide[sample_ids] = wide[sample_ids].fillna(0.0)
long = wide.melt(
    id_vars=["stable_id", "species_id", "gene_id", "subfamily"],
    value_vars=sample_ids,
    var_name="sample_id",
    value_name="expression_value",
)
long["detected"] = long["expression_value"] >= float(snakemake.params.min_tpm_detected)
summary = long.groupby(["stable_id", "species_id", "gene_id"], as_index=False).agg(
    samples_available=("expression_value", "count"),
    expression_detected_samples=("detected", "sum"),
    expression_detected_fraction=("detected", "mean"),
    median_expression=("expression_value", "median"),
    max_expression=("expression_value", "max"),
)
save_table(wide, snakemake.output.matrix)
save_table(long, snakemake.output.long)
save_table(summary, snakemake.output.summary)
save_workbook({"matrix": wide, "long": long, "summary": summary}, snakemake.output.xlsx)

values = wide[sample_ids].to_numpy(dtype=float)
transformed = np.log2(values + 1.0)
if str(snakemake.params.heatmap_transform) == "log2_tpm1_zscore":
    means = transformed.mean(axis=1, keepdims=True)
    sd = transformed.std(axis=1, keepdims=True)
    sd[sd == 0] = 1.0
    transformed = (transformed - means) / sd
fig_height = max(4.8, min(18.0, 0.18 * max(1, wide.shape[0])))
fig, axis = plt.subplots(figsize=(max(7.0, 0.35 * len(sample_ids)), fig_height))
image = axis.imshow(transformed, aspect="auto", interpolation="nearest", cmap="viridis")
axis.set_xticks(range(len(sample_ids)), sample_ids, rotation=45, ha="right")
if wide.shape[0] <= 80:
    axis.set_yticks(range(wide.shape[0]), wide["stable_id"].astype(str), fontsize=6)
else:
    axis.set_yticks([])
fig.colorbar(image, ax=axis, label=str(snakemake.params.heatmap_transform))
fig.tight_layout()
fig.savefig(snakemake.output.plot_pdf)
fig.savefig(snakemake.output.plot_png, dpi=int(snakemake.params.png_dpi))
plt.close(fig)
