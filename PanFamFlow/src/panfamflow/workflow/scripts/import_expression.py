import sys
from pathlib import Path as _ScriptPath

sys.path.insert(0, str(_ScriptPath(snakemake.scriptdir)))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from workflow_utils import resolve_column, save_table, save_workbook

members = pd.read_csv(snakemake.input.members, sep="\t")
source = pd.read_csv(snakemake.input.matrix, sep=None, engine="python")
separator = str(snakemake.params.separator)

stable_column = resolve_column(source, ["stable_id", "protein_id"], required=False)
if stable_column is None:
    species_column = resolve_column(source, ["species_id", "species"], required=False)
    gene_column = resolve_column(source, ["gene_id", "gene"], required=False)
    if species_column and gene_column:
        source["stable_id"] = [
            f"{species}{separator}{gene}"
            for species, gene in zip(
                source[species_column].astype(str), source[gene_column].astype(str), strict=True
            )
        ]
    else:
        first_column = source.columns[0]
        raw_ids = source[first_column].astype(str)
        gene_to_stable = members.groupby("gene_id")["stable_id"].agg(list).to_dict()
        ambiguous = [gene for gene in raw_ids if len(gene_to_stable.get(gene, [])) != 1]
        if ambiguous:
            raise ValueError(
                "Expression row IDs are not stable IDs and are not uniquely mappable gene IDs. "
                f"Examples: {ambiguous[:10]}"
            )
        source["stable_id"] = [gene_to_stable[gene][0] for gene in raw_ids]
else:
    source = source.rename(columns={stable_column: "stable_id"})

metadata_columns = {
    "stable_id",
    "species_id",
    "species",
    "gene_id",
    "gene",
    "transcript_id",
    "protein_id",
}
sample_columns = [
    column for column in source.columns if str(column).lower() not in metadata_columns
]
if not sample_columns:
    raise ValueError("Expression matrix contains no sample columns.")
for column in sample_columns:
    source[column] = pd.to_numeric(source[column], errors="raise")
if source["stable_id"].duplicated().any():
    duplicates = source.loc[source["stable_id"].duplicated(), "stable_id"].astype(str).tolist()
    raise ValueError(f"Expression matrix has duplicate stable IDs: {duplicates[:10]}")

family_ids = set(members["stable_id"].astype(str))
matrix = source.loc[
    source["stable_id"].astype(str).isin(family_ids), ["stable_id", *sample_columns]
]
matrix = members[["stable_id", "species_id", "gene_id", "subfamily"]].merge(
    matrix, on="stable_id", how="left", validate="one_to_one"
)
missing_rows = matrix[sample_columns].isna().all(axis=1)
matrix["expression_data_status"] = np.where(missing_rows, "MISSING", "AVAILABLE")
long = matrix.melt(
    id_vars=["stable_id", "species_id", "gene_id", "subfamily", "expression_data_status"],
    value_vars=sample_columns,
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
save_table(matrix, snakemake.output.matrix)
save_table(long, snakemake.output.long)
save_table(summary, snakemake.output.summary)
save_workbook({"matrix": matrix, "long": long, "summary": summary}, snakemake.output.xlsx)

values = matrix[sample_columns].astype(float).to_numpy()
transformed = np.log2(values + 1.0)
if str(snakemake.params.heatmap_transform) == "log2_tpm1_zscore":
    means = np.nanmean(transformed, axis=1, keepdims=True)
    standard_deviations = np.nanstd(transformed, axis=1, keepdims=True)
    standard_deviations[standard_deviations == 0] = 1.0
    transformed = (transformed - means) / standard_deviations
fig_height = max(4.8, min(18.0, 0.18 * max(1, matrix.shape[0])))
fig, axis = plt.subplots(figsize=(max(7.0, 0.35 * len(sample_columns)), fig_height))
image = axis.imshow(transformed, aspect="auto", interpolation="nearest", cmap="viridis")
axis.set_xticks(range(len(sample_columns)), sample_columns, rotation=45, ha="right")
if matrix.shape[0] <= 80:
    axis.set_yticks(range(matrix.shape[0]), matrix["stable_id"].astype(str), fontsize=6)
else:
    axis.set_yticks([])
axis.set_xlabel("Sample")
axis.set_ylabel("Family gene")
fig.colorbar(image, ax=axis, label=str(snakemake.params.heatmap_transform))
fig.tight_layout()
fig.savefig(snakemake.output.plot_pdf)
fig.savefig(snakemake.output.plot_png, dpi=int(snakemake.params.png_dpi))
plt.close(fig)
