import sys
from pathlib import Path as _ScriptPath

sys.path.insert(0, str(_ScriptPath(snakemake.scriptdir)))

import itertools
from collections import defaultdict
from typing import Any

import pandas as pd
from workflow_utils import first_parent, iter_gff, save_table, save_workbook

members = pd.read_csv(snakemake.input.members, sep="\t")
member_ids = set(members["stable_id"].astype(str))
rows: list[dict[str, Any]] = []

for gff_path, mapping_path in zip(snakemake.input.gff3s, snakemake.input.maps, strict=True):
    mapping = pd.read_csv(mapping_path, sep="\t")
    mapping = mapping.loc[mapping["stable_id"].astype(str).isin(member_ids)].copy()
    if mapping.empty:
        continue
    species = str(mapping["species_id"].iloc[0])
    transcript_to_stable = dict(
        zip(mapping["transcript_id"].astype(str), mapping["stable_id"].astype(str), strict=True)
    )
    gene_to_stable = dict(
        zip(mapping["gene_id"].astype(str), mapping["stable_id"].astype(str), strict=True)
    )
    feature_by_stable: defaultdict[str, defaultdict[str, list[tuple[int, int]]]] = defaultdict(
        lambda: defaultdict(list)
    )
    gene_coords: dict[str, tuple[str, int, int, str]] = {}
    transcript_coords: dict[str, tuple[int, int]] = {}

    for feature in iter_gff(gff_path):
        attrs = feature["attributes"]
        feature_type = str(feature["feature"]).lower()
        if feature_type == "gene":
            gene_id = attrs.get("ID") or attrs.get("gene_id") or attrs.get("Name")
            if gene_id in gene_to_stable:
                stable_id = gene_to_stable[gene_id]
                gene_coords[stable_id] = (
                    str(feature["seqid"]),
                    int(feature["start"]),
                    int(feature["end"]),
                    str(feature["strand"]),
                )
            continue
        if feature_type in {"mrna", "transcript", "ncrna", "trna", "rrna"}:
            transcript_id = attrs.get("ID") or attrs.get("transcript_id")
            if transcript_id in transcript_to_stable:
                transcript_coords[transcript_to_stable[transcript_id]] = (
                    int(feature["start"]),
                    int(feature["end"]),
                )
            continue
        parents = (attrs.get("Parent") or attrs.get("transcript_id") or "").split(",")
        for parent in parents:
            parent = first_parent(parent)
            if parent in transcript_to_stable:
                stable_id = transcript_to_stable[parent]
                feature_by_stable[stable_id][feature_type].append(
                    (int(feature["start"]), int(feature["end"]))
                )

    for record in mapping.to_dict(orient="records"):
        stable_id = str(record["stable_id"])
        chromosome, gene_start, gene_end, strand = gene_coords.get(
            stable_id,
            (
                record.get("chromosome", pd.NA),
                record.get("gene_start", pd.NA),
                record.get("gene_end", pd.NA),
                record.get("strand", pd.NA),
            ),
        )
        transcript_start, transcript_end = transcript_coords.get(
            stable_id,
            (record.get("transcript_start", pd.NA), record.get("transcript_end", pd.NA)),
        )
        features = feature_by_stable[stable_id]
        exons = sorted(features.get("exon", []))
        cds_parts = sorted(features.get("cds", []))
        five_utrs = features.get("five_prime_utr", []) + features.get("5utr", [])
        three_utrs = features.get("three_prime_utr", []) + features.get("3utr", [])
        generic_utrs = features.get("utr", [])
        intron_lengths = [
            max(0, right[0] - left[1] - 1) for left, right in itertools.pairwise(exons)
        ]

        def total_length(intervals: list[tuple[int, int]]) -> int:
            return sum(end - start + 1 for start, end in intervals)

        gene_length = (
            int(gene_end) - int(gene_start) + 1
            if pd.notna(gene_start) and pd.notna(gene_end)
            else pd.NA
        )
        transcript_length = (
            int(transcript_end) - int(transcript_start) + 1
            if pd.notna(transcript_start) and pd.notna(transcript_end)
            else pd.NA
        )
        rows.append(
            {
                "species_id": species,
                "gene_id": record["gene_id"],
                "transcript_id": record["transcript_id"],
                "stable_id": stable_id,
                "subfamily": members.set_index("stable_id").loc[stable_id].get("subfamily", pd.NA),
                "group": record.get("group", pd.NA),
                "chromosome": chromosome,
                "strand": strand,
                "gene_length": gene_length,
                "transcript_span": transcript_length,
                "protein_length": record.get("protein_length", pd.NA),
                "cds_length": total_length(cds_parts),
                "exon_count": len(exons),
                "total_exon_length": total_length(exons),
                "mean_exon_length": total_length(exons) / len(exons) if exons else pd.NA,
                "intron_count": len(intron_lengths),
                "total_intron_length": sum(intron_lengths),
                "mean_intron_length": (
                    sum(intron_lengths) / len(intron_lengths) if intron_lengths else 0.0
                ),
                "five_prime_utr_length": total_length(five_utrs),
                "three_prime_utr_length": total_length(three_utrs),
                "generic_utr_length": total_length(generic_utrs),
                "structure_qc": "PASS" if exons and cds_parts else "MISSING_EXON_OR_CDS",
            }
        )

metrics = pd.DataFrame(rows)
if metrics.empty:
    raise RuntimeError("No gene-structure records were produced for family members.")
summary_columns = [
    "gene_length",
    "protein_length",
    "cds_length",
    "exon_count",
    "intron_count",
    "total_intron_length",
]
summary_rows: list[dict[str, object]] = []
for group_field in ("subfamily", "group"):
    for group_value, subset in metrics.groupby(group_field, dropna=False):
        for metric in summary_columns:
            values = pd.to_numeric(subset[metric], errors="coerce").dropna()
            if values.empty:
                continue
            summary_rows.append(
                {
                    "group_field": group_field,
                    "group_value": group_value,
                    "metric": metric,
                    "n": int(values.shape[0]),
                    "mean": float(values.mean()),
                    "sd": float(values.std(ddof=1)) if values.shape[0] > 1 else pd.NA,
                    "median": float(values.median()),
                    "q1": float(values.quantile(0.25)),
                    "q3": float(values.quantile(0.75)),
                    "min": float(values.min()),
                    "max": float(values.max()),
                }
            )
summary = pd.DataFrame(summary_rows)
save_table(metrics, snakemake.output.metrics)
save_table(summary, snakemake.output.summary)
save_workbook({"metrics": metrics, "summary": summary}, snakemake.output.xlsx)
