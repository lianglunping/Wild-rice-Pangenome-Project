import sys
from pathlib import Path as _ScriptPath

sys.path.insert(0, str(_ScriptPath(snakemake.scriptdir)))

from collections import defaultdict
from typing import Any

import pandas as pd
from workflow_utils import iter_gff, read_fasta, reverse_complement, save_table, write_fasta

members = pd.read_csv(snakemake.input.members, sep="\t")
family_ids = set(members["stable_id"].astype(str))
upstream = int(snakemake.params.upstream_bp)
downstream = int(snakemake.params.downstream_bp)
separator = str(snakemake.params.separator)
rows: list[dict[str, Any]] = []
promoters: dict[str, str] = {}

for species, genome_path, gff_path, mapping_path in zip(
    snakemake.params.species_ids,
    snakemake.input.genomes,
    snakemake.input.gff3s,
    snakemake.input.maps,
    strict=True,
):
    species = str(species)
    genome = read_fasta(genome_path)
    mapping = pd.read_csv(mapping_path, sep="\t")
    mapping = mapping.loc[mapping["stable_id"].astype(str).isin(family_ids)]
    if mapping.empty:
        continue
    gene_intervals: defaultdict[str, list[tuple[int, int, str]]] = defaultdict(list)
    for feature in iter_gff(gff_path):
        if str(feature["feature"]).lower() != "gene":
            continue
        attrs = feature["attributes"]
        gene_id = attrs.get("ID") or attrs.get("gene_id") or attrs.get("Name") or ""
        gene_intervals[str(feature["seqid"])].append(
            (int(feature["start"]), int(feature["end"]), gene_id)
        )
    for record in mapping.to_dict(orient="records"):
        stable_id = str(record["stable_id"])
        chromosome = str(record["chromosome"])
        if chromosome not in genome:
            rows.append(
                {
                    "species_id": species,
                    "gene_id": record["gene_id"],
                    "stable_id": stable_id,
                    "chromosome": chromosome,
                    "promoter_qc": "CHROMOSOME_NOT_FOUND",
                }
            )
            continue
        gene_start = int(record["gene_start"])
        gene_end = int(record["gene_end"])
        strand = str(record["strand"])
        chromosome_length = len(genome[chromosome])
        if strand == "+":
            tss0 = gene_start - 1
            raw_start0 = tss0 - upstream
            raw_end0 = tss0 + downstream
        elif strand == "-":
            tss0 = gene_end
            raw_start0 = tss0 - downstream
            raw_end0 = tss0 + upstream
        else:
            rows.append(
                {
                    "species_id": species,
                    "gene_id": record["gene_id"],
                    "stable_id": stable_id,
                    "chromosome": chromosome,
                    "promoter_qc": "INVALID_STRAND",
                }
            )
            continue
        start0 = max(0, raw_start0)
        end0 = min(chromosome_length, raw_end0)
        sequence = genome[chromosome][start0:end0]
        if strand == "-":
            sequence = reverse_complement(sequence)
        overlap_genes = [
            gene_id
            for other_start, other_end, gene_id in gene_intervals[chromosome]
            if gene_id != record["gene_id"] and not (other_end < start0 + 1 or other_start > end0)
        ]
        truncated = start0 != raw_start0 or end0 != raw_end0
        status = "TRUNCATED_AT_CHROMOSOME_BOUNDARY" if truncated else "PASS"
        promoters[stable_id] = sequence
        rows.append(
            {
                "species_id": species,
                "gene_id": record["gene_id"],
                "stable_id": stable_id,
                "chromosome": chromosome,
                "strand": strand,
                "tss_1based": gene_start if strand == "+" else gene_end,
                "promoter_start_1based": start0 + 1,
                "promoter_end_1based": end0,
                "promoter_length": len(sequence),
                "requested_upstream_bp": upstream,
                "requested_downstream_bp": downstream,
                "truncated_at_chr_boundary": truncated,
                "overlap_gene_count": len(overlap_genes),
                "overlap_gene_ids": ";".join(sorted(overlap_genes)) or pd.NA,
                "promoter_qc": status,
            }
        )

coordinates = pd.DataFrame(rows)
if not promoters:
    raise RuntimeError("No promoter sequence could be extracted for family members.")
write_fasta(promoters, snakemake.output.fasta)
save_table(coordinates, snakemake.output.coordinates, snakemake.output.coordinates_xlsx)
