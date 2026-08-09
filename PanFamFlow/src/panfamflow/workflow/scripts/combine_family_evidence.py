import sys
from pathlib import Path as _ScriptPath

sys.path.insert(0, str(_ScriptPath(snakemake.scriptdir)))

from collections import OrderedDict
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from workflow_utils import (
    iter_fasta_records,
    resolve_column,
    save_table,
    save_workbook,
    write_fasta,
)


def load_fasta_subset(
    paths: Iterable[str | Path], wanted: set[str], label: str
) -> OrderedDict[str, str]:
    """Load only requested records from one or more FASTA files."""

    records: OrderedDict[str, str] = OrderedDict()
    for path in paths:
        for identifier, sequence in iter_fasta_records(path):
            if identifier not in wanted:
                continue
            if identifier in records:
                raise ValueError(f"Duplicate {label} FASTA identifier across inputs: {identifier}")
            records[identifier] = sequence
    missing = sorted(wanted.difference(records))
    if missing:
        raise ValueError(
            f"Missing {label} sequences for {len(missing)} candidate IDs. Examples: "
            + ", ".join(missing[:10])
        )
    return records


hmm_hits: dict[str, dict[str, Any]] = {}
hmm_path = Path(snakemake.input.hmm)
if hmm_path.exists() and hmm_path.stat().st_size:
    with hmm_path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.split(maxsplit=22)
            if len(fields) < 22:
                continue
            target = fields[0]
            full_evalue = float(fields[6])
            domain_ievalue = float(fields[12])
            if not bool(snakemake.params.hmm_cut_ga):
                if full_evalue > float(snakemake.params.hmm_evalue):
                    continue
                if domain_ievalue > float(snakemake.params.hmm_domain_evalue):
                    continue
            hit = {
                "hmm_query": fields[3],
                "hmm_full_evalue": full_evalue,
                "hmm_domain_ievalue": domain_ievalue,
                "hmm_domain_score": float(fields[13]),
                "hmm_ali_from": int(fields[17]),
                "hmm_ali_to": int(fields[18]),
            }
            current = hmm_hits.get(target)
            if current is None or domain_ievalue < float(current["hmm_domain_ievalue"]):
                hmm_hits[target] = hit

blast_hits: dict[str, dict[str, Any]] = {}
blast_path = Path(snakemake.input.blast)
if blast_path.exists() and blast_path.stat().st_size:
    columns = [
        "query_id",
        "stable_id",
        "identity",
        "alignment_length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
        "query_length",
        "subject_length",
    ]
    blast = pd.read_csv(blast_path, sep="\t", names=columns, comment="#")
    if not blast.empty:
        blast["query_coverage"] = blast["alignment_length"] / blast["query_length"] * 100.0
        blast = blast.loc[
            (blast["evalue"] <= float(snakemake.params.blast_evalue))
            & (blast["identity"] >= float(snakemake.params.blast_min_identity))
            & (blast["query_coverage"] >= float(snakemake.params.blast_min_query_coverage))
        ]
        blast = blast.sort_values(
            ["stable_id", "evalue", "bitscore"], ascending=[True, True, False]
        ).drop_duplicates("stable_id")
        for row in blast.to_dict(orient="records"):
            blast_hits[str(row["stable_id"])] = {
                "blast_query_id": row["query_id"],
                "blast_best_evalue": float(row["evalue"]),
                "blast_best_identity": float(row["identity"]),
                "blast_best_query_coverage": float(row["query_coverage"]),
                "blast_best_bitscore": float(row["bitscore"]),
            }

precomputed_ids: set[str] = set()
precomputed_path = str(snakemake.params.precomputed_members or "").strip()
if precomputed_path:
    precomputed = pd.read_csv(precomputed_path, sep=None, engine="python")
    stable_column = resolve_column(precomputed, ["stable_id", "protein_id"], required=False)
    if stable_column:
        precomputed_ids.update(precomputed[stable_column].dropna().astype(str))
    else:
        species_column = resolve_column(precomputed, ["species_id", "species"])
        gene_column = resolve_column(precomputed, ["gene_id", "gene"])
        separator = str(snakemake.params.separator)
        precomputed_ids.update(
            f"{species}{separator}{gene}"
            for species, gene in zip(
                precomputed[species_column].astype(str),
                precomputed[gene_column].astype(str),
                strict=True,
            )
        )

hmm_ids = set(hmm_hits)
blast_ids = set(blast_hits)
if precomputed_path:
    selected_ids = set(precomputed_ids)
    candidate_ids = set(precomputed_ids)
    selection_basis = "precomputed"
else:
    mode = str(snakemake.params.combine_evidence)
    candidate_ids = hmm_ids | blast_ids
    if mode == "intersection":
        selected_ids = hmm_ids & blast_ids
    elif mode == "hmm_only":
        selected_ids = hmm_ids
    elif mode == "blast_only":
        selected_ids = blast_ids
    else:
        selected_ids = candidate_ids
    selection_basis = mode

if not selected_ids:
    raise RuntimeError("No family member passed the configured evidence rule.")

# Canonical mapping files can contain millions of rows.  Filter each file before
# concatenation so memory use scales with the candidate family rather than all genes.
map_frames: list[pd.DataFrame] = []
for path in snakemake.input.maps:
    frame = pd.read_csv(path, sep="\t")
    if "stable_id" not in frame.columns:
        raise ValueError(f"Canonical mapping lacks stable_id: {path}")
    subset = frame.loc[frame["stable_id"].astype(str).isin(candidate_ids)].copy()
    if not subset.empty:
        map_frames.append(subset)
if not map_frames:
    raise ValueError("No evidence IDs were found in the canonical mapping tables.")
maps = pd.concat(map_frames, ignore_index=True)
if maps["stable_id"].duplicated().any():
    duplicates = maps.loc[maps["stable_id"].duplicated(keep=False), "stable_id"].astype(str)
    raise ValueError(f"Duplicate stable IDs in canonical maps: {sorted(duplicates.unique())[:10]}")
map_by_id = maps.set_index("stable_id", drop=False)
unknown = sorted(candidate_ids.difference(map_by_id.index.astype(str)))
if unknown:
    raise ValueError(
        "Evidence contains IDs absent from canonical mapping. First examples: "
        + ", ".join(unknown[:10])
    )

proteins = load_fasta_subset(snakemake.input.proteins, candidate_ids, "protein")
cds = load_fasta_subset(snakemake.input.cds, selected_ids, "CDS")

subfamily_by_id: dict[str, str] = {}
subfamily_path = str(snakemake.params.subfamily_assignments or "").strip()
if subfamily_path:
    table = pd.read_csv(subfamily_path, sep=None, engine="python")
    stable_column = resolve_column(table, ["stable_id", "protein_id"], required=False)
    subfamily_column = resolve_column(table, ["subfamily", "clade", "ogg"])
    if stable_column:
        subfamily_by_id.update(
            zip(table[stable_column].astype(str), table[subfamily_column].astype(str), strict=True)
        )
    else:
        species_column = resolve_column(table, ["species_id", "species"])
        gene_column = resolve_column(table, ["gene_id", "gene"])
        separator = str(snakemake.params.separator)
        for species, gene, subfamily in zip(
            table[species_column].astype(str),
            table[gene_column].astype(str),
            table[subfamily_column].astype(str),
            strict=True,
        ):
            subfamily_by_id[f"{species}{separator}{gene}"] = subfamily


def optional_annotations(path_value: str, prefix: str) -> pd.DataFrame | None:
    if not path_value:
        return None
    table = pd.read_csv(path_value, sep=None, engine="python")
    stable_column = resolve_column(table, ["stable_id", "protein_id"], required=False)
    if stable_column is None:
        species_column = resolve_column(table, ["species_id", "species"])
        gene_column = resolve_column(table, ["gene_id", "gene"])
        separator = str(snakemake.params.separator)
        table["stable_id"] = [
            f"{species}{separator}{gene}"
            for species, gene in zip(
                table[species_column].astype(str), table[gene_column].astype(str), strict=True
            )
        ]
    else:
        table = table.rename(columns={stable_column: "stable_id"})
    if table["stable_id"].duplicated().any():
        duplicates = table.loc[table["stable_id"].duplicated(keep=False), "stable_id"].astype(str)
        raise ValueError(
            f"{prefix} annotation contains duplicate stable IDs: "
            + ", ".join(sorted(duplicates.unique())[:10])
        )
    renamed = {
        column: f"{prefix}_{column}"
        for column in table.columns
        if column != "stable_id" and not column.startswith(f"{prefix}_")
    }
    return table.rename(columns=renamed)


def protein_properties(sequence: str) -> tuple[object, object, str]:
    raw = sequence.upper().rstrip("*")
    clean = "".join(amino_acid for amino_acid in raw if amino_acid in "ACDEFGHIKLMNPQRSTVWY")
    status = "NON_STANDARD_RESIDUES_REMOVED" if len(clean) != len(raw) else "PASS"
    if not clean:
        return pd.NA, pd.NA, "NO_STANDARD_RESIDUES"
    analysis = ProteinAnalysis(clean)
    return analysis.molecular_weight(), analysis.isoelectric_point(), status


rows: list[dict[str, Any]] = []
rejected_rows: list[dict[str, Any]] = []
domain_sequences: OrderedDict[str, str] = OrderedDict()
for stable_id in sorted(candidate_ids):
    base = map_by_id.loc[stable_id].to_dict()
    row: dict[str, Any] = {
        **base,
        "family_name": snakemake.params.family_name,
        "evidence_hmm": stable_id in hmm_ids,
        "evidence_blast": stable_id in blast_ids,
        "evidence_precomputed": stable_id in precomputed_ids,
        "selection_basis": selection_basis,
        "subfamily": subfamily_by_id.get(stable_id, pd.NA),
    }
    row.update(hmm_hits.get(stable_id, {}))
    row.update(blast_hits.get(stable_id, {}))
    if bool(snakemake.params.calculate_properties):
        mw, pi, property_qc = protein_properties(proteins[stable_id])
        row.update(
            {
                "molecular_weight_da": mw,
                "theoretical_pi": pi,
                "protein_property_qc": property_qc,
            }
        )
    if stable_id in selected_ids:
        row["decision"] = "PASS"
        row["rejection_reason"] = pd.NA
        rows.append(row)
        hit = hmm_hits.get(stable_id)
        if hit:
            start = max(1, int(hit["hmm_ali_from"])) - 1
            end = min(len(proteins[stable_id]), int(hit["hmm_ali_to"]))
            if end > start:
                domain_sequences[stable_id] = proteins[stable_id][start:end]
    else:
        row["decision"] = "REJECT"
        missing: list[str] = []
        mode = str(snakemake.params.combine_evidence)
        if mode == "intersection":
            if stable_id not in hmm_ids:
                missing.append("missing_hmm_evidence")
            if stable_id not in blast_ids:
                missing.append("missing_blast_evidence")
        elif mode == "hmm_only":
            missing.append("missing_hmm_evidence")
        elif mode == "blast_only":
            missing.append("missing_blast_evidence")
        row["rejection_reason"] = ";".join(missing) or "did_not_pass_selection_rule"
        rejected_rows.append(row)

members = pd.DataFrame(rows)
for annotation, prefix in (
    (str(snakemake.params.domain_validation_table or ""), "domain_validation"),
    (str(snakemake.params.subcellular_localization_table or ""), "localization"),
):
    annotation_table = optional_annotations(annotation, prefix)
    if annotation_table is not None:
        members = members.merge(annotation_table, on="stable_id", how="left", validate="one_to_one")

rejected = pd.DataFrame(rejected_rows)
columns_front = [
    "species_id",
    "gene_id",
    "transcript_id",
    "stable_id",
    "family_name",
    "subfamily",
    "group",
    "species_subfamily",
    "chromosome",
    "gene_start",
    "gene_end",
    "strand",
]
members = members[
    [column for column in columns_front if column in members.columns]
    + [column for column in members.columns if column not in columns_front]
]
if rejected.empty:
    rejected = members.iloc[0:0].copy()

save_table(members, snakemake.output.members)
save_table(rejected, snakemake.output.rejected)
save_workbook({"members": members, "rejected": rejected}, snakemake.output.xlsx)
write_fasta(
    OrderedDict((stable_id, proteins[stable_id]) for stable_id in sorted(selected_ids)),
    snakemake.output.proteins,
)
write_fasta(
    OrderedDict((stable_id, cds[stable_id]) for stable_id in sorted(selected_ids)),
    snakemake.output.cds,
)
write_fasta(domain_sequences, snakemake.output.domains)
