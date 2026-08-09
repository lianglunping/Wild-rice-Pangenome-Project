import sys
from pathlib import Path as _ScriptPath

sys.path.insert(0, str(_ScriptPath(snakemake.scriptdir)))

from collections import defaultdict
from pathlib import Path

import pandas as pd
from workflow_utils import (
    copy_atomic,
    first_parent,
    iter_gff,
    read_fasta,
    run_command,
    save_table,
    write_fasta,
)

species = str(snakemake.params.species)
separator = str(snakemake.params.separator)
work_dir = Path(snakemake.params.work_dir)
work_dir.mkdir(parents=True, exist_ok=True)
raw_gff = work_dir / "canonical.raw.gff3"
raw_transcripts = work_dir / "transcripts.raw.fa"
raw_cds = work_dir / "cds.raw.fa"
raw_proteins = work_dir / "proteins.raw.fa"

run_command(
    [
        "agat_sp_keep_longest_isoform.pl",
        "--gff",
        str(snakemake.input.gff3),
        "--output",
        str(raw_gff),
    ],
    stdout_path=snakemake.log.agat_stdout,
    stderr_path=snakemake.log.agat_stderr,
)
run_command(
    [
        "gffread",
        str(raw_gff),
        "-g",
        str(snakemake.input.genome),
        "-w",
        str(raw_transcripts),
        "-x",
        str(raw_cds),
        "-y",
        str(raw_proteins),
    ],
    stdout_path=snakemake.log.gffread_stdout,
    stderr_path=snakemake.log.gffread_stderr,
)

transcript_types = {"mrna", "transcript", "ncrna", "trna", "rrna"}
genes: dict[str, dict[str, object]] = {}
transcript_to_gene: dict[str, str] = {}
transcript_coords: dict[str, tuple[int, int]] = {}
children: defaultdict[str, list[dict[str, object]]] = defaultdict(list)

for feature in iter_gff(raw_gff):
    attrs = feature["attributes"]
    feature_name = str(feature["feature"]).lower()
    if feature_name == "gene":
        gene_id = attrs.get("ID") or attrs.get("gene_id") or attrs.get("Name")
        if gene_id:
            genes[gene_id] = {
                "chromosome": feature["seqid"],
                "gene_start": feature["start"],
                "gene_end": feature["end"],
                "strand": feature["strand"],
            }
    elif feature_name in transcript_types:
        transcript_id = attrs.get("ID") or attrs.get("transcript_id")
        gene_id = first_parent(attrs.get("Parent")) or attrs.get("gene_id")
        if transcript_id and gene_id:
            transcript_to_gene[transcript_id] = gene_id
            transcript_coords[transcript_id] = (int(feature["start"]), int(feature["end"]))
    else:
        for parent in (attrs.get("Parent") or attrs.get("transcript_id") or "").split(","):
            if parent:
                children[parent].append(feature)

proteins = read_fasta(raw_proteins)
cds = read_fasta(raw_cds)
transcripts = read_fasta(raw_transcripts)


def sequence_for(records: dict[str, str], transcript_id: str) -> str:
    candidates = [
        transcript_id,
        transcript_id.removeprefix("transcript:"),
        f"transcript:{transcript_id}",
    ]
    for candidate in candidates:
        if candidate in records:
            return records[candidate]
    suffix_matches = [value for key, value in records.items() if key.split("|")[0] == transcript_id]
    if len(suffix_matches) == 1:
        return suffix_matches[0]
    raise KeyError(f"No unique sequence found for transcript {transcript_id!r}")


rows: list[dict[str, object]] = []
protein_output: dict[str, str] = {}
cds_output: dict[str, str] = {}
transcript_output: dict[str, str] = {}
for transcript_id, gene_id in sorted(
    transcript_to_gene.items(), key=lambda item: (item[1], item[0])
):
    protein_sequence = sequence_for(proteins, transcript_id)
    cds_sequence = sequence_for(cds, transcript_id)
    transcript_sequence = sequence_for(transcripts, transcript_id)
    stable_id = f"{species}{separator}{gene_id}"
    if stable_id in protein_output:
        raise ValueError(
            f"More than one canonical transcript remained for gene {gene_id!r} in {species}."
        )
    protein_output[stable_id] = protein_sequence.rstrip("*")
    cds_output[stable_id] = cds_sequence
    transcript_output[stable_id] = transcript_sequence
    gene = genes.get(gene_id, {})
    transcript_start, transcript_end = transcript_coords.get(transcript_id, (pd.NA, pd.NA))
    rows.append(
        {
            "species_id": species,
            "species_name": snakemake.params.species_name,
            "group": snakemake.params.group or pd.NA,
            "species_subfamily": snakemake.params.subfamily or pd.NA,
            "gene_id": gene_id,
            "transcript_id": transcript_id,
            "stable_id": stable_id,
            "chromosome": gene.get("chromosome", pd.NA),
            "gene_start": gene.get("gene_start", pd.NA),
            "gene_end": gene.get("gene_end", pd.NA),
            "strand": gene.get("strand", pd.NA),
            "transcript_start": transcript_start,
            "transcript_end": transcript_end,
            "protein_length": len(protein_output[stable_id]),
            "cds_length": len(cds_sequence),
            "transcript_length": len(transcript_sequence),
            "cds_modulo_3": len(cds_sequence) % 3,
        }
    )

if not rows:
    raise RuntimeError(f"No canonical gene/transcript mapping could be derived for {species}.")

Path(snakemake.output.gff3).parent.mkdir(parents=True, exist_ok=True)
copy_atomic(raw_gff, snakemake.output.gff3)
write_fasta(protein_output, snakemake.output.proteins)
write_fasta(cds_output, snakemake.output.cds)
write_fasta(transcript_output, snakemake.output.transcripts)
save_table(pd.DataFrame(rows), snakemake.output.mapping, snakemake.output.mapping_xlsx)
