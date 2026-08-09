import sys
from pathlib import Path as _ScriptPath

sys.path.insert(0, str(_ScriptPath(snakemake.scriptdir)))

from collections import Counter
from pathlib import Path

import pandas as pd
from workflow_utils import fasta_lengths, iter_gff, save_table, sha256_file, write_json

records = list(snakemake.params.records)
calculate_sha256 = bool(snakemake.params.calculate_sha256)
rows: list[dict[str, object]] = []

for record in records:
    path = Path(record["path"])
    row: dict[str, object] = {
        "species_id": record.get("species_id", ""),
        "role": record["role"],
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": sha256_file(path) if calculate_sha256 else "NOT_CALCULATED",
        "record_count": pd.NA,
        "total_length": pd.NA,
        "min_length": pd.NA,
        "max_length": pd.NA,
        "mean_length": pd.NA,
        "gff_seqid_count": pd.NA,
        "gff_gene_count": pd.NA,
        "gff_transcript_count": pd.NA,
        "gff_cds_count": pd.NA,
        "status": "PASS",
        "notes": "",
    }
    role = str(record["role"])
    try:
        if role in {"genome", "protein", "cds", "reference_proteins", "hmm"} and role != "hmm":
            lengths = list(fasta_lengths(path).values())
            total_length = sum(lengths)
            row.update(
                {
                    "record_count": len(lengths),
                    "total_length": total_length,
                    "min_length": min(lengths) if lengths else 0,
                    "max_length": max(lengths) if lengths else 0,
                    "mean_length": total_length / len(lengths) if lengths else 0.0,
                }
            )
            if not lengths:
                row["status"] = "FAIL"
                row["notes"] = "No FASTA records"
        elif role == "gff3":
            counts: Counter[str] = Counter()
            seqids: set[str] = set()
            for feature in iter_gff(path):
                counts[str(feature["feature"]).lower()] += 1
                seqids.add(str(feature["seqid"]))
            transcript_count = sum(
                counts[name] for name in ("mrna", "transcript", "ncrna", "trna", "rrna")
            )
            row.update(
                {
                    "record_count": sum(counts.values()),
                    "gff_seqid_count": len(seqids),
                    "gff_gene_count": counts["gene"],
                    "gff_transcript_count": transcript_count,
                    "gff_cds_count": counts["cds"],
                }
            )
            if counts["gene"] == 0:
                row["status"] = "WARN"
                row["notes"] = "No explicit gene features; downstream mapping may need review"
    except Exception as error:  # input audit records failures instead of hiding them
        row["status"] = "FAIL"
        row["notes"] = f"{type(error).__name__}: {error}"
    rows.append(row)

table = pd.DataFrame(rows)
save_table(table, snakemake.output.tsv, snakemake.output.xlsx)
manifest = {
    "project": snakemake.config.get("project", {}),
    "config_path": snakemake.config.get("panfamflow_config_path"),
    "selected_modules": snakemake.config.get("panfamflow_selected_modules"),
    "audit_records": len(rows),
    "failed_records": int((table["status"] == "FAIL").sum()),
}
write_json(manifest, snakemake.output.manifest)
if manifest["failed_records"]:
    raise RuntimeError(
        f"Input audit detected {manifest['failed_records']} failed record(s); inspect {snakemake.output.tsv}."
    )
