import sys
from pathlib import Path as _ScriptPath

sys.path.insert(0, str(_ScriptPath(snakemake.scriptdir)))

import concurrent.futures
import hashlib
import itertools
import json
import math
import shutil
from collections import OrderedDict
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import pandas as pd
from workflow_utils import (
    commit_partial,
    partial_path,
    read_fasta,
    run_command,
    save_table,
    save_workbook,
    sha256_file,
    write_fasta,
    write_json,
    write_text_atomic,
)

proteins = read_fasta(snakemake.input.proteins)
cds = read_fasta(snakemake.input.cds)
family_ids = set(proteins)
pair_source = str(snakemake.params.pair_source)
reference_species = str(snakemake.params.reference_species)
separator = str(snakemake.params.separator)
pairs: dict[tuple[str, str], dict[str, Any]] = {}


def add_pair(gene1: str, gene2: str, pair_type: str, group_id: str) -> None:
    if gene1 == gene2 or gene1 not in family_ids or gene2 not in family_ids:
        return
    key = tuple(sorted((gene1, gene2)))
    if key not in pairs:
        pairs[key] = {
            "stable_id_1": key[0],
            "stable_id_2": key[1],
            "pair_type": pair_type,
            "group_id": group_id,
        }
    elif pair_type not in str(pairs[key]["pair_type"]).split(";"):
        pairs[key]["pair_type"] = f"{pairs[key]['pair_type']};{pair_type}"
        pairs[key]["group_id"] = f"{pairs[key]['group_id']};{group_id}"


if pair_source in {"orthology", "both"}:
    membership_value = snakemake.input.membership
    membership_path = Path(
        str(
            membership_value[0] if isinstance(membership_value, (list, tuple)) else membership_value
        )
    )
    if not membership_path.is_file() or membership_path.stat().st_size == 0:
        raise FileNotFoundError(
            "Orthology pair source selected but target-family HOG membership is unavailable."
        )
    membership = pd.read_csv(membership_path, sep="\t")
    for hog_id, group in membership.groupby("HOG_ID", sort=True):
        copy_counts = group.groupby("species_id")["stable_id"].nunique()
        if reference_species not in copy_counts.index or copy_counts.loc[reference_species] != 1:
            continue
        reference_gene = str(
            group.loc[group["species_id"].astype(str) == reference_species, "stable_id"].iloc[0]
        )
        candidates = group.loc[
            (group["species_id"].astype(str) != reference_species)
            & group["species_id"].map(copy_counts).eq(1)
        ]
        for stable_id in candidates["stable_id"].astype(str):
            add_pair(reference_gene, stable_id, "orthology_single_copy_to_reference", str(hog_id))

if pair_source in {"duplication", "both"}:
    duplication_value = snakemake.input.duplication_pairs
    duplication_path = Path(
        str(
            duplication_value[0]
            if isinstance(duplication_value, (list, tuple))
            else duplication_value
        )
    )
    if not duplication_path.is_file() or duplication_path.stat().st_size == 0:
        raise FileNotFoundError(
            "Duplication pair source selected but duplication pairs are unavailable."
        )
    duplication = pd.read_csv(duplication_path, sep="\t")
    if not duplication.empty:
        for row in duplication.to_dict(orient="records"):
            add_pair(
                str(row["stable_id_1"]),
                str(row["stable_id_2"]),
                f"duplication_{row.get('duplication_mode', 'unknown')}",
                str(row.get("species_id", "")),
            )

pair_records = list(pairs.values())
limit = snakemake.params.max_pairs_per_group
if limit not in (None, "", "None"):
    limited: list[dict[str, Any]] = []
    for _, group in itertools.groupby(
        sorted(pair_records, key=lambda row: (row["pair_type"], row["group_id"])),
        key=lambda row: (row["pair_type"], row["group_id"]),
    ):
        limited.extend(list(group)[: int(limit)])
    pair_records = limited
if not pair_records:
    raise RuntimeError(
        "No constrained orthology/duplication sequence pairs were available for Ka/Ks."
    )

kaks_executable = next(
    (candidate for candidate in ("KaKs_Calculator", "KaKs") if shutil.which(candidate)), None
)
if kaks_executable is None:
    raise FileNotFoundError("Neither KaKs_Calculator nor KaKs executable is available in PATH.")
work_root = Path(snakemake.params.work_dir)
work_root.mkdir(parents=True, exist_ok=True)
method = str(snakemake.params.method)
saturation = float(snakemake.params.saturation_ks)
global_signature = hashlib.sha256(
    json.dumps(
        {
            "proteins_sha256": sha256_file(snakemake.input.proteins),
            "cds_sha256": sha256_file(snakemake.input.cds),
            "method": method,
            "saturation_ks": saturation,
        },
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
).hexdigest()


def run_pair(index_record: tuple[int, dict[str, Any]]) -> dict[str, Any]:
    _, record = index_record
    gene1 = str(record["stable_id_1"])
    gene2 = str(record["stable_id_2"])
    pair_signature = hashlib.sha256(
        json.dumps(
            {
                "global_signature": global_signature,
                "gene1": gene1,
                "gene2": gene2,
                "pair_type": record["pair_type"],
                "group_id": record["group_id"],
            },
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
    ).hexdigest()
    pair_id = f"pair_{pair_signature[:16]}"
    pair_dir = work_root / pair_id
    pair_dir.mkdir(parents=True, exist_ok=True)
    state_path = pair_dir / "result.json"
    if state_path.is_file():
        try:
            cached = json.loads(state_path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            cached = {}
        if cached.get("pair_signature") == pair_signature and cached.get("qc_status") != "FAILED":
            cached_result = dict(cached["result"])
            cached_result["resumed_from_cache"] = True
            return cached_result

    result: dict[str, Any] = {
        **record,
        "pair_id": pair_id,
        "species_id_1": gene1.split(separator, 1)[0] if separator in gene1 else pd.NA,
        "gene_id_1": gene1.split(separator, 1)[1] if separator in gene1 else gene1,
        "species_id_2": gene2.split(separator, 1)[0] if separator in gene2 else pd.NA,
        "gene_id_2": gene2.split(separator, 1)[1] if separator in gene2 else gene2,
        "Ka": pd.NA,
        "Ks": pd.NA,
        "Ka_Ks": pd.NA,
        "method": method,
        "qc_status": "FAILED",
        "qc_message": "",
        "resumed_from_cache": False,
    }
    try:
        if len(cds[gene1].replace("-", "")) % 3 or len(cds[gene2].replace("-", "")) % 3:
            raise ValueError("CDS length is not divisible by three")
        protein_fasta = pair_dir / "pair.pep.fa"
        cds_fasta = pair_dir / "pair.cds.fa"
        protein_alignment = pair_dir / "pair.pep.aln.fa"
        codon_alignment = pair_dir / "pair.codon.aln.fa"
        axt = pair_dir / "pair.axt"
        output = pair_dir / "kaks.tsv"
        write_fasta(
            OrderedDict(((gene1, proteins[gene1]), (gene2, proteins[gene2]))), protein_fasta
        )
        write_fasta(OrderedDict(((gene1, cds[gene1]), (gene2, cds[gene2]))), cds_fasta)
        run_command(
            ["mafft", "--auto", "--quiet", str(protein_fasta)],
            stdout_path=protein_alignment,
            stderr_path=pair_dir / "mafft.stderr.log",
        )
        run_command(
            ["pal2nal.pl", str(protein_alignment), str(cds_fasta), "-output", "fasta"],
            stdout_path=codon_alignment,
            stderr_path=pair_dir / "pal2nal.stderr.log",
        )
        aligned = read_fasta(codon_alignment)
        if set(aligned) != {gene1, gene2}:
            raise ValueError("PAL2NAL output identifiers do not match input pair")
        write_text_atomic(f"{pair_id}\n{aligned[gene1]}\n{aligned[gene2]}\n", axt)
        output_temporary = partial_path(output)
        run_command(
            [kaks_executable, "-i", str(axt), "-o", str(output_temporary), "-m", method],
            stdout_path=pair_dir / "kaks.stdout.log",
            stderr_path=pair_dir / "kaks.stderr.log",
        )
        if not output_temporary.is_file():
            raise FileNotFoundError(f"Ka/Ks calculator did not produce {output_temporary}")
        commit_partial(output_temporary, output)
        table = pd.read_csv(output, sep="\t")
        if table.empty:
            raise ValueError("Ka/Ks output table is empty")
        normalized = {
            column.lower().replace("/", "_").replace(" ", "_"): column for column in table.columns
        }
        ka_column = normalized.get("ka")
        ks_column = normalized.get("ks")
        ratio_column = normalized.get("ka_ks") or normalized.get("ka_ks_ratio")
        if ka_column is None or ks_column is None:
            raise ValueError(f"Cannot find Ka/Ks columns in: {list(table.columns)}")
        ka = float(table.iloc[0][ka_column])
        ks = float(table.iloc[0][ks_column])
        ratio = (
            float(table.iloc[0][ratio_column])
            if ratio_column is not None and pd.notna(table.iloc[0][ratio_column])
            else (ka / ks if ks != 0 else math.inf)
        )
        flags: list[str] = []
        if ks == 0:
            flags.append("KS_ZERO")
        if ks >= saturation:
            flags.append("POTENTIAL_SATURATION")
        if not math.isfinite(ratio):
            flags.append("RATIO_NON_FINITE")
        result.update(
            {
                "Ka": ka,
                "Ks": ks,
                "Ka_Ks": ratio,
                "qc_status": ";".join(flags) if flags else "PASS",
                "qc_message": "",
            }
        )
    except Exception as error:
        result["qc_message"] = f"{type(error).__name__}: {error}"
    serializable_result = {
        key: (None if pd.isna(value) else value) for key, value in result.items()
    }
    write_json(
        {
            "pair_signature": pair_signature,
            "qc_status": result["qc_status"],
            "result": serializable_result,
        },
        state_path,
    )
    return result


with concurrent.futures.ThreadPoolExecutor(max_workers=int(snakemake.params.workers)) as executor:
    rows = list(executor.map(run_pair, enumerate(pair_records, start=1)))
results = pd.DataFrame(rows).sort_values("pair_id").reset_index(drop=True)
save_table(results, snakemake.output.tsv)
save_workbook({"kaks_pairs": results}, snakemake.output.xlsx)

valid = results.loc[pd.to_numeric(results["Ka_Ks"], errors="coerce").notna()].copy()
fig, axis = plt.subplots(figsize=(8.0, 4.8))
if valid.empty:
    axis.text(0.5, 0.5, "No valid Ka/Ks estimates", ha="center", va="center")
    axis.set_axis_off()
else:
    groups = [group["Ka_Ks"].astype(float).values for _, group in valid.groupby("pair_type")]
    labels = [str(name) for name, _ in valid.groupby("pair_type")]
    axis.boxplot(groups, tick_labels=labels, showfliers=False)
    axis.set_ylabel("Ka/Ks")
    axis.tick_params(axis="x", rotation=30)
    axis.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(snakemake.output.plot_pdf)
fig.savefig(snakemake.output.plot_png, dpi=int(snakemake.params.png_dpi))
plt.close(fig)
