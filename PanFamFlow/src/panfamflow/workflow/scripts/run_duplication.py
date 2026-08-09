import sys
from pathlib import Path as _ScriptPath

sys.path.insert(0, str(_ScriptPath(snakemake.scriptdir)))

import shutil
from collections import defaultdict
from itertools import chain
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import pandas as pd
from workflow_utils import (
    read_fasta,
    resolve_column,
    run_command,
    save_table,
    save_workbook,
    split_multi_value,
    write_fasta,
)

members = pd.read_csv(snakemake.input.members, sep="\t")
family_ids = set(members["stable_id"].astype(str))
separator = str(snakemake.params.separator)
backend = str(snakemake.params.backend)
mode_rows: list[dict[str, Any]] = []
pair_rows: list[dict[str, Any]] = []

if backend == "precomputed":
    table = pd.read_csv(snakemake.params.precomputed_table, sep=None, engine="python")
    stable_column = resolve_column(table, ["stable_id", "protein_id"], required=False)
    if stable_column is None:
        species_column = resolve_column(table, ["species_id", "species"])
        gene_column = resolve_column(table, ["gene_id", "gene"])
        table["stable_id"] = [
            f"{species}{separator}{gene}"
            for species, gene in zip(
                table[species_column].astype(str), table[gene_column].astype(str), strict=True
            )
        ]
    else:
        table = table.rename(columns={stable_column: "stable_id"})
    mode_column = resolve_column(table, ["duplication_mode", "mode", "type"])
    partner_column = resolve_column(
        table, ["partner_stable_id", "partner_gene_id", "partner"], required=False
    )
    table = table.loc[table["stable_id"].astype(str).isin(family_ids)].copy()
    table["duplication_mode"] = table[mode_column].astype(str)
    if partner_column:
        normalized_partners: list[str | object] = []
        for row in table.to_dict(orient="records"):
            stable_id = str(row["stable_id"])
            species_id = stable_id.split(separator, 1)[0] if separator in stable_id else ""
            partners = []
            for partner in split_multi_value(row.get(partner_column)):
                partner_stable = (
                    partner if separator in partner else f"{species_id}{separator}{partner}"
                )
                partners.append(partner_stable)
                pair_rows.append(
                    {
                        "species_id": species_id,
                        "stable_id_1": stable_id,
                        "stable_id_2": partner_stable,
                        "duplication_mode": str(row[mode_column]),
                        "source_file": str(Path(snakemake.params.precomputed_table).resolve()),
                    }
                )
            normalized_partners.append(";".join(sorted(set(partners))) or pd.NA)
        table["partner_stable_ids"] = normalized_partners
    else:
        table["partner_stable_ids"] = pd.NA
    if table.groupby("stable_id")["duplication_mode"].nunique().gt(1).any():
        conflicting = table.groupby("stable_id")["duplication_mode"].nunique()
        conflicting = conflicting.loc[conflicting > 1].index.tolist()
        raise ValueError(
            "Precomputed duplication table is not uniquely assigned; examples: "
            + ", ".join(conflicting[:10])
        )
    table = table.drop_duplicates(["stable_id", "duplication_mode"])
    mode_rows = table.to_dict(orient="records")
else:
    species_records = {str(record["id"]): record for record in snakemake.params.species_records}
    target_ids = [str(item) for item in snakemake.params.targets]
    species_ids = [str(item) for item in snakemake.params.species_ids]
    map_paths = dict(zip(species_ids, snakemake.input.maps, strict=True))
    protein_paths = dict(zip(species_ids, snakemake.input.proteins, strict=True))
    work_root = Path(snakemake.params.work_dir)
    work_root.mkdir(parents=True, exist_ok=True)
    executable = str(snakemake.params.dupgen_executable)
    if shutil.which(executable) is None and not Path(executable).is_file():
        raise FileNotFoundError(
            f"DupGen_finder executable not found: {executable}. Run scripts/install_dupgen.sh first."
        )

    def write_dupgen_gff(table: pd.DataFrame, path: Path) -> None:
        with path.open("w", encoding="utf-8") as handle:
            for row in table.sort_values(["chromosome", "gene_start", "gene_end"]).to_dict(
                orient="records"
            ):
                if pd.isna(row.get("chromosome")) or pd.isna(row.get("gene_start")):
                    continue
                handle.write(
                    f"{row['species_id']}-{row['chromosome']}\t{row['stable_id']}\t"
                    f"{int(row['gene_start'])}\t{int(row['gene_end'])}\n"
                )

    for target in target_ids:
        record = species_records[target]
        outgroup = str(record.get("outgroup") or "")
        if not outgroup:
            raise ValueError(f"DupGen_finder target {target} has no outgroup")
        target_map = pd.read_csv(map_paths[target], sep="\t")
        outgroup_map = pd.read_csv(map_paths[outgroup], sep="\t")
        target_proteins = read_fasta(protein_paths[target])
        outgroup_proteins = read_fasta(protein_paths[outgroup])
        run_dir = work_root / target
        if run_dir.exists():
            shutil.rmtree(run_dir)
        data_dir = run_dir / "data"
        output_dir = run_dir / "output"
        data_dir.mkdir(parents=True)
        output_dir.mkdir(parents=True)

        write_dupgen_gff(target_map, data_dir / f"{target}.gff")
        combined_map = pd.concat([target_map, outgroup_map], ignore_index=True)
        write_dupgen_gff(combined_map, data_dir / f"{target}_{outgroup}.gff")
        target_fasta = data_dir / f"{target}.pep.fa"
        combined_fasta = data_dir / f"{target}_{outgroup}.pep.fa"
        write_fasta(target_proteins, target_fasta)
        overlap = set(target_proteins).intersection(outgroup_proteins)
        if overlap:
            raise ValueError(
                f"Stable IDs overlap between {target} and {outgroup}: {sorted(overlap)[:5]}"
            )
        write_fasta(chain(target_proteins.items(), outgroup_proteins.items()), combined_fasta)

        target_db = data_dir / target
        combined_db = data_dir / f"{target}_{outgroup}"
        run_command(
            ["diamond", "makedb", "--in", str(target_fasta), "-d", str(target_db)],
            stdout_path=run_dir / "diamond_target_makedb.stdout.log",
            stderr_path=run_dir / "diamond_target_makedb.stderr.log",
        )
        run_command(
            [
                "diamond",
                "blastp",
                "-d",
                str(target_db),
                "-q",
                str(target_fasta),
                "-o",
                str(data_dir / f"{target}.blast"),
                "--threads",
                str(snakemake.threads),
                "--max-target-seqs",
                str(snakemake.params.max_target_seqs),
                "--evalue",
                str(snakemake.params.diamond_evalue),
                "--outfmt",
                "6",
            ],
            stdout_path=run_dir / "diamond_target.stdout.log",
            stderr_path=run_dir / "diamond_target.stderr.log",
        )
        run_command(
            ["diamond", "makedb", "--in", str(combined_fasta), "-d", str(combined_db)],
            stdout_path=run_dir / "diamond_cross_makedb.stdout.log",
            stderr_path=run_dir / "diamond_cross_makedb.stderr.log",
        )
        run_command(
            [
                "diamond",
                "blastp",
                "-d",
                str(combined_db),
                "-q",
                str(combined_fasta),
                "-o",
                str(data_dir / f"{target}_{outgroup}.blast"),
                "--threads",
                str(snakemake.threads),
                "--max-target-seqs",
                str(snakemake.params.max_target_seqs),
                "--evalue",
                str(snakemake.params.diamond_evalue),
                "--outfmt",
                "6",
            ],
            stdout_path=run_dir / "diamond_cross.stdout.log",
            stderr_path=run_dir / "diamond_cross.stderr.log",
        )
        command = [
            executable,
            "-i",
            str(data_dir.resolve()),
            "-t",
            target,
            "-c",
            outgroup,
            "-o",
            str(output_dir.resolve()),
            "-d",
            str(snakemake.params.proximal_max_gene_distance),
        ]
        command.extend(str(value) for value in snakemake.params.extra_args)
        run_command(
            command,
            stdout_path=run_dir / "dupgen.stdout.log",
            stderr_path=run_dir / "dupgen.stderr.log",
        )

        file_modes = {
            "wgd": "WGD",
            "tandem": "Tandem",
            "proximal": "Proximal",
            "transposed": "Transposed",
            "dispersed": "Dispersed",
        }
        assignments: defaultdict[str, set[str]] = defaultdict(set)
        partners: defaultdict[str, set[str]] = defaultdict(set)
        for suffix, mode in file_modes.items():
            candidates = list(output_dir.rglob(f"{target}.{suffix}.pairs"))
            if not candidates:
                continue
            with candidates[0].open("r", encoding="utf-8") as handle:
                for raw in handle:
                    fields = raw.rstrip("\n").split("\t")
                    if len(fields) < 3 or fields[0].lower().startswith("duplicate"):
                        continue
                    gene1, gene2 = fields[0], fields[2]
                    assignments[gene1].add(mode)
                    assignments[gene2].add(mode)
                    partners[gene1].add(gene2)
                    partners[gene2].add(gene1)
                    pair_rows.append(
                        {
                            "species_id": target,
                            "stable_id_1": gene1,
                            "stable_id_2": gene2,
                            "duplication_mode": mode,
                            "source_file": str(candidates[0].resolve()),
                        }
                    )
        singleton_candidates = list(output_dir.rglob(f"{target}.singletons"))
        if singleton_candidates:
            with singleton_candidates[0].open("r", encoding="utf-8") as handle:
                for raw in handle:
                    fields = raw.rstrip("\n").split("\t")
                    if not fields or fields[0].lower().startswith("gene"):
                        continue
                    assignments[fields[0]].add("Singleton")

        target_family_ids = sorted(
            stable_id for stable_id in family_ids if stable_id.startswith(f"{target}{separator}")
        )
        for stable_id in target_family_ids:
            modes = assignments.get(stable_id, set())
            if len(modes) > 1:
                raise RuntimeError(
                    f"DupGen_finder-unique returned multiple modes for {stable_id}: {sorted(modes)}"
                )
            mode = next(iter(modes), "Unclassified")
            gene_id = stable_id.split(separator, 1)[1]
            mode_rows.append(
                {
                    "species_id": target,
                    "gene_id": gene_id,
                    "stable_id": stable_id,
                    "duplication_mode": mode,
                    "partner_stable_ids": ";".join(sorted(partners.get(stable_id, set()))) or pd.NA,
                    "outgroup": outgroup,
                    "backend": "dupgen_finder_unique",
                }
            )

modes = pd.DataFrame(mode_rows)
pairs = pd.DataFrame(pair_rows)
if modes.empty:
    raise RuntimeError("No family duplication assignments were produced.")
if "species_id" not in modes.columns:
    modes = modes.merge(members[["stable_id", "species_id", "gene_id"]], on="stable_id", how="left")
columns = [
    "species_id",
    "gene_id",
    "stable_id",
    "duplication_mode",
    "partner_stable_ids",
    "outgroup",
    "backend",
]
for column in columns:
    if column not in modes.columns:
        modes[column] = pd.NA
modes = modes[columns + [column for column in modes.columns if column not in columns]]
modes = modes.sort_values(["species_id", "stable_id"]).reset_index(drop=True)
save_table(modes, snakemake.output.modes)
save_table(pairs, snakemake.output.pairs)
save_workbook({"duplication_mode": modes, "pairs": pairs}, snakemake.output.xlsx)

counts = modes["duplication_mode"].value_counts().sort_values(ascending=False)
fig, axis = plt.subplots(figsize=(7.2, 4.8))
axis.bar(counts.index, counts.values)
axis.set_xlabel("Duplication mode")
axis.set_ylabel("Number of family genes")
axis.tick_params(axis="x", rotation=35)
axis.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(snakemake.output.plot_pdf)
fig.savefig(snakemake.output.plot_png, dpi=int(snakemake.params.png_dpi))
plt.close(fig)
