import sys
from pathlib import Path as _ScriptPath

sys.path.insert(0, str(_ScriptPath(snakemake.scriptdir)))

import hashlib
import json
import re
import shutil
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd
from workflow_utils import executable_version, run_command, save_table, sha256_file

species = str(snakemake.params.species)
lineage = str(snakemake.params.lineage)
work_dir = Path(snakemake.params.work_dir)
work_dir.mkdir(parents=True, exist_ok=True)
_, busco_version = executable_version(["busco"], ["--version"])
signature_payload = {
    "genome_sha256": sha256_file(snakemake.input.genome),
    "lineage": lineage,
    "mode": str(snakemake.params.mode),
    "offline": bool(snakemake.params.offline),
    "extra_args": [str(value) for value in snakemake.params.extra_args],
    "busco_version": busco_version,
}
signature = hashlib.sha256(
    json.dumps(signature_payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
).hexdigest()
run_name = f"{species}_busco_{signature[:12]}"
run_dir = work_dir / run_name
summary_candidates = sorted(run_dir.glob("short_summary*.txt")) if run_dir.exists() else []
if not summary_candidates:
    if run_dir.exists() and any(run_dir.iterdir()):
        timestamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")
        failed_dir = work_dir / "failed" / f"{run_name}_{timestamp}"
        failed_dir.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(run_dir), failed_dir)
    command = [
        "busco",
        "-i",
        str(snakemake.input.genome),
        "-m",
        "genome",
        "-l",
        lineage,
        "-o",
        run_name,
        "--out_path",
        str(work_dir),
        "-c",
        str(snakemake.threads),
    ]
    if bool(snakemake.params.offline):
        command.append("--offline")
        command.extend(["--download_path", str(snakemake.params.download_path)])
    command.extend(str(value) for value in snakemake.params.extra_args)
    run_command(command, stdout_path=snakemake.log.stdout, stderr_path=snakemake.log.stderr)
    summary_candidates = sorted(run_dir.glob("short_summary*.txt"))

if not summary_candidates:
    raise FileNotFoundError(f"BUSCO summary was not produced under {run_dir}")
summary_text = summary_candidates[0].read_text(encoding="utf-8")
pattern = re.compile(
    r"C:(?P<complete>[0-9.]+)%\[S:(?P<single>[0-9.]+)%,D:(?P<duplicated>[0-9.]+)%\],"
    r"F:(?P<fragmented>[0-9.]+)%,M:(?P<missing>[0-9.]+)%,n:(?P<n>\d+)"
)
match = pattern.search(summary_text.replace(" ", ""))
row: dict[str, object] = {
    "species_id": species,
    "lineage": lineage,
    "complete_pct": pd.NA,
    "single_copy_pct": pd.NA,
    "duplicated_pct": pd.NA,
    "fragmented_pct": pd.NA,
    "missing_pct": pd.NA,
    "busco_n": pd.NA,
    "summary_path": str(summary_candidates[0].resolve()),
}
if match:
    row.update(
        {
            "complete_pct": float(match.group("complete")),
            "single_copy_pct": float(match.group("single")),
            "duplicated_pct": float(match.group("duplicated")),
            "fragmented_pct": float(match.group("fragmented")),
            "missing_pct": float(match.group("missing")),
            "busco_n": int(match.group("n")),
        }
    )
save_table(pd.DataFrame([row]), snakemake.output.tsv)
