import sys
from pathlib import Path as _ScriptPath

sys.path.insert(0, str(_ScriptPath(snakemake.scriptdir)))

import hashlib
import json
import shutil
from datetime import UTC, datetime
from pathlib import Path

from workflow_utils import (
    executable_version,
    run_command,
    sha256_file,
    write_json,
    write_text_atomic,
)

species_ids = [str(item) for item in snakemake.params.species_ids]
proteomes = [Path(path) for path in snakemake.input.proteomes]
if len(species_ids) != len(proteomes):
    raise ValueError("species_ids and proteome paths have different lengths")

executable, version = executable_version(["orthofinder"], ["--version"])
signature_payload = {
    "species": [
        {
            "species_id": species,
            "proteome": str(source.resolve()),
            "sha256": sha256_file(source),
        }
        for species, source in zip(species_ids, proteomes, strict=True)
    ],
    "orthofinder_version": version,
    "search_threads": int(snakemake.params.search_threads),
    "analysis_threads": int(snakemake.params.analysis_threads),
    "extra_args": [str(value) for value in snakemake.params.extra_args],
    "preserve_stable_ids": True,
}
signature = hashlib.sha256(
    json.dumps(signature_payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
).hexdigest()

work_dir = Path(snakemake.params.work_dir)
run_root = work_dir / "runs" / signature[:16]
input_dir = run_root / "proteomes"
result_dir = run_root / "results"
state_path = run_root / "orthofinder.resume.json"
run_root.mkdir(parents=True, exist_ok=True)
input_dir.mkdir(parents=True, exist_ok=True)
for species, source in zip(species_ids, proteomes, strict=True):
    destination = input_dir / f"{species}.fa"
    if destination.exists() or destination.is_symlink():
        continue
    try:
        destination.symlink_to(source.resolve())
    except OSError:
        shutil.copyfile(source, destination)


def find_result() -> Path | None:
    hog_directories = sorted(result_dir.rglob("Phylogenetic_Hierarchical_Orthogroups"))
    return hog_directories[0].parent if hog_directories else None


actual_result_dir = find_result()
resumed_from: str | None = None
command: list[str]
if actual_result_dir is None:
    working_directories = (
        sorted(result_dir.rglob("WorkingDirectory")) if result_dir.exists() else []
    )
    if working_directories:
        resumed_from = str(working_directories[0].resolve())
        command = [
            executable,
            "-b",
            resumed_from,
            "-t",
            str(snakemake.params.search_threads),
            "-a",
            str(snakemake.params.analysis_threads),
            "-X",
        ]
    else:
        if result_dir.exists() and any(result_dir.iterdir()):
            timestamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")
            failed_dir = run_root / "failed" / timestamp
            failed_dir.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(result_dir), failed_dir)
        result_dir.mkdir(parents=True, exist_ok=True)
        command = [
            executable,
            "-f",
            str(input_dir),
            "-t",
            str(snakemake.params.search_threads),
            "-a",
            str(snakemake.params.analysis_threads),
            "-o",
            str(result_dir),
            # Canonical proteins already have globally unique species prefixes.
            "-X",
        ]
    command.extend(str(value) for value in snakemake.params.extra_args)
    write_json(
        {
            "status": "running",
            "signature": signature,
            "signature_payload": signature_payload,
            "command": command,
            "resumed_from": resumed_from,
            "updated_at_utc": datetime.now(UTC).isoformat(),
        },
        state_path,
    )
    run_command(
        command,
        stdout_path=snakemake.log.stdout,
        stderr_path=snakemake.log.stderr,
        cwd=run_root,
    )
    actual_result_dir = find_result()
else:
    command = ["reuse-complete-result", str(actual_result_dir.resolve())]

if actual_result_dir is None:
    working_directories = (
        sorted(result_dir.rglob("WorkingDirectory")) if result_dir.exists() else []
    )
    working_note = (
        f" WorkingDirectory retained at {working_directories[0]} for a retry."
        if working_directories
        else ""
    )
    raise FileNotFoundError(
        "OrthoFinder did not produce a Phylogenetic_Hierarchical_Orthogroups directory."
        + working_note
    )

marker = {
    "status": "complete",
    "signature": signature,
    "signature_payload": signature_payload,
    "orthofinder_executable": executable,
    "orthofinder_version": version,
    "command": command,
    "resumed_from": resumed_from,
    "result_dir": str(actual_result_dir.resolve()),
    "run_root": str(run_root.resolve()),
    "species_ids": species_ids,
    "updated_at_utc": datetime.now(UTC).isoformat(),
}
write_json(marker, state_path)
write_json(marker, snakemake.output.done)
write_text_atomic(str(actual_result_dir.resolve()) + "\n", snakemake.output.result_dir)
write_text_atomic(str(run_root.resolve()) + "\n", work_dir / "latest_run.txt")
