import sys
from pathlib import Path as _ScriptPath

sys.path.insert(0, str(_ScriptPath(snakemake.scriptdir)))

import json
import shutil
from datetime import UTC, datetime
from pathlib import Path

from workflow_utils import (
    copy_atomic,
    executable_version,
    read_fasta,
    run_command,
    sha256_file,
    write_json,
    write_text_atomic,
)

sequences = read_fasta(snakemake.input.proteins)
minimum = int(snakemake.params.min_sequences)
if len(sequences) < minimum:
    raise RuntimeError(
        f"Phylogeny requires at least {minimum} family sequences; observed {len(sequences)}."
    )

work_dir = Path(snakemake.params.work_dir)
work_dir.mkdir(parents=True, exist_ok=True)
alignment = Path(snakemake.output.alignment)
trimmed = Path(snakemake.output.trimmed)
prefix = work_dir / "family"
state_path = work_dir / "family.resume.json"

mode = str(snakemake.params.mafft_mode)
if mode == "linsi":
    mafft_command = ["mafft", "--localpair", "--maxiterate", "1000"]
elif mode == "ginsi":
    mafft_command = ["mafft", "--globalpair", "--maxiterate", "1000"]
elif mode == "einsi":
    mafft_command = ["mafft", "--genafpair", "--maxiterate", "1000"]
else:
    mafft_command = ["mafft", "--auto"]
mafft_command.extend(["--thread", str(snakemake.threads), str(snakemake.input.proteins)])
run_command(mafft_command, stdout_path=alignment, stderr_path=snakemake.log.mafft)

run_command(
    [
        "clipkit",
        str(alignment),
        "-m",
        str(snakemake.params.trim_mode),
        "-o",
        str(trimmed),
    ],
    stdout_path=snakemake.log.clipkit_stdout,
    stderr_path=snakemake.log.clipkit_stderr,
)

bootstrap = int(snakemake.params.ultrafast_bootstrap)
sh_alrt = int(snakemake.params.sh_alrt)
resume_signature = {
    "trimmed_alignment_sha256": sha256_file(trimmed),
    "model": str(snakemake.params.model),
    "ultrafast_bootstrap": bootstrap,
    "sh_alrt": sh_alrt,
    "seed": int(snakemake.params.seed),
}
previous_signature: dict[str, object] | None = None
if state_path.is_file():
    try:
        previous_signature = json.loads(state_path.read_text(encoding="utf-8")).get("signature")
    except (json.JSONDecodeError, OSError):
        previous_signature = None


def archive_stale_prefix() -> None:
    stale_files = sorted(work_dir.glob("family.*"))
    if not stale_files:
        return
    timestamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")
    stale_dir = work_dir / "stale" / timestamp
    stale_dir.mkdir(parents=True, exist_ok=True)
    for source in stale_files:
        if source == state_path:
            continue
        shutil.move(str(source), stale_dir / source.name)


if previous_signature != resume_signature:
    archive_stale_prefix()
    if state_path.exists():
        state_path.unlink()

write_json(
    {
        "status": "prepared",
        "signature": resume_signature,
        "updated_at_utc": datetime.now(UTC).isoformat(),
    },
    state_path,
)

iqtree_executable, version = executable_version(["iqtree3", "iqtree2", "iqtree"], ["--version"])
command = [
    iqtree_executable,
    "-s",
    str(trimmed),
    "-m",
    str(snakemake.params.model),
    "-T",
    str(snakemake.threads),
    "--seed",
    str(snakemake.params.seed),
    "--prefix",
    str(prefix),
]
if bootstrap:
    command.extend(["-B", str(bootstrap)])
if sh_alrt:
    command.extend(["--alrt", str(sh_alrt)])

tree_source = Path(f"{prefix}.treefile")
report_source = Path(f"{prefix}.iqtree")
checkpoint = Path(f"{prefix}.ckp.gz")
if not tree_source.is_file():
    # IQ-TREE reads a compatible .ckp.gz automatically when the command,
    # alignment and prefix are unchanged.  Do not use --redo here.
    run_command(
        command,
        stdout_path=snakemake.log.iqtree_stdout,
        stderr_path=snakemake.log.iqtree_stderr,
    )
if not tree_source.is_file():
    checkpoint_note = f"; checkpoint retained at {checkpoint}" if checkpoint.exists() else ""
    raise FileNotFoundError(
        f"IQ-TREE did not produce {tree_source}; executable was {version}{checkpoint_note}"
    )

copy_atomic(tree_source, snakemake.output.tree)
if report_source.is_file():
    copy_atomic(report_source, snakemake.output.report)
else:
    write_text_atomic(
        f"IQ-TREE executable: {iqtree_executable}\nVersion: {version}\n",
        snakemake.output.report,
    )
write_json(
    {
        "status": "complete",
        "signature": resume_signature,
        "treefile": str(tree_source.resolve()),
        "checkpoint_present": checkpoint.exists(),
        "updated_at_utc": datetime.now(UTC).isoformat(),
    },
    state_path,
)
