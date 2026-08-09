import sys
from pathlib import Path as _ScriptPath

sys.path.insert(0, str(_ScriptPath(snakemake.scriptdir)))

from pathlib import Path

from workflow_utils import run_command

paired = bool(snakemake.params.paired)
for output_path in (
    snakemake.output.r1,
    snakemake.output.r2,
    snakemake.output.json,
    snakemake.output.html,
):
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
command = [
    "fastp",
    "--in1",
    str(snakemake.input.r1),
    "--out1",
    str(snakemake.output.r1),
    "--thread",
    str(snakemake.threads),
    "--json",
    str(snakemake.output.json),
    "--html",
    str(snakemake.output.html),
]
if paired:
    command.extend(["--in2", str(snakemake.input.r2), "--out2", str(snakemake.output.r2)])
command.extend(str(value) for value in snakemake.params.extra_args)
run_command(command, stdout_path=snakemake.log.stdout, stderr_path=snakemake.log.stderr)
if not paired:
    Path(snakemake.output.r2).parent.mkdir(parents=True, exist_ok=True)
    Path(snakemake.output.r2).write_bytes(b"")
