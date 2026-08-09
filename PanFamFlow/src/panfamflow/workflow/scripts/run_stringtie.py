import sys
from pathlib import Path as _ScriptPath

sys.path.insert(0, str(_ScriptPath(snakemake.scriptdir)))

from pathlib import Path

from workflow_utils import run_command

Path(snakemake.output.gtf).parent.mkdir(parents=True, exist_ok=True)
Path(snakemake.output.abundance).parent.mkdir(parents=True, exist_ok=True)
command = [
    "stringtie",
    str(snakemake.input.bam),
    "-G",
    str(snakemake.input.gff3),
    "-e",
    "-o",
    str(snakemake.output.gtf),
    "-A",
    str(snakemake.output.abundance),
    "-p",
    str(snakemake.threads),
]
strandedness = str(snakemake.params.strandedness)
if strandedness == "forward":
    command.append("--fr")
elif strandedness == "reverse":
    command.append("--rf")
command.extend(str(value) for value in snakemake.params.extra_args)
run_command(command, stdout_path=snakemake.log.stdout, stderr_path=snakemake.log.stderr)
