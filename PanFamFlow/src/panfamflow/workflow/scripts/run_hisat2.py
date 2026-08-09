import subprocess
from pathlib import Path

paired = bool(snakemake.params.paired)
command = [
    "hisat2",
    "-x",
    str(snakemake.params.index_prefix),
    "-p",
    str(snakemake.threads),
]
if paired:
    command.extend(["-1", str(snakemake.input.r1), "-2", str(snakemake.input.r2)])
else:
    command.extend(["-U", str(snakemake.input.r1)])
command.extend(str(value) for value in snakemake.params.extra_args)
Path(snakemake.output.bam).parent.mkdir(parents=True, exist_ok=True)
Path(snakemake.log.hisat2).parent.mkdir(parents=True, exist_ok=True)
Path(snakemake.log.samtools).parent.mkdir(parents=True, exist_ok=True)
with (
    Path(snakemake.log.hisat2).open("w", encoding="utf-8") as hisat2_log,
    Path(snakemake.log.samtools).open("w", encoding="utf-8") as samtools_log,
):
    mapper = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=hisat2_log)
    if mapper.stdout is None:
        raise RuntimeError("Failed to open HISAT2 stdout pipe")
    sorter = subprocess.run(
        ["samtools", "sort", "-@", str(snakemake.threads), "-o", str(snakemake.output.bam), "-"],
        stdin=mapper.stdout,
        stderr=samtools_log,
        check=False,
    )
    mapper.stdout.close()
    mapper_return = mapper.wait()
if mapper_return != 0:
    raise RuntimeError(f"HISAT2 failed with exit code {mapper_return}")
if sorter.returncode != 0:
    raise RuntimeError(f"samtools sort failed with exit code {sorter.returncode}")
subprocess.run(
    ["samtools", "index", str(snakemake.output.bam), str(snakemake.output.bai)], check=True
)
