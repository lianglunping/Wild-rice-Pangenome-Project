import sys
from pathlib import Path as _ScriptPath

sys.path.insert(0, str(_ScriptPath(snakemake.scriptdir)))

from collections.abc import Iterator

from workflow_utils import iter_fasta_records, write_fasta


def combined_records() -> Iterator[tuple[str, str]]:
    seen: set[str] = set()
    for path in snakemake.input:
        for identifier, sequence in iter_fasta_records(path):
            if identifier in seen:
                raise ValueError(f"Duplicate stable FASTA ID across species: {identifier}")
            seen.add(identifier)
            yield identifier, sequence


write_fasta(combined_records(), snakemake.output[0])
