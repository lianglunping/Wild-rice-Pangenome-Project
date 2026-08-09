import sys
from pathlib import Path as _ScriptPath

sys.path.insert(0, str(_ScriptPath(snakemake.scriptdir)))

from datetime import UTC, datetime

from workflow_utils import write_json

write_json(
    {
        "status": "complete",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "species_count": len(snakemake.input.maps),
        "mapping_files": [str(path) for path in snakemake.input.maps],
    },
    snakemake.output.done,
)
