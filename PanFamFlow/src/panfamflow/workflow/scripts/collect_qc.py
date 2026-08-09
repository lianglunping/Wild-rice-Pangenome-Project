import sys
from pathlib import Path as _ScriptPath

sys.path.insert(0, str(_ScriptPath(snakemake.scriptdir)))

from datetime import UTC, datetime
from pathlib import Path

import pandas as pd
from workflow_utils import save_table, write_json

busco_paths = [Path(path) for path in snakemake.input.busco]
frames = [
    pd.read_csv(path, sep="\t") for path in busco_paths if path.exists() and path.stat().st_size
]
if frames:
    busco = pd.concat(frames, ignore_index=True)
else:
    busco = pd.DataFrame(
        columns=[
            "species_id",
            "lineage",
            "complete_pct",
            "single_copy_pct",
            "duplicated_pct",
            "fragmented_pct",
            "missing_pct",
            "busco_n",
            "summary_path",
        ]
    )
save_table(busco, snakemake.output.busco_tsv, snakemake.output.busco_xlsx)
write_json(
    {
        "status": "complete",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "input_audit": str(Path(snakemake.input.audit).resolve()),
        "busco_enabled": bool(snakemake.params.busco_enabled),
        "busco_species": int(busco.shape[0]),
    },
    snakemake.output.done,
)
