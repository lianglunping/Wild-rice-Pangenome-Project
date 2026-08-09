from __future__ import annotations

import json
from pathlib import Path

from panfamflow.config import load_config
from panfamflow.runner import (
    build_snakemake_command,
    finalize_run_provenance,
    parse_snakemake_summary,
    write_run_provenance,
)

TOY_CONFIG = Path(__file__).parents[1] / "examples" / "toy" / "config.yaml"


def test_build_snakemake_command_is_list_based() -> None:
    config = load_config(TOY_CONFIG)
    command, stack = build_snakemake_command(config, TOY_CONFIG, ("qc",), dry_run=True)
    with stack:
        assert command[0] == "snakemake"
        assert "--dry-run" in command
        assert any(item == "panfamflow_selected_modules=qc" for item in command)
        assert command[-1].endswith("results/00_qc/qc.done")


def test_smart_resume_flags_are_enabled_by_default() -> None:
    config = load_config(TOY_CONFIG)
    command, stack = build_snakemake_command(config, TOY_CONFIG, ("qc",))
    with stack:
        assert "--rerun-incomplete" in command
        assert "--keep-going" in command
        assert command[command.index("--retries") + 1] == "1"
        trigger_index = command.index("--rerun-triggers")
        assert command[trigger_index + 1 : trigger_index + 6] == [
            "mtime",
            "input",
            "params",
            "code",
            "software-env",
        ]
        assert "--printshellcmds" in command
        assert "--show-failed-logs" in command


def test_resume_mode_off_can_be_overridden_explicitly() -> None:
    config = load_config(TOY_CONFIG)
    disabled = config.model_copy(
        update={"run": config.run.model_copy(update={"resume_mode": "off"})}
    )
    command, stack = build_snakemake_command(disabled, TOY_CONFIG, ("qc",))
    with stack:
        assert "--rerun-incomplete" not in command
        assert "--rerun-triggers" not in command
    resumed, stack = build_snakemake_command(
        disabled,
        TOY_CONFIG,
        ("qc",),
        force_resume=True,
    )
    with stack:
        assert "--rerun-incomplete" in resumed
        assert "--rerun-triggers" in resumed


def test_mtime_only_resume_uses_single_trigger() -> None:
    config = load_config(TOY_CONFIG)
    mtime_only = config.model_copy(
        update={"run": config.run.model_copy(update={"resume_mode": "mtime_only"})}
    )
    command, stack = build_snakemake_command(mtime_only, TOY_CONFIG, ("qc",))
    with stack:
        index = command.index("--rerun-triggers")
        assert command[index + 1] == "mtime"
        assert command[index + 2] == "--"
        assert command[index + 3].startswith("results/")


def test_parse_snakemake_summary() -> None:
    text = (
        "Building DAG...\n"
        "output_file\tdate\trule\tversion\tlog-file(s)\tstatus\tplan\n"
        "results/a.tsv\tSat Aug 9\taudit\t-\tlogs/a.log\tok\tno update\n"
    )
    rows = parse_snakemake_summary(text)
    assert len(rows) == 1
    assert rows[0].output_file == "results/a.tsv"
    assert rows[0].status == "ok"
    assert rows[0].plan == "no update"


def test_run_provenance_is_written_atomically(tmp_path: Path) -> None:
    config = load_config(TOY_CONFIG)
    project = config.project.model_copy(update={"root": tmp_path})
    relocated = config.model_copy(update={"project": project})
    path = write_run_provenance(
        relocated,
        TOY_CONFIG,
        ("qc",),
        ["snakemake", "--summary"],
        mode="test",
    )
    assert path.is_file()
    record = json.loads(path.read_text(encoding="utf-8"))
    assert record["mode"] == "test"
    assert record["modules"] == ["qc"]
    assert record["status"] == "RUNNING"
    finalize_run_provenance(path, 0)
    finalized = json.loads(path.read_text(encoding="utf-8"))
    assert finalized["status"] == "COMPLETED"
    assert finalized["exit_code"] == 0
    assert finalized["finished_at_utc"]
    assert (tmp_path / ".panfamflow" / "provenance" / "resolved_config.yaml").is_file()
    assert not list(tmp_path.rglob("*.partial.*"))
