"""Build, inspect and execute Snakemake commands without shell interpolation."""

from __future__ import annotations

import csv
import io
import json
import os
import shutil
import subprocess
from collections.abc import Sequence
from contextlib import ExitStack
from dataclasses import dataclass
from datetime import UTC, datetime
from importlib import resources
from pathlib import Path
from typing import Any

import yaml

from panfamflow.config import (
    WorkflowConfig,
    analysis_config_hash,
    execution_config_hash,
    project_root,
    resolve_project_path,
)
from panfamflow.modules import targets_for_modules


@dataclass(frozen=True, slots=True)
class SummaryRow:
    """One normalized row from ``snakemake --summary``."""

    output_file: str
    date: str
    rule: str
    log_files: str
    status: str
    plan: str


def _engine_prefix(config: WorkflowConfig) -> list[str]:
    runner = config.run.engine_runner
    environment = config.run.engine_env
    if runner == "current" or environment is None:
        return []
    candidates = [runner] if runner != "auto" else ["mamba", "conda"]
    executable = next((name for name in candidates if shutil.which(name)), None)
    if executable is None:
        expected = " or ".join(candidates)
        raise RuntimeError(f"Cannot find {expected} in PATH.")
    return [executable, "run", "-n", environment]


def _resume_arguments(config: WorkflowConfig, *, force_resume: bool) -> list[str]:
    mode = config.run.resume_mode
    if mode == "off" and not force_resume:
        return []

    arguments: list[str] = []
    if config.run.rerun_incomplete:
        arguments.append("--rerun-incomplete")
    if config.run.keep_going:
        arguments.append("--keep-going")
    arguments.extend(["--retries", str(config.run.retries)])
    triggers = ["mtime"] if mode == "mtime_only" and not force_resume else config.run.rerun_triggers
    arguments.append("--rerun-triggers")
    arguments.extend(triggers)
    return arguments


def build_snakemake_command(
    config: WorkflowConfig,
    config_path: Path,
    modules: Sequence[str],
    *,
    cores: int | None = None,
    profile: Path | None = None,
    dry_run: bool = False,
    unlock: bool = False,
    summary: bool = False,
    dag: bool = False,
    force_resume: bool = False,
) -> tuple[list[str], ExitStack]:
    """Return a subprocess argument list and an ExitStack owning package resources."""

    stack = ExitStack()
    snakefile_resource = resources.files("panfamflow.workflow").joinpath("Snakefile")
    snakefile = stack.enter_context(resources.as_file(snakefile_resource))
    resolved_config = config_path.expanduser().resolve()
    root = project_root(config, config_path)
    targets = list(targets_for_modules(modules, str(config.project.results_dir)))

    command = [*_engine_prefix(config), "snakemake"]
    command.extend(["--snakefile", str(snakefile.resolve())])
    command.extend(["--configfile", str(resolved_config)])
    command.extend(["--directory", str(root)])
    command.extend(["--cores", str(cores or config.run.cores)])
    command.extend(["--jobs", str(config.run.jobs)])
    command.extend(["--latency-wait", str(config.run.latency_wait)])
    command.extend(
        [
            "--config",
            f"panfamflow_selected_modules={','.join(modules)}",
            f"panfamflow_config_path={resolved_config}",
        ]
    )

    selected_profile = profile or config.run.profile
    if selected_profile is not None:
        resolved_profile = resolve_project_path(selected_profile, root)
        assert resolved_profile is not None
        command.extend(["--profile", str(resolved_profile)])
    if config.run.use_conda:
        command.extend(["--software-deployment-method", "conda"])
    if config.run.printshellcmds:
        command.append("--printshellcmds")
    if config.run.show_failed_logs:
        command.append("--show-failed-logs")
    command.extend(_resume_arguments(config, force_resume=force_resume))

    if dry_run:
        command.append("--dry-run")
    if unlock:
        command.append("--unlock")
    if summary:
        command.append("--summary")
    if dag:
        command.append("--dag")
    command.extend(config.run.extra_snakemake_args)
    # ``--rerun-triggers`` accepts multiple values. Without the standard
    # option terminator, argparse can consume the first workflow target as
    # another trigger. Delimit all positional targets explicitly.
    command.append("--")
    command.extend(targets)
    return command, stack


def _atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.partial.{os.getpid()}")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)


def write_run_provenance(
    config: WorkflowConfig,
    config_path: Path,
    modules: Sequence[str],
    command: Sequence[str],
    *,
    mode: str,
) -> Path:
    """Persist resolved configuration and immutable run metadata atomically."""

    root = project_root(config, config_path)
    provenance_root = root / ".panfamflow" / "provenance"
    runs_dir = provenance_root / "runs"
    timestamp = datetime.now(UTC)
    run_id = timestamp.strftime("%Y%m%dT%H%M%S.%fZ")
    resolved = config.model_dump(mode="json", exclude_none=False)
    analysis_hash = analysis_config_hash(config)
    execution_hash = execution_config_hash(config)

    _atomic_write_text(
        provenance_root / "resolved_config.yaml",
        yaml.safe_dump(resolved, sort_keys=False, allow_unicode=True),
    )
    fingerprints: dict[str, Any] = {
        "analysis_config_sha256": analysis_hash,
        "execution_config_sha256": execution_hash,
        "config_path": str(config_path.expanduser().resolve()),
        "updated_at_utc": timestamp.isoformat(),
    }
    _atomic_write_text(
        provenance_root / "fingerprints.json",
        json.dumps(fingerprints, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
    )
    run_record: dict[str, Any] = {
        **fingerprints,
        "run_id": run_id,
        "status": "RUNNING",
        "started_at_utc": timestamp.isoformat(),
        "finished_at_utc": None,
        "exit_code": None,
        "mode": mode,
        "modules": list(modules),
        "command": list(command),
        "resume_mode": config.run.resume_mode,
        "rerun_incomplete": config.run.rerun_incomplete,
        "keep_going": config.run.keep_going,
        "retries": config.run.retries,
        "rerun_triggers": list(config.run.rerun_triggers),
    }
    run_path = runs_dir / f"{run_id}.json"
    _atomic_write_text(
        run_path,
        json.dumps(run_record, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
    )
    return run_path


def finalize_run_provenance(
    run_path: Path,
    exit_code: int,
    *,
    error: str | None = None,
) -> None:
    """Atomically mark a run record as completed or failed.

    The launch record is written before Snakemake starts.  This function makes
    the provenance truthful after the subprocess exits, including launch
    failures that occur before Snakemake can create its own logs.
    """

    try:
        record = json.loads(run_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"Cannot finalize provenance record {run_path}: {exc}") from exc
    finished = datetime.now(UTC)
    record["finished_at_utc"] = finished.isoformat()
    record["exit_code"] = int(exit_code)
    record["status"] = "COMPLETED" if exit_code == 0 else "FAILED"
    if error is not None:
        record["error"] = error
    _atomic_write_text(
        run_path,
        json.dumps(record, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
    )


def execute(command: Sequence[str]) -> int:
    """Execute a command and return its exit code."""

    completed = subprocess.run(list(command), check=False)
    return completed.returncode


def execute_capture(command: Sequence[str]) -> subprocess.CompletedProcess[str]:
    """Execute a command while capturing UTF-8 text output."""

    return subprocess.run(
        list(command),
        check=False,
        capture_output=True,
        text=True,
    )


def parse_snakemake_summary(text: str) -> tuple[SummaryRow, ...]:
    """Parse Snakemake's tab-separated summary while tolerating preamble lines."""

    lines = [line for line in text.splitlines() if line.strip()]
    header_index = next(
        (
            index
            for index, line in enumerate(lines)
            if line.lower().startswith("output_file\t") or line.lower().startswith("output file\t")
        ),
        None,
    )
    if header_index is None:
        return ()
    reader = csv.DictReader(io.StringIO("\n".join(lines[header_index:])), delimiter="\t")
    rows: list[SummaryRow] = []
    for row in reader:
        normalized = {
            str(key).strip().lower().replace(" ", "_"): value or "" for key, value in row.items()
        }
        rows.append(
            SummaryRow(
                output_file=normalized.get("output_file", ""),
                date=normalized.get("date", ""),
                rule=normalized.get("rule", ""),
                log_files=normalized.get("log-file(s)", normalized.get("log_files", "")),
                status=normalized.get("status", ""),
                plan=normalized.get("plan", ""),
            )
        )
    return tuple(rows)
