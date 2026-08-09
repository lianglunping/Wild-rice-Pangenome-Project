"""Command-line interface for PanFamFlow."""

from __future__ import annotations

import json
import shlex
import shutil
from importlib import resources
from pathlib import Path
from typing import Annotated

import typer
from pydantic import ValidationError
from rich.console import Console
from rich.table import Table

from panfamflow import __version__
from panfamflow.benchmark import (
    BenchmarkManifest,
    audit_benchmark,
    default_audit_output_dir,
    initialize_benchmark,
    write_benchmark_audit,
)
from panfamflow.config import (
    WorkflowConfig,
    load_config,
    project_root,
    validate_input_paths,
)
from panfamflow.modules import MODULES, resolve_modules
from panfamflow.runner import (
    build_snakemake_command,
    execute,
    execute_capture,
    finalize_run_provenance,
    parse_snakemake_summary,
    write_run_provenance,
)

app = typer.Typer(
    name="panfamflow",
    no_args_is_help=True,
    add_completion=False,
    help="Configuration-driven target pan-gene-family analysis workflow.",
)
benchmark_app = typer.Typer(
    name="benchmark",
    no_args_is_help=True,
    add_completion=False,
    help="Initialize and audit fail-closed biological benchmark inputs.",
)
app.add_typer(benchmark_app, name="benchmark")
console = Console()


def _version_callback(value: bool) -> None:
    if value:
        console.print(f"PanFamFlow {__version__}")
        raise typer.Exit()


@app.callback()
def main(
    version: Annotated[
        bool | None,
        typer.Option("--version", callback=_version_callback, is_eager=True, help="Show version."),
    ] = None,
) -> None:
    """PanFamFlow command group."""


def _load_or_exit(path: Path) -> WorkflowConfig:
    try:
        return load_config(path)
    except (FileNotFoundError, ValueError, ValidationError) as error:
        console.print(f"[bold red]Configuration error:[/bold red] {error}")
        raise typer.Exit(code=2) from error


def _selected(config: WorkflowConfig, modules: list[str] | None) -> tuple[str, ...]:
    try:
        return resolve_modules(modules or config.run.modules, config)
    except ValueError as error:
        console.print(f"[bold red]Module selection error:[/bold red] {error}")
        raise typer.Exit(code=2) from error


@app.command("init")
def init_project(
    destination: Annotated[Path, typer.Argument(help="New project directory.")],
    force: Annotated[
        bool,
        typer.Option(help="Allow an existing empty destination; never overwrite config.yaml."),
    ] = False,
) -> None:
    """Create a non-destructive project skeleton and config.yaml."""

    target = destination.expanduser().resolve()
    if target.exists() and not target.is_dir():
        console.print(f"[red]Destination is not a directory:[/red] {target}")
        raise typer.Exit(code=2)
    if target.exists() and any(target.iterdir()) and not force:
        console.print(f"[red]Destination is not empty:[/red] {target}")
        raise typer.Exit(code=2)
    target.mkdir(parents=True, exist_ok=True)
    config_target = target / "config.yaml"
    if config_target.exists():
        console.print(f"[red]Refusing to overwrite:[/red] {config_target}")
        raise typer.Exit(code=2)
    schema_target = target / ".panfamflow" / "config.schema.json"
    if schema_target.exists():
        console.print(f"[red]Refusing to overwrite:[/red] {schema_target}")
        raise typer.Exit(code=2)
    template = (
        resources.files("panfamflow.templates").joinpath("config.yaml").read_text(encoding="utf-8")
    )
    config_target.write_text(template, encoding="utf-8")
    schema_target.parent.mkdir(parents=True, exist_ok=True)
    schema_target.write_text(
        json.dumps(WorkflowConfig.model_json_schema(), ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    for directory in ("data", "references", "results", "work", "logs"):
        (target / directory).mkdir(exist_ok=True)
    console.print(f"Created project: [bold]{target}[/bold]")
    console.print("Edit config.yaml, then run: panfamflow validate -c config.yaml")


@app.command("list-modules")
def list_modules() -> None:
    """List selectable modules, direct dependencies and principal tools."""

    table = Table(title="PanFamFlow modules")
    table.add_column("Module", style="bold")
    table.add_column("Direct dependencies")
    table.add_column("Primary tools")
    table.add_column("Description")
    for spec in MODULES.values():
        table.add_row(
            spec.name,
            ", ".join(spec.direct_dependencies) or "—",
            ", ".join(spec.primary_tools),
            spec.description,
        )
    console.print(table)


@app.command("validate")
def validate(
    config_path: Annotated[
        Path, typer.Option("--config", "-c", help="Path to config.yaml.")
    ] = Path("config.yaml"),
    module: Annotated[
        list[str] | None,
        typer.Option("--module", "-m", help="Validate only this module closure; repeatable."),
    ] = None,
) -> None:
    """Validate schema, dependency closure and module-specific input paths."""

    config = _load_or_exit(config_path)
    modules = _selected(config, module)
    issues = validate_input_paths(config, config_path, modules)
    table = Table(title="Validation")
    table.add_column("Severity")
    table.add_column("Field")
    table.add_column("Message")
    for issue in issues:
        style = "red" if issue.severity == "ERROR" else "yellow"
        table.add_row(f"[{style}]{issue.severity}[/{style}]", issue.field, issue.message)
    if issues:
        console.print(table)
    errors = [issue for issue in issues if issue.severity == "ERROR"]
    console.print(f"Project root: {project_root(config, config_path)}")
    console.print(f"Resolved modules: {', '.join(modules)}")
    if errors:
        console.print(f"[bold red]Validation failed with {len(errors)} error(s).[/bold red]")
        raise typer.Exit(code=2)
    console.print("[bold green]Validation passed.[/bold green]")


@app.command("plan")
def plan(
    config_path: Annotated[
        Path, typer.Option("--config", "-c", help="Path to config.yaml.")
    ] = Path("config.yaml"),
    module: Annotated[
        list[str] | None,
        typer.Option("--module", "-m", help="Plan this module; repeatable."),
    ] = None,
    command_only: Annotated[
        bool, typer.Option(help="Print the exact Snakemake command after the module table.")
    ] = False,
) -> None:
    """Show dependency-expanded execution order and final targets."""

    config = _load_or_exit(config_path)
    modules = _selected(config, module)
    table = Table(title="Execution plan")
    table.add_column("Order", justify="right")
    table.add_column("Module", style="bold")
    table.add_column("Final target")
    table.add_column("Tools")
    for index, name in enumerate(modules, start=1):
        spec = MODULES[name]
        target = config.project.results_dir / spec.target
        table.add_row(str(index), name, str(target), ", ".join(spec.primary_tools))
    console.print(table)
    console.print(
        f"Resume mode: {config.run.resume_mode}; "
        f"rerun_incomplete={config.run.rerun_incomplete}; "
        f"keep_going={config.run.keep_going}; retries={config.run.retries}"
    )
    if command_only:
        command, stack = build_snakemake_command(config, config_path, modules, dry_run=True)
        with stack:
            console.print(shlex.join(command))


def _run_selected_workflow(
    *,
    config_path: Path,
    module: list[str] | None,
    cores: int | None,
    profile: Path | None,
    dry_run: bool,
    unlock: bool,
    summary: bool,
    dag: bool,
    skip_validation: bool,
    force_resume: bool,
    invocation_mode: str,
) -> None:
    config = _load_or_exit(config_path)
    modules = _selected(config, module)
    if not skip_validation and not (unlock or summary or dag):
        issues = validate_input_paths(config, config_path, modules)
        errors = [issue for issue in issues if issue.severity == "ERROR"]
        for issue in issues:
            style = "red" if issue.severity == "ERROR" else "yellow"
            console.print(f"[{style}]{issue.severity}[/{style}] {issue.field}: {issue.message}")
        if errors:
            raise typer.Exit(code=2)
    try:
        command, stack = build_snakemake_command(
            config,
            config_path,
            modules,
            cores=cores,
            profile=profile,
            dry_run=dry_run,
            unlock=unlock,
            summary=summary,
            dag=dag,
            force_resume=force_resume,
        )
    except RuntimeError as error:
        console.print(f"[bold red]Engine error:[/bold red] {error}")
        raise typer.Exit(code=2) from error
    console.print(f"Resolved modules: {', '.join(modules)}")
    console.print(f"Command: {shlex.join(command)}")
    provenance: Path | None = None
    with stack:
        if not (unlock or summary or dag):
            provenance = write_run_provenance(
                config,
                config_path,
                modules,
                command,
                mode=("dry-run" if dry_run else invocation_mode),
            )
            console.print(f"Provenance: {provenance}")
        try:
            return_code = execute(command)
        except OSError as error:
            if provenance is not None:
                finalize_run_provenance(
                    provenance,
                    127,
                    error=f"{type(error).__name__}: {error}",
                )
            console.print(f"[bold red]Launch error:[/bold red] {error}")
            raise typer.Exit(code=127) from error
    if provenance is not None:
        finalize_run_provenance(provenance, return_code)
    if return_code != 0:
        raise typer.Exit(code=return_code)


@app.command("run")
def run(
    config_path: Annotated[
        Path, typer.Option("--config", "-c", help="Path to config.yaml.")
    ] = Path("config.yaml"),
    module: Annotated[
        list[str] | None,
        typer.Option("--module", "-m", help="Run this module; repeatable."),
    ] = None,
    cores: Annotated[int | None, typer.Option(min=1, help="Override run.cores.")] = None,
    profile: Annotated[Path | None, typer.Option(help="Override run.profile.")] = None,
    dry_run: Annotated[bool, typer.Option(help="Build the DAG without executing jobs.")] = False,
    unlock: Annotated[bool, typer.Option(help="Remove a confirmed stale Snakemake lock.")] = False,
    summary: Annotated[bool, typer.Option(help="Print Snakemake output-file status.")] = False,
    dag: Annotated[bool, typer.Option(help="Emit the DAG in Graphviz DOT format.")] = False,
    skip_validation: Annotated[
        bool, typer.Option(help="Skip input-path checks; schema checks still apply.")
    ] = False,
) -> None:
    """Run selected modules; completed and still-valid jobs are skipped automatically."""

    _run_selected_workflow(
        config_path=config_path,
        module=module,
        cores=cores,
        profile=profile,
        dry_run=dry_run,
        unlock=unlock,
        summary=summary,
        dag=dag,
        skip_validation=skip_validation,
        force_resume=False,
        invocation_mode="run",
    )


@app.command("resume")
def resume(
    config_path: Annotated[
        Path, typer.Option("--config", "-c", help="Path to config.yaml.")
    ] = Path("config.yaml"),
    module: Annotated[
        list[str] | None,
        typer.Option("--module", "-m", help="Resume this module closure; repeatable."),
    ] = None,
    cores: Annotated[int | None, typer.Option(min=1, help="Override run.cores.")] = None,
    profile: Annotated[Path | None, typer.Option(help="Override run.profile.")] = None,
    dry_run: Annotated[bool, typer.Option(help="Show jobs that would be resumed.")] = False,
    skip_validation: Annotated[
        bool, typer.Option(help="Skip input-path checks; schema checks still apply.")
    ] = False,
) -> None:
    """Force smart resume even when ``run.resume_mode`` is configured as ``off``."""

    _run_selected_workflow(
        config_path=config_path,
        module=module,
        cores=cores,
        profile=profile,
        dry_run=dry_run,
        unlock=False,
        summary=False,
        dag=False,
        skip_validation=skip_validation,
        force_resume=True,
        invocation_mode="resume",
    )


@app.command("retry")
def retry(
    config_path: Annotated[
        Path, typer.Option("--config", "-c", help="Path to config.yaml.")
    ] = Path("config.yaml"),
    module: Annotated[
        list[str] | None,
        typer.Option("--module", "-m", help="Retry failed/incomplete module closure; repeatable."),
    ] = None,
    cores: Annotated[int | None, typer.Option(min=1, help="Override run.cores.")] = None,
    profile: Annotated[Path | None, typer.Option(help="Override run.profile.")] = None,
    dry_run: Annotated[bool, typer.Option(help="Show jobs that would be retried.")] = False,
) -> None:
    """Alias of resume, emphasizing failed or incomplete jobs."""

    _run_selected_workflow(
        config_path=config_path,
        module=module,
        cores=cores,
        profile=profile,
        dry_run=dry_run,
        unlock=False,
        summary=False,
        dag=False,
        skip_validation=False,
        force_resume=True,
        invocation_mode="retry",
    )


@app.command("status")
def status(
    config_path: Annotated[
        Path, typer.Option("--config", "-c", help="Path to config.yaml.")
    ] = Path("config.yaml"),
    module: Annotated[
        list[str] | None,
        typer.Option("--module", "-m", help="Inspect this module closure; repeatable."),
    ] = None,
) -> None:
    """Show Snakemake output status without modifying workflow state."""

    config = _load_or_exit(config_path)
    modules = _selected(config, module)
    try:
        command, stack = build_snakemake_command(
            config,
            config_path,
            modules,
            summary=True,
        )
    except RuntimeError as error:
        console.print(f"[bold red]Engine error:[/bold red] {error}")
        raise typer.Exit(code=2) from error
    with stack:
        completed = execute_capture(command)
    if completed.returncode != 0:
        console.print(completed.stderr or completed.stdout)
        raise typer.Exit(code=completed.returncode)
    rows = parse_snakemake_summary(completed.stdout)
    if not rows:
        console.print(completed.stdout.rstrip())
        return
    table = Table(title="PanFamFlow workflow status")
    table.add_column("Output")
    table.add_column("Rule")
    table.add_column("Status")
    table.add_column("Plan")
    table.add_column("Date")
    for row in rows:
        table.add_row(row.output_file, row.rule, row.status, row.plan, row.date)
    console.print(table)


@app.command("schema")
def schema(
    output: Annotated[Path, typer.Option("--output", "-o", help="Output JSON Schema path.")] = Path(
        "panfamflow-config.schema.json"
    ),
) -> None:
    """Export the machine-readable configuration JSON Schema."""

    target = output.expanduser().resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(
        json.dumps(WorkflowConfig.model_json_schema(), ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    console.print(f"Wrote schema: {target}")


@app.command("doctor")
def doctor() -> None:
    """Check launcher-level executables without assuming rule environments are active."""

    table = Table(title="PanFamFlow launcher doctor")
    table.add_column("Executable")
    table.add_column("Path")
    table.add_column("Status")
    for executable in ("uv", "mamba", "conda", "snakemake"):
        path = shutil.which(executable)
        table.add_row(executable, path or "—", "FOUND" if path else "NOT FOUND")
    console.print(table)
    console.print(
        "Rule-specific tools are installed lazily by Snakemake Conda environments and are not "
        "expected in the launcher PATH."
    )


@benchmark_app.command("init")
def benchmark_init(
    destination: Annotated[Path, typer.Argument(help="New benchmark workspace directory.")],
) -> None:
    """Create a non-destructive biological benchmark intake workspace."""

    try:
        written = initialize_benchmark(destination)
    except (FileExistsError, NotADirectoryError, OSError) as error:
        console.print(f"[bold red]Benchmark initialization error:[/bold red] {error}")
        raise typer.Exit(code=2) from error
    target = destination.expanduser().resolve()
    console.print(f"Created benchmark workspace: [bold]{target}[/bold]")
    for path in written:
        console.print(f"  - {path.relative_to(target)}")
    console.print(
        "Next: freeze benchmark.yaml/species.tsv/manual truth set, then run "
        "panfamflow benchmark audit."
    )


@benchmark_app.command("audit")
def benchmark_audit(
    manifest_path: Annotated[
        Path,
        typer.Option("--manifest", "-m", help="Path to benchmark.yaml."),
    ] = Path("benchmark.yaml"),
    output_dir: Annotated[
        Path | None,
        typer.Option("--output", "-o", help="New audit output directory."),
    ] = None,
    allow_blocked: Annotated[
        bool,
        typer.Option(
            help="Write planning outputs and return success even when the gate is BLOCKED."
        ),
    ] = False,
) -> None:
    """Audit biological benchmark readiness and emit Chinese HTML plus machine outputs."""

    try:
        audit = audit_benchmark(manifest_path)
        target = output_dir or default_audit_output_dir(manifest_path)
        paths = write_benchmark_audit(audit, target)
    except (FileNotFoundError, FileExistsError, ValueError, ValidationError, OSError) as error:
        console.print(f"[bold red]Benchmark audit error:[/bold red] {error}")
        raise typer.Exit(code=2) from error

    table = Table(title="PanFamFlow biological benchmark gate")
    table.add_column("Metric")
    table.add_column("Value")
    table.add_row("Overall status", audit.overall_status)
    table.add_row("Blocking failures", str(audit.blocking_failures))
    table.add_row("Warnings", str(audit.warnings))
    table.add_row("Passed checks", str(audit.passed))
    table.add_row("Manifest SHA256", audit.manifest_sha256)
    console.print(table)
    console.print(f"Chinese HTML: {paths['html']}")
    console.print(f"Machine JSON: {paths['json']}")
    console.print(f"TSV/XLSX: {paths['checks_tsv']} | {paths['xlsx']}")
    if audit.overall_status == "BLOCKED" and not allow_blocked:
        raise typer.Exit(code=2)


@benchmark_app.command("schema")
def benchmark_schema(
    output: Annotated[
        Path,
        typer.Option("--output", "-o", help="Output benchmark JSON Schema path."),
    ] = Path("panfamflow-benchmark.schema.json"),
) -> None:
    """Export the strict biological benchmark manifest JSON Schema."""

    target = output.expanduser().resolve()
    if target.exists():
        console.print(f"[red]Refusing to overwrite:[/red] {target}")
        raise typer.Exit(code=2)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(
        json.dumps(BenchmarkManifest.model_json_schema(), ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    console.print(f"Wrote benchmark schema: {target}")


if __name__ == "__main__":
    app()
