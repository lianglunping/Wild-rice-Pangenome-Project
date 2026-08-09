"""Small, dependency-light helpers shared by Snakemake scripts."""

import gzip
import hashlib
import json
import os
import subprocess
import uuid
from collections import OrderedDict
from collections.abc import Iterable, Iterator, Mapping, Sequence
from pathlib import Path
from typing import Any, TextIO

import pandas as pd


def partial_path(path: str | Path) -> Path:
    """Return a unique sibling path that is never mistaken for a final output."""

    target = Path(path)
    return target.with_name(f".{target.name}.partial.{os.getpid()}.{uuid.uuid4().hex}")


def commit_partial(temporary: str | Path, target: str | Path) -> None:
    """Atomically publish a completed temporary file on the same filesystem."""

    source = Path(temporary)
    destination = Path(target)
    destination.parent.mkdir(parents=True, exist_ok=True)
    os.replace(source, destination)


def write_text_atomic(text: str, path: str | Path) -> None:
    """Write text through a uniquely named partial file and atomic rename."""

    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = partial_path(target)
    temporary.write_text(text, encoding="utf-8")
    commit_partial(temporary, target)


def copy_atomic(source: str | Path, target: str | Path) -> None:
    """Copy a file without exposing a partially written destination."""

    import shutil

    destination = Path(target)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = partial_path(destination)
    shutil.copy2(source, temporary)
    commit_partial(temporary, destination)


def ensure_nonempty(path: str | Path) -> None:
    """Reject missing or zero-byte outputs before completion metadata is written."""

    target = Path(path)
    if not target.is_file() or target.stat().st_size == 0:
        raise RuntimeError(f"Expected non-empty output was not produced: {target}")


def open_text(path: str | Path, mode: str = "rt") -> TextIO:
    file_path = Path(path)
    if file_path.suffix == ".gz":
        return gzip.open(file_path, mode, encoding="utf-8")  # type: ignore[return-value]
    return file_path.open(mode, encoding="utf-8")


def sha256_file(path: str | Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def iter_fasta_records(path: str | Path) -> Iterator[tuple[str, str]]:
    """Yield FASTA records one at a time while rejecting duplicate identifiers."""

    seen: set[str] = set()
    identifier: str | None = None
    chunks: list[str] = []
    with open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    yield identifier, "".join(chunks)
                identifier = line[1:].split()[0]
                if identifier in seen:
                    raise ValueError(f"Duplicate FASTA identifier {identifier!r} in {path}")
                seen.add(identifier)
                chunks = []
            else:
                if identifier is None:
                    raise ValueError(f"FASTA sequence encountered before header in {path}")
                chunks.append(line)
    if identifier is not None:
        yield identifier, "".join(chunks)


def read_fasta(path: str | Path) -> OrderedDict[str, str]:
    sequences: OrderedDict[str, str] = OrderedDict()
    for identifier, sequence in iter_fasta_records(path):
        sequences[identifier] = sequence
    return sequences


def fasta_lengths(path: str | Path) -> OrderedDict[str, int]:
    """Read only FASTA identifiers and sequence lengths without storing sequences."""

    lengths: OrderedDict[str, int] = OrderedDict()
    identifier: str | None = None
    sequence_length = 0
    with open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    lengths[identifier] = sequence_length
                identifier = line[1:].split()[0]
                if identifier in lengths:
                    raise ValueError(f"Duplicate FASTA identifier {identifier!r} in {path}")
                sequence_length = 0
            else:
                if identifier is None:
                    raise ValueError(f"FASTA sequence encountered before header in {path}")
                sequence_length += len(line)
    if identifier is not None:
        lengths[identifier] = sequence_length
    return lengths


def write_fasta(
    records: Mapping[str, str] | Iterable[tuple[str, str]],
    path: str | Path,
    width: int = 60,
) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = partial_path(target)
    items = records.items() if isinstance(records, Mapping) else records
    with temporary.open("w", encoding="utf-8") as handle:
        for identifier, sequence in items:
            handle.write(f">{identifier}\n")
            for start in range(0, len(sequence), width):
                handle.write(sequence[start : start + width] + "\n")
    commit_partial(temporary, target)


def parse_gff_attributes(text: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for item in text.strip().strip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        elif " " in item:
            key, value = item.split(" ", 1)
            value = value.strip().strip('"')
        else:
            continue
        attributes[key.strip()] = value.strip()
    return attributes


def first_parent(value: str | None) -> str | None:
    if value is None:
        return None
    return value.split(",", 1)[0]


def iter_gff(path: str | Path) -> Iterator[dict[str, Any]]:
    with open_text(path) as handle:
        for line_number, raw in enumerate(handle, start=1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"Expected 9 GFF/GTF columns at {path}:{line_number}")
            seqid, source, feature, start, end, score, strand, phase, attributes = fields
            yield {
                "seqid": seqid,
                "source": source,
                "feature": feature,
                "start": int(start),
                "end": int(end),
                "score": score,
                "strand": strand,
                "phase": phase,
                "attributes": parse_gff_attributes(attributes),
                "raw": raw,
            }


def save_table(df: pd.DataFrame, tsv: str | Path, xlsx: str | Path | None = None) -> None:
    tsv_path = Path(tsv)
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    tsv_temporary = partial_path(tsv_path)
    xlsx_path = Path(xlsx) if xlsx is not None else None
    xlsx_temporary = partial_path(xlsx_path) if xlsx_path is not None else None
    df.to_csv(tsv_temporary, sep="\t", index=False, na_rep="NA")
    if xlsx_path is not None and xlsx_temporary is not None:
        xlsx_path.parent.mkdir(parents=True, exist_ok=True)
        with pd.ExcelWriter(xlsx_temporary, engine="openpyxl") as writer:
            df.to_excel(writer, sheet_name="data", index=False)
    commit_partial(tsv_temporary, tsv_path)
    if xlsx_path is not None and xlsx_temporary is not None:
        commit_partial(xlsx_temporary, xlsx_path)


def save_workbook(tables: Mapping[str, pd.DataFrame], path: str | Path) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = partial_path(target)
    with pd.ExcelWriter(temporary, engine="openpyxl") as writer:
        for name, table in tables.items():
            safe_name = name[:31] or "data"
            table.to_excel(writer, sheet_name=safe_name, index=False)
    commit_partial(temporary, target)


def write_json(data: Any, path: str | Path) -> None:
    write_text_atomic(
        json.dumps(data, indent=2, ensure_ascii=False) + "\n",
        path,
    )


def run_command(
    command: Sequence[str],
    *,
    stdout_path: str | Path | None = None,
    stderr_path: str | Path | None = None,
    cwd: str | Path | None = None,
    env: Mapping[str, str] | None = None,
) -> None:
    stdout_target = Path(stdout_path) if stdout_path is not None else None
    stdout_temporary = partial_path(stdout_target) if stdout_target is not None else None
    if stdout_target is not None:
        stdout_target.parent.mkdir(parents=True, exist_ok=True)
    if stderr_path is not None:
        Path(stderr_path).parent.mkdir(parents=True, exist_ok=True)
    stdout_handle = (
        stdout_temporary.open("w", encoding="utf-8") if stdout_temporary is not None else None
    )
    stderr_handle = Path(stderr_path).open("w", encoding="utf-8") if stderr_path else None  # noqa: SIM115
    try:
        completed = subprocess.run(
            list(command),
            check=False,
            cwd=Path(cwd) if cwd is not None else None,
            env={**os.environ, **dict(env or {})},
            stdout=stdout_handle,
            stderr=stderr_handle,
            text=True,
        )
    finally:
        if stdout_handle is not None:
            stdout_handle.close()
        if stderr_handle is not None:
            stderr_handle.close()
    if completed.returncode != 0:
        rendered = " ".join(command)
        partial_note = (
            f" Partial stdout retained at {stdout_temporary}."
            if stdout_temporary is not None and stdout_temporary.exists()
            else ""
        )
        raise RuntimeError(
            f"Command failed with exit code {completed.returncode}: {rendered}.{partial_note}"
        )
    if stdout_target is not None and stdout_temporary is not None:
        commit_partial(stdout_temporary, stdout_target)


def executable_version(candidates: Sequence[str], arguments: Sequence[str]) -> tuple[str, str]:
    for executable in candidates:
        try:
            completed = subprocess.run(
                [executable, *arguments],
                check=False,
                capture_output=True,
                text=True,
                timeout=30,
            )
        except (FileNotFoundError, subprocess.TimeoutExpired):
            continue
        text = (completed.stdout or completed.stderr).strip().splitlines()
        return executable, text[0] if text else "version unavailable"
    raise FileNotFoundError(f"None of the executables were found: {', '.join(candidates)}")


def reverse_complement(sequence: str) -> str:
    table = str.maketrans("ACGTRYMKBDHVNacgtrymkbdhvn", "TGCAYRKMVHDBNtgcayrkmvhdbn")
    return sequence.translate(table)[::-1]


def split_multi_value(value: Any) -> list[str]:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return []
    text = str(value).strip()
    if not text or text.upper() == "NA":
        return []
    normalized = text.replace(";", ",")
    return [item.strip() for item in normalized.split(",") if item.strip()]


def resolve_column(
    df: pd.DataFrame, candidates: Sequence[str], required: bool = True
) -> str | None:
    lookup = {column.strip().lower().replace(" ", "_"): column for column in df.columns}
    for candidate in candidates:
        key = candidate.strip().lower().replace(" ", "_")
        if key in lookup:
            return lookup[key]
    if required:
        raise ValueError(f"None of the required columns were found: {', '.join(candidates)}")
    return None
