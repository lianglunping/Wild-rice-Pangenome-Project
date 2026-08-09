from __future__ import annotations

from pathlib import Path

import pytest

from panfamflow.workflow.scripts.workflow_utils import fasta_lengths, iter_fasta_records


def test_fasta_lengths_streaming_summary(tmp_path: Path) -> None:
    fasta = tmp_path / "sequences.fa"
    fasta.write_text(">a\nAAAA\nAA\n>b description\nTTT\n", encoding="utf-8")
    assert dict(fasta_lengths(fasta)) == {"a": 6, "b": 3}
    assert list(iter_fasta_records(fasta)) == [("a", "AAAAAA"), ("b", "TTT")]


def test_fasta_duplicate_identifier_is_rejected(tmp_path: Path) -> None:
    fasta = tmp_path / "duplicate.fa"
    fasta.write_text(">a\nAA\n>a\nTT\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Duplicate FASTA identifier"):
        list(iter_fasta_records(fasta))


def test_run_command_does_not_publish_partial_stdout(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import subprocess

    from panfamflow.workflow.scripts.workflow_utils import run_command

    output = tmp_path / "result.tsv"

    def fake_run(*args: object, **kwargs: object) -> subprocess.CompletedProcess[str]:
        stdout = kwargs["stdout"]
        assert hasattr(stdout, "write")
        stdout.write("partial\n")  # type: ignore[union-attr]
        return subprocess.CompletedProcess(args=[], returncode=1)

    monkeypatch.setattr(subprocess, "run", fake_run)
    with pytest.raises(RuntimeError, match="exit code 1"):
        run_command(["fake"], stdout_path=output)
    assert not output.exists()
    partials = list(tmp_path.glob(".result.tsv.partial.*"))
    assert len(partials) == 1
    assert partials[0].read_text(encoding="utf-8") == "partial\n"


def test_run_command_atomically_publishes_success(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import subprocess

    from panfamflow.workflow.scripts.workflow_utils import run_command

    output = tmp_path / "result.tsv"

    def fake_run(*args: object, **kwargs: object) -> subprocess.CompletedProcess[str]:
        stdout = kwargs["stdout"]
        assert hasattr(stdout, "write")
        stdout.write("complete\n")  # type: ignore[union-attr]
        return subprocess.CompletedProcess(args=[], returncode=0)

    monkeypatch.setattr(subprocess, "run", fake_run)
    run_command(["fake"], stdout_path=output)
    assert output.read_text(encoding="utf-8") == "complete\n"
    assert not list(tmp_path.glob(".result.tsv.partial.*"))
