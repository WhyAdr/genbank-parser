"""CLI orchestration, stdout, and error-surface tests for gbparse mobilome."""

import json
from pathlib import Path

import pytest

from genbank_parser.cli import main


def test_mobilome_cli_formats_and_output_behavior(tmp_path: Path, capsys) -> None:
    fixture = Path("tests/fixtures/mobilome_evidence.gb")

    assert main(["mobilome", str(fixture), "--format", "json"]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["schema_version"] == "gbparse.mobilome.v2"

    output = tmp_path / "report.tsv"
    assert (
        main(
            [
                "mobilome",
                str(fixture),
                "--format",
                "tsv",
                "--output",
                str(output),
            ]
        )
        == 0
    )
    assert capsys.readouterr().out == ""
    assert output.read_text(encoding="utf-8").startswith("schema_version\trow_type")


def test_mobilome_cli_reports_protected_output_without_traceback(capsys) -> None:
    fixture = Path("tests/fixtures/mobilome_evidence.gb")

    with pytest.raises(SystemExit):
        main(["mobilome", str(fixture), "--output", str(fixture)])
    captured = capsys.readouterr()
    assert "must not overwrite" in captured.err
    assert "Traceback" not in captured.err


def test_mobilome_cli_help_is_available(capsys) -> None:
    with pytest.raises(SystemExit) as excinfo:
        main(["mobilome", "--help"])
    assert excinfo.value.code == 0
    assert "--min-evidence" in capsys.readouterr().out
