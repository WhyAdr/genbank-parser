"""Deterministic report serialization and protected-output contracts."""

import hashlib
import json
from pathlib import Path

import pytest

from genbank_parser.mobilome.models import (
    AnalysisParameters,
    MobilomeOutputError,
    MobilomeReport,
)
from genbank_parser.mobilome.report import (
    serialize_mobilome_report,
    write_mobilome_report,
)


def _report() -> MobilomeReport:
    return MobilomeReport(
        source_file="fixture.gb",
        source_sha256=hashlib.sha256(b"fixture\n").hexdigest(),
        tool_version="0.5.0",
        catalog_version="1.0.0",
        inference_version="1.0.0",
        resources=(),
        parameters=AnalysisParameters(
            include="all", min_evidence=1, database_source="packaged"
        ),
        annotation_provenance=(),
        inventory=(),
        replicons=(),
        cross_record_hypotheses=(),
        disabled_aggregate_rules=(),
        handoffs=(),
        limitations=(),
    )


def test_empty_report_serializes_deterministically_in_all_formats() -> None:
    report = _report()

    rendered_json = serialize_mobilome_report(report, "json")

    assert rendered_json == serialize_mobilome_report(report, "json")
    assert json.loads(rendered_json)["schema_version"] == "gbparse.mobilome.v1"
    assert serialize_mobilome_report(report, "tsv").startswith(
        "schema_version\trow_type"
    )
    assert "Mobilome evidence report" in serialize_mobilome_report(report, "text")


def test_writer_is_atomic_and_protects_source_path(tmp_path: Path) -> None:
    source = tmp_path / "fixture.gb"
    source.write_bytes(b"fixture\n")
    report = _report()
    output = tmp_path / "nested" / "report.json"

    assert (
        write_mobilome_report(
            report,
            source_path=source,
            database_paths=(),
            format_type="json",
            output_path=output,
        )
        == ""
    )
    assert output.exists()
    with pytest.raises(MobilomeOutputError, match="must not overwrite"):
        write_mobilome_report(
            report,
            source_path=source,
            database_paths=(),
            format_type="json",
            output_path=source,
        )
