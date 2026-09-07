"""Serializer, schema, normalized TSV, and protected-output regressions."""

import csv
import io
import json
from importlib import resources
from pathlib import Path

import jsonschema
import pytest

from genbank_parser.mobilome import analyze_mobilome
from genbank_parser.mobilome.models import MobilomeOutputError
from genbank_parser.mobilome.report import (
    TSV_COLUMNS,
    mobilome_report_to_dict,
    serialize_mobilome_report,
    write_mobilome_report,
)


def test_json_and_tsv_are_schema_valid_and_lossless_at_reason_level() -> None:
    report = analyze_mobilome(Path("tests/fixtures/mobilome_evidence.gb"))
    payload = mobilome_report_to_dict(report)
    schema = json.loads(
        resources.files("genbank_parser")
        .joinpath("data", "mobilome", "report.schema.json")
        .read_text(encoding="utf-8")
    )
    jsonschema.Draft202012Validator.check_schema(schema)
    jsonschema.validate(payload, schema)

    rendered_json = serialize_mobilome_report(report, "json")
    rendered_tsv = serialize_mobilome_report(report, "tsv")
    rows = list(csv.DictReader(io.StringIO(rendered_tsv), delimiter="\t"))
    raw_reasons = sum(
        len(hit.reasons) for assessment in report.replicons for hit in assessment.hits
    )

    assert rendered_json == serialize_mobilome_report(report, "json")
    assert json.loads(rendered_json)["summary"]["raw_reasons"] == raw_reasons
    assert tuple(rows[0]) == TSV_COLUMNS
    assert len([row for row in rows if row["row_type"] == "inventory"]) == len(
        report.inventory
    )
    assert len([row for row in rows if row["row_type"] == "evidence"]) == raw_reasons
    assert all(
        json.loads(row["segments_json"])
        for row in rows
        if row["row_type"] == "evidence"
    )


def test_mobilome_serializers_match_versioned_goldens() -> None:
    report = analyze_mobilome(Path("tests/fixtures/mobilome_evidence.gb"))
    golden_dir = Path("tests/golden")

    assert serialize_mobilome_report(report, "json") == (
        golden_dir / "mobilome_v2.json"
    ).read_bytes().decode("utf-8")
    assert serialize_mobilome_report(report, "tsv") == (
        golden_dir / "mobilome_v2.tsv"
    ).read_bytes().decode("utf-8")
    assert serialize_mobilome_report(report, "text") == (
        golden_dir / "mobilome_v2.txt"
    ).read_bytes().decode("utf-8")


def test_include_filters_details_but_not_inventory_or_scanned_summary() -> None:
    report = analyze_mobilome(
        Path("tests/fixtures/mobilome_inference.gb"), include="plasmid"
    )
    payload = mobilome_report_to_dict(report)

    assert len(report.replicons) == 1
    assert len(report.inventory) == 3
    assert payload["summary"]["records_scanned"] == 3
    assert payload["summary"]["raw_hits"] > sum(
        len(assessment.hits) for assessment in report.replicons
    )


def test_writer_is_atomic_and_protects_input_and_custom_database_paths(
    tmp_path: Path,
) -> None:
    source = Path("tests/fixtures/mobilome_evidence.gb")
    before = source.read_bytes()
    report = analyze_mobilome(source)
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
    assert source.read_bytes() == before
    with pytest.raises(MobilomeOutputError, match="must not overwrite"):
        write_mobilome_report(
            report,
            source_path=source,
            database_paths=(),
            format_type="json",
            output_path=source,
        )
    protected = tmp_path / "markers.yaml"
    protected.write_text("markers", encoding="utf-8")
    with pytest.raises(MobilomeOutputError, match="must not overwrite"):
        write_mobilome_report(
            report,
            source_path=source,
            database_paths=(protected,),
            format_type="json",
            output_path=protected,
        )
