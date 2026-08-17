"""Stable JSON, TSV, and text reporting contracts."""
from __future__ import annotations

import json

from genbank_parser.meor import analyze_meor
from genbank_parser.meor.report import TSV_COLUMNS, render_json, render_text, render_tsv


def test_report_contracts() -> None:
    report = analyze_meor("tests/fixtures/meor_parity.gb", window_size=1000)
    payload = json.loads(render_json(report))
    assert list(payload) == [
        "schema_version",
        "tool",
        "analysis",
        "tool_version",
        "catalog_version",
        "file",
        "total_features",
        "total_hits",
        "total_clusters",
        "parameters",
        "hits",
        "clusters",
        "pathways",
        "karyograms",
    ]
    assert payload["schema_version"] == "gbparse.meor.v1"
    assert payload["tool_version"] == "0.5.0"
    assert payload["catalog_version"] == "1.1"
    assert payload["parameters"] == {
        "min_weight": 1,
        "max_gap": 200,
        "window_size": 1000,
    }
    assert render_tsv(report).splitlines()[0].split("\t") == list(TSV_COLUMNS)
    text = render_text(report)
    for section in (
        "[1] CATEGORY SUMMARY & HIT COUNTS",
        "[2] MEOR PATHWAY COMPLETENESS PROFILE",
        "[3] CO-LOCALIZED MEOR OPERONS & BGC CLUSTERS",
        "[4] GENOMIC DENSITY KARYOGRAMS",
    ):
        assert section in text
    assert "Low (W=1)" in text
    assert "1 kb windows" in text
    assert "gbparse Version: 0.5.0" in text
    assert "MEOR Catalog Version: 1.1" in text


def test_min_weight_filters_hits_and_pathways() -> None:
    reports = {
        weight: analyze_meor("tests/fixtures/meor_parity.gb", min_weight=weight)
        for weight in (1, 2, 3)
    }
    assert [reports[weight].total_hits for weight in (1, 2, 3)] == [11, 9, 7]
    long_chain = {
        weight: reports[weight].pathways[1].completeness_pct for weight in (1, 2, 3)
    }
    assert long_chain == {1: 100.0, 2: 50.0, 3: 0.0}
