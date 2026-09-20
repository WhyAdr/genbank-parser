"""v0.9.1 hardening regressions for cohort identity, publication, and resume."""

from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

import jsonschema
import pytest

from genbank_parser.bakta import run_batch_summary
from genbank_parser.batch import execute_batch
from genbank_parser.cli import main
from genbank_parser.index import inspect_index, query_index
from genbank_parser.query import QueryExpressionError, query_features


def test_index_identity_and_all_failed_rebuild_preserves_known_good(
    simple_cds_gbff: Path, tmp_path: Path
) -> None:
    index = tmp_path / "cohort.gbidx"
    assert main(["index", "build", str(simple_cds_gbff), str(index)]) == 0
    assert main(["index", "update", str(index), str(simple_cds_gbff.resolve())]) == 0
    assert inspect_index(index)["counts"] == {"sources": 1, "records": 1, "features": 5}

    before = hashlib.sha256(index.read_bytes()).hexdigest()
    bad = tmp_path / "bad.gb"
    bad.write_text("not a GenBank file\n", encoding="utf-8")
    assert main(["index", "build", str(bad), str(index), "--force", "--on-error", "skip"]) == 3
    assert hashlib.sha256(index.read_bytes()).hexdigest() == before


def test_index_query_type_parity_and_complete_projection(
    simple_cds_gbff: Path, tmp_path: Path
) -> None:
    index = tmp_path / "cohort.gbidx"
    assert main(["index", "build", str(simple_cds_gbff), str(index)]) == 0
    expressions = (
        'start contains "1"',
        'length contains "3"',
        'start == "1"',
        'start != "1"',
        'casefold(type) == "cds"',
    )
    for expression in expressions:
        direct = [match.row.feature_index for match in query_features(simple_cds_gbff, expression)]
        indexed = [row.feature_index for row in query_index(index, expression)]
        assert indexed == direct

    for expression in ('gene > 3', 'gene < "zzz"', 'ko > 3'):
        with pytest.raises(QueryExpressionError):
            query_features(simple_cds_gbff, expression)
        with pytest.raises(QueryExpressionError):
            query_index(index, expression)

    direct_row = query_features(simple_cds_gbff.resolve(), 'feature_index == 3')[0].row.to_dict()
    indexed_row = query_index(index, 'feature_index == 3')[0].to_dict()
    assert set(direct_row["xrefs"]) == {
        "cog_ids",
        "db_xrefs",
        "ec_numbers",
        "go_terms",
        "kegg_kos",
        "pfam",
        "rfam",
    }
    assert set(indexed_row["xrefs"]) == set(direct_row["xrefs"])
    assert indexed_row["sample_key"] == "simple_cds"
    for key in direct_row:
        assert indexed_row[key] == direct_row[key]


def test_index_query_multivalue_rhs_preserves_direct_semantics(
    simple_cds_gbff: Path, tmp_path: Path
) -> None:
    source = tmp_path / "multivalue.gb"
    source.write_text(
        simple_cds_gbff.read_text(encoding="utf-8").replace(
            '                     /translation="MKVLWAGLIT"',
            '                     /note="first"\n'
            '                     /note="KEGG:K00844"\n'
            '                     /translation="MKVLWAGLIT"',
            1,
        ),
        encoding="utf-8",
    )
    index = tmp_path / "cohort.gbidx"
    assert main(["index", "build", str(source), str(index)]) == 0
    expression = 'qualifier("db_xref") == qualifier("note")'
    direct = [match.row.feature_index for match in query_features(source, expression)]
    indexed = [row.feature_index for row in query_index(index, expression)]
    assert indexed == direct


def test_index_report_is_preflighted_before_database_mutation(
    simple_cds_gbff: Path, multi_record_circular_gbff: Path, tmp_path: Path
) -> None:
    index = tmp_path / "cohort.gbidx"
    report = tmp_path / "report.json"
    assert main(["index", "build", str(simple_cds_gbff), str(index), "--report", str(report)]) == 0
    before = hashlib.sha256(index.read_bytes()).hexdigest()
    assert main(["index", "update", str(index), str(multi_record_circular_gbff), "--report", str(report)]) == 4
    assert hashlib.sha256(index.read_bytes()).hexdigest() == before
    assert main(
        [
            "index",
            "update",
            str(index),
            str(multi_record_circular_gbff),
            "--report",
            str(report),
            "--report-force",
        ]
    ) == 0


def test_batch_resume_preserves_threshold_outcome_and_checks_log(
    tmp_path: Path,
) -> None:
    run_dir = tmp_path / "run"
    first = execute_batch(
        ["tests/fixtures/duplicate_locus.gb"],
        command="validate",
        output_dir=run_dir,
        tail=("--format", "json", "--fail-on", "warning"),
    )
    assert first.exit_code == 1
    second = execute_batch(
        ["tests/fixtures/duplicate_locus.gb"],
        command="validate",
        output_dir=run_dir,
        tail=("--format", "json", "--fail-on", "warning"),
        resume=True,
    )
    assert second.exit_code == 1
    assert second.manifest["jobs"][0]["status"] == "skipped_unchanged"
    log = run_dir / str(second.manifest["jobs"][0]["stderr_log"])
    log.write_text("tampered\n", encoding="utf-8")
    third = execute_batch(
        ["tests/fixtures/duplicate_locus.gb"],
        command="validate",
        output_dir=run_dir,
        tail=("--format", "json", "--fail-on", "warning"),
        resume=True,
    )
    assert third.exit_code == 1
    assert third.manifest["jobs"][0]["status"] == "threshold_failed"


def test_batch_adapter_contracts_use_parsed_options(simple_cds_gbff: Path, tmp_path: Path) -> None:
    cases = (
        ("validate", (), "result.txt"),
        ("validate", ("--format=json",), "result.json"),
        ("validate", ("--json",), "result.json"),
        ("query", ("--where", 'type == "CDS"', "--emit", "faa"), "result.faa"),
        ("phylo", (), "report.txt"),
        ("neighborhood", ("TEST_001",), "result.txt"),
    )
    for number, (command, tail, expected) in enumerate(cases):
        result = execute_batch(
            [simple_cds_gbff],
            command=command,
            output_dir=tmp_path / f"run-{number}",
            tail=tail,
        )
        assert result.exit_code == 0
        assert expected in {item["relative_path"] for item in result.manifest["jobs"][0]["outputs"]}


def test_every_registered_batch_adapter_executes(simple_cds_gbff: Path, tmp_path: Path) -> None:
    cases = (
        ("validate", ()),
        ("summary", ()),
        ("metadata", ()),
        ("extract", ()),
        ("search", ()),
        ("locus", ("TEST_001",)),
        ("neighborhood", ("TEST_001",)),
        ("region", ("--record", "CONTIG_1", "--start", "1", "--end", "30")),
        ("fasta", ()),
        ("sequence", ()),
        ("codon", ()),
        ("functional", ()),
        ("discover", ()),
        ("phylo", ()),
        ("crispr", ()),
        ("gff", ()),
        ("meor", ()),
        ("mobilome", ()),
        ("query", ("--where", 'type == "CDS"')),
        ("operons", ()),
        ("export", ("--format", "bed12")),
    )
    for number, (command, tail) in enumerate(cases):
        result = execute_batch(
            [simple_cds_gbff],
            command=command,
            output_dir=tmp_path / f"adapter-{number}-{command}",
            tail=tail,
        )
        assert result.exit_code == 0, command
        assert result.manifest["jobs"][0]["status"] == "succeeded"


def test_records_reject_malformed_zero_record_input(tmp_path: Path) -> None:
    malformed = tmp_path / "bad.gb"
    malformed.write_text("not a GenBank file\n", encoding="utf-8")
    assert main(["records", "list", str(malformed), "--format", "json"]) == 3


def test_batch_summary_fail_is_fail_closed_and_skip_publishes_failures(
    simple_cds_gbff: Path, tmp_path: Path
) -> None:
    sources = tmp_path / "sources"
    sources.mkdir()
    shutil.copyfile(simple_cds_gbff, sources / "good.gb")
    (sources / "bad.gb").write_text("not a GenBank file\n", encoding="utf-8")
    csv_path = tmp_path / "summary.csv"
    tsv_path = tmp_path / "summary.tsv"
    md_path = tmp_path / "summary.md"
    failed = run_batch_summary(
        [sources], csv_path, tsv_path, md_path, on_error="fail"
    )
    assert failed.exit_code == 3
    assert not csv_path.exists() and not tsv_path.exists() and not md_path.exists()
    skipped = run_batch_summary(
        [sources], csv_path, tsv_path, md_path, on_error="skip", force=True
    )
    assert skipped.exit_code == 3
    failure_path = tmp_path / "summary.failures.json"
    payload = json.loads(failure_path.read_text(encoding="utf-8"))
    assert payload["schema_version"] == "gbparse.batch-summary.v1"
    assert payload["failures"]


def test_batch_manifest_schema_and_index_report_schema(
    simple_cds_gbff: Path, tmp_path: Path
) -> None:
    result = execute_batch(
        [simple_cds_gbff], command="validate", output_dir=tmp_path / "run", tail=("--format", "json")
    )
    from importlib import resources

    batch_schema = json.loads(
        resources.files("genbank_parser").joinpath("data", "schemas", "batch-v1.schema.json").read_text()
    )
    jsonschema.validate(result.manifest, batch_schema)
    report = tmp_path / "index-report.json"
    index = tmp_path / "cohort.gbidx"
    assert main(["index", "build", str(simple_cds_gbff), str(index), "--report", str(report)]) == 0
    index_schema = json.loads(
        resources.files("genbank_parser").joinpath("data", "schemas", "index-report-v1.schema.json").read_text()
    )
    jsonschema.validate(json.loads(report.read_text(encoding="utf-8")), index_schema)
