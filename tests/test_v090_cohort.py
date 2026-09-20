"""Focused v0.9.0 cohort, record, index, and batch contracts."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

from genbank_parser.batch import execute_batch
from genbank_parser.cli import main
from genbank_parser.discovery import (
    assign_sample_keys,
    discover_inputs,
    fingerprint_file,
)
from genbank_parser.index import inspect_index, query_index
from genbank_parser.io import read_genbank
from genbank_parser.query import query_features
from genbank_parser.records import render_selected_records


def test_discovery_deduplicates_overlapping_roots_and_assigns_collision_safe_keys(
    simple_cds_gbff: Path, tmp_path: Path
) -> None:
    root = tmp_path / "cohort"
    root.mkdir()
    first = root / "one.gbff"
    second_dir = root / "nested"
    second_dir.mkdir()
    second = second_dir / "one.gbff"
    first.write_bytes(simple_cds_gbff.read_bytes())
    second.write_bytes(simple_cds_gbff.read_bytes())
    discovered = discover_inputs([root, first])
    assert [item.display_path for item in discovered] == sorted(
        [str(first), str(second)], key=lambda value: value.casefold()
    )
    keys = assign_sample_keys(discovered)
    assert len(set(keys.values())) == 2
    assert all(key.startswith("one-") for key in keys.values())
    assert fingerprint_file(first).sha256 == fingerprint_file(second).sha256


def test_records_round_trip_and_split_manifest(
    multi_record_circular_gbff: Path, tmp_path: Path
) -> None:
    document = read_genbank(multi_record_circular_gbff)
    payload = render_selected_records(tuple(document.records), "genbank")
    round_tripped = read_genbank(__import__("io").StringIO(payload.decode("utf-8")))
    assert [(record.id, record.length, len(record.features)) for record in round_tripped.records] == [
        (record.id, record.length, len(record.features)) for record in document.records
    ]
    split_dir = tmp_path / "split"
    assert main(
        [
            "records",
            "split",
            str(multi_record_circular_gbff),
            "--output-dir",
            str(split_dir),
        ]
    ) == 0
    manifest = (split_dir / "records.manifest.tsv").read_text(encoding="utf-8").splitlines()
    assert manifest[0].split("\t") == [
        "schema_version",
        "source",
        "record_index",
        "record_id",
        "record_name",
        "length",
        "filename",
        "sha256",
    ]
    assert len(manifest) == 3


def test_index_query_matches_direct_query_and_status(
    simple_cds_gbff: Path, tmp_path: Path
) -> None:
    index = tmp_path / "cohort.gbidx"
    assert main(["index", "build", str(simple_cds_gbff), str(index)]) == 0
    expressions = (
        'type == "CDS" and (gene == "testA" or ko == "K00844")',
        '"CDS" == type',
        'type != ["CDS", "source"]',
        'null == topology',
    )
    for expression in expressions:
        direct = [match.row.locus_tag for match in query_features(simple_cds_gbff, expression)]
        indexed = [row.locus_tag for row in query_index(index, expression)]
        assert indexed == direct
    status = inspect_index(index, check=True, verify_sources=True)
    assert status["schema_version"] == "gbparse.index-status.v1"
    assert status["counts"] == {"sources": 1, "records": 1, "features": 5}
    assert status["checks"]["quick_check"] == "ok"
    assert status["source_checks"][0]["sha256_matches"] is True


def test_index_update_prune_and_query_output_alias_rejection(
    simple_cds_gbff: Path, multi_record_circular_gbff: Path, tmp_path: Path
) -> None:
    source_dir = tmp_path / "sources"
    source_dir.mkdir()
    first = source_dir / "first.gb"
    second = source_dir / "second.gb"
    shutil.copyfile(simple_cds_gbff, first)
    shutil.copyfile(multi_record_circular_gbff, second)
    index = tmp_path / "cohort.gbidx"
    assert main(["index", "build", str(source_dir), str(index)]) == 0
    assert inspect_index(index)["counts"]["sources"] == 2
    second.unlink()
    assert main(["index", "update", str(index), str(source_dir), "--prune"]) == 0
    assert inspect_index(index)["counts"]["sources"] == 1
    assert main(
        [
            "index",
            "query",
            str(index),
            "--where",
            'type == "CDS"',
            "--output",
            str(index),
        ]
    ) == 4


def test_batch_manifest_resume_and_failed_job_has_no_artifact(
    simple_cds_gbff: Path, tmp_path: Path
) -> None:
    source_dir = tmp_path / "sources"
    source_dir.mkdir()
    shutil.copyfile(simple_cds_gbff, source_dir / "good.gb")
    run_dir = tmp_path / "run"
    result = execute_batch(
        [source_dir],
        command="validate",
        output_dir=run_dir,
        tail=("--format", "json"),
    )
    assert result.exit_code == 0
    manifest = json.loads((run_dir / "batch-manifest.json").read_text(encoding="utf-8"))
    assert manifest["schema_version"] == "gbparse.batch.v1"
    assert manifest["jobs"][0]["status"] == "succeeded"
    resumed = execute_batch(
        [source_dir],
        command="validate",
        output_dir=run_dir,
        tail=("--format", "json"),
        resume=True,
    )
    assert resumed.exit_code == 0
    assert resumed.manifest["jobs"][0]["status"] == "skipped_unchanged"
    (run_dir / "outputs" / "good" / "result.json").write_text("edited", encoding="utf-8")
    rerun = execute_batch(
        [source_dir],
        command="validate",
        output_dir=run_dir,
        tail=("--format", "json"),
        resume=True,
    )
    assert rerun.manifest["jobs"][0]["status"] == "succeeded"
    (source_dir / "bad.gb").write_text("not a GenBank file\n", encoding="utf-8")
    bad_run = tmp_path / "bad-run"
    failed = execute_batch(
        [source_dir],
        command="validate",
        output_dir=bad_run,
        tail=("--format", "json"),
    )
    assert failed.exit_code == 3
    bad_job = next(job for job in failed.manifest["jobs"] if job["sample_key"] == "bad")
    assert bad_job["status"] == "failed"
    assert not (bad_run / "outputs" / "bad").exists()


def test_batch_jobs_gt_one_preserves_output_artifacts(
    simple_cds_gbff: Path, tmp_path: Path
) -> None:
    source_dir = tmp_path / "sources"
    source_dir.mkdir()
    shutil.copyfile(simple_cds_gbff, source_dir / "first.gb")
    shutil.copyfile(simple_cds_gbff, source_dir / "second.gb")
    serial_dir = tmp_path / "serial"
    parallel_dir = tmp_path / "parallel"
    serial = execute_batch(
        [source_dir],
        command="validate",
        output_dir=serial_dir,
        tail=("--format", "json"),
        jobs=1,
    )
    parallel = execute_batch(
        [source_dir],
        command="validate",
        output_dir=parallel_dir,
        tail=("--format", "json"),
        jobs=2,
    )
    assert serial.exit_code == parallel.exit_code == 0
    assert {job["sample_key"] for job in serial.manifest["jobs"]} == {
        job["sample_key"] for job in parallel.manifest["jobs"]
    }
    for job in serial.manifest["jobs"]:
        serial_output = serial_dir / "outputs" / str(job["sample_key"]) / str(job["outputs"][0]["relative_path"])
        parallel_output = parallel_dir / "outputs" / str(job["sample_key"]) / str(job["outputs"][0]["relative_path"])
        assert serial_output.read_bytes() == parallel_output.read_bytes()


def test_batch_rejects_unsupported_and_runner_owned_options(simple_cds_gbff: Path, tmp_path: Path) -> None:
    assert main(
        [
            "batch",
            str(simple_cds_gbff),
            "--command",
            "compare",
            "--output-dir",
            str(tmp_path / "run"),
        ]
    ) == 2
    assert main(
        [
            "batch",
            str(simple_cds_gbff),
            "--command",
            "validate",
            "--output-dir",
            str(tmp_path / "run2"),
            "--",
            "--output",
            "bad.json",
        ]
    ) == 2
