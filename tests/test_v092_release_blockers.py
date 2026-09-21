"""Release-blocker regressions for the v0.9.2 hardening plan."""

from __future__ import annotations

import copy
import gzip
import hashlib
import json
import os
import shutil
import sqlite3
import subprocess
from importlib import resources
from pathlib import Path

import jsonschema
import pytest

from genbank_parser import read_genbank
from genbank_parser.batch import _validate_manifest, execute_batch
from genbank_parser.cli import main
from genbank_parser.cli_io import InputError, OutputRecoveryError, publish_file_set
from genbank_parser.index import build_index, inspect_index, query_index, update_index
from genbank_parser.io import GenBankInputError
from genbank_parser.query import QueryExpressionError, parse_select_fields


def _copy_source(source: Path, directory: Path, name: str) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    target = directory / name
    shutil.copyfile(source, target)
    return target


def _adapter_cases() -> tuple[tuple[str, tuple[str, ...]], ...]:
    return (
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


def test_plain_and_gzip_sources_bind_the_raw_digest_once(
    simple_cds_gbff: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    gzip_source = tmp_path / "simple_cds.gb.gz"
    raw = simple_cds_gbff.read_bytes()
    gzip_source.write_bytes(gzip.compress(raw))

    for source in (simple_cds_gbff, gzip_source):
        monkeypatch.setenv("GBPARSE_SOURCE_LABEL", f"logical/{source.name}")
        monkeypatch.setenv("GBPARSE_SOURCE_SHA256", hashlib.sha256(source.read_bytes()).hexdigest())
        document = read_genbank(source)
        assert document.source_label == f"logical/{source.name}"
        assert document.all_features

    monkeypatch.delenv("GBPARSE_SOURCE_SHA256")
    with pytest.raises(GenBankInputError, match="provided together"):
        read_genbank(simple_cds_gbff)

    monkeypatch.setenv("GBPARSE_SOURCE_SHA256", "0" * 64)
    with pytest.raises(GenBankInputError, match="digest mismatch"):
        read_genbank(simple_cds_gbff)

    monkeypatch.delenv("GBPARSE_SOURCE_LABEL")
    monkeypatch.delenv("GBPARSE_SOURCE_SHA256")
    truncated = tmp_path / "truncated.gb.gz"
    truncated.write_bytes(b"\x1f\x8b\x08\x00not-a-complete-gzip")
    with pytest.raises(GenBankInputError, match="truncated or invalid gzip"):
        read_genbank(truncated)


def test_gzip_batch_preserves_logical_source_and_succeeds(
    simple_cds_gbff: Path, tmp_path: Path
) -> None:
    gzip_source = tmp_path / "simple_cds.gb.gz"
    with gzip_source.open("wb") as handle:
        handle.write(simple_cds_gbff.read_bytes())
    source = gzip_source.resolve()
    result = execute_batch(
        [source],
        command="metadata",
        output_dir=tmp_path / "gzip-run",
        tail=("--format", "json"),
    )
    assert result.exit_code == 0
    output = tmp_path / "gzip-run" / "outputs" / "simple_cds" / "result.json"
    assert json.loads(output.read_text(encoding="utf-8"))["source"] == str(source)


def test_all_batch_adapters_scrub_transient_paths(
    simple_cds_gbff: Path, tmp_path: Path
) -> None:
    forbidden = (".gbparse-input.", ".gbparse-inprogress", ".work.")
    for number, (command, tail) in enumerate(_adapter_cases()):
        run = tmp_path / f"adapter-{number}-{command}"
        result = execute_batch(
            [simple_cds_gbff],
            command=command,
            output_dir=run,
            tail=tail,
        )
        assert result.exit_code == 0, command
        for path in run.rglob("*"):
            if not path.is_file():
                continue
            text = path.read_text(encoding="utf-8", errors="ignore")
            assert not any(token in text for token in forbidden), path


def test_child_interrupt_cleans_snapshot_and_leaves_resume_manifest(
    simple_cds_gbff: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    real_run = subprocess.run

    def interrupt_child(*args: object, **kwargs: object) -> subprocess.CompletedProcess[str]:
        argv = args[0] if args else kwargs.get("args", ())
        if isinstance(argv, (list, tuple)) and "genbank_parser.cli" in argv:
            raise KeyboardInterrupt
        return real_run(*args, **kwargs)  # type: ignore[arg-type]

    monkeypatch.setattr("genbank_parser.batch.subprocess.run", interrupt_child)
    run = tmp_path / "interrupted"
    with pytest.raises(KeyboardInterrupt):
        execute_batch([simple_cds_gbff], command="validate", output_dir=run)

    progress = tmp_path / ".interrupted.gbparse-inprogress"
    manifest = json.loads((progress / "batch-manifest.json").read_text(encoding="utf-8"))
    assert manifest["jobs"][0]["status"] == "running"
    assert not list(tmp_path.glob(".gbparse-input.*"))
    assert not list(progress.glob(".gbparse-input.*"))


def _assert_recovery_paths_exist(error: OutputRecoveryError) -> None:
    assert error.backups or error.staging
    assert all(path.exists() for path in (*error.backups, *error.staging))


def test_publish_file_set_retains_every_reported_recovery_path(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    first = tmp_path / "first.txt"
    second = tmp_path / "second.txt"
    first.write_text("old-first", encoding="utf-8")
    second.write_text("old-second", encoding="utf-8")
    real_replace = os.replace

    def fail_install(source: str | Path, target: str | Path) -> None:
        source_path = Path(source)
        target_path = Path(target)
        if target_path == second and source_path.suffix == ".tmp":
            raise OSError("install failed")
        if target_path == first and ".backup" in source_path.name:
            raise OSError("restore failed")
        real_replace(source, target)

    monkeypatch.setattr("genbank_parser.cli_io.os.replace", fail_install)
    with pytest.raises(OutputRecoveryError) as caught:
        publish_file_set({first: "new-first", second: "new-second"}, force=True)
    _assert_recovery_paths_exist(caught.value)


@pytest.mark.parametrize("operation", ("build", "update"))
def test_index_publication_retains_recovery_paths(
    simple_cds_gbff: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    operation: str,
) -> None:
    index = tmp_path / f"{operation}.gbidx"
    report = tmp_path / f"{operation}-report.json"
    assert build_index([simple_cds_gbff], index, report_path=report).report
    old_index = index.read_bytes()
    old_report = report.read_bytes()
    real_replace = os.replace

    def fail_install(source: str | Path, target: str | Path) -> None:
        source_path = Path(source)
        target_path = Path(target)
        if target_path == report and source_path.suffix == ".tmp":
            raise OSError("report install failed")
        if target_path == index and ".backup" in source_path.name:
            raise OSError("database restore failed")
        real_replace(source, target)

    monkeypatch.setattr("genbank_parser.cli_io.os.replace", fail_install)
    with pytest.raises(OutputRecoveryError) as caught:
        if operation == "build":
            build_index(
                [simple_cds_gbff],
                index,
                force=True,
                report_path=report,
                report_force=True,
            )
        else:
            update_index(
                index,
                [simple_cds_gbff],
                report_path=report,
                report_force=True,
            )
    _assert_recovery_paths_exist(caught.value)
    assert report.read_bytes() == old_report
    if index.exists():
        assert index.read_bytes() == old_index


def test_stop_latch_is_preserved_and_changed_success_clears_it(
    duplicate_locus_gbff: Path,
    simple_cds_gbff: Path,
    tmp_path: Path,
) -> None:
    sources = tmp_path / "sources"
    first = _copy_source(duplicate_locus_gbff, sources, "a.gb")
    _copy_source(duplicate_locus_gbff, sources, "b.gb")
    run = tmp_path / "run"
    first_result = execute_batch(
        [sources],
        command="validate",
        output_dir=run,
        tail=("--format", "json", "--fail-on", "warning"),
        on_error="stop",
    )
    assert first_result.manifest["jobs"][0]["status"] == "threshold_failed"
    assert first_result.manifest["jobs"][1]["status"] == "pending"

    unchanged = execute_batch(
        [sources],
        command="validate",
        output_dir=run,
        tail=("--format", "json", "--fail-on", "warning"),
        on_error="stop",
        resume=True,
    )
    assert unchanged.manifest["jobs"][0]["resume_action"] == "reused_unchanged"
    assert unchanged.manifest["jobs"][1]["status"] == "pending"
    assert unchanged.manifest["jobs"][1]["resume_action"] == "deferred_after_stop"

    shutil.copyfile(simple_cds_gbff, first)
    resumed = execute_batch(
        [sources],
        command="validate",
        output_dir=run,
        tail=("--format", "json", "--fail-on", "warning"),
        on_error="stop",
        resume=True,
    )
    assert resumed.manifest["jobs"][0]["status"] == "succeeded"
    assert resumed.manifest["jobs"][1]["status"] == "threshold_failed"
    assert resumed.manifest["jobs"][0]["resume_action"] == "rerun_changed"


def test_unchanged_ordinary_failure_is_reused_by_stop_resume(
    simple_cds_gbff: Path, tmp_path: Path
) -> None:
    sources = tmp_path / "sources"
    bad = sources / "a.gb"
    sources.mkdir()
    bad.write_text("not GenBank", encoding="utf-8")
    _copy_source(simple_cds_gbff, sources, "b.gb")
    run = tmp_path / "run"
    first = execute_batch([sources], command="validate", output_dir=run, on_error="stop")
    assert first.manifest["jobs"][0]["status"] == "failed"
    assert first.manifest["jobs"][1]["status"] == "pending"
    second = execute_batch(
        [sources], command="validate", output_dir=run, on_error="stop", resume=True
    )
    assert second.manifest["jobs"][0]["resume_action"] == "reused_unchanged"
    assert second.manifest["jobs"][1]["status"] == "pending"
    assert second.manifest["jobs"][1]["resume_action"] == "deferred_after_stop"


def test_continue_mode_executes_pending_after_reused_failure(
    simple_cds_gbff: Path, tmp_path: Path
) -> None:
    sources = tmp_path / "sources"
    sources.mkdir()
    (sources / "a.gb").write_text("not GenBank", encoding="utf-8")
    _copy_source(simple_cds_gbff, sources, "b.gb")
    run = tmp_path / "run"
    execute_batch([sources], command="validate", output_dir=run, on_error="continue")
    result = execute_batch(
        [sources], command="validate", output_dir=run, on_error="continue", resume=True
    )
    assert result.manifest["jobs"][0]["resume_action"] == "reused_unchanged"
    assert result.manifest["jobs"][1]["status"] == "succeeded"


@pytest.mark.parametrize("field", ("sample", "sample_key"))
def test_direct_query_rejects_cohort_projection_fields(
    simple_cds_gbff: Path, field: str, capsys: pytest.CaptureFixture[str]
) -> None:
    with pytest.raises(QueryExpressionError, match="cohort index"):
        parse_select_fields(field, allow_cohort_fields=False)
    assert main(
        [
            "query",
            str(simple_cds_gbff),
            "--where",
            'type == "CDS"',
            "--select",
            field,
        ]
    ) == 2
    assert "cohort index" in capsys.readouterr().err


def test_index_cohort_projection_remains_valid(
    simple_cds_gbff: Path, tmp_path: Path
) -> None:
    index = tmp_path / "cohort.gbidx"
    assert main(["index", "build", str(simple_cds_gbff), str(index)]) == 0
    rows = query_index(index, 'type == "CDS"')
    assert rows and all(row.sample_key == "simple_cds" for row in rows)
    assert main(
        [
            "index",
            "query",
            str(index),
            "--where",
            'type == "CDS"',
            "--select",
            "sample_key",
        ]
    ) == 0


def _mutate_sample_key_index(index: Path, variant: str) -> None:
    connection = sqlite3.connect(index)
    connection.execute("DROP INDEX uq_sources_sample_key_nocase")
    if variant == "nonunique":
        connection.execute(
            "CREATE INDEX uq_sources_sample_key_nocase ON sources(sample_key)"
        )
    elif variant == "binary":
        connection.execute(
            "CREATE UNIQUE INDEX uq_sources_sample_key_nocase ON sources(sample_key)"
        )
    elif variant == "extra-column":
        connection.execute(
            "CREATE UNIQUE INDEX uq_sources_sample_key_nocase "
            "ON sources(sample_key COLLATE NOCASE, source_pk)"
        )
    elif variant == "partial":
        connection.execute(
            "CREATE UNIQUE INDEX uq_sources_sample_key_nocase "
            "ON sources(sample_key COLLATE NOCASE) WHERE sample_key IS NOT NULL"
        )
    else:
        raise AssertionError(variant)
    connection.commit()
    connection.close()


@pytest.mark.parametrize("variant", ("nonunique", "binary", "extra-column", "partial"))
def test_index_status_and_queries_validate_the_actual_definition(
    simple_cds_gbff: Path, tmp_path: Path, variant: str
) -> None:
    valid = tmp_path / "valid.gbidx"
    invalid = tmp_path / f"{variant}.gbidx"
    assert main(["index", "build", str(simple_cds_gbff), str(valid)]) == 0
    shutil.copyfile(valid, invalid)
    _mutate_sample_key_index(invalid, variant)
    with pytest.raises(ValueError, match="sample-key"):
        inspect_index(invalid)
    with pytest.raises(ValueError, match="sample-key"):
        query_index(invalid, 'type == "CDS"')
    assert main(["index", "status", str(invalid)]) == 3


def test_replay_working_directory_is_recorded_and_legacy_manifests_remain_readable(
    simple_cds_gbff: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)
    relative_source = Path("inputs") / "sample.gb"
    _copy_source(simple_cds_gbff, relative_source.parent, relative_source.name)
    result = execute_batch(
        [relative_source], command="validate", output_dir=Path("run")
    )
    assert Path(result.manifest["runner"]["working_directory"]) == tmp_path.resolve()

    legacy = copy.deepcopy(result.manifest)
    legacy["gbparse_version"] = "0.9.1"
    legacy["runner"].pop("working_directory")
    assert _validate_manifest(legacy) is legacy

    current_without_cwd = copy.deepcopy(result.manifest)
    current_without_cwd["runner"].pop("working_directory")
    with pytest.raises(InputError, match="working_directory"):
        _validate_manifest(current_without_cwd)

    relative_cwd = copy.deepcopy(result.manifest)
    relative_cwd["runner"]["working_directory"] = "."
    with pytest.raises(InputError, match="absolute"):
        _validate_manifest(relative_cwd)


def test_index_migrate_cli_is_exposed(simple_cds_gbff: Path, tmp_path: Path) -> None:
    index = tmp_path / "cohort.gbidx"
    assert main(["index", "build", str(simple_cds_gbff), str(index)]) == 0
    assert main(["index", "migrate", str(index)]) == 0


def test_index_status_instance_matches_packaged_schema(
    simple_cds_gbff: Path, tmp_path: Path
) -> None:
    index = tmp_path / "cohort.gbidx"
    status = tmp_path / "status.json"
    assert main(["index", "build", str(simple_cds_gbff), str(index)]) == 0
    assert main(
        ["index", "status", str(index), "--format", "json", "--output", str(status)]
    ) == 0
    schema = json.loads(
        resources.files("genbank_parser")
        .joinpath("data", "schemas", "index-status-v1.schema.json")
        .read_text(encoding="utf-8")
    )
    jsonschema.validate(json.loads(status.read_text(encoding="utf-8")), schema)
