"""Release-blocker regressions for the v0.9.1 hardening plan."""

from __future__ import annotations

import json
import shutil
import sqlite3
from pathlib import Path
from unittest.mock import patch

import pytest

from genbank_parser.batch import (
    ExpectedOutput,
    _verified_outputs_payload,
    execute_batch,
)
from genbank_parser.cli import main
from genbank_parser.cli_io import (
    OutputError,
    OutputRecoveryError,
    publish_directory_tree,
    publish_file_set,
    publish_staged_files,
)
from genbank_parser.discovery import discover_inputs
from genbank_parser.index import inspect_index, migrate_index, query_index
from genbank_parser.io import GenBankInputError, read_genbank
from genbank_parser.query import QueryExpressionError, query_features


def _copy_source(source: Path, directory: Path, name: str) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    target = directory / name
    shutil.copyfile(source, target)
    return target


def test_batch_manifest_keeps_logical_provenance_and_no_transient_paths(
    simple_cds_gbff: Path, tmp_path: Path
) -> None:
    result = execute_batch(
        [simple_cds_gbff],
        command="validate",
        output_dir=tmp_path / "run",
        tail=("--format", "json"),
    )
    job = result.manifest["jobs"][0]
    text = json.dumps(result.manifest)
    assert ".gbparse-input." not in text
    assert ".gbparse-inprogress" not in text
    assert ".work." not in text
    assert any("outputs\\" in value or "outputs/" in value for value in job["argv"])
    assert str(simple_cds_gbff) in job["argv"]


def test_direct_and_batch_metadata_share_logical_source_label(
    simple_cds_gbff: Path, tmp_path: Path
) -> None:
    from genbank_parser.metadata import build_metadata_report

    direct = build_metadata_report(simple_cds_gbff.resolve())
    result = execute_batch(
        [simple_cds_gbff.resolve()],
        command="metadata",
        output_dir=tmp_path / "run",
        tail=("--format", "json"),
    )
    rendered = json.loads(
        (tmp_path / "run" / "outputs" / "simple_cds" / "result.json").read_text(
            encoding="utf-8"
        )
    )
    assert result.exit_code == 0
    assert rendered["source"] == direct["source"]


def test_batch_requires_all_declared_multi_file_outputs(tmp_path: Path) -> None:
    sample = tmp_path / "sample"
    sample.mkdir()
    (sample / "one.txt").write_text("ok", encoding="utf-8")
    (sample / "extra.txt").write_text("unexpected", encoding="utf-8")
    with pytest.raises(OutputError, match="missing"):
        _verified_outputs_payload(
            sample,
            (ExpectedOutput("one.txt"), ExpectedOutput("two.txt")),
            capture_extra_outputs=False,
        )


def test_batch_requires_all_declared_directory_outputs(tmp_path: Path) -> None:
    sample = tmp_path / "sample"
    sample.mkdir()
    (sample / "one.txt").write_text("ok", encoding="utf-8")
    (sample / "extra.txt").write_text("unexpected", encoding="utf-8")
    with pytest.raises(OutputError, match="undeclared"):
        _verified_outputs_payload(
            sample,
            (ExpectedOutput("one.txt"),),
            capture_extra_outputs=False,
        )
    (sample / "extra.txt").unlink()
    assert _verified_outputs_payload(
        sample,
        (ExpectedOutput("one.txt"),),
        capture_extra_outputs=False,
    )[0]["relative_path"] == "one.txt"


def test_index_update_allocates_unique_keys_from_final_cohort(
    simple_cds_gbff: Path, tmp_path: Path
) -> None:
    sources = tmp_path / "sources"
    _copy_source(simple_cds_gbff, sources, "one.gb")
    index = tmp_path / "cohort.gbidx"
    assert main(["index", "build", str(sources), str(index)]) == 0
    _copy_source(simple_cds_gbff, sources, "one.gbk")
    assert main(["index", "update", str(index), str(sources)]) == 0
    keys = [row[0] for row in sqlite3.connect(index).execute("SELECT sample_key FROM sources")]
    assert len(keys) == len({key.casefold() for key in keys}) == 2
    assert all(key.startswith("one-") for key in keys)


def test_revision_two_index_migrates_and_rekeys_deterministically(
    simple_cds_gbff: Path, tmp_path: Path
) -> None:
    sources = tmp_path / "sources"
    _copy_source(simple_cds_gbff, sources, "one.gb")
    _copy_source(simple_cds_gbff, sources, "one.gbk")
    index = tmp_path / "cohort.gbidx"
    assert main(["index", "build", str(sources), str(index)]) == 0
    connection = sqlite3.connect(index)
    connection.execute("DROP INDEX uq_sources_sample_key_nocase")
    connection.execute("UPDATE sources SET sample_key = 'same'")
    connection.execute("PRAGMA user_version = 2")
    connection.execute("UPDATE metadata SET value = '2' WHERE key = 'schema_revision'")
    connection.commit()
    connection.close()
    (sources / "one.gb").unlink()
    (sources / "one.gbk").unlink()
    assert main(["index", "status", str(index)]) == 3
    migrated = migrate_index(index)
    assert migrated["schema_revision"] == "3"
    keys = [row[0] for row in sqlite3.connect(index).execute("SELECT sample_key FROM sources")]
    assert len(keys) == len({key.casefold() for key in keys}) == 2


def test_index_skipped_source_does_not_perturb_keys(
    simple_cds_gbff: Path, tmp_path: Path
) -> None:
    good_only = tmp_path / "good-only"
    _copy_source(simple_cds_gbff, good_only, "same.gb")
    first = tmp_path / "first.gbidx"
    assert main(["index", "build", str(good_only), str(first)]) == 0
    expected = sqlite3.connect(first).execute("SELECT sample_key FROM sources").fetchone()[0]
    mixed = tmp_path / "mixed"
    _copy_source(simple_cds_gbff, mixed, "same.gb")
    (mixed / "same.gbk").write_text("not GenBank", encoding="utf-8")
    second = tmp_path / "second.gbidx"
    assert main(["index", "build", str(mixed), str(second), "--on-error", "skip"]) == 3
    assert sqlite3.connect(second).execute("SELECT sample_key FROM sources").fetchone()[0] == expected


def test_batch_resume_rekeys_and_prunes_old_same_stem_tree(
    simple_cds_gbff: Path, tmp_path: Path
) -> None:
    sources = tmp_path / "sources"
    _copy_source(simple_cds_gbff, sources, "one.gb")
    run = tmp_path / "run"
    assert execute_batch([sources], command="validate", output_dir=run, tail=("--format", "json")).exit_code == 0
    _copy_source(simple_cds_gbff, sources, "one.gbk")
    result = execute_batch([sources], command="validate", output_dir=run, tail=("--format", "json"), resume=True)
    keys = {job["sample_key"] for job in result.manifest["jobs"] if job["status"] != "not_requested"}
    assert len(keys) == 2 and all(key.startswith("one-") for key in keys)
    assert not (run / "outputs" / "one").exists()


def test_batch_resume_is_idempotent_and_prunes_removed_artifacts(
    simple_cds_gbff: Path, tmp_path: Path
) -> None:
    sources = tmp_path / "sources"
    _copy_source(simple_cds_gbff, sources, "one.gb")
    _copy_source(simple_cds_gbff, sources, "two.gb")
    run = tmp_path / "run"
    execute_batch([sources], command="validate", output_dir=run, tail=("--format", "json"))
    (sources / "two.gb").unlink()
    first = execute_batch([sources], command="validate", output_dir=run, tail=("--format", "json"), resume=True)
    second = execute_batch([sources], command="validate", output_dir=run, tail=("--format", "json"), resume=True)
    removed = next(job for job in second.manifest["jobs"] if job["input"]["display_path"].endswith("two.gb"))
    assert removed["status"] == "not_requested"
    assert not (run / "outputs" / "two").exists()
    assert next(job for job in second.manifest["jobs"] if job["sample_key"] == "one")["resume_action"] == "reused_unchanged"
    assert first.manifest["jobs"][0]["status"] in {"succeeded", "not_requested"}


def test_batch_recovers_interrupted_first_execution(
    simple_cds_gbff: Path, tmp_path: Path
) -> None:
    run = tmp_path / "run"
    with patch("genbank_parser.batch._run_one", side_effect=KeyboardInterrupt), pytest.raises(KeyboardInterrupt):
        execute_batch([simple_cds_gbff], command="validate", output_dir=run, tail=("--format", "json"))
    progress = tmp_path / ".run.gbparse-inprogress"
    assert json.loads((progress / "batch-manifest.json").read_text())["jobs"][0]["status"] == "running"
    result = execute_batch([simple_cds_gbff], command="validate", output_dir=run, tail=("--format", "json"), resume=True)
    assert result.exit_code == 0
    assert result.manifest["jobs"][0]["resume_action"] == "replayed_running"
    assert not progress.exists()


def test_direct_query_rejects_cohort_only_sample_fields(
    simple_cds_gbff: Path, tmp_path: Path
) -> None:
    with pytest.raises(QueryExpressionError, match="cohort index"):
        query_features(simple_cds_gbff, 'sample_key == "simple_cds"')
    index = tmp_path / "cohort.gbidx"
    assert main(["index", "build", str(simple_cds_gbff), str(index)]) == 0
    assert query_index(index, 'sample_key == "simple_cds"')


def test_index_status_preserves_wal_mode(simple_cds_gbff: Path, tmp_path: Path) -> None:
    index = tmp_path / "cohort.gbidx"
    assert main(["index", "build", str(simple_cds_gbff), str(index)]) == 0
    connection = sqlite3.connect(index)
    assert connection.execute("PRAGMA journal_mode = WAL").fetchone()[0].casefold() == "wal"
    connection.close()
    assert inspect_index(index)["counts"]["sources"] == 1
    check = sqlite3.connect(index).execute("PRAGMA journal_mode").fetchone()[0]
    assert check.casefold() == "wal"


def test_index_reads_work_with_a_concurrent_wal_reader(simple_cds_gbff: Path, tmp_path: Path) -> None:
    index = tmp_path / "cohort.gbidx"
    assert main(["index", "build", str(simple_cds_gbff), str(index)]) == 0
    writer = sqlite3.connect(index)
    assert writer.execute("PRAGMA journal_mode = WAL").fetchone()[0].casefold() == "wal"
    writer.execute("BEGIN")
    assert inspect_index(index, check=True)["checks"]["quick_check"] == "ok"
    writer.rollback()
    writer.close()


def test_incomplete_index_fails_as_cli_input_error(tmp_path: Path) -> None:
    index = tmp_path / "incomplete.db"
    connection = sqlite3.connect(index)
    connection.execute("PRAGMA user_version = 3")
    connection.execute("CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL)")
    connection.execute("INSERT INTO metadata VALUES ('schema_version', 'gbparse.index.v1')")
    connection.execute("INSERT INTO metadata VALUES ('schema_revision', '3')")
    connection.commit()
    connection.close()
    assert main(["index", "status", str(index)]) == 3


def test_publication_cleanup_failure_is_post_commit_warning(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    destination = tmp_path / "result.txt"
    destination.write_text("old", encoding="utf-8")
    original_unlink = Path.unlink
    calls = {"backup": 0}

    def flaky_unlink(path: Path, *args: object, **kwargs: object) -> None:
        if path.name.endswith(".backup"):
            calls["backup"] += 1
            if calls["backup"] > 1:
                raise OSError("locked backup")
        original_unlink(path, *args, **kwargs)

    monkeypatch.setattr(Path, "unlink", flaky_unlink)
    with pytest.warns(RuntimeWarning, match="retained backup"):
        publish_file_set({destination: "new"}, force=True)
    assert destination.read_text(encoding="utf-8") == "new"
    assert list(tmp_path.glob(".*.backup"))


def test_publication_rollback_failure_preserves_recovery_artifacts(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    destination = tmp_path / "result.txt"
    second = tmp_path / "second.txt"
    destination.write_text("old", encoding="utf-8")
    second.write_text("old2", encoding="utf-8")
    staged_one = tmp_path / "one.tmp"
    staged_two = tmp_path / "two.tmp"
    staged_one.write_text("new", encoding="utf-8")
    staged_two.write_text("new2", encoding="utf-8")
    real_replace = __import__("os").replace

    def fail_install(source: str | Path, target: str | Path) -> None:
        if Path(target) == second and Path(source).suffix == ".tmp":
            raise OSError("install failed")
        if Path(target) == destination and ".backup" in Path(source).name:
            raise OSError("restore failed")
        real_replace(source, target)

    monkeypatch.setattr("genbank_parser.cli_io.os.replace", fail_install)
    with pytest.raises(OutputRecoveryError) as caught:
        publish_staged_files(((staged_one, destination), (staged_two, second)), force=True)
    assert str(destination) in str(caught.value)
    assert list(tmp_path.glob(".*.backup"))
    assert staged_two.exists()


def test_directory_tree_rollback_failure_preserves_recovery_artifacts(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    destination = tmp_path / "run"
    destination.mkdir()
    (destination / "old.txt").write_text("old", encoding="utf-8")
    staging = tmp_path / "staging"
    staging.mkdir()
    (staging / "new.txt").write_text("new", encoding="utf-8")
    real_replace = __import__("os").replace

    def fail_restore(source: str | Path, target: str | Path) -> None:
        if Path(target) == destination and Path(source) == staging:
            raise OSError("install failed")
        if Path(target) == destination and ".backup." in Path(source).name:
            raise OSError("restore failed")
        real_replace(source, target)

    monkeypatch.setattr("genbank_parser.cli_io.os.replace", fail_restore)
    with pytest.raises(OutputRecoveryError):
        publish_directory_tree(staging, destination, force=True)
    assert staging.exists()
    assert list(tmp_path.glob(".run.backup.*"))


def test_batch_unknown_child_exit_is_nonzero(simple_cds_gbff: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    completed = __import__("subprocess").CompletedProcess([], 5, "", "")
    monkeypatch.setattr("genbank_parser.batch.subprocess.run", lambda *args, **kwargs: completed)
    result = execute_batch([simple_cds_gbff], command="validate", output_dir=tmp_path / "run")
    assert result.exit_code == 3
    assert result.failed_count == 1
    assert result.manifest["jobs"][0]["exit_code"] == 5


def test_batch_stop_treats_threshold_failure_as_stop(
    duplicate_locus_gbff: Path, tmp_path: Path
) -> None:
    sources = tmp_path / "sources"
    _copy_source(duplicate_locus_gbff, sources, "a.gb")
    _copy_source(duplicate_locus_gbff, sources, "b.gb")
    result = execute_batch(
        [sources],
        command="validate",
        output_dir=tmp_path / "run",
        tail=("--format", "json", "--fail-on", "warning"),
        on_error="stop",
    )
    assert result.exit_code == 1
    assert result.failed_count == 1
    assert sum(job["status"] == "threshold_failed" for job in result.manifest["jobs"]) == 1
    assert sum(job["status"] == "pending" for job in result.manifest["jobs"]) == 1


def test_batch_manifest_schema_separates_status_and_resume_action(
    simple_cds_gbff: Path, tmp_path: Path
) -> None:
    run = tmp_path / "run"
    execute_batch([simple_cds_gbff], command="validate", output_dir=run)
    result = execute_batch([simple_cds_gbff], command="validate", output_dir=run, resume=True)
    job = result.manifest["jobs"][0]
    assert job["status"] == "succeeded"
    assert job["resume_action"] == "reused_unchanged"
    assert job["status"] != "skipped_unchanged"


def test_source_label_digest_binding_rejects_mutated_snapshot(simple_cds_gbff: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("GBPARSE_SOURCE_LABEL", "logical.gb")
    monkeypatch.setenv("GBPARSE_SOURCE_SHA256", "0" * 64)
    with pytest.raises(GenBankInputError, match="digest mismatch"):
        read_genbank(simple_cds_gbff)


def test_discovery_metadata_is_invariant_to_input_permutation(simple_cds_gbff: Path, tmp_path: Path) -> None:
    root = tmp_path / "cohort"
    source = _copy_source(simple_cds_gbff, root, "one.gb")
    first = discover_inputs([root, source])
    second = discover_inputs([source, root])
    assert first == second
