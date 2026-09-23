"""Release-blocker regressions for the v0.9.3 hardening plan."""

from __future__ import annotations

import hashlib
import os
import shutil
import sqlite3
import subprocess
from importlib import resources
from pathlib import Path

import pytest

from genbank_parser import read_genbank
from genbank_parser.batch import BatchUsageError, execute_batch
from genbank_parser.cli import main
from genbank_parser.cli_io import InputError, OutputError, OutputRecoveryError
from genbank_parser.index import build_index, migrate_index
from genbank_parser.index.schema import SCHEMA_SQL
from genbank_parser.io import iter_genbank

# Pinned historical DDL from 9db07e0c7c9adee4675f564f6bf281e5ad7d457c.
_LEGACY_V1_SQL = """
CREATE TABLE metadata (
    key TEXT PRIMARY KEY,
    value TEXT NOT NULL
);
CREATE TABLE sources (
    source_pk INTEGER PRIMARY KEY,
    display_path TEXT NOT NULL UNIQUE,
    resolved_path_at_index TEXT NOT NULL,
    sample_key TEXT NOT NULL,
    sha256 TEXT NOT NULL,
    size_bytes INTEGER NOT NULL,
    mtime_ns INTEGER NOT NULL,
    compressed INTEGER NOT NULL CHECK (compressed IN (0, 1)),
    record_count INTEGER NOT NULL,
    feature_count INTEGER NOT NULL
);
CREATE TABLE records (
    record_pk INTEGER PRIMARY KEY,
    source_pk INTEGER NOT NULL REFERENCES sources(source_pk) ON DELETE CASCADE,
    record_index INTEGER NOT NULL,
    record_id TEXT NOT NULL,
    record_name TEXT NOT NULL,
    description TEXT NOT NULL,
    length INTEGER NOT NULL,
    topology TEXT,
    molecule_type TEXT,
    organism TEXT,
    strain TEXT,
    gc_percent REAL,
    UNIQUE(source_pk, record_index)
);
CREATE TABLE features (
    feature_pk INTEGER PRIMARY KEY,
    record_pk INTEGER NOT NULL REFERENCES records(record_pk) ON DELETE CASCADE,
    feature_index INTEGER NOT NULL,
    type TEXT NOT NULL,
    locus_tag TEXT,
    gene TEXT,
    product TEXT,
    protein_id TEXT,
    start INTEGER NOT NULL,
    end INTEGER NOT NULL,
    strand INTEGER,
    biological_length INTEGER NOT NULL,
    partial INTEGER NOT NULL CHECK (partial IN (0, 1)),
    pseudo INTEGER NOT NULL CHECK (pseudo IN (0, 1)),
    UNIQUE(record_pk, feature_index)
);
CREATE TABLE segments (
    feature_pk INTEGER NOT NULL REFERENCES features(feature_pk) ON DELETE CASCADE,
    ordinal INTEGER NOT NULL,
    start INTEGER NOT NULL,
    end INTEGER NOT NULL,
    PRIMARY KEY(feature_pk, ordinal)
);
CREATE TABLE qualifiers (
    feature_pk INTEGER NOT NULL REFERENCES features(feature_pk) ON DELETE CASCADE,
    key TEXT NOT NULL,
    ordinal INTEGER NOT NULL,
    value TEXT NOT NULL,
    PRIMARY KEY(feature_pk, key, ordinal)
);
CREATE TABLE xrefs (
    feature_pk INTEGER NOT NULL REFERENCES features(feature_pk) ON DELETE CASCADE,
    namespace TEXT NOT NULL,
    value TEXT NOT NULL,
    source_field TEXT NOT NULL,
    PRIMARY KEY(feature_pk, namespace, value, source_field)
);
CREATE INDEX idx_sources_sha256 ON sources(sha256);
CREATE INDEX idx_records_id ON records(record_id);
CREATE INDEX idx_records_name ON records(record_name);
CREATE INDEX idx_records_organism ON records(organism);
CREATE INDEX idx_features_type ON features(type);
CREATE INDEX idx_features_locus_tag ON features(locus_tag);
CREATE INDEX idx_features_gene ON features(gene);
CREATE INDEX idx_features_protein_id ON features(protein_id);
CREATE INDEX idx_features_coordinates ON features(record_pk, start, end);
CREATE INDEX idx_xrefs_namespace_value ON xrefs(namespace, value);
CREATE INDEX idx_qualifiers_key_value ON qualifiers(key, value);
"""


def _copy_source(source: Path, directory: Path, name: str) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    target = directory / name
    shutil.copyfile(source, target)
    return target


def _make_authentic_legacy_index(
    source: Path, directory: Path, revision: int
) -> tuple[Path, Path]:
    """Build a pinned historical-D DL index with real parsed fixture rows."""

    fresh = directory / "fresh.gbidx"
    assert build_index([source], fresh).exit_code == 0
    legacy = directory / f"legacy-v{revision}.gbidx"
    connection = sqlite3.connect(legacy)
    try:
        if revision == 1:
            connection.executescript(_LEGACY_V1_SQL)
        elif revision == 2:
            # The v0.9.1 table definitions are the v0.9.2 definitions before
            # the explicit v3 source/sample-key indexes were added.
            connection.executescript(SCHEMA_SQL)
            connection.execute("DROP INDEX uq_sources_source_identity")
            connection.execute("DROP INDEX uq_sources_sample_key_nocase")
        else:
            raise AssertionError(revision)
        connection.execute("ATTACH DATABASE ? AS source_db", (str(fresh),))
        connection.execute("INSERT INTO metadata SELECT * FROM source_db.metadata")
        connection.execute(
            "UPDATE metadata SET value = ? WHERE key = 'schema_revision'",
            (str(revision),),
        )
        connection.execute(
            "UPDATE metadata SET value = ? WHERE key = 'gbparse_version'",
            ("0.9.0" if revision == 1 else "0.9.1",),
        )
        connection.execute(
            "UPDATE metadata SET value = ? WHERE key = 'query_semantics_version'",
            ("gbparse.query.v1" if revision == 1 else "gbparse.query.v2",),
        )
        if revision == 1:
            connection.execute(
                "INSERT INTO sources SELECT source_pk, display_path, resolved_path_at_index, sample_key, sha256, size_bytes, mtime_ns, compressed, record_count, feature_count FROM source_db.sources"
            )
        else:
            connection.execute("INSERT INTO sources SELECT * FROM source_db.sources")
        for table in ("records", "features", "segments", "qualifiers"):
            connection.execute(f"INSERT INTO {table} SELECT * FROM source_db.{table}")
        if revision == 1:
            connection.execute(
                "INSERT INTO xrefs SELECT feature_pk, namespace, value, source_field FROM source_db.xrefs"
            )
        else:
            connection.execute("INSERT INTO xrefs SELECT * FROM source_db.xrefs")
        connection.execute(f"PRAGMA user_version = {revision}")
        connection.commit()
        connection.execute("DETACH DATABASE source_db")
    finally:
        connection.close()
    return fresh, legacy


def _index_signature(path: Path) -> dict[str, object]:
    connection = sqlite3.connect(path)
    try:
        tables = tuple(
            str(row[0])
            for row in connection.execute(
                "SELECT name FROM sqlite_master WHERE type = 'table' ORDER BY name"
            )
        )
        master = tuple(
            tuple(row)
            for row in connection.execute(
                "SELECT type, name, tbl_name, sql FROM sqlite_master "
                "WHERE type IN ('table', 'index') ORDER BY type, name"
            )
        )
        table_info = {
            table: tuple(tuple(row) for row in connection.execute(f"PRAGMA table_info({table})"))
            for table in tables
        }
        index_list = {
            table: tuple(tuple(row) for row in connection.execute(f"PRAGMA index_list({table})"))
            for table in tables
        }
        index_xinfo = {
            str(row[1]): tuple(tuple(item) for item in connection.execute(f"PRAGMA index_xinfo('{row[1]}')"))
            for table in tables
            for row in connection.execute(f"PRAGMA index_list({table})")
        }
        representative = tuple(
            connection.execute(
                "SELECT s.source_identity, s.sample_key, r.record_id, f.feature_index, "
                "f.locus_tag, x.namespace, x.ordinal, x.value "
                "FROM sources AS s JOIN records AS r USING (source_pk) "
                "JOIN features AS f USING (record_pk) "
                "LEFT JOIN xrefs AS x USING (feature_pk) "
                "ORDER BY s.source_identity, r.record_index, f.feature_index, x.namespace, x.ordinal"
            ).fetchall()
        )
        counts = tuple(
            (table, int(connection.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]))
            for table in ("sources", "records", "features", "segments", "qualifiers", "xrefs")
        )
        return {
            "master": master,
            "table_info": table_info,
            "index_list": index_list,
            "index_xinfo": index_xinfo,
            "foreign_key_check": tuple(connection.execute("PRAGMA foreign_key_check").fetchall()),
            "integrity_check": connection.execute("PRAGMA integrity_check").fetchone()[0],
            "representative": representative,
            "counts": counts,
        }
    finally:
        connection.close()


def _index_revision_metadata(path: Path) -> tuple[int, str]:
    connection = sqlite3.connect(path)
    try:
        return (
            int(connection.execute("PRAGMA user_version").fetchone()[0]),
            str(
                connection.execute(
                    "SELECT value FROM metadata WHERE key = 'schema_revision'"
                ).fetchone()[0]
            ),
        )
    finally:
        connection.close()


@pytest.mark.parametrize(
    ("command", "flag"),
    (
        ("discover", "--rules"),
        ("meor", "--markers"),
        ("meor", "--pathways"),
        ("mobilome", "--database-dir"),
    ),
)
@pytest.mark.parametrize("equals_form", (False, True))
def test_batch_rejects_custom_resources_before_recovery_tree(
    simple_cds_gbff: Path,
    tmp_path: Path,
    command: str,
    flag: str,
    equals_form: bool,
) -> None:
    resource = tmp_path / "custom-resource"
    if command == "mobilome":
        resource.mkdir()
    else:
        resource.write_text("custom", encoding="utf-8")
    output = tmp_path / "run"
    token = f"{flag}={resource}" if equals_form else flag
    tail = (token,) if equals_form else (flag, str(resource))

    with pytest.raises(BatchUsageError, match="auxiliary-resource snapshot"):
        execute_batch(
            [simple_cds_gbff],
            command=command,
            output_dir=output,
            tail=tail,
        )
    assert not output.exists()
    assert not (tmp_path / ".run.gbparse-inprogress").exists()


def test_batch_force_cannot_replace_auxiliary_input_tree(
    simple_cds_gbff: Path, tmp_path: Path
) -> None:
    resource = tmp_path / "custom-database"
    resource.mkdir()
    marker = resource / "markers.yaml"
    marker.write_text("preserve", encoding="utf-8")

    with pytest.raises(BatchUsageError, match="--database-dir"):
        execute_batch(
            [simple_cds_gbff],
            command="mobilome",
            output_dir=resource,
            tail=(f"--database-dir={resource}",),
            force=True,
        )
    assert marker.read_text(encoding="utf-8") == "preserve"


def test_direct_commands_retain_custom_resource_support(
    simple_cds_gbff: Path, tmp_path: Path
) -> None:
    package = resources.files("genbank_parser")
    rules = tmp_path / "rules.yaml"
    rules.write_bytes(package.joinpath("rulesets", "mobilome.yaml").read_bytes())
    assert (
        main(
            [
                "discover",
                str(simple_cds_gbff),
                "--rules",
                str(rules),
                "--format",
                "json",
                "--output",
                str(tmp_path / "discover.json"),
            ]
        )
        == 0
    )

    meor_dir = tmp_path / "meor"
    meor_dir.mkdir()
    meor_markers = meor_dir / "markers.yaml"
    meor_pathways = meor_dir / "pathways.yaml"
    meor_markers.write_bytes(
        package.joinpath("data", "meor", "markers.yaml").read_bytes()
    )
    meor_pathways.write_bytes(
        package.joinpath("data", "meor", "pathways.yaml").read_bytes()
    )
    assert (
        main(
            [
                "meor",
                str(simple_cds_gbff),
                "--markers",
                str(meor_markers),
                "--pathways",
                str(meor_pathways),
                "--format",
                "json",
                "--output",
                str(tmp_path / "meor.json"),
            ]
        )
        == 0
    )

    mobilome_dir = tmp_path / "mobilome"
    mobilome_dir.mkdir()
    for name in ("markers.yaml", "provenance.yaml", "inference.yaml"):
        target = mobilome_dir / name
        target.write_bytes(
            package.joinpath("data", "mobilome", name).read_bytes()
        )
    assert (
        main(
            [
                "mobilome",
                str(simple_cds_gbff),
                "--database-dir",
                str(mobilome_dir),
                "--format",
                "json",
                "--output",
                str(tmp_path / "mobilome.json"),
            ]
        )
        == 0
    )


def test_stop_barrier_runs_changed_job_before_reused_failure(
    simple_cds_gbff: Path, duplicate_locus_gbff: Path, tmp_path: Path
) -> None:
    sources = tmp_path / "sources"
    first = _copy_source(simple_cds_gbff, sources, "a.gb")
    (sources / "b.gb").write_text("not GenBank", encoding="utf-8")
    _copy_source(simple_cds_gbff, sources, "c.gb")
    run = tmp_path / "run"

    first_result = execute_batch(
        [sources], command="validate", output_dir=run, on_error="stop"
    )
    assert [job["status"] for job in first_result.manifest["jobs"]] == [
        "succeeded",
        "failed",
        "pending",
    ]

    shutil.copyfile(duplicate_locus_gbff, first)
    resumed = execute_batch(
        [sources], command="validate", output_dir=run, on_error="stop", resume=True
    )
    jobs = {job["sample_key"]: job for job in resumed.manifest["jobs"]}
    assert jobs["a"]["resume_action"] == "rerun_changed"
    assert jobs["a"]["status"] == "succeeded"
    assert jobs["b"]["resume_action"] == "reused_unchanged"
    assert jobs["b"]["status"] == "failed"
    assert jobs["c"]["resume_action"] == "deferred_after_stop"
    assert jobs["c"]["status"] == "pending"


def test_stop_barrier_handles_threshold_failure_before_changed_job(
    simple_cds_gbff: Path, duplicate_locus_gbff: Path, tmp_path: Path
) -> None:
    sources = tmp_path / "sources"
    first = _copy_source(simple_cds_gbff, sources, "a.gb")
    _copy_source(duplicate_locus_gbff, sources, "b.gb")
    _copy_source(simple_cds_gbff, sources, "c.gb")
    run = tmp_path / "run"
    tail = ("--format", "json", "--fail-on", "warning")

    first_result = execute_batch(
        [sources],
        command="validate",
        output_dir=run,
        tail=tail,
        on_error="stop",
    )
    assert [job["status"] for job in first_result.manifest["jobs"]] == [
        "succeeded",
        "threshold_failed",
        "pending",
    ]

    shutil.copyfile(duplicate_locus_gbff, first)
    resumed = execute_batch(
        [sources],
        command="validate",
        output_dir=run,
        tail=tail,
        on_error="stop",
        resume=True,
    )
    jobs = {job["sample_key"]: job for job in resumed.manifest["jobs"]}
    assert jobs["a"]["resume_action"] == "rerun_changed"
    assert jobs["a"]["status"] == "threshold_failed"
    assert jobs["b"]["resume_action"] == "reused_unchanged"
    assert jobs["b"]["status"] == "threshold_failed"
    assert jobs["c"]["resume_action"] == "deferred_after_stop"
    assert jobs["c"]["status"] == "pending"


def test_batch_publication_failure_preserves_inprogress_tree(
    simple_cds_gbff: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    run = tmp_path / "run"
    execute_batch([simple_cds_gbff], command="validate", output_dir=run)
    before = (run / "batch-manifest.json").read_bytes()
    progress = tmp_path / ".run.gbparse-inprogress"
    real_replace = os.replace

    def fail_final(source: str | Path, target: str | Path) -> None:
        if Path(source) == progress and Path(target) == run:
            raise OSError("final publication failed")
        real_replace(source, target)

    monkeypatch.setattr("genbank_parser.cli_io.os.replace", fail_final)
    with pytest.raises(OutputRecoveryError, match="durable staging retained") as caught:
        execute_batch([simple_cds_gbff], command="validate", output_dir=run, force=True)

    assert run.joinpath("batch-manifest.json").read_bytes() == before
    assert progress.joinpath("batch-manifest.json").is_file()
    assert progress.joinpath("logs").is_dir()
    assert list(progress.joinpath("outputs").rglob("*"))
    assert str(progress) in str(caught.value)


def test_resume_copy_interruption_is_self_recoverable(
    simple_cds_gbff: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    run = tmp_path / "run"
    execute_batch([simple_cds_gbff], command="validate", output_dir=run)
    before = {path.relative_to(run): path.read_bytes() for path in run.rglob("*") if path.is_file()}
    real_copy2 = shutil.copy2
    real_copytree = shutil.copytree

    def interrupt_copytree(source: str | Path, target: str | Path, *args: object, **kwargs: object) -> None:
        source_path = Path(source)
        target_path = Path(target)
        first = next(path for path in source_path.rglob("*") if path.is_file())
        copied = target_path / first.relative_to(source_path)
        copied.parent.mkdir(parents=True, exist_ok=True)
        real_copy2(first, copied)
        raise OSError("resume copy interrupted")

    monkeypatch.setattr("genbank_parser.batch.shutil.copytree", interrupt_copytree)
    with pytest.raises(OutputError, match="durable batch resume tree"):
        execute_batch([simple_cds_gbff], command="validate", output_dir=run, resume=True)

    assert not (tmp_path / ".run.gbparse-inprogress").exists()
    assert not list(tmp_path.glob(".run.gbparse-resume-copy.*"))
    assert {path.relative_to(run): path.read_bytes() for path in run.rglob("*") if path.is_file()} == before

    monkeypatch.setattr("genbank_parser.batch.shutil.copytree", real_copytree)
    resumed = execute_batch([simple_cds_gbff], command="validate", output_dir=run, resume=True)
    assert resumed.exit_code == 0
    assert not (tmp_path / ".run.gbparse-inprogress").exists()


def test_iter_genbank_feature_indices_are_document_global(
    multi_record_circular_gbff: Path,
) -> None:
    streamed = [
        feature.feature_index
        for record in iter_genbank(multi_record_circular_gbff)
        for feature in record.features
    ]
    materialized = [
        feature.feature_index
        for record in read_genbank(multi_record_circular_gbff).records
        for feature in record.features
    ]
    assert streamed == materialized == list(range(1, len(streamed) + 1))


def _assert_authentic_legacy_index_migration(
    simple_cds_gbff: Path, tmp_path: Path, revision: int
) -> None:
    fresh, legacy = _make_authentic_legacy_index(simple_cds_gbff, tmp_path, revision)
    immutable = tmp_path / f"legacy-v{revision}-immutable.gbidx"
    target = tmp_path / f"legacy-v{revision}-target.gbidx"
    shutil.copyfile(legacy, immutable)
    shutil.copyfile(immutable, target)
    before_digest = hashlib.sha256(immutable.read_bytes()).hexdigest()
    before_metadata = _index_revision_metadata(immutable)

    migrated = migrate_index(target)
    assert migrated["schema_revision"] == "3"
    assert _index_signature(target) == _index_signature(fresh)
    assert hashlib.sha256(immutable.read_bytes()).hexdigest() == before_digest
    assert _index_revision_metadata(immutable) == before_metadata


def test_authentic_revision1_to_revision3_schema_equivalence(
    simple_cds_gbff: Path, tmp_path: Path
) -> None:
    _assert_authentic_legacy_index_migration(simple_cds_gbff, tmp_path, 1)


def test_authentic_revision2_to_revision3_schema_equivalence(
    simple_cds_gbff: Path, tmp_path: Path
) -> None:
    _assert_authentic_legacy_index_migration(simple_cds_gbff, tmp_path, 2)


def test_failed_migration_does_not_advance_schema_revision(
    simple_cds_gbff: Path, tmp_path: Path
) -> None:
    _fresh, legacy = _make_authentic_legacy_index(simple_cds_gbff, tmp_path, 2)
    target = tmp_path / "tampered.gbidx"
    shutil.copyfile(legacy, target)
    connection = sqlite3.connect(target)
    connection.execute(
        "CREATE UNIQUE INDEX uq_sources_sample_key_nocase ON sources(sample_key)"
    )
    connection.commit()
    connection.close()

    with pytest.raises(InputError, match="invalid existing index"):
        migrate_index(target)
    connection = sqlite3.connect(target)
    try:
        assert connection.execute("PRAGMA user_version").fetchone()[0] == 2
        assert connection.execute(
            "SELECT value FROM metadata WHERE key = 'schema_revision'"
        ).fetchone()[0] == "2"
        index_sql = connection.execute(
            "SELECT sql FROM sqlite_master WHERE name = 'uq_sources_sample_key_nocase'"
        ).fetchone()[0]
        assert "COLLATE NOCASE" not in index_sql.upper()
    finally:
        connection.close()


def test_migrate_current_wal_preserves_journal_mode(
    simple_cds_gbff: Path, tmp_path: Path
) -> None:
    index = tmp_path / "wal.gbidx"
    assert build_index([simple_cds_gbff], index).exit_code == 0
    connection = sqlite3.connect(index)
    try:
        assert connection.execute("PRAGMA journal_mode = WAL").fetchone()[0].casefold() == "wal"
    finally:
        connection.close()

    assert main(["index", "migrate", str(index)]) == 0
    connection = sqlite3.connect(index)
    try:
        assert connection.execute("PRAGMA journal_mode").fetchone()[0].casefold() == "wal"
    finally:
        connection.close()


def test_migrate_with_active_wal_reader_is_noop(
    simple_cds_gbff: Path, tmp_path: Path
) -> None:
    index = tmp_path / "active-wal.gbidx"
    assert build_index([simple_cds_gbff], index).exit_code == 0
    reader = sqlite3.connect(index)
    try:
        assert reader.execute("PRAGMA journal_mode = WAL").fetchone()[0].casefold() == "wal"
        reader.execute("BEGIN")
        assert main(["index", "migrate", str(index)]) == 0
    finally:
        reader.rollback()
        reader.close()


def test_migrate_corrupt_database_returns_exit_3(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    corrupt = tmp_path / "corrupt.gbidx"
    corrupt.write_bytes(b"not a sqlite database")
    assert main(["index", "migrate", str(corrupt)]) == 3
    assert "Traceback" not in capsys.readouterr().err


def test_update_corrupt_database_returns_exit_3(
    simple_cds_gbff: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    corrupt = tmp_path / "corrupt.gbidx"
    corrupt.write_bytes(b"not a sqlite database")
    assert main(["index", "update", str(corrupt), str(simple_cds_gbff)]) == 3
    assert "Traceback" not in capsys.readouterr().err


@pytest.mark.parametrize("relative_source", (True, False))
def test_manifest_argv_is_replayable_from_recorded_working_directory(
    simple_cds_gbff: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    relative_source: bool,
) -> None:
    monkeypatch.chdir(tmp_path)
    source = Path("inputs") / "sample.gb"
    _copy_source(simple_cds_gbff, source.parent, source.name)
    source_argument = source if relative_source else source.resolve()
    run = Path("run")

    result = execute_batch(
        [source_argument],
        command="validate",
        output_dir=run,
        tail=("--format", "json"),
    )
    job = result.manifest["jobs"][0]
    argv = job["argv"]
    assert str(source.resolve()) in argv
    assert str((tmp_path / "run" / "outputs" / "sample" / "result.json").resolve()) in argv
    assert ".gbparse-input." not in " ".join(argv)
    assert ".gbparse-inprogress" not in " ".join(argv)
    assert ".work." not in " ".join(argv)

    completed = subprocess.run(
        argv,
        cwd=result.manifest["runner"]["working_directory"],
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    assert (tmp_path / "run" / "outputs" / "sample" / "result.json").is_file()
    assert not (tmp_path / "outputs").exists()
