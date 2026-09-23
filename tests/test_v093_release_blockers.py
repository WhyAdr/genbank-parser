"""Release-blocker regressions for the v0.9.3 hardening plan."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from genbank_parser.batch import BatchUsageError, execute_batch


def _copy_source(source: Path, directory: Path, name: str) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    target = directory / name
    shutil.copyfile(source, target)
    return target


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
