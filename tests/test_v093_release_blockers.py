"""Release-blocker regressions for the v0.9.3 hardening plan."""

from __future__ import annotations

from pathlib import Path

import pytest

from genbank_parser.batch import BatchUsageError, execute_batch


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
