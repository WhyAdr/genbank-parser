"""Tests for fail-closed mobilome evidence-resource loading."""

import json
from importlib import resources
from pathlib import Path

import pytest
import yaml

from genbank_parser.mobilome import load_mobilome_database
from genbank_parser.mobilome.database import MobilomeDatabaseError


def _copy_custom_database(destination: Path) -> None:
    source = resources.files("genbank_parser").joinpath("data", "mobilome")
    for name in ("markers.yaml", "provenance.yaml", "inference.yaml"):
        (destination / name).write_bytes(source.joinpath(name).read_bytes())


def test_packaged_database_is_versioned_hashed_and_schema_valid() -> None:
    database = load_mobilome_database()
    schema_path = resources.files("genbank_parser").joinpath(
        "data", "mobilome", "report.schema.json"
    )
    schema = json.loads(schema_path.read_text(encoding="utf-8"))

    assert database.catalog_version == "1.0.0"
    assert database.inference_version == "1.0.0"
    assert database.database_source == "packaged"
    assert database.source_paths == ()
    assert [resource.name for resource in database.resources] == [
        "markers.yaml",
        "provenance.yaml",
        "inference.yaml",
        "report.schema.json",
    ]
    assert all(len(resource.sha256) == 64 for resource in database.resources)
    assert schema["$schema"] == "https://json-schema.org/draft/2020-12/schema"


def test_complete_custom_database_is_loaded_as_one_resource_set(tmp_path: Path) -> None:
    _copy_custom_database(tmp_path)

    database = load_mobilome_database(tmp_path)

    assert database.database_source == "custom"
    assert tuple(path.name for path in database.source_paths) == (
        "markers.yaml",
        "provenance.yaml",
        "inference.yaml",
    )


def test_database_rejects_incomplete_or_invalid_custom_resources(
    tmp_path: Path,
) -> None:
    _copy_custom_database(tmp_path)
    (tmp_path / "extra.txt").write_text("unexpected", encoding="utf-8")
    with pytest.raises(MobilomeDatabaseError, match="exactly"):
        load_mobilome_database(tmp_path)

    (tmp_path / "extra.txt").unlink()
    markers_path = tmp_path / "markers.yaml"
    payload = yaml.safe_load(markers_path.read_text(encoding="utf-8"))
    payload["markers"][0]["matchers"][0]["strength"] = 4
    markers_path.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")
    with pytest.raises(MobilomeDatabaseError, match="strength"):
        load_mobilome_database(tmp_path)
