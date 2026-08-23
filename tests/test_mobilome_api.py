"""Public native mobilome API contracts."""

from pathlib import Path

import pytest

import genbank_parser
from genbank_parser.mobilome import analyze_mobilome, load_mobilome_database
from genbank_parser.mobilome.models import MobilomeParameterError
from genbank_parser.mobilome.report import (
    serialize_mobilome_report,
    write_mobilome_report,
)


def test_public_api_is_pure_and_root_api_remains_small() -> None:
    report = analyze_mobilome(Path("tests/fixtures/mobilome_inventory.gb"))

    assert report.source_file == "mobilome_inventory.gb"
    assert report.catalog_version == load_mobilome_database().catalog_version
    assert not hasattr(genbank_parser, "analyze_mobilome")


def test_invalid_public_options_fail_before_input_work() -> None:
    missing = Path("does-not-exist.gbff")
    with pytest.raises(MobilomeParameterError, match="include"):
        analyze_mobilome(missing, include="invalid")
    with pytest.raises(MobilomeParameterError, match="min_evidence"):
        analyze_mobilome(missing, min_evidence=4)
    with pytest.raises(MobilomeParameterError, match="format_type"):
        serialize_mobilome_report(
            analyze_mobilome(Path("tests/fixtures/mobilome_inventory.gb")), "invalid"
        )
    report = analyze_mobilome(Path("tests/fixtures/mobilome_inventory.gb"))
    with pytest.raises(MobilomeParameterError, match="format_type"):
        write_mobilome_report(
            report,
            source_path=missing,
            database_paths=(),
            format_type="invalid",
        )
