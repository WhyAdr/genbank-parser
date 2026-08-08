"""Pytest configuration and shared fixtures for genbank_parser."""
from pathlib import Path
import pytest

FIXTURES_DIR = Path(__file__).resolve().parent / "fixtures"


@pytest.fixture
def fixtures_dir() -> Path:
    return FIXTURES_DIR


@pytest.fixture
def simple_cds_gbff() -> Path:
    return FIXTURES_DIR / "simple_cds.gb"


@pytest.fixture
def compound_joined_gbff() -> Path:
    return FIXTURES_DIR / "compound_joined.gb"


@pytest.fixture
def special_cds_gbff() -> Path:
    return FIXTURES_DIR / "special_cds.gb"


@pytest.fixture
def multi_record_circular_gbff() -> Path:
    return FIXTURES_DIR / "multi_record_circular.gb"


@pytest.fixture
def duplicate_locus_gbff() -> Path:
    return FIXTURES_DIR / "duplicate_locus.gb"


@pytest.fixture
def real_c14_gbff() -> Path:
    p = Path(__file__).resolve().parent.parent / "C14-NMZ.gbff"
    return p if p.exists() else None
