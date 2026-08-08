"""Test validation findings, pseudogene tolerance, codon start, and duplicate locus tags."""
from pathlib import Path

from genbank_parser.validate import validate


def test_validate_clean_fixture(simple_cds_gbff: Path) -> None:
    findings = validate(simple_cds_gbff, json_mode=True)
    # Should have no ERROR findings
    errors = [f for f in findings if f.severity == 'ERROR']
    assert len(errors) == 0


def test_validate_duplicate_locus_tag(duplicate_locus_gbff: Path) -> None:
    findings = validate(duplicate_locus_gbff, json_mode=True)
    dup_findings = [f for f in findings if f.code == 'DUPLICATE_LOCUS_TAG']
    assert len(dup_findings) >= 1
    assert dup_findings[0].locus_tag == 'DUP_TAG_001'
    assert dup_findings[0].severity == 'WARNING'


def test_validate_special_features_tolerance(special_cds_gbff: Path) -> None:
    findings = validate(special_cds_gbff, json_mode=True)
    # Pseudogene should NOT produce missing translation warning
    pseudo_findings = [f for f in findings if f.locus_tag == 'PSEUDO_001']
    assert not any(f.code == 'MISSING_TRANSLATION' for f in pseudo_findings)

    # Partial feature should not fail translation mismatch
    part_findings = [f for f in findings if f.locus_tag == 'PARTIAL_001']
    assert not any(f.code == 'CDS_TRANSLATION_MISMATCH' for f in part_findings)
