"""Test query search, region extraction with rebase, and annotation diffing."""
from pathlib import Path

from genbank_parser.diff import diff_annotations
from genbank_parser.query import search_features
from genbank_parser.region import extract_region


def test_search_features(simple_cds_gbff: Path) -> None:
    results_gene = search_features(simple_cds_gbff, gene="testA", ftype="CDS", format_type="json")
    assert len(results_gene) == 1
    assert results_gene[0]['locus_tag'] == 'TEST_001'

    results_ko = search_features(simple_cds_gbff, ko="K00844", format_type="json")
    assert len(results_ko) == 1
    assert results_ko[0]['locus_tag'] == 'TEST_001'


def test_region_extraction_and_rebase(simple_cds_gbff: Path) -> None:
    sub_rec = extract_region(simple_cds_gbff, locus_tag="TEST_002", rebase=True)
    assert len(sub_rec.seq) == 30
    assert len(sub_rec.features) >= 1

    rebased_feat = sub_rec.features[0]
    assert int(rebased_feat.location.start) == 0
    assert int(rebased_feat.location.end) == 30


def test_diff_annotations(simple_cds_gbff: Path, special_cds_gbff: Path) -> None:
    diff_res = diff_annotations(simple_cds_gbff, special_cds_gbff, format_type="json")
    assert diff_res['old_cds_count'] == 2
    assert diff_res['new_cds_count'] >= 3
    assert diff_res['added_cds'] >= 1
