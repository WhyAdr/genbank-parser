"""Test typed data models and biological length/span invariants."""
from pathlib import Path

from genbank_parser.io import read_genbank
from genbank_parser.model import GenBankDocument, GenBankFeature, GenBankRecord


def test_simple_feature_length_and_span(simple_cds_gbff: Path) -> None:
    doc = read_genbank(simple_cds_gbff)
    assert len(doc.records) == 1
    rec = doc.records[0]

    cdss = rec.cds_features
    assert len(cdss) == 2

    # Forward simple CDS
    fwd = cdss[0]
    assert fwd.start == 1
    assert fwd.end == 30
    assert fwd.strand == 1
    assert fwd.strand_symbol == '+'
    assert fwd.length == 30
    assert fwd.genomic_span == 30
    assert fwd.is_compound is False
    assert fwd.locus_tag == "TEST_001"
    assert fwd.gene == "testA"

    # Reverse simple CDS
    rev = cdss[1]
    assert rev.start == 51
    assert rev.end == 80
    assert rev.strand == -1
    assert rev.strand_symbol == '-'
    assert rev.length == 30
    assert rev.genomic_span == 30
    assert rev.is_compound is False
    assert rev.locus_tag == "TEST_002"


def test_compound_feature_length_vs_span(compound_joined_gbff: Path) -> None:
    doc = read_genbank(compound_joined_gbff)
    rec = doc.records[0]
    cdss = rec.cds_features
    assert len(cdss) == 3

    ordinary = cdss[0]
    assert ordinary.is_compound is False
    assert ordinary.length == 9
    assert ordinary.genomic_span == 9

    fwd_join = cdss[1]
    assert fwd_join.is_compound is True
    assert fwd_join.length == 30
    assert fwd_join.genomic_span == 39
    assert fwd_join.length < fwd_join.genomic_span
    assert len(fwd_join.join_segments) == 2
    assert fwd_join.join_segments == [(10, 24), (34, 48)]

    rev_join = cdss[2]
    assert rev_join.is_compound is True
    assert rev_join.strand == -1
    assert rev_join.strand_symbol == '-'
    assert rev_join.length == 30
    assert rev_join.genomic_span == 39


def test_feature_dict_compatibility(simple_cds_gbff: Path) -> None:
    doc = read_genbank(simple_cds_gbff)
    feat = doc.records[0].cds_features[0]

    assert feat['type'] == 'CDS'
    assert feat['start'] == 1
    assert feat['end'] == 30
    assert feat['strand'] == '+'
    assert feat['contig'] == 'CONTIG_1'
    assert feat['feature_index'] == 3
    assert 'qualifiers' in feat
    assert feat['join_segments'] == []

    d = feat.to_dict()
    assert isinstance(d, dict)
    assert d['start'] == 1
    assert d['end'] == 30


def test_document_queries(simple_cds_gbff: Path, multi_record_circular_gbff: Path) -> None:
    doc = read_genbank(simple_cds_gbff)
    match = doc.find_locus("TEST_001")
    assert match is not None
    r, f = match
    assert r.id == "CONTIG_1"
    assert f.gene == "testA"

    doc_multi = read_genbank(multi_record_circular_gbff)
    assert len(doc_multi.records) == 2
    assert doc_multi.records[0].id == "CHR_CIRCULAR"
    assert doc_multi.records[0].topology == "circular"
    assert doc_multi.records[1].id == "PLASMID_1"
    assert doc_multi.records[1].topology == "circular"
