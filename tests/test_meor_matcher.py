"""Confidence and precedence tests for MEOR feature matching."""
from __future__ import annotations

from Bio.SeqFeature import SeqFeature, SimpleLocation

from genbank_parser.io import extract_xrefs
from genbank_parser.meor.database import load_meor_database
from genbank_parser.meor.matcher import match_feature_to_marker
from genbank_parser.model import GenBankFeature


def make_feature(**qualifiers: list[str]) -> GenBankFeature:
    raw = SeqFeature(SimpleLocation(0, 100, strand=1), type="CDS", qualifiers=qualifiers)
    return GenBankFeature(
        record_id="contig",
        record_index=1,
        feature_index=1,
        type="CDS",
        location=raw.location,
        qualifiers=qualifiers,
        record_length=1000,
        raw_feature=raw,
    )


def test_extract_xrefs_note_selection_is_backward_compatible() -> None:
    feature = make_feature(note=["KEGG:K00496", "EC:1.14.15.3"])
    assert extract_xrefs(feature)["kegg_kos"] == ["K00496"]
    assert extract_xrefs(feature)["ec_numbers"] == ["1.14.15.3"]
    assert extract_xrefs(feature, include_notes=False)["kegg_kos"] == []
    assert extract_xrefs(feature, include_notes=False)["ec_numbers"] == []


def test_matcher_preserves_confidence_tiers_and_precedence() -> None:
    marker = load_meor_database().marker_map["alkB"]

    structured_ko = make_feature(
        db_xref=["KEGG:K00496"], product=["alkane 1-monooxygenase"]
    )
    assert match_feature_to_marker(structured_ko, marker).reason == "KEGG:K00496"
    assert match_feature_to_marker(structured_ko, marker).weight == 3

    structured_ec = make_feature(EC_number=["1.14.15.3"])
    assert match_feature_to_marker(structured_ec, marker).evidence_type == "ec"
    assert match_feature_to_marker(structured_ec, marker).weight == 3

    gene = make_feature(gene=["alkB2"])
    assert match_feature_to_marker(gene, marker).evidence_type == "gene"
    assert match_feature_to_marker(gene, marker).weight == 3

    product = make_feature(product=["putative alkane hydroxylase"])
    assert match_feature_to_marker(product, marker).evidence_type == "product"
    assert match_feature_to_marker(product, marker).weight == 2

    note_ko = make_feature(note=["KEGG:K00496"])
    assert match_feature_to_marker(note_ko, marker).evidence_type == "note_kegg"
    assert match_feature_to_marker(note_ko, marker).weight == 1

    note_ec = make_feature(note=["EC:1.14.15.3"])
    assert match_feature_to_marker(note_ec, marker).evidence_type == "note_ec"
    assert match_feature_to_marker(note_ec, marker).weight == 1

    note_product = make_feature(note=["possible alkane monooxygenase"])
    assert match_feature_to_marker(note_product, marker).evidence_type == "note"
    assert match_feature_to_marker(note_product, marker).weight == 1
