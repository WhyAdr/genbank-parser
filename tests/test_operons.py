"""Operon-pair spatial semantics."""

from __future__ import annotations

from Bio.SeqFeature import FeatureLocation

from genbank_parser.model import GenBankFeature
from genbank_parser.operons import build_operon_result, find_operon_pairs


def _cds(index: int, start: int, end: int, strand: int = 1) -> GenBankFeature:
    return GenBankFeature(
        record_id="R",
        record_index=0,
        feature_index=index,
        type="CDS",
        location=FeatureLocation(start - 1, end, strand=strand),
        qualifiers={"locus_tag": [f"CDS_{index}"]},
    )


def test_operon_gap_boundaries_and_strand_filter() -> None:
    features = [
        _cds(1, 1, 100),
        _cds(2, 51, 150),
        _cds(3, 301, 350),
        _cds(4, 401, 450, strand=-1),
    ]
    pairs = find_operon_pairs(features, min_gap=-50, max_gap=150)
    assert [(a.feature_index, b.feature_index, gap) for a, b, gap in pairs] == [
        (1, 2, -50),
        (2, 3, 150),
    ]


def test_circular_operon_pair_connects_last_to_first() -> None:
    features = [_cds(1, 20, 80), _cds(2, 900, 960)]
    pairs = find_operon_pairs(
        features,
        min_gap=0,
        max_gap=59,
        circular=True,
        record_length=1000,
    )
    assert [(a.feature_index, b.feature_index, gap) for a, b, gap in pairs] == [
        (2, 1, 59)
    ]


def test_single_cds_is_not_paired_with_itself() -> None:
    feature = _cds(1, 20, 80)
    assert not find_operon_pairs(
        [feature], circular=True, record_length=1000
    )


def test_structured_three_gene_cluster_and_break() -> None:
    features = [
        _cds(1, 1, 100),
        _cds(2, 120, 200),
        _cds(3, 220, 300),
        _cds(4, 600, 700),
    ]
    result = build_operon_result(features, max_gap=50)
    assert len(result.pairs) == 2
    assert len(result.clusters) == 1
    assert [feature.feature_index for feature in result.clusters[0].features] == [
        1,
        2,
        3,
    ]
    assert result.clusters[0].gaps == (19, 19)


def test_structured_circular_cluster_orders_across_origin() -> None:
    features = [
        _cds(1, 20, 80),
        _cds(2, 100, 160),
        _cds(3, 900, 960),
    ]
    result = build_operon_result(
        features,
        min_gap=0,
        max_gap=59,
        circular=True,
        record_length=1000,
    )
    assert len(result.clusters) == 1
    assert [feature.feature_index for feature in result.clusters[0].features] == [
        3,
        1,
        2,
    ]
    assert result.clusters[0].gaps == (59, 19)
