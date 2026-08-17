"""Core genomic-neighborhood selection regressions."""

from __future__ import annotations

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation

from genbank_parser.model import GenBankFeature, GenBankRecord
from genbank_parser.spatial import select_cds_window


def _feature(index: int, start: int, end: int, *, strand: int = 1) -> GenBankFeature:
    return GenBankFeature(
        record_id="CIRC",
        record_index=0,
        feature_index=index,
        type="CDS",
        location=FeatureLocation(start - 1, end, strand=strand),
        qualifiers={"locus_tag": [f"CDS_{index}"]},
        record_length=1000,
        topology="circular",
    )


def _record(*, circular: bool = True) -> GenBankRecord:
    features = [
        _feature(1, 20, 80),
        _feature(2, 200, 260),
        _feature(3, 500, 560),
        _feature(4, 900, 960),
    ]
    return GenBankRecord(
        id="CIRC",
        name="CIRC",
        description="",
        seq=Seq("A" * 1000),
        length=1000,
        topology="circular" if circular else "linear",
        features=features,
    )


def test_circular_window_crosses_origin_once() -> None:
    record = _record()
    result = select_cds_window(record, record.features[0], 1)
    assert [feature.locus_tag for feature in result.features] == [
        "CDS_4",
        "CDS_1",
        "CDS_2",
    ]
    assert result.wraps_origin
    assert sum(feature.feature_index == 1 for feature in result.features) == 1


def test_oversized_circular_window_has_no_duplicates() -> None:
    record = _record()
    result = select_cds_window(record, record.features[1], 50)
    indices = [feature.feature_index for feature in result.features]
    assert len(indices) == len(set(indices)) == 4
    assert indices == [1, 2, 3, 4]


def test_linear_boundary_and_zero_windows() -> None:
    record = _record(circular=False)
    assert [f.feature_index for f in select_cds_window(record, record.features[0], 2).features] == [1, 2, 3]
    assert [f.feature_index for f in select_cds_window(record, record.features[-1], 2).features] == [2, 3, 4]
    assert [f.feature_index for f in select_cds_window(record, record.features[2], 0).features] == [3]
