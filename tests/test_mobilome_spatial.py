"""One-based, segment-aware linear and circular distance behavior."""

import pytest

from genbank_parser.spatial import (
    circular_feature_distance_bp,
    intervening_gap_bp,
    minimum_feature_covering_arc,
)


def test_intervening_gap_uses_one_based_inclusive_coordinates() -> None:
    assert intervening_gap_bp(((1, 10),), ((11, 20),)) == 0
    assert intervening_gap_bp(((1, 10),), ((20, 30),)) == 9
    assert intervening_gap_bp(((1, 10),), ((5, 15),)) == 0
    assert intervening_gap_bp(((1, 10), (90, 100)), ((20, 30),)) == 9


def test_circular_distance_handles_origin_spanning_segments() -> None:
    assert (
        circular_feature_distance_bp(((900, 950),), ((20, 50),), record_length=1000)
        == 69
    )
    assert (
        circular_feature_distance_bp(((990, 1000),), ((1, 5),), record_length=1000) == 0
    )
    assert (
        circular_feature_distance_bp(
            ((900, 950), (20, 50)), ((60, 70),), record_length=1000
        )
        == 9
    )


def test_spatial_helpers_reject_unlocatable_or_out_of_range_segments() -> None:
    with pytest.raises(ValueError):
        intervening_gap_bp((), ((1, 2),))
    with pytest.raises(ValueError):
        circular_feature_distance_bp(((1, 10),), ((20, 30),), record_length=0)
    with pytest.raises(ValueError):
        circular_feature_distance_bp(((1, 10),), ((20, 101),), record_length=100)


def test_minimum_covering_arc_is_globally_coherent_not_pairwise() -> None:
    arc = minimum_feature_covering_arc(
        ((100, 110), (400, 410), (700, 710)), record_length=1000
    )

    assert arc.span_bp == 611
    assert (arc.start, arc.end, arc.wraps_origin) == (100, 710, False)


def test_minimum_covering_arc_handles_origin_and_inclusive_adjacency() -> None:
    origin = minimum_feature_covering_arc(
        ((990, 1000), (1, 10), (20, 30)), record_length=1000
    )
    compound = minimum_feature_covering_arc(((990, 1000), (1, 10)), record_length=1000)
    adjacent = minimum_feature_covering_arc(((100, 100), (101, 101)), record_length=300)

    assert origin.span_bp == 41
    assert (origin.start, origin.end, origin.wraps_origin) == (990, 30, True)
    assert compound.span_bp == 21
    assert (compound.start, compound.end, compound.wraps_origin) == (990, 10, True)
    assert adjacent.span_bp == 2
    assert (adjacent.start, adjacent.end, adjacent.wraps_origin) == (100, 101, False)


def test_minimum_covering_arc_uses_deterministic_tie_breaking() -> None:
    arc = minimum_feature_covering_arc(
        ((100, 100), (200, 200), (300, 300)), record_length=300
    )

    assert (arc.span_bp, arc.start, arc.end) == (201, 100, 300)
