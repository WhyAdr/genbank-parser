"""One-based, segment-aware linear and circular distance behavior."""

import pytest

from genbank_parser.spatial import circular_feature_distance_bp, intervening_gap_bp


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
