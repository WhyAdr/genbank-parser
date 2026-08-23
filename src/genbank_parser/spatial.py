"""Shared target resolution and circular genomic-window helpers."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from itertools import pairwise

from .model import GenBankDocument, GenBankFeature, GenBankRecord


class TargetNotFoundError(LookupError):
    """Raised when neither a locus tag nor a gene name resolves a target."""


class NonCDSTargetError(ValueError):
    """Raised when a target resolves only to a non-CDS feature."""


@dataclass(frozen=True)
class FeatureWindow:
    """A deterministic CDS window in display order."""

    features: tuple[GenBankFeature, ...]
    wraps_origin: bool


def _validated_segments(
    segments: Sequence[tuple[int, int]],
    *,
    owner: str,
) -> tuple[tuple[int, int], ...]:
    if not segments:
        raise ValueError(f"{owner} must contain at least one located segment")
    normalized: list[tuple[int, int]] = []
    for start, end in segments:
        if isinstance(start, bool) or isinstance(end, bool) or start < 1 or end < start:
            raise ValueError(f"{owner} has an invalid one-based inclusive segment")
        normalized.append((start, end))
    return tuple(normalized)


def intervening_gap_bp(
    left_segments: Sequence[tuple[int, int]],
    right_segments: Sequence[tuple[int, int]],
) -> int:
    """Return the minimum non-negative linear gap strictly between segments.

    Coordinates are one-based and inclusive.  Adjacent or overlapping segments
    have a gap of zero.  Compound locations retain their individual segments,
    so a broad min/max envelope cannot inflate a gap through an internal join.
    """

    left = _validated_segments(left_segments, owner="left_segments")
    right = _validated_segments(right_segments, owner="right_segments")
    gaps: list[int] = []
    for left_start, left_end in left:
        for right_start, right_end in right:
            if left_end < right_start:
                gaps.append(right_start - left_end - 1)
            elif right_end < left_start:
                gaps.append(left_start - right_end - 1)
            else:
                gaps.append(0)
    return min(gaps)


def circular_feature_distance_bp(
    first_segments: Sequence[tuple[int, int]],
    second_segments: Sequence[tuple[int, int]],
    *,
    record_length: int,
) -> int:
    """Return the minimum intervening gap in either direction on a circle."""

    if isinstance(record_length, bool) or record_length <= 0:
        raise ValueError("record_length must be a positive integer")
    first = _validated_segments(first_segments, owner="first_segments")
    second = _validated_segments(second_segments, owner="second_segments")
    if any(end > record_length for _, end in (*first, *second)):
        raise ValueError("segments must lie within record_length")
    if intervening_gap_bp(first, second) == 0:
        return 0
    gaps: list[int] = []
    for first_start, first_end in first:
        for second_start, second_end in second:
            clockwise = (second_start - first_end - 1) % record_length
            anticlockwise = (first_start - second_end - 1) % record_length
            gaps.extend((clockwise, anticlockwise))
    return min(gaps)


def feature_display_bounds(
    feature: GenBankFeature,
    *,
    circular: bool,
    record_length: int,
) -> tuple[int, int]:
    """Return one unwrapped span, including compound features crossing origin."""
    if circular and feature.is_compound:
        segments = feature.join_segments
        touches_start = any(start == 1 for start, _ in segments)
        touches_end = any(end == record_length for _, end in segments)
        if touches_start and touches_end:
            high_start = max(start for start, _ in segments)
            low_end = min(end for _, end in segments)
            return high_start, low_end + record_length
    return feature.start, feature.end


def resolve_target(
    document: GenBankDocument,
    target: str,
) -> tuple[GenBankRecord, GenBankFeature]:
    """Resolve an exact locus tag, then an exact gene name, preferring CDSs."""
    locus_match = document.find_locus(target)
    if locus_match is not None:
        record, feature = locus_match
        if feature.type.casefold() != "cds":
            raise NonCDSTargetError(
                f"Resolved target {target!r} is a {feature.type!r}, not a CDS"
            )
        return record, feature

    folded = target.casefold()
    non_cds_match: tuple[GenBankRecord, GenBankFeature] | None = None
    for record in document.records:
        for feature in record.features:
            if not feature.gene or feature.gene.casefold() != folded:
                continue
            if feature.type.casefold() == "cds":
                return record, feature
            if non_cds_match is None:
                non_cds_match = (record, feature)

    if non_cds_match is not None:
        _, feature = non_cds_match
        raise NonCDSTargetError(
            f"Resolved target {target!r} is a {feature.type!r}, not a CDS"
        )
    raise TargetNotFoundError(f"Locus tag or gene {target!r} was not found")


def select_cds_window(
    record: GenBankRecord,
    target: GenBankFeature,
    window: int,
) -> FeatureWindow:
    """Select at most one copy of each physical CDS around ``target``."""
    if window < 0:
        raise ValueError("window must be non-negative")

    cdss = sorted(
        record.cds_features,
        key=lambda feature: (feature.start, feature.end, feature.feature_index),
    )
    if not cdss:
        raise ValueError(f"Record {record.id!r} has no CDS features")
    try:
        target_index = next(
            index
            for index, feature in enumerate(cdss)
            if feature.feature_index == target.feature_index
        )
    except StopIteration:
        raise NonCDSTargetError(
            f"Resolved target {target.locus_tag or target.gene!r} is not in CDS ordering"
        ) from None

    if record.topology != "circular":
        start = max(0, target_index - window)
        end = min(len(cdss), target_index + window + 1)
        return FeatureWindow(tuple(cdss[start:end]), False)

    count = min(len(cdss), 2 * window + 1)
    left_count = min(window, (count - 1) // 2)
    right_count = count - left_count - 1
    indices = [
        (target_index + offset) % len(cdss)
        for offset in range(-left_count, right_count + 1)
    ]
    wraps_origin = any(current < previous for previous, current in pairwise(indices))
    return FeatureWindow(tuple(cdss[index] for index in indices), wraps_origin)
