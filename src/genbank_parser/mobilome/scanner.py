"""Field-aware evidence scanning over canonical typed GenBank features."""

from __future__ import annotations

import re
from collections.abc import Iterable
from itertools import pairwise

from ..model import GenBankFeature
from .models import (
    EvidenceReason,
    FeatureRef,
    MobilomeDatabase,
    MobilomeHit,
    MobilomeMarker,
    MobilomeParameterError,
)


def validate_min_evidence(min_evidence: int) -> None:
    """Validate the closed evidence-threshold set shared by API entry points."""

    if isinstance(min_evidence, bool) or min_evidence not in (1, 2, 3):
        raise MobilomeParameterError("min_evidence must be 1, 2, or 3")


def _optional_text(value: str) -> str | None:
    text = value.strip()
    return text or None


def _location_segments(feature: GenBankFeature) -> tuple[tuple[int, int], ...]:
    """Return real one-based inclusive segments without model fallback values."""

    location = feature.location
    try:
        parts = tuple(getattr(location, "parts", ()) or (location,))
    except (AttributeError, TypeError):
        return ()
    segments: list[tuple[int, int]] = []
    try:
        for part in parts:
            start = int(part.start) + 1
            end = int(part.end)
            if start < 1 or end < start:
                return ()
            segments.append((start, end))
    except (AttributeError, TypeError, ValueError):
        return ()
    return tuple(segments)


def _wraps_origin(
    segments: tuple[tuple[int, int], ...],
    *,
    strand: int | None,
    topology: str | None,
    record_length: int,
) -> bool:
    """Detect an origin-spanning ordered compound location on a circle.

    GenBank compound parts are ordered in biological extraction order.  Forward
    parts that step from high to low coordinates, and reverse parts that step
    from low to high coordinates, therefore cross coordinate one.  A simple
    min/max span is deliberately not used here.
    """

    if topology is None or topology.casefold() != "circular" or len(segments) < 2:
        return False
    if record_length <= 0 or any(end > record_length for _, end in segments):
        return False
    if strand == -1:
        anchors = [end for _, end in segments]
        return any(current > previous for previous, current in pairwise(anchors))
    anchors = [start for start, _ in segments]
    return any(current < previous for previous, current in pairwise(anchors))


def feature_ref(feature: GenBankFeature) -> FeatureRef:
    """Build a public reference solely from canonical model semantics."""

    segments = _location_segments(feature)
    location = feature.location
    if segments:
        start = (
            segments[0][0] if len(segments) == 1 else min(item[0] for item in segments)
        )
        end = (
            segments[0][1] if len(segments) == 1 else max(item[1] for item in segments)
        )
        try:
            biological_length = len(location)
        except (AttributeError, TypeError, ValueError):
            biological_length = None
    else:
        start = None
        end = None
        biological_length = None
    try:
        strand = location.strand
    except AttributeError:
        strand = None
    if strand not in (-1, 0, 1, None):
        strand = None
    operator = getattr(location, "operator", None)
    if not isinstance(operator, str):
        operator = None
    return FeatureRef(
        record_index=feature.record_index,
        record_id=feature.record_id,
        feature_index=feature.feature_index,
        feature_type=feature.type,
        locus_tag=_optional_text(feature.locus_tag),
        gene=_optional_text(feature.gene),
        product=_optional_text(feature.product),
        start=start,
        end=end,
        strand=strand,
        segments=segments,
        location_operator=operator,
        biological_length=biological_length,
        wraps_origin=_wraps_origin(
            segments,
            strand=strand,
            topology=feature.topology,
            record_length=feature.record_length,
        ),
        is_partial=feature.is_partial,
        is_pseudo=feature.is_pseudo,
    )


def values_for_declared_field(feature: GenBankFeature, field: str) -> tuple[str, ...]:
    """Return values for one catalog-declared field without cross-field joining."""

    if field == "feature_type":
        return (feature.type,)
    return tuple(str(value) for value in feature.get_quals(field))


def _match_value(
    value: str,
    *,
    mode: str,
    pattern: str,
    negative_patterns: tuple[str, ...],
) -> str | None:
    trimmed = value.strip()
    if not trimmed:
        return None
    if any(re.search(negative, trimmed) for negative in negative_patterns):
        return None
    if mode == "exact":
        return trimmed if trimmed == pattern else None
    if mode == "casefold_exact":
        return trimmed if trimmed.casefold() == pattern.strip().casefold() else None
    match = re.search(pattern, trimmed)
    return match.group(0) if match is not None else None


def _reason_key(
    marker: MobilomeMarker,
    reference: FeatureRef,
    reason: EvidenceReason,
) -> tuple[object, ...]:
    return (
        marker.id,
        reference.record_index,
        reference.feature_index,
        reason.field,
        reason.match_mode,
        reason.matched_value,
        reason.matched_text,
        reason.pattern,
        reason.strength,
        reason.source_ids,
    )


def scan_mobilome_features(
    features: Iterable[GenBankFeature],
    database: MobilomeDatabase,
    *,
    min_evidence: int = 1,
) -> tuple[MobilomeHit, ...]:
    """Retain every declared match, including pseudo and below-threshold hits."""

    validate_min_evidence(min_evidence)
    hits: list[MobilomeHit] = []
    for feature in sorted(
        features, key=lambda item: (item.record_index, item.feature_index)
    ):
        reference = feature_ref(feature)
        for marker in database.markers:
            if feature.type.casefold() not in marker.feature_types_casefold:
                continue
            seen: set[tuple[object, ...]] = set()
            reasons: list[tuple[int, int, EvidenceReason]] = []
            for matcher_index, matcher in enumerate(marker.matchers):
                for value in values_for_declared_field(feature, matcher.field):
                    for pattern_index, pattern in enumerate(matcher.patterns):
                        matched_text = _match_value(
                            value,
                            mode=matcher.mode,
                            pattern=pattern,
                            negative_patterns=matcher.negative_patterns,
                        )
                        if matched_text is None:
                            continue
                        eligible = (
                            matcher.strength >= min_evidence
                            and (
                                not marker.requires_non_pseudo
                                or not reference.is_pseudo
                            )
                            and (
                                not marker.requires_complete or not reference.is_partial
                            )
                        )
                        reason = EvidenceReason(
                            reason_id=(
                                f"r{reference.record_index}:f{reference.feature_index}:"
                                f"m{marker.id}:{matcher.field}:{matcher_index}:{pattern_index}"
                            ),
                            field=matcher.field,
                            match_mode=matcher.mode,
                            matched_value=value,
                            matched_text=matched_text,
                            pattern=pattern,
                            strength=matcher.strength,
                            eligible=eligible,
                            source_ids=marker.sources,
                        )
                        key = _reason_key(marker, reference, reason)
                        if key in seen:
                            continue
                        seen.add(key)
                        reasons.append((matcher_index, pattern_index, reason))
            if not reasons:
                continue
            reasons.sort(
                key=lambda item: (
                    item[2].field,
                    item[0],
                    item[1],
                    item[2].matched_value,
                )
            )
            ordered_reasons = tuple(item[2] for item in reasons)
            hits.append(
                MobilomeHit(
                    hit_id=f"r{reference.record_index}:f{reference.feature_index}:m{marker.id}",
                    marker_id=marker.id,
                    marker_label=marker.label,
                    facets=marker.facets,
                    feature=reference,
                    reasons=ordered_reasons,
                    max_strength=max(reason.strength for reason in ordered_reasons),
                    eligible_for_inference=any(
                        reason.eligible for reason in ordered_reasons
                    ),
                )
            )
    return tuple(
        sorted(
            hits,
            key=lambda item: (
                item.feature.record_index,
                item.feature.feature_index,
                item.marker_id,
            ),
        )
    )


__all__ = [
    "feature_ref",
    "scan_mobilome_features",
    "validate_min_evidence",
    "values_for_declared_field",
]
