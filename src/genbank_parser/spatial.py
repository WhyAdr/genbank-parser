"""Shared target resolution and circular genomic-window helpers."""

from __future__ import annotations

from dataclasses import dataclass

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
    wraps_origin = any(
        current < previous for previous, current in zip(indices, indices[1:])
    )
    return FeatureWindow(tuple(cdss[index] for index in indices), wraps_origin)
