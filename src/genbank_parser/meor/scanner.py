"""MEOR feature scanning over the canonical typed parser model."""
from __future__ import annotations

from collections.abc import Iterable

from ..model import GenBankFeature
from .matcher import match_feature_to_marker
from .models import MeorDatabase, MeorHit


def scan_meor_features(
    features: Iterable[GenBankFeature],
    database: MeorDatabase,
    *,
    min_weight: int = 1,
) -> list[MeorHit]:
    if min_weight not in (1, 2, 3):
        raise ValueError("min_weight must be 1, 2, or 3")
    categories = database.category_map
    markers_by_category = {
        category.id: [
            marker for marker in database.markers if marker.category_id == category.id
        ]
        for category in database.categories
    }
    hits: list[MeorHit] = []
    for feature in features:
        if feature.type not in {"CDS", "misc_feature"}:
            continue
        for category in database.categories:
            for marker in markers_by_category[category.id]:
                match = match_feature_to_marker(feature, marker)
                if match is None or match.weight < min_weight:
                    continue
                hits.append(
                    MeorHit(
                        feature_index=feature.feature_index,
                        contig=feature.record_id,
                        locus_tag=feature.locus_tag or "NO_LOCUS",
                        gene=feature.gene,
                        product=feature.product,
                        start=feature.start,
                        end=feature.end,
                        strand=feature.strand_symbol,
                        category_id=category.id,
                        category_title=category.title,
                        marker_id=marker.id,
                        marker_name=marker.name,
                        weight=match.weight,
                        evidence_type=match.evidence_type,
                        evidence_value=match.evidence_value,
                    )
                )
                # Legacy v1.1.0 retains only the first marker match per category.
                break
    return hits
