"""Conservative, configured mobilome inference over retained evidence hits."""

from __future__ import annotations

from collections.abc import Sequence

from ..spatial import circular_feature_distance_bp, intervening_gap_bp
from .models import (
    HypothesisParticipant,
    InferenceComponent,
    InferenceRule,
    MissingComponent,
    MobilomeDatabase,
    MobilomeHit,
    MobilomeHypothesis,
    RepliconAssessment,
    RepliconInventory,
)


def _functional_hits(hits: Sequence[MobilomeHit]) -> tuple[MobilomeHit, ...]:
    """Retain eligible, non-pseudo, complete candidates for component rules."""

    return tuple(
        hit
        for hit in hits
        if hit.eligible_for_inference
        and not hit.feature.is_pseudo
        and not hit.feature.is_partial
    )


def _feature_key(hit: MobilomeHit) -> tuple[int, int]:
    return hit.feature.record_index, hit.feature.feature_index


def _facet_hits(hits: Sequence[MobilomeHit], facet: str) -> tuple[MobilomeHit, ...]:
    return tuple(hit for hit in _functional_hits(hits) if facet in hit.facets)


def _choose_hits(
    candidates: Sequence[MobilomeHit],
    count: int,
    *,
    forbidden_features: set[tuple[int, int]],
) -> tuple[MobilomeHit, ...] | None:
    """Choose deterministic distinct-feature support for one required facet."""

    selected: list[MobilomeHit] = []
    seen_features = set(forbidden_features)
    for hit in candidates:
        key = _feature_key(hit)
        if key in seen_features:
            continue
        selected.append(hit)
        seen_features.add(key)
        if len(selected) == count:
            return tuple(selected)
    return None


def _component_support(
    hits: Sequence[MobilomeHit],
    component: InferenceComponent,
    *,
    forbidden_features: set[tuple[int, int]] | None = None,
    distinct_features: bool,
) -> tuple[tuple[MobilomeHit, ...] | None, tuple[str, ...]]:
    """Evaluate one configured component and return selected hits plus misses."""

    forbidden = set() if forbidden_features is None else set(forbidden_features)
    minimums = component.minimums
    candidates_by_facet = {
        facet: _facet_hits(hits, facet) for facet in component.facets
    }
    missing = tuple(
        facet
        for facet in component.facets
        if len({_feature_key(hit) for hit in candidates_by_facet[facet]})
        < minimums.get(facet, 1)
    )
    if component.operator == "any":
        for facet in component.facets:
            chosen_any = _choose_hits(
                candidates_by_facet[facet],
                minimums.get(facet, 1),
                forbidden_features=forbidden,
            )
            if chosen_any is not None:
                return chosen_any, ()
        return None, missing or component.facets
    if missing:
        return None, missing
    missing_facets: list[str] = []
    selected_hits: list[MobilomeHit] = []
    used = set(forbidden)
    for facet in component.facets:
        chosen = _choose_hits(
            candidates_by_facet[facet],
            minimums.get(facet, 1),
            forbidden_features=used if distinct_features else set(),
        )
        if chosen is None:
            missing_facets.append(facet)
            continue
        selected_hits.extend(chosen)
        if distinct_features:
            used.update(_feature_key(hit) for hit in chosen)
    if missing_facets:
        return None, tuple(missing_facets)
    return tuple(selected_hits), ()


def _rule_component_support(
    hits: Sequence[MobilomeHit],
    rule: InferenceRule,
    database: MobilomeDatabase,
    *,
    role: str,
) -> tuple[tuple[MobilomeHit, ...] | None, tuple[MissingComponent, ...]]:
    """Evaluate every component required for one role in a rule."""

    component_ids = rule.components_by_role.get(role, ())
    selected: list[MobilomeHit] = []
    used: set[tuple[int, int]] = set()
    missing: list[MissingComponent] = []
    for component_id in component_ids:
        component = database.component_map[component_id]
        component_hits, component_missing = _component_support(
            hits,
            component,
            forbidden_features=used,
            distinct_features=rule.distinct_features,
        )
        if component_hits is None:
            missing.append(
                MissingComponent(
                    component_id=component_id,
                    missing_facets=component_missing,
                )
            )
            continue
        selected.extend(component_hits)
        if rule.distinct_features:
            used.update(_feature_key(hit) for hit in component_hits)
    if missing:
        return None, tuple(missing)
    return tuple(selected), ()


def _participant(
    inventory: RepliconInventory, role: str = "subject"
) -> HypothesisParticipant:
    return HypothesisParticipant(
        role=role,  # type: ignore[arg-type]
        record_index=inventory.record_index,
        record_id=inventory.record_id,
    )


def _hypothesis(
    *,
    hypothesis_id: str,
    rule: InferenceRule,
    inventory: RepliconInventory,
    supporting: Sequence[MobilomeHit],
    missing: Sequence[MissingComponent] = (),
    summary: str | None = None,
    limitations: Sequence[str] | None = None,
) -> MobilomeHypothesis:
    return MobilomeHypothesis(
        hypothesis_id=hypothesis_id,
        rule_id=rule.id,
        kind=rule.kind,
        status=rule.status if not missing else "insufficient",
        summary=summary or rule.wording,
        participants=(_participant(inventory),),
        supporting_hit_ids=tuple(sorted({hit.hit_id for hit in supporting})),
        missing_components=tuple(missing),
        conflicting_hit_ids=(),
        limitations=tuple(limitations if limitations is not None else rule.limitations),
        source_ids=tuple(sorted(rule.sources)),
    )


def _raw_marker_hits(
    hits: Sequence[MobilomeHit], component: InferenceComponent
) -> tuple[MobilomeHit, ...]:
    return tuple(
        hit for hit in hits if any(facet in hit.facets for facet in component.facets)
    )


def _infer_component_rule(
    inventory: RepliconInventory,
    hits: Sequence[MobilomeHit],
    database: MobilomeDatabase,
    rule_id: str,
    role: str = "subject_record",
) -> MobilomeHypothesis | None:
    rule = database.inference_rule_map.get(rule_id)
    if rule is None:
        return None
    selected, missing = _rule_component_support(hits, rule, database, role=role)
    if selected is not None:
        return _hypothesis(
            hypothesis_id=f"h:{rule.id}:{inventory.record_index}",
            rule=rule,
            inventory=inventory,
            supporting=selected,
        )
    component_ids = rule.components_by_role.get(role, ())
    raw = tuple(
        hit
        for component_id in component_ids
        for hit in _raw_marker_hits(hits, database.component_map[component_id])
    )
    if not raw and rule_id != "replication_annotation_candidate":
        return None
    extra_limitations: list[str] = list(rule.limitations)
    if raw and all(hit.feature.is_pseudo for hit in raw):
        extra_limitations.append(
            "Observed candidates are pseudogenized; autonomy and function are unresolved."
        )
    elif raw:
        extra_limitations.append(
            "Observed candidates do not meet the configured functional component policy."
        )
    else:
        extra_limitations.append(
            "No annotated replication candidate was observed; autonomy is unresolved."
        )
    default_missing = tuple(
        MissingComponent(
            component_id=cid,
            missing_facets=database.component_map[cid].facets,
        )
        for cid in component_ids
    )
    return _hypothesis(
        hypothesis_id=f"h:{rule.id}:{inventory.record_index}",
        rule=rule,
        inventory=inventory,
        supporting=raw,
        missing=missing or default_missing,
        summary="Replication evidence insufficient; mechanism unresolved"
        if rule_id == "replication_annotation_candidate"
        else f"{rule.wording} component pattern is incomplete",
        limitations=tuple(extra_limitations),
    )


def _distance_between_hits(
    first: MobilomeHit,
    second: MobilomeHit,
    inventory: RepliconInventory,
) -> int | None:
    if not first.feature.segments or not second.feature.segments:
        return None
    try:
        if inventory.topology == "circular":
            return circular_feature_distance_bp(
                first.feature.segments,
                second.feature.segments,
                record_length=inventory.length,
            )
        return intervening_gap_bp(first.feature.segments, second.feature.segments)
    except ValueError:
        return None


def _infer_toxin_antitoxin_pairs(
    inventory: RepliconInventory,
    hits: Sequence[MobilomeHit],
    database: MobilomeDatabase,
) -> tuple[MobilomeHypothesis, ...]:
    hypotheses: list[MobilomeHypothesis] = []
    for rule in database.inference_rules:
        if rule.kind != "toxin_antitoxin" or len(rule.required_marker_ids) != 2:
            continue
        first_marker, second_marker = rule.required_marker_ids
        first_hits = [
            hit for hit in _functional_hits(hits) if hit.marker_id == first_marker
        ]
        second_hits = [
            hit for hit in _functional_hits(hits) if hit.marker_id == second_marker
        ]
        for first in first_hits:
            for second in second_hits:
                if rule.distinct_features and _feature_key(first) == _feature_key(
                    second
                ):
                    continue
                gap = _distance_between_hits(first, second, inventory)
                if (
                    gap is None
                    or rule.max_circular_gap_bp is None
                    or gap > rule.max_circular_gap_bp
                ):
                    continue
                hypotheses.append(
                    _hypothesis(
                        hypothesis_id=(
                            f"h:{rule.id}:{inventory.record_index}:"
                            f"{min(first.feature.feature_index, second.feature.feature_index)}:"
                            f"{max(first.feature.feature_index, second.feature.feature_index)}"
                        ),
                        rule=rule,
                        inventory=inventory,
                        supporting=(first, second),
                        limitations=tuple(rule.limitations)
                        + (
                            (
                                f"Configured maximum gap ({inventory.topology}): "
                                f"{rule.max_circular_gap_bp} bp; observed gap: {gap} bp."
                            ),
                        ),
                    )
                )
    return tuple(hypotheses)


def record_inference_limitations(
    inventory: RepliconInventory,
    hits: Sequence[MobilomeHit],
    database: MobilomeDatabase,
) -> tuple[str, ...]:
    """Return record-level limits that cannot be represented by a pair hit."""

    limitations: list[str] = []
    local_hits = tuple(
        hit for hit in hits if hit.feature.record_index == inventory.record_index
    )
    for rule in database.inference_rules:
        if rule.kind != "toxin_antitoxin" or len(rule.required_marker_ids) != 2:
            continue
        first_marker, second_marker = rule.required_marker_ids
        first_hits = tuple(
            hit for hit in _functional_hits(local_hits) if hit.marker_id == first_marker
        )
        second_hits = tuple(
            hit
            for hit in _functional_hits(local_hits)
            if hit.marker_id == second_marker
        )
        if any(
            not first.feature.segments or not second.feature.segments
            for first in first_hits
            for second in second_hits
        ):
            limitations.append(
                f"{rule.id}: at least one same-record toxin/antitoxin candidate is unlocatable; spatial pairing was skipped."
            )
    return tuple(sorted(set(limitations)))


def infer_replicon_hypotheses(
    replicon: RepliconInventory,
    hits: Sequence[MobilomeHit],
    database: MobilomeDatabase,
) -> tuple[MobilomeHypothesis, ...]:
    """Infer cautious record-local assessments from one record's retained hits."""

    local_hits = tuple(
        hit for hit in hits if hit.feature.record_index == replicon.record_index
    )
    hypotheses: list[MobilomeHypothesis] = []
    for rule_id in (
        "replication_annotation_candidate",
        "mobilizable_core_annotation_candidate",
        "helper_machinery_annotation_candidate",
    ):
        hypothesis = _infer_component_rule(replicon, local_hits, database, rule_id)
        if hypothesis is not None:
            hypotheses.append(hypothesis)
    hypotheses.extend(_infer_toxin_antitoxin_pairs(replicon, local_hits, database))
    return tuple(
        sorted(
            hypotheses,
            key=lambda item: (
                item.participants[0].record_index,
                item.rule_id,
                item.hypothesis_id,
            ),
        )
    )


def infer_cross_record_hypotheses(
    assessments: Sequence[RepliconAssessment],
    database: MobilomeDatabase,
) -> tuple[MobilomeHypothesis, ...]:
    """Evaluate the single v1 helper-dependent mobility rule across records."""

    rule = database.inference_rule_map.get("possible_helper_dependent_mobilization")
    if rule is None:
        return ()
    hypotheses: list[MobilomeHypothesis] = []
    for target in assessments:
        target_support, _ = _rule_component_support(
            target.hits, rule, database, role="target_record"
        )
        if target_support is None:
            continue
        for helper in assessments:
            if helper.inventory.record_index == target.inventory.record_index:
                continue
            helper_support, _ = _rule_component_support(
                helper.hits, rule, database, role="distinct_helper_record"
            )
            if helper_support is None:
                continue
            hypotheses.append(
                MobilomeHypothesis(
                    hypothesis_id=(
                        f"h:{rule.id}:{target.inventory.record_index}:"
                        f"{helper.inventory.record_index}"
                    ),
                    rule_id=rule.id,
                    kind=rule.kind,
                    status=rule.status,
                    summary=rule.wording,
                    participants=(
                        _participant(target.inventory, "target"),
                        _participant(helper.inventory, "helper"),
                    ),
                    supporting_hit_ids=tuple(
                        sorted(
                            {hit.hit_id for hit in (*target_support, *helper_support)}
                        )
                    ),
                    missing_components=(),
                    conflicting_hit_ids=(),
                    limitations=rule.limitations,
                    source_ids=tuple(sorted(rule.sources)),
                )
            )
    return tuple(
        sorted(
            hypotheses,
            key=lambda item: (
                item.participants[0].record_index,
                item.rule_id,
                item.hypothesis_id,
            ),
        )
    )


__all__ = [
    "infer_cross_record_hypotheses",
    "infer_replicon_hypotheses",
    "record_inference_limitations",
]
