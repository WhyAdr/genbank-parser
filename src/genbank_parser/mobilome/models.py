"""Immutable contracts for the annotation-supported mobilome report.

The objects in this module deliberately describe annotations and configured
inference policy.  They do not encode phenotype, transfer, replicon ancestry,
or sequence-homology conclusions.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

SCHEMA_VERSION = "gbparse.mobilome.v1"

RepliconClass = Literal["chromosome", "plasmid", "unknown"]
RepliconTopology = Literal["circular", "linear", "unknown"]
AssessmentStatus = Literal["supported", "tentative", "insufficient", "conflicting"]
MatchMode = Literal["exact", "casefold_exact", "regex"]
HandoffStatus = Literal["not_run"]
ParticipantRole = Literal["subject", "target", "helper"]
DatabaseSource = Literal["packaged", "custom"]


class MobilomeError(Exception):
    """Base class for expected, user-facing mobilome failures."""


class MobilomeInputError(MobilomeError):
    """Raised when an input cannot produce a typed GenBank document."""


class MobilomeDatabaseError(MobilomeError):
    """Raised when a packaged or custom evidence catalog is invalid."""


class MobilomeParameterError(MobilomeError):
    """Raised when a public mobilome API parameter is outside its contract."""


class MobilomeOutputError(MobilomeError):
    """Raised when a report cannot be rendered or written safely."""


@dataclass(frozen=True)
class FeatureRef:
    """Stable, location-aware identity for one canonical parser feature."""

    record_index: int
    record_id: str
    feature_index: int
    feature_type: str
    locus_tag: str | None
    gene: str | None
    product: str | None
    start: int | None
    end: int | None
    strand: int | None
    segments: tuple[tuple[int, int], ...]
    location_operator: str | None
    biological_length: int | None
    wraps_origin: bool
    is_partial: bool
    is_pseudo: bool


@dataclass(frozen=True)
class RecordEvidence:
    """One record-level classification reason with exact annotation provenance."""

    rule_id: str
    source: Literal[
        "source_qualifier",
        "record_annotation",
        "record_description",
        "record_name",
    ]
    field: str
    value: str
    matched_text: str
    pattern: str
    source_feature_index: int | None
    interpreted_as: RepliconClass
    strength: int
    source_ids: tuple[str, ...]


@dataclass(frozen=True)
class SourceQualifier:
    """Normalized qualifiers from a typed ``source`` feature."""

    key: str
    values: tuple[str, ...]


@dataclass(frozen=True)
class CatalogResource:
    """Exact-byte provenance for a resource used by an analysis."""

    name: str
    sha256: str


@dataclass(frozen=True)
class AnalysisParameters:
    include: Literal["all", "chromosome", "plasmid", "unknown"]
    min_evidence: int
    database_source: DatabaseSource


@dataclass(frozen=True)
class AnnotationProvenance:
    record_indices: tuple[int, ...]
    name: str | None
    version: str | None
    database_version: str | None
    evidence_fields: tuple[str, ...]


@dataclass(frozen=True)
class EvidenceReason:
    """One exact field/value/pattern match retained within a marker hit."""

    reason_id: str
    field: str
    match_mode: MatchMode
    matched_value: str
    matched_text: str
    pattern: str
    strength: int
    eligible: bool
    source_ids: tuple[str, ...]


@dataclass(frozen=True)
class MobilomeHit:
    """All evidence reasons for one marker on one feature."""

    hit_id: str
    marker_id: str
    marker_label: str
    facets: tuple[str, ...]
    feature: FeatureRef
    reasons: tuple[EvidenceReason, ...]
    max_strength: int
    eligible_for_inference: bool


@dataclass(frozen=True)
class RepliconClassification:
    label: RepliconClass
    status: AssessmentStatus
    supporting_evidence: tuple[RecordEvidence, ...]
    conflicting_evidence: tuple[RecordEvidence, ...]
    limitations: tuple[str, ...]


@dataclass(frozen=True)
class HypothesisParticipant:
    role: ParticipantRole
    record_index: int
    record_id: str


@dataclass(frozen=True)
class MissingComponent:
    """One unsatisfied component requirement with its missing facets."""

    component_id: str
    missing_facets: tuple[str, ...] = ()


@dataclass(frozen=True)
class MobilomeHypothesis:
    """A configured, cautious hypothesis or an explicit insufficiency result."""

    hypothesis_id: str
    rule_id: str
    kind: str
    status: AssessmentStatus
    summary: str
    participants: tuple[HypothesisParticipant, ...]
    supporting_hit_ids: tuple[str, ...]
    missing_components: tuple[MissingComponent, ...]
    conflicting_hit_ids: tuple[str, ...]
    limitations: tuple[str, ...]
    source_ids: tuple[str, ...]


@dataclass(frozen=True)
class RepliconInventory:
    record_index: int
    record_id: str
    name: str
    description: str
    length: int
    topology: RepliconTopology
    source_qualifiers: tuple[SourceQualifier, ...]
    classification: RepliconClassification
    feature_count: int
    cds_count: int
    pseudo_feature_count: int
    pseudo_cds_count: int
    pseudogene_feature_count: int
    spatial_limitations: tuple[str, ...] = ()


@dataclass(frozen=True)
class RepliconAssessment:
    inventory: RepliconInventory
    hits: tuple[MobilomeHit, ...]
    hypotheses: tuple[MobilomeHypothesis, ...]


@dataclass(frozen=True)
class DisabledAggregateRule:
    rule_id: str
    reason: str
    evidence_retained_as: tuple[str, ...]
    enablement_requirements: tuple[str, ...]
    source_ids: tuple[str, ...]


@dataclass(frozen=True)
class ExternalHandoff:
    handoff_id: str
    tool: str
    purpose: str
    status: HandoffStatus
    reason: str
    required_input: tuple[str, ...]
    required_provenance_fields: tuple[str, ...]
    limitations: tuple[str, ...]
    source_ids: tuple[str, ...] = ()


@dataclass(frozen=True)
class MobilomeReport:
    source_file: str
    source_sha256: str
    tool_version: str
    catalog_version: str
    inference_version: str
    resources: tuple[CatalogResource, ...]
    parameters: AnalysisParameters
    annotation_provenance: tuple[AnnotationProvenance, ...]
    inventory: tuple[RepliconInventory, ...]
    replicons: tuple[RepliconAssessment, ...]
    cross_record_hypotheses: tuple[MobilomeHypothesis, ...]
    disabled_aggregate_rules: tuple[DisabledAggregateRule, ...]
    handoffs: tuple[ExternalHandoff, ...]
    limitations: tuple[str, ...]
    scanned_replicons: tuple[RepliconAssessment, ...] = field(
        default_factory=tuple, repr=False, compare=False
    )


@dataclass(frozen=True)
class Facet:
    id: str
    title: str


@dataclass(frozen=True)
class MarkerMatcher:
    field: str
    mode: MatchMode
    patterns: tuple[str, ...]
    negative_patterns: tuple[str, ...]
    strength: int


@dataclass(frozen=True)
class MobilomeMarker:
    id: str
    label: str
    facets: tuple[str, ...]
    feature_types: tuple[str, ...]
    matchers: tuple[MarkerMatcher, ...]
    requires_non_pseudo: bool
    requires_complete: bool
    sources: tuple[str, ...]
    interpretation: str
    limitations: tuple[str, ...]

    @property
    def feature_types_casefold(self) -> frozenset[str]:
        return frozenset(feature_type.casefold() for feature_type in self.feature_types)


@dataclass(frozen=True)
class ProvenanceSource:
    id: str
    source_type: str
    payload: dict[str, object] = field(repr=False, compare=False)


@dataclass(frozen=True)
class InferenceComponent:
    id: str
    operator: Literal["any", "all"]
    facets: tuple[str, ...]
    min_distinct_features: tuple[tuple[str, int], ...]

    @property
    def minimums(self) -> dict[str, int]:
        return dict(self.min_distinct_features)


@dataclass(frozen=True)
class InferenceRule:
    id: str
    kind: str
    status: AssessmentStatus
    wording: str
    sources: tuple[str, ...]
    limitations: tuple[str, ...]
    required_components: tuple[tuple[str, tuple[str, ...]], ...] = ()
    required_marker_ids: tuple[str, ...] = ()
    same_record: bool = False
    distinct_features: bool = False
    max_circular_gap_bp: int | None = None

    @property
    def components_by_role(self) -> dict[str, tuple[str, ...]]:
        return dict(self.required_components)


@dataclass(frozen=True)
class MobilomeDatabase:
    """Validated evidence and inference resources used for one analysis."""

    facets: tuple[Facet, ...]
    markers: tuple[MobilomeMarker, ...]
    provenance_sources: tuple[ProvenanceSource, ...]
    components: tuple[InferenceComponent, ...]
    inference_rules: tuple[InferenceRule, ...]
    disabled_aggregate_rules: tuple[DisabledAggregateRule, ...]
    handoffs: tuple[ExternalHandoff, ...]
    catalog_version: str
    inference_version: str
    provenance_version: str
    resources: tuple[CatalogResource, ...]
    source_paths: tuple[Path, ...]
    database_source: DatabaseSource

    @property
    def facet_map(self) -> dict[str, Facet]:
        return {facet.id: facet for facet in self.facets}

    @property
    def marker_map(self) -> dict[str, MobilomeMarker]:
        return {marker.id: marker for marker in self.markers}

    @property
    def component_map(self) -> dict[str, InferenceComponent]:
        return {component.id: component for component in self.components}

    @property
    def inference_rule_map(self) -> dict[str, InferenceRule]:
        return {rule.id: rule for rule in self.inference_rules}


__all__ = [
    "SCHEMA_VERSION",
    "AnalysisParameters",
    "AnnotationProvenance",
    "AssessmentStatus",
    "CatalogResource",
    "DatabaseSource",
    "DisabledAggregateRule",
    "EvidenceReason",
    "ExternalHandoff",
    "Facet",
    "FeatureRef",
    "HandoffStatus",
    "HypothesisParticipant",
    "InferenceComponent",
    "InferenceRule",
    "MarkerMatcher",
    "MatchMode",
    "MissingComponent",
    "MobilomeDatabase",
    "MobilomeDatabaseError",
    "MobilomeError",
    "MobilomeHit",
    "MobilomeHypothesis",
    "MobilomeInputError",
    "MobilomeMarker",
    "MobilomeOutputError",
    "MobilomeParameterError",
    "MobilomeReport",
    "ParticipantRole",
    "ProvenanceSource",
    "RecordEvidence",
    "RepliconAssessment",
    "RepliconClass",
    "RepliconClassification",
    "RepliconInventory",
    "RepliconTopology",
    "SourceQualifier",
]
