"""Typed models for MEOR marker discovery and reporting."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal


@dataclass(frozen=True)
class MeorCategory:
    id: str
    title: str
    role: str


@dataclass(frozen=True)
class MeorMarker:
    id: str
    name: str
    category_id: str
    gene_patterns: tuple[str, ...] = ()
    ecs: tuple[str, ...] = ()
    kos: tuple[str, ...] = ()
    product_patterns: tuple[str, ...] = ()
    note_patterns: tuple[str, ...] = ()
    sources: tuple[str, ...] = ()


@dataclass(frozen=True)
class MeorPathwayStep:
    name: str
    any_of: tuple[str, ...]


@dataclass(frozen=True)
class MeorPathway:
    id: str
    name: str
    category_id: str
    steps: tuple[MeorPathwayStep, ...]


@dataclass(frozen=True)
class MeorDatabase:
    categories: tuple[MeorCategory, ...]
    markers: tuple[MeorMarker, ...]
    pathways: tuple[MeorPathway, ...]
    provenance: dict[str, Any] = field(default_factory=dict, compare=False)

    @property
    def category_map(self) -> dict[str, MeorCategory]:
        return {category.id: category for category in self.categories}

    @property
    def marker_map(self) -> dict[str, MeorMarker]:
        return {marker.id: marker for marker in self.markers}


EvidenceType = Literal[
    "kegg", "ec", "gene", "product", "note", "note_kegg", "note_ec"
]


@dataclass(frozen=True)
class MarkerMatch:
    weight: int
    evidence_type: EvidenceType
    evidence_value: str

    @property
    def reason(self) -> str:
        prefixes = {
            "kegg": "KEGG:",
            "ec": "EC:",
            "gene": "gene:",
            "product": "product:",
            "note": "note:",
            "note_kegg": "note:KEGG:",
            "note_ec": "note:EC:",
        }
        return f"{prefixes[self.evidence_type]}{self.evidence_value}"


@dataclass(frozen=True)
class MeorHit:
    feature_index: int
    contig: str
    locus_tag: str
    gene: str
    product: str
    start: int
    end: int
    strand: str
    category_id: str
    category_title: str
    marker_id: str
    marker_name: str
    weight: int
    evidence_type: EvidenceType
    evidence_value: str

    @property
    def reason(self) -> str:
        return MarkerMatch(
            self.weight, self.evidence_type, self.evidence_value
        ).reason

    def to_dict(self) -> dict[str, Any]:
        return {
            "contig": self.contig,
            "locus_tag": self.locus_tag,
            "gene": self.gene,
            "product": self.product,
            "start": self.start,
            "end": self.end,
            "strand": self.strand,
            "category_key": self.category_id,
            "category_title": self.category_title,
            "marker_id": self.marker_id,
            "marker_name": self.marker_name,
            "weight": self.weight,
            "reason": self.reason,
            "feature_index": self.feature_index,
            "evidence_type": self.evidence_type,
            "evidence_value": self.evidence_value,
        }


@dataclass(frozen=True)
class MeorCluster:
    cluster_id: str
    contig: str
    start: int
    end: int
    span_bp: int
    gene_count: int
    hit_count: int
    categories: tuple[str, ...]
    locus_tags: tuple[str, ...]
    genes: tuple[str, ...]
    members: tuple[MeorHit, ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "cluster_id": self.cluster_id,
            "contig": self.contig,
            "start": self.start,
            "end": self.end,
            "span_bp": self.span_bp,
            "gene_count": self.gene_count,
            "hit_count": self.hit_count,
            "categories": list(self.categories),
            "locus_tags": list(self.locus_tags),
            "genes": list(self.genes),
            "members": [member.to_dict() for member in self.members],
        }


@dataclass(frozen=True)
class MeorPathwayResult:
    pathway_id: str
    pathway: str
    category_id: str
    completeness_pct: float
    steps_present: int
    total_steps: int
    steps: tuple[dict[str, Any], ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "pathway": self.pathway,
            "completeness_pct": self.completeness_pct,
            "steps_present": self.steps_present,
            "total_steps": self.total_steps,
            "steps": [dict(step) for step in self.steps],
            "pathway_id": self.pathway_id,
            "category_key": self.category_id,
        }


@dataclass(frozen=True)
class MeorReport:
    source_file: str
    total_features: int
    parameters: dict[str, int]
    hits: tuple[MeorHit, ...]
    clusters: tuple[MeorCluster, ...]
    pathways: tuple[MeorPathwayResult, ...]
    karyograms: dict[str, dict[str, Any]]
    categories: tuple[MeorCategory, ...] = field(repr=False)

    @property
    def total_hits(self) -> int:
        return len(self.hits)

    @property
    def total_clusters(self) -> int:
        return len(self.clusters)

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": "gbparse.meor.v1",
            "tool": "gbparse",
            "analysis": "meor",
            "file": self.source_file,
            "total_features": self.total_features,
            "total_hits": self.total_hits,
            "total_clusters": self.total_clusters,
            "parameters": dict(self.parameters),
            "hits": [hit.to_dict() for hit in self.hits],
            "clusters": [cluster.to_dict() for cluster in self.clusters],
            "pathways": [pathway.to_dict() for pathway in self.pathways],
            "karyograms": self.karyograms,
        }
