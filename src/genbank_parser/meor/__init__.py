"""MEOR, hydrocarbon-degradation, and biosurfactant analysis API."""
from __future__ import annotations

from pathlib import Path

from ..io import read_genbank
from .database import MeorDatabaseError, load_meor_database
from .karyogram import generate_meor_karyograms
from .models import (
    MarkerMatch,
    MeorCategory,
    MeorCluster,
    MeorDatabase,
    MeorHit,
    MeorMarker,
    MeorPathway,
    MeorPathwayResult,
    MeorReport,
)
from .pathways import evaluate_meor_pathways
from .scanner import scan_meor_features
from .spatial import cluster_meor_hits


def analyze_meor(
    filepath: str | Path,
    *,
    min_weight: int = 1,
    max_gap: int = 200,
    window_size: int = 50_000,
    marker_database: MeorDatabase | None = None,
) -> MeorReport:
    """Analyze annotation-supported MEOR genomic potential without printing."""
    if min_weight not in (1, 2, 3):
        raise ValueError("min_weight must be 1, 2, or 3")
    if max_gap < 0:
        raise ValueError("max_gap must be non-negative")
    if window_size <= 0:
        raise ValueError("window_size must be positive")
    database = marker_database or load_meor_database()
    document = read_genbank(filepath)
    hits = scan_meor_features(
        document.all_features, database, min_weight=min_weight
    )
    clusters = cluster_meor_hits(hits, max_gap=max_gap)
    pathways = evaluate_meor_pathways(hits, database)
    karyograms = generate_meor_karyograms(
        document.records, hits, window_size=window_size
    )
    return MeorReport(
        source_file=Path(filepath).name,
        total_features=document.total_features,
        parameters={
            "min_weight": min_weight,
            "max_gap": max_gap,
            "window_size": window_size,
        },
        hits=tuple(hits),
        clusters=tuple(clusters),
        pathways=tuple(pathways),
        karyograms=karyograms,
        categories=database.categories,
    )


__all__ = [
    "MarkerMatch",
    "MeorCategory",
    "MeorCluster",
    "MeorDatabase",
    "MeorDatabaseError",
    "MeorHit",
    "MeorMarker",
    "MeorPathway",
    "MeorPathwayResult",
    "MeorReport",
    "analyze_meor",
    "cluster_meor_hits",
    "evaluate_meor_pathways",
    "generate_meor_karyograms",
    "load_meor_database",
    "scan_meor_features",
]
