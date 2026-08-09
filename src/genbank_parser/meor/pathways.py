"""Genome-level MEOR pathway completeness evaluation."""
from __future__ import annotations

from .models import MeorDatabase, MeorHit, MeorPathwayResult


def evaluate_meor_pathways(
    hits: list[MeorHit] | tuple[MeorHit, ...], database: MeorDatabase
) -> list[MeorPathwayResult]:
    marker_ids = {hit.marker_id for hit in hits}
    results: list[MeorPathwayResult] = []
    for pathway in database.pathways:
        present = 0
        steps: list[dict[str, object]] = []
        for step in pathway.steps:
            is_present = any(marker_id in marker_ids for marker_id in step.any_of)
            present += int(is_present)
            steps.append(
                {
                    "step": step.name,
                    "status": "PRESENT" if is_present else "ABSENT",
                    "target_markers": list(step.any_of),
                }
            )
        total = len(pathway.steps)
        completeness = round(100.0 * present / total, 1) if total else 0.0
        results.append(
            MeorPathwayResult(
                pathway_id=pathway.id,
                pathway=pathway.name,
                category_id=pathway.category_id,
                completeness_pct=completeness,
                steps_present=present,
                total_steps=total,
                steps=tuple(steps),
            )
        )
    return results
