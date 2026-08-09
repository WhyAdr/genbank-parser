"""Source-aware, parity-preserving MEOR marker matching."""
from __future__ import annotations

import re

from ..io import extract_xrefs
from ..model import GenBankFeature
from .models import MarkerMatch, MeorMarker


def match_feature_to_marker(
    feature: GenBankFeature, marker: MeorMarker
) -> MarkerMatch | None:
    """Return the highest-precedence legacy-compatible match for one marker."""
    structured = extract_xrefs(feature, include_notes=False)
    gene = feature.gene.strip()
    product = feature.product.lower().strip()
    notes = " ".join(feature.get_quals("note")).lower()

    for ko in structured["kegg_kos"]:
        if ko in marker.kos:
            return MarkerMatch(3, "kegg", ko)
    for ec in structured["ec_numbers"]:
        if ec in marker.ecs:
            return MarkerMatch(3, "ec", ec)
    if gene:
        for pattern in marker.gene_patterns:
            if re.search(pattern, gene, re.IGNORECASE):
                return MarkerMatch(3, "gene", gene)
    if product:
        for pattern in marker.product_patterns:
            if re.search(pattern, product, re.IGNORECASE):
                return MarkerMatch(2, "product", pattern)
    if notes:
        for pattern in marker.product_patterns:
            if re.search(pattern, notes, re.IGNORECASE):
                return MarkerMatch(1, "note", pattern)
        for pattern in marker.note_patterns:
            if re.search(pattern, notes, re.IGNORECASE):
                return MarkerMatch(1, "note", pattern)
        for ec in marker.ecs:
            ec_lower = ec.lower()
            if f"ec:{ec_lower}" in notes or f"ec_number:{ec_lower}" in notes:
                return MarkerMatch(1, "note_ec", ec)
        for ko in marker.kos:
            if ko.lower() in notes:
                return MarkerMatch(1, "note_kegg", ko)
    return None
