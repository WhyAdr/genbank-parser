"""Native, typed, annotation-supported mobilome evidence analysis."""

from __future__ import annotations

import hashlib
import re
from pathlib import Path

from .. import __version__
from ..io import read_genbank
from ..model import GenBankDocument
from .database import load_mobilome_database
from .inference import infer_cross_record_hypotheses, infer_replicon_hypotheses
from .models import (
    AnalysisParameters,
    AnnotationProvenance,
    MobilomeDatabase,
    MobilomeDatabaseError,
    MobilomeError,
    MobilomeInputError,
    MobilomeOutputError,
    MobilomeParameterError,
    MobilomeReport,
    RepliconAssessment,
    RepliconInventory,
)
from .replicons import inventory_replicons
from .scanner import scan_mobilome_features, validate_min_evidence

_VALID_INCLUDE = frozenset({"all", "chromosome", "plasmid", "unknown"})


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
    except OSError as exc:
        raise MobilomeInputError(
            f"Could not read input GenBank file {path}: {exc}"
        ) from exc
    return digest.hexdigest()


def _validate_include(include: str) -> None:
    if include not in _VALID_INCLUDE:
        raise MobilomeParameterError(
            "include must be one of 'all', 'chromosome', 'plasmid', or 'unknown'"
        )


def _annotation_provenance(
    document: GenBankDocument,
) -> tuple[AnnotationProvenance, ...]:
    """Capture only input-declared annotation-pipeline context when available."""

    grouped: dict[tuple[str | None, str | None, str | None], list[int]] = {}
    for record_index, record in enumerate(document.records, 1):
        values = [
            str(value)
            for value in record.annotations.values()
            if isinstance(value, str)
        ]
        comments = "\n".join(values)
        name: str | None = None
        version: str | None = None
        database_version: str | None = None
        if re.search(r"\bBakta\b", comments, re.IGNORECASE):
            name = "Bakta"
            version_match = re.search(
                r"\bBakta\s+(?:v(?:ersion)?\s*)?([0-9][A-Za-z0-9._-]*)",
                comments,
                re.IGNORECASE,
            )
            if version_match is not None:
                version = version_match.group(1)
            database_match = re.search(
                r"\b(?:database|db)\s*(?:version)?\s*[:= ]\s*([A-Za-z0-9._-]+)",
                comments,
                re.IGNORECASE,
            )
            if database_match is not None:
                database_version = database_match.group(1)
        key = (name, version, database_version)
        grouped.setdefault(key, []).append(record_index)
    return tuple(
        AnnotationProvenance(
            record_indices=tuple(record_indices),
            name=name,
            version=version,
            database_version=database_version,
            evidence_fields=("record_annotations",),
        )
        for (name, version, database_version), record_indices in sorted(
            grouped.items(), key=lambda item: item[1]
        )
    )


def analyze_mobilome(
    filepath: str | Path,
    *,
    include: str = "all",
    min_evidence: int = 1,
    database: MobilomeDatabase | None = None,
) -> MobilomeReport:
    """Return a pure annotation-supported mobilome evidence report.

    The input is parsed exactly once through :func:`read_genbank`.  All records
    are inventoried and scanned; ``include`` only filters detailed report
    presentation after cross-record evaluation is complete.
    """

    _validate_include(include)
    validate_min_evidence(min_evidence)
    path = Path(filepath)
    if not path.is_file():
        raise MobilomeInputError(
            f"Input GenBank file does not exist or is not a file: {path}"
        )
    source_before = _sha256(path)
    try:
        document = read_genbank(path)
    except Exception as exc:  # Biopython parser exceptions vary by malformed input.
        raise MobilomeInputError(f"Could not parse input GenBank file: {exc}") from exc
    source_after = _sha256(path)
    if source_after != source_before:
        raise MobilomeInputError(
            "Input changed while parsing; analysis was not retained"
        )
    if not document.records:
        raise MobilomeInputError("Input GenBank file contains no parsed records")
    active_database = database or load_mobilome_database()
    inventory = inventory_replicons(document)
    hits = scan_mobilome_features(
        document.all_features, active_database, min_evidence=min_evidence
    )
    all_assessments: list[RepliconAssessment] = []
    for member in inventory:
        record_hits = tuple(
            hit for hit in hits if hit.feature.record_index == member.record_index
        )
        all_assessments.append(
            RepliconAssessment(
                inventory=member,
                hits=record_hits,
                hypotheses=infer_replicon_hypotheses(
                    member, record_hits, active_database
                ),
            )
        )
    scanned_replicons = tuple(all_assessments)
    cross_record_hypotheses = infer_cross_record_hypotheses(
        scanned_replicons, active_database
    )
    detailed_replicons = tuple(
        assessment
        for assessment in scanned_replicons
        if include == "all" or assessment.inventory.classification.label == include
    )
    return MobilomeReport(
        source_file=path.name,
        source_sha256=source_after,
        tool_version=__version__,
        catalog_version=active_database.catalog_version,
        inference_version=active_database.inference_version,
        resources=active_database.resources,
        parameters=AnalysisParameters(
            include=include,  # type: ignore[arg-type]
            min_evidence=min_evidence,
            database_source=active_database.database_source,
        ),
        annotation_provenance=_annotation_provenance(document),
        inventory=inventory,
        replicons=detailed_replicons,
        cross_record_hypotheses=cross_record_hypotheses,
        disabled_aggregate_rules=active_database.disabled_aggregate_rules,
        handoffs=active_database.handoffs,
        limitations=(
            "This report summarizes input annotations and configured evidence; it does not establish phenotype, expression, transfer, replicon ancestry, replication mechanism, or element identity.",
            "Absent annotations are not evidence that a biological function is absent.",
        ),
        scanned_replicons=scanned_replicons,
    )


__all__ = [
    "MobilomeDatabase",
    "MobilomeDatabaseError",
    "MobilomeError",
    "MobilomeInputError",
    "MobilomeOutputError",
    "MobilomeParameterError",
    "MobilomeReport",
    "RepliconInventory",
    "analyze_mobilome",
    "infer_cross_record_hypotheses",
    "infer_replicon_hypotheses",
    "inventory_replicons",
    "load_mobilome_database",
    "scan_mobilome_features",
]
