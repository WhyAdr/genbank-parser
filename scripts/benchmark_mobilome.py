#!/usr/bin/env python3
"""Optional real-genome calibration for the mobilome analysis (non-CI).

Runs the packaged mobilome analysis over a local real GenBank file (default:
``PF_NNT_reoriented.gbff`` at the repository root) and reports runtime,
rendered output volume, and the calibrated expectations from the mobilome
integration plan (section 16.9). The file is intentionally gitignored, so the
calibration must skip clearly when the file is absent and must never stage it.
"""

from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

from genbank_parser.mobilome import analyze_mobilome
from genbank_parser.mobilome.report import serialize_mobilome_report

FORBIDDEN_CLAIMS = (
    "obligate co-transfer",
    "confirmed conjugative",
    "satellite plasmid",
    "theta replication",
    "rolling-circle replication",
    "confirmed phagemid",
    "oriv",
)

DEFAULT_INPUT = Path("PF_NNT_reoriented.gbff")
DEFAULT_MANIFEST = Path("tests/data/PF_NNT_reoriented.manifest.json")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_calibration_manifest(path: Path = DEFAULT_MANIFEST) -> dict[str, object]:
    """Load and minimally validate the committed calibration manifest."""

    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"Could not load calibration manifest {path}: {exc}") from exc
    if not isinstance(payload, dict) or payload.get("schema_version") != 1:
        raise ValueError(f"Calibration manifest {path} must use schema_version 1")
    fixture_id = payload.get("fixture_id")
    expected = payload.get("sha256")
    source_record_id = payload.get("source_record_id")
    transformation_provenance = payload.get("transformation_provenance")
    if not isinstance(fixture_id, str) or not fixture_id:
        raise ValueError("Calibration manifest fixture_id must be a non-empty string")
    if not isinstance(source_record_id, str) or not source_record_id:
        raise ValueError(
            "Calibration manifest source_record_id must be a non-empty string"
        )
    if transformation_provenance not in {"complete", "incomplete"}:
        raise ValueError(
            "Calibration manifest transformation_provenance must be "
            "'complete' or 'incomplete'"
        )
    if not isinstance(expected, str) or len(expected) != 64:
        raise ValueError(
            "Calibration manifest sha256 must be a 64-character hex string"
        )
    try:
        int(expected, 16)
    except ValueError as exc:
        raise ValueError("Calibration manifest sha256 must be hexadecimal") from exc
    return payload


def verify_calibration_input(
    path: Path,
    manifest_path: Path = DEFAULT_MANIFEST,
) -> dict[str, object]:
    """Verify exact input bytes before running real-genome calibration."""

    manifest = load_calibration_manifest(manifest_path)
    observed = _sha256(path)
    expected = str(manifest["sha256"]).lower()
    if observed != expected:
        raise ValueError(
            "Calibration fixture hash mismatch: "
            f"expected {expected}, observed {observed}"
        )
    return {**manifest, "observed_sha256": observed}


def calibrate(
    path: Path,
    *,
    manifest_path: Path = DEFAULT_MANIFEST,
) -> int:
    metadata = verify_calibration_input(path, manifest_path)
    started = time.perf_counter()
    report = analyze_mobilome(path)
    elapsed = time.perf_counter() - started

    inventories = list(report.inventory)
    plasmid_evidence = [
        item
        for item in inventories
        if any(
            qualifier.key.casefold() == "plasmid"
            for qualifier in item.source_qualifiers
        )
    ]
    pseudo_cds = sum(item.pseudo_cds_count for item in inventories)
    summaries = [
        hypothesis.summary
        for assessment in report.replicons
        for hypothesis in assessment.hypotheses
    ]
    all_text = " ".join(summaries).casefold()
    scanned = report.scanned_replicons or report.replicons
    hits = tuple(hit for assessment in scanned for hit in assessment.hits)
    reasons = tuple(reason for hit in hits for reason in hit.reasons)
    summary = {
        "raw_hits": len(hits),
        "eligible_hits": sum(1 for hit in hits if hit.eligible_for_inference),
        "raw_reasons": len(reasons),
        "eligible_reasons": sum(1 for reason in reasons if reason.eligible),
        "hypotheses": len(
            tuple(
                hypothesis
                for assessment in scanned
                for hypothesis in assessment.hypotheses
            )
            + report.cross_record_hypotheses
        ),
    }
    volume = {
        name: len(payload)
        for name, payload in (
            ("json", serialize_mobilome_report(report, "json")),
            ("tsv", serialize_mobilome_report(report, "tsv")),
            ("text", serialize_mobilome_report(report, "text")),
        )
    }

    print(f"calibration fixture:    {metadata['fixture_id']}")
    print(f"input:                  {path}")
    print("hash verification:      verified")
    print(f"source record ID:       {metadata['source_record_id']}")
    print(
        "transformation provenance: "
        f"{metadata['transformation_provenance']}"
    )
    print(f"input SHA-256:          {metadata['observed_sha256']}")
    print(f"records inventoried:    {len(inventories)}")
    print(f"source-level plasmids:  {len(plasmid_evidence)}")
    print(f"pseudo CDS retained:    {pseudo_cds}")
    print(
        f"raw hits / eligible:    {summary['raw_hits']} / {summary['eligible_hits']} "
        f"(reasons {summary['raw_reasons']} / {summary['eligible_reasons']})"
    )
    print(f"hypotheses:             {summary['hypotheses']}")
    print(f"runtime:                {elapsed:.3f} s")
    for name, size in volume.items():
        print(f"rendered {name:<14}  {size} bytes")

    failures: list[str] = []
    if len(inventories) != 5:
        failures.append(f"expected 5 records, observed {len(inventories)}")
    if len(plasmid_evidence) != 4:
        failures.append(
            f"expected 4 records with explicit source-level plasmid evidence, "
            f"observed {len(plasmid_evidence)}"
        )
    if pseudo_cds != 21:
        failures.append(f"expected 21 canonical pseudo CDS, observed {pseudo_cds}")
    for claim in FORBIDDEN_CLAIMS:
        if claim in all_text:
            failures.append(f"forbidden categorical claim present: {claim!r}")
    for item in inventories:
        label = item.classification.label
        if label in ("chromosome", "plasmid") and not (
            item.classification.supporting_evidence
        ):
            failures.append(
                f"record {item.record_id} classified as {label} without "
                "explicit annotation evidence (circularity alone is not evidence)"
            )
    for assessment in report.replicons:
        for hypothesis in assessment.hypotheses:
            if hypothesis.rule_id == "replication_annotation_candidate" and (
                hypothesis.summary
                not in (
                    "Replication annotation candidate with unresolved mechanism",
                    "Replication evidence insufficient; mechanism unresolved",
                )
            ):
                failures.append(
                    "replication hypothesis selects a mechanism: "
                    f"{hypothesis.summary!r}"
                )

    if failures:
        print("\nCALIBRATION FAILURES:")
        for failure in failures:
            print(f"  - {failure}")
        return 1
    print("\ncalibration: all section 16.9 expectations met")
    return 0


def benchmark(path: Path) -> int:
    """Compatibility wrapper for callers of the pre-0.8.2 script API."""

    return calibrate(path)


def main(argv: list[str]) -> int:
    path = Path(argv[1]) if len(argv) > 1 else DEFAULT_INPUT
    if not path.is_file():
        print(
            f"skip: real-genome calibration input {path} is not present "
            "(it is gitignored by design and never staged)"
        )
        return 0
    try:
        return calibrate(path)
    except ValueError as exc:
        print(f"calibration error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
