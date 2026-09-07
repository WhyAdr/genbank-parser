#!/usr/bin/env python3
"""Optional real-genome benchmark for the mobilome analysis (non-CI).

Runs the packaged mobilome analysis over a local real GenBank file (default:
``PF_NNT_reoriented.gbff`` at the repository root) and reports runtime,
rendered output volume, and the calibrated expectations from the mobilome
integration plan (section 16.9).  The file is intentionally gitignored, so the
benchmark must skip clearly when the file is absent and must never stage it.
"""

from __future__ import annotations

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


def benchmark(path: Path) -> int:
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

    print(f"input:                  {path}")
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


def main(argv: list[str]) -> int:
    path = Path(argv[1]) if len(argv) > 1 else DEFAULT_INPUT
    if not path.is_file():
        print(
            f"skip: real-genome input {path} is not present "
            "(it is gitignored by design and never staged)"
        )
        return 0
    return benchmark(path)


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
