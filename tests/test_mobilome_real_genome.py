"""Opt-in real-genome calibration test for the mobilome analysis.

The real input ``PF_NNT_reoriented.gbff`` is intentionally gitignored (plan
section 16.9): the test skips clearly when the file is absent and must never
stage it.  When present, it asserts the calibrated expectations documented in
the mobilome integration plan.
"""

from pathlib import Path

import pytest

from genbank_parser.mobilome import analyze_mobilome

REAL_GENOME = Path("PF_NNT_reoriented.gbff")

pytestmark = pytest.mark.skipif(
    not REAL_GENOME.is_file(),
    reason="real-genome input PF_NNT_reoriented.gbff is not present (gitignored by design)",
)

FORBIDDEN_CLAIMS = (
    "obligate co-transfer",
    "confirmed conjugative",
    "satellite plasmid",
    "theta replication",
    "rolling-circle replication",
    "confirmed phagemid",
    "oriv",
)


def test_real_genome_meets_section_16_9_calibration() -> None:
    report = analyze_mobilome(REAL_GENOME)

    inventories = list(report.inventory)
    assert len(inventories) == 5

    plasmid_evidence = [
        item
        for item in inventories
        if any(
            qualifier.key.casefold() == "plasmid"
            for qualifier in item.source_qualifiers
        )
    ]
    assert len(plasmid_evidence) == 4

    assert sum(item.pseudo_cds_count for item in inventories) == 21

    for item in inventories:
        label = item.classification.label
        if label in ("chromosome", "plasmid"):
            # No record is classified from circularity alone.
            assert item.classification.supporting_evidence, item.record_id

    summaries = " ".join(
        hypothesis.summary
        for assessment in report.replicons
        for hypothesis in assessment.hypotheses
    ).casefold()
    for claim in FORBIDDEN_CLAIMS:
        # No generic product produces a categorical AMR, VF, phagemid,
        # satellite, or transfer claim, and the chromosome's generic
        # rep_origin is never called an oriV or a replication mechanism.
        assert claim not in summaries, claim
    for assessment in report.replicons:
        for hypothesis in assessment.hypotheses:
            if hypothesis.rule_id == "replication_annotation_candidate":
                assert hypothesis.summary in (
                    "Replication annotation candidate with unresolved mechanism",
                    "Replication evidence insufficient; mechanism unresolved",
                )
