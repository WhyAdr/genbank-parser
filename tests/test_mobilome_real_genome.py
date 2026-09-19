"""Opt-in real-genome calibration tests for the mobilome analysis.

The real input ``PF_NNT_reoriented.gbff`` is intentionally gitignored (plan
section 16.9): the test skips clearly when the file is absent and must never
stage it. When present, its manifest-locked bytes are verified before the
calibrated expectations documented in the mobilome integration plan run.
"""

import importlib.util
from pathlib import Path

import pytest

from genbank_parser.mobilome import analyze_mobilome

_SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "benchmark_mobilome.py"
_SCRIPT_SPEC = importlib.util.spec_from_file_location(
    "gbparse_benchmark_mobilome", _SCRIPT_PATH
)
if _SCRIPT_SPEC is None or _SCRIPT_SPEC.loader is None:
    raise ImportError(f"Could not load calibration script from {_SCRIPT_PATH}")
_SCRIPT_MODULE = importlib.util.module_from_spec(_SCRIPT_SPEC)
_SCRIPT_SPEC.loader.exec_module(_SCRIPT_MODULE)
main = _SCRIPT_MODULE.main
verify_calibration_input = _SCRIPT_MODULE.verify_calibration_input

REAL_GENOME = Path("PF_NNT_reoriented.gbff")
MANIFEST = Path("tests/data/PF_NNT_reoriented.manifest.json")

REAL_GENOME_SKIP = pytest.mark.skipif(
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


@REAL_GENOME_SKIP
def test_real_genome_meets_section_16_9_calibration() -> None:
    metadata = verify_calibration_input(REAL_GENOME, MANIFEST)
    assert metadata["fixture_id"] == "PF_NNT_reoriented"
    assert metadata["observed_sha256"] == metadata["sha256"]

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


@REAL_GENOME_SKIP
def test_altered_calibration_bytes_fail_before_analysis(tmp_path: Path) -> None:
    altered = tmp_path / "PF_NNT_reoriented.gbff"
    altered.write_bytes(REAL_GENOME.read_bytes() + b"\n")

    with pytest.raises(ValueError, match="hash mismatch"):
        verify_calibration_input(altered, MANIFEST)


def test_missing_calibration_input_skips_cleanly(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    missing = tmp_path / "missing.gbff"

    assert main(["benchmark_mobilome.py", str(missing)]) == 0
    assert "skip: real-genome calibration input" in capsys.readouterr().out
