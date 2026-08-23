"""Cautious record-local and cross-record mobilome inference tests."""

from pathlib import Path

from genbank_parser.mobilome import analyze_mobilome

FORBIDDEN = (
    "obligate co-transfer",
    "confirmed conjugative",
    "satellite plasmid",
    "theta replication",
    "rolling-circle replication",
    "confirmed phagemid",
)


def test_cross_record_helper_hypothesis_is_tentative_and_structured() -> None:
    report = analyze_mobilome(Path("tests/fixtures/mobilome_inference.gb"))

    assert len(report.cross_record_hypotheses) == 1
    hypothesis = report.cross_record_hypotheses[0]
    assert hypothesis.summary == "Possible helper-dependent mobilization"
    assert [participant.role for participant in hypothesis.participants] == [
        "target",
        "helper",
    ]
    assert (
        hypothesis.participants[0].record_index
        != hypothesis.participants[1].record_index
    )
    assert "Co-transfer of both records is not implied." in hypothesis.limitations


def test_toxin_antitoxin_pairs_and_empty_record_are_calibrated() -> None:
    report = analyze_mobilome(Path("tests/fixtures/mobilome_evidence.gb"))
    summaries = [item.summary for item in report.replicons[0].hypotheses]

    assert "Type III ToxIN annotation-pair candidate" in summaries
    assert "HEPN/MNT annotation-pair candidate" in summaries
    assert all(
        forbidden not in "\n".join(summaries).casefold() for forbidden in FORBIDDEN
    )

    empty = analyze_mobilome(
        Path("tests/fixtures/mobilome_inference.gb")
    ).scanned_replicons[2]
    empty_summaries = [item.summary for item in empty.hypotheses]
    assert empty_summaries == [
        "Replication evidence insufficient; mechanism unresolved"
    ]
    assert "satellite" not in " ".join(empty_summaries).casefold()
