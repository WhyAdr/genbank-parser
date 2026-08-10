"""Regression tests for the offline MEOR identifier audit."""

from __future__ import annotations

import hashlib
from pathlib import Path

import yaml

from genbank_parser.meor.audit import (
    Finding,
    ReferenceData,
    audit_markers,
    load_explorenz,
    load_fallback_ec,
    load_pfams,
    write_reference_manifest,
    write_reports,
)


def _markers_by_id() -> dict[str, dict[str, object]]:
    payload = yaml.safe_load(
        Path("src/genbank_parser/data/meor/markers.yaml").read_text(encoding="utf-8")
    )
    return {marker["id"]: marker for marker in payload["markers"]}


def test_identifier_remediation_assignments_are_regression_protected() -> None:
    markers = _markers_by_id()
    assert markers["prmABCD"]["kos"] == []
    assert markers["alkB"]["kos"] == ["K00496"]
    assert markers["CYP153"]["kos"] == []
    assert markers["rubAB"]["kos"] == []
    assert markers["almA"]["kos"] == []
    assert markers["dmpO_tomA"]["kos"] == []
    assert markers["assA_masD"]["kos"] == []
    assert markers["bssABC"]["kos"] == ["K07540"]
    assert markers["nmsA"]["kos"] == []
    assert markers["cmdA_ebdA"]["kos"] == []
    assert markers["rhlA"]["kos"] == []
    assert markers["rhlB"]["kos"] == []
    assert markers["rhlC"]["kos"] == []
    assert markers["srfA"]["kos"] == []
    assert markers["licA_D"]["kos"] == []
    assert markers["fengycin_iturin"]["kos"] == ["K15666", "K15667", "K15668"]
    assert markers["sfp"]["kos"] == []
    assert markers["wza_wzb_wzc"]["kos"] == ["K01991"]
    assert markers["pcaGH"]["ecs"] == ["1.13.11.3"]
    assert markers["dszC"]["ecs"] == ["1.14.14.21"]
    assert markers["alkJ"]["ecs"] == ["1.1.1.1"]
    assert markers["ncrA"]["kos"] == []
    assert markers["rhlRI"]["kos"] == ["K13061"]
    assert markers["rhlRI"]["ecs"] == []
    assert markers["etnE"]["kos"] == []


def test_explorenz_hist_action_and_class_heading_are_authoritative(
    tmp_path: Path,
) -> None:
    source = tmp_path / "explorenz.xml"
    source.write_text(
        """<mysqldump>
<table_data name="class"><row>
<field name="class">1</field><field name="subclass">14</field>
<field name="subsubclass">13</field><field name="heading">With &lt;small&gt;oxygen&lt;/small&gt;</field>
</row></table_data>
<table_data name="entry"><row>
<field name="ec_num">1.1.1.1</field><field name="accepted_name">alcohol dehydrogenase</field>
</row></table_data>
<table_data name="hist"><row>
<field name="ec_num">1.1.1.1</field><field name="action">transferred</field>
<field name="note">Use &lt;a href="query.php?ec=1.1.1.2"&gt;successor&lt;/a&gt;</field>
</row></table_data>
</mysqldump>""",
        encoding="utf-8",
    )
    entries, classes, history = load_explorenz(source)
    assert entries["1.1.1.1"]["accepted_name"] == "alcohol dehydrogenase"
    assert classes[("1", "14", "13")] == "With oxygen"
    assert history["1.1.1.1"]["action"] == "transferred"
    assert "1.1.1.2" in history["1.1.1.1"]["note"]


def test_ec_fallback_and_versioned_pfam_aliases(tmp_path: Path) -> None:
    ec_file = tmp_path / "enzyme.tsv"
    ec_file.write_text("1.1.1.1\talcohol dehydrogenase\n", encoding="utf-8")
    entries, history = load_fallback_ec(ec_file)
    assert entries["1.1.1.1"]["accepted_name"] == "alcohol dehydrogenase"
    assert history["1.1.1.1"]["action"] == "created"

    pfam_file = tmp_path / "Pfam-A.hmm"
    pfam_file.write_text(
        "HMMER3/f\nNAME  ABC_tran\nACC   PF00005.31\nDESC  ABC transporter\n//\n",
        encoding="utf-8",
    )
    pfams = load_pfams(pfam_file, {"PF00005"})
    assert pfams["PF00005"]["full_accession"] == "PF00005.31"
    assert pfams["PF00005.31"]["name"] == "ABC_tran"


def test_ko_ec_contradictions_are_reported_without_catalog_mutation() -> None:
    marker = {
        "id": "example",
        "name": "Example marker",
        "kos": ["K00001"],
        "ecs": ["2.1.1.-"],
        "cogs": [],
        "pfams": [],
    }
    references = ReferenceData(
        kegg={
            "K00001": {"definition": "example enzyme [EC:1.1.1.1]", "ecs": ["1.1.1.1"]}
        },
        cogs={},
        ec_entries={},
        ec_classes={},
        ec_history={},
        pfams={},
        ec_source="test",
    )
    findings = audit_markers([marker], references)
    ko = next(item for item in findings if item.identifier_type == "KO")
    consistency = next(
        item for item in findings if item.identifier_type == "CONSISTENCY"
    )
    assert ko.audit_status == "REVIEW"
    assert consistency.audit_status == "REVIEW"
    assert consistency.reference_status == "KO_EC_NO_COMPATIBLE_EC"
    assert marker["kos"] == ["K00001"]


def test_ko_ec_additional_scope_is_broader_not_a_contradiction() -> None:
    marker = {
        "id": "example",
        "name": "Example marker",
        "kos": ["K00001"],
        "ecs": ["1.1.1.1"],
        "cogs": [],
        "pfams": [],
    }
    references = ReferenceData(
        kegg={
            "K00001": {
                "definition": "example enzyme [EC:1.1.1.1 1.1.1.2]",
                "ecs": ["1.1.1.1", "1.1.1.2"],
            }
        },
        cogs={},
        ec_entries={},
        ec_classes={},
        ec_history={},
        pfams={},
        ec_source="test",
    )
    findings = audit_markers([marker], references)
    ko = next(item for item in findings if item.identifier_type == "KO")
    consistency = next(
        item for item in findings if item.identifier_type == "CONSISTENCY"
    )
    assert ko.audit_status == "PRESENT_UNREVIEWED"
    assert consistency.audit_status == "PRESENT_UNREVIEWED"
    assert consistency.reference_status == "KO_EC_COMPATIBLE_BUT_BROADER"


def test_reference_manifest_is_hashed_and_posix_normalized(tmp_path: Path) -> None:
    references: dict[str, Path] = {}
    for label in ("kegg", "cog", "ec", "pfam"):
        source = tmp_path / f"{label}.dat"
        source.write_bytes(label.encode("ascii"))
        references[label] = source
    manifest = tmp_path / "reference_manifest.tsv"
    write_reference_manifest(manifest, references)
    text = manifest.read_text(encoding="utf-8")
    assert "\\" not in text
    assert hashlib.sha256(b"kegg").hexdigest() in text
    assert text == manifest.read_text(encoding="utf-8")


def test_report_exposes_broader_ko_ec_qualifications(tmp_path: Path) -> None:
    reference_paths: dict[str, Path] = {}
    for label in ("kegg", "cog", "ec", "pfam"):
        source = tmp_path / f"{label}.dat"
        source.write_bytes(label.encode("ascii"))
        reference_paths[label] = source
    markers = tmp_path / "markers.yaml"
    markers.write_text("schema_version: 1\nmarkers: []\n", encoding="utf-8")
    references = ReferenceData(
        kegg={},
        cogs={},
        ec_entries={},
        ec_classes={},
        ec_history={},
        pfams={},
        ec_source=str(reference_paths["ec"]),
    )
    findings = [
        Finding(
            "tmo_tod",
            "CONSISTENCY",
            "tmo_tod",
            "KO_EC_COMPATIBLE_BUT_BROADER",
            "K00001: 1.1.1.2",
            "1.1.1.-",
            "PRESENT_UNREVIEWED",
            "additional EC scope is broader",
        ),
        Finding(
            "alkB",
            "CONSISTENCY",
            "alkB",
            "KO_EC_EXACT_COMPATIBLE",
            "K00496: 1.14.15.3",
            "1.14.15.3",
            "PRESENT_UNREVIEWED",
            "all embedded KO ECs are compatible",
        ),
    ]
    write_reports(tmp_path / "audit", findings, references, markers, reference_paths)
    summary = (tmp_path / "audit" / "audit_summary.md").read_text(encoding="utf-8")
    assert "KO_EC_COMPATIBLE_BUT_BROADER" not in summary
    assert "`tmo_tod` | `CONSISTENCY` | `tmo_tod`" in summary
    assert "additional EC scope is broader" in summary
