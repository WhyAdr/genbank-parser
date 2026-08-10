"""Regression tests for the offline MEOR identifier audit."""

from __future__ import annotations

from pathlib import Path

import yaml

from scripts.audit_meor_identifiers import (
    ReferenceData,
    audit_markers,
    load_explorenz,
    load_fallback_ec,
    load_pfams,
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
    assert marker["kos"] == ["K00001"]
