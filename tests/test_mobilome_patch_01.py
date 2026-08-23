"""Post-0.6.0 mobilome hardening and adversarial contract regressions."""

from __future__ import annotations

import json
from importlib import resources
from pathlib import Path

import pytest
import yaml
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation

from genbank_parser import read_genbank
from genbank_parser.cli import main
from genbank_parser.mobilome import analyze_mobilome, load_mobilome_database
from genbank_parser.mobilome.inference import (
    infer_cross_record_hypotheses,
    infer_replicon_hypotheses,
    record_inference_limitations,
)
from genbank_parser.mobilome.models import (
    MobilomeDatabaseError,
    MobilomeInputError,
    MobilomeOutputError,
    RepliconAssessment,
)
from genbank_parser.mobilome.replicons import inventory_replicons
from genbank_parser.mobilome.report import (
    render_text,
    serialize_mobilome_report,
    write_mobilome_report,
)
from genbank_parser.mobilome.scanner import scan_mobilome_features
from genbank_parser.model import GenBankDocument, GenBankFeature, GenBankRecord

FORBIDDEN_CLAIMS = (
    "obligate co-mobilisation",
    "obligate co-mobilization",
    "obligate co-transfer",
    "guaranteed transfer",
    "confirmed conjugative",
    "diagnostic pxo",
    "proves pxo ancestry",
    "confirmed phagemid",
    "confirmed resistance",
    "confirmed virulence",
    "active communication",
    "satellite plasmid",
    "evolutionary transition",
    "theta replication",
    "rolling-circle replication",
)


def _genbank(
    tmp_path: Path,
    *,
    name: str,
    products: tuple[str, ...] = (),
    source: str | None = None,
    comment: str | None = None,
    pseudo_rep: bool = False,
    crispr: bool = False,
    cas: bool = False,
) -> Path:
    lines = [
        f"LOCUS       {name:<16} 800 bp    DNA     circular BCT 23-AUG-2026".ljust(79),
        f"DEFINITION  Synthetic {name} mobilome patch record.",
        f"ACCESSION   {name}",
        f"VERSION     {name}.1",
    ]
    if comment is not None:
        lines.extend(
            [
                "COMMENT     ##Patch-START##",
                f"            Note :: {comment}",
                "            ##Patch-END##",
            ]
        )
    lines.extend(
        ["FEATURES             Location/Qualifiers", "     source          1..800"]
    )
    lines.append('                     /organism="Synthetic bacterium"')
    if source is not None:
        lines.append(f'                     /{source}="patch"')
    if crispr:
        lines.extend(
            [
                "     repeat_region   100..120",
                '                     /rpt_type="CRISPR"',
            ]
        )
    if cas:
        lines.extend(
            [
                "     CDS             140..190",
                '                     /gene="cas9"',
                '                     /product="CRISPR-associated protein Cas9"',
                '                     /translation="MMMM"',
            ]
        )
    for index, product in enumerate(products, 1):
        start = 220 + (index - 1) * 55
        end = start + 40
        gene = "repA" if "replication initiation" in product.casefold() else "patch"
        lines.extend(
            [
                f"     CDS             {start}..{end}",
                f'                     /locus_tag="PATCH_{index:04d}"',
                f'                     /gene="{gene}"',
                f'                     /product="{product}"',
                '                     /translation="MMMM"',
            ]
        )
        if pseudo_rep and index == 1:
            lines.append("                     /pseudo")
    if pseudo_rep and not products:
        lines.extend(
            [
                "     CDS             220..260",
                '                     /gene="repA"',
                '                     /product="replication initiation protein"',
                "                     /pseudo",
            ]
        )
    lines.extend(
        [
            "ORIGIN",
            "        1 " + "a" * 60,
            "       61 " + "a" * 60,
            "      121 " + "a" * 60,
            "      181 " + "a" * 60,
            "      241 " + "a" * 60,
            "      301 " + "a" * 60,
            "      361 " + "a" * 60,
            "      421 " + "a" * 60,
            "      481 " + "a" * 60,
            "      541 " + "a" * 60,
            "      601 " + "a" * 60,
            "      661 " + "a" * 60,
            "      721 " + "a" * 60,
            "      781 " + "a" * 20,
            "//",
        ]
    )
    path = tmp_path / f"{name}.gb"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    return path


def _hypotheses(report):
    return (
        tuple(
            hypothesis
            for assessment in report.scanned_replicons
            for hypothesis in assessment.hypotheses
        )
        + report.cross_record_hypotheses
    )


def _synthetic_features(
    specs: tuple[tuple[str, int, int, dict[str, list[str]]], ...],
    *,
    topology: str = "linear",
    length: int = 10_000,
):
    features = [
        GenBankFeature(
            record_id="synthetic.1",
            record_index=1,
            feature_index=index,
            type=feature_type,
            location=FeatureLocation(start - 1, end, strand=1),
            qualifiers=qualifiers,
            record_length=length,
            topology=topology,
        )
        for index, (feature_type, start, end, qualifiers) in enumerate(specs, 1)
    ]
    record = GenBankRecord(
        id="synthetic.1",
        name="synthetic",
        description="Synthetic direct mobilome contract record.",
        seq=Seq("A" * length),
        length=length,
        topology=topology,
        features=features,
    )
    document = GenBankDocument(path=Path("synthetic.gb"), records=[record])
    inventory = inventory_replicons(document)[0]
    database = load_mobilome_database()
    hits = scan_mobilome_features(features, database)
    return inventory, hits, database


def test_empty_and_malformed_inputs_fail_as_mobilome_errors(tmp_path: Path) -> None:
    empty = tmp_path / "empty.gb"
    empty.write_text("", encoding="utf-8")
    malformed = tmp_path / "malformed.gb"
    malformed.write_text("LOCUS broken\nnot a GenBank record\n", encoding="utf-8")

    with pytest.raises(MobilomeInputError):
        analyze_mobilome(empty)
    with pytest.raises(MobilomeInputError):
        analyze_mobilome(malformed)


def test_forbidden_claims_are_absent_from_all_outputs_and_catalog_wording() -> None:
    fixture_dir = Path("tests/fixtures")
    for fixture in (
        fixture_dir / "mobilome_inventory.gb",
        fixture_dir / "mobilome_evidence.gb",
        fixture_dir / "mobilome_inference.gb",
    ):
        report = analyze_mobilome(fixture)
        rendered = "\n".join(
            (
                render_text(report),
                serialize_mobilome_report(report, "json"),
                serialize_mobilome_report(report, "tsv"),
            )
        ).casefold()
    assert not any(claim in rendered for claim in FORBIDDEN_CLAIMS)

    data_dir = resources.files("genbank_parser").joinpath("data", "mobilome")
    catalog_text = "\n".join(
        yaml.safe_dump(yaml.safe_load(data_dir.joinpath(name).read_text()))
        for name in ("markers.yaml", "provenance.yaml", "inference.yaml")
    ).casefold()
    assert not any(claim in catalog_text for claim in FORBIDDEN_CLAIMS)


@pytest.mark.parametrize(
    "products",
    [
        ("pXO1-01 hypothetical protein",),
        ("pXO1-02 hypothetical protein", "pXO2-01 hypothetical protein"),
        (
            "pXO2-02 hypothetical protein",
            "pXO1-02 hypothetical protein",
            "pXO1-03 hypothetical protein",
        ),
        ("pXO1-01 hypothetical protein", "pXO1-01 hypothetical protein"),
    ],
)
def test_pxo_marker_counts_and_order_never_enable_an_aggregate(
    tmp_path: Path, products: tuple[str, ...]
) -> None:
    report = analyze_mobilome(
        _genbank(tmp_path, name="PXO_PATCH", products=products, source="plasmid")
    )

    assert sum(
        hit.marker_id == "pxo_numbered_product_annotation"
        for hit in report.scanned_replicons[0].hits
    ) == len(products)
    assert "pxo_like_annotation_pattern_candidate" not in {
        item.rule_id for item in _hypotheses(report)
    }


def test_single_and_multiple_phage_modules_never_enable_phagemid_claims(
    tmp_path: Path,
) -> None:
    products = (
        "phage integrase",
        "terminase",
        "capsid protein",
        "holin",
        "phage replication protein",
        "AimR-like regulator",
        "replication initiation protein",
    )
    report = analyze_mobilome(
        _genbank(tmp_path, name="PHAGE_PATCH", products=products, source="plasmid")
    )

    assert {hit.marker_id for hit in report.scanned_replicons[0].hits} >= {
        "phage_integration_candidate",
        "phage_packaging_candidate",
        "phage_structure_candidate",
        "phage_lysis_candidate",
        "phage_replication_candidate",
        "aimr_like_candidate",
    }
    assert "phage_module_bearing_plasmid_candidate" not in {
        item.rule_id for item in _hypotheses(report)
    }


def test_crispr_comment_isolation_and_separate_observations(tmp_path: Path) -> None:
    comment_only = analyze_mobilome(
        _genbank(tmp_path, name="CRISPR_COMMENT", comment="CRISPR is mentioned here")
    )
    assert not comment_only.scanned_replicons[0].hits

    annotated = analyze_mobilome(
        _genbank(
            tmp_path,
            name="CRISPR_TYPED",
            crispr=True,
            cas=True,
            comment="CRISPR is mentioned here",
        )
    )
    assert {hit.marker_id for hit in annotated.scanned_replicons[0].hits} >= {
        "crispr_array",
        "crispr_cas",
    }
    assert not any(item.kind == "crispr" for item in _hypotheses(annotated))


def test_pseudogenized_rep_only_is_insufficient_and_mechanism_unresolved(
    tmp_path: Path,
) -> None:
    report = analyze_mobilome(_genbank(tmp_path, name="PSEUDO_REP", pseudo_rep=True))
    hypotheses = _hypotheses(report)

    assert len(hypotheses) == 1
    assert hypotheses[0].status == "insufficient"
    assert (
        hypotheses[0].summary
        == "Replication evidence insufficient; mechanism unresolved"
    )
    assert "satellite" not in render_text(report).casefold()
    assert not any(
        term in render_text(report).casefold()
        for term in (
            "theta replication",
            "rolling-circle replication",
            "strand-displacement",
        )
    )


def test_replication_and_pxo_observations_do_not_select_a_mechanism(
    tmp_path: Path,
) -> None:
    report = analyze_mobilome(
        _genbank(
            tmp_path,
            name="MECHANISM_PATCH",
            products=(
                "replication initiation protein",
                "pXO1-01 hypothetical protein",
                "pXO1-02 hypothetical protein",
            ),
        )
    )
    text = render_text(report).casefold()
    assert "unresolved mechanism" in text
    assert not any(
        term in text for term in ("theta", "rolling-circle", "strand-displacement")
    )


@pytest.mark.parametrize("topology", ["linear", "circular"])
@pytest.mark.parametrize("gap, expected", [(1999, True), (2000, True), (2001, False)])
def test_toxin_antitoxin_gap_boundary_is_inclusive(
    topology: str, gap: int, expected: bool
) -> None:
    inventory, hits, database = _synthetic_features(
        (
            ("CDS", 100, 110, {"gene": ["toxN"], "product": ["ToxN-family toxin"]}),
            (
                "ncRNA",
                111 + gap,
                121 + gap,
                {"gene": ["toxI"], "product": ["ToxI RNA antitoxin"]},
            ),
        ),
        topology=topology,
    )
    hypotheses = infer_replicon_hypotheses(inventory, hits, database)

    assert any(item.kind == "toxin_antitoxin" for item in hypotheses) is expected


@pytest.mark.parametrize(
    "spec",
    [
        (("CDS", 100, 110, {"gene": ["toxN"], "product": ["ToxN-family toxin"]}),),
        (
            (
                "ncRNA",
                100,
                110,
                {"gene": ["toxI"], "product": ["ToxI RNA antitoxin"]},
            ),
        ),
    ],
)
def test_lone_toxin_or_antitoxin_is_observation_only(spec) -> None:
    inventory, hits, database = _synthetic_features(spec)

    assert hits
    assert not any(
        item.kind == "toxin_antitoxin"
        for item in infer_replicon_hypotheses(inventory, hits, database)
    )


def test_mobility_component_subsets_report_exact_missing_facets() -> None:
    cases = {
        "relaxase_only": (
            (("CDS", 100, 140, {"gene": ["mobA"], "product": ["relaxase"]}),),
            "mobilizable_core_annotation_candidate",
            ("mobility_orit",),
        ),
        "t4cp_only": (
            (
                (
                    "CDS",
                    100,
                    140,
                    {"gene": ["traD"], "product": ["type IV coupling protein"]},
                ),
            ),
            "helper_machinery_annotation_candidate",
            ("mobility_mpf",),
        ),
        "single_mpf": (
            (
                (
                    "CDS",
                    100,
                    140,
                    {"gene": ["traA"], "product": ["conjugal transfer protein"]},
                ),
            ),
            "helper_machinery_annotation_candidate",
            ("mobility_mpf", "mobility_t4cp"),
        ),
        "t4cp_plus_one_mpf": (
            (
                (
                    "CDS",
                    100,
                    140,
                    {"gene": ["traD"], "product": ["type IV coupling protein"]},
                ),
                (
                    "CDS",
                    200,
                    240,
                    {"gene": ["traA"], "product": ["conjugal transfer protein"]},
                ),
            ),
            "helper_machinery_annotation_candidate",
            ("mobility_mpf",),
        ),
    }
    for specs, rule_id, missing in cases.values():
        inventory, hits, database = _synthetic_features(specs)
        hypotheses = infer_replicon_hypotheses(inventory, hits, database)
        helper = next(item for item in hypotheses if item.rule_id == rule_id)
        assert helper.status == "insufficient"
        assert helper.missing_components == missing


def test_complete_mobilizable_core_and_same_record_helper_do_not_cross_pair() -> None:
    specs = (
        ("regulatory", 100, 110, {"regulatory_class": ["oriT"]}),
        ("CDS", 200, 240, {"gene": ["mobA"], "product": ["relaxase"]}),
        (
            "CDS",
            300,
            340,
            {"gene": ["traD"], "product": ["type IV coupling protein"]},
        ),
        (
            "CDS",
            400,
            440,
            {"gene": ["traA"], "product": ["conjugal transfer protein"]},
        ),
        (
            "CDS",
            500,
            540,
            {"gene": ["traB"], "product": ["type IV secretion system protein"]},
        ),
    )
    inventory, hits, database = _synthetic_features(specs)
    local = infer_replicon_hypotheses(inventory, hits, database)
    assessment = RepliconAssessment(
        inventory=inventory,
        hits=hits,
        hypotheses=local,
    )

    assert any(
        item.rule_id == "mobilizable_core_annotation_candidate"
        and item.status == "tentative"
        for item in local
    )
    assert any(
        item.rule_id == "helper_machinery_annotation_candidate"
        and item.status == "tentative"
        for item in local
    )
    assert not infer_cross_record_hypotheses((assessment,), database)


def test_unlocatable_toxin_antitoxin_pair_emits_a_record_limitation() -> None:
    inventory = inventory_replicons(
        read_genbank(Path("tests/fixtures/mobilome_evidence.gb"))
    )[0]
    features = [
        GenBankFeature(
            record_id=inventory.record_id,
            record_index=inventory.record_index,
            feature_index=101,
            type="CDS",
            location=None,
            qualifiers={"gene": ["toxN"], "product": ["ToxN-family toxin"]},
            record_length=inventory.length,
            topology=inventory.topology,
        ),
        GenBankFeature(
            record_id=inventory.record_id,
            record_index=inventory.record_index,
            feature_index=102,
            type="ncRNA",
            location=None,
            qualifiers={"gene": ["toxI"], "product": ["ToxI RNA antitoxin"]},
            record_length=inventory.length,
            topology=inventory.topology,
        ),
    ]
    hits = scan_mobilome_features(features, load_mobilome_database())

    limitations = record_inference_limitations(
        inventory, hits, load_mobilome_database()
    )

    assert any("unlocatable" in limitation for limitation in limitations)


def test_cli_matrix_covers_include_database_and_input_errors(
    tmp_path: Path, capsys
) -> None:
    fixture = Path("tests/fixtures/mobilome_inventory.gb")
    for include in ("all", "chromosome", "plasmid", "unknown"):
        assert (
            main(["mobilome", str(fixture), "--include", include, "--format", "json"])
            == 0
        )
        assert (
            json.loads(capsys.readouterr().out)["schema_version"]
            == "gbparse.mobilome.v1"
        )

    database_dir = tmp_path / "database"
    database_dir.mkdir()
    source = resources.files("genbank_parser").joinpath("data", "mobilome")
    for name in ("markers.yaml", "provenance.yaml", "inference.yaml"):
        (database_dir / name).write_bytes(source.joinpath(name).read_bytes())
    assert (
        main(
            [
                "mobilome",
                str(fixture),
                "--database-dir",
                str(database_dir),
                "--format",
                "json",
            ]
        )
        == 0
    )
    capsys.readouterr()

    (database_dir / "inference.yaml").unlink()
    with pytest.raises(SystemExit):
        main(["mobilome", str(fixture), "--database-dir", str(database_dir)])
    missing_database_error = capsys.readouterr().err
    assert "exactly" in missing_database_error
    assert "Traceback" not in missing_database_error

    with pytest.raises(SystemExit) as invalid_threshold:
        main(["mobilome", str(fixture), "--min-evidence", "4"])
    threshold_error = capsys.readouterr().err
    assert invalid_threshold.value.code == 2
    assert "invalid choice" in threshold_error
    assert "Traceback" not in threshold_error

    for input_path in (tmp_path / "empty.gb", tmp_path / "malformed.gb"):
        input_path.write_text(
            "" if input_path.name == "empty.gb" else "broken", encoding="utf-8"
        )
        with pytest.raises(SystemExit):
            main(["mobilome", str(input_path)])
        assert "Traceback" not in capsys.readouterr().err


def test_database_rejects_curated_sources_without_supports(tmp_path: Path) -> None:
    source = resources.files("genbank_parser").joinpath("data", "mobilome")
    for name in ("markers.yaml", "provenance.yaml", "inference.yaml"):
        (tmp_path / name).write_bytes(source.joinpath(name).read_bytes())
    provenance_path = tmp_path / "provenance.yaml"
    payload = yaml.safe_load(provenance_path.read_text(encoding="utf-8"))
    payload["sources"].append(
        {
            "id": "bad-curated",
            "type": "curated_database",
            "name": "Synthetic database",
            "release": "1",
            "url": "https://example.invalid",
            "taxonomic_scope": "synthetic",
            "limitations": ["test"],
        }
    )
    provenance_path.write_text(
        yaml.safe_dump(payload, sort_keys=False), encoding="utf-8"
    )
    with pytest.raises(MobilomeDatabaseError, match="supports"):
        load_mobilome_database(tmp_path)


def test_writer_rejects_aliases_cleans_interrupted_temp_and_detects_mutation(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import genbank_parser.mobilome.report as report_module

    source = tmp_path / "input.gb"
    source.write_bytes(Path("tests/fixtures/mobilome_evidence.gb").read_bytes())
    report = analyze_mobilome(source)

    symlink = tmp_path / "input-alias.gb"
    try:
        symlink.symlink_to(source)
    except OSError:
        symlink = None
    if symlink is not None:
        with pytest.raises(MobilomeOutputError, match="must not overwrite"):
            write_mobilome_report(
                report,
                source_path=source,
                database_paths=(),
                format_type="json",
                output_path=symlink,
            )

    hardlink = tmp_path / "input-hardlink.gb"
    try:
        hardlink.hardlink_to(source)
    except OSError:
        hardlink = None
    if hardlink is not None:
        with pytest.raises(MobilomeOutputError, match="must not overwrite"):
            write_mobilome_report(
                report,
                source_path=source,
                database_paths=(),
                format_type="json",
                output_path=hardlink,
            )

    output = tmp_path / "report.json"
    real_replace = report_module.os.replace

    def interrupted(*args, **kwargs):
        raise OSError("simulated interruption")

    monkeypatch.setattr(report_module.os, "replace", interrupted)
    with pytest.raises(MobilomeOutputError, match="safely"):
        write_mobilome_report(
            report,
            source_path=source,
            database_paths=(),
            format_type="json",
            output_path=output,
        )
    assert not output.exists()
    assert not list(tmp_path.glob(".report.json.*.tmp"))
    monkeypatch.setattr(report_module.os, "replace", real_replace)

    source.write_bytes(source.read_bytes() + b"\n")
    with pytest.raises(MobilomeOutputError, match="changed after analysis"):
        write_mobilome_report(
            report,
            source_path=source,
            database_paths=(),
            format_type="json",
            output_path=tmp_path / "mutated.json",
        )
