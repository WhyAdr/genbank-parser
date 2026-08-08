"""Adversarial regressions for the 0.2.1 correctness patch."""

from pathlib import Path

from Bio.SeqFeature import FeatureLocation

from genbank_parser.codon import analyze_codon_usage
from genbank_parser.compare import _word_in
from genbank_parser.crispr import detect_crispr, interval_distance
from genbank_parser.diff import diff_annotations
from genbank_parser.discover import discover_clusters
from genbank_parser.gff import convert_to_gff3
from genbank_parser.io import read_genbank
from genbank_parser.model import GenBankFeature, GenBankRecord
from genbank_parser.neighborhood import extract_neighborhood
from genbank_parser.phylo import extract_phylogenomic_markers
from genbank_parser.region import extract_region
from genbank_parser.validate import validate


def test_find_locus_prefers_cds_and_exposes_all_matches(simple_cds_gbff: Path) -> None:
    doc = read_genbank(simple_cds_gbff)
    matches = doc.records[0].find_locus_features("TEST_002")
    assert [feature.type for feature in matches] == ["gene", "CDS"]
    assert doc.records[0].find_locus("TEST_002").type == "CDS"
    assert doc.find_cds("TEST_002")[1].type == "CDS"


def test_neighborhood_does_not_fall_back_to_first_cds(simple_cds_gbff: Path) -> None:
    neighborhood = extract_neighborhood(simple_cds_gbff, "TEST_002", window=0)
    assert [feature.locus_tag for feature in neighborhood] == ["TEST_002"]


def test_region_is_local_and_preserves_joined_locations(
    compound_joined_gbff: Path,
) -> None:
    forward = extract_region(compound_joined_gbff, locus_tag="JOIN_FWD_001")
    reverse = extract_region(compound_joined_gbff, locus_tag="JOIN_REV_001")
    forward_cds = next(feature for feature in forward.features if feature.type == "CDS")
    reverse_cds = next(feature for feature in reverse.features if feature.type == "CDS")
    assert forward_cds.location.__class__.__name__ == "CompoundLocation"
    assert reverse_cds.location.__class__.__name__ == "CompoundLocation"
    assert all(
        int(feature.location.end) <= len(forward.seq) for feature in forward.features
    )
    assert all(
        int(feature.location.end) <= len(reverse.seq) for feature in reverse.features
    )
    assert forward.annotations["parent_record"] == "CONTIG_JOIN"


def test_circular_region_wraps_origin(multi_record_circular_gbff: Path) -> None:
    region = extract_region(
        multi_record_circular_gbff,
        record_id="CHR_CIRCULAR",
        start=90,
        end=10,
    )
    assert len(region.seq) == 21
    assert all(
        int(feature.location.end) <= len(region.seq) for feature in region.features
    )


def test_legacy_dict_preserves_unknown_strand_states() -> None:
    for strand, symbol in ((0, "?"), (None, ".")):
        feature = GenBankFeature(
            record_id="TEST",
            record_index=1,
            feature_index=1,
            type="misc_feature",
            location=FeatureLocation(0, 5, strand=strand),
        )
        assert feature["strand"] == symbol
        assert feature.to_dict()["strand"] == symbol


def test_coding_density_uses_nonredundant_union() -> None:
    features = [
        GenBankFeature("R", 1, 1, "CDS", FeatureLocation(0, 10, strand=1)),
        GenBankFeature("R", 1, 2, "CDS", FeatureLocation(5, 15, strand=1)),
    ]
    record = GenBankRecord("R", "R", "", "A" * 20, 20, features=features)
    assert record.cds_feature_bp_sum == 20
    assert record.nonredundant_coding_bp == 15
    assert record.coding_density == 75.0


def test_gff_ids_are_unique_and_phase_respects_codon_start(
    simple_cds_gbff: Path, special_cds_gbff: Path, tmp_path: Path
) -> None:
    simple_gff = convert_to_gff3(simple_cds_gbff)
    ids = [
        line.split("ID=", 1)[1].split(";", 1)[0]
        for line in simple_gff.splitlines()
        if "\t" in line
    ]
    assert len(ids) == len(set(ids))
    assert simple_gff.count("ID=gene:TEST_001") == 1

    special_gff = convert_to_gff3(special_cds_gbff)
    codon_start_lines = [
        line for line in special_gff.splitlines() if "CODON_START_001" in line
    ]
    cds_line = next(line for line in codon_start_lines if "\tCDS\t" in line)
    assert cds_line.split("\t")[7] == "1"

    table3 = tmp_path / "codon-start-3.gb"
    table3.write_text(
        special_cds_gbff.read_text(encoding="utf-8").replace(
            "/codon_start=2", "/codon_start=3", 1
        ),
        encoding="utf-8",
    )
    table3_line = next(
        line
        for line in convert_to_gff3(table3).splitlines()
        if "CODON_START_001" in line and "\tCDS\t" in line
    )
    assert table3_line.split("\t")[7] == "2"


def test_validator_reports_total_multirecord_length(
    multi_record_circular_gbff: Path, capsys
) -> None:
    validate(multi_record_circular_gbff)
    assert "Total genome length   : 150 bp" in capsys.readouterr().out


def test_validator_reports_exception_and_unknown_translation_table(
    simple_cds_gbff: Path, tmp_path: Path
) -> None:
    exceptional = tmp_path / "exceptional.gb"
    exceptional.write_text(
        simple_cds_gbff.read_text(encoding="utf-8").replace(
            '                     /translation="MKVLWAGLIT"',
            '                     /exception="ribosomal frameshifting"\n'
            '                     /translation="MKVLWAGLIT"',
            1,
        ),
        encoding="utf-8",
    )
    exceptional_findings = validate(exceptional, json_mode=True)
    assert any(
        finding.code == "EXCEPTIONAL_TRANSLATION_ANNOTATION"
        for finding in exceptional_findings
    )

    unknown = tmp_path / "unknown-table.gb"
    unknown.write_text(
        simple_cds_gbff.read_text(encoding="utf-8").replace(
            '                     /translation="MKVLWAGLIT"',
            "                     /transl_table=99\n"
            '                     /translation="MKVLWAGLIT"',
            1,
        ),
        encoding="utf-8",
    )
    unknown_findings = validate(unknown, json_mode=True)
    assert any(
        finding.code == "UNKNOWN_TRANSLATION_TABLE" for finding in unknown_findings
    )


def test_compare_product_matching_is_not_substring_matching() -> None:
    assert _word_in("ladA", "long-chain ladA enzyme")
    assert not _word_in("ladA", "ladA-like unrelated protein")


def test_codon_usage_uses_alternate_translation_table(tmp_path: Path) -> None:
    table4 = tmp_path / "table4.gb"
    table4.write_text(
        """LOCUS       TABLE4                   9 bp    DNA     linear   BCT 08-AUG-2026
DEFINITION  Synthetic translation-table test.
ACCESSION   TABLE4
FEATURES             Location/Qualifiers
     CDS             1..9
                     /locus_tag="TABLE4_001"
                     /transl_table=4
ORIGIN
        1 atgtgataa
//
""",
        encoding="utf-8",
    )
    result = analyze_codon_usage(table4, min_len_aa=0)
    assert result["translation_tables_encountered"] == [4]
    assert result["codon_counts"]["TGA"] == 1
    assert result["terminal_stop_codons"] == 1


def test_discovery_packaged_rules_emit_tsv(
    simple_cds_gbff: Path, tmp_path: Path, capsys
) -> None:
    transposase = tmp_path / "transposase.gb"
    transposase.write_text(
        simple_cds_gbff.read_text(encoding="utf-8").replace(
            '/product="test protein A"', '/product="transposase"'
        ),
        encoding="utf-8",
    )
    result = discover_clusters(transposase, format_type="tsv")
    output = capsys.readouterr().out
    assert result["total_hits"] == 1
    assert output.splitlines()[0].startswith("record\tcluster_id\t")


def test_phylo_reports_multi_copy_candidates(
    simple_cds_gbff: Path, tmp_path: Path
) -> None:
    rec_a = tmp_path / "recA.gb"
    rec_a.write_text(
        simple_cds_gbff.read_text(encoding="utf-8")
        .replace('/gene="testA"', '/gene="recA"')
        .replace('/gene="testB"', '/gene="recA"')
        .replace('/product="test protein A"', '/product="DNA recombinase A"')
        .replace('/product="test protein B"', '/product="DNA recombinase A"'),
        encoding="utf-8",
    )
    candidates = extract_phylogenomic_markers(
        rec_a, marker_set="housekeeping", min_length=1
    )
    assert len(candidates["recA"]) == 2


def test_crispr_interval_distance_and_zero_schema(simple_cds_gbff: Path) -> None:
    assert interval_distance(10, 20, 21, 30) == 1
    assert interval_distance(10, 20, 15, 30) == 0
    result = detect_crispr(simple_cds_gbff)
    assert result["colocalized"] == []
    assert result["colocalized_count"] == 0


def test_crispr_positive_schema(tmp_path: Path) -> None:
    crispr = tmp_path / "crispr.gb"
    crispr.write_text(
        """LOCUS       CRISPR                   100 bp    DNA     linear   BCT 08-AUG-2026
DEFINITION  Synthetic CRISPR annotation test.
ACCESSION   CRISPR
FEATURES             Location/Qualifiers
     repeat_region   10..20
                     /note="CRISPR direct repeat"
     CDS             30..50
                     /locus_tag="CAS1_001"
                     /gene="cas1"
                     /product="CRISPR-associated protein Cas1"
ORIGIN
        1 aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
//
""",
        encoding="utf-8",
    )
    result = detect_crispr(crispr, window=10)
    assert len(result["arrays"]) == 1
    assert len(result["cas_genes"]) == 1
    assert result["colocalized_count"] == 1
    assert isinstance(result["colocalized"], list)


def test_diff_recognizes_boundary_shift_by_locus(
    simple_cds_gbff: Path, tmp_path: Path
) -> None:
    shifted = tmp_path / "shifted.gb"
    text = simple_cds_gbff.read_text(encoding="utf-8")
    shifted.write_text(
        text.replace("complement(51..80)", "complement(52..81)"), encoding="utf-8"
    )
    result = diff_annotations(simple_cds_gbff, shifted, format_type="json")
    assert result["boundary_shifted_cds"] == 1
    assert result["removed_cds"] == 0
    assert result["added_cds"] == 0


def test_diff_reports_ec_changes(simple_cds_gbff: Path, tmp_path: Path) -> None:
    updated = tmp_path / "ec-updated.gb"
    updated.write_text(
        simple_cds_gbff.read_text(encoding="utf-8").replace('"1.1.1.1"', '"2.2.2.2"'),
        encoding="utf-8",
    )
    result = diff_annotations(simple_cds_gbff, updated, format_type="json")
    assert result["xref_changes"] == 1
    assert result["details"]["xref_changes"]
