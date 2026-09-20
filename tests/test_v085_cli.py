"""Focused v0.8.5 interoperability and CLI-contract tests."""

from __future__ import annotations

import gzip
import io
import json
from pathlib import Path

import pytest
from Bio.SeqFeature import FeatureLocation

from genbank_parser.cli import main
from genbank_parser.cli_io import OutputError, write_text
from genbank_parser.crispr import render_crispr_report
from genbank_parser.diff import diff_annotations
from genbank_parser.io import GenBankInputError, read_genbank
from genbank_parser.model import GenBankFeature
from genbank_parser.operons import (
    build_operon_result,
    find_operon_pairs,
    render_operons,
)
from genbank_parser.query import parse_query, query_features


def test_plain_gzip_stream_and_stdin_parse_to_same_features(
    simple_cds_gbff: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    compressed = tmp_path / "input.unusual"
    compressed.write_bytes(gzip.compress(simple_cds_gbff.read_bytes()))

    def identities(document):
        return [
            (feature.record_id, feature.feature_index, feature.type, feature.locus_tag)
            for feature in document.all_features
        ]

    plain = read_genbank(simple_cds_gbff)
    assert identities(read_genbank(compressed)) == identities(plain)
    assert identities(read_genbank(io.StringIO(simple_cds_gbff.read_text(encoding="utf-8")))) == identities(plain)

    monkeypatch.setattr("sys.stdin", io.StringIO(simple_cds_gbff.read_text(encoding="utf-8")))
    assert identities(read_genbank("-")) == identities(plain)

    with pytest.raises(GenBankInputError, match="truncated or invalid gzip"):
        bad = tmp_path / "bad.gbff"
        bad.write_bytes(b"\x1f\x8b\x08truncated")
        read_genbank(bad)


def test_query_grammar_is_bounded_and_multivalue_aware(simple_cds_gbff: Path) -> None:
    matches = query_features(
        simple_cds_gbff,
        'type == "CDS" and gene in ["testA", "testB"] and not pseudo',
    )
    assert [match.row.locus_tag for match in matches] == ["TEST_001", "TEST_002"]
    assert len(query_features(simple_cds_gbff, 'type == "CDS" and qualifier("gene") ~= "^test"')) == 2
    assert len(query_features(simple_cds_gbff, 'type == "CDS" and length > 29.5')) == 2
    assert parse_query("not not pseudo")
    with pytest.raises(ValueError, match="unknown query field"):
        parse_query("__import__ == \"os\"")
    with pytest.raises(ValueError, match="nesting"):
        parse_query("(" * 40 + "type == \"CDS\"" + ")" * 40)


def test_query_jsonl_is_data_only_and_output_is_atomic(
    simple_cds_gbff: Path, tmp_path: Path, capsys
) -> None:
    output = tmp_path / "query.jsonl"
    assert main(
        [
            "query",
            str(simple_cds_gbff),
            "--where",
            'type == "CDS"',
            "--select",
            "record,locus_tag,gene,ko",
            "--format",
            "jsonl",
            "--output",
            str(output),
        ]
    ) == 0
    assert capsys.readouterr().out == ""
    rows = [json.loads(line) for line in output.read_text(encoding="utf-8").splitlines()]
    assert [row["locus_tag"] for row in rows] == ["TEST_001", "TEST_002"]
    assert main(
        [
            "query",
            str(simple_cds_gbff),
            "--where",
            'type == "CDS"',
            "--format",
            "tsv",
            "--output",
            str(output),
        ]
    ) == 4
    assert len(rows) == 2

    write_text("old\n", output, force=True)
    with pytest.raises(OutputError):
        write_text("new\n", output)
    assert output.read_text(encoding="utf-8") == "old\n"
    blocked_parent = tmp_path / "blocked-parent"
    blocked_parent.write_text("not a directory\n", encoding="utf-8")
    with pytest.raises(OutputError):
        write_text("data\n", blocked_parent / "out.txt")


def test_validate_gate_writes_report_before_exit(
    duplicate_locus_gbff: Path, tmp_path: Path
) -> None:
    output = tmp_path / "validation.json"
    assert main(
        [
            "validate",
            str(duplicate_locus_gbff),
            "--format",
            "json",
            "--output",
            str(output),
            "--fail-on",
            "warning",
        ]
    ) == 1
    payload = json.loads(output.read_text(encoding="utf-8"))
    assert any(item["code"] == "DUPLICATE_LOCUS_TAG" for item in payload)


def test_unified_export_formats_and_ncbi_ids(
    simple_cds_gbff: Path, tmp_path: Path
) -> None:
    jsonl = tmp_path / "features.jsonl"
    assert main(["export", str(simple_cds_gbff), "--format", "jsonl", "--output", str(jsonl)]) == 0
    assert all(json.loads(line)["schema_version"] == "gbparse.feature.v1" for line in jsonl.read_text(encoding="utf-8").splitlines())

    bed = tmp_path / "features.bed"
    assert main(["export", str(simple_cds_gbff), "--format", "bed12", "--output", str(bed)]) == 0
    for line in bed.read_text(encoding="utf-8").splitlines():
        fields = line.split("\t")
        assert len(fields) == 12
        assert 0 <= int(fields[1]) < int(fields[2]) <= 100

    output_dir = tmp_path / "ncbi"
    assert main(["export", str(simple_cds_gbff), "--format", "ncbi-table", "--output-dir", str(output_dir)]) == 0
    fasta_ids = {
        line[1:].split()[0]
        for line in (output_dir / "simple_cds.fsa").read_text(encoding="utf-8").splitlines()
        if line.startswith(">")
    }
    table_ids = {
        line.removeprefix(">Feature ")
        for line in (output_dir / "simple_cds.tbl").read_text(encoding="utf-8").splitlines()
        if line.startswith(">Feature ")
    }
    assert fasta_ids == table_ids
    assert "table2asn" in (output_dir / "simple_cds.export.json").read_text(encoding="utf-8")


def test_legacy_sequence_region_and_codon_outputs_are_data_only_and_atomic(
    simple_cds_gbff: Path, tmp_path: Path, capsys
) -> None:
    fna = tmp_path / "sequence.fna"
    ffn = tmp_path / "sequence.ffn"
    assert main(
        [
            "sequence",
            str(simple_cds_gbff),
            "--fna",
            str(fna),
            "--ffn",
            str(ffn),
        ]
    ) == 0
    assert capsys.readouterr().out == ""
    assert fna.read_text(encoding="utf-8").startswith(">CONTIG_1")
    assert ">TEST_002" in ffn.read_text(encoding="utf-8")
    fna.write_text("keep\n", encoding="utf-8")
    assert main(
        [
            "sequence",
            str(simple_cds_gbff),
            "--fna",
            str(fna),
            "--ffn",
            str(tmp_path / "sequence-2.ffn"),
        ]
    ) == 4
    assert fna.read_text(encoding="utf-8") == "keep\n"

    region = tmp_path / "region.gb"
    assert main(
        [
            "region",
            str(simple_cds_gbff),
            "--record",
            "CONTIG_1",
            "--start",
            "1",
            "--end",
            "30",
            "--output",
            str(region),
        ]
    ) == 0
    assert capsys.readouterr().out == ""
    assert read_genbank(region).records[0].length == 30

    codon = tmp_path / "codon.tsv"
    assert main(
        [
            "codon",
            str(simple_cds_gbff),
            "--min-len",
            "1",
            "--output",
            str(codon),
        ]
    ) == 0
    assert capsys.readouterr().out == ""
    assert codon.read_text(encoding="utf-8").startswith("TranslationTable\t")


def test_compare_evidence_reconciles_with_matrix(simple_cds_gbff: Path, tmp_path: Path) -> None:
    matrix = tmp_path / "matrix.tsv"
    evidence = tmp_path / "evidence.tsv"
    assert main(
        [
            "compare",
            str(simple_cds_gbff),
            "--targets",
            "testA,K00844",
            "--mode",
            "presence",
            "--format",
            "tsv",
            "--output",
            str(matrix),
            "--evidence-output",
            str(evidence),
        ]
    ) == 0
    rows = matrix.read_text(encoding="utf-8").splitlines()
    assert rows[0].endswith("testA\tK00844")
    assert rows[1].endswith("1\t1")
    assert len(evidence.read_text(encoding="utf-8").splitlines()) == 3


def test_compare_rejects_shared_matrix_and_evidence_output(
    simple_cds_gbff: Path, tmp_path: Path
) -> None:
    output = tmp_path / "same.tsv"
    assert main(
        [
            "compare",
            str(simple_cds_gbff),
            "--targets",
            "testA",
            "--output",
            str(output),
            "--evidence-output",
            str(output),
        ]
    ) == 4
    assert not output.exists()


def test_operons_json_contract(simple_cds_gbff: Path, capsys) -> None:
    assert main(["operons", str(simple_cds_gbff), "--format", "json"]) == 0
    report = json.loads(capsys.readouterr().out)
    assert report["schema_version"] == "gbparse.operons.v1"
    assert "candidate" in report["interpretation"]


def test_reverse_strand_operon_gap_and_biological_order() -> None:
    def feature(index: int, start: int, end: int) -> GenBankFeature:
        return GenBankFeature(
            record_id="contig",
            record_index=1,
            feature_index=index,
            type="CDS",
            location=FeatureLocation(start - 1, end, strand=-1),
            qualifiers={"locus_tag": [f"REV_{index}"]},
            record_length=100,
        )

    features = [feature(1, 1, 10), feature(2, 30, 40), feature(3, 60, 70)]
    pairs = find_operon_pairs(features, max_gap=25)
    assert [(pair[0].locus_tag, pair[1].locus_tag, pair[2]) for pair in pairs] == [
        ("REV_2", "REV_1", 19),
        ("REV_3", "REV_2", 19),
    ]
    result = build_operon_result(features, max_gap=25)
    assert [item.locus_tag for item in result.clusters[0].features] == [
        "REV_3",
        "REV_2",
        "REV_1",
    ]
    tsv = render_operons(
        {
            "records": [
                {
                    "record": "contig",
                    "pairs": [
                        {
                            "first": {"feature_index": 2, "locus_tag": "REV_2", "gene": None, "start": 30, "end": 40, "strand": "-"},
                            "second": {"feature_index": 1, "locus_tag": "REV_1", "gene": None, "start": 1, "end": 10, "strand": "-"},
                            "gap": 19,
                        }
                    ],
                    "clusters": [],
                }
            ],
            "pair_count": 1,
            "cluster_count": 0,
        },
        "tsv",
    )
    assert "\t1\t40\t19\ttrue" in tsv


def test_structured_tsv_uses_empty_fields_for_missing_crispr_qualifiers() -> None:
    rendered = render_crispr_report(
        {
            "arrays": [
                {
                    "record": "contig",
                    "type": "repeat_region",
                    "locus_tag": None,
                    "gene": None,
                    "product": None,
                    "start": 1,
                    "end": 10,
                    "strand": ".",
                }
            ],
            "cas_genes": [],
            "colocalized": [],
            "colocalized_count": 0,
        },
        "tsv",
    )
    assert "None" not in rendered


def test_diff_paf_preserves_mapped_coordinates(
    simple_cds_gbff: Path, tmp_path: Path, capsys
) -> None:
    renamed = tmp_path / "renamed.gb"
    renamed.write_text(
        simple_cds_gbff.read_text(encoding="utf-8").replace("CONTIG_1", "RENAMED"),
        encoding="utf-8",
    )
    paf = tmp_path / "old-to-new.paf"
    paf.write_text(
        "CONTIG_1\t100\t0\t100\t+\tRENAMED\t100\t0\t100\t100\t100\t60\n",
        encoding="utf-8",
    )

    result = diff_annotations(
        simple_cds_gbff,
        renamed,
        format_type="json",
        coordinate_map=paf,
    )
    capsys.readouterr()
    mapped = next(change for change in result["changes"] if change["match"] == "locus_tag")
    assert mapped["mapped_old"] == {
        "record": "RENAMED",
        "start": 1,
        "end": 30,
        "strand": "+",
        "strand_value": 1,
    }


def test_diff_jsonl_summary_only_retains_explicit_metadata(
    simple_cds_gbff: Path, capsys
) -> None:
    result = diff_annotations(
        simple_cds_gbff,
        simple_cds_gbff,
        format_type="jsonl",
        summary_only=True,
    )
    line = capsys.readouterr().out.strip()
    payload = json.loads(line)
    assert payload["record_type"] == "summary"
    assert payload["summary_only"] is True
    assert payload["result_complete"] is True
    assert payload["change_count"] == result["change_count"]
