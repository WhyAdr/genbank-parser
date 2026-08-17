"""Test unified gbparse CLI subcommands."""
from pathlib import Path

import pytest

from genbank_parser.cli import main


def test_cli_subcommands_dispatch(simple_cds_gbff: Path, tmp_path: Path) -> None:
    # validate
    assert main(["validate", str(simple_cds_gbff)]) == 0

    # summary
    assert main(["summary", str(simple_cds_gbff)]) == 0

    # search
    assert main(["search", str(simple_cds_gbff), "--gene", "testA"]) == 0

    # locus
    assert main(["locus", str(simple_cds_gbff), "TEST_001"]) == 0

    # neighborhood
    assert main(["neighborhood", str(simple_cds_gbff), "TEST_001", "2"]) == 0

    # fasta
    fa_out = tmp_path / "out_cli.faa"
    assert main(["fasta", str(simple_cds_gbff), str(fa_out)]) == 0
    assert fa_out.exists()

    # extract
    tsv_out = tmp_path / "out_cli.tsv"
    assert main(["extract", str(simple_cds_gbff), str(tsv_out)]) == 0
    assert tsv_out.exists()

    # codon
    assert main(["codon", str(simple_cds_gbff), "--min-len", "5"]) == 0

    # functional
    assert main(["functional", str(simple_cds_gbff)]) == 0

    # discover
    assert main(["discover", str(simple_cds_gbff)]) == 0

    # gff
    gff_out = tmp_path / "out_cli.gff3"
    assert main(["gff", str(simple_cds_gbff), str(gff_out)]) == 0
    assert gff_out.exists()


def test_neighborhood_cli_machine_outputs(simple_cds_gbff: Path, tmp_path: Path, capsys) -> None:
    assert main([
        "neighborhood",
        str(simple_cds_gbff),
        "TEST_001",
        "1",
        "--format",
        "json",
    ]) == 0
    assert capsys.readouterr().out.startswith('{\n  "schema_version"')

    output = tmp_path / "neighborhood.tsv"
    assert main([
        "neighborhood",
        str(simple_cds_gbff),
        "TEST_001",
        "1",
        "--format",
        "tsv",
        "--output",
        str(output),
    ]) == 0
    assert capsys.readouterr().out == ""
    assert output.read_text(encoding="utf-8").startswith("record\torder\t")

    with pytest.raises(SystemExit):
        main([
            "neighborhood",
            str(simple_cds_gbff),
            "TEST_001",
            "--output",
            str(output),
            "--viz-output",
            str(output),
        ])


def test_meor_cli_formats_output_and_validation(tmp_path: Path, capsys) -> None:
    fixture = Path("tests/fixtures/meor_parity.gb")
    assert main(["meor", str(fixture)]) == 0
    assert "MEOR & BIOSURFACTANT DISCOVERY REPORT" in capsys.readouterr().out

    json_out = tmp_path / "meor.json"
    assert main([
        "meor",
        str(fixture),
        "--min-weight",
        "2",
        "--format",
        "json",
        "--output",
        str(json_out),
    ]) == 0
    assert '"schema_version": "gbparse.meor.v1"' in json_out.read_text(encoding="utf-8")

    tsv_out = tmp_path / "meor.tsv"
    assert main([
        "meor",
        str(fixture),
        "--format",
        "tsv",
        "--max-gap",
        "300",
        "--window-size",
        "1000",
        "--output",
        str(tsv_out),
    ]) == 0
    assert tsv_out.read_text(encoding="utf-8").startswith("contig\tlocus_tag")

    for args in (
        ["--min-weight", "0"],
        ["--min-weight", "4"],
        ["--max-gap", "-1"],
        ["--window-size", "0"],
    ):
        with pytest.raises(SystemExit):
            main(["meor", str(fixture), *args])

    with pytest.raises(SystemExit):
        main(["meor", str(fixture), "--output", str(fixture)])
