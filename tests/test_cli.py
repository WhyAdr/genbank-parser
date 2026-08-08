"""Test unified gbparse CLI and backwards-compatible scripts."""
from pathlib import Path
import subprocess
import sys

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


def test_legacy_scripts_execution(simple_cds_gbff: Path) -> None:
    scripts_dir = Path(__file__).resolve().parent.parent / "scripts"
    val_script = scripts_dir / "genbank_validate.py"

    proc = subprocess.run([sys.executable, str(val_script), str(simple_cds_gbff)], capture_output=True, text=True)
    assert proc.returncode == 0
    assert "GENBANK FEATURE TABLE -- STRUCTURAL & BIOLOGICAL REPORT" in proc.stdout
