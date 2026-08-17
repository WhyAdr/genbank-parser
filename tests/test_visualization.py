"""Optional dna_features_viewer neighborhood rendering."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

from genbank_parser.neighborhood import build_neighborhood
from genbank_parser.visualization.neighborhood import (
    VisualizationDependencyError,
    default_visualization_path,
    feature_label,
    render_neighborhood,
)


def test_core_neighborhood_does_not_import_plotting(simple_cds_gbff: Path) -> None:
    build_neighborhood(simple_cds_gbff, "TEST_001", 1)
    assert "dna_features_viewer" not in sys.modules
    assert "matplotlib.pyplot" not in sys.modules


def test_missing_visualization_dependency_is_actionable(
    simple_cds_gbff: Path, tmp_path: Path, monkeypatch
) -> None:
    result = build_neighborhood(simple_cds_gbff, "TEST_001", 1)
    monkeypatch.setitem(sys.modules, "dna_features_viewer", None)
    with pytest.raises(VisualizationDependencyError, match="genbank-parser\\[viz\\]"):
        render_neighborhood(result, tmp_path / "missing.svg")


def test_visualization_names_and_label_priority(simple_cds_gbff: Path) -> None:
    result = build_neighborhood(simple_cds_gbff, "TEST_001", 1)
    assert default_visualization_path(result).name == (
        "simple_cds-TEST_001-neighborhood.svg"
    )
    assert feature_label(result.target) == "testA"


@pytest.mark.parametrize("suffix", [".svg", ".png", ".pdf"])
def test_render_neighborhood_smoke(
    simple_cds_gbff: Path, tmp_path: Path, suffix: str
) -> None:
    pytest.importorskip("dna_features_viewer")
    result = build_neighborhood(simple_cds_gbff, "TEST_001", 1)
    output = render_neighborhood(result, tmp_path / f"neighborhood{suffix}")
    assert output.is_file()
    assert output.stat().st_size > 100
    if suffix == ".svg":
        assert "<svg" in output.read_text(encoding="utf-8")


def test_cli_visualization_smoke(simple_cds_gbff: Path, tmp_path: Path) -> None:
    pytest.importorskip("dna_features_viewer")
    from genbank_parser.cli import main

    output = tmp_path / "cli-neighborhood.svg"
    assert main([
        "neighborhood",
        str(simple_cds_gbff),
        "TEST_001",
        "1",
        "--visualize",
        "--viz-output",
        str(output),
        "--output",
        str(tmp_path / "cli-neighborhood.txt"),
    ]) == 0
    assert "<svg" in output.read_text(encoding="utf-8")


def test_ruleset_context_and_operon_render(
    context_neighborhood_gbff: Path, tmp_path: Path
) -> None:
    pytest.importorskip("dna_features_viewer")
    result = build_neighborhood(
        context_neighborhood_gbff,
        "CTX_001",
        1,
        include_feature_types=("CDS", "tRNA"),
        ruleset="mobilome",
        show_operons=True,
        operon_gap=60,
    )
    output = render_neighborhood(
        result,
        tmp_path / "context.svg",
        color_mode="ruleset",
    )
    svg = output.read_text(encoding="utf-8")
    assert "annotation-rule match" in svg
    assert "Same-strand proximity link" in svg


def test_cli_ruleset_context_overlay(
    context_neighborhood_gbff: Path, tmp_path: Path
) -> None:
    pytest.importorskip("dna_features_viewer")
    from genbank_parser.cli import main

    data_output = tmp_path / "context.json"
    figure_output = tmp_path / "context.svg"
    assert main([
        "neighborhood",
        str(context_neighborhood_gbff),
        "CTX_001",
        "1",
        "--format",
        "json",
        "--output",
        str(data_output),
        "--visualize",
        "--viz-output",
        str(figure_output),
        "--include-feature-types",
        "CDS,tRNA",
        "--color-by",
        "ruleset",
        "--ruleset",
        "mobilome",
        "--show-operons",
        "--operon-gap",
        "60",
    ]) == 0
    assert '"ruleset": "mobilome"' in data_output.read_text(encoding="utf-8")
    assert "annotation-rule match" in figure_output.read_text(encoding="utf-8")
