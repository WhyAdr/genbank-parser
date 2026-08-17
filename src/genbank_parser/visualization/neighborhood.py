"""Static neighborhood rendering with dna_features_viewer."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

from ..neighborhood import NeighborhoodFeature, NeighborhoodResult

SUPPORTED_SUFFIXES = {".svg", ".png", ".pdf"}
RULE_PALETTE = (
    "#4c78a8",
    "#f58518",
    "#54a24b",
    "#b279a2",
    "#e45756",
    "#72b7b2",
    "#ff9da6",
)
DEPENDENCY_MESSAGE = (
    "Neighborhood visualization requires the optional 'viz' dependencies.\n"
    "Install with: pip install 'genbank-parser[viz]'"
)


class VisualizationDependencyError(RuntimeError):
    """Raised when the optional plotting stack is unavailable."""


def _load_plotting() -> tuple[type[Any], type[Any], Any]:
    try:
        from dna_features_viewer import GraphicFeature, GraphicRecord
        import matplotlib.pyplot as pyplot
    except (ImportError, ModuleNotFoundError) as exc:
        raise VisualizationDependencyError(DEPENDENCY_MESSAGE) from exc
    return GraphicFeature, GraphicRecord, pyplot


def _shorten(text: str, limit: int = 32) -> str:
    normalized = " ".join(text.split())
    if len(normalized) <= limit:
        return normalized
    return normalized[: limit - 1].rstrip() + "…"


def feature_label(item: NeighborhoodFeature, label_mode: str = "auto") -> str | None:
    """Select a compact, deterministic label from authoritative qualifiers."""
    feature = item.feature
    if label_mode == "none":
        return None
    if label_mode == "gene":
        return feature.gene or None
    if label_mode == "locus_tag":
        return feature.locus_tag or None
    if label_mode == "product":
        return _shorten(feature.product) if feature.product else None
    if label_mode != "auto":
        raise ValueError(f"Unsupported label mode: {label_mode}")
    return (
        feature.gene
        or feature.locus_tag
        or (_shorten(feature.product) if feature.product else "")
        or feature.type
    )


def default_visualization_path(result: NeighborhoodResult) -> Path:
    """Return a deterministic SVG name in the current directory."""
    slug = re.sub(r"[^A-Za-z0-9._-]+", "-", result.target_query).strip("-._")
    if not slug:
        slug = "target"
    return Path(f"{result.input_path.stem}-{slug}-neighborhood.svg")


def render_neighborhood(
    result: NeighborhoodResult,
    output_path: str | Path | None = None,
    *,
    label_mode: str = "auto",
    color_mode: str = "default",
) -> Path:
    """Render one linear-unwrapped neighborhood and return the written path."""
    if color_mode not in {"default", "ruleset"}:
        raise ValueError(f"Unsupported color mode: {color_mode}")
    if color_mode == "ruleset" and result.ruleset is None:
        raise ValueError("Ruleset coloring requires canonical rule matches in the result")

    output = Path(output_path) if output_path is not None else default_visualization_path(result)
    if not output.suffix:
        output = output.with_suffix(".svg")
    if output.suffix.casefold() not in SUPPORTED_SUFFIXES:
        supported = ", ".join(sorted(SUPPORTED_SUFFIXES))
        raise ValueError(f"Unsupported visualization format {output.suffix!r}; use {supported}")
    if output.resolve() == result.input_path.resolve():
        raise ValueError("Visualization output cannot overwrite the input GenBank file")

    GraphicFeature, GraphicRecord, pyplot = _load_plotting()
    graphic_features = []
    rule_ids = sorted(
        {
            match.rule_id
            for item in result.features
            for match in item.rule_matches
        }
    )
    rule_colors = {
        rule_id: RULE_PALETTE[index % len(RULE_PALETTE)]
        for index, rule_id in enumerate(rule_ids)
    }
    for item in result.features:
        feature = item.feature
        if color_mode == "ruleset" and item.rule_matches:
            color = rule_colors[item.rule_matches[0].rule_id]
        elif item.is_target:
            color = "#d1495b"
        elif feature.is_pseudo:
            color = "#bdbdbd"
        else:
            color = "#9ecae1"
        graphic_features.append(
            GraphicFeature(
                start=item.local_start - 1,
                end=item.local_end,
                strand=feature.strand if feature.strand in (-1, 1) else 0,
                color=color,
                linecolor=(
                    "#d1495b"
                    if item.is_target
                    else "#7f0000"
                    if feature.is_partial
                    else "#355c7d"
                ),
                label=feature_label(item, label_mode),
            )
        )

    sequence_length = max(item.local_end for item in result.features)
    record = GraphicRecord(
        sequence_length=sequence_length,
        features=graphic_features,
        first_index=1,
        plots_indexing="biopython",
    )
    figure_width = max(8.0, min(16.0, 1.35 * len(graphic_features)))
    axis, _ = record.plot(figure_width=figure_width, with_ruler=True)
    axis.set_title(
        f"{result.record_id}: {result.target_query} neighborhood",
        loc="left",
        weight="bold",
    )
    item_by_index = {
        item.feature.feature_index: item for item in result.features
    }
    if result.operon_links:
        ymin, ymax = axis.get_ylim()
        bracket_y = ymin + 0.08 * (ymax - ymin)
        cap = 0.025 * (ymax - ymin)
        for link in result.operon_links:
            first = item_by_index[link.first_feature_index]
            second = item_by_index[link.second_feature_index]
            left = min(first.local_end, second.local_end)
            right = max(first.local_start, second.local_start)
            axis.plot(
                [left, left, right, right],
                [bracket_y + cap, bracket_y, bracket_y, bracket_y + cap],
                color="#4d4d4d",
                linewidth=1.5,
            )

    if color_mode == "ruleset" or result.operon_links:
        for rule_id in rule_ids:
            axis.plot(
                [],
                [],
                color=rule_colors[rule_id],
                linewidth=8,
                label=f"{rule_id} annotation-rule match",
            )
        if result.operon_links:
            axis.plot(
                [],
                [],
                color="#4d4d4d",
                linewidth=1.5,
                label="Same-strand proximity link",
            )
        axis.legend(loc="upper right", frameon=False, fontsize="small")
    output.parent.mkdir(parents=True, exist_ok=True)
    axis.figure.savefig(output, bbox_inches="tight", dpi=200)
    pyplot.close(axis.figure)
    return output
