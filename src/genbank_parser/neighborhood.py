"""Structured genomic neighborhoods and deterministic serializers."""

from __future__ import annotations

import argparse
import csv
import io
import json
from dataclasses import dataclass
from pathlib import Path
from collections.abc import Sequence

from .discover import RuleMatch, load_ruleset, match_feature_rules
from .io import read_genbank
from .model import GenBankFeature
from .operons import build_operon_result
from .spatial import feature_display_bounds, resolve_target, select_cds_window

SCHEMA_VERSION = "gbparse.neighborhood.v1"
TSV_COLUMNS = (
    "record",
    "order",
    "is_target",
    "feature_index",
    "type",
    "locus_tag",
    "gene",
    "product",
    "strand",
    "genomic_start",
    "genomic_end",
    "local_start",
    "local_end",
    "partial",
    "pseudo",
    "rule_matches",
    "operon_neighbors",
)


@dataclass(frozen=True)
class NeighborhoodFeature:
    """One selected feature with monotonic local display coordinates."""

    feature: GenBankFeature
    order: int
    is_target: bool
    local_start: int
    local_end: int
    rule_matches: tuple[RuleMatch, ...] = ()
    operon_neighbors: tuple[int, ...] = ()

    def to_dict(self) -> dict[str, object]:
        feature = self.feature
        return {
            "record": feature.record_id,
            "order": self.order,
            "is_target": self.is_target,
            "feature_index": feature.feature_index,
            "type": feature.type,
            "locus_tag": feature.locus_tag,
            "gene": feature.gene,
            "product": feature.product,
            "strand": feature.strand_symbol,
            "genomic_start": feature.start,
            "genomic_end": feature.end,
            "local_start": self.local_start,
            "local_end": self.local_end,
            "partial": feature.is_partial,
            "pseudo": feature.is_pseudo,
            "rule_matches": [match.to_dict() for match in self.rule_matches],
            "operon_neighbors": list(self.operon_neighbors),
        }

    def to_tsv_dict(self) -> dict[str, object]:
        row = self.to_dict()
        row["rule_matches"] = ";".join(
            f"{match.rule_id}:{match.term}:{match.weight}"
            for match in self.rule_matches
        )
        row["operon_neighbors"] = ";".join(
            str(index) for index in self.operon_neighbors
        )
        return row


@dataclass(frozen=True)
class NeighborhoodOperonLink:
    """A tested same-strand proximity link between two displayed CDSs."""

    first_feature_index: int
    second_feature_index: int
    gap: int

    def to_dict(self) -> dict[str, int]:
        return {
            "first_feature_index": self.first_feature_index,
            "second_feature_index": self.second_feature_index,
            "gap": self.gap,
        }


@dataclass(frozen=True)
class NeighborhoodResult:
    """Renderer-independent genomic-neighborhood data contract."""

    input_path: Path
    target_query: str
    record_id: str
    record_length: int
    topology: str
    window: int
    wraps_origin: bool
    features: tuple[NeighborhoodFeature, ...]
    ruleset: str | None = None
    operon_links: tuple[NeighborhoodOperonLink, ...] = ()

    @property
    def target(self) -> NeighborhoodFeature:
        return next(feature for feature in self.features if feature.is_target)

    def to_dict(self) -> dict[str, object]:
        return {
            "schema_version": SCHEMA_VERSION,
            "input": str(self.input_path),
            "target_query": self.target_query,
            "record": {
                "id": self.record_id,
                "length": self.record_length,
                "topology": self.topology,
            },
            "window": self.window,
            "wraps_origin": self.wraps_origin,
            "ruleset": self.ruleset,
            "operon_links": [link.to_dict() for link in self.operon_links],
            "features": [feature.to_dict() for feature in self.features],
        }


def _local_coordinates(
    features: tuple[GenBankFeature, ...],
    *,
    circular: bool,
    record_length: int,
) -> tuple[tuple[int, int], ...]:
    """Map genomic coordinates onto one monotonic, linear display interval."""
    if not features:
        return ()

    unwrapped: list[tuple[int, int]] = []
    shift = 0
    previous_start: int | None = None
    for feature in features:
        raw_start, raw_end = feature_display_bounds(
            feature,
            circular=circular,
            record_length=record_length,
        )
        start = raw_start + shift
        if circular and previous_start is not None and start < previous_start:
            shift += record_length
            start = raw_start + shift
        end = raw_end + shift
        if end < start:
            end += record_length
        unwrapped.append((start, end))
        previous_start = start

    origin = unwrapped[0][0]
    return tuple((start - origin + 1, end - origin + 1) for start, end in unwrapped)


def _context_coordinates(
    feature: GenBankFeature,
    *,
    anchor_start: int,
    span_end: int,
    circular: bool,
    record_length: int,
) -> tuple[int, int] | None:
    """Map a context feature to the selected unwrapped span when it overlaps."""
    shifts = (-record_length, 0, record_length) if circular else (0,)
    for shift in shifts:
        start = feature.start + shift
        end = feature.end + shift
        if end < start:
            end += record_length
        local_start = start - anchor_start + 1
        local_end = end - anchor_start + 1
        if local_end >= 1 and local_start <= span_end:
            return max(1, local_start), min(span_end, local_end)
    return None


def build_neighborhood(
    filepath: str | Path,
    target: str,
    window: int = 5,
    *,
    include_feature_types: Sequence[str] = ("CDS",),
    ruleset: str | None = None,
    show_operons: bool = False,
    operon_gap: int = 150,
) -> NeighborhoodResult:
    """Build a neighborhood without printing or terminating the process."""
    input_path = Path(filepath)
    document = read_genbank(input_path)
    record, target_feature = resolve_target(document, target)
    selected = select_cds_window(record, target_feature, window)
    coordinates = _local_coordinates(
        selected.features,
        circular=record.topology == "circular",
        record_length=record.length,
    )
    selected_coordinates = {
        feature.feature_index: coordinate
        for feature, coordinate in zip(selected.features, coordinates)
    }
    allowed_types = {feature_type.casefold() for feature_type in include_feature_types}
    allowed_types.add("cds")
    anchor_start = selected.features[0].start
    span_end = max(end for _, end in coordinates)

    displayed: list[tuple[GenBankFeature, int, int]] = [
        (feature, *selected_coordinates[feature.feature_index])
        for feature in selected.features
    ]
    selected_indices = set(selected_coordinates)
    for feature in record.features:
        if feature.feature_index in selected_indices or feature.type.casefold() == "cds":
            continue
        if feature.type.casefold() not in allowed_types:
            continue
        context_coordinates = _context_coordinates(
            feature,
            anchor_start=anchor_start,
            span_end=span_end,
            circular=record.topology == "circular",
            record_length=record.length,
        )
        if context_coordinates is not None:
            displayed.append((feature, *context_coordinates))

    displayed.sort(key=lambda item: (item[1], item[2], item[0].feature_index))
    rules = load_ruleset(ruleset) if ruleset is not None else None
    links: tuple[NeighborhoodOperonLink, ...] = ()
    neighbor_indices: dict[int, set[int]] = {}
    if show_operons:
        operon_result = build_operon_result(
            record.features,
            max_gap=operon_gap,
            circular=record.topology == "circular",
            record_length=record.length,
        )
        links = tuple(
            NeighborhoodOperonLink(
                pair.first.feature_index,
                pair.second.feature_index,
                pair.gap,
            )
            for pair in operon_result.pairs
            if pair.first.feature_index in selected_indices
            and pair.second.feature_index in selected_indices
        )
        for link in links:
            neighbor_indices.setdefault(link.first_feature_index, set()).add(
                link.second_feature_index
            )
            neighbor_indices.setdefault(link.second_feature_index, set()).add(
                link.first_feature_index
            )

    features = tuple(
        NeighborhoodFeature(
            feature=feature,
            order=order,
            is_target=feature.feature_index == target_feature.feature_index,
            local_start=local_start,
            local_end=local_end,
            rule_matches=match_feature_rules(feature, rules) if rules is not None else (),
            operon_neighbors=tuple(sorted(neighbor_indices.get(feature.feature_index, ()))),
        )
        for order, (feature, local_start, local_end) in enumerate(displayed, 1)
    )
    return NeighborhoodResult(
        input_path=input_path,
        target_query=target,
        record_id=record.id,
        record_length=record.length,
        topology=record.topology or "linear",
        window=window,
        wraps_origin=selected.wraps_origin,
        features=features,
        ruleset=ruleset,
        operon_links=links,
    )


def serialize_neighborhood(
    result: NeighborhoodResult,
    format_type: str = "text",
) -> str:
    """Serialize a structured neighborhood as text, TSV, or JSON."""
    if format_type == "json":
        return json.dumps(result.to_dict(), indent=2, ensure_ascii=False) + "\n"
    if format_type == "tsv":
        stream = io.StringIO()
        writer = csv.DictWriter(
            stream,
            fieldnames=TSV_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(feature.to_tsv_dict() for feature in result.features)
        return stream.getvalue()
    if format_type != "text":
        raise ValueError(f"Unsupported neighborhood format: {format_type}")

    lines = [
        "=" * 90,
        f"  GENOMIC NEIGHBORHOOD AROUND {result.target_query}",
        "=" * 90,
        (
            f"  Contig       : {result.record_id} (topology: {result.topology}, "
            f"length: {result.record_length:,} bp)"
        ),
        f"  Window size  : +/- {result.window} CDSs",
        "",
        (
            f"{'Strand':6s}  {'Locus Tag':18s}  {'Gene':8s}  "
            f"{'Start':>10s}  {'End':>10s}  Product"
        ),
        "-" * 90,
    ]
    for item in result.features:
        feature = item.feature
        pointer = ">> " if item.is_target else "   "
        lines.append(
            f"{pointer}[{feature.strand_symbol}]  "
            f"{(feature.locus_tag or '-'):18s}  "
            f"{(feature.gene or '-'):8s}  "
            f"{feature.start:>10,d}  {feature.end:>10,d}  "
            f"{(feature.product or '-')[:40]}"
        )
    lines.append("=" * 90)
    return "\n".join(lines) + "\n"


def write_neighborhood(
    result: NeighborhoodResult,
    *,
    format_type: str = "text",
    output_path: str | Path | None = None,
) -> str:
    """Return serialized data and optionally write it to a safe output path."""
    rendered = serialize_neighborhood(result, format_type)
    if output_path is not None:
        output = Path(output_path)
        if output.resolve() == result.input_path.resolve():
            raise ValueError("Neighborhood output cannot overwrite the input GenBank file")
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(rendered, encoding="utf-8", newline="")
    return rendered


def extract_neighborhood(
    filepath: str | Path,
    locus_tag: str,
    window: int = 5,
) -> list[GenBankFeature]:
    """Compatibility wrapper that prints text and returns selected CDS features."""
    result = build_neighborhood(filepath, locus_tag, window)
    print(serialize_neighborhood(result, "text"), end="")
    return [item.feature for item in result.features]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="View genomic neighborhood (+/- N genes) around a target locus."
    )
    parser.add_argument("input", help="Input GenBank file")
    parser.add_argument("locus_tag", help="Target locus tag or gene name")
    parser.add_argument("window", nargs="?", type=int, default=5)
    parser.add_argument("--format", choices=["text", "tsv", "json"], default="text")
    parser.add_argument("--output", help="Data output path (default: stdout)")
    args = parser.parse_args()

    result = build_neighborhood(args.input, args.locus_tag, args.window)
    rendered = write_neighborhood(
        result, format_type=args.format, output_path=args.output
    )
    if args.output is None:
        print(rendered, end="")


if __name__ == "__main__":
    main()
