"""Structured genomic neighborhoods and deterministic serializers."""

from __future__ import annotations

import argparse
import csv
import io
import json
from dataclasses import dataclass
from pathlib import Path

from .io import read_genbank
from .model import GenBankFeature
from .spatial import resolve_target, select_cds_window

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
)


@dataclass(frozen=True)
class NeighborhoodFeature:
    """One selected feature with monotonic local display coordinates."""

    feature: GenBankFeature
    order: int
    is_target: bool
    local_start: int
    local_end: int

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
        start = feature.start + shift
        if circular and previous_start is not None and start < previous_start:
            shift += record_length
            start = feature.start + shift
        end = feature.end + shift
        if end < start:
            end += record_length
        unwrapped.append((start, end))
        previous_start = start

    origin = unwrapped[0][0]
    return tuple((start - origin + 1, end - origin + 1) for start, end in unwrapped)


def build_neighborhood(
    filepath: str | Path,
    target: str,
    window: int = 5,
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
    features = tuple(
        NeighborhoodFeature(
            feature=feature,
            order=order,
            is_target=feature.feature_index == target_feature.feature_index,
            local_start=local_start,
            local_end=local_end,
        )
        for order, (feature, (local_start, local_end)) in enumerate(
            zip(selected.features, coordinates), 1
        )
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
        writer.writerows(feature.to_dict() for feature in result.features)
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
