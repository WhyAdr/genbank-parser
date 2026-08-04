#!/usr/bin/env python3
"""Build comparable genome and annotation summaries from Bakta GenBank output.

GenBank files (``.gbff``, ``.gbk``, or ``.gb``) supply the feature-level
annotations. A same-named Bakta ``.txt`` summary, when present, supplies the
published headline counts and sequence statistics.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Iterable

from genbank_parser import extract_xrefs, get_qual, parse_features
from genbank_sequence import parse_sequences


GENBANK_SUFFIXES = {".gbff", ".gbk", ".gb"}
FIELDNAMES = [
    "Sample / Isolate",
    "Genome Size (Mbp)",
    "Total Gene/Feature Count",
    "G + C content (%)",
    "CDS",
    "rRNA",
    "tRNA / tmRNA",
    "Genes assigned to UniRef",
    "Genes assigned to RefSeq",
    "Genes assigned to COGs",
    "Genes assigned to KEGG/KO",
    "Genes without function prediction",
    "Pseudogenes",
    "Regulatory ncRNAs",
]


def discover_genbank_files(inputs: Iterable[Path]) -> list[Path]:
    """Return unique compatible GenBank files from files and directories."""
    found: set[Path] = set()
    for input_path in inputs:
        if input_path.is_file() and input_path.suffix.lower() in GENBANK_SUFFIXES:
            found.add(input_path.resolve())
        elif input_path.is_dir():
            found.update(
                path.resolve()
                for path in input_path.rglob("*")
                if path.is_file() and path.suffix.lower() in GENBANK_SUFFIXES
            )
        else:
            print(f"Warning: skipping unavailable or unsupported input: {input_path}", file=sys.stderr)
    return sorted(found, key=lambda path: (path.stem.casefold(), str(path).casefold()))


def parse_bakta_summary(gbff_path: Path) -> dict[str, int | float]:
    """Read stable headline values from a same-named Bakta text report."""
    txt_path = gbff_path.with_suffix(".txt")
    if not txt_path.exists():
        return {}

    key_map = {
        "size": "genome_size_bp", "length": "genome_size_bp", "gc": "gc_content",
        "gc content": "gc_content", "cdss": "cds", "trnas": "trna",
        "tmrnas": "tmrna", "rrnas": "rrna", "ncrnas": "ncrna",
        "ncrna regions": "regulatory", "crispr arrays": "crispr", "pseudogenes": "pseudogene",
        "sorfs": "sorf", "orics": "oric", "orivs": "oriv", "orits": "orit",
    }
    summary: dict[str, int | float] = {}
    line_re = re.compile(r"^\s*([^:]+):\s*([\d,.]+)\s*(?:bp|%)?\s*$", re.IGNORECASE)
    with txt_path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            match = line_re.match(line)
            if not match:
                continue
            raw_key = re.sub(r"\s+", " ", match.group(1).strip().casefold())
            key = key_map.get(raw_key)
            if not key:
                continue
            raw_value = match.group(2).replace(",", "")
            summary[key] = float(raw_value) if key == "gc_content" else int(float(raw_value))
    return summary


def _feature_identity(feature: dict) -> tuple[str, str, int, int]:
    """Give paired gene/CDS pseudogene records a stable shared identity."""
    tag = get_qual(feature, "locus_tag")
    if tag:
        return ("tag", tag, 0, 0)
    return (feature["contig"], "", feature["start"], feature["end"])


def _is_pseudogene(feature: dict) -> bool:
    qualifiers = feature["qualifiers"]
    return feature["type"].casefold() == "pseudogene" or bool(
        qualifiers.get("pseudo") or qualifiers.get("pseudogene")
    )


def _xref_values(feature: dict) -> list[str]:
    qualifiers = feature["qualifiers"]
    return qualifiers.get("db_xref", []) + qualifiers.get("note", [])


def _has_prefixed_xref(feature: dict, prefix: str) -> bool:
    return any(value.casefold().startswith(prefix.casefold()) for value in _xref_values(feature))


def _genome_stats(
    gbff_path: Path, features: list[dict], summary: dict[str, int | float]
) -> tuple[int | None, float | None]:
    """Compute genome size and GC from sequence, with feature/text fallbacks."""
    sequences = parse_sequences(gbff_path)
    if sequences:
        total_bp = sum(len(sequence) for sequence in sequences.values())
        total_acgt = sum(
            sequence.count(base)
            for sequence in sequences.values()
            for base in "ACGT"
        )
        total_gc = sum(
            sequence.count(base)
            for sequence in sequences.values()
            for base in "GC"
        )
        genome_size = int(summary.get("genome_size_bp", total_bp))
        gc_content = summary.get("gc_content", 100 * total_gc / total_acgt if total_acgt else None)
        return genome_size, float(gc_content) if gc_content is not None else None

    source_spans = [
        feature["end"] - feature["start"] + 1
        for feature in features
        if feature["type"].casefold() == "source"
    ]
    source_size = sum(source_spans) if source_spans else None
    return (
        int(summary.get("genome_size_bp", source_size)) if summary.get("genome_size_bp", source_size) is not None else None,
        float(summary["gc_content"]) if "gc_content" in summary else None,
    )


def _summary_count(summary: dict[str, int | float], key: str, fallback: int) -> int:
    """Prefer an integer Bakta headline count when it is available."""
    return int(summary[key]) if key in summary else fallback


def summarise_file(gbff_path: Path) -> dict[str, object]:
    """Create one output row from a Bakta GenBank annotation file."""
    features = parse_features(gbff_path)
    if not features:
        raise ValueError("no GenBank features were parsed")

    summary = parse_bakta_summary(gbff_path)
    counts = Counter(feature["type"].casefold() for feature in features)
    pseudogene_ids = {
        _feature_identity(feature)
        for feature in features
        if _is_pseudogene(feature)
    }
    non_gene_ids = {
        _feature_identity(feature)
        for feature in features
        if feature["type"].casefold() not in {"source", "gene"}
    }
    pseudogene_only_genes = sum(identity not in non_gene_ids for identity in pseudogene_ids)
    gbff_total_features = sum(
        1
        for feature in features
        if feature["type"].casefold() not in {"source", "gene"}
    ) + pseudogene_only_genes

    cdss = [feature for feature in features if feature["type"].casefold() == "cds"]
    genome_size_bp, gc_content = _genome_stats(gbff_path, features, summary)
    headline_counts = {
        "cds": _summary_count(summary, "cds", len(cdss)),
        "trna": _summary_count(summary, "trna", counts["trna"]),
        "tmrna": _summary_count(summary, "tmrna", counts["tmrna"]),
        "rrna": _summary_count(summary, "rrna", counts["rrna"]),
        "ncrna": _summary_count(summary, "ncrna", counts["ncrna"]),
        "regulatory": _summary_count(summary, "regulatory", counts["regulatory"]),
        "crispr": _summary_count(summary, "crispr", 0),
        "pseudogene": _summary_count(summary, "pseudogene", len(pseudogene_ids)),
        "sorf": _summary_count(summary, "sorf", 0),
        "oric": _summary_count(summary, "oric", 0),
        "oriv": _summary_count(summary, "oriv", 0),
        "orit": _summary_count(summary, "orit", 0),
    }
    total_features = (
        sum(headline_counts.values()) if "cds" in summary else gbff_total_features
    )

    return {
        "Sample / Isolate": gbff_path.stem,
        "Genome Size (Mbp)": _format_number(
            genome_size_bp / 1_000_000 if genome_size_bp is not None else None, 3
        ),
        "Total Gene/Feature Count": total_features,
        "G + C content (%)": _format_number(gc_content, 2),
        "CDS": headline_counts["cds"],
        "rRNA": headline_counts["rrna"],
        "tRNA / tmRNA": f"{headline_counts['trna']} / {headline_counts['tmrna']}",
        "Genes assigned to UniRef": sum(_has_prefixed_xref(cds, "uniref") for cds in cdss),
        "Genes assigned to RefSeq": sum(_has_prefixed_xref(cds, "refseq") for cds in cdss),
        "Genes assigned to COGs": sum(bool(extract_xrefs(cds)["cog_ids"]) for cds in cdss),
        "Genes assigned to KEGG/KO": sum(
            bool(extract_xrefs(cds)["kegg_kos"]) or _has_prefixed_xref(cds, "ko:")
            for cds in cdss
        ),
        "Genes without function prediction": sum(
            get_qual(cds, "product").strip().casefold() == "hypothetical protein"
            for cds in cdss
        ),
        "Pseudogenes": headline_counts["pseudogene"],
        "Regulatory ncRNAs": headline_counts["regulatory"],
    }


def _format_number(value: float | None, places: int) -> str:
    return "" if value is None else f"{value:.{places}f}"


def _markdown_escape(value: object) -> str:
    return str(value).replace("\\", "\\\\").replace("|", "\\|").replace("\n", "<br>")


def write_outputs(rows: list[dict[str, object]], csv_path: Path, tsv_path: Path, md_path: Path) -> None:
    """Write the same report in CSV, TSV, and Markdown forms."""
    for path, delimiter in ((csv_path, ","), (tsv_path, "\t")):
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, delimiter=delimiter)
            writer.writeheader()
            writer.writerows(rows)

    md_path.parent.mkdir(parents=True, exist_ok=True)
    markdown_lines = ["# Bakta Isolates Summary Table", ""]
    markdown_lines.append("| " + " | ".join(FIELDNAMES) + " |")
    markdown_lines.append("| " + " | ".join("---" for _ in FIELDNAMES) + " |")
    markdown_lines.extend(
        "| " + " | ".join(_markdown_escape(row[column]) for column in FIELDNAMES) + " |"
        for row in rows
    )
    md_path.write_text("\n".join(markdown_lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Summarise Bakta GenBank annotations into CSV, TSV, and Markdown tables."
    )
    parser.add_argument(
        "-i", "--input", "--input-dir", nargs="+", required=True, type=Path,
        help="One or more GenBank files or directories to search recursively.",
    )
    parser.add_argument("--output-csv", type=Path, default=Path("bakta_genome_summary.csv"))
    parser.add_argument("--output-tsv", type=Path, default=Path("bakta_genome_summary.tsv"))
    parser.add_argument("--output-md", type=Path, default=Path("bakta_genome_summary.md"))
    args = parser.parse_args()

    files = discover_genbank_files(args.input)
    if not files:
        parser.error("no .gbff, .gbk, or .gb files were found")

    rows: list[dict[str, object]] = []
    for gbff_path in files:
        try:
            rows.append(summarise_file(gbff_path))
        except (OSError, ValueError) as error:
            print(f"Warning: unable to summarise {gbff_path}: {error}", file=sys.stderr)

    if not rows:
        print("Error: no summary records could be created.", file=sys.stderr)
        return 1

    rows.sort(key=lambda row: str(row["Sample / Isolate"]).casefold())
    write_outputs(rows, args.output_csv, args.output_tsv, args.output_md)
    print(f"Summarised {len(rows)} isolate(s).")
    print(f"CSV: {args.output_csv}")
    print(f"TSV: {args.output_tsv}")
    print(f"Markdown: {args.output_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
