"""Build comparable genome and annotation summaries from Bakta GenBank output."""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from pathlib import Path
import re
import sys
from typing import Any, Iterable

from .io import extract_xrefs, get_qual, read_genbank
from .model import GenBankFeature

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


def discover_genbank_files(inputs: Iterable[str | Path]) -> list[Path]:
    """Return unique compatible GenBank files from files and directories."""
    found: set[Path] = set()
    for input_item in inputs:
        input_path = Path(input_item)
        if input_path.is_file() and input_path.suffix.lower() in GENBANK_SUFFIXES:
            found.add(input_path.resolve())
        elif input_path.is_dir():
            found.update(
                p.resolve()
                for p in input_path.rglob("*")
                if p.is_file() and p.suffix.lower() in GENBANK_SUFFIXES
            )
        else:
            print(f"Warning: skipping unavailable or unsupported input: {input_path}", file=sys.stderr)
    return sorted(found, key=lambda p: (p.stem.casefold(), str(p).casefold()))


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


def _feature_identity(feature: GenBankFeature | dict[str, Any]) -> tuple[str, str, int, int]:
    tag = get_qual(feature, "locus_tag")
    if tag:
        return ("tag", tag, 0, 0)
    contig = getattr(feature, "record_id", None) or feature.get("contig", "")
    st = getattr(feature, "start", None) or feature.get("start", 0)
    en = getattr(feature, "end", None) or feature.get("end", 0)
    return (contig, "", st, en)


def _is_pseudogene(feature: GenBankFeature | dict[str, Any]) -> bool:
    if hasattr(feature, "is_pseudo"):
        return feature.is_pseudo
    quals = feature.get("qualifiers", {})
    ftype = feature.get("type", "")
    return ftype.casefold() == "pseudogene" or bool(quals.get("pseudo") or quals.get("pseudogene"))


def _xref_values(feature: GenBankFeature | dict[str, Any]) -> list[str]:
    if hasattr(feature, "qualifiers"):
        quals = feature.qualifiers
    else:
        quals = feature.get("qualifiers", {})
    return quals.get("db_xref", []) + quals.get("note", [])


def _has_prefixed_xref(feature: GenBankFeature | dict[str, Any], prefix: str) -> bool:
    return any(v.casefold().startswith(prefix.casefold()) for v in _xref_values(feature))


def _format_number(value: float | None, places: int) -> str:
    return "" if value is None else f"{value:.{places}f}"


def _markdown_escape(value: object) -> str:
    return str(value).replace("\\", "\\\\").replace("|", "\\|").replace("\n", "<br>")


def summarise_file(gbff_path: Path) -> dict[str, object]:
    """Create one output row from a Bakta GenBank annotation file."""
    doc = read_genbank(gbff_path)
    features = doc.all_features
    if not features:
        raise ValueError("no GenBank features were parsed")

    summary = parse_bakta_summary(gbff_path)
    counts = Counter(feature.type.casefold() for feature in features)
    pseudogene_ids = {
        _feature_identity(feature)
        for feature in features
        if _is_pseudogene(feature)
    }
    non_gene_ids = {
        _feature_identity(feature)
        for feature in features
        if feature.type.casefold() not in {"source", "gene"}
    }
    pseudogene_only_genes = sum(identity not in non_gene_ids for identity in pseudogene_ids)
    gbff_total_features = sum(
        1
        for feature in features
        if feature.type.casefold() not in {"source", "gene"}
    ) + pseudogene_only_genes

    cdss = [f for f in features if f.type.casefold() == "cds"]
    total_bp = doc.total_length
    genome_size_bp = int(summary.get("genome_size_bp", total_bp)) if total_bp > 0 or "genome_size_bp" in summary else None

    # GC content
    if "gc_content" in summary:
        gc_content: float | None = float(summary["gc_content"])
    else:
        gc_content = sum(r.gc_content * r.length for r in doc.records) / total_bp if total_bp > 0 else None

    def _summary_count(key: str, fallback: int) -> int:
        return int(summary[key]) if key in summary else fallback

    headline_counts = {
        "cds": _summary_count("cds", len(cdss)),
        "trna": _summary_count("trna", counts["trna"]),
        "tmrna": _summary_count("tmrna", counts["tmrna"]),
        "rrna": _summary_count("rrna", counts["rrna"]),
        "ncrna": _summary_count("ncrna", counts["ncrna"]),
        "regulatory": _summary_count("regulatory", counts["regulatory"]),
        "crispr": _summary_count("crispr", 0),
        "pseudogene": _summary_count("pseudogene", len(pseudogene_ids)),
        "sorf": _summary_count("sorf", 0),
        "oric": _summary_count("oric", 0),
        "oriv": _summary_count("oriv", 0),
        "orit": _summary_count("orit", 0),
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


def write_outputs(rows: list[dict[str, object]], csv_path: Path, tsv_path: Path, md_path: Path) -> None:
    """Write the report in CSV, TSV, and Markdown formats."""
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


def batch_summary(
    inputs: list[str | Path],
    csv_out: str | Path = "bakta_summary.csv",
    tsv_out: str | Path = "bakta_summary.tsv",
    md_out: str | Path = "bakta_summary.md",
) -> list[dict[str, object]]:
    files = discover_genbank_files(inputs)
    if not files:
        print("ERROR: No GenBank files found in specified inputs.", file=sys.stderr)
        return []

    rows: list[dict[str, object]] = []
    for f in files:
        try:
            row = summarise_file(f)
            rows.append(row)
        except Exception as err:
            print(f"Warning: Failed to summarise {f}: {err}", file=sys.stderr)

    if rows:
        write_outputs(rows, Path(csv_out), Path(tsv_out), Path(md_out))
        print(f"Generated summary across {len(rows)} isolate(s):")
        print(f"  CSV: {csv_out}")
        print(f"  TSV: {tsv_out}")
        print(f"  MD : {md_out}")

    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description="Parse Bakta GenBank annotations and write comparison tables.")
    parser.add_argument("inputs", nargs="+", help="GenBank files or directories containing .gbff/.gbk/.gb files")
    parser.add_argument("--csv", default="bakta_summary.csv", help="Output CSV path (default: bakta_summary.csv)")
    parser.add_argument("--tsv", default="bakta_summary.tsv", help="Output TSV path (default: bakta_summary.tsv)")
    parser.add_argument("--md", default="bakta_summary.md", help="Output Markdown path (default: bakta_summary.md)")
    args = parser.parse_args()

    rows = batch_summary(args.inputs, csv_out=args.csv, tsv_out=args.tsv, md_out=args.md)
    return 0 if rows else 1


if __name__ == "__main__":
    sys.exit(main())
