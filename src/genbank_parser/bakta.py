"""Build comparable genome and annotation summaries from Bakta GenBank output."""
from __future__ import annotations

import argparse
import csv
import io
import json
import re
import sys
from collections import Counter
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .cli_io import (
    OutputError,
    eprint,
    paths_same,
    publish_directory,
    publish_file_set,
)
from .discovery import DISCOVERABLE_SUFFIXES, discover_inputs
from .io import extract_xrefs, get_qual, read_genbank
from .model import GenBankFeature

GENBANK_SUFFIXES = set(DISCOVERABLE_SUFFIXES)
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
    return [item.resolved_path for item in discover_inputs(inputs)]


def parse_bakta_summary(gbff_path: Path) -> dict[str, int | float]:
    """Read stable headline values from a same-named Bakta text report."""
    source_path = gbff_path.with_suffix("") if gbff_path.suffix.casefold() == ".gz" else gbff_path
    txt_path = source_path.with_suffix(".txt")
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

    source_name = gbff_path.name[:-3] if gbff_path.name.casefold().endswith(".gz") else gbff_path.name
    return {
        "Sample / Isolate": Path(source_name).stem,
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


def _render_delimited(rows: list[dict[str, object]], delimiter: str) -> str:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=FIELDNAMES, delimiter=delimiter, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue()


def _render_markdown(rows: list[dict[str, object]]) -> str:
    markdown_lines = ["# Bakta Isolates Summary Table", ""]
    markdown_lines.append("| " + " | ".join(FIELDNAMES) + " |")
    markdown_lines.append("| " + " | ".join("---" for _ in FIELDNAMES) + " |")
    markdown_lines.extend(
        "| " + " | ".join(_markdown_escape(row[column]) for column in FIELDNAMES) + " |"
        for row in rows
    )
    return "\n".join(markdown_lines) + "\n"


@dataclass(frozen=True)
class BatchSummaryResult:
    rows: tuple[dict[str, object], ...]
    failed_sources: tuple[dict[str, str], ...] = ()

    @property
    def exit_code(self) -> int:
        return 0 if self.rows and not self.failed_sources else 3


def write_outputs(
    rows: list[dict[str, object]],
    csv_path: Path,
    tsv_path: Path,
    md_path: Path,
    *,
    force: bool = False,
    inputs: Iterable[str | Path] = (),
    failures: Iterable[dict[str, str]] = (),
) -> None:
    """Render and publish the legacy report set as one transaction."""

    inputs = tuple(inputs)
    rendered = {
        csv_path: _render_delimited(rows, ","),
        tsv_path: _render_delimited(rows, "\t"),
        md_path: _render_markdown(rows),
    }
    failure_rows = tuple(failures)
    if failure_rows:
        rendered[md_path.with_name(f"{md_path.stem}.failures.json")] = (
            json.dumps(
                {"schema_version": "gbparse.batch-summary.v1", "failures": list(failure_rows)},
                ensure_ascii=False,
                allow_nan=False,
                sort_keys=True,
                indent=2,
            )
            + "\n"
        )
    output_paths = tuple(rendered)
    for index, path in enumerate(output_paths):
        if any(paths_same(path, other) for other in output_paths[index + 1 :]):
            raise OutputError(f"batch-summary outputs must be distinct: {path}")
    publish_file_set(rendered, force=force, inputs=inputs)


def run_batch_summary(
    inputs: list[str | Path],
    csv_out: str | Path = "bakta_summary.csv",
    tsv_out: str | Path = "bakta_summary.tsv",
    md_out: str | Path = "bakta_summary.md",
    *,
    output_dir: str | Path | None = None,
    force: bool = False,
    on_error: str = "fail",
) -> BatchSummaryResult:
    if on_error not in {"fail", "skip"}:
        raise ValueError("on_error must be fail or skip")
    if output_dir:
        output_path = Path(output_dir)
        staging_siblings = tuple(output_path.parent.glob(f".{output_path.name}.*")) if output_path.parent.exists() else ()
        excluded = (output_path, *staging_siblings)
    else:
        excluded = (csv_out, tsv_out, md_out)
    files = [item.resolved_path for item in discover_inputs(inputs, exclude=excluded)]
    if not files:
        eprint("ERROR: no GenBank files found in specified inputs")
        return BatchSummaryResult(())

    rows: list[dict[str, object]] = []
    failures: list[dict[str, str]] = []
    for file_path in files:
        try:
            rows.append(summarise_file(file_path))
        except (OSError, ValueError, KeyError, TypeError) as exc:
            failures.append({"source": str(file_path), "error": str(exc)})
            eprint(f"WARNING: failed to summarise {file_path}: {exc}")
    if not rows:
        eprint("ERROR: no usable GenBank summaries were produced")
        return BatchSummaryResult((), tuple(failures))
    if failures and on_error == "fail":
        eprint("ERROR: one or more GenBank summaries failed; no batch-summary outputs were published")
        return BatchSummaryResult(tuple(rows), tuple(failures))
    if output_dir is not None:
        files_payload = {
            "bakta_summary.csv": _render_delimited(rows, ","),
            "bakta_summary.tsv": _render_delimited(rows, "\t"),
            "bakta_summary.md": _render_markdown(rows),
        }
        if failures:
            files_payload["bakta_summary.failures.json"] = (
                json.dumps(
                    {"schema_version": "gbparse.batch-summary.v1", "failures": failures},
                    ensure_ascii=False,
                    allow_nan=False,
                    sort_keys=True,
                    indent=2,
                )
                + "\n"
            )
        publish_directory(files_payload, output_dir, force=force, inputs=inputs)
    else:
        write_outputs(
            rows,
            Path(csv_out),
            Path(tsv_out),
            Path(md_out),
            force=force,
            inputs=inputs,
            failures=failures,
        )
    eprint(f"Generated summary across {len(rows)} isolate(s)")
    return BatchSummaryResult(tuple(rows), tuple(failures))


def batch_summary(
    inputs: list[str | Path],
    csv_out: str | Path = "bakta_summary.csv",
    tsv_out: str | Path = "bakta_summary.tsv",
    md_out: str | Path = "bakta_summary.md",
) -> list[dict[str, object]]:
    return list(run_batch_summary(inputs, csv_out, tsv_out, md_out).rows)


def main() -> int:
    parser = argparse.ArgumentParser(description="Parse Bakta GenBank annotations and write comparison tables.")
    parser.add_argument("inputs", nargs="+", help="GenBank files or directories containing .gbff/.gbk/.gb files")
    parser.add_argument("--csv", default="bakta_summary.csv", help="Output CSV path (default: bakta_summary.csv)")
    parser.add_argument("--tsv", default="bakta_summary.tsv", help="Output TSV path (default: bakta_summary.tsv)")
    parser.add_argument("--md", default="bakta_summary.md", help="Output Markdown path (default: bakta_summary.md)")
    parser.add_argument("--output-dir", help="Publish all three tables as one atomic directory")
    parser.add_argument("--force", action="store_true", help="Replace existing output files or directory")
    parser.add_argument("--on-error", choices=("fail", "skip"), default="fail")
    args = parser.parse_args()

    result = run_batch_summary(
        args.inputs,
        csv_out=args.csv,
        tsv_out=args.tsv,
        md_out=args.md,
        output_dir=args.output_dir,
        force=args.force,
        on_error=args.on_error,
    )
    return result.exit_code


if __name__ == "__main__":
    sys.exit(main())
