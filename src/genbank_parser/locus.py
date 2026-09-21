"""Single-locus deep-dive: display all qualifiers and cross-references for a target locus."""
from __future__ import annotations

import argparse
import io
import json
import sys
from pathlib import Path
from typing import Any

from .io import _source_label, extract_xrefs, read_genbank
from .model import GenBankFeature


def build_locus_report(filepath: str | Path, locus_tag: str) -> dict[str, Any]:
    """Return a JSON-safe locus report without printing or exiting."""

    doc = read_genbank(filepath)
    source_label = _source_label(doc, filepath)
    match = doc.find_locus(locus_tag)
    if match is None:
        for rec in doc.records:
            for feature in rec.features:
                if feature.gene and feature.gene.casefold() == locus_tag.casefold():
                    match = (rec, feature)
                    break
            if match:
                break
    if match is None:
        raise ValueError(f"Locus tag or gene '{locus_tag}' not found in {source_label}")
    record, feature = match
    xrefs = extract_xrefs(feature)
    return {
        "query": locus_tag,
        "record": record.id,
        "record_length": record.length,
        "type": feature.type,
        "feature_index": feature.feature_index,
        "locus_tag": feature.locus_tag or None,
        "gene": feature.gene or None,
        "product": feature.product or None,
        "start": feature.start,
        "end": feature.end,
        "strand": feature.strand_symbol,
        "strand_value": feature.strand,
        "length": feature.length,
        "genomic_span": feature.genomic_span,
        "compound": feature.is_compound,
        "segments": [{"start": start, "end": end} for start, end in (feature.join_segments if feature.is_compound else [(feature.start, feature.end)])],
        "partial_start": feature.is_partial_start,
        "partial_end": feature.is_partial_end,
        "pseudo": feature.is_pseudo,
        "codon_start": feature.codon_start,
        "transl_table": feature.transl_table,
        "xrefs": xrefs,
        "qualifiers": {key: list(values) for key, values in sorted(feature.qualifiers.items())},
    }


def render_locus_report(report: dict[str, Any], format_type: str = "text") -> str:
    if format_type == "json":
        return json.dumps(report, ensure_ascii=False, sort_keys=True, indent=2) + "\n"
    if format_type == "tsv":
        output = io.StringIO(newline="")
        output.write("field\tvalue\n")
        for key in sorted(report):
            value = report[key]
            if isinstance(value, (dict, list)):
                value = json.dumps(value, ensure_ascii=False, sort_keys=True)
            output.write(f"{key}\t{str(value if value is not None else '') .replace(chr(9), ' ')}\n")
        return output.getvalue()
    lines = [
        "=" * 70,
        f"  FEATURE DEEP-DIVE: {report['query']}",
        "=" * 70,
        f"  Record / Contig   : {report['record']} (length: {report['record_length']:,} bp)",
        f"  Feature type      : {report['type']}",
        f"  Coordinates       : {report['start']:,} .. {report['end']:,} ({report['strand']})",
        f"  Biological length : {report['length']:,} bp",
        f"  Genomic span      : {report['genomic_span']:,} bp",
        f"  Compound / Join   : {report['compound']}",
        f"  Partial coords    : start={report['partial_start']}, end={report['partial_end']}",
        f"  Pseudogene        : {report['pseudo']}",
        "",
        "-- Standard Qualifiers --",
    ]
    for key in ("gene", "product", "locus_tag"):
        if report[key]:
            lines.append(f"  {key:16s} : {report[key]}")
    lines.extend(["", "-- Cross-References --"])
    for key, values in report["xrefs"].items():
        if values:
            lines.append(f"  {key:16s} : {', '.join(values)}")
    lines.extend(["", "-- All Qualifiers --"])
    for key, values in report["qualifiers"].items():
        for value in values:
            lines.append(f"  /{key:16s} : {value}")
    lines.append("=" * 70)
    return "\n".join(lines) + "\n"


def inspect_locus(filepath: str | Path, locus_tag: str) -> GenBankFeature | None:
    doc = read_genbank(filepath)
    source_label = _source_label(doc, filepath)
    match = doc.find_locus(locus_tag)

    if match is None:
        # Check by gene name
        for rec in doc.records:
            for f in rec.features:
                if f.gene and f.gene.casefold() == locus_tag.casefold():
                    match = (rec, f)
                    break
            if match:
                break

    if match is None:
        print(f"ERROR: Locus tag or gene '{locus_tag}' not found in {source_label}", file=sys.stderr)
        sys.exit(1)

    rec, f = match
    xrefs = extract_xrefs(f)

    print("=" * 70)
    print(f"  FEATURE DEEP-DIVE: {locus_tag}")
    print("=" * 70)
    print(f"  Record / Contig   : {rec.id} (length: {rec.length:,} bp)")
    print(f"  Feature type      : {f.type}")
    print(f"  Coordinates       : {f.start:,} .. {f.end:,} ({f.strand_symbol})")
    print(f"  Biological length : {f.length:,} bp ({f.length // 3} aa)")
    print(f"  Genomic span      : {f.genomic_span:,} bp")
    print(f"  Compound / Join   : {f.is_compound}")
    if f.is_compound:
        print(f"  Join segments     : {f.join_segments}")
    print(f"  Partial coords    : start={f.is_partial_start}, end={f.is_partial_end}")
    print(f"  Pseudogene        : {f.is_pseudo}")
    print()
    print("-- Standard Qualifiers --")
    if f.gene:
        print(f"  Gene symbol       : {f.gene}")
    if f.product:
        print(f"  Product           : {f.product}")
    if f.protein_id:
        print(f"  Protein ID        : {f.protein_id}")
    print(f"  Codon start       : {f.codon_start}")
    print(f"  Transl table      : {f.transl_table}")
    print()
    print("-- Cross-References --")
    if xrefs['kegg_kos']:
        print(f"  KEGG KO           : {', '.join(xrefs['kegg_kos'])}")
    if xrefs['ec_numbers']:
        print(f"  EC Number         : {', '.join(xrefs['ec_numbers'])}")
    if xrefs['cog_ids']:
        print(f"  COG IDs           : {', '.join(xrefs['cog_ids'])}")
    if xrefs['pfam']:
        print(f"  Pfam              : {', '.join(xrefs['pfam'])}")
    if xrefs['rfam']:
        print(f"  Rfam              : {', '.join(xrefs['rfam'])}")
    if xrefs['go_terms']:
        print(f"  GO Terms          : {', '.join(xrefs['go_terms'])}")
    if xrefs['db_xrefs']:
        print(f"  Other db_xrefs    : {', '.join(xrefs['db_xrefs'][:10])}")
    print()
    print("-- All Qualifiers --")
    for k, v in sorted(f.qualifiers.items()):
        if k == 'translation':
            trans = v[0]
            print(f"  /{k:16s} : {trans[:40]}... (length: {len(trans)} aa)")
        else:
            for item in v:
                print(f"  /{k:16s} : {item}")
    print("=" * 70)

    return f


def main() -> None:
    parser = argparse.ArgumentParser(description="Deep-dive inspection of a single locus tag or gene.")
    parser.add_argument('input', help="Input GenBank file")
    parser.add_argument('locus_tag', help="Target locus tag or gene name")
    args = parser.parse_args()

    inspect_locus(args.input, args.locus_tag)


if __name__ == '__main__':
    main()
