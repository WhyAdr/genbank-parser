"""CRISPR/Cas annotation scanner (not a sequence-level array detector)."""

from __future__ import annotations

import argparse
import io
import json
from pathlib import Path
from typing import Any

from .io import _source_label, read_genbank
from .model import GenBankFeature

CAS_KEYWORDS = [
    "cas1",
    "cas2",
    "cas4",
    "cas3",
    "cas5",
    "cas6",
    "cas7",
    "cas8",
    "cas9",
    "cpf1",
    "cas12",
    "cas10",
    "csm",
    "cmr",
    "cas13",
    "c2c2",
    "crispr-associated",
    "crispr associated",
]

REPEAT_KEYWORDS = [
    "crispr",
    "direct repeat",
    "palindromic repeat",
]


def interval_distance(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    """Return zero for overlap, otherwise the gap between two intervals."""
    if a_end < b_start:
        return b_start - a_end
    if b_end < a_start:
        return a_start - b_end
    return 0


def _feature_summary(feature: GenBankFeature) -> dict[str, Any]:
    """Return a stable, JSON-serializable annotation summary."""
    return {
        "record": feature.record_id,
        "type": feature.type,
        "locus_tag": feature.locus_tag or None,
        "gene": feature.gene or None,
        "product": feature.product or None,
        "start": feature.start,
        "end": feature.end,
        "strand": feature.strand_symbol,
    }


def _is_cas(f: GenBankFeature) -> bool:
    gene = (f.gene or "").lower()
    product = (f.product or "").lower()
    text = f"{gene} {product}"
    return any(kw in text for kw in CAS_KEYWORDS)


def _is_crispr_array(f: GenBankFeature) -> bool:
    if f.type != "repeat_region":
        return False
    product = (f.product or "").lower()
    notes = " ".join(f.qualifiers.get("note", [])).lower()
    text = f"{product} {notes}"
    return any(kw in text for kw in REPEAT_KEYWORDS)


def build_crispr_report(filepath: str | Path, window: int = 15000) -> dict[str, Any]:
    """Build a CRISPR/Cas annotation report without printing."""

    doc = read_genbank(filepath)
    arrays = [feature for feature in doc.all_features if _is_crispr_array(feature)]
    cas_cdss = [feature for feature in doc.all_features if feature.type == "CDS" and _is_cas(feature)]
    colocalized: list[dict[str, Any]] = []
    for array in arrays:
        nearby = [
            cas
            for cas in cas_cdss
            if cas.record_id == array.record_id
            and interval_distance(array.start, array.end, cas.start, cas.end) <= window
        ]
        if nearby:
            colocalized.append(
                {
                    "array": _feature_summary(array),
                    "cas_genes": [_feature_summary(cas) for cas in nearby],
                }
            )
    return {
        "schema_version": "gbparse.crispr.v1",
        "source": _source_label(doc, filepath),
        "window": window,
        "arrays": [_feature_summary(array) for array in arrays],
        "cas_genes": [_feature_summary(cas) for cas in cas_cdss],
        "colocalized": colocalized,
        "colocalized_count": len(colocalized),
    }


def render_crispr_report(report: dict[str, Any], format_type: str = "text") -> str:
    if format_type == "json":
        return json.dumps(report, ensure_ascii=False, sort_keys=True, indent=2) + "\n"
    if format_type == "tsv":
        output = io.StringIO(newline="")
        output.write("row_type\trecord\ttype\tlocus_tag\tgene\tproduct\tstart\tend\tstrand\tarray_record\tarray_start\tarray_end\n")
        def write_row(values: tuple[object, ...]) -> None:
            output.write("\t".join("" if value is None else str(value) for value in values) + "\n")

        for array in report["arrays"]:
            write_row(("array", array["record"], array["type"], array["locus_tag"], array["gene"], array["product"], array["start"], array["end"], array["strand"], "", "", ""))
        for cas in report["cas_genes"]:
            write_row(("cas", cas["record"], cas["type"], cas["locus_tag"], cas["gene"], cas["product"], cas["start"], cas["end"], cas["strand"], "", "", ""))
        for pair in report["colocalized"]:
            array = pair["array"]
            for cas in pair["cas_genes"]:
                write_row(("link", cas["record"], cas["type"], cas["locus_tag"], cas["gene"], cas["product"], cas["start"], cas["end"], cas["strand"], array["record"], array["start"], array["end"]))
        return output.getvalue()
    lines = [
        "=" * 70,
        "  CRISPR/Cas ANNOTATION SCANNER",
        "=" * 70,
        f"  File          : {report['source']}",
        f"  CRISPR arrays : {len(report['arrays'])}",
        f"  Cas CDSs      : {len(report['cas_genes'])}",
        f"  Co-localized  : {report['colocalized_count']}",
        "=" * 70,
    ]
    return "\n".join(lines) + "\n"


def detect_crispr(filepath: str | Path, window: int = 15000) -> dict[str, Any]:
    doc = read_genbank(filepath)
    all_features = doc.all_features

    arrays = [f for f in all_features if _is_crispr_array(f)]
    cas_cdss = [f for f in all_features if f.type == "CDS" and _is_cas(f)]

    print("=" * 70)
    print("  CRISPR/Cas ANNOTATION SCANNER")
    print("=" * 70)
    print(f"  File          : {_source_label(doc, filepath)}")
    print(f"  CRISPR arrays : {len(arrays)}")
    print(f"  Cas CDSs      : {len(cas_cdss)}")
    print()

    if not arrays and not cas_cdss:
        print("  No CRISPR-associated features detected.")
        print("=" * 70)
        return {
            "arrays": [],
            "cas_genes": [],
            "colocalized": [],
            "colocalized_count": 0,
        }

    if arrays:
        print("-- CRISPR Arrays --")
        for a in arrays:
            note = "; ".join(a.qualifiers.get("note", [])) or "-"
            print(
                f"  {a.record_id:20s}  {a.start:>10,d}..{a.end:>10,d}  {a.strand_symbol}  note: {note[:60]}"
            )
        print()

    if cas_cdss:
        print("-- Cas Genes --")
        for c in cas_cdss:
            tag = c.locus_tag or "-"
            gene = c.gene or "-"
            prod = (c.product or "-")[:40]
            print(
                f"  {c.record_id:20s}  {c.start:>10,d}..{c.end:>10,d}  {c.strand_symbol}  {tag:16s}  {gene:8s}  {prod}"
            )
        print()

    # Spatial co-localization
    colocalized: list[dict[str, Any]] = []
    for a in arrays:
        nearby_cas = [
            c
            for c in cas_cdss
            if c.record_id == a.record_id
            and interval_distance(a.start, a.end, c.start, c.end) <= window
        ]
        if nearby_cas:
            colocalized.append({"array": a, "cas_genes": nearby_cas})

    if colocalized:
        print(f"-- Spatial Co-localization (within {window:,} bp) --")
        for pair in colocalized:
            arr = pair["array"]
            print(
                f"  Array at {arr.record_id}:{arr.start}..{arr.end} linked with {len(pair['cas_genes'])} Cas gene(s):"
            )
            for c in pair["cas_genes"]:
                dist = interval_distance(arr.start, arr.end, c.start, c.end)
                print(
                    f"    -> {c.locus_tag or '-'} ({c.gene or '-'}) {c.product} [dist: {dist:,} bp]"
                )
        print()

    print("=" * 70)
    array_dicts = [_feature_summary(a) for a in arrays]
    cas_dicts = [_feature_summary(c) for c in cas_cdss]
    array_index = {id(feature): index for index, feature in enumerate(arrays)}
    colocalized_dicts = [
        {
            "array": array_dicts[array_index[id(pair["array"])]],
            "cas_genes": [_feature_summary(c) for c in pair["cas_genes"]],
        }
        for pair in colocalized
    ]
    return {
        "arrays": array_dicts,
        "cas_genes": cas_dicts,
        "colocalized": colocalized_dicts,
        "colocalized_count": len(colocalized_dicts),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Scan CRISPR/Cas annotations and spatial proximity."
    )
    parser.add_argument("input", help="Input GenBank file")
    parser.add_argument(
        "--window",
        type=int,
        default=15000,
        help="Window size in bp for array-Cas linking (default: 15000)",
    )
    args = parser.parse_args()

    detect_crispr(args.input, window=args.window)


if __name__ == "__main__":
    main()
