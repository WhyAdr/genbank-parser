"""Extract record metadata (LOCUS, DEFINITION, SOURCE, ORGANISM, TOPOLOGY) via typed parser."""
from __future__ import annotations

import argparse
import io
from pathlib import Path
from typing import Any

from .io import read_genbank


def build_metadata_report(filepath: str | Path) -> dict[str, Any]:
    """Return record metadata without printing."""

    doc = read_genbank(filepath)
    records_info: list[dict[str, Any]] = []
    for rec in doc.records:
        strain = ""
        for feature in rec.features:
            if feature.type == "source":
                strain = feature.get_qual("strain") or feature.get_qual("isolate")
                if strain:
                    break
        records_info.append(
            {
                "locus": rec.id,
                "length": rec.length,
                "mol_type": rec.molecule_type or rec.annotations.get("molecule_type", "DNA"),
                "topology": rec.topology or "linear",
                "division": rec.division or rec.annotations.get("data_file_division", ""),
                "date": rec.date or rec.annotations.get("date", ""),
                "definition": rec.description,
                "accession": rec.annotations.get("accessions", [""])[0]
                if rec.annotations.get("accessions")
                else "",
                "organism": rec.annotations.get("organism", ""),
                "source": rec.annotations.get("source", ""),
                "strain": strain,
            }
        )
    return {
        "source": str(filepath),
        "record_count": len(records_info),
        "total_length": sum(int(record["length"]) for record in records_info),
        "records": records_info,
    }


def render_metadata(report: dict[str, Any], format_type: str = "text") -> str:
    if format_type == "json":
        import json

        return json.dumps(report, ensure_ascii=False, sort_keys=True, indent=2) + "\n"
    if format_type == "tsv":
        columns = ("locus", "length", "mol_type", "topology", "division", "date", "definition", "accession", "organism", "source", "strain")
        output = io.StringIO(newline="")
        output.write("\t".join(columns) + "\n")
        for record in report["records"]:
            output.write("\t".join(str(record.get(column, "")).replace("\t", " ") for column in columns) + "\n")
        return output.getvalue()
    lines = [
        "=" * 70,
        "  GENBANK METADATA REPORT",
        "=" * 70,
        f"  File    : {report['source']}",
        f"  Records : {report['record_count']}",
        "",
    ]
    for index, record in enumerate(report["records"], 1):
        lines.extend(
            [
                f"  Record {index}: {record['locus']}",
                f"    Length     : {record['length']:,} bp",
                f"    Topology   : {record['topology']}",
                f"    Mol. type  : {record['mol_type']}",
                f"    Definition : {str(record['definition'])[:80]}",
            ]
        )
        if record["organism"]:
            lines.append(f"    Organism   : {record['organism']}")
        if record["strain"]:
            lines.append(f"    Strain     : {record['strain']}")
        lines.append("")
    lines.extend([f"  Total: {report['record_count']} record(s), {report['total_length']:,} bp", "=" * 70])
    return "\n".join(lines) + "\n"


def extract_metadata(filepath: str | Path) -> list[dict[str, Any]]:
    doc = read_genbank(filepath)
    records_info: list[dict[str, Any]] = []

    print("=" * 70)
    print("  GENBANK METADATA REPORT")
    print("=" * 70)
    print(f"  File    : {filepath}")
    print(f"  Records : {len(doc.records)}")
    print()

    total_len = 0
    for i, rec in enumerate(doc.records, 1):
        rec_len = rec.length
        total_len += rec_len

        if rec_len >= 1_000_000:
            size_str = f"{rec_len / 1_000_000:.2f} Mb"
        elif rec_len >= 1_000:
            size_str = f"{rec_len / 1_000:.1f} kb"
        else:
            size_str = f"{rec_len} bp"

        organism = rec.annotations.get('organism', '')
        source = rec.annotations.get('source', '')
        strain = ''
        for f in rec.features:
            if f.type == 'source':
                strain = f.get_qual('strain') or f.get_qual('isolate')
                if strain:
                    break

        info = {
            'locus': rec.id,
            'length': rec_len,
            'mol_type': rec.molecule_type or rec.annotations.get('molecule_type', 'DNA'),
            'topology': rec.topology or 'linear',
            'division': rec.division or rec.annotations.get('data_file_division', ''),
            'date': rec.date or rec.annotations.get('date', ''),
            'definition': rec.description,
            'accession': rec.annotations.get('accessions', [''])[0] if rec.annotations.get('accessions') else '',
            'organism': organism,
            'source': source,
            'strain': strain,
        }
        records_info.append(info)

        print(f"  Record {i}: {rec.id}")
        print(f"    Length     : {size_str} ({rec_len:,} bp)")
        print(f"    Topology   : {info['topology']}")
        print(f"    Mol. type  : {info['mol_type']}")
        print(f"    Definition : {rec.description[:80]}")
        if organism:
            print(f"    Organism   : {organism}")
        if strain:
            print(f"    Strain     : {strain}")
        print()

    if total_len >= 1_000_000:
        total_str = f"{total_len / 1_000_000:.2f} Mb"
    else:
        total_str = f"{total_len:,} bp"

    print(f"  Total: {len(doc.records)} record(s), {total_str}")
    if len(doc.records) > 1:
        print("  WARNING: Multi-contig file. Downstream scripts must respect contig boundaries.")
    print("=" * 70)

    return records_info


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract top-level GenBank metadata and summary statistics.")
    parser.add_argument('input', help="Input GenBank file")
    args = parser.parse_args()
    extract_metadata(args.input)


if __name__ == '__main__':
    main()
