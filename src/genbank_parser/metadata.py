"""Extract record metadata (LOCUS, DEFINITION, SOURCE, ORGANISM, TOPOLOGY) via typed parser."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Any

from .io import read_genbank
from .model import GenBankDocument, GenBankRecord


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
