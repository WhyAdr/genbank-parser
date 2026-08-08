"""Export GenBank annotations to a clean tabular TSV format."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys

from .io import extract_xrefs, get_qual, read_genbank
from .model import GenBankDocument

FIELDNAMES = [
    'contig', 'start', 'end', 'strand', 'type',
    'locus_tag', 'gene', 'product', 'ec_number',
    'cog', 'kegg_ko', 'pfam', 'rfam', 'db_xref', 'protein_id'
]


def export_annotations_tsv(filepath: str | Path, output_path: str | Path | None = None) -> list[dict[str, str]]:
    doc = read_genbank(filepath)
    rows: list[dict[str, str]] = []

    for f in doc.all_features:
        xrefs = extract_xrefs(f)
        row = {
            'contig': f.record_id,
            'start': str(f.start),
            'end': str(f.end),
            'strand': f.strand_symbol,
            'type': f.type,
            'locus_tag': f.locus_tag or '',
            'gene': f.gene or '',
            'product': f.product or '',
            'ec_number': ';'.join(xrefs['ec_numbers']),
            'cog': ';'.join(xrefs['cog_ids']),
            'kegg_ko': ';'.join(xrefs['kegg_kos']),
            'pfam': ';'.join(xrefs['pfam']),
            'rfam': ';'.join(xrefs['rfam']),
            'db_xref': ';'.join(xrefs['db_xrefs']),
            'protein_id': f.protein_id or '',
        }
        rows.append(row)

    if output_path:
        out_p = Path(output_path)
        with out_p.open('w', newline='', encoding='utf-8') as fh:
            writer = csv.DictWriter(fh, fieldnames=FIELDNAMES, delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)
        print(f"Exported {len(rows):,} annotations to {output_path}")
    else:
        writer = csv.DictWriter(sys.stdout, fieldnames=FIELDNAMES, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)

    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Export GenBank feature table to a tab-delimited TSV file.")
    parser.add_argument('input', help="Input GenBank file")
    parser.add_argument('output', nargs='?', help="Output TSV path (default: stdout)")
    args = parser.parse_args()

    export_annotations_tsv(args.input, args.output)


if __name__ == '__main__':
    main()
