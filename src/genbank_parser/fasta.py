"""Export annotated protein sequences from GenBank translation qualifiers to FASTA (.faa)."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

from .io import read_genbank
from .model import GenBankDocument


def export_protein_fasta(filepath: str | Path, output_path: str | Path | None = None) -> int:
    doc = read_genbank(filepath)
    cdss = [f for f in doc.all_features if f.type == 'CDS']

    lines: list[str] = []
    written = 0

    for f in cdss:
        trans = f.translation
        if not trans:
            continue

        tag = f.locus_tag or f"cds_{f.feature_index}"
        gene = f.gene
        prod = f.product or 'hypothetical protein'

        header = f">{tag}"
        if gene:
            header += f" gene={gene}"
        header += f" product={prod} [{f.record_id}:{f.start}..{f.end}({f.strand_symbol})]"

        lines.append(header)
        for i in range(0, len(trans), 60):
            lines.append(trans[i:i + 60])
        written += 1

    content = '\n'.join(lines) + '\n' if lines else ''

    if output_path:
        out_p = Path(output_path)
        out_p.write_text(content, encoding='utf-8')
        print(f"Exported {written:,} protein sequences to {output_path}")
    else:
        sys.stdout.write(content)

    return written


def main() -> None:
    parser = argparse.ArgumentParser(description="Export all CDS translations as protein FASTA.")
    parser.add_argument('input', help="Input GenBank file")
    parser.add_argument('output', nargs='?', help="Output FASTA path (.faa, default: stdout)")
    args = parser.parse_args()

    export_protein_fasta(args.input, args.output)


if __name__ == '__main__':
    main()
