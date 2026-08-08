"""Multi-genome comparative presence/absence matrix for target marker genes/KOs."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys
from typing import Any, Sequence

from .io import extract_xrefs, get_qual, read_genbank
from .model import GenBankFeature


def compare_genomes(
    genome_inputs: Sequence[str | Path],
    targets: Sequence[str],
    output_path: str | Path | None = None,
) -> list[dict[str, Any]]:
    target_set = [t.strip() for t in targets if t.strip()]

    # Collect genomes
    genome_files: list[Path] = []
    for inp in genome_inputs:
        p = Path(inp)
        if p.is_file() and p.suffix.lower() in ('.gbff', '.gbk', '.gb'):
            genome_files.append(p)
        elif p.is_dir():
            for f in sorted(p.rglob('*')):
                if f.is_file() and f.suffix.lower() in ('.gbff', '.gbk', '.gb'):
                    genome_files.append(f)

    if not genome_files:
        print("ERROR: No compatible GenBank files found in input path(s).", file=sys.stderr)
        sys.exit(1)

    print("=" * 80)
    print("  MULTI-GENOME COMPARATIVE MARKER SCANNER")
    print("=" * 80)
    print(f"  Genomes scanned : {len(genome_files)}")
    print(f"  Target markers  : {', '.join(target_set)}")
    print()

    matrix_rows: list[dict[str, Any]] = []

    header = f"{'Genome / Isolate':30s}  " + "  ".join(f"{t[:10]:>10s}" for t in target_set)
    print(header)
    print("-" * len(header))

    for g_path in genome_files:
        doc = read_genbank(g_path)
        genome_name = g_path.stem

        marker_counts: dict[str, int] = {t: 0 for t in target_set}
        marker_tags: dict[str, list[str]] = {t: [] for t in target_set}

        for f in doc.all_features:
            if f.type != 'CDS':
                continue
            gene = (f.gene or '').casefold()
            prod = (f.product or '').casefold()
            xrefs = extract_xrefs(f)
            kos = [k.casefold() for k in xrefs['kegg_kos']]
            ecs = [e.casefold() for e in xrefs['ec_numbers']]
            cogs = [c.casefold() for c in xrefs['cog_ids']]

            for t in target_set:
                t_low = t.casefold()
                if (
                    t_low == gene
                    or t_low in prod
                    or t_low in kos
                    or t_low in ecs
                    or t_low in cogs
                ):
                    marker_counts[t] += 1
                    marker_tags[t].append(f.locus_tag or f"{f.record_id}:{f.start}")

        row: dict[str, Any] = {'Genome': genome_name, 'Path': str(g_path)}
        row.update({t: marker_counts[t] for t in target_set})
        matrix_rows.append(row)

        counts_str = "  ".join(f"{marker_counts[t]:>10d}" for t in target_set)
        print(f"{genome_name:30s}  {counts_str}")

    print("=" * 80)

    if output_path:
        out_p = Path(output_path)
        with out_p.open('w', newline='', encoding='utf-8') as fh:
            writer = csv.DictWriter(fh, fieldnames=['Genome', 'Path'] + target_set, delimiter='\t')
            writer.writeheader()
            writer.writerows(matrix_rows)
        print(f"\nWrote comparative matrix TSV to {output_path}")

    return matrix_rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Multi-genome marker presence/absence matrix.")
    parser.add_argument('genomes', nargs='+', help="GenBank files or directory containing .gbff/.gbk files")
    parser.add_argument('--targets', required=True, help="Comma-separated marker names, genes, KOs, or ECs (e.g. 'ladA,ssuD,K20938')")
    parser.add_argument('--output', help="Output matrix TSV path")
    args = parser.parse_args()

    targets = [t.strip() for t in args.targets.split(',') if t.strip()]
    compare_genomes(args.genomes, targets, output_path=args.output)


if __name__ == '__main__':
    main()
