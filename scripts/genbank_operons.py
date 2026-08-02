#!/usr/bin/env python3
"""Flag consecutive same-strand genes within max_gap bp -- operon candidates.

Contig boundaries are respected: pairs from different contigs are never linked.
"""
import sys, argparse
from genbank_parser import parse_features, get_qual


def operon_candidates(filepath, max_gap=150):
    features = parse_features(filepath)
    cdss = sorted(
        [f for f in features if f['type'] == 'CDS'],
        key=lambda f: (f['contig'], f['start']),
    )

    pairs = []
    for i in range(len(cdss) - 1):
        a, b = cdss[i], cdss[i + 1]
        if a['contig'] == b['contig'] and a['strand'] == b['strand']:
            gap = b['start'] - a['end'] - 1
            if 0 <= gap <= max_gap:
                pairs.append((a, b, gap))

    print(f"-- Operon candidates (same-strand, gap <= {max_gap} bp) --")
    print(f"{'Locus A':18s}  {'Gene':6s}  ->  {'Locus B':18s}  {'Gene':6s}  "
          f"{'Gap':>5s}  Strand  Products")
    print("-" * 110)
    for a, b, gap in pairs:
        ta = get_qual(a, 'locus_tag', '?')
        tb = get_qual(b, 'locus_tag', '?')
        ga = get_qual(a, 'gene') or '-'
        gb = get_qual(b, 'gene') or '-'
        pa = get_qual(a, 'product')[:30]
        pb = get_qual(b, 'product')[:30]
        print(f"  {ta:18s}  {ga:6s}  ->  {tb:18s}  {gb:6s}  "
              f"{gap:>4d}bp  {a['strand']:>3s}    {pa} | {pb}")

    print(f"\nTotal candidate pairs: {len(pairs)}")

    # Tight clusters (3+ consecutive co-directional genes)
    if pairs:
        print("\n-- Tight clusters (>=3 consecutive co-directional genes) --")
        clusters = []
        current_cluster = [pairs[0][0], pairs[0][1]]
        for a, b, gap in pairs[1:]:
            # Unique feature index from parse_features() ensures exact identity match
            if current_cluster[-1]['line'] == a['line']:
                current_cluster.append(b)
            else:
                if len(current_cluster) >= 3:
                    clusters.append(current_cluster)
                current_cluster = [a, b]
        if len(current_cluster) >= 3:
            clusters.append(current_cluster)

        if clusters:
            for cl in clusters:
                tags  = [get_qual(f, 'locus_tag', '?') for f in cl]
                genes = [get_qual(f, 'gene') or '-' for f in cl]
                span  = f"{cl[0]['start']:,}..{cl[-1]['end']:,}  (contig: {cl[0]['contig']})"
                print(f"  [{cl[0]['strand']}] {' -> '.join(genes)}  ({span})")
                print(f"       Tags: {', '.join(tags)}")
        else:
            print("  (none found)")


def main():
    parser = argparse.ArgumentParser(
        description="Identify operon candidates: consecutive same-strand genes within a gap threshold."
    )
    parser.add_argument('input',   help="Input GenBank file")
    parser.add_argument('max_gap', nargs='?', type=int, default=150,
                        help="Maximum intergenic gap in bp (default: 150)")
    args = parser.parse_args()
    operon_candidates(args.input, args.max_gap)


if __name__ == '__main__':
    main()
