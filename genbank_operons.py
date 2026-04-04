#!/usr/bin/env python3
"""Flag consecutive same-strand genes within max_gap bp -- operon candidates."""
import sys
from genbank_parser import parse_features, get_qual


def operon_candidates(filepath, max_gap=150):
    features = parse_features(filepath)
    cdss = sorted(
        [f for f in features if f['type'] == 'CDS'],
        key=lambda f: f['start']
    )

    pairs = []
    for i in range(len(cdss) - 1):
        a, b = cdss[i], cdss[i + 1]
        if a['contig'] == b['contig'] and a['strand'] == b['strand']:
            gap = b['start'] - a['end'] - 1
            if 0 <= gap <= max_gap:
                pairs.append((a, b, gap))

    print(f"-- Operon candidates (same-strand, gap <= {max_gap} bp) --")
    print(f"{'Locus A':18s}  {'Gene':6s}  ->  {'Locus B':18s}  {'Gene':6s}  {'Gap':>5s}  Strand  Products")
    print("-" * 110)
    for a, b, gap in pairs:
        ta = get_qual(a, 'locus_tag', '?')
        tb = get_qual(b, 'locus_tag', '?')
        ga = get_qual(a, 'gene') or '-'
        gb = get_qual(b, 'gene') or '-'
        pa = get_qual(a, 'product')[:30]
        pb = get_qual(b, 'product')[:30]
        print(f"  {ta:18s}  {ga:6s}  ->  {tb:18s}  {gb:6s}  {gap:>4d}bp  {a['strand']:>3s}    {pa} | {pb}")

    print(f"\nTotal candidate pairs: {len(pairs)}")

    # Also flag tight clusters (3+ consecutive genes within threshold)
    if pairs:
        print("\n-- Tight clusters (>=3 consecutive co-directional genes) --")
        clusters = []
        current_cluster = [pairs[0][0]]
        for a, b, gap in pairs:
            if current_cluster[-1] is a:
                current_cluster.append(b)
            else:
                if len(current_cluster) >= 3:
                    clusters.append(current_cluster)
                current_cluster = [a, b]
        if len(current_cluster) >= 3:
            clusters.append(current_cluster)

        if clusters:
            for cl in clusters:
                tags = [get_qual(f, 'locus_tag', '?') for f in cl]
                genes = [get_qual(f, 'gene') or '-' for f in cl]
                span = f"{cl[0]['start']}..{cl[-1]['end']}"
                print(f"  [{cl[0]['strand']}] {' -> '.join(genes)}  ({span})")
                print(f"       Tags: {', '.join(tags)}")
        else:
            print("  (none found)")


if __name__ == '__main__':
    filepath = sys.argv[1]
    max_gap = int(sys.argv[2]) if len(sys.argv) > 2 else 150
    operon_candidates(filepath, max_gap)
