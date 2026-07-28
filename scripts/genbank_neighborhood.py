#!/usr/bin/env python3
"""Show gene neighborhood for a target locus tag.

Strictly respects contig boundaries: the neighborhood window is built only
from genes on the same contig as the target (Principle #5).
"""
import sys, argparse
from genbank_parser import parse_features, get_qual


def neighborhood(filepath, target_tag, window=5):
    features = parse_features(filepath)
    genes = [f for f in features if f['type'] in ('gene', 'CDS')]

    # Deduplicate by locus_tag, preferring CDS (richer annotation)
    seen = {}
    for f in genes:
        tag = get_qual(f, 'locus_tag')
        if tag and (tag not in seen or f['type'] == 'CDS'):
            seen[tag] = f

    # Locate the target to identify its contig BEFORE any sorting
    target_feature = None
    for f in seen.values():
        if get_qual(f, 'locus_tag') == target_tag:
            target_feature = f
            break

    if target_feature is None:
        print(f"Locus tag '{target_tag}' not found.")
        sys.exit(1)

    target_contig = target_feature['contig']

    # Sort only genes on the SAME contig to avoid cross-contig contamination
    ordered = sorted(
        [f for f in seen.values() if f['contig'] == target_contig],
        key=lambda f: f['start'],
    )

    target_idx = next(
        i for i, f in enumerate(ordered) if get_qual(f, 'locus_tag') == target_tag
    )

    lo = max(0, target_idx - window)
    hi = min(len(ordered), target_idx + window + 1)

    print(f"-- Gene neighborhood for {target_tag} (contig: {target_contig}, +/-{window}) --")
    print(f"{'#':>3}  {'Locus Tag':18s}  {'Gene':8s}  {'Strand':6s}  "
          f"{'Start':>10s}  {'End':>10s}  Product")
    print("-" * 100)
    for i in range(lo, hi):
        f = ordered[i]
        marker = " >>>" if i == target_idx else "    "
        print(f"{marker}{i - target_idx:+d}  "
              f"{get_qual(f, 'locus_tag'):18s}  "
              f"{get_qual(f, 'gene'):8s}  "
              f"{f['strand']:6s}  "
              f"{f['start']:>10,d}  "
              f"{f['end']:>10,d}  "
              f"{get_qual(f, 'product')[:55]}")


def main():
    parser = argparse.ArgumentParser(
        description="Display gene neighborhood for a target locus tag."
    )
    parser.add_argument('input',      help="Input GenBank file")
    parser.add_argument('locus_tag',  help="Target locus tag")
    parser.add_argument('window',     nargs='?', type=int, default=5,
                        help="Number of flanking genes to show (default: 5)")
    args = parser.parse_args()
    neighborhood(args.input, args.locus_tag, args.window)


if __name__ == '__main__':
    main()
