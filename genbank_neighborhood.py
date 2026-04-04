#!/usr/bin/env python3
"""Show gene neighborhood for a target locus tag."""
import sys
from genbank_parser import parse_features, get_qual


def neighborhood(filepath, target_tag, window=5):
    features = parse_features(filepath)
    genes = [f for f in features if f['type'] in ('gene', 'CDS')]
    # Deduplicate by locus_tag, keeping CDS (richer annotation)
    seen = {}
    for f in genes:
        tag = get_qual(f, 'locus_tag')
        if tag and (tag not in seen or f['type'] == 'CDS'):
            seen[tag] = f
    ordered = sorted(seen.values(), key=lambda f: f['start'])

    target_idx = None
    for i, f in enumerate(ordered):
        if get_qual(f, 'locus_tag') == target_tag:
            target_idx = i
            break

    if target_idx is None:
        print(f"Locus tag '{target_tag}' not found.")
        sys.exit(1)

    lo = max(0, target_idx - window)
    hi = min(len(ordered), target_idx + window + 1)

    print(f"-- Gene neighborhood for {target_tag} (+/-{window}) --")
    print(f"{'#':>3}  {'Locus Tag':18s}  {'Gene':8s}  {'Strand':6s}  {'Start':>10s}  {'End':>10s}  Product")
    print("-" * 95)
    for i in range(lo, hi):
        f = ordered[i]
        marker = " >>>" if i == target_idx else "    "
        print(f"{marker}{i-target_idx:+d}  "
              f"{get_qual(f, 'locus_tag'):18s}  "
              f"{get_qual(f, 'gene'):8s}  "
              f"{f['strand']:6s}  "
              f"{f['start']:>10,d}  "
              f"{f['end']:>10,d}  "
              f"{get_qual(f, 'product')[:50]}")


if __name__ == '__main__':
    neighborhood(sys.argv[1], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 5)
