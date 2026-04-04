#!/usr/bin/env python3
"""
Discovery Mode: Scans a GenBank feature table for high-value biological
features using keyword heuristics, spatial clustering, and dark-matter
operon detection. Designed to help researchers identify interesting loci
in an uncharacterized isolate before targeted analysis.

Usage:
    python genbank_discover.py <file> [--cluster-gap 5000] [--operon-gap 150]
"""
import sys
from genbank_parser import parse_features, get_qual


# ============================================================
# Keyword dictionaries -- curated for Bakta/Prokka annotations
# ============================================================
# Each entry: (keyword, weight)
# Weight: 3 = highly specific, 2 = moderately specific, 1 = broad

CATEGORIES = {
    "Xenobiotic Degradation & Bioremediation": [
        ("alkane monooxygenase", 3), ("alkane hydrolase", 3), ("alkane hydroxylase", 3),
        ("alkb", 3), ("lada", 3), ("alma", 3),
        ("catechol dioxygenase", 3), ("protocatechuate", 3), ("gentisate", 3),
        ("ring-hydroxylating", 3), ("ring-cleaving", 3), ("benzoate", 2),
        ("biphenyl", 3), ("phenol hydroxylase", 3), ("toluene", 3),
        ("alkanesulfonate", 3), ("dibenzothiophene", 3), ("desulfurization", 3),
        ("dehalogenase", 3), ("haloalkane", 3), ("haloacid", 2),
        ("petase", 3), ("mhetase", 3), ("cutinase", 2),
        ("monooxygenase", 1), ("dioxygenase", 1),
    ],
    "Heavy Metal & Stress Resistance": [
        ("mercuric reductase", 3), ("merr", 3), ("mera", 3), ("merb", 3),
        ("organomercurial", 3),
        ("arsenate reductase", 3), ("arsenite", 3), ("arsb", 3), ("arsc", 3),
        ("arsr", 3),
        ("copper resistance", 3), ("copper efflux", 3), ("copr", 3), ("copa", 3),
        ("cusab", 2), ("multicopper oxidase", 2),
        ("cadmium", 3), ("czca", 3), ("czcb", 3), ("czcd", 3),
        ("chromate", 3), ("tellurite", 3), ("terb", 3), ("chrr", 3),
        ("cold-shock", 2), ("heat-shock", 2), ("universal stress", 2),
        ("oxidative stress", 2), ("superoxide dismutase", 2), ("catalase", 1),
    ],
    "Secondary Metabolites & BGCs": [
        ("non-ribosomal peptide", 3), ("nrps", 3), ("polyketide synthase", 3),
        ("pks", 3), ("condensation domain", 3), ("acyl carrier protein", 2),
        ("bacteriocin", 3), ("lantibiotic", 3), ("lanthipeptide", 3),
        ("microcin", 3), ("colicin", 3),
        ("siderophore", 3), ("enterobactin", 3), ("pyoverdine", 3),
        ("aerobactin", 3), ("yersiniabactin", 3), ("vibriobactin", 3),
        ("desferrioxamine", 3), ("dihydroxybenzoate", 2),
        ("phenazine", 3), ("indole", 2), ("terpene cyclase", 3),
        ("violacein", 3), ("prodigiosin", 3),
    ],
    "Mobilome & HGT Markers": [
        ("transposase", 3), ("insertion sequence", 3),
        ("is element", 3), ("tnp", 2),
        ("integrase", 3), ("site-specific recombinase", 2),
        ("recombinase", 1),
        ("phage portal", 3), ("phage tail", 3), ("phage capsid", 3),
        ("phage terminase", 3), ("phage baseplate", 3),
        ("phage head", 3), ("phage lysin", 3),
        ("conjugative", 3), ("type iv secretion", 3), ("relaxase", 3),
        ("tra protein", 2), ("mob protein", 2),
        ("toxin-antitoxin", 3), ("addiction module", 3),
        ("restriction endonuclease", 2), ("methyltransferase", 1),
    ],
}


def is_hypothetical(product):
    """Check if a protein product is functionally uncharacterized."""
    p = product.lower()
    return (
        "hypothetical" in p
        or "uncharacterized" in p
        or "domain-containing protein" in p
        or p.startswith("duf")
        or "domain of unknown function" in p
        or p == ""
    )


def classify_cds(f):
    """Classify a CDS against all categories. Returns list of (category, keyword, weight)."""
    gene = get_qual(f, 'gene').lower()
    product = get_qual(f, 'product').lower()
    search_text = gene + " " + product
    hits = []
    for cat, keywords in CATEGORIES.items():
        for kw, weight in keywords:
            if kw in search_text:
                hits.append((cat, kw, weight))
                break  # One hit per category max
    return hits


def spatial_cluster(items, max_gap):
    """Group items into spatial clusters based on max_gap between consecutive items."""
    if not items:
        return []
    sorted_items = sorted(items, key=lambda x: x['start'])
    clusters = [[sorted_items[0]]]
    for item in sorted_items[1:]:
        if item['start'] - clusters[-1][-1]['end'] <= max_gap:
            clusters[-1].append(item)
        else:
            clusters.append([item])
    return clusters


def discover(filepath, cluster_gap=5000, operon_gap=150):
    features = parse_features(filepath)
    cdss = sorted([f for f in features if f['type'] == 'CDS'], key=lambda f: f['start'])

    if not cdss:
        print("No CDS features found.")
        sys.exit(1)

    total_cds = len(cdss)
    hyp_count = sum(1 for f in cdss if is_hypothetical(get_qual(f, 'product')))

    # -- Classify every CDS --
    cat_hits = {cat: [] for cat in CATEGORIES}
    for f in cdss:
        hits = classify_cds(f)
        for cat, kw, weight in hits:
            cat_hits[cat].append({
                'feature': f, 'start': f['start'], 'end': f['end'],
                'tag': get_qual(f, 'locus_tag', '?'),
                'gene': get_qual(f, 'gene') or '-',
                'product': get_qual(f, 'product'),
                'keyword': kw, 'weight': weight,
            })

    # ============================================================
    # REPORT
    # ============================================================
    print("=" * 70)
    print("  GENOME DISCOVERY REPORT")
    print("=" * 70)
    print(f"  File           : {filepath}")
    print(f"  Total CDS      : {total_cds}")
    print(f"  Hypothetical   : {hyp_count} ({100*hyp_count//total_cds}%)")
    print()

    # -- Summary tally --
    print("-- Summary --")
    for cat in CATEGORIES:
        n = len(cat_hits[cat])
        high = sum(1 for h in cat_hits[cat] if h['weight'] >= 3)
        bar = "#" * min(n, 40)
        print(f"  {cat:42s}  {n:3d} hits ({high} high-confidence)")
        if bar:
            print(f"  {'':42s}  {bar}")
    print()

    # -- Per-category detail with spatial clustering --
    for cat in CATEGORIES:
        items = cat_hits[cat]
        if not items:
            print(f"-- {cat} --")
            print("  (no hits)")
            print()
            continue

        clusters = spatial_cluster(items, cluster_gap)
        isolated = [c for c in clusters if len(c) == 1]
        grouped = [c for c in clusters if len(c) >= 2]

        print(f"-- {cat} ({len(items)} hits, {len(grouped)} clusters) --")
        print()

        for ci, cluster in enumerate(grouped, 1):
            span_start = cluster[0]['start']
            span_end = cluster[-1]['end']
            span_kb = (span_end - span_start + 1) / 1000
            print(f"  ** Cluster {ci}: {len(cluster)} genes spanning {span_kb:.1f} kb "
                  f"({span_start:,}..{span_end:,}) **")
            for h in cluster:
                conf = "***" if h['weight'] >= 3 else " * " if h['weight'] >= 2 else " . "
                print(f"    {conf} {h['tag']:18s} {h['gene']:8s}  "
                      f"{h['product'][:45]:45s}  [{h['keyword']}]")
            print()

        if isolated:
            print(f"  Isolated hits ({len(isolated)}):")
            for cl in isolated:
                h = cl[0]
                conf = "***" if h['weight'] >= 3 else " * " if h['weight'] >= 2 else " . "
                print(f"    {conf} {h['tag']:18s} {h['gene']:8s}  "
                      f"{h['product'][:45]:45s}  [{h['keyword']}]")
            print()

    # ============================================================
    # DARK MATTER OPERONS
    # ============================================================
    print("=" * 70)
    print("  DARK MATTER OPERONS")
    print("=" * 70)
    print()

    all_clusters = []
    current = [cdss[0]]
    for i in range(1, len(cdss)):
        a, b = cdss[i - 1], cdss[i]
        gap = b['start'] - a['end'] - 1
        if a['contig'] == b['contig'] and a['strand'] == b['strand'] and 0 <= gap <= operon_gap:
            current.append(b)
        else:
            if len(current) >= 3:
                all_clusters.append(current)
            current = [b]
    if len(current) >= 3:
        all_clusters.append(current)

    pure_dark = []
    mixed_dark = []

    for cl in all_clusters:
        products = [get_qual(f, 'product') for f in cl]
        hyp_frac = sum(1 for p in products if is_hypothetical(p)) / len(products)
        if hyp_frac == 1.0:
            pure_dark.append(cl)
        elif hyp_frac >= 0.6 and len(cl) >= 4:
            mixed_dark.append((cl, hyp_frac))

    if pure_dark:
        print(f"  -- Pure dark matter ({len(pure_dark)} clusters, ALL genes unknown) --")
        for i, cl in enumerate(pure_dark, 1):
            span = f"{cl[0]['start']:,}..{cl[-1]['end']:,}"
            print(f"  Cluster {i}: {len(cl)} genes, {cl[0]['strand']} strand, span {span}")
            for f in cl:
                tag = get_qual(f, 'locus_tag', '?')
                prod = get_qual(f, 'product') or '(hypothetical)'
                print(f"    {tag:18s}  {prod[:55]}")
            print()
    else:
        print("  -- Pure dark matter --")
        print("  (none found)")
        print()

    if mixed_dark:
        print(f"  -- Mixed dark matter ({len(mixed_dark)} clusters, >=60% unknown + anchor genes) --")
        for i, (cl, frac) in enumerate(mixed_dark, 1):
            span = f"{cl[0]['start']:,}..{cl[-1]['end']:,}"
            print(f"  Cluster {i}: {len(cl)} genes ({frac:.0%} hypothetical), "
                  f"{cl[0]['strand']} strand, span {span}")
            for f in cl:
                tag = get_qual(f, 'locus_tag', '?')
                prod = get_qual(f, 'product')
                marker = "  ???" if is_hypothetical(prod) else "  <--"
                print(f"    {tag:18s}  {(prod or '(hypothetical)')[:50]:50s} {marker}")
            print()
    else:
        print("  -- Mixed dark matter --")
        print("  (none found)")
        print()

    # ============================================================
    # CROSS-ANALYSIS: Mobilome-adjacent functional islands
    # ============================================================
    mob_items = cat_hits.get("Mobilome & HGT Markers", [])
    func_items = []
    for cat in ["Xenobiotic Degradation & Bioremediation",
                "Heavy Metal & Stress Resistance",
                "Secondary Metabolites & BGCs"]:
        func_items.extend(cat_hits.get(cat, []))

    if mob_items and func_items:
        print("=" * 70)
        print("  MOBILOME-ADJACENT FUNCTIONAL ISLANDS")
        print("=" * 70)
        print()
        print("  (Functional genes within 10 kb of a mobilome marker)")
        print()

        found_islands = False
        for mob in mob_items:
            nearby = [fi for fi in func_items
                      if abs(fi['start'] - mob['start']) <= 10000]
            if nearby:
                found_islands = True
                print(f"  Near {mob['tag']} ({mob['gene']}, {mob['keyword']}):")
                for fi in nearby:
                    dist = fi['start'] - mob['start']
                    direction = f"+{dist//1000}kb" if dist >= 0 else f"{dist//1000}kb"
                    print(f"    {direction:>6s}  {fi['tag']:18s}  {fi['product'][:45]}")
                print()

        if not found_islands:
            print("  (No functional genes found within 10 kb of mobilome markers)")
            print()

    print("=" * 70)
    print("  END OF DISCOVERY REPORT")
    print("=" * 70)


if __name__ == '__main__':
    filepath = sys.argv[1] if len(sys.argv) > 1 else input("File path: ")
    cluster_gap = 5000
    operon_gap = 150
    for i, arg in enumerate(sys.argv):
        if arg == '--cluster-gap' and i + 1 < len(sys.argv):
            cluster_gap = int(sys.argv[i + 1])
        if arg == '--operon-gap' and i + 1 < len(sys.argv):
            operon_gap = int(sys.argv[i + 1])
    discover(filepath, cluster_gap, operon_gap)
