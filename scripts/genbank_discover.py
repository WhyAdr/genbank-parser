#!/usr/bin/env python3
"""
Discovery Mode: Scans a GenBank feature table for high-value biological
features using keyword heuristics, spatial clustering, and dark-matter
operon detection. Designed to help researchers identify interesting loci
in an uncharacterized isolate before targeted analysis.

Usage:
    python genbank_discover.py <file> [--cluster-gap 5000] [--operon-gap 150]
                                      [--min-weight 1] [--format text|json|tsv]

Contig boundary correctness:
    All spatial operations (spatial_cluster, dark-matter operon detection,
    mobilome-adjacent island detection) group genes by contig first.
    Genes from different contigs are never considered spatially adjacent.
"""
import sys, argparse, collections, re, os, json
from genbank_parser import parse_features, get_qual


# Category-specific exclusion patterns — genes matching these substrings
# are skipped for the given category to prevent housekeeping false positives.
EXCLUDE_LISTS = {
    "Xenobiotic Degradation & Bioremediation": [
        "menaquinone", "o-succinylbenzoate",
    ],
    "Mobilome & HGT Markers": [
        "rrna methyltransferase", "trna methyltransferase",
        "dna methyltransferase",
    ],
}

# Pre-compiled DUF pattern for is_hypothetical()
_DUF_PAT = re.compile(r'\bduf\d+\b', re.IGNORECASE)


# ============================================================
# Keyword dictionaries -- curated for Bakta/Prokka annotations
# ============================================================
CATEGORIES = {
    "Xenobiotic Degradation & Bioremediation": [
        ("alkane monooxygenase", 3), ("alkane hydrolase", 3), ("alkane hydroxylase", 3),
        ("alkb", 3), ("lada", 3), ("alma", 3),
        ("catechol dioxygenase", 3), ("protocatechuate", 3), ("gentisate", 3),
        ("ring-hydroxylating", 3), ("ring-cleaving", 3),
        ("benzoate 1,2-dioxygenase", 3), ("benzoate catabolism", 3), ("benzoyl-coa", 2),
        ("biphenyl", 3), ("phenol hydroxylase", 3), ("toluene", 3),
        ("alkanesulfonate", 3), ("dibenzothiophene", 3), ("desulfurization", 3),
        ("dehalogenase", 3), ("haloalkane", 3), ("haloacid", 2),
        ("petase", 3), ("mhetase", 3), ("cutinase", 2),
        ("monooxygenase", 1), ("dioxygenase", 1),
    ],
    "Heavy Metal & Stress Resistance": [
        ("mercuric reductase", 3), ("merr", 3), ("merb", 3),
        ("mert", 3), ("merc", 3), ("merd", 2), ("merp", 2),
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
    """Classify a product annotation as functionally unknown.

    Tiered rules:
      1. Empty / bare "protein" — no annotation at all
      2. Explicit hypothetical/uncharacterized labels
      3. DUF#### domains (domain of unknown function)
    Does NOT flag "X-type domain-containing protein" — those describe
    characterised folds with known regulatory or enzymatic function.
    """
    p = product.lower().strip()
    return (
        not p                                    # empty
        or p == "protein"                        # bare "protein"
        or "hypothetical" in p
        or "uncharacterized" in p
        or _DUF_PAT.search(p) is not None        # DUF1234 etc.
        or "domain of unknown function" in p
    )


def classify_cds(f, min_weight=1):
    """Classify a CDS feature against the keyword dictionary.

    Uses word-boundary regex matching to avoid substring collisions
    (e.g. 'merA' must not match 'polymerase'). Applies per-category
    exclude-lists and minimum weight filtering.
    """
    gene        = get_qual(f, 'gene').lower()
    product     = get_qual(f, 'product').lower()
    search_text = gene + " " + product
    hits = []
    for cat, keywords in CATEGORIES.items():
        # Skip category if any exclude-list term matches
        excludes = EXCLUDE_LISTS.get(cat, [])
        if any(exc in search_text for exc in excludes):
            continue
        for kw, weight in keywords:
            if weight < min_weight:
                continue
            if re.search(r'\b' + re.escape(kw) + r'\b', search_text):
                hits.append((cat, kw, weight))
                break
    return hits


def spatial_cluster(items, max_gap):
    """Group items into spatial clusters, strictly per-contig.

    Items must carry a 'feature' key pointing to the parsed feature dict
    (which carries a 'contig' key).  Items from different contigs are
    never merged into the same cluster.
    """
    if not items:
        return []

    by_contig = collections.defaultdict(list)
    for item in items:
        by_contig[item['feature']['contig']].append(item)

    all_clusters = []
    for contig_items in by_contig.values():
        sorted_items = sorted(contig_items, key=lambda x: x['start'])
        clusters = [[sorted_items[0]]]
        for item in sorted_items[1:]:
            if item['start'] - clusters[-1][-1]['end'] <= max_gap:
                clusters[-1].append(item)
            else:
                clusters.append([item])
        all_clusters.extend(clusters)
    return all_clusters


def interval_dist(a_start, a_end, b_start, b_end):
    """Minimum distance between two genomic intervals.

    Returns 0 for overlapping intervals, otherwise the gap in bp.
    """
    if a_end < b_start:
        return b_start - a_end
    if b_end < a_start:
        return a_start - b_end
    return 0  # overlapping


def _trunc(text, width=45):
    """Truncate text with ellipsis if it exceeds width."""
    if len(text) > width:
        return text[:width - 2] + '..'
    return text


def _dark_matter_clusters(cdss, operon_gap):
    """Detect operon-like clusters per contig; never crosses contig boundaries."""
    by_contig = collections.defaultdict(list)
    for f in cdss:
        by_contig[f['contig']].append(f)

    all_clusters = []
    for contig_cdss in by_contig.values():
        sorted_cdss = sorted(contig_cdss, key=lambda f: f['start'])
        if not sorted_cdss:
            continue
        current = [sorted_cdss[0]]
        for i in range(1, len(sorted_cdss)):
            a, b = sorted_cdss[i - 1], sorted_cdss[i]
            gap  = b['start'] - a['end'] - 1
            if a['strand'] == b['strand'] and gap <= operon_gap:
                current.append(b)
            else:
                if len(current) >= 3:
                    all_clusters.append(current)
                current = [b]
        if len(current) >= 3:
            all_clusters.append(current)
    return all_clusters


def discover(filepath, cluster_gap=5000, operon_gap=150,
             min_weight=1, out_format='text'):
    """Run the full discovery pipeline on a GenBank file."""
    if not os.path.isfile(filepath):
        sys.exit(f"Error: file not found: {filepath}")

    try:
        features = parse_features(filepath)
    except Exception as e:
        sys.exit(f"Error parsing {filepath}: {e}")

    cdss = [f for f in features if f['type'] == 'CDS']

    if not cdss:
        print("No CDS features found.")
        sys.exit(1)

    total_cds = len(cdss)
    hyp_count = sum(1 for f in cdss if is_hypothetical(get_qual(f, 'product')))

    # Classify every CDS
    cat_hits = {cat: [] for cat in CATEGORIES}
    for f in cdss:
        hits = classify_cds(f, min_weight)
        for cat, kw, weight in hits:
            cat_hits[cat].append({
                'feature': f,
                'start':   f['start'],
                'end':     f['end'],
                'tag':     get_qual(f, 'locus_tag', '?'),
                'gene':    get_qual(f, 'gene') or '-',
                'product': get_qual(f, 'product'),
                'keyword': kw,
                'weight':  weight,
            })

    # Tag multi-category hits
    tag_cats = collections.defaultdict(set)
    for cat, items in cat_hits.items():
        for item in items:
            tag_cats[item['tag']].add(cat)
    for cat, items in cat_hits.items():
        for item in items:
            item['is_multi'] = len(tag_cats[item['tag']]) > 1

    # Contig summary statistics
    contig_stats = collections.defaultdict(lambda: {'cds': 0, 'min': float('inf'),
                                                     'max': 0, 'hits': 0})
    for f in cdss:
        cs = contig_stats[f['contig']]
        cs['cds'] += 1
        cs['min'] = min(cs['min'], f['start'])
        cs['max'] = max(cs['max'], f['end'])
    for cat_items_list in cat_hits.values():
        for item in cat_items_list:
            contig_stats[item['feature']['contig']]['hits'] += 1

    # ---- JSON output mode ----
    if out_format == 'json':
        _emit_json(filepath, total_cds, hyp_count, contig_stats,
                   cat_hits, cdss, operon_gap, cluster_gap)
        return

    # ---- TSV output mode ----
    if out_format == 'tsv':
        _emit_tsv(cat_hits)
        return

    # ============================================================
    # REPORT
    # ============================================================
    print("=" * 70)
    print("  GENOME DISCOVERY REPORT")
    print("=" * 70)
    print(f"  File           : {filepath}")
    print(f"  Total CDS      : {total_cds}")
    print(f"  Hypothetical   : {hyp_count} ({100 * hyp_count // total_cds}%)")
    print()

    # ---- Contig Summary ----
    print("-- Contig Summary --")
    for contig, cs in sorted(contig_stats.items()):
        span_kb = (cs['max'] - cs['min'] + 1) / 1000
        density = cs['hits'] / cs['cds'] * 100 if cs['cds'] else 0
        print(f"  {contig:30s}  {span_kb:8.1f} kb  {cs['cds']:5d} CDS  "
              f"{cs['hits']:3d} hits ({density:.1f}%)")
    print()

    # ---- Genomic Density Karyogram ----
    WINDOW = 50000
    print("-- Genomic Density Karyogram (50 kb windows) --")
    all_hit_positions = []
    for cat_items_list in cat_hits.values():
        for item in cat_items_list:
            all_hit_positions.append((item['feature']['contig'], item['start']))
    for contig, cs in sorted(contig_stats.items()):
        n_windows = max(1, (cs['max'] - cs['min']) // WINDOW + 1)
        bins = [0] * n_windows
        for c, pos in all_hit_positions:
            if c == contig:
                idx = min((pos - cs['min']) // WINDOW, n_windows - 1)
                bins[idx] += 1
        max_bin = max(bins) if bins else 1
        bar = ''.join('#' if b > max_bin * 0.66
                      else '+' if b > max_bin * 0.33
                      else '.' if b > 0
                      else ' ' for b in bins)
        print(f"  {contig:20s} |{bar}| ({n_windows} windows)")
    print()

    # ---- Summary ----
    for cat in CATEGORIES:
        n    = len(cat_hits[cat])
        high = sum(1 for h in cat_hits[cat] if h['weight'] >= 3)
        bar  = "#" * min(n, 40)
        print(f"  {cat:42s}  {n:3d} hits ({high} high-confidence)")
        if bar:
            print(f"  {'':42s}  {bar}")
    print()

    for cat in CATEGORIES:
        items = cat_hits[cat]
        if not items:
            print(f"-- {cat} --")
            print("  (no hits)")
            print()
            continue

        clusters = spatial_cluster(items, cluster_gap)
        isolated = [c for c in clusters if len(c) == 1]
        grouped  = [c for c in clusters if len(c) >= 2]

        print(f"-- {cat} ({len(items)} hits, {len(grouped)} clusters) --")
        print()

        for ci, cluster in enumerate(grouped, 1):
            span_start = cluster[0]['start']
            span_end   = cluster[-1]['end']
            span_kb    = (span_end - span_start + 1) / 1000
            contig     = cluster[0]['feature']['contig']
            print(f"  ** Cluster {ci}: {len(cluster)} genes spanning {span_kb:.1f} kb "
                  f"({span_start:,}..{span_end:,}, contig: {contig}) **")
            for h in cluster:
                conf = "***" if h['weight'] >= 3 else " * " if h['weight'] >= 2 else " . "
                strand_sym = ">" if h['feature']['strand'] == '+' else "<"
                multi_flag = "*" if h.get('is_multi') else " "
                print(f"    {conf}{multi_flag} {h['tag']:18s} {strand_sym} {h['gene']:8s}  "
                      f"{_trunc(h['product']):45s}  [{h['keyword']}]")
            print()

        if isolated:
            print(f"  Isolated hits ({len(isolated)}):")
            for cl in isolated:
                h    = cl[0]
                conf = "***" if h['weight'] >= 3 else " * " if h['weight'] >= 2 else " . "
                strand_sym = ">" if h['feature']['strand'] == '+' else "<"
                multi_flag = "*" if h.get('is_multi') else " "
                print(f"    {conf}{multi_flag} {h['tag']:18s} {strand_sym} {h['gene']:8s}  "
                      f"{_trunc(h['product']):45s}  [{h['keyword']}]")
            print()

    # ============================================================
    # DARK MATTER OPERONS
    # ============================================================
    print("=" * 70)
    print("  DARK MATTER OPERONS")
    print("=" * 70)
    print()

    all_clusters = _dark_matter_clusters(cdss, operon_gap)

    pure_dark  = []
    mixed_dark = []
    for cl in all_clusters:
        products  = [get_qual(f, 'product') for f in cl]
        hyp_frac  = sum(1 for p in products if is_hypothetical(p)) / len(products)
        if hyp_frac == 1.0:
            pure_dark.append(cl)
        elif hyp_frac >= 0.6 and len(cl) >= 3:
            mixed_dark.append((cl, hyp_frac))

    if pure_dark:
        print(f"  -- Pure dark matter ({len(pure_dark)} clusters, ALL genes unknown) --")
        for i, cl in enumerate(pure_dark, 1):
            span = f"{cl[0]['start']:,}..{cl[-1]['end']:,}"
            print(f"  Cluster {i}: {len(cl)} genes, {cl[0]['strand']} strand, "
                  f"span {span}, contig: {cl[0]['contig']}")
            for f in cl:
                tag  = get_qual(f, 'locus_tag', '?')
                prod = get_qual(f, 'product') or '(hypothetical)'
                strand_sym = ">" if f['strand'] == '+' else "<"
                print(f"    {strand_sym} {tag:18s}  {_trunc(prod, 55)}")
            print()
    else:
        print("  -- Pure dark matter --\n  (none found)\n")

    if mixed_dark:
        print(f"  -- Mixed dark matter ({len(mixed_dark)} clusters, >=60% unknown + anchor genes) --")
        for i, (cl, frac) in enumerate(mixed_dark, 1):
            span = f"{cl[0]['start']:,}..{cl[-1]['end']:,}"
            print(f"  Cluster {i}: {len(cl)} genes ({frac:.0%} hypothetical), "
                  f"{cl[0]['strand']} strand, span {span}, contig: {cl[0]['contig']}")
            for f in cl:
                tag    = get_qual(f, 'locus_tag', '?')
                prod   = get_qual(f, 'product')
                marker = "  ???" if is_hypothetical(prod) else "  <--"
                strand_sym = ">" if f['strand'] == '+' else "<"
                print(f"    {strand_sym} {tag:18s}  {_trunc(prod or '(hypothetical)', 50):50s} {marker}")
            print()
    else:
        print("  -- Mixed dark matter --\n  (none found)\n")

    # ============================================================
    # CROSS-ANALYSIS: Mobilome-adjacent functional islands
    # ============================================================
    mob_items  = cat_hits.get("Mobilome & HGT Markers", [])
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
        print("  (Functional genes within 10 kb of a mobilome marker, same contig)")
        print()

        found_islands = False
        for mob in mob_items:
            mob_contig = mob['feature']['contig']
            nearby = [
                fi for fi in func_items
                if fi['feature']['contig'] == mob_contig
                and interval_dist(mob['start'], mob['end'],
                                  fi['start'], fi['end']) <= 10000
            ]
            if nearby:
                found_islands = True
                print(f"  Near {mob['tag']} ({mob['gene']}, {mob['keyword']}, "
                      f"contig: {mob_contig}):")
                for fi in nearby:
                    dist = interval_dist(mob['start'], mob['end'],
                                        fi['start'], fi['end'])
                    kb_dist = round(dist / 1000)
                    direction = f"+{kb_dist}kb" if fi['start'] >= mob['start'] else f"-{kb_dist}kb"
                    strand_sym = ">" if fi['feature']['strand'] == '+' else "<"
                    print(f"    {direction:>6s}  {strand_sym} {fi['tag']:18s}  {_trunc(fi['product'])}")
                print()

        if not found_islands:
            print("  (No functional genes found within 10 kb of mobilome markers)\n")

    print("=" * 70)
    print("  END OF DISCOVERY REPORT")
    print("=" * 70)


def _emit_json(filepath, total_cds, hyp_count, contig_stats,
               cat_hits, cdss, operon_gap, cluster_gap):
    """Emit full discovery report as structured JSON."""
    result = {
        'file': filepath,
        'total_cds': total_cds,
        'hypothetical_count': hyp_count,
        'contig_summary': {},
        'categories': {},
        'dark_matter': {'pure': [], 'mixed': []},
        'mobilome_islands': [],
    }

    for contig, cs in sorted(contig_stats.items()):
        result['contig_summary'][contig] = {
            'span_kb': round((cs['max'] - cs['min'] + 1) / 1000, 1),
            'cds': cs['cds'], 'hits': cs['hits'],
        }

    for cat in CATEGORIES:
        items = cat_hits[cat]
        clusters = spatial_cluster(items, cluster_gap)
        result['categories'][cat] = {
            'total_hits': len(items),
            'high_confidence': sum(1 for h in items if h['weight'] >= 3),
            'clusters': [
                {'genes': [{'tag': h['tag'], 'gene': h['gene'],
                            'product': h['product'], 'keyword': h['keyword'],
                            'weight': h['weight'], 'strand': h['feature']['strand'],
                            'start': h['start'], 'end': h['end'],
                            'contig': h['feature']['contig'],
                            'is_multi': h.get('is_multi', False)}
                           for h in c],
                 'span_kb': round((c[-1]['end'] - c[0]['start'] + 1) / 1000, 1)}
                for c in clusters if len(c) >= 2
            ],
            'isolated': [
                {'tag': c[0]['tag'], 'gene': c[0]['gene'],
                 'product': c[0]['product'], 'keyword': c[0]['keyword'],
                 'weight': c[0]['weight'], 'strand': c[0]['feature']['strand'],
                 'start': c[0]['start'], 'end': c[0]['end'],
                 'contig': c[0]['feature']['contig'],
                 'is_multi': c[0].get('is_multi', False)}
                for c in clusters if len(c) == 1
            ],
        }

    # Dark matter
    all_clusters = _dark_matter_clusters(cdss, operon_gap)
    for cl in all_clusters:
        products = [get_qual(f, 'product') for f in cl]
        hyp_frac = sum(1 for p in products if is_hypothetical(p)) / len(products)
        genes = [{'tag': get_qual(f, 'locus_tag', '?'),
                  'product': get_qual(f, 'product'),
                  'strand': f['strand'], 'start': f['start'], 'end': f['end'],
                  'contig': f['contig']} for f in cl]
        if hyp_frac == 1.0:
            result['dark_matter']['pure'].append({'genes': genes, 'hyp_frac': 1.0})
        elif hyp_frac >= 0.6 and len(cl) >= 3:
            result['dark_matter']['mixed'].append({'genes': genes, 'hyp_frac': round(hyp_frac, 2)})

    # Mobilome islands
    mob_items = cat_hits.get("Mobilome & HGT Markers", [])
    func_items = []
    for cat in ["Xenobiotic Degradation & Bioremediation",
                "Heavy Metal & Stress Resistance",
                "Secondary Metabolites & BGCs"]:
        func_items.extend(cat_hits.get(cat, []))
    for mob in mob_items:
        mob_contig = mob['feature']['contig']
        nearby = [fi for fi in func_items
                  if fi['feature']['contig'] == mob_contig
                  and interval_dist(mob['start'], mob['end'],
                                    fi['start'], fi['end']) <= 10000]
        if nearby:
            result['mobilome_islands'].append({
                'mobilome': {'tag': mob['tag'], 'gene': mob['gene'],
                             'keyword': mob['keyword'], 'contig': mob_contig},
                'nearby': [{'tag': fi['tag'], 'product': fi['product'],
                            'dist_bp': interval_dist(mob['start'], mob['end'],
                                                     fi['start'], fi['end'])}
                           for fi in nearby],
            })

    print(json.dumps(result, indent=2, ensure_ascii=False))


def _emit_tsv(cat_hits):
    """Emit flat TSV table of all category hits."""
    print('\t'.join(['contig', 'start', 'end', 'strand', 'locus_tag', 'gene',
                     'product', 'category', 'keyword', 'weight', 'is_multi']))
    for cat in CATEGORIES:
        for h in sorted(cat_hits[cat], key=lambda x: (x['feature']['contig'], x['start'])):
            print('\t'.join([
                h['feature']['contig'], str(h['start']), str(h['end']),
                h['feature']['strand'], h['tag'], h['gene'], h['product'],
                cat, h['keyword'], str(h['weight']),
                str(h.get('is_multi', False)),
            ]))


def main():
    parser = argparse.ArgumentParser(
        description="Discovery mode: keyword scan, spatial clustering, and dark-matter operon detection."
    )
    parser.add_argument('input', help="Input GenBank file")
    parser.add_argument('--cluster-gap', type=int, default=5000,
                        help="Max intergenic gap for spatial clustering (default: 5000 bp)")
    parser.add_argument('--operon-gap',  type=int, default=150,
                        help="Max intergenic gap for operon clustering (default: 150 bp)")
    parser.add_argument('--min-weight', type=int, default=1, choices=[1, 2, 3],
                        help="Minimum keyword weight to report (default: 1)")
    parser.add_argument('--format', choices=['text', 'json', 'tsv'], default='text',
                        help="Output format (default: text)")
    args = parser.parse_args()
    discover(args.input, args.cluster_gap, args.operon_gap,
             args.min_weight, args.format)


if __name__ == '__main__':
    main()
