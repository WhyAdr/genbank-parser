#!/usr/bin/env python3
"""Phylogenetic marker gene extractor.

Scans a GenBank file for universally conserved single-copy genes commonly
used in prokaryotic phylogenetics, then exports each marker as a separate
protein FASTA or produces a concatenated multi-marker alignment-ready file.

Marker sets supported:
  --markers core         : 30 ribosomal proteins (rpsB-L, rplA-L subsets)
  --markers housekeeping : RecA, GyrB, RpoB, DnaK, GroEL + EF series
  --markers all          : union of both sets (default)

Usage:
    python genbank_phylo.py INPUT.gbff [--markers all] [--output-dir ./phylo_markers/]

Matching priority (per CDS, per marker):
  1. /gene qualifier exact match                 (most reliable)
  2. /product qualifier exact match
  3. /product regex match (only when /gene is absent)

The regex fallback (step 3) uses a product-description pattern anchored with
a negative lookahead (?![\\w\\-]) so that e.g. 'L1' does NOT match 'L11',
'L14', 'L18' etc., and 'Elongation factor G' does NOT match 'Elongation
factor G-binding protein'.  This resolves the two known failure modes of an
unanchored .search():
  * Numeric suffix bleeding  : 'S2' into S20/S21, 'L1' into L10/L11/L14-L18
  * Hyphen-linked derivatives: 'factor G' into 'factor G-binding protein'

When /gene IS present and does not match the target, the regex fallback is
suppressed entirely -- trusting the annotator's explicit gene assignment over
a fuzzy product string match.  This eliminates false positives from genes like
dksA ('DnaK suppressor protein'), rpoC ('RNA polymerase subunit beta-prime'),
and greA ('transcription elongation factor GreA').
"""

import sys
import re
import os
import argparse
from genbank_parser import parse_features, get_qual

# ---------------------------------------------------------------------------
# Curated marker dictionaries: gene_name -> canonical product description
# ---------------------------------------------------------------------------

RIBOSOMAL = {
    # Small subunit (30S)
    'rpsB': '30S ribosomal protein S2',
    'rpsC': '30S ribosomal protein S3',
    'rpsD': '30S ribosomal protein S4',
    'rpsE': '30S ribosomal protein S5',
    'rpsG': '30S ribosomal protein S7',
    'rpsH': '30S ribosomal protein S8',
    'rpsI': '30S ribosomal protein S9',
    'rpsJ': '30S ribosomal protein S10',
    'rpsK': '30S ribosomal protein S11',
    'rpsL': '30S ribosomal protein S12',
    'rpsM': '30S ribosomal protein S13',
    # Large subunit (50S)
    'rplA': '50S ribosomal protein L1',
    'rplB': '50S ribosomal protein L2',
    'rplC': '50S ribosomal protein L3',
    'rplD': '50S ribosomal protein L4',
    'rplE': '50S ribosomal protein L5',
    'rplF': '50S ribosomal protein L6',
    'rplK': '50S ribosomal protein L11',
    'rplL': '50S ribosomal protein L12',
    'rplN': '50S ribosomal protein L14',
    'rplP': '50S ribosomal protein L16',
    'rplR': '50S ribosomal protein L18',
    'rplV': '50S ribosomal protein L22',
    'rplW': '50S ribosomal protein L23',
}

HOUSEKEEPING = {
    'recA':  'DNA recombinase A',
    'gyrB':  'DNA gyrase subunit B',
    'rpoB':  'RNA polymerase subunit beta',
    'dnaK':  'DnaK chaperone',
    'groEL': 'GroEL chaperonin (Hsp60)',
    'atpD':  'ATP synthase subunit beta',
    'fusA':  'Elongation factor G',
    'tsf':   'Elongation factor Ts',
    'tuf':   'Elongation factor Tu',
}

ALL_MARKERS = {**RIBOSOMAL, **HOUSEKEEPING}
MARKER_SETS = {
    'core':         RIBOSOMAL,
    'housekeeping': HOUSEKEEPING,
    'all':          ALL_MARKERS,
}

# ---------------------------------------------------------------------------
# Product regex patterns -- anchored with negative lookahead.
#
# (?![\w\-]) fails if the pattern is immediately followed by a word character
# OR a hyphen, preventing:
#   * 'L1'  matching 'L11', 'L14', 'L16', 'L18'  (word chars after digit)
#   * 'G'   matching 'G-binding'                  (hyphen after letter)
#   * 'S2'  matching 'S20', 'S21'                 (word chars after digit)
#
# The gene-name-as-bare-substring alternative (e.g. '|dnak') from the
# original version has been REMOVED. It caused dksA ('DnaK suppressor') to
# be pulled in as a dnaK hit. The /gene exact-match handles legitimate
# gene-name matching without any regex needed.
# ---------------------------------------------------------------------------
# The description pattern anchored with (?![\w\-]) prevents suffix bleed
# (e.g. 'L1' into 'L11/L14/L18', 'G' into 'G-binding').
#
# The \b-anchored gene-name alternative handles products that name the gene
# AFTER a modifier, e.g. 'chaperonin GroEL' or 'co-chaperonin GroES'.
# This alternative is ONLY exercised when /gene is absent (guarded in
# _matches_marker), so it cannot reach proteins like dksA that do carry
# a /gene qualifier.
_PRODUCT_PATTERNS = {
    marker: re.compile(
        re.escape(desc.lower()) + r'(?![\w\-])'
        + r'|' +
        r'\b' + re.escape(marker.lower()) + r'\b',
        re.IGNORECASE
    )
    for marker, desc in ALL_MARKERS.items()
}

# Minimum translation length (aa) -- shorter sequences are likely pseudogene
# fragments or assembly frameshifts, not suitable for phylogenetic inference.
_MIN_AA_LENGTH = 50


def _matches_marker(f, marker, desc):
    """Return True if feature f is a credible hit for (marker, desc).

    Matching priority:
      1. /gene exact match            -- most authoritative
      2. /product exact match
      3. /product regex match         -- ONLY when /gene is absent

    When /gene IS present but does not equal the target marker, the regex
    fallback is suppressed to avoid grabbing proteins that merely contain
    a marker-related substring in their product description (e.g. dksA,
    rpoC, greA, rpmD, rplO ...).
    """
    gene    = get_qual(f, 'gene').lower()
    product = get_qual(f, 'product').lower()

    # Priority 1: explicit gene annotation
    if gene == marker.lower():
        return True

    # Priority 2: exact product match
    if product == desc.lower():
        return True

    # Priority 3: fuzzy product regex -- only when /gene is absent.
    # If a gene IS annotated (even a different one), trust that over any
    # substring match in the product description.
    if not gene:
        return bool(_PRODUCT_PATTERNS[marker].search(product))

    return False


def _product_word_overlap(desc, product):
    """Count shared significant words (>2 chars) between desc and product.
    Low overlap (<1) indicates a genuine gene-product mismatch worth flagging.
    """
    desc_words = set(re.findall(r'\b\w{3,}\b', desc.lower()))
    prod_words = set(re.findall(r'\b\w{3,}\b', product.lower()))
    return len(desc_words & prod_words)


def _warn_product_mismatch(marker, desc, tag, product):
    """Return a warning string for gene-match hits with unexpected products."""
    return (
        f"  [WARN] {tag}: /gene='{marker}' matched but product is unexpected.\n"
        f"         Expected : '{desc}'\n"
        f"         Observed : '{product}'\n"
        f"         Possible Bakta/Prokka misannotation -- verify manually.\n"
    )


def extract_phylo_markers(filepath, marker_set='all', output_dir=None,
                          min_length=None):
    if min_length is None:
        min_length = _MIN_AA_LENGTH

    markers  = MARKER_SETS.get(marker_set, ALL_MARKERS)
    features = parse_features(filepath)
    cdss     = [f for f in features if f['type'] == 'CDS']

    # Scan -- one independent pass per marker over all CDSs
    found    = {}   # marker -> list of (locus_tag, translation, contig)
    warnings = []   # gene-product mismatch warnings

    for marker, desc in markers.items():
        hits = []
        for f in cdss:
            if not _matches_marker(f, marker, desc):
                continue
            trans = get_qual(f, 'translation')
            tag   = get_qual(f, 'locus_tag', '?')
            if not trans or len(trans) < min_length:
                continue
            # Flag gene-product mismatches for manual review.
            # Use word-overlap heuristic: if the expected description shares
            # fewer than 1 significant word (>2 chars) with the observed product,
            # this is likely a genuine mismatch (e.g. gene=fusA but product is
            # 'Tetracycline resistance protein').  Legitimate alternative names
            # like 'recombinase RecA', 'chaperone protein DnaK', or
            # '50S ribosomal protein L7/L12' all share >=1 key word and are
            # therefore not flagged.
            product = get_qual(f, 'product', '').lower()
            if product and _product_word_overlap(desc, product) < 1:
                warnings.append(_warn_product_mismatch(marker, desc, tag, product))
            hits.append((tag, trans, f['contig']))
        found[marker] = hits

    # -----------------------------------------------------------------------
    # Report
    # -----------------------------------------------------------------------
    print("=" * 65)
    print("  PHYLOGENETIC MARKER EXTRACTION REPORT")
    print("=" * 65)
    print(f"  File        : {filepath}")
    print(f"  Marker set  : {marker_set} ({len(markers)} markers)")
    print()

    recovered   = sum(1 for v in found.values() if v)
    single_copy = {m: v for m, v in found.items() if len(v) == 1}
    multi_copy  = {m: v for m, v in found.items() if len(v) > 1}
    absent_list = [m for m, v in found.items() if not v]

    print(f"  Recovered   : {recovered}/{len(markers)} markers")
    print(f"  Single-copy : {len(single_copy)}")
    print(f"  Multi-copy  : {len(multi_copy)}  (require manual verification)")
    print(f"  Absent      : {len(absent_list)}")
    print()

    cw = (8, 10, 18, 20)
    print(f"  {'Marker':{cw[0]}s}  {'Status':{cw[1]}s}  "
          f"{'Locus Tag':{cw[2]}s}  {'Contig':{cw[3]}s}  Description")
    print("  " + "-" * 88)

    for marker, hits in found.items():
        if hits:
            flag = "  *** PARALOGS -- verify ***" if len(hits) > 1 else ""
            for tag, _, contig in hits:
                print(f"  {marker:{cw[0]}s}  {'FOUND':{cw[1]}s}  "
                      f"{tag:{cw[2]}s}  {contig:{cw[3]}s}  "
                      f"{markers[marker]}{flag}")
        else:
            print(f"  {marker:{cw[0]}s}  {'ABSENT':{cw[1]}s}  "
                  f"{'':18s}  {'':20s}  {markers[marker]}")

    if warnings:
        print()
        print("  GENE-PRODUCT MISMATCH WARNINGS")
        print("  " + "-" * 55)
        for w in warnings:
            print(w, end='')

    # -----------------------------------------------------------------------
    # Write per-marker FASTA files
    # -----------------------------------------------------------------------
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        for marker, hits in found.items():
            if not hits:
                continue
            fpath = os.path.join(output_dir, f"{marker}.faa")
            with open(fpath, 'w') as fh:
                for tag, trans, contig in hits:
                    fh.write(f">{tag} [{marker}] [{contig}]\n")
                    for i in range(0, len(trans), 60):
                        fh.write(trans[i:i + 60] + '\n')

        print()
        print(f"  Per-marker FASTA files written to: {output_dir}")

        # Concatenated alignment skeleton -- single-copy markers ONLY
        sc_hits = {m: v[0] for m, v in found.items() if len(v) == 1}
        if sc_hits:
            concat_path = os.path.join(output_dir, 'concatenated_markers.faa')
            with open(concat_path, 'w') as fh:
                genome_name = os.path.basename(filepath).rsplit('.', 1)[0]
                concat_seq  = ''.join(trans for _, trans, _ in sc_hits.values())
                fh.write(f">{genome_name}\n")
                for i in range(0, len(concat_seq), 60):
                    fh.write(concat_seq[i:i + 60] + '\n')
            print(f"  Concatenated ({len(sc_hits)} single-copy markers): {concat_path}")
    else:
        for marker, hits in found.items():
            for tag, trans, contig in hits:
                print(f">{tag} [{marker}] [{contig}]")
                for i in range(0, len(trans), 60):
                    print(trans[i:i + 60])


def main():
    parser = argparse.ArgumentParser(
        description="Extract phylogenetic marker genes from a GenBank file.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('input',
                        help="Input GenBank (.gbff/.gbk) file")
    parser.add_argument('--markers',
                        choices=['core', 'housekeeping', 'all'],
                        default='all',
                        help="Marker set to use (default: all)")
    parser.add_argument('--output-dir',
                        help="Directory for per-marker FASTA files and "
                             "concatenated alignment skeleton")
    parser.add_argument('--min-length',
                        type=int, default=_MIN_AA_LENGTH,
                        metavar='AA',
                        help=f"Minimum translation length in amino acids "
                             f"(default: {_MIN_AA_LENGTH})")
    args = parser.parse_args()
    extract_phylo_markers(args.input, args.markers, args.output_dir,
                          args.min_length)


if __name__ == '__main__':
    main()
