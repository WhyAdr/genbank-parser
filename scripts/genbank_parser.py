#!/usr/bin/env python3
"""
Canonical GenBank feature table parser using Biopython (Bio.SeqIO).
Single source of truth -- all analysis scripts import from here.

Handles:
  - Multi-contig/multi-record files (LOCUS splitting via SeqRecord)
  - complement() and strand handling
  - Partial features (<start, >end) via BeforePosition/AfterPosition
  - Multi-line qualifier continuation with proper spacing and quote stripping
  - Simple and complex join()/order() compound locations
  - All feature keys parsed structurally
"""
import collections
import re
import sys
from Bio import SeqIO
from Bio.SeqFeature import BeforePosition, AfterPosition, CompoundLocation


# ---------------------------------------------------------------------------
# Compiled regexes for cross-reference extraction
# ---------------------------------------------------------------------------
_known_xref_prefixes = ('GO:', 'COG', 'PF', 'RFAM', 'EC:')
_kegg_ko_re = re.compile(r'^(?:KEGG:)?(K\d{5})')
_cog_re = re.compile(r'^(?:COG:)?(COG\d+)')
_pfam_re = re.compile(r'^(?:Pfam:)?(PF\d+)')


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def parse_features(filepath):
    """Parse a GenBank-style feature table into structured records via Biopython.

    Returns a list of dicts, each with keys:
        type (str), start (int), end (int), strand ('+'/'-'),
        contig (str), qualifiers (defaultdict(list)), line (int),
        partial_start (bool), partial_end (bool),
        join_segments (list[tuple[int,int]])

    Each feature includes a 'contig' key from its parent LOCUS record,
    preventing cross-contig clustering in downstream analyses.
    """
    features = []
    feat_counter = 0
    for record in SeqIO.parse(filepath, "genbank"):
        contig = record.id if (record.id and record.id != '.') else record.name
        for feat in record.features:
            feat_counter += 1
            strand = '-' if feat.location.strand == -1 else '+'

            # Partiality detection
            start_pos = feat.location.start
            end_pos = feat.location.end
            partial_start = isinstance(start_pos, BeforePosition) or '<' in str(start_pos)
            partial_end = isinstance(end_pos, AfterPosition) or '>' in str(end_pos)

            # Qualifiers: Biopython returns dict[str, list[str]]
            quals = collections.defaultdict(list)
            for k, v_list in feat.qualifiers.items():
                if isinstance(v_list, list):
                    quals[k] = list(v_list)
                else:
                    quals[k] = [str(v_list)]

            # Compound location segments (join/order)
            segments = []
            if isinstance(feat.location, CompoundLocation):
                for part in feat.location.parts:
                    segments.append((int(part.start) + 1, int(part.end)))

            features.append({
                'type':          feat.type,
                'start':         int(start_pos) + 1,  # Biopython 0-based -> 1-based
                'end':           int(end_pos),
                'strand':        strand,
                'contig':        contig,
                'qualifiers':    quals,
                'line':          feat_counter,
                'partial_start': partial_start,
                'partial_end':   partial_end,
                'join_segments': segments,
            })
    return features


def get_qual(feature, key, default=''):
    """Return first qualifier value for *key*, or *default*."""
    vals = feature['qualifiers'].get(key, [])
    return vals[0] if vals else default


def get_notes(feature, prefix):
    """Return all /note values that start with *prefix* (e.g. 'COG:', 'GO:')."""
    return [n for n in feature['qualifiers'].get('note', []) if n.startswith(prefix)]


def extract_xrefs(feature):
    """Extract semantically typed cross-references from a feature.

    Searches both /db_xref and /note for prefixed identifiers, then adds
    /EC_number with higher priority than EC: in /note (INSDC vs Bakta style).

    Returns a dict:
        go_terms   : list[str]   GO:0001234
        cog_ids    : list[str]   COG1234
        kegg_kos   : list[str]   K01234
        pfam       : list[str]   PF01234
        rfam       : list[str]   RFAM:RF01234
        ec_numbers : list[str]   1.2.3.4  (from /EC_number or EC: in /note)
        db_xrefs   : list[str]   everything else in /db_xref
    """
    all_xrefs = (
        feature['qualifiers'].get('db_xref', []) +
        feature['qualifiers'].get('note', [])
    )

    go_terms = [x for x in all_xrefs if x.startswith('GO:')]
    cog_ids  = []
    for x in all_xrefs:
        m = _cog_re.match(x)
        if m:
            cog_ids.append(m.group(1))
    kegg_kos = []
    for x in all_xrefs:
        m = _kegg_ko_re.match(x)
        if m:
            kegg_kos.append(m.group(1))
    pfam     = []
    for x in all_xrefs:
        m = _pfam_re.match(x)
        if m:
            pfam.append(m.group(1))
    rfam     = [x for x in all_xrefs if x.startswith('RFAM')]

    # EC_number: /EC_number qualifier takes priority (INSDC submission style).
    # Fall back to EC:-prefixed values in /note or /db_xref (Bakta style).
    ec_qual    = feature['qualifiers'].get('EC_number', [])
    ec_note    = [x[3:] for x in all_xrefs if x.startswith('EC:')]
    ec_numbers = ec_qual if ec_qual else ec_note

    db_xrefs = [
        x for x in feature['qualifiers'].get('db_xref', [])
        if not any(x.startswith(p) for p in _known_xref_prefixes)
        and not _kegg_ko_re.match(x)
        and not _cog_re.match(x)
        and not _pfam_re.match(x)
    ]

    return {
        'go_terms':   go_terms,
        'cog_ids':    cog_ids,
        'kegg_kos':   kegg_kos,
        'pfam':       pfam,
        'rfam':       rfam,
        'ec_numbers': ec_numbers,
        'db_xrefs':   db_xrefs,
    }
