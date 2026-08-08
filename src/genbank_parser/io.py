"""Input/output handling and canonical GenBank parsing using Biopython."""
from __future__ import annotations

import collections
from pathlib import Path
import re
from typing import Any, Sequence

from Bio import SeqIO
from Bio.Seq import Seq

from .model import GenBankDocument, GenBankFeature, GenBankRecord

# ---------------------------------------------------------------------------
# Cross-reference regexes and prefix registries
# ---------------------------------------------------------------------------
_known_xref_prefixes = ('GO:', 'COG', 'PF', 'RFAM', 'Rfam:', 'EC:')
_kegg_ko_re = re.compile(r'^(?:KEGG:)?(K\d{5})')
_cog_re = re.compile(r'^(?:COG:)?(COG\d+)')
_pfam_re = re.compile(r'^(?:Pfam:)?(PF\d+)')
_rfam_re = re.compile(r'^(?:Rfam:)?(RF\d+)', re.IGNORECASE)


def read_genbank(filepath: str | Path) -> GenBankDocument:
    """Read a GenBank flatfile into a fully typed GenBankDocument."""
    path = Path(filepath)
    records: list[GenBankRecord] = []
    global_feat_counter = 0

    with path.open('r', encoding='utf-8', errors='replace') as handle:
        for rec_idx, rec in enumerate(SeqIO.parse(handle, "genbank"), 1):
            contig = rec.id if (rec.id and rec.id != '.') else rec.name
            topology = rec.annotations.get('topology')
            mol_type = rec.annotations.get('molecule_type')
            division = rec.annotations.get('data_file_division')
            date = rec.annotations.get('date')
            seq = rec.seq if rec.seq is not None else Seq("")
            rec_len = len(seq)
            if rec_len == 0 and 'length' in rec.annotations:
                rec_len = int(rec.annotations['length'])

            features: list[GenBankFeature] = []
            for feat in rec.features:
                global_feat_counter += 1

                # Normalise qualifiers to dict[str, list[str]]
                quals: dict[str, list[str]] = collections.defaultdict(list)
                for k, v_list in feat.qualifiers.items():
                    if isinstance(v_list, list):
                        quals[k] = [str(x) for x in v_list]
                    else:
                        quals[k] = [str(v_list)]

                gb_feat = GenBankFeature(
                    record_id=contig,
                    record_index=rec_idx,
                    feature_index=global_feat_counter,
                    type=feat.type,
                    location=feat.location,
                    qualifiers=dict(quals),
                    record_length=rec_len,
                    topology=topology,
                    raw_feature=feat,
                )
                features.append(gb_feat)

            gb_rec = GenBankRecord(
                id=contig,
                name=rec.name,
                description=rec.description,
                seq=seq,
                length=rec_len,
                topology=topology,
                molecule_type=mol_type,
                division=division,
                date=date,
                annotations=dict(rec.annotations),
                features=features,
            )
            records.append(gb_rec)

    return GenBankDocument(path=path, records=records)


def parse_features(filepath: str | Path) -> list[GenBankFeature]:
    """Parse GenBank file and return all features as a flat list."""
    doc = read_genbank(filepath)
    return doc.all_features


def get_qual(feature: Any, key: str, default: str = '') -> str:
    """Return first qualifier value for key, or default."""
    if hasattr(feature, 'get_qual'):
        return feature.get_qual(key, default)
    if isinstance(feature, dict):
        vals = feature.get('qualifiers', {}).get(key, [])
        return vals[0] if vals else default
    return default


def get_notes(feature: Any, prefix: str) -> list[str]:
    """Return all /note values that start with prefix (e.g. 'COG:', 'GO:')."""
    if hasattr(feature, 'qualifiers'):
        notes = feature.qualifiers.get('note', [])
    elif isinstance(feature, dict):
        notes = feature.get('qualifiers', {}).get('note', [])
    else:
        notes = []
    return [n for n in notes if n.startswith(prefix)]


def extract_xrefs(feature: Any) -> dict[str, list[str]]:
    """Extract semantically typed cross-references from a feature.

    Searches both /db_xref and /note for prefixed identifiers, then adds
    /EC_number with higher priority than EC: in /note (INSDC vs Bakta style).

    Returns a dict:
        go_terms   : list[str]   GO:0001234
        cog_ids    : list[str]   COG1234
        kegg_kos   : list[str]   K01234
        pfam       : list[str]   PF01234
        rfam       : list[str]   RF01234
        ec_numbers : list[str]   1.2.3.4  (from /EC_number or EC: in /note)
        db_xrefs   : list[str]   everything else in /db_xref
    """
    if hasattr(feature, 'qualifiers'):
        quals = feature.qualifiers
    elif isinstance(feature, dict):
        quals = feature.get('qualifiers', {})
    else:
        quals = {}

    all_xrefs = quals.get('db_xref', []) + quals.get('note', [])

    go_terms = [x for x in all_xrefs if x.startswith('GO:')]

    cog_ids: list[str] = []
    for x in all_xrefs:
        m = _cog_re.match(x)
        if m:
            cog_ids.append(m.group(1))

    kegg_kos: list[str] = []
    for x in all_xrefs:
        m = _kegg_ko_re.match(x)
        if m:
            kegg_kos.append(m.group(1))

    pfam: list[str] = []
    for x in all_xrefs:
        m = _pfam_re.match(x)
        if m:
            pfam.append(m.group(1))

    rfam: list[str] = []
    for x in all_xrefs:
        m = _rfam_re.match(x)
        if m:
            rfam.append(m.group(1))

    # EC_number: /EC_number qualifier takes priority (INSDC submission style).
    # Fall back to EC:-prefixed values in /note or /db_xref (Bakta style).
    ec_qual = quals.get('EC_number', [])
    ec_note = [x[3:] for x in all_xrefs if x.startswith('EC:')]
    ec_numbers = ec_qual if ec_qual else ec_note

    db_xrefs = [
        x for x in quals.get('db_xref', [])
        if not any(x.startswith(p) for p in _known_xref_prefixes)
        and not _kegg_ko_re.match(x)
        and not _cog_re.match(x)
        and not _pfam_re.match(x)
        and not _rfam_re.match(x)
    ]

    return {
        'go_terms': list(dict.fromkeys(go_terms)),
        'cog_ids': list(dict.fromkeys(cog_ids)),
        'kegg_kos': list(dict.fromkeys(kegg_kos)),
        'pfam': list(dict.fromkeys(pfam)),
        'rfam': list(dict.fromkeys(rfam)),
        'ec_numbers': list(dict.fromkeys(ec_numbers)),
        'db_xrefs': list(dict.fromkeys(db_xrefs)),
    }
