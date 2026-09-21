"""Input/output handling and canonical GenBank parsing using Biopython."""
from __future__ import annotations

import collections
import gzip
import hashlib
import io
import os
import re
import sys
from collections.abc import Iterator
from pathlib import Path
from typing import Any, TextIO

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


class GenBankInputError(ValueError):
    """Raised when a GenBank source cannot be decoded or parsed."""


def _read_source_text(source: str | Path | TextIO) -> tuple[str, Path | None, str, str, str]:
    """Read a source once and return text plus source metadata."""

    if hasattr(source, "read"):
        stream = source  # type: ignore[assignment]
        label = str(getattr(stream, "name", "<stream>"))
        raw_stream = getattr(stream, "buffer", None)
        try:
            raw = raw_stream.read() if raw_stream is not None else stream.read()
        except (OSError, UnicodeError) as exc:
            raise GenBankInputError(f"could not read input stream {label}: {exc}") from exc
        source_kind = "stdin" if label in {"<stdin>", "-"} else "stream"
        path = None
    else:
        source_text = str(source)
        if source_text == "-":
            stream = getattr(sys.stdin, "buffer", sys.stdin)
            try:
                raw = stream.read()
            except (OSError, UnicodeError) as exc:
                raise GenBankInputError(f"could not read stdin: {exc}") from exc
            label = "-"
            source_kind = "stdin"
            path = None
        else:
            path = Path(source)
            label = str(path)
            source_kind = "path"
            try:
                raw = path.read_bytes()
            except OSError as exc:
                raise GenBankInputError(f"could not read input {path}: {exc}") from exc

    if isinstance(raw, str):
        raw_bytes = raw.encode("utf-8")
    elif isinstance(raw, (bytes, bytearray)):
        raw_bytes = bytes(raw)
    else:
        raise GenBankInputError(f"input {label} did not produce text or bytes")
    raw_digest = hashlib.sha256(raw_bytes).hexdigest()
    expected_digest = os.environ.get("GBPARSE_SOURCE_SHA256")
    logical_label = os.environ.get("GBPARSE_SOURCE_LABEL")
    if expected_digest is not None or logical_label is not None:
        if not expected_digest or not logical_label:
            raise GenBankInputError(
                "GBPARSE_SOURCE_SHA256 and GBPARSE_SOURCE_LABEL must be provided together"
            )
        if expected_digest.casefold() != raw_digest:
            raise GenBankInputError(
                f"source snapshot digest mismatch for {label}: expected {expected_digest}, got {raw_digest}"
            )
        label = logical_label
    if raw_bytes.startswith(b"\x1f\x8b"):
        try:
            raw_bytes = gzip.decompress(raw_bytes)
        except (OSError, EOFError) as exc:
            raise GenBankInputError(f"truncated or invalid gzip input: {label}") from exc
    return raw_bytes.decode("utf-8", errors="replace"), path, label, source_kind, raw_digest


def iter_genbank(source: str | Path | TextIO) -> Iterator[GenBankRecord]:
    """Yield typed records from a GenBank path, gzip path, stdin, or stream."""

    text, _path, _label, _source_kind, _raw_digest = _read_source_text(source)
    global_feature_index = 0
    for rec_idx, rec in enumerate(SeqIO.parse(io.StringIO(text), "genbank"), 1):
        contig = rec.id if (rec.id and rec.id != ".") else rec.name
        topology = rec.annotations.get("topology")
        mol_type = rec.annotations.get("molecule_type")
        division = rec.annotations.get("data_file_division")
        date = rec.annotations.get("date")
        seq = rec.seq if rec.seq is not None else Seq("")
        rec_len = len(seq)
        if rec_len == 0 and "length" in rec.annotations:
            try:
                rec_len = int(rec.annotations["length"])
            except (TypeError, ValueError):
                rec_len = 0

        features: list[GenBankFeature] = []
        for feat in rec.features:
            global_feature_index += 1
            quals: dict[str, list[str]] = collections.defaultdict(list)
            for key, values in feat.qualifiers.items():
                if isinstance(values, list):
                    quals[key] = [str(value) for value in values]
                else:
                    quals[key] = [str(values)]
            features.append(
                GenBankFeature(
                    record_id=contig,
                    record_index=rec_idx,
                    feature_index=global_feature_index,
                    type=feat.type,
                    location=feat.location,
                    qualifiers=dict(quals),
                    record_length=rec_len,
                    topology=topology,
                    raw_feature=feat,
                )
            )
        yield GenBankRecord(
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
            raw_record=rec,
            record_index=rec_idx,
        )


def read_genbank(source: str | Path | TextIO) -> GenBankDocument:
    """Read a GenBank flatfile into a fully typed GenBankDocument.

    The source may be a path, a gzip-compressed path regardless of suffix,
    ``-`` for stdin, or a text/binary-backed stream.  The returned document
    remains the canonical in-memory parser model used by every analyzer.
    """

    text, path, label, source_kind, _raw_digest = _read_source_text(source)
    records: list[GenBankRecord] = []
    try:
        records.extend(iter_genbank(io.StringIO(text)))
    except (AttributeError, IndexError, KeyError, OSError, TypeError, ValueError) as exc:
        raise GenBankInputError(f"could not parse GenBank input {label}: {exc}") from exc
    # Feature indices historically were global across a document.  Keep that
    # compatibility even though iter_genbank exposes record-local indices.
    global_index = 0
    for record in records:
        for feature in record.features:
            global_index += 1
            feature.feature_index = global_index
    return GenBankDocument(
        path=path,
        records=records,
        source_label=label,
        source_kind=source_kind,
    )


def parse_features(filepath: str | Path | TextIO) -> list[GenBankFeature]:
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


def extract_xrefs(
    feature: Any,
    *,
    include_notes: bool = True,
) -> dict[str, list[str]]:
    """Extract semantically typed cross-references from a feature.

    Searches /db_xref and, by default, /note for prefixed identifiers, then
    adds /EC_number with higher priority than EC-prefixed values. Set
    ``include_notes=False`` when free-text note evidence must remain distinct
    from structured qualifiers.

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

    all_xrefs = list(quals.get('db_xref', []))
    if include_notes:
        all_xrefs.extend(quals.get('note', []))

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


def extract_xref_sources(feature: Any) -> dict[str, list[dict[str, str]]]:
    """Return typed xrefs with the qualifier field that supplied each value."""

    if hasattr(feature, "qualifiers"):
        qualifiers = feature.qualifiers
    elif isinstance(feature, dict):
        qualifiers = feature.get("qualifiers", {})
    else:
        qualifiers = {}
    typed = extract_xrefs(feature, include_notes=True)
    result: dict[str, list[dict[str, str]]] = {key: [] for key in typed}
    fields = ("db_xref", "note", "EC_number")
    raw_values = [(field, str(value)) for field in fields for value in qualifiers.get(field, [])]
    for category, values in typed.items():
        for value in values:
            source_field = "EC_number" if category == "ec_numbers" and value in qualifiers.get("EC_number", []) else "db_xref"
            if any(raw == value and field == "note" for field, raw in raw_values):
                source_field = "note"
            result[category].append({"value": value, "source_field": source_field})
    return result


__all__ = [
    "GenBankInputError",
    "extract_xref_sources",
    "extract_xrefs",
    "get_notes",
    "get_qual",
    "iter_genbank",
    "parse_features",
    "read_genbank",
]
