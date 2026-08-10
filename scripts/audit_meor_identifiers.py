#!/usr/bin/env python3
"""Audit the MEOR identifier catalog against offline reference snapshots.

The audit is deliberately read-only.  It validates the identifiers recorded in
``markers.yaml``, resolves them against the supplied KEGG, COG, ExplorEnz (or
ExPASy/TSV), and Pfam snapshots, and writes deterministic TSV/Markdown reports.
It does not infer new mappings or rewrite the catalog.

Example::

    python scripts/audit_meor_identifiers.py \
        --markers src/genbank_parser/data/meor/markers.yaml \
        --kegg audit-reference/kegg_ko_list_2026-08-02.tsv \
        --cog audit-reference/cog-24.def.tab \
        --ec-xml audit-reference/explorenz-enzyme-data.xml \
        --pfam audit-reference/Pfam-A.hmm \
        --outdir audit
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from collections.abc import Iterable
from dataclasses import dataclass
from html import unescape
from pathlib import Path
from typing import Any

import yaml

KO_RE = re.compile(r"K\d{5}\Z")
EC_RE = re.compile(r"\d+\.(?:\d+|-)\.(?:\d+|-)\.(?:\d+|-)\Z")
COG_RE = re.compile(r"COG\d+\Z")
PFAM_RE = re.compile(r"PF\d{5}(?:\.\d+)?\Z")
EC_BRACKET_RE = re.compile(r"\[EC:([^\]]+)\]")
EC_VALUE_RE = re.compile(r"\d+\.\d+\.\d+(?:\.\d+|\.-)\Z")
SUCCESSOR_RE = re.compile(r"ec=(\d+\.\d+\.\d+\.\d+)")
HTML_TAG_RE = re.compile(r"<[^>]*>")
GENERIC_TERMS_RE = re.compile(
    r"\b(?:hypothetical|uncharacterized|unknown function|general function prediction only|DUF)\b",
    re.IGNORECASE,
)

AUDIT_HEADER = [
    "marker_id",
    "identifier_type",
    "identifier",
    "reference_status",
    "reference_name",
    "catalog_interpretation",
    "audit_status",
    "notes",
]
STATUSES = {
    "PASS",
    "MISSING_FROM_REFERENCE",
    "DEFINITION_MISMATCH",
    "NONSPECIFIC",
    "DUPLICATED",
    "MALFORMED",
    "REVIEW",
    "REMOVE",
    "REPLACE",
}


class AuditInputError(ValueError):
    """Raised when an audit input is missing or has an unusable structure."""


@dataclass
class Finding:
    marker_id: str
    identifier_type: str
    identifier: str
    reference_status: str
    reference_name: str
    catalog_interpretation: str
    audit_status: str
    notes: str = ""

    def row(self) -> list[str]:
        return [
            self.marker_id,
            self.identifier_type,
            self.identifier,
            self.reference_status,
            self.reference_name,
            self.catalog_interpretation,
            self.audit_status,
            self.notes or "-",
        ]


@dataclass
class ReferenceData:
    kegg: dict[str, dict[str, Any]]
    cogs: dict[str, dict[str, str]]
    ec_entries: dict[str, dict[str, str]]
    ec_classes: dict[tuple[str, str, str], str]
    ec_history: dict[str, dict[str, str]]
    pfams: dict[str, dict[str, str]]
    ec_source: str


def _require_file(path: Path, label: str) -> Path:
    if not path.is_file():
        raise AuditInputError(f"{label} not found: {path}")
    return path


def load_markers(path: Path) -> list[dict[str, Any]]:
    _require_file(path, "markers file")
    try:
        payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise AuditInputError(f"could not read markers file {path}: {exc}") from exc
    if not isinstance(payload, dict) or not isinstance(payload.get("markers"), list):
        raise AuditInputError("markers file must contain a top-level 'markers' list")
    return [
        item if isinstance(item, dict) else {"_malformed_marker": item}
        for item in payload["markers"]
    ]


def _string_list(marker: dict[str, Any], key: str) -> tuple[list[str], list[str]]:
    """Return scalar values and human-readable malformed values for a field."""
    value = marker.get(key, [])
    if value is None:
        return [], []
    if not isinstance(value, list):
        return [], [f"{key} must be a list, got {type(value).__name__}"]
    values: list[str] = []
    malformed: list[str] = []
    for item in value:
        if isinstance(item, str):
            values.append(item.strip())
        else:
            malformed.append(repr(item))
    return values, malformed


def _clean_html(value: str) -> str:
    return " ".join(HTML_TAG_RE.sub(" ", unescape(value)).split())


def load_kegg(path: Path) -> dict[str, dict[str, Any]]:
    _require_file(path, "KEGG reference")
    result: dict[str, dict[str, Any]] = {}
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, 1):
            line = line.rstrip("\r\n")
            if not line.strip():
                continue
            fields = line.split("\t", 1)
            if len(fields) != 2 or not fields[0].strip():
                raise AuditInputError(f"malformed KEGG row at {path}:{line_number}")
            identifier, definition = fields[0].strip(), fields[1].strip()
            embedded = [
                value
                for value in EC_BRACKET_RE.findall(definition)
                for value in value.split()
            ]
            result[identifier] = {"definition": definition, "ecs": embedded}
    return result


def load_cogs(path: Path) -> dict[str, dict[str, str]]:
    _require_file(path, "COG reference")
    result: dict[str, dict[str, str]] = {}
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, 1):
            line = line.rstrip("\r\n")
            if not line.strip():
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                raise AuditInputError(f"malformed COG row at {path}:{line_number}")
            identifier = fields[0].strip()
            result[identifier] = {
                "category": fields[1].strip(),
                "description": fields[2].strip(),
                "gene": fields[3].strip() if len(fields) > 3 else "",
                "pathway": fields[4].strip() if len(fields) > 4 else "",
            }
    return result


def load_explorenz(
    path: Path,
) -> tuple[
    dict[str, dict[str, str]],
    dict[tuple[str, str, str], str],
    dict[str, dict[str, str]],
]:
    """Stream the entry, class, and hist tables from the ExplorEnz dump."""
    _require_file(path, "ExplorEnz XML reference")
    entries: dict[str, dict[str, str]] = {}
    classes: dict[tuple[str, str, str], str] = {}
    history: dict[str, dict[str, str]] = {}
    current_table: str | None = None
    current_row: dict[str, str] | None = None
    current_field: str | None = None

    try:
        context = ET.iterparse(path, events=("start", "end"))
        for event, element in context:
            tag = element.tag.rsplit("}", 1)[-1]
            if event == "start":
                if tag == "table_data":
                    current_table = element.attrib.get("name")
                elif tag == "row" and current_table in {"entry", "class", "hist"}:
                    current_row = {}
                elif tag == "field" and current_table in {"entry", "class", "hist"}:
                    current_field = element.attrib.get("name") or element.attrib.get(
                        "Field"
                    )
            else:
                if tag == "field" and current_row is not None and current_field:
                    current_row[current_field] = element.text or ""
                    current_field = None
                    element.clear()
                elif tag == "field":
                    # Discard fields from cite/refs/html and table schemas.
                    element.clear()
                elif tag == "row":
                    if current_row is not None:
                        if current_table == "entry" and current_row.get("ec_num"):
                            entries[current_row["ec_num"]] = dict(current_row)
                        elif current_table == "class" and current_row.get("class"):
                            key = (
                                current_row.get("class", ""),
                                current_row.get("subclass", ""),
                                current_row.get("subsubclass", ""),
                            )
                            classes[key] = _clean_html(current_row.get("heading", ""))
                        elif current_table == "hist" and current_row.get("ec_num"):
                            history[current_row["ec_num"]] = dict(current_row)
                        current_row = None
                    # Rows in cite/html/refs are intentionally discarded too;
                    # leaving them attached to the XML root defeats streaming.
                    element.clear()
                elif tag == "table_data":
                    current_table = None
                    element.clear()
                elif tag not in {"mysqldump", "database"}:
                    element.clear()
    except ET.ParseError as exc:
        raise AuditInputError(f"malformed ExplorEnz XML {path}: {exc}") from exc
    return entries, classes, history


def load_fallback_ec(
    path: Path,
) -> tuple[dict[str, dict[str, str]], dict[str, dict[str, str]]]:
    """Load ExPASy enzyme.dat or a simple ``EC<TAB>name`` fallback."""
    _require_file(path, "EC fallback reference")
    entries: dict[str, dict[str, str]] = {}
    history: dict[str, dict[str, str]] = {}
    first_nonempty = ""
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.strip():
                first_nonempty = line.rstrip("\r\n")
                break
    is_tsv = "\t" in first_nonempty
    if is_tsv:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            for line_number, line in enumerate(handle, 1):
                line = line.rstrip("\r\n")
                if not line.strip():
                    continue
                fields = line.split("\t", 1)
                if len(fields) != 2:
                    raise AuditInputError(
                        f"malformed EC TSV row at {path}:{line_number}"
                    )
                ec, name = fields[0].strip(), fields[1].strip()
                entries[ec] = {"ec_num": ec, "accepted_name": name}
                history[ec] = {"ec_num": ec, "action": "created", "note": ""}
        return entries, history

    current_ec: str | None = None
    current_name: list[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("ID   "):
                current_ec = line[5:].strip()
                current_name = []
            elif line.startswith("DE   ") and current_ec:
                current_name.append(line[5:].strip())
            elif line.startswith("//") and current_ec:
                entries[current_ec] = {
                    "ec_num": current_ec,
                    "accepted_name": " ".join(current_name),
                }
                history[current_ec] = {
                    "ec_num": current_ec,
                    "action": "created",
                    "note": "",
                }
                current_ec = None
                current_name = []
    if not entries:
        raise AuditInputError(
            f"EC fallback is neither usable TSV nor enzyme.dat: {path}"
        )
    return entries, history


def load_pfams(
    path: Path, requested: set[str] | None = None
) -> dict[str, dict[str, str]]:
    """Index Pfam aliases while streaming only model header lines."""
    _require_file(path, "Pfam reference")
    if requested is not None and not requested:
        return {}
    result: dict[str, dict[str, str]] = {}
    name: str | None = None
    accession: str | None = None
    description: str | None = None
    # The HMM payload is ASCII-compatible and very large.  Reading bytes keeps
    # Python's text decoder from retaining large transient Unicode buffers on
    # Windows while still decoding the three short metadata fields.
    with path.open("rb") as handle:
        for line in handle:
            if line.startswith(b"NAME  "):
                name = line[6:].decode("utf-8", errors="replace").strip()
            elif line.startswith(b"ACC   "):
                accession = line[6:].decode("utf-8", errors="replace").strip()
            elif line.startswith(b"DESC  "):
                description = line[6:].decode("utf-8", errors="replace").strip()
            elif line.startswith(b"//"):
                if accession:
                    unversioned = accession.split(".", 1)[0]
                    if (
                        requested is None
                        or unversioned in requested
                        or accession in requested
                    ):
                        record = {
                            "accession": unversioned,
                            "full_accession": accession,
                            "name": name or "",
                            "description": description or "",
                        }
                        result[unversioned] = record
                        result[accession] = record
                name = None
                accession = None
                description = None
    return result


def load_references(
    kegg_path: Path,
    cog_path: Path,
    ec_xml_path: Path | None,
    ec_fallback_path: Path | None,
    pfam_path: Path,
    pfam_targets: set[str] | None = None,
) -> ReferenceData:
    kegg = load_kegg(kegg_path)
    cogs = load_cogs(cog_path)
    ec_source: str
    if ec_xml_path is not None and ec_xml_path.is_file():
        ec_entries, ec_classes, ec_history = load_explorenz(ec_xml_path)
        ec_source = str(ec_xml_path)
    elif ec_fallback_path is not None:
        ec_entries, ec_history = load_fallback_ec(ec_fallback_path)
        ec_classes = {}
        ec_source = str(ec_fallback_path)
    elif ec_xml_path is not None:
        raise AuditInputError(
            f"ExplorEnz XML not found and no --ec fallback supplied: {ec_xml_path}"
        )
    else:
        raise AuditInputError("one of --ec-xml or --ec is required")
    pfams = load_pfams(pfam_path, pfam_targets)
    return ReferenceData(
        kegg, cogs, ec_entries, ec_classes, ec_history, pfams, ec_source
    )


def _marker_id(marker: dict[str, Any], index: int) -> str:
    value = marker.get("id")
    return value if isinstance(value, str) and value.strip() else f"<marker-{index}>"


def _catalog_name(marker: dict[str, Any]) -> str:
    value = marker.get("name", "")
    return value if isinstance(value, str) else str(value)


def _is_wildcard_ec(value: str) -> bool:
    return "-" in value


def _ec_class_key(value: str) -> tuple[str, str, str]:
    parts = value.split(".")
    return tuple(part if part != "-" else "0" for part in parts[:3])  # type: ignore[return-value]


def _ec_compatible(left: str, right: str) -> bool:
    """Return whether two ECs overlap at their declared specificity."""
    left_parts, right_parts = left.split("."), right.split(".")
    return all(a == b or a == "-" or b == "-" for a, b in zip(left_parts, right_parts))


def _definition_mismatch(marker_id: str, definition: str) -> str | None:
    """Flag only high-confidence, catalog-specific definition conflicts.

    A general natural-language similarity score would create false positives:
    KO definitions often describe one subunit while a marker names a whole
    complex.  These rules therefore cover only the incompatible assignments
    exposed by this snapshot and leave uncertain mappings as REVIEW candidates.
    """
    text = definition.casefold()
    if marker_id == "CYP153" and re.search(
        r"\balkb(?:1_2|1|2)?\b|alkane 1-monooxygenase", text
    ):
        return "KEGG definition is AlkB-family, not a CYP153 cytochrome P450"
    direct_conflicts = {
        "prmABCD": (
            r"\balcohol dehydrogenase\b",
            "an alcohol dehydrogenase, not a propane monooxygenase",
        ),
        "rubAB": (
            r"\byidH\b|apolipoprotein N-acyltransferase",
            "an unrelated membrane protein/acyltransferase",
        ),
        "almA": (
            r"\bargO\b|N-acetylglutamate synthase",
            "ArgO/N-acetylglutamate synthase, not AlmA",
        ),
        "dmpO_tomA": (
            r"dihydrobenzene|dihydroxynaphthalene|dihydrodiol dehydrogenase",
            "a downstream dihydrodiol dehydrogenase",
        ),
        "assA_masD": (r"\bbestrophin\b", "bestrophin, not alkylsuccinate synthase"),
        "bssABC": (
            r"GPI mannosyltransferase",
            "a GPI mannosyltransferase, not benzylsuccinate synthase",
        ),
        "nmsA": (
            r"S-adenosyl-L-methionine hydrolase",
            "an S-adenosyl-L-methionine hydrolase",
        ),
        "cmdA_ebdA": (
            r"XRE family transcriptional regulator|transcriptional regulator",
            "a transcriptional regulator, not ethylbenzene dehydrogenase",
        ),
        "rhlA": (r"trehalose synthase", "trehalose synthase, not HAA synthase"),
        "rhlB": (
            r"mannosylfructose-phosphate synthase",
            "mannosylfructose-phosphate synthase, not rhamnosyltransferase I",
        ),
        "rhlC": (
            r"N-acetylhexosamine 1-kinase",
            "N-acetylhexosamine kinase, not rhamnosyltransferase II",
        ),
        "sfp": (
            r"translation initiation factor",
            "translation initiation factor IF-1, not a phosphopantetheinyl transferase",
        ),
        "wza_wzb_wzc": (
            r"ABC-2 type transport system ATP-binding protein",
            "an ABC ATPase, not Wza/Wzb/Wzc export",
        ),
    }
    conflict = direct_conflicts.get(marker_id)
    if conflict and re.search(conflict[0], text, re.IGNORECASE):
        return "KEGG definition is " + conflict[1]
    if marker_id in {"srfA", "licA_D"} and re.search(
        r"\b(?:pps|fen|plipastatin|fengycin|rifamycin|heptose)\b", text
    ):
        return "KEGG definition is a different lipopeptide/PKS family"
    if marker_id == "fengycin_iturin" and re.search(r"\b(?:rifamycin|heptose)\b", text):
        return "KEGG definition is not a fengycin/iturin lipopeptide synthetase"
    return None


def _generic_reference(description: str) -> bool:
    return bool(GENERIC_TERMS_RE.search(description))


def _malformed_finding(
    marker_id: str, kind: str, value: str, name: str, note: str
) -> Finding:
    return Finding(marker_id, kind, value, "MALFORMED", "N/A", name, "MALFORMED", note)


def audit_markers(
    markers: list[dict[str, Any]], references: ReferenceData
) -> list[Finding]:
    findings: list[Finding] = []
    by_identifier: dict[tuple[str, str], list[Finding]] = defaultdict(list)
    marker_data: dict[str, dict[str, Any]] = {}

    for index, marker in enumerate(markers, 1):
        marker_id = _marker_id(marker, index)
        name = _catalog_name(marker)
        marker_data[marker_id] = {"marker": marker, "name": name, "kos": [], "ecs": []}

        if "_malformed_marker" in marker:
            findings.append(
                _malformed_finding(
                    marker_id,
                    "MARKER",
                    repr(marker["_malformed_marker"]),
                    name,
                    "marker entry must be a mapping",
                )
            )
            continue

        fields = (
            ("kos", "KO", KO_RE),
            ("ecs", "EC", EC_RE),
            ("cogs", "COG", COG_RE),
            ("pfams", "Pfam", PFAM_RE),
        )
        for field, kind, pattern in fields:
            values, malformed = _string_list(marker, field)
            for note in malformed:
                findings.append(
                    _malformed_finding(marker_id, kind, f"<{field}>", name, note)
                )
            for value in values:
                if not value or not pattern.fullmatch(value):
                    finding = _malformed_finding(
                        marker_id,
                        kind,
                        value,
                        name,
                        f"invalid {kind} identifier format",
                    )
                    findings.append(finding)
                    continue
                if kind == "KO":
                    reference = references.kegg.get(value)
                    if reference is None:
                        finding = Finding(
                            marker_id,
                            kind,
                            value,
                            "ABSENT",
                            "N/A",
                            name,
                            "MISSING_FROM_REFERENCE",
                            "KO identifier not present in KEGG snapshot",
                        )
                    else:
                        definition = str(reference["definition"])
                        mismatch = _definition_mismatch(marker_id, definition)
                        status = "DEFINITION_MISMATCH" if mismatch else "PASS"
                        note = mismatch or ""
                        if reference.get("ecs"):
                            note = f"{note}; " if note else ""
                            note += "embedded ECs: " + ", ".join(reference["ecs"])
                        finding = Finding(
                            marker_id,
                            kind,
                            value,
                            "EXISTS",
                            definition,
                            name,
                            status,
                            note,
                        )
                        marker_data[marker_id]["kos"].append(value)
                elif kind == "EC":
                    marker_data[marker_id]["ecs"].append(value)
                    if _is_wildcard_ec(value):
                        heading = references.ec_classes.get(
                            _ec_class_key(value), "Unresolved EC class heading"
                        )
                        finding = Finding(
                            marker_id,
                            kind,
                            value,
                            "WILDCARD",
                            heading,
                            name,
                            "NONSPECIFIC",
                            "wildcard EC is broad contextual evidence",
                        )
                    else:
                        entry = references.ec_entries.get(value)
                        history = references.ec_history.get(value, {})
                        action = history.get("action", "").casefold()
                        ref_name = entry.get("accepted_name", "") if entry else ""
                        if action in {"deleted", "transferred"}:
                            successors = sorted(
                                set(SUCCESSOR_RE.findall(history.get("note", "")))
                            )
                            successor_note = (
                                ", ".join(successors)
                                if successors
                                else "no successor EC linked"
                            )
                            finding = Finding(
                                marker_id,
                                kind,
                                value,
                                f"HISTORICAL_{action.upper()}",
                                ref_name or history.get("note", ""),
                                name,
                                "REPLACE",
                                f"{action} EC; successor(s): {successor_note}",
                            )
                        elif action == "created":
                            finding = Finding(
                                marker_id,
                                kind,
                                value,
                                "EXISTS",
                                ref_name,
                                name,
                                "PASS",
                                "",
                            )
                        elif entry is None and not history:
                            finding = Finding(
                                marker_id,
                                kind,
                                value,
                                "ABSENT",
                                "N/A",
                                name,
                                "MISSING_FROM_REFERENCE",
                                "EC not present in EC snapshot",
                            )
                        else:
                            finding = Finding(
                                marker_id,
                                kind,
                                value,
                                "UNCLASSIFIED",
                                ref_name,
                                name,
                                "REVIEW",
                                "EC entry has no authoritative hist.action",
                            )
                elif kind == "COG":
                    reference = references.cogs.get(value)
                    if reference is None:
                        finding = Finding(
                            marker_id,
                            kind,
                            value,
                            "ABSENT",
                            "N/A",
                            name,
                            "MISSING_FROM_REFERENCE",
                            "COG identifier not present in COG snapshot",
                        )
                    else:
                        ref_name = (
                            f"[{reference['category']}] {reference['description']}"
                        )
                        status = (
                            "NONSPECIFIC"
                            if _generic_reference(reference["description"])
                            else "PASS"
                        )
                        note = (
                            "broad/generic COG definition"
                            if status == "NONSPECIFIC"
                            else ""
                        )
                        finding = Finding(
                            marker_id,
                            kind,
                            value,
                            "EXISTS",
                            ref_name,
                            name,
                            status,
                            note,
                        )
                else:
                    reference = references.pfams.get(value) or references.pfams.get(
                        value.split(".", 1)[0]
                    )
                    if reference is None:
                        finding = Finding(
                            marker_id,
                            kind,
                            value,
                            "ABSENT",
                            "N/A",
                            name,
                            "MISSING_FROM_REFERENCE",
                            "Pfam accession not present in Pfam-A snapshot",
                        )
                    else:
                        ref_name = f"{reference['name']}: {reference['description']}"
                        status = (
                            "NONSPECIFIC"
                            if _generic_reference(reference["description"])
                            else "PASS"
                        )
                        note = (
                            "broad/generic Pfam definition"
                            if status == "NONSPECIFIC"
                            else ""
                        )
                        finding = Finding(
                            marker_id,
                            kind,
                            value,
                            "EXISTS",
                            ref_name,
                            name,
                            status,
                            note,
                        )

                findings.append(finding)
                by_identifier[(kind, value)].append(finding)

    # Mark shared identifiers, retaining a more specific definition mismatch
    # where one is already established.
    for (kind, identifier), occurrences in by_identifier.items():
        occurrences = [item for item in occurrences if item is not None]
        marker_ids = list(dict.fromkeys(item.marker_id for item in occurrences))
        if len(marker_ids) < 2:
            continue
        shared_note = f"shared across markers: {', '.join(marker_ids)}"
        for finding in occurrences:
            if shared_note not in finding.notes:
                finding.notes = f"{finding.notes}; {shared_note}".lstrip("; ")
            if finding.audit_status == "PASS":
                finding.audit_status = "DUPLICATED"

    # KO-definition ECs are a useful consistency check, but only at matching
    # EC specificity.  A wildcard on either side is compatible with its child.
    for marker_id, data in marker_data.items():
        marker_ecs = [value for value in data.get("ecs", []) if EC_RE.fullmatch(value)]
        embedded: list[tuple[str, str]] = []
        for ko in data.get("kos", []):
            for ec in references.kegg.get(ko, {}).get("ecs", []):
                if EC_VALUE_RE.fullmatch(ec):
                    embedded.append((ko, ec))
        if not embedded:
            continue
        extras = [
            (ko, ec)
            for ko, ec in embedded
            if not any(_ec_compatible(ec, own) for own in marker_ecs)
        ]
        if not extras:
            continue
        reference_name = "; ".join(f"{ko}: {ec}" for ko, ec in extras)
        catalog = ", ".join(marker_ecs) or "<no marker EC>"
        findings.append(
            Finding(
                marker_id,
                "CONSISTENCY",
                marker_id,
                "KO_EC_CONTRADICTION",
                reference_name,
                catalog,
                "REVIEW",
                "KO definition embeds EC(s) not compatible with marker ecs",
            )
        )
        for finding in findings:
            if (
                finding.marker_id == marker_id
                and finding.identifier_type == "KO"
                and finding.identifier in {ko for ko, _ in extras}
            ):
                if "KO/EC contradiction" not in finding.notes:
                    finding.notes = f"{finding.notes}; KO/EC contradiction".lstrip("; ")
                if finding.audit_status == "PASS":
                    finding.audit_status = "REVIEW"

    for finding in findings:
        if finding.audit_status not in STATUSES:
            raise AssertionError(f"unexpected audit status: {finding.audit_status}")
    return findings


def _write_tsv(path: Path, rows: Iterable[Iterable[str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows(rows)


def write_reports(
    outdir: Path,
    findings: list[Finding],
    references: ReferenceData,
    input_paths: dict[str, Path],
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    _write_tsv(
        outdir / "meor_identifier_audit.tsv",
        [AUDIT_HEADER, *(finding.row() for finding in findings)],
    )
    for kind, filename in (
        ("KO", "ko_diff.tsv"),
        ("COG", "cog_diff.tsv"),
        ("EC", "ec_diff.tsv"),
        ("Pfam", "pfam_diff.tsv"),
    ):
        rows = [
            AUDIT_HEADER,
            *(finding.row() for finding in findings if finding.identifier_type == kind),
        ]
        _write_tsv(outdir / filename, rows)

    status_counts = Counter(finding.audit_status for finding in findings)
    catalog_counts = Counter(
        finding.identifier_type
        for finding in findings
        if finding.identifier_type in {"KO", "COG", "EC", "Pfam"}
    )
    lines = [
        "# MEOR identifier audit",
        "",
        "This report was generated offline from the supplied snapshots. The audit is read-only and does not modify the marker catalog.",
        "",
        "## Inputs",
        "",
    ]
    for label in ("markers", "kegg", "cog", "ec_xml", "pfam"):
        if label in input_paths:
            lines.append(f"- `{label}`: `{input_paths[label]}`")
    lines.extend(
        [
            f"- `ec_source_used`: `{references.ec_source}`",
            "",
            "## Coverage",
            "",
            "| Identifier type | Catalog assignments | Reference entries indexed |",
            "| --- | ---: | ---: |",
            f"| KO | {catalog_counts['KO']} | {len(references.kegg):,} |",
            f"| COG | {catalog_counts['COG']} | {len(references.cogs):,} |",
            f"| EC | {catalog_counts['EC']} | {len(references.ec_entries):,} entries; {len(references.ec_history):,} history rows |",
            f"| Pfam | {catalog_counts['Pfam']} | {len(references.pfams):,} aliases |",
            "",
            "## Status counts",
            "",
            "| Status | Count |",
            "| --- | ---: |",
        ]
    )
    for status in sorted(status_counts):
        lines.append(f"| `{status}` | {status_counts[status]} |")
    lines.extend(
        [
            "",
            "## Findings requiring review or remediation",
            "",
            "| Marker | Type | Identifier | Status | Notes |",
            "| --- | --- | --- | --- | --- |",
        ]
    )
    for finding in findings:
        if finding.audit_status != "PASS":
            note = finding.notes.replace("|", "\\|")
            lines.append(
                f"| `{finding.marker_id}` | `{finding.identifier_type}` | `{finding.identifier}` | `{finding.audit_status}` | {note} |"
            )
    lines.extend(
        [
            "",
            "## Reference index counts",
            "",
            f"- KEGG definitions: {len(references.kegg):,} (embedded EC annotations are parsed from `[EC:...]` brackets).",
            f"- COG definitions: {len(references.cogs):,}; categories are preserved as complete strings.",
            f"- EC entries: {len(references.ec_entries):,}; history rows: {len(references.ec_history):,}; class headings: {len(references.ec_classes):,}.",
            f"- Pfam aliases: {len(references.pfams):,}; aliases include versioned and unversioned accessions.",
            "- Pfam indexing retains only catalog-requested accessions; the current catalog has no Pfam assignments, so the supplied Pfam snapshot was validated for presence but no model aliases were retained.",
            "",
            "## Notes",
            "",
            "`NONSPECIFIC` wildcard ECs are retained as broad contextual evidence. `DUPLICATED` identifies reuse across markers; a high-confidence incompatible definition is reported as `DEFINITION_MISMATCH`. `REVIEW` consistency findings require curator judgment.",
            "",
        ]
    )
    _write_tsv(
        outdir / "audit_status_counts.tsv",
        [
            ["audit_status", "count"],
            *([status, str(status_counts[status])] for status in sorted(status_counts)),
        ],
    )
    (outdir / "audit_summary.md").write_text("\n".join(lines), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Audit MEOR identifiers against offline reference snapshots"
    )
    parser.add_argument(
        "--markers",
        type=Path,
        default=Path("src/genbank_parser/data/meor/markers.yaml"),
    )
    parser.add_argument(
        "--kegg", type=Path, default=Path("audit-reference/kegg_ko_list_2026-08-02.tsv")
    )
    parser.add_argument(
        "--cog", type=Path, default=Path("audit-reference/cog-24.def.tab")
    )
    parser.add_argument(
        "--ec-xml", type=Path, default=Path("audit-reference/explorenz-enzyme-data.xml")
    )
    parser.add_argument(
        "--ec", type=Path, default=None, help="enzyme.dat or EC<TAB>name fallback"
    )
    parser.add_argument("--pfam", type=Path, default=Path("audit-reference/Pfam-A.hmm"))
    parser.add_argument("--outdir", type=Path, default=Path("audit"))
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        markers = load_markers(args.markers)
        pfam_targets: set[str] = set()
        for marker in markers:
            values, _ = _string_list(marker, "pfams")
            pfam_targets.update(values)
        references = load_references(
            args.kegg, args.cog, args.ec_xml, args.ec, args.pfam, pfam_targets
        )
        findings = audit_markers(markers, references)
        write_reports(
            args.outdir,
            findings,
            references,
            {
                "markers": args.markers,
                "kegg": args.kegg,
                "cog": args.cog,
                "ec_xml": args.ec_xml,
                "pfam": args.pfam,
            },
        )
    except AuditInputError as exc:
        print(f"audit error: {exc}", file=sys.stderr)
        return 2
    print(f"Audited {len(markers)} markers and wrote reports to {args.outdir}")
    counts = Counter(finding.audit_status for finding in findings)
    for status in sorted(counts):
        print(f"  {status}: {counts[status]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
