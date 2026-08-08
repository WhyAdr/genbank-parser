"""Typed data models representing GenBank documents, records, and features."""

from __future__ import annotations

import collections
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from Bio.Seq import Seq
from Bio.SeqFeature import (
    AfterPosition,
    BeforePosition,
    CompoundLocation,
    FeatureLocation,
)


@dataclass
class GenBankFeature:
    """Represents an individual annotated genomic feature with location semantics."""

    record_id: str
    record_index: int
    feature_index: int
    type: str
    location: FeatureLocation | CompoundLocation | Any
    qualifiers: dict[str, list[str]] = field(default_factory=dict)
    record_length: int = 0
    topology: str | None = None
    raw_feature: Any = None

    @property
    def start(self) -> int:
        """1-based inclusive start coordinate."""
        try:
            return int(self.location.start) + 1
        except (AttributeError, TypeError, ValueError):
            return 1

    @property
    def end(self) -> int:
        """1-based inclusive end coordinate."""
        try:
            return int(self.location.end)
        except (AttributeError, TypeError, ValueError):
            return 0

    @property
    def strand(self) -> int | None:
        """Raw strand state: +1, -1, 0, or None."""
        try:
            return self.location.strand
        except AttributeError:
            return None

    @property
    def strand_symbol(self) -> str:
        """Strand symbol: '+' (+1), '-' (-1), '?' (0), '.' (None)."""
        s = self.strand
        if s == 1:
            return "+"
        if s == -1:
            return "-"
        if s == 0:
            return "?"
        return "."

    @property
    def length(self) -> int:
        """Biological sequence length (excludes introns/gaps in compound features)."""
        try:
            return len(self.location)
        except (AttributeError, TypeError, ValueError):
            return max(0, self.end - self.start + 1)

    @property
    def genomic_span(self) -> int:
        """Total genomic span from start to end (max_end - min_start + 1)."""
        return max(0, self.end - self.start + 1)

    @property
    def is_compound(self) -> bool:
        """True if the feature has a compound location (join/order)."""
        return isinstance(self.location, CompoundLocation)

    @property
    def is_partial_start(self) -> bool:
        """True if feature start coordinate is partial (<start)."""
        try:
            sp = self.location.start
            return isinstance(sp, BeforePosition) or "<" in str(sp)
        except (AttributeError, TypeError, ValueError):
            return False

    @property
    def is_partial_end(self) -> bool:
        """True if feature end coordinate is partial (>end)."""
        try:
            ep = self.location.end
            return isinstance(ep, AfterPosition) or ">" in str(ep)
        except (AttributeError, TypeError, ValueError):
            return False

    @property
    def is_partial(self) -> bool:
        """True if either start or end coordinate is partial."""
        return self.is_partial_start or self.is_partial_end

    @property
    def join_segments(self) -> list[tuple[int, int]]:
        """List of 1-based (start, end) tuples for compound location parts."""
        segments: list[tuple[int, int]] = []
        if isinstance(self.location, CompoundLocation):
            for part in self.location.parts:
                segments.append((int(part.start) + 1, int(part.end)))
        return segments

    @property
    def locus_tag(self) -> str:
        return self.get_qual("locus_tag")

    @property
    def gene(self) -> str:
        return self.get_qual("gene")

    @property
    def product(self) -> str:
        return self.get_qual("product")

    @property
    def protein_id(self) -> str:
        return self.get_qual("protein_id")

    @property
    def translation(self) -> str:
        return self.get_qual("translation")

    @property
    def codon_start(self) -> int:
        val = self.get_qual("codon_start", "1")
        try:
            iv = int(val)
            return iv if iv in (1, 2, 3) else 1
        except ValueError:
            return 1

    @property
    def transl_table(self) -> int:
        val = self.get_qual("transl_table", "11")
        try:
            return int(val)
        except ValueError:
            return 11

    @property
    def is_pseudo(self) -> bool:
        if self.type.lower() == "pseudogene":
            return True
        return bool(self.qualifiers.get("pseudo") or self.qualifiers.get("pseudogene"))

    def get_qual(self, key: str, default: str = "") -> str:
        """Return first qualifier value for key, or default."""
        vals = self.qualifiers.get(key, [])
        return vals[0] if vals else default

    def get_quals(self, key: str) -> list[str]:
        """Return all qualifier values for key."""
        return list(self.qualifiers.get(key, []))

    def extract(self, record_seq: Seq | str) -> Seq:
        """Extract the biological sequence for this feature from parent record sequence."""
        if isinstance(record_seq, str):
            record_seq = Seq(record_seq)
        return self.location.extract(record_seq)

    def __getitem__(self, key: str) -> Any:
        """Dictionary-like access for backwards compatibility."""
        if key == "type":
            return self.type
        if key == "start":
            return self.start
        if key == "end":
            return self.end
        if key == "strand":
            return self.strand_symbol
        if key == "contig":
            return self.record_id
        if key == "qualifiers":
            return self.qualifiers
        if key in ("line", "feature_index"):
            return self.feature_index
        if key == "partial_start":
            return self.is_partial_start
        if key == "partial_end":
            return self.is_partial_end
        if key == "join_segments":
            return self.join_segments
        if key == "location":
            return self.location
        if key == "record_length":
            return self.record_length
        if key == "topology":
            return self.topology
        if key in self.qualifiers:
            return self.qualifiers[key]
        raise KeyError(key)

    def __contains__(self, key: str) -> bool:
        known_keys = {
            "type",
            "start",
            "end",
            "strand",
            "contig",
            "qualifiers",
            "line",
            "feature_index",
            "partial_start",
            "partial_end",
            "join_segments",
            "location",
            "record_length",
            "topology",
        }
        return key in known_keys or key in self.qualifiers

    def get(self, key: str, default: Any = None) -> Any:
        try:
            return self[key]
        except KeyError:
            return default

    def to_dict(self) -> dict[str, Any]:
        """Convert to legacy feature dictionary format."""
        return {
            "type": self.type,
            "start": self.start,
            "end": self.end,
            "strand": self.strand_symbol,
            "contig": self.record_id,
            "qualifiers": collections.defaultdict(list, self.qualifiers),
            "line": self.feature_index,
            "feature_index": self.feature_index,
            "partial_start": self.is_partial_start,
            "partial_end": self.is_partial_end,
            "join_segments": self.join_segments,
            "location": self.location,
            "record_length": self.record_length,
            "topology": self.topology,
        }


@dataclass
class GenBankRecord:
    """Represents a single contig, chromosome, or plasmid record."""

    id: str
    name: str
    description: str
    seq: Seq
    length: int
    topology: str | None = None
    molecule_type: str | None = None
    division: str | None = None
    date: str | None = None
    annotations: dict[str, Any] = field(default_factory=dict)
    features: list[GenBankFeature] = field(default_factory=list)

    @property
    def cds_features(self) -> list[GenBankFeature]:
        return [f for f in self.features if f.type == "CDS"]

    @property
    def gc_content(self) -> float:
        """Percentage of G and C nucleotides in sequence."""
        if not self.seq:
            return 0.0
        s_upper = str(self.seq).upper()
        gc = s_upper.count("G") + s_upper.count("C")
        total_acgt = gc + s_upper.count("A") + s_upper.count("T")
        return (100.0 * gc / total_acgt) if total_acgt > 0 else 0.0

    @property
    def coding_density(self) -> float:
        """Percentage of record nucleotides covered by the CDS union."""
        if self.length == 0:
            return 0.0
        return 100.0 * self.nonredundant_coding_bp / self.length

    @property
    def cds_feature_bp_sum(self) -> int:
        """Sum of biological CDS lengths, including overlaps."""
        return sum(f.length for f in self.cds_features)

    @property
    def nonredundant_coding_bp(self) -> int:
        """Number of record bases covered by the union of CDS locations."""
        intervals: list[tuple[int, int]] = []
        for feature in self.cds_features:
            location = feature.location
            parts = getattr(location, "parts", (location,))
            for part in parts:
                start = max(0, int(part.start))
                end = min(self.length, int(part.end)) if self.length else int(part.end)
                if end > start:
                    intervals.append((start, end))

        if not intervals:
            return 0

        intervals.sort()
        covered = 0
        current_start, current_end = intervals[0]
        for start, end in intervals[1:]:
            if start <= current_end:
                current_end = max(current_end, end)
            else:
                covered += current_end - current_start
                current_start, current_end = start, end
        return covered + current_end - current_start

    def find_locus_features(self, locus_tag: str) -> list[GenBankFeature]:
        """Return every feature carrying ``locus_tag`` in source order."""
        return [f for f in self.features if f.locus_tag == locus_tag]

    def find_locus(
        self,
        locus_tag: str,
        prefer_type: str | None = "CDS",
    ) -> GenBankFeature | None:
        """Find a locus, preferring the requested feature type when present."""
        matches = self.find_locus_features(locus_tag)
        if not matches:
            return None
        if prefer_type is not None:
            preferred = prefer_type.casefold()
            for feat in matches:
                if feat.type.casefold() == preferred:
                    return feat
        return matches[0]

    def search_gene(self, gene_name: str) -> list[GenBankFeature]:
        """Find features with matching /gene qualifier."""
        return [f for f in self.features if f.gene.casefold() == gene_name.casefold()]

    def features_by_type(self, ftype: str) -> list[GenBankFeature]:
        """Filter features by type (e.g. 'CDS', 'tRNA', 'rRNA')."""
        return [f for f in self.features if f.type == ftype]


@dataclass
class GenBankDocument:
    """Represents an entire parsed GenBank file with one or more records."""

    path: Path
    records: list[GenBankRecord] = field(default_factory=list)

    def __len__(self) -> int:
        return len(self.records)

    def __iter__(self) -> Iterator[GenBankRecord]:
        return iter(self.records)

    def __getitem__(self, idx: int) -> GenBankRecord:
        return self.records[idx]

    @property
    def total_length(self) -> int:
        return sum(r.length for r in self.records)

    @property
    def total_features(self) -> int:
        return sum(len(r.features) for r in self.records)

    @property
    def all_features(self) -> list[GenBankFeature]:
        res: list[GenBankFeature] = []
        for r in self.records:
            res.extend(r.features)
        return res

    def get_record(self, record_id: str) -> GenBankRecord | None:
        """Retrieve record by ID or name."""
        for rec in self.records:
            if rec.id == record_id or rec.name == record_id:
                return rec
        return None

    def find_locus(
        self,
        locus_tag: str,
        prefer_type: str | None = "CDS",
    ) -> tuple[GenBankRecord, GenBankFeature] | None:
        """Search records, preferring a matching feature type (normally CDS)."""
        for rec in self.records:
            feat = rec.find_locus(locus_tag, prefer_type=prefer_type)
            if feat is not None:
                return (rec, feat)
        return None

    def find_cds(self, locus_tag: str) -> tuple[GenBankRecord, GenBankFeature] | None:
        """Find a CDS by locus tag across all records."""
        return self.find_locus(locus_tag, prefer_type="CDS")
