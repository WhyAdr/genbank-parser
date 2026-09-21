"""GenBank structural and biological semantics validator."""

from __future__ import annotations

import argparse
import collections
import io
import json
from dataclasses import asdict, dataclass
from pathlib import Path

from Bio.Data import CodonTable

from .io import _source_label, read_genbank
from .model import GenBankFeature


@dataclass
class ValidationFinding:
    severity: str  # 'ERROR', 'WARNING', 'INFO'
    code: str
    record_id: str
    locus_tag: str
    feature_type: str
    coordinates: str
    message: str


@dataclass(frozen=True)
class ValidationReport:
    """Renderer-independent validation result and summary statistics."""

    source: str
    record_count: int
    feature_count: int
    total_length: int
    type_counts: dict[str, int]
    strand_counts: dict[str, int]
    cds_count: int
    named_gene_count: int
    hypothetical_count: int
    locus_tag_count: int
    unique_locus_tag_count: int
    findings: tuple[ValidationFinding, ...]

    @property
    def severity_counts(self) -> dict[str, int]:
        return {
            severity: sum(1 for finding in self.findings if finding.severity == severity)
            for severity in ("ERROR", "WARNING", "INFO")
        }


def _translation_table_id(feature: GenBankFeature) -> int | None:
    """Return a known NCBI translation-table ID, or ``None`` if invalid."""
    raw = feature.get_qual("transl_table")
    if not raw:
        return 11
    try:
        table_id = int(raw)
    except ValueError:
        return None
    return table_id if table_id in CodonTable.unambiguous_dna_by_id else None


def build_validation_report(filepath: str | Path) -> ValidationReport:
    doc = read_genbank(filepath)
    if len(doc.records) == 0:
        raise ValueError("No records parsed. Check file format.")

    all_features = doc.all_features
    if not all_features:
        raise ValueError("No features parsed. Check file format.")

    type_counts = collections.Counter(f.type for f in all_features)
    strand_counts = collections.Counter(f.strand_symbol for f in all_features)
    all_quals: set[str] = set()
    for f in all_features:
        all_quals.update(f.qualifiers.keys())

    cdss = [f for f in all_features if f.type == "CDS"]

    # Locus tag indexing
    locus_tags: list[str] = []
    for f in all_features:
        lt = f.locus_tag
        if lt:
            locus_tags.append(lt)

    unique_tags = set(locus_tags)
    named_genes = [f for f in cdss if f.gene]
    products = [f.product for f in cdss if f.product]
    hypothetical = [
        p
        for p in products
        if "hypothetical" in p.lower()
        or "domain-containing" in p.lower()
        or "duf" in p.lower()
    ]

    findings: list[ValidationFinding] = []

    # Check for duplicate locus tags across genuinely different genes
    locus_tag_map: dict[str, list[GenBankFeature]] = collections.defaultdict(list)
    for f in all_features:
        if f.locus_tag:
            locus_tag_map[f.locus_tag].append(f)

    for tag, feats in locus_tag_map.items():
        # Gene + CDS sharing tag is expected. Multiple distinct CDSs with same tag is suspicious.
        cds_in_tag = [f for f in feats if f.type == "CDS"]
        if len(cds_in_tag) > 1:
            coords = ", ".join(f"{f.record_id}:{f.start}..{f.end}" for f in cds_in_tag)
            findings.append(
                ValidationFinding(
                    severity="WARNING",
                    code="DUPLICATE_LOCUS_TAG",
                    record_id=cds_in_tag[0].record_id,
                    locus_tag=tag,
                    feature_type="CDS",
                    coordinates=coords,
                    message=f"Locus tag '{tag}' assigned to {len(cds_in_tag)} distinct CDS features ({coords})",
                )
            )

    # Feature-level validation
    for rec in doc.records:
        rec_seq = rec.seq
        has_seq = len(rec_seq) > 0

        for f in rec.features:
            coords_str = f"{rec.id}:{f.start}..{f.end}({f.strand_symbol})"
            tag = f.locus_tag or "-"

            if f.start < 1:
                findings.append(
                    ValidationFinding(
                        severity="ERROR",
                        code="FEATURE_START_BELOW_ONE",
                        record_id=rec.id,
                        locus_tag=tag,
                        feature_type=f.type,
                        coordinates=coords_str,
                        message=f"Feature start ({f.start}) is below one",
                    )
                )

            parts = list(getattr(f.location, "parts", (f.location,))) if f.location is not None else []
            if not parts:
                findings.append(
                    ValidationFinding(
                        severity="ERROR",
                        code="EMPTY_LOCATION",
                        record_id=rec.id,
                        locus_tag=tag,
                        feature_type=f.type,
                        coordinates=coords_str,
                        message="Feature has no location segments",
                    )
                )
            for part in parts:
                try:
                    part_start = int(part.start) + 1
                    part_end = int(part.end)
                except (AttributeError, TypeError, ValueError):
                    findings.append(
                        ValidationFinding(
                            severity="ERROR",
                            code="INVALID_LOCATION_SEGMENT",
                            record_id=rec.id,
                            locus_tag=tag,
                            feature_type=f.type,
                            coordinates=coords_str,
                            message="Feature contains an invalid location segment",
                        )
                    )
                    continue
                if rec.length > 0 and (part_start < 1 or part_end > rec.length):
                    findings.append(
                        ValidationFinding(
                            severity="ERROR",
                            code="COMPOUND_SEGMENT_OUT_OF_BOUNDS",
                            record_id=rec.id,
                            locus_tag=tag,
                            feature_type=f.type,
                            coordinates=coords_str,
                            message=(
                                f"Location segment ({part_start}..{part_end}) exceeds "
                                f"record bounds (1..{rec.length})"
                            ),
                        )
                    )

            # Out of bounds check
            if rec.length > 0 and f.end > rec.length:
                findings.append(
                    ValidationFinding(
                        severity="ERROR",
                        code="COORDINATE_OUT_OF_BOUNDS",
                        record_id=rec.id,
                        locus_tag=tag,
                        feature_type=f.type,
                        coordinates=coords_str,
                        message=f"Feature end ({f.end}) exceeds record length ({rec.length})",
                    )
                )

            if f.length <= 0:
                findings.append(
                    ValidationFinding(
                        severity="ERROR",
                        code="ZERO_LENGTH_FEATURE",
                        record_id=rec.id,
                        locus_tag=tag,
                        feature_type=f.type,
                        coordinates=coords_str,
                        message=f"Feature has non-positive biological length ({f.length})",
                    )
                )

            if f.type == "CDS":
                is_pseudo = f.is_pseudo
                is_partial = f.is_partial

                raw_codon_start = f.get_qual("codon_start")
                if raw_codon_start:
                    try:
                        codon_start = int(raw_codon_start)
                    except ValueError:
                        codon_start = None
                    if codon_start not in (1, 2, 3):
                        findings.append(
                            ValidationFinding(
                                severity="ERROR",
                                code="INVALID_CODON_START",
                                record_id=rec.id,
                                locus_tag=tag,
                                feature_type="CDS",
                                coordinates=coords_str,
                                message=f"/codon_start must be 1, 2, or 3, got {raw_codon_start!r}",
                            )
                        )

                if not f.product:
                    findings.append(
                        ValidationFinding(
                            severity="WARNING",
                            code="MISSING_PRODUCT",
                            record_id=rec.id,
                            locus_tag=tag,
                            feature_type="CDS",
                            coordinates=coords_str,
                            message="CDS is missing /product qualifier",
                        )
                    )

                if not is_pseudo and not f.translation:
                    findings.append(
                        ValidationFinding(
                            severity="WARNING",
                            code="MISSING_TRANSLATION",
                            record_id=rec.id,
                            locus_tag=tag,
                            feature_type="CDS",
                            coordinates=coords_str,
                            message="Non-pseudogene CDS is missing /translation qualifier",
                        )
                    )

                # Translation verification if genome sequence is present
                if has_seq and not is_pseudo:
                    exceptional = any(
                        f.qualifiers.get(key)
                        for key in (
                            "transl_except",
                            "exception",
                            "ribosomal_slippage",
                            "frameshift",
                        )
                    )
                    if exceptional:
                        findings.append(
                            ValidationFinding(
                                severity="INFO",
                                code="EXCEPTIONAL_TRANSLATION_ANNOTATION",
                                record_id=rec.id,
                                locus_tag=tag,
                                feature_type="CDS",
                                coordinates=coords_str,
                                message="CDS carries an exceptional translation annotation; strict comparison is relaxed",
                            )
                        )

                    try:
                        extracted_nt = f.extract(rec_seq)
                        offset = f.codon_start - 1
                        coding_nt = extracted_nt[offset:]

                        if not is_partial and len(coding_nt) % 3 != 0:
                            findings.append(
                                ValidationFinding(
                                    severity="WARNING",
                                    code="LENGTH_NOT_DIVISIBLE_BY_THREE",
                                    record_id=rec.id,
                                    locus_tag=tag,
                                    feature_type="CDS",
                                    coordinates=coords_str,
                                    message=f"Coding sequence length ({len(coding_nt)} bp) not divisible by 3",
                                )
                            )

                        # Translate only with a known table.  A silent fallback
                        # to table 11 can manufacture a false validation result.
                        table_id = _translation_table_id(f)
                        if table_id is None:
                            findings.append(
                                ValidationFinding(
                                    severity="WARNING",
                                    code="UNKNOWN_TRANSLATION_TABLE",
                                    record_id=rec.id,
                                    locus_tag=tag,
                                    feature_type="CDS",
                                    coordinates=coords_str,
                                    message=f"Unknown /transl_table value {f.get_qual('transl_table')!r}; strict comparison skipped",
                                )
                            )
                            continue
                        computed_aa = str(
                            coding_nt.translate(table=table_id, to_stop=False)
                        )

                        # Strip terminal stop for comparison
                        computed_trimmed = computed_aa.rstrip("*")
                        annotated_trans = f.translation.strip()

                        if (
                            annotated_trans
                            and computed_trimmed != annotated_trans
                            and not exceptional
                        ):
                            # If first amino acid is alternative start codon (e.g. TTG/GTG translated as M)
                            if (
                                len(computed_trimmed) == len(annotated_trans)
                                and annotated_trans[0] == "M"
                                and computed_trimmed[1:] == annotated_trans[1:]
                            ):
                                pass  # Standard alternative bacterial start codon translation
                            elif not is_partial:
                                findings.append(
                                    ValidationFinding(
                                        severity="WARNING",
                                        code="CDS_TRANSLATION_MISMATCH",
                                        record_id=rec.id,
                                        locus_tag=tag,
                                        feature_type="CDS",
                                        coordinates=coords_str,
                                        message=(
                                            f"Computed translation ({len(computed_trimmed)} aa) differs from "
                                            f"/translation ({len(annotated_trans)} aa)"
                                        ),
                                    )
                                )

                        # Internal stop codon check
                        if (
                            "*" in computed_trimmed
                            and not is_partial
                            and not exceptional
                        ):
                            findings.append(
                                ValidationFinding(
                                    severity="ERROR",
                                    code="INTERNAL_STOP_CODON",
                                    record_id=rec.id,
                                    locus_tag=tag,
                                    feature_type="CDS",
                                    coordinates=coords_str,
                                    message="Computed translation contains internal stop codon(s)",
                                )
                            )

                    except (IndexError, KeyError, TypeError, ValueError) as err:
                        findings.append(
                            ValidationFinding(
                                severity="INFO",
                                code="TRANSLATION_EXTRACTION_ERROR",
                                record_id=rec.id,
                                locus_tag=tag,
                                feature_type="CDS",
                                coordinates=coords_str,
                                message=f"Could not verify translation: {err}",
                            )
                        )

    severity_order = {"ERROR": 0, "WARNING": 1, "INFO": 2}
    findings.sort(
        key=lambda finding: (
            finding.record_id,
            finding.coordinates,
            severity_order.get(finding.severity, 99),
            finding.code,
            finding.message,
        )
    )
    return ValidationReport(
        source=_source_label(doc, filepath),
        record_count=len(doc.records),
        feature_count=len(all_features),
        total_length=doc.total_length,
        type_counts=dict(sorted(type_counts.items())),
        strand_counts=dict(sorted(strand_counts.items())),
        cds_count=len(cdss),
        named_gene_count=len(named_genes),
        hypothetical_count=len(hypothetical),
        locus_tag_count=len(locus_tags),
        unique_locus_tag_count=len(unique_tags),
        findings=tuple(findings),
    )


def render_validation_json(report: ValidationReport) -> str:
    # Keep the legacy JSON top-level array in v0.8.5.
    return json.dumps(
        [asdict(finding) for finding in report.findings],
        ensure_ascii=False,
        allow_nan=False,
        sort_keys=True,
        indent=2,
    ) + "\n"


def render_validation_tsv(report: ValidationReport) -> str:
    output = io.StringIO(newline="")
    output.write("severity\tcode\trecord_id\tlocus_tag\tfeature_type\tcoordinates\tmessage\n")
    for finding in report.findings:
        values = (
            finding.severity,
            finding.code,
            finding.record_id,
            finding.locus_tag,
            finding.feature_type,
            finding.coordinates,
            finding.message,
        )
        output.write("\t".join(str(value).replace("\t", " ") for value in values) + "\n")
    return output.getvalue()


def render_validation_text(report: ValidationReport, *, max_display: int = 25) -> str:
    counts = report.severity_counts
    lines = [
        "=" * 70,
        "  GENBANK FEATURE TABLE -- STRUCTURAL & BIOLOGICAL REPORT",
        "=" * 70,
        f"  File                  : {report.source}",
        f"  Total records         : {report.record_count}",
        f"  Total features        : {report.feature_count}",
        f"  Total genome length   : {report.total_length:,} bp",
        "",
        "-- Feature type counts --",
    ]
    lines.extend(f"  {feature_type:20s}  {count}" for feature_type, count in report.type_counts.items())
    lines.extend(["", "-- Strand distribution --"])
    lines.extend(f"  {strand}  {count}" for strand, count in report.strand_counts.items())
    lines.extend(
        [
            "",
            "-- CDS statistics --",
            f"  Count                 : {report.cds_count}",
            f"  Named genes           : {report.named_gene_count}",
            f"  Hypothetical / DUF    : {report.hypothetical_count}",
            "",
            "-- Locus tags --",
            f"  Total occurrences     : {report.locus_tag_count}",
            f"  Unique tags           : {report.unique_locus_tag_count}",
            "",
            "-- Validation findings --",
            f"  Errors: {counts['ERROR']} | Warnings: {counts['WARNING']} | Info: {counts['INFO']}",
        ]
    )
    if report.findings:
        for finding in report.findings[:max_display]:
            lines.append(
                f"  [{finding.severity}] {finding.code} ({finding.coordinates}) - {finding.message}"
            )
        if len(report.findings) > max_display:
            lines.append(f"  ... and {len(report.findings) - max_display} more findings.")
    else:
        lines.append("  [OK] Clean: No structural or translation abnormalities detected.")
    lines.append("=" * 70)
    return "\n".join(lines) + "\n"


def validate(filepath: str | Path, json_mode: bool = False) -> list[ValidationFinding]:
    """Compatibility wrapper retaining the historical printing behavior."""

    report = build_validation_report(filepath)
    print(render_validation_json(report) if json_mode else render_validation_text(report), end="")
    return list(report.findings)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Validate GenBank feature table structure and biological semantics."
    )
    parser.add_argument("input", help="Input GenBank file")
    parser.add_argument(
        "--json", action="store_true", help="Output findings in JSON format"
    )
    args = parser.parse_args()
    validate(args.input, json_mode=args.json)


if __name__ == "__main__":
    main()
