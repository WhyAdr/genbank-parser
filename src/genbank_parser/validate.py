"""GenBank structural and biological semantics validator."""

from __future__ import annotations

import argparse
import collections
import json
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

from Bio.Data import CodonTable

from .io import read_genbank
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


def validate(filepath: str | Path, json_mode: bool = False) -> list[ValidationFinding]:
    doc = read_genbank(filepath)
    if len(doc.records) == 0:
        print("ERROR: No records parsed. Check file format.", file=sys.stderr)
        sys.exit(1)

    all_features = doc.all_features
    if not all_features:
        print("ERROR: No features parsed. Check file format.", file=sys.stderr)
        sys.exit(1)

    type_counts = collections.Counter(f.type for f in all_features)
    strand_counts = collections.Counter(f.strand_symbol for f in all_features)
    all_quals: set[str] = set()
    for f in all_features:
        all_quals.update(f.qualifiers.keys())

    cdss = [f for f in all_features if f.type == "CDS"]

    # CDS lengths
    cds_lens = [f.length for f in cdss]

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

    if json_mode:
        print(json.dumps([asdict(f) for f in findings], indent=2))
        return findings

    # Pretty-print report
    print("=" * 70)
    print("  GENBANK FEATURE TABLE -- STRUCTURAL & BIOLOGICAL REPORT")
    print("=" * 70)
    print(f"  File                  : {filepath}")
    print(f"  Total records         : {len(doc.records)}")
    print(f"  Total features        : {len(all_features)}")
    print(f"  Total genome length   : {doc.total_length:,} bp")
    print()
    print("-- Feature type counts --")
    for ft, c in type_counts.most_common():
        print(f"  {ft:20s}  {c}")
    print()
    print("-- Strand distribution --")
    for s, c in strand_counts.items():
        print(f"  {s}  {c}")
    print()
    print("-- CDS statistics --")
    print(f"  Count                 : {len(cdss)}")
    if cds_lens:
        print(f"  Min length            : {min(cds_lens):,} bp")
        print(f"  Max length            : {max(cds_lens):,} bp")
        print(f"  Mean length           : {sum(cds_lens) / len(cds_lens):,.0f} bp")
        print(f"  Median length         : {sorted(cds_lens)[len(cds_lens) // 2]:,} bp")
    print(f"  Named genes           : {len(named_genes)}")
    print(f"  Hypothetical / DUF    : {len(hypothetical)}")
    print()
    print("-- Locus tags --")
    print(f"  Total occurrences     : {len(locus_tags)}")
    print(f"  Unique tags           : {len(unique_tags)}")
    if unique_tags:
        tags_sorted = sorted(unique_tags)
        print(f"  First tag             : {tags_sorted[0]}")
        print(f"  Last tag              : {tags_sorted[-1]}")
    print()
    print("-- Validation findings --")
    error_cnt = sum(1 for finding in findings if finding.severity == "ERROR")
    warn_cnt = sum(1 for finding in findings if finding.severity == "WARNING")
    info_cnt = sum(1 for finding in findings if finding.severity == "INFO")
    print(f"  Errors: {error_cnt} | Warnings: {warn_cnt} | Info: {info_cnt}")
    if findings:
        for finding in findings[:25]:
            print(
                f"  [{finding.severity}] {finding.code} ({finding.coordinates}) - {finding.message}"
            )
        if len(findings) > 25:
            print(f"  ... and {len(findings) - 25} more findings.")
    else:
        print("  [OK] Clean: No structural or translation abnormalities detected.")
    print("=" * 70)

    return findings


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
