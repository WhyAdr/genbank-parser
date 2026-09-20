"""Evidence-preserving multi-genome annotation comparison."""

from __future__ import annotations

import argparse
import csv
import io
import json
import os
import re
from pathlib import Path
from typing import Any

import yaml  # type: ignore[import-untyped]

from .io import extract_xrefs, read_genbank

EVIDENCE_COLUMNS = (
    "genome",
    "path",
    "target",
    "target_kind",
    "record",
    "feature_index",
    "locus_tag",
    "gene",
    "product",
    "start",
    "end",
    "strand",
    "matched_field",
    "matched_value",
    "match_mode",
)


def _word_in(needle: str, haystack: str) -> bool:
    """Match a product marker as a token, avoiding substring false positives."""

    return bool(re.search(r"(?<![\w-])" + re.escape(needle) + r"(?![\w-])", haystack, re.IGNORECASE))


def discover_genome_files(genome_inputs: list[str | Path] | tuple[str | Path, ...]) -> list[Path]:
    paths: list[Path] = []
    extensions = {".gb", ".gbk", ".gbff"}
    for raw in genome_inputs:
        path = Path(raw)
        if path.is_file() and (path.suffix.casefold() in extensions or path.name.casefold().endswith(".gz")):
            paths.append(path)
        elif path.is_dir():
            paths.extend(
                candidate
                for candidate in sorted(path.rglob("*"), key=lambda item: item.as_posix().casefold())
                if candidate.is_file()
                and (
                    candidate.suffix.casefold() in extensions
                    or candidate.name.casefold().endswith((".gb.gz", ".gbk.gz", ".gbff.gz"))
                )
            )
    unique: dict[str, Path] = {}
    for path in paths:
        unique[str(path.resolve()).casefold()] = path
    return [unique[key] for key in sorted(unique)]


def _sample_ids(paths: list[Path]) -> dict[Path, str]:
    if not paths:
        return {}
    try:
        # ``commonpath`` is more reliable than Path.parents for mixed roots.
        common_text = os.path.commonpath([str(path.resolve()) for path in paths])
        common = Path(common_text)
        if common.is_file():
            common = common.parent
    except (ValueError, OSError):
        common = paths[0].parent
    result: dict[Path, str] = {}
    seen: set[str] = set()
    for path in paths:
        relative = path.resolve().relative_to(common) if path.resolve().is_relative_to(common) else path.name
        sample = relative.as_posix()
        for suffix in (".gbff.gz", ".gbff", ".gbk", ".gb", ".gz"):
            if sample.casefold().endswith(suffix):
                sample = sample[: -len(suffix)]
                break
        sample = sample.replace("/", "__") or path.stem
        if sample in seen:
            raise ValueError(f"duplicate derived genome sample ID: {sample}")
        seen.add(sample)
        result[path] = sample
    return result


def _target_definitions(targets: list[str] | None, targets_file: str | Path | None) -> list[dict[str, Any]]:
    if (targets is None) == (targets_file is None):
        raise ValueError("provide exactly one of --targets and --targets-file")
    if targets is not None:
        values = [value.strip() for value in targets if value.strip()]
        if not values:
            raise ValueError("--targets must contain at least one target")
        if len(set(values)) != len(values):
            raise ValueError("--targets must not contain duplicate target IDs")
        return [{"id": value, "legacy": value} for value in values]
    path = Path(targets_file)  # type: ignore[arg-type]
    try:
        with path.open(encoding="utf-8") as handle:
            payload = json.load(handle) if path.suffix.casefold() == ".json" else yaml.safe_load(handle)
    except (OSError, ValueError, yaml.YAMLError) as exc:
        raise ValueError(f"could not read targets file {path}: {exc}") from exc
    if not isinstance(payload, dict) or payload.get("schema_version") != 1 or not isinstance(payload.get("targets"), list):
        raise ValueError("targets file must contain schema_version: 1 and a targets list")
    definitions: list[dict[str, Any]] = []
    target_ids: set[str] = set()
    for target in payload["targets"]:
        if not isinstance(target, dict) or not isinstance(target.get("id"), str):
            raise ValueError("each target must contain a string id")  # noqa: TRY004
        target_id = target["id"].strip()
        if not target_id:
            raise ValueError("target ids must not be empty")
        if target_id in target_ids:
            raise ValueError(f"duplicate target id {target_id!r}")
        target_ids.add(target_id)
        clauses = target.get("any")
        if not isinstance(clauses, list) or not clauses:
            raise ValueError(f"target {target_id!r} must contain a non-empty any list")
        normalized: list[dict[str, str]] = []
        for clause in clauses:
            if not isinstance(clause, dict) or not isinstance(clause.get("field"), str) or "value" not in clause:
                raise ValueError(f"target {target_id!r} contains an invalid clause")
            mode = str(clause.get("mode", "exact"))
            if mode not in {"exact", "casefold_exact", "contains", "token", "regex"}:
                raise ValueError(f"unsupported target matching mode {mode!r}")
            if clause["field"] not in {
                "record",
                "type",
                "locus_tag",
                "gene",
                "product",
                "protein_id",
                "ko",
                "kegg_ko",
                "ec",
                "ec_number",
                "cog",
                "pfam",
                "rfam",
                "go",
                "db_xref",
            }:
                raise ValueError(f"unsupported target field {clause['field']!r}")
            value = str(clause["value"])
            if mode == "regex":
                try:
                    re.compile(value)
                except re.error as exc:
                    raise ValueError(f"invalid target regex {value!r}: {exc}") from exc
            normalized.append({"field": clause["field"], "mode": mode, "value": value})
        definitions.append({"id": target_id, "clauses": normalized})
    return definitions


def _feature_values(feature: Any) -> dict[str, list[str]]:
    xrefs = extract_xrefs(feature)
    return {
        "record": [feature.record_id],
        "type": [feature.type],
        "locus_tag": [feature.locus_tag] if feature.locus_tag else [],
        "gene": [feature.gene] if feature.gene else [],
        "product": [feature.product] if feature.product else [],
        "protein_id": [feature.protein_id] if feature.protein_id else [],
        "ko": list(xrefs["kegg_kos"]),
        "kegg_ko": list(xrefs["kegg_kos"]),
        "ec": list(xrefs["ec_numbers"]),
        "ec_number": list(xrefs["ec_numbers"]),
        "cog": list(xrefs["cog_ids"]),
        "pfam": list(xrefs["pfam"]),
        "rfam": list(xrefs["rfam"]),
        "go": list(xrefs["go_terms"]),
        "db_xref": list(feature.qualifiers.get("db_xref", [])),
    }


def _match_clause(values: list[str], mode: str, expected: str) -> list[str]:
    matches: list[str] = []
    for value in values:
        if mode == "exact" and value == expected or mode == "casefold_exact" and value.casefold() == expected.casefold() or mode == "contains" and expected.casefold() in value.casefold() or mode == "token" and _word_in(expected, value) or mode == "regex" and re.search(expected, value) is not None:
            matches.append(value)
    return matches


def _legacy_clauses(value: str) -> list[dict[str, str]]:
    return [
        {"field": "gene", "mode": "casefold_exact", "value": value},
        {"field": "product", "mode": "token", "value": value},
        {"field": "ko", "mode": "casefold_exact", "value": value},
        {"field": "ec", "mode": "casefold_exact", "value": value},
        {"field": "cog", "mode": "casefold_exact", "value": value},
    ]


def compare_data(
    genome_inputs: list[str | Path] | tuple[str | Path, ...],
    *,
    targets: list[str] | None = None,
    targets_file: str | Path | None = None,
    mode: str = "count",
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    if mode not in {"count", "presence", "status"}:
        raise ValueError("mode must be count, presence, or status")
    definitions = _target_definitions(targets, targets_file)
    paths = discover_genome_files(genome_inputs)
    if not paths:
        raise ValueError("no compatible GenBank files found in input path(s)")
    samples = _sample_ids(paths)
    matrix: list[dict[str, Any]] = []
    evidence: list[dict[str, Any]] = []
    for path in paths:
        document = read_genbank(path)
        counts = {definition["id"]: 0 for definition in definitions}
        for feature in document.all_features:
            if feature.type != "CDS":
                continue
            values = _feature_values(feature)
            for definition in definitions:
                clauses = definition.get("clauses") or _legacy_clauses(definition["legacy"])
                matched: list[tuple[str, str, str]] = []
                for clause in clauses:
                    field = clause["field"]
                    clause_matches = _match_clause(values.get(field, []), clause["mode"], clause["value"])
                    matched.extend((field, value, clause["mode"]) for value in clause_matches)
                if not matched:
                    continue
                counts[definition["id"]] += 1
                evidence.append(
                    {
                        "genome": samples[path],
                        "path": str(path),
                        "target": definition["id"],
                        "target_kind": "typed" if definition.get("clauses") else "legacy",
                        "record": feature.record_id,
                        "feature_index": feature.feature_index,
                        "locus_tag": feature.locus_tag or "",
                        "gene": feature.gene or "",
                        "product": feature.product or "",
                        "start": feature.start,
                        "end": feature.end,
                        "strand": feature.strand_symbol,
                        "matched_field": ";".join(dict.fromkeys(item[0] for item in matched)),
                        "matched_value": ";".join(dict.fromkeys(item[1] for item in matched)),
                        "match_mode": ";".join(dict.fromkeys(item[2] for item in matched)),
                    }
                )
        row: dict[str, Any] = {"Genome": samples[path], "Path": str(path)}
        for definition in definitions:
            count = counts[definition["id"]]
            row[definition["id"]] = (
                count
                if mode == "count"
                else int(bool(count))
                if mode == "presence"
                else "ABSENT"
                if count == 0
                else "SINGLE_COPY"
                if count == 1
                else "MULTI_COPY"
            )
        matrix.append(row)
    return matrix, evidence


def render_compare(
    matrix: list[dict[str, Any]],
    evidence: list[dict[str, Any]],
    *,
    format_type: str = "text",
    mode: str = "count",
) -> str:
    targets = [key for key in matrix[0] if key not in {"Genome", "Path"}] if matrix else []
    if format_type == "json":
        return json.dumps({"mode": mode, "targets": targets, "matrix": matrix, "evidence": evidence}, ensure_ascii=False, sort_keys=True, indent=2) + "\n"
    if format_type == "tsv":
        output = io.StringIO(newline="")
        writer = csv.DictWriter(output, fieldnames=["Genome", "Path", *targets], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(matrix)
        return output.getvalue()
    lines = [
        "=" * 80,
        "  MULTI-GENOME ANNOTATION/XREF MARKER SCANNER",
        "=" * 80,
        f"  Mode: {mode}",
        "  Genome / Isolate  " + "  ".join(targets),
    ]
    for row in matrix:
        lines.append("  " + row["Genome"] + "  " + "  ".join(str(row[target]) for target in targets))
    lines.append("=" * 80)
    return "\n".join(lines) + "\n"


def render_evidence_tsv(evidence: list[dict[str, Any]]) -> str:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=EVIDENCE_COLUMNS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(evidence)
    return output.getvalue()


def compare_genomes(
    genome_inputs: list[str | Path] | tuple[str | Path, ...],
    targets: list[str] | tuple[str, ...] | None = None,
    output_path: str | Path | None = None,
    *,
    targets_file: str | Path | None = None,
    mode: str = "count",
    evidence_output: str | Path | None = None,
    format_type: str = "text",
) -> list[dict[str, Any]]:
    """Compatibility wrapper returning matrix rows and optionally printing."""

    matrix, evidence = compare_data(
        genome_inputs,
        targets=list(targets) if targets is not None else None,
        targets_file=targets_file,
        mode=mode,
    )
    rendered = render_compare(matrix, evidence, format_type=format_type, mode=mode)
    if output_path:
        Path(output_path).write_text(rendered, encoding="utf-8", newline="")
    else:
        print(rendered, end="")
    if evidence_output:
        Path(evidence_output).write_text(render_evidence_tsv(evidence), encoding="utf-8", newline="")
    return matrix


def main() -> None:
    parser = argparse.ArgumentParser(description="Multi-genome marker presence/absence matrix.")
    parser.add_argument("genomes", nargs="+", help="GenBank files or directories")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--targets", help="Comma-separated marker names")
    group.add_argument("--targets-file", help="YAML/JSON target definition")
    parser.add_argument("--mode", choices=["count", "presence", "status"], default="count")
    parser.add_argument("--format", choices=["text", "tsv", "json"], default="text")
    parser.add_argument("--output", help="Output matrix path")
    parser.add_argument("--evidence-output", help="Evidence TSV path")
    args = parser.parse_args()
    targets = [item.strip() for item in args.targets.split(",") if item.strip()] if args.targets else None
    compare_genomes(args.genomes, targets, args.output, targets_file=args.targets_file, mode=args.mode, evidence_output=args.evidence_output, format_type=args.format)


if __name__ == "__main__":
    main()
