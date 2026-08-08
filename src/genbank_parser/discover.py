"""Annotation-based mobilome and xenobiotic discovery with declarative rules."""

from __future__ import annotations

import argparse
import collections
import csv
import json
import sys
from importlib import resources
from pathlib import Path
from typing import Any

import yaml  # type: ignore[import-untyped]

from .io import read_genbank
from .operons import find_operon_pairs


def _load_packaged_rules(ruleset: str) -> list[dict[str, Any]]:
    if ruleset not in {"mobilome", "xenobiotics"}:
        raise ValueError(f"Unknown packaged ruleset {ruleset!r}")
    resource = resources.files("genbank_parser.rulesets").joinpath(f"{ruleset}.yaml")
    loaded = yaml.safe_load(resource.read_text(encoding="utf-8"))
    if not isinstance(loaded, list):
        raise TypeError(f"Packaged ruleset {ruleset!r} must contain a list")
    return loaded


DEFAULT_MOBILOME_RULES = _load_packaged_rules("mobilome")


def _load_rules(
    ruleset: str,
    rules_file: str | Path | None,
) -> list[dict[str, Any]]:
    if rules_file is None:
        return _load_packaged_rules(ruleset)
    path = Path(rules_file)
    if not path.exists():
        raise FileNotFoundError(f"Rules file not found: {path}")
    with path.open(encoding="utf-8") as handle:
        loaded = (
            yaml.safe_load(handle)
            if path.suffix.lower() in {".yaml", ".yml"}
            else json.load(handle)
        )
    if not isinstance(loaded, list):
        raise TypeError("Rules file must contain a list of rule objects")
    return loaded


def discover_clusters(
    filepath: str | Path,
    cluster_gap: int = 5000,
    operon_gap: int = 150,
    min_weight: int = 1,
    format_type: str = "text",
    rules_file: str | Path | None = None,
    ruleset: str = "mobilome",
) -> dict[str, Any]:
    """Find annotation-rule hits and cluster them by contig proximity.

    This is intentionally an annotation scanner.  It does not infer dark
    matter, orthology, or sequence-level mobilome evidence.
    """
    if cluster_gap < 0 or operon_gap < 0:
        raise ValueError("cluster_gap and operon_gap must be non-negative")
    if format_type not in {"text", "json", "tsv"}:
        raise ValueError(f"Unsupported discovery format: {format_type}")

    doc = read_genbank(filepath)
    rules = _load_rules(ruleset, rules_file)
    hits: list[dict[str, Any]] = []
    operon_links: list[dict[str, Any]] = []

    for rec in doc.records:
        for feature in rec.cds_features:
            gene = (feature.gene or "").casefold()
            product = (feature.product or "").casefold()
            notes = " ".join(feature.qualifiers.get("note", [])).casefold()
            matched_rules: list[dict[str, Any]] = []
            total_weight = 0
            for rule in rules:
                for term in rule.get("terms", []):
                    term_text = str(term)
                    term_lower = term_text.casefold()
                    if (
                        term_lower in gene
                        or term_lower in product
                        or term_lower in notes
                    ):
                        weight = int(rule.get("weight", 1))
                        matched_rules.append(
                            {"rule": rule["id"], "term": term_text, "weight": weight}
                        )
                        total_weight += weight
                        break
            if total_weight >= min_weight:
                hits.append(
                    {
                        "contig": rec.id,
                        "locus_tag": feature.locus_tag or "-",
                        "gene": feature.gene or "-",
                        "product": feature.product or "-",
                        "start": feature.start,
                        "end": feature.end,
                        "strand": feature.strand_symbol,
                        "weight": total_weight,
                        "matches": matched_rules,
                        "feature_index": feature.feature_index,
                    }
                )

        # ``operon_gap`` now has an observable effect: report same-strand
        # proximity links among the discovered hits.
        hit_index = {
            hit["feature_index"]: hit for hit in hits if hit["contig"] == rec.id
        }
        for first, second, gap in find_operon_pairs(rec.features, max_gap=operon_gap):
            if first.feature_index in hit_index and second.feature_index in hit_index:
                operon_links.append(
                    {
                        "contig": rec.id,
                        "first_locus_tag": first.locus_tag or "-",
                        "second_locus_tag": second.locus_tag or "-",
                        "first_feature_index": first.feature_index,
                        "second_feature_index": second.feature_index,
                        "gap": gap,
                    }
                )

    clusters: list[list[dict[str, Any]]] = []
    hits_by_contig: dict[str, list[dict[str, Any]]] = collections.defaultdict(list)
    for hit in hits:
        hits_by_contig[hit["contig"]].append(hit)
    for contig_hits in hits_by_contig.values():
        sorted_hits = sorted(contig_hits, key=lambda hit: hit["start"])
        current = [sorted_hits[0]]
        for hit in sorted_hits[1:]:
            if hit["start"] - current[-1]["end"] <= cluster_gap:
                current.append(hit)
            else:
                if len(current) >= 2 or sum(item["weight"] for item in current) >= 4:
                    clusters.append(current)
                current = [hit]
        if len(current) >= 2 or sum(item["weight"] for item in current) >= 4:
            clusters.append(current)

    islands: list[dict[str, Any]] = []
    for cluster_id, cluster in enumerate(clusters, 1):
        cluster_indices = {feature["feature_index"] for feature in cluster}
        links = [
            link
            for link in operon_links
            if link["first_feature_index"] in cluster_indices
            and link["second_feature_index"] in cluster_indices
        ]
        islands.append(
            {
                "cluster_id": f"island_{cluster_id}",
                "contig": cluster[0]["contig"],
                "start": cluster[0]["start"],
                "end": cluster[-1]["end"],
                "span": cluster[-1]["end"] - cluster[0]["start"] + 1,
                "genes": len(cluster),
                "total_weight": sum(item["weight"] for item in cluster),
                "operon_links": links,
                "features": cluster,
            }
        )

    result: dict[str, Any] = {
        "file": str(filepath),
        "ruleset": ruleset if rules_file is None else str(rules_file),
        "total_hits": len(hits),
        "total_islands": len(islands),
        "operon_gap": operon_gap,
        "hits": hits,
        "operon_links": operon_links,
        "islands": islands,
    }

    if format_type == "json":
        print(json.dumps(result, indent=2))
        return result
    if format_type == "tsv":
        writer = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "record",
                "cluster_id",
                "cluster_start",
                "cluster_end",
                "locus_tag",
                "gene",
                "product",
                "strand",
                "rule_id",
                "matched_term",
                "weight",
                "total_feature_score",
            ]
        )
        cluster_by_feature = {
            feature["feature_index"]: island
            for island in islands
            for feature in island["features"]
        }
        for hit in hits:
            island = cluster_by_feature.get(hit["feature_index"])
            matches = hit["matches"] or [{"rule": "", "term": "", "weight": ""}]
            for match in matches:
                writer.writerow(
                    [
                        hit["contig"],
                        island["cluster_id"] if island else "",
                        island["start"] if island else "",
                        island["end"] if island else "",
                        hit["locus_tag"],
                        hit["gene"],
                        hit["product"],
                        hit["strand"],
                        match["rule"],
                        match["term"],
                        match["weight"],
                        hit["weight"],
                    ]
                )
        return result

    print("=" * 70)
    print("  GENOMIC ANNOTATION DISCOVERY SCANNER")
    print("=" * 70)
    print(f"  File             : {filepath}")
    print(f"  Ruleset          : {result['ruleset']}")
    print(f"  Individual hits  : {len(hits)}")
    print(f"  Clustered islands: {len(islands)} (gap <= {cluster_gap:,} bp)")
    print(f"  Operon links     : {len(operon_links)} (gap <= {operon_gap:,} bp)")
    if islands:
        print("\n-- Discovered Annotation Islands --")
        for island in islands:
            print(
                f"  {island['cluster_id']}: {island['contig']}:{island['start']:,}..{island['end']:,} "
                f"({island['span']:,} bp, {island['genes']} genes, score: {island['total_weight']})"
            )
            for feature in island["features"]:
                matched = ", ".join(
                    f"{item['rule']}:{item['term']}" for item in feature["matches"]
                )
                print(
                    f"    [{feature['strand']}] {feature['locus_tag']:16s} {feature['gene']:8s} "
                    f"{feature['product'][:35]:35s} (matches: {matched})"
                )
    else:
        print("\n  (No annotation islands detected at chosen threshold)")
    print("=" * 70)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Mine annotation-supported mobilome or xenobiotic islands."
    )
    parser.add_argument("input", help="Input GenBank file")
    parser.add_argument("--cluster-gap", type=int, default=5000)
    parser.add_argument("--operon-gap", type=int, default=150)
    parser.add_argument("--min-weight", type=int, default=1)
    parser.add_argument("--format", choices=["text", "json", "tsv"], default="text")
    parser.add_argument(
        "--ruleset", choices=["mobilome", "xenobiotics"], default="mobilome"
    )
    parser.add_argument(
        "--rules", dest="rules_file", help="Custom YAML/JSON ruleset file"
    )
    args = parser.parse_args()
    discover_clusters(
        args.input,
        cluster_gap=args.cluster_gap,
        operon_gap=args.operon_gap,
        min_weight=args.min_weight,
        format_type=args.format,
        rules_file=args.rules_file,
        ruleset=args.ruleset,
    )


if __name__ == "__main__":
    main()
