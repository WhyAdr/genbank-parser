#!/usr/bin/env python3
"""MEOR & Biosurfactant Biosynthesis Discovery Engine.

Parses GenBank feature table files (.gbff, .gbk, .gb) computationally and evaluates
isolate genomes for Microbial Enhanced Oil Recovery (MEOR) potential, hydrocarbon
degradation pathways (aliphatic, aromatic, anaerobic), and biosurfactant/bio-emulsifier
biosynthetic gene clusters (glycolipids, lipopeptides, polymeric EPS).

Incorporates marker definitions and cross-references synthesized from CANT-HYD,
HMDB (Hydrocarbon Monooxygenase Gene Database), and HADEG databases.

Usage:
    python scripts/genbank_meor.py INPUT.gbff
    python scripts/genbank_meor.py INPUT.gbff --format json
    python scripts/genbank_meor.py INPUT.gbff --format tsv
    python scripts/genbank_meor.py INPUT.gbff --min-weight 2 --max-gap 300
"""

import sys
import os
import argparse
import collections
import json
import csv
import re
from genbank_parser import parse_features, get_qual, extract_xrefs

# ---------------------------------------------------------------------------
# MEOR & Biosurfactant Marker Database (CANT-HYD + HMDB + HADEG)
# ---------------------------------------------------------------------------

MEOR_CATEGORIES = {
    "cat1_short_alkane": {
        "title": "Short-Chain Alkanes & Gaseous Hydrocarbons (C1-C4)",
        "role": "Methane/propane/butane oxidation & gas utilization",
        "markers": [
            {
                "id": "pmoABC",
                "name": "Particulate Methane Monooxygenase (pmoABC)",
                "genes": [r"\bpmo[ABC]\b"],
                "ecs": ["1.14.18.3"],
                "kos": ["K10968", "K10969", "K10970"],
                "products": [r"particulate methane monooxygenase", r"methane monooxygenase subunit"]
            },
            {
                "id": "mmoXYZ",
                "name": "Soluble Methane Monooxygenase (mmoXYZ)",
                "genes": [r"\bmmo[X-Z]\b", r"\bmmo[B-D]\b"],
                "ecs": ["1.14.13.25"],
                "kos": ["K16157", "K16158", "K16159", "K16160", "K16161", "K16162"],
                "products": [r"soluble methane monooxygenase", r"methane monooxygenase hydroxylase"]
            },
            {
                "id": "prmABCD",
                "name": "Propane Monooxygenase (prmABCD / bmo)",
                "genes": [r"\bprm[A-D]\b", r"\bbmo[X-Z]\b"],
                "ecs": ["1.14.13.-"],
                "kos": ["K22473", "K22474", "K22475"],
                "products": [r"propane monooxygenase", r"butane monooxygenase"]
            }
        ]
    },
    "cat2_medium_alkane": {
        "title": "Medium-Chain n-Alkanes (C5-C16)",
        "role": "Liquid alkane degradation & oil viscosity reduction",
        "markers": [
            {
                "id": "alkB",
                "name": "Alkane 1-monooxygenase (alkB / alkB1 / alkB2)",
                "genes": [r"\balkB\d?\b"],
                "ecs": ["1.14.15.3"],
                "kos": ["K00496"],
                "products": [r"alkane 1-monooxygenase", r"alkane hydroxylase", r"alkane monooxygenase"]
            },
            {
                "id": "CYP153",
                "name": "Cytochrome P450 Alkane Hydroxylase (CYP153)",
                "genes": [r"\bcyp153\b", r"\bcyp153A\d?\b"],
                "ecs": ["1.14.15.-"],
                "kos": ["K00496"],
                "products": [r"cytochrome P450.*alkane", r"CYP153"]
            },
            {
                "id": "rubAB",
                "name": "Rubredoxin & Rubredoxin Reductase (rubA/rubB)",
                "genes": [r"\brubA\d?\b", r"\brubB\b"],
                "ecs": ["1.18.1.1"],
                "kos": ["K00389", "K03820"],
                "products": [r"rubredoxin", r"rubredoxin--NAD\+ reductase", r"rubredoxin reductase"]
            },
            {
                "id": "alkJ",
                "name": "Alcohol Dehydrogenase (alkJ / Alkane pathway)",
                "genes": [r"\balkJ\b"],
                "ecs": ["1.1.1.1", "1.1.99.8"],
                "kos": ["K00001", "K13953"],
                "products": [r"alkane alcohol dehydrogenase", r"fatty alcohol dehydrogenase"]
            },
            {
                "id": "alkH",
                "name": "Aldehyde Dehydrogenase (alkH / Alkane pathway)",
                "genes": [r"\balkH\b"],
                "ecs": ["1.2.1.3"],
                "kos": ["K00128"],
                "products": [r"fatty aldehyde dehydrogenase", r"alkane aldehyde dehydrogenase"]
            },
            {
                "id": "alkK",
                "name": "Fatty Acyl-CoA Synthetase (alkK)",
                "genes": [r"\balkK\b"],
                "ecs": ["6.2.1.3"],
                "kos": ["K01897"],
                "products": [r"acyl-CoA synthetase", r"medium-chain fatty-acid--CoA ligase"]
            }
        ]
    },
    "cat3_long_alkane": {
        "title": "Long-Chain n-Alkanes & Heavy Paraffins (C18-C36+)",
        "role": "Heavy paraffin wax degradation & pour point reduction",
        "markers": [
            {
                "id": "ladA",
                "name": "Long-chain Alkane Monooxygenase LadA (Flavoprotein C15-C36)",
                "genes": [r"\bladA\d?\b"],
                "ecs": ["1.14.14.28"],
                "kos": ["K22476"],
                "products": [r"long-chain alkane monooxygenase", r"flavoprotein monooxygenase ladA"]
            },
            {
                "id": "almA",
                "name": "Long-chain Alkane Hydroxylase AlmA (Flavin-binding C20-C32)",
                "genes": [r"\balmA\b"],
                "ecs": ["1.14.14.-"],
                "kos": ["K22477"],
                "products": [r"flavin-binding monooxygenase almA", r"long-chain alkane hydroxylase almA"]
            }
        ]
    },
    "cat4_aromatic_pah": {
        "title": "Aromatic & Polycyclic Aromatic Hydrocarbons (BTEX & PAHs)",
        "role": "Aromatic fraction breakdown & heavy crude liquefaction",
        "markers": [
            {
                "id": "tmo_tod",
                "name": "Toluene / BTEX Monooxygenase (tmoABCDEF / todC1C2BA)",
                "genes": [r"\btmo[A-F]\b", r"\btod[A-D]\b", r"\bxyl[MA]\b"],
                "ecs": ["1.14.13.-"],
                "kos": ["K15760", "K15761", "K15762"],
                "products": [r"toluene monooxygenase", r"xylene monooxygenase", r"toluene dioxygenase"]
            },
            {
                "id": "ndo_nah",
                "name": "Naphthalene / PAH Dioxygenase (ndoABC / nahAcAb / phnAc)",
                "genes": [r"\bndo[A-C]\b", r"\bnahA[cb]?\b", r"\bphnA[cb]?\b", r"\bbphA[1-4]\b"],
                "ecs": ["1.14.12.12", "1.14.12.18"],
                "kos": ["K14579", "K14580", "K14581"],
                "products": [r"naphthalene 1,2-dioxygenase", r"biphenyl 2,3-dioxygenase", r"polycyclic aromatic hydrocarbon dioxygenase"]
            },
            {
                "id": "catA_xylE",
                "name": "Central Ring-Cleavage Dioxygenases (catA / xylE)",
                "genes": [r"\bcatA\b", r"\bxylE\b"],
                "ecs": ["1.13.11.1", "1.13.11.2"],
                "kos": ["K03381", "K00446"],
                "products": [r"catechol 1,2-dioxygenase", r"catechol 2,3-dioxygenase"]
            }
        ]
    },
    "cat5_anaerobic_activation": {
        "title": "Anaerobic Hydrocarbon Activation (Reservoir Anoxic Conditions)",
        "role": "Anoxic in-situ hydrocarbon activation via fumarate addition",
        "markers": [
            {
                "id": "assA_masD",
                "name": "Alkylsuccinate Synthase / MasD (assA / masD)",
                "genes": [r"\bassA\b", r"\bmasD\b"],
                "ecs": ["4.1.99.16"],
                "kos": ["K22204"],
                "products": [r"alkylsuccinate synthase", r"1-nnet-succinate synthase"]
            },
            {
                "id": "bssABC",
                "name": "Benzylsuccinate Synthase (bssA / bssB / bssC)",
                "genes": [r"\bbss[A-C]\b"],
                "ecs": ["4.1.99.11"],
                "kos": ["K07540", "K07541", "K07542"],
                "products": [r"benzylsuccinate synthase"]
            },
            {
                "id": "nmsA",
                "name": "Naphthylmethylsuccinate Synthase (nmsA)",
                "genes": [r"\bnmsA\b"],
                "ecs": ["4.1.99.-"],
                "kos": ["K22205"],
                "products": [r"naphthylmethylsuccinate synthase"]
            }
        ]
    },
    "cat6_glycolipid_biosurfactant": {
        "title": "Glycolipid Biosurfactants (Rhamnolipids & Trehalolipids)",
        "role": "Surface tension reduction & bio-emulsification (IFT reduction)",
        "markers": [
            {
                "id": "rhlA",
                "name": "HAA Synthase (rhlA / Rhamnolipid precursor)",
                "genes": [r"\brhlA\b"],
                "ecs": ["2.3.1.266"],
                "kos": ["K13057"],
                "products": [r"rhamnolipid biosynthesis 3-hydroxyacyl-acyl carrier protein", r"HAA synthase", r"rhlA"]
            },
            {
                "id": "rhlB",
                "name": "Rhamnosyltransferase I (rhlB / Mono-rhamnolipid)",
                "genes": [r"\brhlB\b"],
                "ecs": ["2.4.1.189"],
                "kos": ["K13058"],
                "products": [r"rhamnosyltransferase chain B", r"rhamnosyltransferase 1", r"rhlB"]
            },
            {
                "id": "rhlC",
                "name": "Rhamnosyltransferase II (rhlC / Di-rhamnolipid)",
                "genes": [r"\brhlC\b"],
                "ecs": ["2.4.1.298"],
                "kos": ["K13059"],
                "products": [r"rhamnosyltransferase 2", r"rhlC"]
            },
            {
                "id": "rhlRI",
                "name": "Rhamnolipid Quorum Sensing Regulators (rhlR / rhlI)",
                "genes": [r"\brhlR\b", r"\brhlI\b"],
                "ecs": ["2.3.1.184"],
                "kos": ["K13060", "K13061"],
                "products": [r"regulatory protein RhlR", r"acyl-homoserine-lactone synthase RhlI"]
            },
            {
                "id": "trehalolipid",
                "name": "Trehalolipid Biosynthesis (sdtA / otsA / otsB)",
                "genes": [r"\bsdtA\b", r"\tshA\b", r"\botsA\b", r"\botsB\b"],
                "ecs": ["2.4.1.15", "3.1.3.12"],
                "kos": ["K00697", "K01087"],
                "products": [r"trehalose-6-phosphate synthase", r"trehalose-6-phosphatase", r"sulfetolipid"]
            }
        ]
    },
    "cat7_lipopeptide_biosurfactant": {
        "title": "Lipopeptide Biosurfactant NRPS BGCs (Surfactin & Lichenysin)",
        "role": "Potent lipopeptide biosurfactants & interfacial tension reduction",
        "markers": [
            {
                "id": "srfA",
                "name": "Surfactin Synthetase NRPS (srfA-A / srfA-B / srfA-C / srfAD)",
                "genes": [r"\bsrfA-[A-D]\b", r"\bsrfA[ABCD]\b", r"\bsrfA\b"],
                "ecs": ["2.7.7.-", "6.3.2.-"],
                "kos": ["K15666", "K15667", "K15668"],
                "products": [r"surfactin synthetase", r"surfactin synthase", r"surfactin biosynthesis"]
            },
            {
                "id": "licA_D",
                "name": "Lichenysin Synthetase NRPS (licA / licB / licC / licD)",
                "genes": [r"\blic[A-D]\b"],
                "ecs": ["2.7.7.-"],
                "kos": ["K15666"],
                "products": [r"lichenysin synthetase", r"lichenysin synthase", r"lichenysin biosynthesis"]
            },
            {
                "id": "fengycin_iturin",
                "name": "Fengycin / Iturin / Viscosin NRPS (fenA-E / ituA-C / vsnA-C)",
                "genes": [r"\bfen[A-E]\b", r"\bpps[A-E]\b", r"\bitu[A-C]\b", r"\bvsn[A-C]\b"],
                "ecs": ["2.7.7.-"],
                "kos": ["K15669", "K15670"],
                "products": [r"fengycin synthetase", r"iturin synthetase", r"plipastatin synthase", r"viscosin synthetase"]
            }
        ]
    },
    "cat8_emulsifier_acid_gas": {
        "title": "Polymeric Bio-emulsifiers & Organic Acid / Gas Production",
        "role": "Emulsification, rock dissolution, acid flooding & reservoir repressurization",
        "markers": [
            {
                "id": "wza_wzb_wzc",
                "name": "Polymeric Bio-emulsifier Secretion (wza / wzb / wzc / EPS)",
                "genes": [r"\bwz[abc]\b", r"\balmA\b"],
                "ecs": ["3.1.3.-"],
                "kos": ["K01990", "K01991"],
                "products": [r"polysaccharide export protein Wza", r"tyrosine-protein phosphatase Wzb", r"emulsan", r"alasan"]
            },
            {
                "id": "meor_acid_gas",
                "name": "MEOR Organic Acid & Gas Drivers (pflB / ackA / pta / ca)",
                "genes": [r"\bpflB\b", r"\backA\b", r"\bpta\b", r"\bcan?\b"],
                "ecs": ["2.3.1.54", "2.7.2.1", "2.3.1.8", "4.2.1.1"],
                "kos": ["K00656", "K00925", "K00625", "K01672", "K01673"],
                "products": [r"pyruvate formate-lyase", r"acetate kinase", r"phosphate acetyltransferase", r"carbonic anhydrase"]
            }
        ]
    }
}

# Key MEOR Pathways for completeness evaluation
MEOR_PATHWAYS = {
    "Medium Alkane Degradation Chain (C5-C16)": {
        "category": "cat2_medium_alkane",
        "steps": [
            {"step": "Primary Alkane Hydroxylation", "ids": ["alkB", "CYP153"]},
            {"step": "Electron Transfer (Rubredoxin System)", "ids": ["rubAB"]},
            {"step": "Alcohol Oxidation (Alcohol Dehydrogenase)", "ids": ["alkJ"]},
            {"step": "Aldehyde Oxidation (Aldehyde Dehydrogenase)", "ids": ["alkH"]},
            {"step": "Acyl-CoA Activation (Acyl-CoA Synthetase)", "ids": ["alkK"]}
        ]
    },
    "Long-Chain Paraffin Wax Degradation (C18+)": {
        "category": "cat3_long_alkane",
        "steps": [
            {"step": "Long-Chain Monooxygenase (LadA)", "ids": ["ladA"]},
            {"step": "Long-Chain Hydroxylase (AlmA)", "ids": ["almA"]}
        ]
    },
    "Rhamnolipid Biosurfactant Biosynthesis": {
        "category": "cat6_glycolipid_biosurfactant",
        "steps": [
            {"step": "HAA Fatty Acid Precursor Synthase (RhlA)", "ids": ["rhlA"]},
            {"step": "Mono-rhamnolipid Synthase (RhlB)", "ids": ["rhlB"]},
            {"step": "Di-rhamnolipid Synthase (RhlC)", "ids": ["rhlC"]}
        ]
    },
    "Lipopeptide Biosurfactant BGC (Surfactin/Lichenysin)": {
        "category": "cat7_lipopeptide_biosurfactant",
        "steps": [
            {"step": "Surfactin / Lichenysin NRPS Modules", "ids": ["srfA", "licA_D"]}
        ]
    },
    "Anaerobic Alkane Activation (AssA)": {
        "category": "cat5_anaerobic_activation",
        "steps": [
            {"step": "Alkylsuccinate Synthase catalytic alpha subunit", "ids": ["assA_masD"]}
        ]
    }
}


# ---------------------------------------------------------------------------
# Matching & Scanner Core
# ---------------------------------------------------------------------------

def match_feature_to_marker(feat, marker):
    """Evaluate a feature against a marker definition.

    Returns tuple (matched: bool, weight: int, match_reason: str).
    Weight 3 = Exact KO, EC, or Gene Symbol match (High Confidence)
    Weight 2 = Specific product keyword match (Medium Confidence)
    Weight 1 = Broad note/product keyword match (Low Confidence)
    """
    qualifiers = feat['qualifiers']
    xrefs = extract_xrefs(feat)

    gene_name = get_qual(feat, 'gene', '').strip()
    product = get_qual(feat, 'product', '').lower().strip()
    notes = " ".join(qualifiers.get('note', [])).lower()

    # Tier 1: KEGG KO match
    for ko in xrefs['kegg_kos']:
        if ko in marker.get('kos', []):
            return True, 3, f"KEGG:{ko}"

    # Tier 1: EC number match
    for ec in xrefs['ec_numbers']:
        if ec in marker.get('ecs', []):
            return True, 3, f"EC:{ec}"

    # Tier 1: Gene symbol regex match
    if gene_name:
        for g_pat in marker.get('genes', []):
            if re.search(g_pat, gene_name, re.IGNORECASE):
                return True, 3, f"gene:{gene_name}"

    # Tier 2: Product keyword regex match
    if product:
        for p_pat in marker.get('products', []):
            if re.search(p_pat, product, re.IGNORECASE):
                return True, 2, f"product:{p_pat}"

    # Tier 1/2: Note search for gene or KO/EC
    if notes:
        for ec in marker.get('ecs', []):
            if f"ec:{ec.lower()}" in notes or f"ec_number:{ec.lower()}" in notes:
                return True, 3, f"note:EC:{ec}"
        for ko in marker.get('kos', []):
            if ko.lower() in notes:
                return True, 3, f"note:KEGG:{ko}"

    return False, 0, ""


def scan_meor_features(features, min_weight=1):
    """Scan all features for MEOR & biosurfactant markers.

    Returns list of dicts containing hit details.
    """
    hits = []
    for feat in features:
        if feat['type'] not in ('CDS', 'misc_feature'):
            continue

        locus_tag = get_qual(feat, 'locus_tag', 'NO_LOCUS')
        gene = get_qual(feat, 'gene', '')
        product = get_qual(feat, 'product', '')

        for cat_key, cat_info in MEOR_CATEGORIES.items():
            for marker in cat_info['markers']:
                matched, weight, reason = match_feature_to_marker(feat, marker)
                if matched and weight >= min_weight:
                    hits.append({
                        'contig': feat['contig'],
                        'locus_tag': locus_tag,
                        'gene': gene,
                        'product': product,
                        'start': feat['start'],
                        'end': feat['end'],
                        'strand': feat['strand'],
                        'category_key': cat_key,
                        'category_title': cat_info['title'],
                        'marker_id': marker['id'],
                        'marker_name': marker['name'],
                        'weight': weight,
                        'reason': reason
                    })
                    break  # Avoid double-counting same feature in same marker
    return hits


def find_meor_clusters(hits, max_gap=200):
    """Group co-localized MEOR & biosurfactant genes into clusters/operons.

    Partitioned strictly by contig. Max gap in bp between adjacent genes.
    """
    if not hits:
        return []

    # Sort hits by contig, then start coordinate
    sorted_hits = sorted(hits, key=lambda x: (x['contig'], x['start']))

    clusters = []
    curr_cluster = [sorted_hits[0]]

    for hit in sorted_hits[1:]:
        prev_hit = curr_cluster[-1]
        if hit['contig'] == prev_hit['contig'] and (hit['start'] - prev_hit['end']) <= max_gap:
            curr_cluster.append(hit)
        else:
            if len(curr_cluster) >= 2:
                clusters.append(curr_cluster)
            curr_cluster = [hit]

    if len(curr_cluster) >= 2:
        clusters.append(curr_cluster)

    formatted_clusters = []
    for idx, cl in enumerate(clusters, 1):
        categories = sorted(list(set(h['category_key'] for h in cl)))
        locus_tags = [h['locus_tag'] for h in cl]
        genes = [h['gene'] if h['gene'] else h['marker_id'] for h in cl]
        span_bp = cl[-1]['end'] - cl[0]['start'] + 1

        formatted_clusters.append({
            'cluster_id': f"MEOR_Cluster_{idx:02d}",
            'contig': cl[0]['contig'],
            'start': cl[0]['start'],
            'end': cl[-1]['end'],
            'span_bp': span_bp,
            'gene_count': len(cl),
            'categories': categories,
            'locus_tags': locus_tags,
            'genes': genes,
            'members': cl
        })

    return formatted_clusters


def evaluate_pathways(hits):
    """Evaluate percentage completeness of predefined key MEOR pathways."""
    hit_marker_ids = set(h['marker_id'] for h in hits)
    pathway_results = []

    for pw_name, pw_info in MEOR_PATHWAYS.items():
        total_steps = len(pw_info['steps'])
        present_steps = 0
        step_details = []

        for step in pw_info['steps']:
            step_name = step['step']
            req_ids = step['ids']
            is_present = any(mid in hit_marker_ids for mid in req_ids)
            if is_present:
                present_steps += 1
            step_details.append({
                'step': step_name,
                'status': 'PRESENT' if is_present else 'ABSENT',
                'target_markers': req_ids
            })

        pct = (present_steps / total_steps * 100.0) if total_steps > 0 else 0.0
        pathway_results.append({
            'pathway': pw_name,
            'completeness_pct': round(pct, 1),
            'steps_present': present_steps,
            'total_steps': total_steps,
            'steps': step_details
        })

    return pathway_results


def generate_karyograms(features, hits, window_size=50000):
    """Generate ASCII karyogram visual density overview per contig."""
    contig_spans = collections.defaultdict(int)
    for feat in features:
        contig_spans[feat['contig']] = max(contig_spans[feat['contig']], feat['end'])

    karyograms = {}
    for contig, span in contig_spans.items():
        num_windows = max(1, (span + window_size - 1) // window_size)
        window_counts = [0] * num_windows

        contig_hits = [h for h in hits if h['contig'] == contig]
        for h in contig_hits:
            w_idx = min(num_windows - 1, (h['start'] - 1) // window_size)
            window_counts[w_idx] += 1

        ascii_bar = "".join("#" if c > 0 else "-" for c in window_counts)
        karyograms[contig] = {
            'span_bp': span,
            'window_size': window_size,
            'total_hits': len(contig_hits),
            'ascii_map': ascii_bar
        }
    return karyograms


# ---------------------------------------------------------------------------
# Main Output & CLI
# ---------------------------------------------------------------------------

def print_text_report(filepath, features, hits, clusters, pathways, karyograms):
    """Print clean human-readable console report."""
    print("=" * 80)
    print(" MEOR & BIOSURFACTANT DISCOVERY REPORT")
    print(f" Target File: {os.path.basename(filepath)}")
    print(f" Total Features Parsed: {len(features):,}")
    print(f" Total MEOR & Biosurfactant Hits: {len(hits)}")
    print(f" Co-localized BGC / Operon Candidates: {len(clusters)}")
    print("=" * 80)

    # 1. Summary by Category
    print("\n[1] CATEGORY SUMMARY & HIT COUNTS")
    print("-" * 80)
    print(f"{'Category':<52} | {'High (W=3)':<10} | {'Med (W=2)':<10} | {'Total':<6}")
    print("-" * 80)

    cat_counts = collections.defaultdict(lambda: {'w3': 0, 'w2': 0, 'w1': 0})
    for h in hits:
        w_key = f"w{h['weight']}"
        cat_counts[h['category_key']][w_key] += 1

    for cat_key, cat_info in MEOR_CATEGORIES.items():
        counts = cat_counts[cat_key]
        total = counts['w3'] + counts['w2'] + counts['w1']
        title = cat_info['title'][:50]
        print(f"{title:<52} | {counts['w3']:<10} | {counts['w2']:<10} | {total:<6}")
    print("-" * 80)

    # 2. Pathway Completeness
    print("\n[2] MEOR PATHWAY COMPLETENESS PROFILE")
    print("-" * 80)
    for pw in pathways:
        bar_len = int(pw['completeness_pct'] // 10)
        bar = "=" * bar_len + "-" * (10 - bar_len)
        print(f" {pw['pathway']:<45} [{bar}] {pw['completeness_pct']:>5.1f}% ({pw['steps_present']}/{pw['total_steps']} steps)")
        for st in pw['steps']:
            mark = "+" if st['status'] == 'PRESENT' else "-"
            print(f"    [{mark}] {st['step']}")
    print("-" * 80)

    # 3. Spatial BGC & Operon Clusters
    print("\n[3] CO-LOCALIZED MEOR OPERONS & BGC CLUSTERS (Gap <= 200 bp)")
    print("-" * 80)
    if not clusters:
        print(" No co-localized multi-gene MEOR clusters detected.")
    else:
        for cl in clusters:
            print(f" Cluster ID : {cl['cluster_id']} ({cl['contig']}:{cl['start']:,}-{cl['end']:,}, span: {cl['span_bp']:,} bp)")
            print(f" Genes ({cl['gene_count']}): {', '.join(cl['genes'])}")
            print(f" Locus Tags : {', '.join(cl['locus_tags'])}")
            print(f" Members:")
            for m in cl['members']:
                gene_str = f"({m['gene']})" if m['gene'] else ""
                prod_str = m['product'][:40] + ".." if len(m['product']) > 40 else m['product']
                print(f"   * {m['locus_tag']:<14} {m['strand']} {m['start']:>8d}-{m['end']:<8d} {m['marker_name'][:32]:<32} {gene_str} [{prod_str}]")
            print()
    print("-" * 80)

    # 4. Karyogram Density
    print("\n[4] GENOMIC DENSITY KARYOGRAMS (50 kb windows)")
    print("-" * 80)
    for contig, k_data in karyograms.items():
        print(f" Contig: {contig} ({k_data['span_bp']:,} bp, {k_data['total_hits']} hits)")
        print(f" Map: [{k_data['ascii_map']}] (# = hit present in 50kb window)")
    print("=" * 80)


def main():
    parser = argparse.ArgumentParser(
        description="MEOR & Biosurfactant Biosynthesis Discovery Engine"
    )
    parser.add_argument("input_file", help="Path to GenBank file (.gbff, .gbk, .gb)")
    parser.add_argument("--format", choices=["text", "json", "tsv"], default="text",
                        help="Output format (default: text console report)")
    parser.add_argument("--min-weight", type=int, choices=[1, 2, 3], default=1,
                        help="Minimum confidence weight (1=All, 2=Medium+, 3=High confidence only)")
    parser.add_argument("--max-gap", type=int, default=200,
                        help="Max intergenic gap (bp) for operon/BGC clustering (default: 200)")

    args = parser.parse_args()

    if not os.path.exists(args.input_file):
        print(f"Error: File '{args.input_file}' not found.", file=sys.stderr)
        sys.exit(1)

    try:
        features = parse_features(args.input_file)
    except Exception as e:
        print(f"Error parsing GenBank file '{args.input_file}': {e}", file=sys.stderr)
        sys.exit(1)

    hits = scan_meor_features(features, min_weight=args.min_weight)
    clusters = find_meor_clusters(hits, max_gap=args.max_gap)
    pathways = evaluate_pathways(hits)
    karyograms = generate_karyograms(features, hits)

    if args.format == "json":
        output = {
            "file": os.path.basename(args.input_file),
            "total_features": len(features),
            "total_hits": len(hits),
            "total_clusters": len(clusters),
            "hits": hits,
            "clusters": clusters,
            "pathways": pathways,
            "karyograms": karyograms
        }
        print(json.dumps(output, indent=2))
    elif args.format == "tsv":
        writer = csv.writer(sys.stdout, delimiter='\t')
        writer.writerow(["contig", "locus_tag", "gene", "start", "end", "strand",
                         "category_key", "marker_id", "marker_name", "weight", "reason", "product"])
        for h in hits:
            writer.writerow([
                h['contig'], h['locus_tag'], h['gene'], h['start'], h['end'], h['strand'],
                h['category_key'], h['marker_id'], h['marker_name'], h['weight'], h['reason'], h['product']
            ])
    else:
        print_text_report(args.input_file, features, hits, clusters, pathways, karyograms)


if __name__ == "__main__":
    main()
