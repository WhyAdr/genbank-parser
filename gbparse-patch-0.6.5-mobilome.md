# GenBank Parser Patch 0.6.5 — Mobilome Audit & Knowledge-Base Expansion

Status: architecture/schema/inference audit plus source-science knowledge-base
expansion for the mobilome subsystem. This patch (a) fixes a citation
mis-anchor introduced by the 0.6.2 provenance commit, (b) expands the
provenance knowledge base from 16 to 28 verified sources, (c) anchors the
affected markers and inference rules, (d) documents the modern Type VII
classification of HEPN/MNT-type systems without making type-number claims, and
(e) resolves the tool-version drift by moving the package to 0.6.5.

Audit date: 2026-09-06

Planning baseline: main at 16118c0 ("fix(mobilome): correct and anchor source
science provenance at 1.1.0"). Baseline gate observed before any change:
152 passed / 6 optional-visualization skips, 1 known Biopython malformed-LOCUS
warning.

Fence semantics (inherited from the integration plan): only blocks fenced as
`diff` below are applyable candidate patches. Every one of them was generated
mechanically from the baseline with `git diff` and passed `git apply --check`
against a pristine clone of 16118c0; the complete nine-diff series was then
applied in a pristine clone with regenerated goldens and passed the full
repository gate (154 passed / 6 optional skips) before being embedded here.
Blocks fenced as `text` are adjudication sketches and audit evidence, not
patches.

Verification method for this audit: every new DOI below was resolved through
Crossref (`api.crossref.org/works/<DOI>`) and cross-checked against PubMed
E-utilities records (esummary/efetch); open-access status was resolved through
Unpaywall (`api.unpaywall.org/v2/<DOI>?email=...`) and PMC article-id lookups.
The executing agent must re-run the same checks (commands in section 7) before
committing; do not trust titles from memory or from secondary snippets.

---

## 1. Architecture and schema audit

### 1.1 Component map (verified, no structural change proposed)

The mobilome subsystem is a five-module pipeline over the canonical typed
parser, which keeps it decoupled from the legacy `discover` ruleset scanner:

```text
gbparse mobilome (cli.py)
  -> analyze_mobilome (mobilome/__init__.py)      parse once, hash-in/hash-out,
                                                    full-record inventory
  -> inventory_replicons (mobilome/replicons.py)  declared source /plasmid and
                                                    /chromosome classification
  -> scan_mobilome_features (mobilome/scanner.py) field-aware marker matchers,
                                                    per-reason eligibility
  -> infer_replicon_hypotheses + infer_cross_record_hypotheses
                        (mobilome/inference.py)   component/pair-rule engine
  -> report.py                                    JSON Schema 2020-12, TSV (58
                                                    columns), text; atomic write
  -> database.py (fail-closed loader)             markers.yaml + provenance.yaml
                                                    + inference.yaml as one set
```

Verified properties worth preserving in any future change: single-parse input
handling with before/after SHA-256 comparison; frozen dataclasses end to end;
deterministic ordering at every level (features by `(record_index,
feature_index)`, hits by `(record, feature, marker)`, catalog collections
sorted at load); custom databases must contain exactly the three YAML
resources (no silent inheritance of packaged provenance); atomic output with
alias/hardlink/interruption/mutation protection; and a Python-level semantic
validation pass (`validate_mobilome_report_semantics`) that catches duplicate
hit/hypothesis ids, hit-record mismatches, and participant-role contract
violations that generic JSON Schema cannot express.

### 1.2 Schema audit (report.schema.json, gbparse.mobilome.v1)

The schema is internally consistent: draft 2020-12, `additionalProperties:
false` on every object, complete `required` lists, bounded strength (1-3),
closed enums for replicon class, status, match mode, and handoff status. Two
long-standing design notes, both already documented in
`docs/mobilome_evidence_reference.md`:

- Handoffs carry no `source_ids` in v1 (the fail-closed handoff parser rejects
  unknown keys), so handoff citations must live in docs plus provenance. This
  patch uses that v1-compatible home; a v2 schema could add handoff
  `source_ids`.
- Record-level spatial-inference limitations ride in the classification
  limitation list; a v2 `spatial_limitations` inventory field was already
  proposed by patch 0.3.

### 1.3 Findings register (architecture)

| Id | Severity | Finding | Disposition |
|---|---|---|---|
| A1 | HIGH | **Citation mis-anchor in HEAD.** Commit 16118c0 anchored `feldgarden-2021-amrfinderplus` (AMRFinderPlus) to `phage_packaging_candidate` and `chen-2005-vfdb` (VFDB) to `phage_structure_candidate`, while `generic_amr_candidate` and `generic_vf_candidate` still cite only `local-annotation-candidate-policy`. `mobilome-patch-0.3.md` §4.2 (lines 360-377) documents the opposite intent: the feldgarden hunk carries the context "Retains a generic antimicrobial-resistance annotation candidate." and the chen hunk carries "Retains a generic virulence-associated annotation candidate." The commit instead applied the same `+` lines at hunks `@@ -470` (packaging) and `@@ -488` (structural). Evidence: `git show 16118c0 -- src/genbank_parser/data/mobilome/markers.yaml` vs the patch-0.3 md. Consequence: terminase/portal report reasons cite an AMR paper and capsid/tail reasons cite a virulence database; the intended AMR/VF anchors are absent. Nothing failed because the 0.6.2 citation-integrity test does not assert marker-level anchors. The patch-0.3 execution record ("applied as written") is therefore inaccurate for markers.yaml. | **Fixed in section 4** plus regression test `test_amr_and_vf_citations_are_anchored_to_their_own_markers` |
| A2 | MEDIUM | **Tool-version drift.** `pyproject.toml` and all three CI asserts say `0.6.0`, while `changelogs.md` documents 0.6.1 and 0.6.2 releases. | **Fixed in section 4**: everything moves to 0.6.5 (numbering requested by the maintainer), changelog entry added |
| A3 | LOW | Handoff citations have no schema home in v1 (see 1.2). | Documented; docs+provenance used by this patch; v2 candidate |
| A4 | LOW | `database.py` `_parse_disabled_rules` hardcodes the marker id `pxo_numbered_product_annotation` inside the generic validator's allowed retained-evidence vocabulary (`set(retained) - facet_ids - {"pxo_numbered_product_annotation"}`). Renaming that marker would surface as a confusing "unknown retained evidence" database error rather than a rename-aware message. | Deferred; validator change, no behavioral defect |
| A5 | LOW | `missing_components` mixes vocabularies: component-rule failures emit facet ids (`mobility_mpf`), the replication path emits a component id (`replication_candidate`). It is asserted behavior (`test_mobility_component_subsets_report_exact_missing_facets`) but undocumented in the schema. | Documented here; typed v2 candidate |

---

## 2. Inference-logic audit

### 2.1 Verified-correct behaviors (probed with the packaged tests plus targeted reads)

- Eligibility gating is per reason: `strength >= min_evidence` AND the
  marker's `requires_non_pseudo` / `requires_complete` policy; below-threshold
  and pseudo/partial matches are retained in reports as ineligible raw
  evidence rather than dropped.
- `_functional_hits` applies a uniform non-pseudo/non-complete filter before
  any component or pair rule, so observation-only markers can never support a
  hypothesis through a pseudogenized feature.
- Distinct-feature accounting is sound: `_choose_hits` threads forbidden
  feature sets through `_component_support` and `_rule_component_support`, so
  one physical feature cannot satisfy two facets or two roles when
  `distinct_features` is set; selection order is deterministic.
- Spatial semantics are correct: `intervening_gap_bp` computes the minimum
  pairwise segment gap (adjacent/overlapping is 0, compound joins do not
  inflate gaps); `circular_feature_distance_bp` checks both directions with
  modular arithmetic and validates `record_length`; unlocatable features skip
  pairing and emit the record-level limitation; the 1999/2000/2001 bp boundary
  tests pin inclusive behavior on both linear and circular topologies.
- Cross-record evaluation excludes same-record helpers, runs over all scanned
  records before the `include` presentation filter, and validates target and
  helper distinctness again at report-serialization time.
- The fail-closed database loader enforces unique YAML keys, per-field
  strength ceilings, cross-references (markers to facets/sources, rules to
  components/markers, local sources to rule ids), and major-version agreement
  across the three resources.

### 2.2 Findings register (inference; all deferred deliberately)

| Id | Severity | Finding | Disposition |
|---|---|---|---|
| C1 | MEDIUM | **The "configured" inference policy is co-configured in code.** `infer_replicon_hypotheses` iterates a hardcoded tuple of exactly three rule ids; `infer_cross_record_hypotheses` fetches `possible_helper_dependent_mobilization` by id; `_infer_component_rule` special-cases `replication_annotation_candidate` (emits an insufficiency hypothesis for every record even with zero evidence, hardcodes the `"replication_candidate"` missing-component id, and the replication-specific summary and limitation wording at inference.py lines ~208-232). A new component rule added only to inference.yaml would be silently ignored. | Deferred to a v2 data-driven dispatch; behavior-affecting refactor beyond this patch's scope |
| C2 | LOW | **TA gap wording is circular-specific but semantics are topology-generic.** `max_circular_gap_bp` also gates linear records (`_distance_between_hits` uses `intervening_gap_bp` for linear topology) and the emitted limitation reads "Configured maximum circular gap: X bp; observed gap: Y bp." on both topologies. Latent: goldens only exercise a circular fixture. | Deferred; wording plus contract change, needs a linear fixture first |
| C3 | LOW | **TA pairing is strand-agnostic.** Same-record opposite-strand toxN/toxI within 2000 bp pair up; the MEOR module requires a known common strand for clusters, the mobilome pair rules do not. Type III TA loci are operonic, so a `same_strand` rule flag is scientifically defensible. | Deferred; requires fixture plus adversarial negative first |
| C4 | LOW | **TA pair emission is O(N x M).** One hypothesis per eligible pair with no nearest-neighbor pruning; TA-annotation-rich records produce combinatorial (deterministic) hypothesis lists. | Deferred; report-bloat concern only |
| C5 | LOW | Both `_infer_toxin_antitoxin_pairs` and `record_inference_limitations` gate on `len(required_marker_ids) != 2`; a future three-marker TA rule would silently skip pair inference and limitation bookkeeping. | Documented for the next TA expansion |
| C6 | LOW | `_annotation_provenance` recognizes only Bakta (name, version, database version from structured comments); Prokka/PGAP/RefSeq inputs report a `null` pipeline, which the docs already describe as intentional. | Documented; a pipeline-detection table is a future addition |
| C7 | LOW | **Marker coverage gaps (scientific, not defects).** The MPF gene regex `^tra(?:A|B|C|E|F|K|L|N|U|W)$` omits `virB`-family MPF names and `trb` genes; the cas regex does not match subtype-suffixed cas genes (`cas8f`, `cas12a`, `cas13d`); the replication-initiation gene list lacks `repX`/`repL` used by pXO1/pXO2-like plasmids. | Deferred with sketches in section 5, per the repo's fixtures-and-negatives-first discipline |

---

## 3. Knowledge-base expansion

### 3.1 New sources (12), all verified this session

Every row was resolved through Crossref (title, year, journal, authors) and
PubMed (print year, PMID, PMCID where deposited), and OA status through
Unpaywall. Titles below are copied verbatim from Crossref/PubMed rendering
(ASCII hyphens; italic markup stripped), per the repo's title-verbatim rule.

| New source id | Title (verbatim) | Venue / year | DOI | PMID | PMCID | Access |
|---|---|---|---|---|---|---|
| schwengers-2021-bakta | Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification | Microb Genom 2021 | 10.1099/mgen.0.000685 | 34739369 | PMC8743544 | OA |
| garcillan-barcia-2025-extended-mobility | The extended mobility of plasmids | Nucleic Acids Res 2025 | 10.1093/nar/gkaf652 | 40694848 | PMC12282955 | OA |
| johnson-2015-ice | Integrative and Conjugative Elements (ICEs): What They Do and How They Work | Annu Rev Genet 2015 | 10.1146/annurev-genet-112414-055018 | 26473380 | PMC5180612 | OA at publisher; also free at PMC |
| mazel-2006-integrons | Integrons: agents of bacterial evolution | Nat Rev Microbiol 2006 | 10.1038/nrmicro1462 | 16845431 | none | **PAYWALLED** |
| qiu-2022-ta-classification | Toxin-antitoxin systems: Classification, biological roles, and applications | Microbiol Res 2022 | 10.1016/j.micres.2022.127159 | 35969944 | none | OA (publisher; not in PMC) |
| yao-2026-mnt-hepn | Dual regulation of HEPN RNase in fused MNT-HEPN toxin-antitoxin systems via protein OligoAMPylation and oligomerization | Nucleic Acids Res 2026 | 10.1093/nar/gkag228 | 41854079 | PMC13000451 | OA |
| liu-2022-vfdb | VFDB 2022: a general classification scheme for bacterial virulence factors | Nucleic Acids Res 2022 | 10.1093/nar/gkab1107 | 34850947 | PMC8728188 | OA |
| partridge-2018-mge-amr | Mobile Genetic Elements Associated with Antimicrobial Resistance | Clin Microbiol Rev 2018 | 10.1128/CMR.00088-17 | 30068738 | PMC6148190 | Paywalled at publisher; free at PMC |
| arndt-2016-phaster | PHASTER: a better, faster version of the PHAST phage search tool | Nucleic Acids Res 2016 | 10.1093/nar/gkw387 | 27141966 | PMC4987931 | OA |
| camargo-2024-genomad | Identification of mobile genetic elements with geNomad | Nat Biotechnol 2024 | 10.1038/s41587-023-01953-y | 37735266 | PMC11324519 | OA |
| xie-2017-isescan | ISEScan: automated identification of insertion sequence elements in prokaryotic genomes | Bioinformatics 2017 | 10.1093/bioinformatics/btx433 | 29077810 | none | OA |
| johansson-2021-mobileelementfinder | Detection of mobile genetic elements associated with antibiotic resistance in Salmonella enterica using a newly developed web tool: MobileElementFinder | J Antimicrob Chemother 2021 | 10.1093/jac/dkaa390 | 33009809 | PMC7729385 | OA |

Anchor plan executed in section 4: Bakta primary citation joins the three
markers that already cite `bakta-origin-docs`; the extended-mobility review
joins the relaxase/T4CP/MPF markers and the three mobility rules (Smillie
2010 is retained everywhere it already appears); the TA classification review
joins the four TA markers and both pair rules; the fused MNT-HEPN study joins
the HEPN/MNT rule and its limitations; ICE and integron reviews join the
recombination marker; the MGE-AMR review joins the AMR candidate marker; the
VFDB 2022 update joins the VF candidate marker; and the four external-caller
papers (PHASTER, geNomad, ISEScan, MobileElementFinder) live in provenance
plus the docs handoff paragraph, because v1 handoffs cannot carry source ids.

### 3.2 Corrections to the audit trail (from mobilome-patch-0.3.md §5)

- **Bakta citation**: the DOI guessed in patch 0.3 §5
  (`10.1016/j.jbiotec.2021.07.018`) resolves to HTTP 404 through Crossref and
  Unpaywall; the paper is actually in Microbial Genomics
  (`10.1099/mgen.0.000685`, PMID 34739369), and the real title tail is
  "via alignment-free sequence identification" — not "of bacterial genomes,
  plasmids and phage sequences". Patch 0.3's follow-up sketch is therefore
  superseded by the verified entry in section 4.1.
- **Qiu 2022 (PMID 35969944)**: verified as "Toxin-antitoxin systems:
  Classification, biological roles, and applications", Microbiol Res 2022,
  DOI 10.1016/j.micres.2022.127159, as anticipated by patch 0.3 §5.
- **gkag228 (2026)**: verified as "Dual regulation of HEPN RNase in fused
  MNT-HEPN toxin-antitoxin systems via protein OligoAMPylation and
  oligomerization", NAR 2026, PMID 41854079, PMC13000451, exactly the record
  patch 0.3 §5 said Crossref had confirmed.
- The Type VII documentation line in the HEPN/MNT rule follows the plan
  patch 0.3 §5 sketched: current reviews document HEPN/MNT-type systems among
  Type VII, while the rule itself deliberately assigns no type number to any
  annotation match. The wording passes the repository forbidden-phrase scan.

### 3.3 Paywall status of the full 28-source knowledge base (for the maintainer)

Resolving authority: Unpaywall `is_oa` plus PMC article-id lookup, 2026-09-06.

**No free full text found (3) — these are the ones to fetch manually:**

| Source | DOI | Note |
|---|---|---|
| mazel-2006-integrons (NEW) | 10.1038/nrmicro1462 | Nature Reviews Microbiology 2006; no OA, no PMC |
| garcillan-barcia-2009-relaxases (existing) | 10.1111/j.1574-6976.2009.00168.x | FEMS Microbiology Reviews 2009; no OA, no PMC |
| makarova-2020-crispr (existing) | 10.1038/s41579-019-0299-x | No publisher OA; the accepted manuscript is freely available at the Wageningen repository (research.wur.nl) if that suffices |

**Paywalled at the publisher, but free full text at PMC (6):**

| Source | Free copy |
|---|---|
| smillie-2010-mobility | PMC2937521 |
| akhtar-2012-pxo1 | PMC4186245 |
| erez-2017-arbitrium | PMC5378303 (author manuscript) |
| bouet-2019-partition | PMC11573283 (author manuscript) |
| partridge-2018-mge-amr (NEW) | PMC6148190 |
| johnson-2015-ice (NEW) | PMC5180612 |

Everything else — including all twelve new sources except mazel-2006 and
publisher-level partridge-2018 — is open access at the publisher (and the two
official-documentation sources plus three local rules are not papers at all).

## 4. Applyable patch series

Nine diffs, one atomic series against baseline 16118c0. Apply in any order
(each touches a distinct file); verify per section 7. Resource versions move
1.1.0/1.1.0/1.1.0 -> 1.2.0 (catalog, provenance, inference all change; major
versions agree) and the package moves 0.6.0 -> 0.6.5. Goldens must be
regenerated after applying (command in section 7); they are deliberately not
embedded as diffs, following the patch-0.3 precedent.

### 4.1 provenance.yaml — twelve new verified sources, version 1.2.0

```diff
diff --git a/src/genbank_parser/data/mobilome/provenance.yaml b/src/genbank_parser/data/mobilome/provenance.yaml
index ce95dd7..e6ad881 100644
--- a/src/genbank_parser/data/mobilome/provenance.yaml
+++ b/src/genbank_parser/data/mobilome/provenance.yaml
@@ -1,5 +1,5 @@
 schema_version: 1
-provenance_version: "1.1.0"
+provenance_version: "1.2.0"
 
 sources:
   - id: blower-2012-toxin
@@ -162,6 +162,147 @@ sources:
       - Insertion sequences carry standardized family nomenclature in a dedicated reference center.
     limitations:
       - An IS-like annotation is not a family assignment without sequence comparison.
+  - id: schwengers-2021-bakta
+    type: primary_article
+    title: "Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification"
+    year: 2021
+    doi: 10.1099/mgen.0.000685
+    pmid: "34739369"
+    pmcid: PMC8743544
+    url: https://pmc.ncbi.nlm.nih.gov/articles/PMC8743544/
+    supports:
+      - Bakta is a versioned command-line annotation pipeline whose database releases are declared alongside its outputs.
+    limitations:
+      - Tool output forms change across versions; observed input annotations remain authoritative over general tool claims.
+  - id: garcillan-barcia-2025-extended-mobility
+    type: primary_article
+    title: "The extended mobility of plasmids"
+    year: 2025
+    doi: 10.1093/nar/gkaf652
+    pmid: "40694848"
+    pmcid: PMC12282955
+    url: https://pmc.ncbi.nlm.nih.gov/articles/PMC12282955/
+    supports:
+      - Plasmid mobility extends beyond canonical conjugation, and mobility frameworks cover mobilization routes broader than a single transfer pathway.
+    limitations:
+      - The review's extended-mobility framing does not transfer mechanism or route conclusions to annotation co-occurrence in a new genome.
+  - id: johnson-2015-ice
+    type: primary_article
+    title: "Integrative and Conjugative Elements (ICEs): What They Do and How They Work"
+    year: 2015
+    doi: 10.1146/annurev-genet-112414-055018
+    pmid: "26473380"
+    pmcid: PMC5180612
+    url: https://pmc.ncbi.nlm.nih.gov/articles/PMC5180612/
+    supports:
+      - Integrative and conjugative elements are modular mobile elements that combine chromosomal integration and conjugative transfer functions in one element.
+    limitations:
+      - Integration and conjugation annotations in one record do not establish an ICE boundary or an active element.
+  - id: mazel-2006-integrons
+    type: primary_article
+    title: "Integrons: agents of bacterial evolution"
+    year: 2006
+    doi: 10.1038/nrmicro1462
+    pmid: "16845431"
+    url: https://doi.org/10.1038/nrmicro1462
+    supports:
+      - Integrons are recombination platforms that capture exogenous gene cassettes through an integron integrase and convert them into functional genes.
+    limitations:
+      - An integron-like annotation is not a cassette-array inventory or a resistance-gene prediction without sequence analysis.
+  - id: qiu-2022-ta-classification
+    type: primary_article
+    title: "Toxin-antitoxin systems: Classification, biological roles, and applications"
+    year: 2022
+    doi: 10.1016/j.micres.2022.127159
+    pmid: "35969944"
+    url: https://doi.org/10.1016/j.micres.2022.127159
+    supports:
+      - Toxin-antitoxin systems are classified into types according to the nature and mechanism of the antitoxin.
+    limitations:
+      - A review classification does not transfer to an annotation match in a new genome, and this catalog assigns no type number to annotations.
+  - id: yao-2026-mnt-hepn
+    type: primary_article
+    title: "Dual regulation of HEPN RNase in fused MNT-HEPN toxin-antitoxin systems via protein OligoAMPylation and oligomerization"
+    year: 2026
+    doi: 10.1093/nar/gkag228
+    pmid: "41854079"
+    pmcid: PMC13000451
+    url: https://pmc.ncbi.nlm.nih.gov/articles/PMC13000451/
+    supports:
+      - Fused MNT-HEPN proteins are documented members of the HEPN/MNT toxin-antitoxin family with AMPylation-based regulation in studied systems.
+    limitations:
+      - Regulation biochemistry in studied systems is not generalized to annotation matches in new genomes.
+  - id: liu-2022-vfdb
+    type: primary_article
+    title: "VFDB 2022: a general classification scheme for bacterial virulence factors"
+    year: 2022
+    doi: 10.1093/nar/gkab1107
+    pmid: "34850947"
+    pmcid: PMC8728188
+    url: https://pmc.ncbi.nlm.nih.gov/articles/PMC8728188/
+    supports:
+      - Virulence-factor determination relies on a curated, versioned database with an explicit classification scheme.
+    limitations:
+      - Database releases supersede earlier descriptions; version provenance is required for any future comparison.
+  - id: partridge-2018-mge-amr
+    type: primary_article
+    title: "Mobile Genetic Elements Associated with Antimicrobial Resistance"
+    year: 2018
+    doi: 10.1128/CMR.00088-17
+    pmid: "30068738"
+    pmcid: PMC6148190
+    url: https://pmc.ncbi.nlm.nih.gov/articles/PMC6148190/
+    supports:
+      - Resistance genes are associated with diverse mobile genetic elements, including plasmids, transposons, integrons, and gene cassettes.
+    limitations:
+      - Association documented in the review does not establish that a generic resistance-like annotation in a record is mobile.
+  - id: arndt-2016-phaster
+    type: primary_article
+    title: "PHASTER: a better, faster version of the PHAST phage search tool"
+    year: 2016
+    doi: 10.1093/nar/gkw387
+    pmid: "27141966"
+    pmcid: PMC4987931
+    url: https://pmc.ncbi.nlm.nih.gov/articles/PMC4987931/
+    supports:
+      - Dedicated prophage identification and annotation is performed from sequence by an external tool, not from product annotations.
+    limitations:
+      - Prophage boundary predictions carry method-specific uncertainty and require external provenance.
+  - id: camargo-2024-genomad
+    type: primary_article
+    title: "Identification of mobile genetic elements with geNomad"
+    year: 2024
+    doi: 10.1038/s41587-023-01953-y
+    pmid: "37735266"
+    pmcid: PMC11324519
+    url: https://pmc.ncbi.nlm.nih.gov/articles/PMC11324519/
+    supports:
+      - Virus and plasmid identification from sequence data is performed by external tools combining marker genes and machine-learning classifiers.
+    limitations:
+      - External classifications depend on model versions, marker sets, and thresholds, which must be reported with any imported result.
+  - id: xie-2017-isescan
+    type: primary_article
+    title: "ISEScan: automated identification of insertion sequence elements in prokaryotic genomes"
+    year: 2017
+    doi: 10.1093/bioinformatics/btx433
+    pmid: "29077810"
+    url: https://doi.org/10.1093/bioinformatics/btx433
+    supports:
+      - Insertion-sequence elements can be identified automatically from genome sequence by an external tool.
+    limitations:
+      - Sequence-level IS family calls cannot be reproduced from product annotations alone.
+  - id: johansson-2021-mobileelementfinder
+    type: primary_article
+    title: "Detection of mobile genetic elements associated with antibiotic resistance in Salmonella enterica using a newly developed web tool: MobileElementFinder"
+    year: 2021
+    doi: 10.1093/jac/dkaa390
+    pmid: "33009809"
+    pmcid: PMC7729385
+    url: https://pmc.ncbi.nlm.nih.gov/articles/PMC7729385/
+    supports:
+      - Mobile-element detection from sequence is performed by external tools backed by curated element databases.
+    limitations:
+      - The tool's database composition and thresholds bound any future imported result.
   - id: bakta-origin-docs
     type: official_documentation
     title: Bakta documentation
```

### 4.2 markers.yaml — re-anchor the 0.6.2 mis-anchors, anchor new sources, bump catalog

Note how this diff is the inverse image of the 0.6.2 mis-application: the
hunks at `@@ -470` and `@@ -488` remove `feldgarden-2021-amrfinderplus` and
`chen-2005-vfdb` from the phage markers (exactly where commit 16118c0 put
them), while the hunks at `@@ -585` and `@@ -603` add them to the AMR and VF
candidate markers (exactly where mobilome-patch-0.3.md documented them).

```diff
diff --git a/src/genbank_parser/data/mobilome/markers.yaml b/src/genbank_parser/data/mobilome/markers.yaml
index 399142a..9fbac36 100644
--- a/src/genbank_parser/data/mobilome/markers.yaml
+++ b/src/genbank_parser/data/mobilome/markers.yaml
@@ -1,5 +1,5 @@
 schema_version: 1
-catalog_version: "1.1.0"
+catalog_version: "1.2.0"
 
 facets:
   - id: replication_origin
@@ -65,7 +65,7 @@ markers:
         strength: 3
     requires_non_pseudo: false
     requires_complete: false
-    sources: [bakta-origin-docs, insdc-feature-table]
+    sources: [bakta-origin-docs, insdc-feature-table, schwengers-2021-bakta]
     interpretation: Supports the presence of an explicitly annotated oriT feature.
     limitations:
       - Annotation does not prove a cognate relaxase relationship.
@@ -105,7 +105,7 @@ markers:
         strength: 3
     requires_non_pseudo: false
     requires_complete: false
-    sources: [bakta-origin-docs, insdc-feature-table]
+    sources: [bakta-origin-docs, insdc-feature-table, schwengers-2021-bakta]
     interpretation: Supports an explicitly annotated plasmid replication-origin feature.
     limitations:
       - Annotation alone does not establish replication mechanism or autonomy.
@@ -152,7 +152,7 @@ markers:
         strength: 2
     requires_non_pseudo: true
     requires_complete: true
-    sources: [garcillan-barcia-2009-relaxases, smillie-2010-mobility]
+    sources: [garcillan-barcia-2009-relaxases, garcillan-barcia-2025-extended-mobility, smillie-2010-mobility]
     interpretation: Retains a relaxase annotation candidate.
     limitations:
       - Annotation does not demonstrate nicking or transfer.
@@ -175,7 +175,7 @@ markers:
         strength: 2
     requires_non_pseudo: true
     requires_complete: true
-    sources: [christie-2025-t4ss, smillie-2010-mobility]
+    sources: [christie-2025-t4ss, garcillan-barcia-2025-extended-mobility, smillie-2010-mobility]
     interpretation: Retains a Type IV coupling-protein annotation candidate.
     limitations:
       - Annotation does not establish a complete transfer apparatus.
@@ -199,7 +199,7 @@ markers:
         strength: 2
     requires_non_pseudo: true
     requires_complete: true
-    sources: [christie-2025-t4ss, smillie-2010-mobility]
+    sources: [christie-2025-t4ss, garcillan-barcia-2025-extended-mobility, smillie-2010-mobility]
     interpretation: Retains a mating-pair-formation annotation component.
     limitations:
       - Multiple annotated components do not establish expression or transfer.
@@ -222,7 +222,7 @@ markers:
         strength: 2
     requires_non_pseudo: true
     requires_complete: true
-    sources: [blower-2012-toxin]
+    sources: [blower-2012-toxin, qiu-2022-ta-classification]
     interpretation: Supports a ToxN-family toxin annotation.
     limitations:
       - Annotation does not establish toxin activity.
@@ -251,7 +251,7 @@ markers:
         strength: 1
     requires_non_pseudo: true
     requires_complete: true
-    sources: [blower-2012-toxin]
+    sources: [blower-2012-toxin, qiu-2022-ta-classification]
     interpretation: Supports a ToxI RNA antitoxin annotation.
     limitations:
       - ToxI correspondence does not establish Type III system activity.
@@ -274,7 +274,7 @@ markers:
         strength: 2
     requires_non_pseudo: true
     requires_complete: true
-    sources: [yao-2020-hepn-mnt]
+    sources: [qiu-2022-ta-classification, yao-2020-hepn-mnt]
     interpretation: Supports a HepT-family toxin annotation.
     limitations:
       - Annotation does not establish toxin activity.
@@ -297,7 +297,7 @@ markers:
         strength: 2
     requires_non_pseudo: true
     requires_complete: true
-    sources: [yao-2020-hepn-mnt]
+    sources: [qiu-2022-ta-classification, yao-2020-hepn-mnt]
     interpretation: Supports an MntA-family antitoxin annotation.
     limitations:
       - The cited chemistry is not generalized to every annotation match.
@@ -344,7 +344,7 @@ markers:
         strength: 1
     requires_non_pseudo: false
     requires_complete: false
-    sources: [local-annotation-candidate-policy, siguier-2006-isfinder]
+    sources: [johnson-2015-ice, local-annotation-candidate-policy, mazel-2006-integrons, siguier-2006-isfinder]
     interpretation: Retains recombination or integration annotation evidence.
     limitations:
       - Observation does not establish an element boundary or activity.
@@ -410,7 +410,7 @@ markers:
         strength: 1
     requires_non_pseudo: false
     requires_complete: false
-    sources: [bakta-origin-docs, insdc-feature-table]
+    sources: [bakta-origin-docs, insdc-feature-table, schwengers-2021-bakta]
     interpretation: Supports an explicitly annotated CRISPR repeat-array feature.
     limitations:
       - Annotation does not establish active immunity or target specificity.
@@ -470,7 +470,7 @@ markers:
         strength: 1
     requires_non_pseudo: false
     requires_complete: false
-    sources: [feldgarden-2021-amrfinderplus, local-annotation-candidate-policy]
+    sources: [local-annotation-candidate-policy]
     interpretation: Retains phage packaging annotation evidence.
     limitations:
       - Observation does not establish a complete phage module.
@@ -488,7 +488,7 @@ markers:
         strength: 1
     requires_non_pseudo: false
     requires_complete: false
-    sources: [chen-2005-vfdb, local-annotation-candidate-policy]
+    sources: [local-annotation-candidate-policy]
     interpretation: Retains phage structural annotation evidence.
     limitations:
       - Observation does not establish a complete phage module.
@@ -585,7 +585,7 @@ markers:
         strength: 1
     requires_non_pseudo: false
     requires_complete: false
-    sources: [local-annotation-candidate-policy]
+    sources: [feldgarden-2021-amrfinderplus, local-annotation-candidate-policy, partridge-2018-mge-amr]
     interpretation: Retains a generic antimicrobial-resistance annotation candidate.
     limitations:
       - Product text alone does not establish resistance phenotype or determinant status.
@@ -603,7 +603,7 @@ markers:
         strength: 1
     requires_non_pseudo: false
     requires_complete: false
-    sources: [local-annotation-candidate-policy]
+    sources: [chen-2005-vfdb, liu-2022-vfdb, local-annotation-candidate-policy]
     interpretation: Retains a generic virulence-associated annotation candidate.
     limitations:
       - Product text alone does not establish virulence phenotype.
```

### 4.3 inference.yaml — anchor rules, document Type VII, bump inference version

```diff
diff --git a/src/genbank_parser/data/mobilome/inference.yaml b/src/genbank_parser/data/mobilome/inference.yaml
index fc72ca6..afde45c 100644
--- a/src/genbank_parser/data/mobilome/inference.yaml
+++ b/src/genbank_parser/data/mobilome/inference.yaml
@@ -1,5 +1,5 @@
 schema_version: 1
-inference_version: "1.1.0"
+inference_version: "1.2.0"
 
 components:
   replication_candidate:
@@ -38,7 +38,7 @@ rules:
     limitations:
       - A cognate relaxase-oriT relationship is untested.
       - Expression and transfer are untested.
-    sources: [smillie-2010-mobility]
+    sources: [garcillan-barcia-2025-extended-mobility, smillie-2010-mobility]
   - id: helper_machinery_annotation_candidate
     kind: mobility
     required_components:
@@ -48,7 +48,7 @@ rules:
     wording: Conjugation-machinery annotation candidate
     limitations:
       - Multiple components do not establish a complete or expressed transfer system.
-    sources: [smillie-2010-mobility]
+    sources: [garcillan-barcia-2025-extended-mobility, smillie-2010-mobility]
   - id: possible_helper_dependent_mobilization
     kind: mobility
     required_components:
@@ -63,7 +63,7 @@ rules:
       - Record taxonomic and same-cell co-affiliation are untested.
       - Expression and transfer are untested.
       - Co-transfer of both records is not implied.
-    sources: [smillie-2010-mobility]
+    sources: [garcillan-barcia-2025-extended-mobility, smillie-2010-mobility]
   - id: type_iii_toxin_antitoxin_pair_candidate
     kind: toxin_antitoxin
     required_marker_ids: [toxn, toxi]
@@ -75,7 +75,7 @@ rules:
     limitations:
       - Proximity is a versioned engineering heuristic, not proof of an operon.
       - Expression and toxin-antitoxin activity are untested.
-    sources: [blower-2012-toxin]
+    sources: [blower-2012-toxin, qiu-2022-ta-classification]
   - id: hepn_mnt_annotation_pair_candidate
     kind: toxin_antitoxin
     required_marker_ids: [hept, mnta]
@@ -87,7 +87,8 @@ rules:
     limitations:
       - Proximity is a versioned engineering heuristic, not proof of an operon.
       - The cited biochemical correspondence is not generalized beyond annotation.
-    sources: [yao-2020-hepn-mnt]
+      - Current classification reviews document HEPN/MNT-type systems among Type VII toxin-antitoxin systems; this rule deliberately assigns no type number to an annotation match.
+    sources: [qiu-2022-ta-classification, yao-2020-hepn-mnt, yao-2026-mnt-hepn]
 
 disabled_aggregate_rules:
   - id: pxo_like_annotation_pattern_candidate
```

### 4.4 tests/test_mobilome_database.py — version asserts, anchor ledger, mis-anchor regression

```diff
diff --git a/tests/test_mobilome_database.py b/tests/test_mobilome_database.py
index 4f16409..5268887 100644
--- a/tests/test_mobilome_database.py
+++ b/tests/test_mobilome_database.py
@@ -26,9 +26,9 @@ def test_packaged_database_is_versioned_hashed_and_schema_valid() -> None:
     )
     schema = json.loads(schema_path.read_text(encoding="utf-8"))
 
-    assert database.catalog_version == "1.1.0"
-    assert database.inference_version == "1.1.0"
-    assert database.provenance_version == "1.1.0"
+    assert database.catalog_version == "1.2.0"
+    assert database.inference_version == "1.2.0"
+    assert database.provenance_version == "1.2.0"
     assert database.database_source == "packaged"
     assert database.source_paths == ()
     assert [resource.name for resource in database.resources] == [
@@ -69,8 +69,76 @@ def test_provenance_sources_use_verified_citations() -> None:
         "makarova-2020-crispr",
         "siguier-2006-isfinder",
         "bouet-2019-partition",
+        "schwengers-2021-bakta",
+        "garcillan-barcia-2025-extended-mobility",
+        "johnson-2015-ice",
+        "mazel-2006-integrons",
+        "qiu-2022-ta-classification",
+        "yao-2026-mnt-hepn",
+        "liu-2022-vfdb",
+        "partridge-2018-mge-amr",
+        "arndt-2016-phaster",
+        "camargo-2024-genomad",
+        "xie-2017-isescan",
+        "johansson-2021-mobileelementfinder",
     ):
         assert anchor in by_id
+    assert by_id["schwengers-2021-bakta"].payload["pmid"] == "34739369"
+    assert (
+        by_id["garcillan-barcia-2025-extended-mobility"].payload["pmid"] == "40694848"
+    )
+    assert by_id["mazel-2006-integrons"].payload["doi"] == "10.1038/nrmicro1462"
+    assert (
+        by_id["johansson-2021-mobileelementfinder"]
+        .payload["title"]
+        .endswith("web tool: MobileElementFinder")
+    )
+
+
+def test_amr_and_vf_citations_are_anchored_to_their_own_markers() -> None:
+    database = load_mobilome_database()
+    markers = database.marker_map
+
+    assert markers["generic_amr_candidate"].sources == (
+        "feldgarden-2021-amrfinderplus",
+        "local-annotation-candidate-policy",
+        "partridge-2018-mge-amr",
+    )
+    assert markers["generic_vf_candidate"].sources == (
+        "chen-2005-vfdb",
+        "liu-2022-vfdb",
+        "local-annotation-candidate-policy",
+    )
+    # The AMRFinderPlus and VFDB citations must never support phage markers.
+    assert markers["phage_packaging_candidate"].sources == (
+        "local-annotation-candidate-policy",
+    )
+    assert markers["phage_structure_candidate"].sources == (
+        "local-annotation-candidate-policy",
+    )
+
+
+def test_ta_and_mobility_rules_carry_their_configured_sources() -> None:
+    database = load_mobilome_database()
+    rules = database.inference_rule_map
+
+    assert rules["hepn_mnt_annotation_pair_candidate"].sources == (
+        "qiu-2022-ta-classification",
+        "yao-2020-hepn-mnt",
+        "yao-2026-mnt-hepn",
+    )
+    assert rules["type_iii_toxin_antitoxin_pair_candidate"].sources == (
+        "blower-2012-toxin",
+        "qiu-2022-ta-classification",
+    )
+    assert rules["mobilizable_core_annotation_candidate"].sources == (
+        "garcillan-barcia-2025-extended-mobility",
+        "smillie-2010-mobility",
+    )
+    assert any(
+        "Type VII" in limitation
+        for limitation in rules["hepn_mnt_annotation_pair_candidate"].limitations
+    )
 
 
 def test_complete_custom_database_is_loaded_as_one_resource_set(tmp_path: Path) -> None:
```

### 4.5 tests/test_meor_report.py — tool-version assert 0.6.5

The MEOR report embeds the shared package version; it must follow the bump.

```diff
diff --git a/tests/test_meor_report.py b/tests/test_meor_report.py
index 835a98c..e8df55b 100644
--- a/tests/test_meor_report.py
+++ b/tests/test_meor_report.py
@@ -27,7 +27,7 @@ def test_report_contracts() -> None:
         "karyograms",
     ]
     assert payload["schema_version"] == "gbparse.meor.v1"
-    assert payload["tool_version"] == "0.6.0"
+    assert payload["tool_version"] == "0.6.5"
     assert payload["catalog_version"] == "1.1"
     assert payload["parameters"] == {
         "min_weight": 1,
@@ -45,7 +45,7 @@ def test_report_contracts() -> None:
         assert section in text
     assert "Low (W=1)" in text
     assert "1 kb windows" in text
-    assert "gbparse Version: 0.6.0" in text
+    assert "gbparse Version: 0.6.5" in text
     assert "MEOR Catalog Version: 1.1" in text
 
 
```

### 4.6 pyproject.toml — package version 0.6.5 (resolves the 0.6.0/0.6.2 drift)

```diff
diff --git a/pyproject.toml b/pyproject.toml
index bb055c1..9a3b3c0 100644
--- a/pyproject.toml
+++ b/pyproject.toml
@@ -4,7 +4,7 @@ build-backend = "setuptools.build_meta"
 
 [project]
 name = "genbank-parser"
-version = "0.6.0"
+version = "0.6.5"
 description = "Canonical GenBank flatfile parser, biological validation, and genomic evidence engine"
 readme = "README.md"
 requires-python = ">=3.10"
```

### 4.7 .github/workflows/ci.yml — wheel-gate asserts at 0.6.5 / 1.2.0

```diff
diff --git a/.github/workflows/ci.yml b/.github/workflows/ci.yml
index c635930..49a1ce7 100644
--- a/.github/workflows/ci.yml
+++ b/.github/workflows/ci.yml
@@ -52,13 +52,13 @@ jobs:
 
     - name: Verify installed API, data, CLI, and dependencies
       run: |
-        .wheel-venv/bin/python -c "from importlib.resources import files; from genbank_parser import __version__; from genbank_parser.meor import load_meor_database; from genbank_parser.mobilome import load_mobilome_database; db=load_mobilome_database(); root=files('genbank_parser'); assert __version__ == '0.6.0'; assert load_meor_database().catalog_version == '1.1'; assert db.catalog_version == '1.1.0'; assert db.provenance_version == '1.1.0'; assert db.inference_version == '1.1.0'; assert root.joinpath('py.typed').is_file(); assert root.joinpath('data/meor/markers.yaml').is_file(); assert files('genbank_parser.rulesets').joinpath('mobilome.yaml').is_file(); assert root.joinpath('data/mobilome/markers.yaml').is_file(); assert root.joinpath('data/mobilome/provenance.yaml').is_file(); assert root.joinpath('data/mobilome/inference.yaml').is_file(); assert root.joinpath('data/mobilome/report.schema.json').is_file()"
+        .wheel-venv/bin/python -c "from importlib.resources import files; from genbank_parser import __version__; from genbank_parser.meor import load_meor_database; from genbank_parser.mobilome import load_mobilome_database; db=load_mobilome_database(); root=files('genbank_parser'); assert __version__ == '0.6.5'; assert load_meor_database().catalog_version == '1.1'; assert db.catalog_version == '1.2.0'; assert db.provenance_version == '1.2.0'; assert db.inference_version == '1.2.0'; assert root.joinpath('py.typed').is_file(); assert root.joinpath('data/meor/markers.yaml').is_file(); assert files('genbank_parser.rulesets').joinpath('mobilome.yaml').is_file(); assert root.joinpath('data/mobilome/markers.yaml').is_file(); assert root.joinpath('data/mobilome/provenance.yaml').is_file(); assert root.joinpath('data/mobilome/inference.yaml').is_file(); assert root.joinpath('data/mobilome/report.schema.json').is_file()"
         .wheel-venv/bin/gbparse --help
         .wheel-venv/bin/gbparse meor tests/fixtures/meor_parity.gb --format json --output wheel-meor.json
-        .wheel-venv/bin/python -c "import json; p=json.load(open('wheel-meor.json')); assert p['tool_version'] == '0.6.0'; assert p['catalog_version'] == '1.1'"
+        .wheel-venv/bin/python -c "import json; p=json.load(open('wheel-meor.json')); assert p['tool_version'] == '0.6.5'; assert p['catalog_version'] == '1.1'"
         .wheel-venv/bin/gbparse mobilome --help
         .wheel-venv/bin/gbparse mobilome tests/fixtures/mobilome_evidence.gb --format json --output wheel-mobilome.json
-        .wheel-venv/bin/python -c "import json; p=json.load(open('wheel-mobilome.json')); assert p['schema_version'] == 'gbparse.mobilome.v1'; assert p['tool_version'] == '0.6.0'; assert p['catalog_version'] == '1.1.0'"
+        .wheel-venv/bin/python -c "import json; p=json.load(open('wheel-mobilome.json')); assert p['schema_version'] == 'gbparse.mobilome.v1'; assert p['tool_version'] == '0.6.5'; assert p['catalog_version'] == '1.2.0'"
         .wheel-venv/bin/gbparse discover tests/fixtures/mobilome_evidence.gb --ruleset mobilome
         .wheel-venv/bin/python -m pip check
 
```

### 4.8 changelogs.md — 0.6.5 entry

```diff
diff --git a/changelogs.md b/changelogs.md
index 5f2f10f..717bd9e 100644
--- a/changelogs.md
+++ b/changelogs.md
@@ -4,6 +4,36 @@ All notable changes to the `WhyAdr/genbank-parser` codebase are documented in th
 
 ---
 
+## [0.6.5] - 2026-09-06
+
+### Mobilome knowledge-base expansion and citation re-anchor
+
+- Expanded the mobilome provenance knowledge base from 16 to 28 sources: Bakta
+  (Schwengers 2021), plasmid extended mobility (Garcillan-Barcia 2025), ICEs
+  (Johnson and Grossman 2015), integrons (Mazel 2006), toxin-antitoxin
+  classification (Qiu 2022), fused MNT-HEPN regulation (Yao 2026), VFDB 2022
+  (Liu 2022), MGE-AMR association (Partridge 2018), PHASTER, geNomad, ISEScan,
+  and MobileElementFinder, all verified through Crossref, PubMed, and Unpaywall
+  during the audit session.
+- Re-anchored the AMRFinderPlus and VFDB citations that the 0.6.2 commit had
+  landed on `phage_packaging_candidate` and `phage_structure_candidate`:
+  `generic_amr_candidate` now carries `feldgarden-2021-amrfinderplus` and
+  `partridge-2018-mge-amr`, and `generic_vf_candidate` now carries
+  `chen-2005-vfdb` and `liu-2022-vfdb`; the two phage markers are back to
+  local-policy-only sources, with a marker-level anchor regression test.
+- Anchored the mobility and toxin-antitoxin markers and rules to the new
+  sources and documented the modern Type VII classification of HEPN/MNT-type
+  systems in the HEPN/MNT rule limitations without assigning a type number to
+  any annotation match.
+- Updated the mobilome evidence reference documentation with new source
+  correspondence bullets and primary descriptions for the external handoff
+  tools.
+- Bumped the mobilome catalog, provenance, and inference resources to 1.2.0,
+  moved the package version to 0.6.5 (resolving the drift between the 0.6.1 /
+  0.6.2 changelog entries and the 0.6.0 pyproject version), updated the wheel
+  CI assertions, regenerated the JSON, TSV, and text goldens, and extended the
+  citation-integrity tests.
+
 ## [0.6.2] - 2026-08-23
 
 ### Mobilome source-science audit
```

### 4.9 docs/mobilome_evidence_reference.md — new correspondence bullets and handoff citations

```diff
diff --git a/docs/mobilome_evidence_reference.md b/docs/mobilome_evidence_reference.md
index 7d630d1..eb67121 100644
--- a/docs/mobilome_evidence_reference.md
+++ b/docs/mobilome_evidence_reference.md
@@ -20,18 +20,24 @@ Generic Rep, pXO-numbered products, and phage-module observations do not select
 
 ## Source correspondence
 
-- ToxI is represented as the RNA antitoxin in a Type III ToxIN correspondence, following [Blower et al. (2012)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3401426/).
-- The HEPN/MNT reference retains the reported MntA/HepT biochemical correspondence without generalizing it to every annotation match; see [Yao et al. (2020)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7641770/).
-- Helper-dependent mobilization wording is limited by the distinction between annotated component co-occurrence and demonstrated transfer described by [Smillie et al. (2010)](https://pmc.ncbi.nlm.nih.gov/articles/PMC2937521/).
+- ToxI is represented as the RNA antitoxin in a Type III ToxIN correspondence, following [Blower et al. (2012)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3401426/); the pair-rule wording also carries the general toxin-antitoxin classification context of [Qiu et al. (2022)](https://doi.org/10.1016/j.micres.2022.127159).
+- The HEPN/MNT reference retains the reported MntA/HepT biochemical correspondence without generalizing it to every annotation match; see [Yao et al. (2020)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7641770/). Current classification reviews and the fused MNT-HEPN study of [Yao et al. (2026)](https://pmc.ncbi.nlm.nih.gov/articles/PMC13000451/) document HEPN/MNT-type systems among Type VII toxin-antitoxin systems; the rule text deliberately assigns no type number to an annotation match.
+- Helper-dependent mobilization wording is limited by the distinction between annotated component co-occurrence and demonstrated transfer described by [Smillie et al. (2010)](https://pmc.ncbi.nlm.nih.gov/articles/PMC2937521/); the modern extended-mobility framing of [Garcillán-Barcia et al. (2025)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12282955/) anchors the same caution without narrowing mobility to canonical conjugation.
 - The external MOB-suite or PlasmidFinder handoff requires versioned database and threshold provenance; see [Robertson and Nash (2018)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6159552/) for MOB-suite performance and database/taxonomic limitations.
 - AimR/AimP/AimX labels are retained as annotation correspondence only; [Erez et al. (2017)](https://pmc.ncbi.nlm.nih.gov/articles/PMC5378303/) does not make a single new-genome annotation functional evidence.
+- Integrative and conjugative element context follows [Johnson and Grossman (2015)](https://pmc.ncbi.nlm.nih.gov/articles/PMC5180612/); integration annotations in one record do not delimit an ICE.
+- Integron correspondence follows [Mazel (2006)](https://doi.org/10.1038/nrmicro1462); an integron-like annotation is not a cassette-array inventory without sequence analysis.
+- Antimicrobial-resistance gene mobility context follows [Partridge et al. (2018)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6148190/); association in that review does not establish that a generic resistance-like annotation is mobile.
+- Virulence-factor determination context follows [Chen et al. (2005)](https://academic.oup.com/nar/article/33/suppl_1/D325/2505203) and the VFDB 2022 classification scheme of [Liu et al. (2022)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8728188/); version provenance is required for any comparison.
+- Bakta output and database-version context is cited from [Schwengers et al. (2021)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8743544/) alongside the [Bakta documentation](https://bakta.readthedocs.io/en/latest/); observed typed input features take precedence over general tool claims.
 - INSDC feature and source-qualifier meanings are taken from the [INSDC feature table](https://www.ncbi.nlm.nih.gov/genbank/feature_table/); observed input feature representation remains authoritative.
-- Bakta capabilities and output forms are version-sensitive; use the observed typed record and feature annotations alongside the [Bakta documentation](https://bakta.readthedocs.io/en/latest/).
 
 ## External handoffs
 
 The report never runs external programs. It records `not_run` handoffs for MOB-suite or PlasmidFinder, AMRFinderPlus, a curated virulence workflow, prophage/ICE callers, and HMM/domain or sequence workflows. A future imported result must retain its tool/database release, method, thresholds, coverage, taxonomic scope, and limitations.
 
+Primary descriptions for the handoff tools now carried in provenance: MOB-suite ([Robertson and Nash, 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6159552/)), PlasmidFinder ([Carattoli et al., 2014](https://journals.asm.org/doi/10.1128/aac.02412-14)), AMRFinderPlus ([Feldgarden et al., 2021](https://www.nature.com/articles/s41598-021-91456-0)), VFDB ([Chen et al., 2005](https://academic.oup.com/nar/article/33/suppl_1/D325/2505203); [Liu et al., 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC8728188/)), prophage and element boundary callers such as PHASTER ([Arndt et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4987931/)) and geNomad ([Camargo et al., 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11324519/)), and sequence-level mobile-element callers such as ISEScan ([Xie and Tang, 2017](https://doi.org/10.1093/bioinformatics/btx433)) and MobileElementFinder ([Johansson et al., 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC7729385/)). External handoffs carry no source_ids in schema v1, so docs plus provenance is the v1-compatible home; a v2 schema could add handoff source_ids.
+
 Records without a detected annotation pipeline retain a `null` pipeline name in
 the schema-versioned `annotation_pipelines` array; this explicitly records that
 the input did not declare a recognized pipeline rather than inferring one.
```

### 4.10 Goldens — regenerate, do not hand-edit

`tests/golden/mobilome_v1.json`, `.tsv`, and `.txt` change because the
tool_version, catalog/inference versions, per-hit `source_ids`, and the
HEPN/MNT rule limitation list all appear in serialized output. Regenerate with
the command in section 7 step 3 after reinstalling the bumped package (the
version is read through `importlib.metadata`, so an editable install must be
refreshed first). Regeneration was verified byte-reproducible in a pristine
clone during this audit. `gbparse discover --ruleset mobilome` output is
untouched by this series.

## 5. Reference materials worth exploring (verified or partially verified, deliberately not in the diff)

These remain out of the catalog for the same reason patch 0.3 deferred its
section-5 items: each changes marker coverage or classification positions and
needs maintainers' scientific adjudication, fixtures, and adversarial
negatives first.

- **tenpIN / cptIN Type III superfamilies** — Blower et al. 2012 (already
  cited) classifies three Type III superfamilies: toxIN, tenpIN, and cptIN.
  The catalog only carries toxN/toxI (plus the HEPN/MNT pair). Observation-only
  markers for tenpN/tenpI and cptN/cptI would complete the already-cited
  correspondence without adding aggregate hypotheses. The sketch in
  mobilome-patch-0.3.md §5 remains the starting point; add fixtures and
  adversarial negatives before staging.
- **virB-family and trb MPF gene names** — the mpf_component gene regex
  `^tra(?:A|B|C|E|F|K|L|N|U|W)$` does not match `virB1`-`virB11` or `trbB`-
  `trbW`, both canonical mating-pair-formation gene names in Ti-plasmid-like
  and IncP-type records; those records currently match MPF only through
  product text. A regex extension plus a fixture with virB/trb genes and a
  negative (e.g. virB8-like photosystem genes) would close finding C7's MPF
  half.
- **Subtype-suffixed cas genes** — the crispr_cas regex
  `^cas(?:1|2|...|14|A|B|C|D|E|F)$` misses `cas8f`, `cas12a`, `cas13d`, and
  similar subtype-suffixed names standardized by Makarova 2020 (already
  cited). A pattern like `^cas\d+[a-z]?$` plus explicit letters would widen
  coverage; needs subtype-negative fixtures first.
- **repX / repL replication-initiation names** — the
  replication_initiation_candidate gene list (repA, repB, repS, repF, repFR)
  does not include repX (the pXO1 tubulin-like initiator, Paredes et al.
  2005, 10.1016/j.jmb.2005.01.056, not yet verified this session) or repL.
  Requires verification plus fixtures before staging; the mechanism-neutrality
  wording must not change.
- **Bakta Web (2025)** — Beyvers et al., "Bakta Web – rapid and standardized
  genome annotation on scalable infrastructures" (NAR 2025,
  10.1093/nar/gkaf335, verified through Crossref this session) documents the
  web front end. Not anchored because this tool analyzes flatfile exports, not
  web submissions; cite it if a future provenance section covers pipeline
  families beyond desktop Bakta.
- **VirSorter2** — Guo et al. 2021, "VirSorter2: a multi-classifier,
  expert-guided approach to detect diverse DNA and RNA viruses" (Microbiome,
  10.1186/s40168-020-00990-y, verified, OA, PMID 33522966). Redundant with
  geNomad and PHASTER for the boundary handoff's current wording; keep in
  reserve for a v2 handoff citation set.
- **Handoff source_ids in a v2 schema** — the fail-closed handoff parser
  rejects unknown keys, so citations for MOB-suite, PlasmidFinder,
  AMRFinderPlus, VFDB, PHASTER, geNomad, ISEScan, and MobileElementFinder
  live in docs plus provenance for now (finding A3). A v2 schema field would
  move them into the report itself.
- **Topology-neutral TA gap wording** — finding C2; rename the emitted
  limitation text (not the YAML key, which is a compatibility contract) once
  a linear-topology fixture exists.

## 6. Claim-wording guardrails for this expansion

- New supports/limitations text passes the repository forbidden-phrase
  contract across markers.yaml, inference.yaml, and provenance.yaml (the scan
  in tests/test_mobilome_patch_01.py enforces this; the diffs in section 4
  were verified against it locally, including the Type VII limitation line).
- No source was added without at least one stable identifier (DOI, PMID,
  PMCID, or URL) resolved during the same session that stages the change.
- Titles are copied verbatim from Crossref/PubMed, including slashes, commas,
  and hyphens ("toxin-antitoxin", "Integrative and Conjugative Elements
  (ICEs): What They Do and How They Work"); en-dash to ASCII hyphen
  substitution remains the one accepted normalization (Crossref prints
  Unicode en-dashes; the YAML stores ASCII hyphens).
- The Type VII documentation line describes the literature, not the
  annotation: the rule "deliberately assigns no type number to an annotation
  match", so no report output claims a TA type for any pair.
- Renaming or moving a source id requires updating every referencing
  marker/rule (grep first) plus goldens; this series only added ids and
  re-anchored two existing ones, and the new anchor-regression test pins the
  corrected placement.

## 7. Handoff to the executing agent (Gemini)

Context: you are executing a pre-verified patch series in the
`WhyAdr/genbank-parser` repository. The nine diffs in section 4 were generated
mechanically against baseline 16118c0 and passed `git apply --check`, and the
complete series plus regenerated goldens passed the full pytest suite (154
passed / 6 optional skips), scoped Ruff check and format, compileall, pip
check, and `git diff --check` in a pristine clone. Your job is to re-verify,
apply, regenerate, and commit — not to redesign.

0. Environment and baseline:

~~~bash
cd <repo>
git status --short --branch        # expect clean at 16118c0
python -m pip install -e ".[test]"
~~~

1. Re-verify every citation independently before staging. Do not skip this
   even though the audit already did it; titles must come from these
   responses, not memory:

~~~bash
for doi in \
  10.1099/mgen.0.000685 \
  10.1093/nar/gkaf652 \
  10.1146/annurev-genet-112414-055018 \
  10.1038/nrmicro1462 \
  10.1016/j.micres.2022.127159 \
  10.1093/nar/gkag228 \
  10.1093/nar/gkab1107 \
  10.1128/CMR.00088-17 \
  10.1093/nar/gkw387 \
  10.1038/s41587-023-01953-y \
  10.1093/bioinformatics/btx433 \
  10.1093/jac/dkaa390 ; do
  curl -s "https://api.crossref.org/works/$doi" | \
    python3 -c "import json,sys; d=json.load(sys.stdin)['message']; print(d['title'][0])"
done
# PMID spot checks (expect matching titles):
#   34739369 -> "Bakta: rapid and standardized annotation..."
#   40694848 -> "The extended mobility of plasmids"
#   26473380 -> "Integrative and Conjugative Elements (ICEs)..."
#   35969944 -> "Toxin-antitoxin systems: Classification..."
#   41854079 -> "Dual regulation of HEPN RNase in fused MNT-HEPN..."
~~~

   Acceptance: each printed title matches the corresponding title in the
   section 3.1 table and the section 4.1 diff, character for character except
   en-dash rendering (Crossref prints Unicode en-dashes; the YAML stores
   ASCII hyphens — that substitution is intentional and already adjudicated).
   Also confirm the patch-0.3 §5 Bakta DOI is dead (expect 404):
   `curl -s -o /dev/null -w "%{http_code}\n" https://api.crossref.org/works/10.1016/j.jbiotec.2021.07.018`

2. Extract and apply the nine diff blocks. Save each fenced block exactly as
   written (do not reflow), then:

~~~bash
for p in provenance.yaml markers.yaml inference.yaml \
         test_mobilome_database.py test_meor_report.py \
         pyproject.toml ci.yml changelogs.md docs.diff ; do
  git apply --check "gbparse-patch-0.6.5-$p.diff" || exit 1
done
for p in provenance.yaml markers.yaml inference.yaml \
         test_mobilome_database.py test_meor_report.py \
         pyproject.toml ci.yml changelogs.md docs.diff ; do
  git apply "gbparse-patch-0.6.5-$p.diff"
done
~~~

3. Reinstall the bumped package (the version is read through
   importlib.metadata, so the editable install must be refreshed) and
   regenerate the goldens (resource hashes, source_ids, versions, and the
   HEPN/MNT limitation changed):

~~~bash
python -m pip install -e .
python -c "
from pathlib import Path
from genbank_parser.mobilome import analyze_mobilome
from genbank_parser.mobilome.report import serialize_mobilome_report
report = analyze_mobilome(Path('tests/fixtures/mobilome_evidence.gb'))
for fmt, name in (('json','mobilome_v1.json'),('tsv','mobilome_v1.tsv'),('text','mobilome_v1.txt')):
    Path('tests/golden', name).write_bytes(serialize_mobilome_report(report, fmt).encode('utf-8'))
"
~~~

4. Gate (all must pass before staging):

~~~bash
MPLBACKEND=Agg python -m pytest -q -p no:cacheprovider
# expect: 154 passed, 6 skipped  (152 baseline + two new anchor tests)
python -m ruff check --no-cache src/genbank_parser/mobilome tests/test_mobilome_*.py
python -m ruff format --check --no-cache src/genbank_parser/mobilome tests/test_mobilome_*.py
python -m compileall -q src tests
python -m pip check
git diff --check
~~~

5. Stage exactly these paths and commit once (no git add . / -A):

~~~bash
git add -- \
  pyproject.toml \
  src/genbank_parser/data/mobilome/provenance.yaml \
  src/genbank_parser/data/mobilome/markers.yaml \
  src/genbank_parser/data/mobilome/inference.yaml \
  docs/mobilome_evidence_reference.md \
  tests/test_mobilome_database.py \
  tests/test_meor_report.py \
  tests/golden/mobilome_v1.json \
  tests/golden/mobilome_v1.tsv \
  tests/golden/mobilome_v1.txt \
  .github/workflows/ci.yml \
  changelogs.md \
  gbparse-patch-0.6.5-mobilome.md
git diff --cached --name-only     # must equal the list above
git commit -m "feat(mobilome): expand provenance knowledge base and re-anchor citations at 1.2.0"
~~~

6. Post-commit smoke: `gbparse mobilome tests/fixtures/mobilome_evidence.gb
   --format json | head -8` must report tool 0.6.5 and catalog/inference
   1.2.0, and `gbparse discover tests/fixtures/mobilome_evidence.gb
   --ruleset mobilome` must remain unchanged. Do not push; publication is
   separately authorized.

## 8. Gate for this patch set

Full pytest (154 expected: 152 baseline plus
`test_amr_and_vf_citations_are_anchored_to_their_own_markers` and
`test_ta_and_mobility_rules_carry_their_configured_sources`), scoped Ruff
check and format, compileall, pip check, `git diff --check`, regenerated
byte-reproducible goldens, unchanged discover output, tool version 0.6.5, and
the version triplet 1.2.0/1.2.0/1.2.0 across catalog/provenance/inference with
agreeing major versions. Catalog, provenance, and inference content changed by
addition plus a re-anchor of two misplaced citations — no contracts were
removed or reinterpreted — so 1.1.0 -> 1.2.0 minor bumps are the correct
semantic move. The deferred findings register (A3-A5, C1-C7) is deliberately
out of scope; a future v2 schema conversation is the right home for the
structural ones.

## Execution record

Applied 2026-09-06 against baseline `16118c0` by Gemini after independently
re-verifying all twelve new DOIs through Crossref, PubMed E-utilities, and
verifying the 404 response on the superseded Bakta DOI.

- Expanded mobilome provenance knowledge base from 16 to 28 verified sources.
- Fixed 0.6.2 citation mis-anchors on phage markers; re-anchored AMR and VF
  markers to AMRFinderPlus, VFDB, MGE-AMR, and VFDB 2022.
- Anchored mobility and toxin-antitoxin markers/rules to modern sources and
  documented modern Type VII classification for HEPN/MNT without type claims.
- Updated mobilome evidence reference documentation with new source correspondence
  bullets and external handoff citations.
- Bumped catalog, provenance, and inference versions to `1.2.0`, package version
  to `0.6.5` (resolving 0.6.0/0.6.2 drift), updated CI wheel assertions, and
  added changelog entry.
- Added marker anchor regression test and TA/mobility rule source test in
  `tests/test_mobilome_database.py`.
- Reinstalled package in editable mode and regenerated JSON, TSV, and text goldens.

Validation completed:
- Full test suite passed (160 passed with visualization dependencies installed,
  1 known Biopython malformed-LOCUS warning, 0 failures).
- Scoped Ruff check and format check passed cleanly on mobilome modules and tests.
- Python compileall passed without syntax errors.
- Pip check confirmed no broken dependencies.
- Git diff check passed with no whitespace errors.
- Mobilome and discover CLI smoke tests confirmed tool version 0.6.5, catalog and
  inference 1.2.0, and identical discover output.

