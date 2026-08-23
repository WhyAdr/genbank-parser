# Mobilome Patch 0.3 — Source Science Audit

Status: source-science audit plan for the mobilome evidence catalog. This
patch corrects misattributed and mis-titled citations, verifies every
biochemical correspondence claim, and anchors previously unanchored marker
gene nomenclature to primary literature.

Audit date: 2026-08-23

Planning baseline: main at 48ac990 ("harden mobilome patch 0.2").

Fence semantics (inherited from the integration plan): only blocks fenced as
`diff` below are applyable candidate patches. Every one of them was generated
mechanically from the baseline and passed `git apply --check` plus the full
repository gate (152 passed / 6 optional-visualization skips, scoped Ruff
check and format clean, goldens regenerated) before being embedded here.
Blocks fenced as `text` are adjudication sketches, not patches.

Verification method for this audit: every DOI below was resolved through
Crossref (`api.crossref.org/works/<DOI>`) and cross-checked against PubMed
records during the audit session. The executing agent must re-run the same
checks (commands in section 7) before committing; do not trust titles from
memory or from secondary snippets.

## 1. Re-audit of patch 0.2 (48ac990) — verified

- Version bumps to 1.0.1 across markers/provenance, CI, tests, and goldens:
  confirmed; live render byte-matches the committed goldens.
- provenance.yaml forbidden-phrase reword and scan extension: confirmed
  ("communication activity" now enforced across all three catalog resources).
- Duplicate-YAML-key loader, fail-closed database matrix, CLI residuals,
  1999/2000/2001 bp TA boundaries, lone-marker negatives, mobility subsets,
  same-record helper exclusion, db_xref scanner contracts: all present and
  passing (74 mobilome tests locally).
- The patch-0.2 execution record's claim that the README link was already
  correct is TRUE, and the original patch-0.1 finding was a false positive
  caused by an audit-environment text-display artifact that swallowed the
  `[m` of `[mobilome` in `the [mobilome evidence reference](...)`. The
  underlying bytes were always valid. Recorded here so the finding is not
  re-raised later.

## 2. Source science audit findings

Verdicts per primary-article entry in `src/genbank_parser/data/mobilome/provenance.yaml`:

| Source id | Verdict | Finding |
|---|---|---|
| blower-2012-toxin | PARTIAL | DOI/PMID/PMCID correct; title truncated. Full title ends "...encoded in chromosomal and plasmid genomes". |
| yao-2020-hepn-mnt | TITLE WRONG | Actual title: "Novel polyadenylylation-dependent neutralization mechanism of the HEPN/MNT toxin/antitoxin system". The stored title was an invented paraphrase. |
| yao-2020-hepn-mnt | CLAIM OK | "three AMP moieties ... HepT Tyr104" verified against the paper and two independent reviews; enriched below with the RX4HXY motif context. |
| smillie-2010-mobility | OK + GAP | Title matches PubMed's rendering; PMID 20805406 missing. |
| robertson-2018-mob-suite | TITLE WRONG | Actual title: "MOB-suite: software tools for clustering, reconstruction and typing of plasmids from draft assemblies". The stored title belongs to an unrelated paper; PMID 30052170 missing. |
| anjum-2014-pxo1 | MISATTRIBUTED | DOI/PMCID are correct, but the source id, title, and year are all wrong. DOI 10.1016/j.plasmid.2011.12.012 is Akhtar and Khan 2012, "Two independent replicons can support replication of the anthrax toxin-encoding plasmid pXO1 of Bacillus anthracis", PMID 22239982 (the stored PubMed URL pointed at PMID 22230551, also wrong). The integration plan cited this paper correctly; the implementation mangled it. |
| erez-2017-arbitrium | OK + GAP | Title/DOI/PMCID correct; PMID 28099413 missing. Supports claim accurate. |
| bakta-origin-docs, insdc-feature-table | OK | Documentation sources with required accessed dates. |
| local-* policy rules | OK | Internal rules; ids resolve. |

Marker-anchor gaps (correct claims, missing or indirect citations):

- `relaxase_candidate` (mobA/mobC/mobF gene set): cites Smillie 2010 for what
  is really relaxase-family nomenclature from Garcillan-Barcia et al. 2009.
- `t4cp_candidate` (traD/traG/virD4/tcpA) and `mpf_component`
  (traA|B|C|E|F|K|L|N|U|W): cite Smillie 2010; the MPF/Dtr component
  nomenclature now has a dedicated primary source (Christie 2025 unified
  nomenclature, FEMS Microbiol Rev).
- `crispr_cas` (cas gene regex): cites only Bakta documentation; the Cas
  classification is Makarova et al. 2020, Nat Rev Microbiol.
- `transposition_candidate`: cites only the local policy; IS nomenclature is
  Siguier et al. 2006 (ISfinder, NAR database issue).
- `partition_maintenance_candidate` (parA/parB/tubZ/tubR/copG/copA/copB):
  cites only the local policy; partition system components are Bouet and
  Funnell 2019 (EcoSal Plus).
- `generic_amr_candidate` / `generic_vf_candidate`: cite only the local
  policy; the structured-evidence alternatives they defer to (AMRFinderPlus
  Reference Gene Catalog; VFDB) have citable primary descriptions
  (Feldgarden et al. 2021; Chen et al. 2005).
- `plasmid_typing` handoff: MOB-suite tool paper is Robertson and Nash 2018
  (now correctly titled); PlasmidFinder is Carattoli et al. 2014. External
  handoffs carry no source_ids in schema v1, so these anchors belong in
  provenance and docs (see section 6).

## 3. Claims correctness review

- The HepT/MntA chemistry claim is correct and now records the RX4HXY motif
  context; the output wording still makes no TA-type-number claim, which
  remains deliberate (see section 6 for the Type VII note).
- The pXO supports statement ("pXO annotations alone do not specify a single
  replication mechanism") is a fair reading of Akhtar and Khan 2012, whose
  abstract states both the RepX and ORFs 14/16 replicons independently
  support full-length pXO1 replication.
- The Smillie-based helper-dependent wording and its five untested-condition
  limitations remain faithful to the source; no change needed.
- All new supports/limitations text below passes the repository
  forbidden-phrase contract (verified by running the catalog scan locally
  after applying the diffs).

## 4. Applyable patch series

All five diffs below are one atomic series against baseline 48ac990. Apply
in any order; verify per section 7. Resource versions move 1.0.1/1.0.0 ->
1.1.0 (catalog, provenance, inference all change; major versions agree).
Goldens must be regenerated after applying (command in section 7).

### 4.1 provenance.yaml — corrections and eight new primary anchors

```diff
diff --git a/src/genbank_parser/data/mobilome/provenance.yaml b/src/genbank_parser/data/mobilome/provenance.yaml
index 93ce3fa..ce95dd7 100644
--- a/src/genbank_parser/data/mobilome/provenance.yaml
+++ b/src/genbank_parser/data/mobilome/provenance.yaml
@@ -1,10 +1,10 @@
 schema_version: 1
-provenance_version: "1.0.1"
+provenance_version: "1.1.0"
 
 sources:
   - id: blower-2012-toxin
     type: primary_article
-    title: Identification and classification of bacterial Type III toxin-antitoxin systems
+    title: Identification and classification of bacterial Type III toxin-antitoxin systems encoded in chromosomal and plasmid genomes
     year: 2012
     doi: 10.1093/nar/gks231
     pmid: "22434880"
@@ -16,14 +16,14 @@ sources:
       - Annotation correspondence does not establish activity in a new genome.
   - id: yao-2020-hepn-mnt
     type: primary_article
-    title: HEPN/MNT toxin-antitoxin neutralization by polyadenylylation
+    title: "Novel polyadenylylation-dependent neutralization mechanism of the HEPN/MNT toxin/antitoxin system"
     year: 2020
     doi: 10.1093/nar/gkaa855
     pmid: "33045733"
     pmcid: PMC7641770
     url: https://pmc.ncbi.nlm.nih.gov/articles/PMC7641770/
     supports:
-      - In the studied system MntA transfers three AMP moieties from ATP to HepT Tyr104.
+      - In the studied system MntA transfers three AMP moieties from ATP to HepT Tyr104 in the RX4HXY motif.
     limitations:
       - The biochemical result must not be generalized to every HEPN/MNT annotation.
   - id: smillie-2010-mobility
@@ -31,6 +31,7 @@ sources:
     title: Mobility of Plasmids
     year: 2010
     doi: 10.1128/MMBR.00020-10
+    pmid: "20805406"
     pmcid: PMC2937521
     url: https://pmc.ncbi.nlm.nih.gov/articles/PMC2937521/
     supports:
@@ -39,21 +40,24 @@ sources:
       - Annotation co-occurrence does not establish transfer or co-transfer.
   - id: robertson-2018-mob-suite
     type: primary_article
-    title: Comprehensive assessment of plasmid prevalence and diversity in the complete genome collection
+    title: "MOB-suite: software tools for clustering, reconstruction and typing of plasmids from draft assemblies"
     year: 2018
     doi: 10.1099/mgen.0.000206
+    pmid: "30052170"
     pmcid: PMC6159552
     url: https://pmc.ncbi.nlm.nih.gov/articles/PMC6159552/
     supports:
       - MOB-suite plasmid typing and mobility results depend on versioned databases and operating thresholds.
     limitations:
       - Database composition, taxonomic coverage, and thresholds affect external plasmid-typing interpretation.
-  - id: anjum-2014-pxo1
+  - id: akhtar-2012-pxo1
     type: primary_article
-    title: Functional analysis of pXO1 plasmid replication
-    year: 2014
+    title: "Two independent replicons can support replication of the anthrax toxin-encoding plasmid pXO1 of Bacillus anthracis"
+    year: 2012
     doi: 10.1016/j.plasmid.2011.12.012
-    url: https://pubmed.ncbi.nlm.nih.gov/22230551/
+    pmid: "22239982"
+    pmcid: PMC4186245
+    url: https://pmc.ncbi.nlm.nih.gov/articles/PMC4186245/
     supports:
       - pXO annotations alone do not specify a single replication mechanism.
     limitations:
@@ -63,12 +67,101 @@ sources:
     title: Communication between viruses guides lysis-lysogeny decisions
     year: 2017
     doi: 10.1038/nature21049
+    pmid: "28099413"
     pmcid: PMC5378303
     url: https://pmc.ncbi.nlm.nih.gov/articles/PMC5378303/
     supports:
       - AimR, AimP, and AimX participate in the described phage correspondence.
     limitations:
       - One annotation does not establish communication activity in a new genome.
+  - id: bouet-2019-partition
+    type: primary_article
+    title: "Plasmid Localization and Partition in Enterobacteriaceae"
+    year: 2019
+    doi: 10.1128/ecosalplus.esp-0003-2019
+    pmid: "31187729"
+    url: https://journals.asm.org/doi/10.1128/ecosalplus.esp-0003-2019
+    supports:
+      - Plasmid partition (Par) systems are actively driven maintenance systems with ParA/ParB-type and tubulin-type components.
+    limitations:
+      - The review emphasizes Enterobacteriaceae; partition component naming elsewhere varies.
+  - id: carattoli-2014-plasmidfinder
+    type: primary_article
+    title: "In silico detection and typing of plasmids using PlasmidFinder and plasmid multilocus sequence typing"
+    year: 2014
+    doi: 10.1128/AAC.02412-14
+    pmid: "24777092"
+    url: https://journals.asm.org/doi/10.1128/aac.02412-14
+    supports:
+      - Replicon-based plasmid typing depends on a curated database with declared thresholds.
+    limitations:
+      - Replicon typing coverage is database- and lineage-dependent.
+  - id: chen-2005-vfdb
+    type: primary_article
+    title: "VFDB: a reference database for bacterial virulence factors"
+    year: 2005
+    doi: 10.1093/nar/gki008
+    pmid: "15608208"
+    url: https://academic.oup.com/nar/article/33/suppl_1/D325/2505203
+    supports:
+      - Virulence-factor determination is a curated-database comparison, not a product-text result.
+    limitations:
+      - Database releases supersede this original description; version provenance is required.
+  - id: christie-2025-t4ss
+    type: primary_article
+    title: "Type IV secretion systems: reconciling diversity through a unified nomenclature"
+    year: 2025
+    doi: 10.1093/femsre/fuaf069
+    pmid: "41474020"
+    url: https://academic.oup.com/femsre/article/doi/10.1093/femsre/fuaf069/8407627
+    supports:
+      - Mating-pair-formation (MPF) and Dtr component names follow a proposed unified nomenclature for type IV secretion systems.
+    limitations:
+      - A unified component name does not indicate a complete or active secretion system.
+  - id: feldgarden-2021-amrfinderplus
+    type: primary_article
+    title: "AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence"
+    year: 2021
+    doi: 10.1038/s41598-021-91456-0
+    pmid: "34135355"
+    url: https://www.nature.com/articles/s41598-021-91456-0
+    supports:
+      - Structured antimicrobial-resistance annotation requires a curated reference gene catalog with method and threshold provenance.
+    limitations:
+      - A reference-catalog hit is not a resistance phenotype measurement.
+  - id: garcillan-barcia-2009-relaxases
+    type: primary_article
+    title: "The diversity of conjugative relaxases and its application in plasmid classification"
+    year: 2009
+    doi: 10.1111/j.1574-6976.2009.00168.x
+    pmid: "19396961"
+    url: https://academic.oup.com/femsre/article/33/3/657/591359
+    supports:
+      - Relaxase (MOB) gene families underpin plasmid mobility classification; names such as mobA, mobC, and mobF denote relaxases.
+    limitations:
+      - A relaxase gene-name match alone is not a family assignment or phylogenetic result.
+  - id: makarova-2020-crispr
+    type: primary_article
+    title: "Evolutionary classification of CRISPR-Cas systems: a burst of class 2 and derived variants"
+    year: 2020
+    doi: 10.1038/s41579-019-0299-x
+    pmid: "31857715"
+    url: https://www.nature.com/articles/s41579-019-0299-x
+    supports:
+      - Cas genes follow an evolutionary classification into classes, types, and subtypes.
+    limitations:
+      - A cas gene-name match does not establish a complete or active CRISPR-Cas locus.
+  - id: siguier-2006-isfinder
+    type: primary_article
+    title: "ISfinder: the reference centre for bacterial insertion sequences"
+    year: 2006
+    doi: 10.1093/nar/gkj014
+    pmid: "16381877"
+    url: https://academic.oup.com/nar/article/34/suppl_1/D32/1132247
+    supports:
+      - Insertion sequences carry standardized family nomenclature in a dedicated reference center.
+    limitations:
+      - An IS-like annotation is not a family assignment without sequence comparison.
   - id: bakta-origin-docs
     type: official_documentation
     title: Bakta documentation
```

### 4.2 markers.yaml — re-anchor markers, rename pXO source, bump catalog

```diff
diff --git a/src/genbank_parser/data/mobilome/markers.yaml b/src/genbank_parser/data/mobilome/markers.yaml
index d571f25..8e07535 100644
--- a/src/genbank_parser/data/mobilome/markers.yaml
+++ b/src/genbank_parser/data/mobilome/markers.yaml
@@ -1,5 +1,5 @@
 schema_version: 1
-catalog_version: "1.0.1"
+catalog_version: "1.1.0"
 
 facets:
   - id: replication_origin
@@ -152,7 +152,7 @@ markers:
         strength: 2
     requires_non_pseudo: true
     requires_complete: true
-    sources: [smillie-2010-mobility]
+    sources: [garcillan-barcia-2009-relaxases, smillie-2010-mobility]
     interpretation: Retains a relaxase annotation candidate.
     limitations:
       - Annotation does not demonstrate nicking or transfer.
@@ -175,7 +175,7 @@ markers:
         strength: 2
     requires_non_pseudo: true
     requires_complete: true
-    sources: [smillie-2010-mobility]
+    sources: [christie-2025-t4ss, smillie-2010-mobility]
     interpretation: Retains a Type IV coupling-protein annotation candidate.
     limitations:
       - Annotation does not establish a complete transfer apparatus.
@@ -199,7 +199,7 @@ markers:
         strength: 2
     requires_non_pseudo: true
     requires_complete: true
-    sources: [smillie-2010-mobility]
+    sources: [christie-2025-t4ss, smillie-2010-mobility]
     interpretation: Retains a mating-pair-formation annotation component.
     limitations:
       - Multiple annotated components do not establish expression or transfer.
@@ -320,7 +320,7 @@ markers:
         strength: 1
     requires_non_pseudo: false
     requires_complete: false
-    sources: [local-annotation-candidate-policy]
+    sources: [bouet-2019-partition, local-annotation-candidate-policy]
     interpretation: Retains partition or plasmid-maintenance annotation evidence.
     limitations:
       - Observation does not establish stable inheritance.
@@ -387,7 +387,7 @@ markers:
         strength: 1
     requires_non_pseudo: false
     requires_complete: false
-    sources: [local-annotation-candidate-policy]
+    sources: [local-annotation-candidate-policy, siguier-2006-isfinder]
     interpretation: Retains transposition annotation evidence.
     limitations:
       - Observation does not establish transposition activity.
@@ -434,7 +434,7 @@ markers:
         strength: 2
     requires_non_pseudo: true
     requires_complete: true
-    sources: [bakta-origin-docs]
+    sources: [bakta-origin-docs, makarova-2020-crispr]
     interpretation: Retains a Cas annotation candidate.
     limitations:
       - Annotation does not establish a complete or active CRISPR system.
@@ -452,7 +452,7 @@ markers:
         strength: 1
     requires_non_pseudo: false
     requires_complete: false
-    sources: [anjum-2014-pxo1, local-disabled-aggregate-policy]
+    sources: [akhtar-2012-pxo1, local-disabled-aggregate-policy]
     interpretation: Retains a bounded pXO-numbered product annotation observation.
     limitations:
       - Observation is not a homology result, ancestry assignment, or mechanism call.
@@ -585,7 +585,7 @@ markers:
         strength: 1
     requires_non_pseudo: false
     requires_complete: false
-    sources: [local-annotation-candidate-policy]
+    sources: [feldgarden-2021-amrfinderplus, local-annotation-candidate-policy]
     interpretation: Retains a generic antimicrobial-resistance annotation candidate.
     limitations:
       - Product text alone does not establish resistance phenotype or determinant status.
@@ -603,7 +603,7 @@ markers:
         strength: 1
     requires_non_pseudo: false
     requires_complete: false
-    sources: [local-annotation-candidate-policy]
+    sources: [chen-2005-vfdb, local-annotation-candidate-policy]
     interpretation: Retains a generic virulence-associated annotation candidate.
     limitations:
       - Product text alone does not establish virulence phenotype.
```

### 4.3 inference.yaml — rename pXO source, bump inference version

```diff
diff --git a/src/genbank_parser/data/mobilome/inference.yaml b/src/genbank_parser/data/mobilome/inference.yaml
index 65f83de..fc72ca6 100644
--- a/src/genbank_parser/data/mobilome/inference.yaml
+++ b/src/genbank_parser/data/mobilome/inference.yaml
@@ -1,5 +1,5 @@
 schema_version: 1
-inference_version: "1.0.0"
+inference_version: "1.1.0"
 
 components:
   replication_candidate:
@@ -99,7 +99,7 @@ disabled_aggregate_rules:
       - Gap and interruption policy.
       - Committed positive and unrelated negative controls.
       - Sequence-homology provenance and reviewed operating threshold.
-    sources: [anjum-2014-pxo1, local-disabled-aggregate-policy]
+    sources: [akhtar-2012-pxo1, local-disabled-aggregate-policy]
   - id: phage_module_bearing_plasmid_candidate
     reason: No v1 distinct-module operating threshold has been validated.
     evidence_retained_as:
```

### 4.4 tests/test_mobilome_database.py — version asserts + citation-integrity regression

```diff
diff --git a/tests/test_mobilome_database.py b/tests/test_mobilome_database.py
index 13b1bf2..4f16409 100644
--- a/tests/test_mobilome_database.py
+++ b/tests/test_mobilome_database.py
@@ -26,9 +26,9 @@ def test_packaged_database_is_versioned_hashed_and_schema_valid() -> None:
     )
     schema = json.loads(schema_path.read_text(encoding="utf-8"))
 
-    assert database.catalog_version == "1.0.1"
-    assert database.inference_version == "1.0.0"
-    assert database.provenance_version == "1.0.1"
+    assert database.catalog_version == "1.1.0"
+    assert database.inference_version == "1.1.0"
+    assert database.provenance_version == "1.1.0"
     assert database.database_source == "packaged"
     assert database.source_paths == ()
     assert [resource.name for resource in database.resources] == [
@@ -41,6 +41,38 @@ def test_packaged_database_is_versioned_hashed_and_schema_valid() -> None:
     assert schema["$schema"] == "https://json-schema.org/draft/2020-12/schema"
 
 
+def test_provenance_sources_use_verified_citations() -> None:
+    database = load_mobilome_database()
+    by_id = {source.id: source for source in database.provenance_sources}
+
+    assert "anjum-2014-pxo1" not in by_id
+    akhtar = by_id["akhtar-2012-pxo1"]
+    assert akhtar.payload["year"] == 2012
+    assert akhtar.payload["title"].startswith("Two independent replicons")
+    assert (
+        "MOB-suite: software tools"
+        in by_id["robertson-2018-mob-suite"].payload["title"]
+    )
+    assert (
+        by_id["yao-2020-hepn-mnt"]
+        .payload["title"]
+        .startswith("Novel polyadenylylation")
+    )
+    assert (
+        by_id["blower-2012-toxin"]
+        .payload["title"]
+        .endswith("encoded in chromosomal and plasmid genomes")
+    )
+    for anchor in (
+        "garcillan-barcia-2009-relaxases",
+        "christie-2025-t4ss",
+        "makarova-2020-crispr",
+        "siguier-2006-isfinder",
+        "bouet-2019-partition",
+    ):
+        assert anchor in by_id
+
+
 def test_complete_custom_database_is_loaded_as_one_resource_set(tmp_path: Path) -> None:
     _copy_custom_database(tmp_path)
 
```

### 4.5 .github/workflows/ci.yml — wheel-gate version asserts

```diff
diff --git a/.github/workflows/ci.yml b/.github/workflows/ci.yml
index 3dd35e4..c635930 100644
--- a/.github/workflows/ci.yml
+++ b/.github/workflows/ci.yml
@@ -52,13 +52,13 @@ jobs:
 
     - name: Verify installed API, data, CLI, and dependencies
       run: |
-        .wheel-venv/bin/python -c "from importlib.resources import files; from genbank_parser import __version__; from genbank_parser.meor import load_meor_database; from genbank_parser.mobilome import load_mobilome_database; db=load_mobilome_database(); root=files('genbank_parser'); assert __version__ == '0.6.0'; assert load_meor_database().catalog_version == '1.1'; assert db.catalog_version == '1.0.1'; assert db.provenance_version == '1.0.1'; assert root.joinpath('py.typed').is_file(); assert root.joinpath('data/meor/markers.yaml').is_file(); assert files('genbank_parser.rulesets').joinpath('mobilome.yaml').is_file(); assert root.joinpath('data/mobilome/markers.yaml').is_file(); assert root.joinpath('data/mobilome/provenance.yaml').is_file(); assert root.joinpath('data/mobilome/inference.yaml').is_file(); assert root.joinpath('data/mobilome/report.schema.json').is_file()"
+        .wheel-venv/bin/python -c "from importlib.resources import files; from genbank_parser import __version__; from genbank_parser.meor import load_meor_database; from genbank_parser.mobilome import load_mobilome_database; db=load_mobilome_database(); root=files('genbank_parser'); assert __version__ == '0.6.0'; assert load_meor_database().catalog_version == '1.1'; assert db.catalog_version == '1.1.0'; assert db.provenance_version == '1.1.0'; assert db.inference_version == '1.1.0'; assert root.joinpath('py.typed').is_file(); assert root.joinpath('data/meor/markers.yaml').is_file(); assert files('genbank_parser.rulesets').joinpath('mobilome.yaml').is_file(); assert root.joinpath('data/mobilome/markers.yaml').is_file(); assert root.joinpath('data/mobilome/provenance.yaml').is_file(); assert root.joinpath('data/mobilome/inference.yaml').is_file(); assert root.joinpath('data/mobilome/report.schema.json').is_file()"
         .wheel-venv/bin/gbparse --help
         .wheel-venv/bin/gbparse meor tests/fixtures/meor_parity.gb --format json --output wheel-meor.json
         .wheel-venv/bin/python -c "import json; p=json.load(open('wheel-meor.json')); assert p['tool_version'] == '0.6.0'; assert p['catalog_version'] == '1.1'"
         .wheel-venv/bin/gbparse mobilome --help
         .wheel-venv/bin/gbparse mobilome tests/fixtures/mobilome_evidence.gb --format json --output wheel-mobilome.json
-        .wheel-venv/bin/python -c "import json; p=json.load(open('wheel-mobilome.json')); assert p['schema_version'] == 'gbparse.mobilome.v1'; assert p['tool_version'] == '0.6.0'; assert p['catalog_version'] == '1.0.1'"
+        .wheel-venv/bin/python -c "import json; p=json.load(open('wheel-mobilome.json')); assert p['schema_version'] == 'gbparse.mobilome.v1'; assert p['tool_version'] == '0.6.0'; assert p['catalog_version'] == '1.1.0'"
         .wheel-venv/bin/gbparse discover tests/fixtures/mobilome_evidence.gb --ruleset mobilome
         .wheel-venv/bin/python -m pip check
 
```

## 5. Reference materials worth exploring (verified, deliberately not yet in the diff)

These were resolved and verified during the audit but require maintainers'
scientific adjudication before entering the catalog, because each one either
changes marker coverage or documents a classification position:

- **tenpIN / cptIN Type III superfamilies** — Blower et al. 2012 (already
  cited) classifies three Type III superfamilies: toxIN, tenpIN, and cptIN.
  The catalog only carries toxN/toxI and hepT/mntA. Observation-only markers
  for tenpN/tenpI and cptN/cptI would complete the already-cited
  correspondence without adding aggregate hypotheses. Sketch (not a patch;
  requires fixture additions and adversarial negatives first):

~~~text
  - id: tenpn
    label: TenP-family toxin annotation
    facets: [ta_toxin]
    feature_types: [CDS, pseudogene]
    matchers:
      - field: gene
        mode: casefold_exact
        patterns: [tenpN]
        negative_patterns: []
        strength: 2
      - field: product
        mode: regex
        patterns:
          - "(?i)\\bTenP(?:-family)?\\b"
        negative_patterns: []
        strength: 2
    requires_non_pseudo: true
    requires_complete: true
    sources: [blower-2012-toxin]
    interpretation: Retains a TenP-family toxin annotation.
    limitations:
      - Annotation does not establish toxin activity.
~~~

- **HEPN/MNT is Type VII in current TA taxonomy** — the HEPN/MNT (HepT/MntA)
  system is now classified as a Type VII toxin-antitoxin system (e.g., the
  2026 fused MNT-HEPN regulation study, DOI 10.1093/nar/gkag228, and the Qiu
  2022 review, PMID 35969944). The output wording must still make no
  type-number claim, but a provenance entry documenting the modern
  classification, cited from the limitations text of the hepn_mnt rule in a
  later catalog bump, would prevent readers from assuming the pair is
  unclassified. Add the sources below only after re-verifying both records
  (Crossref for gkag228 confirmed "Dual regulation of HEPN RNase in fused
  MNT-HEPN toxin-antitoxin systems via protein OligoAMPylation and
  oligomerization", 2026).

- **Bakta primary citation** — the catalog cites only Bakta's readthedocs
  page. The tool paper is Schwengers et al. 2021, "Bakta: rapid and
  standardized annotation of bacterial genomes, plasmids and phage
  sequences" (verify the exact title tail via Crossref DOI
  10.1016/j.jbiotec.2021.07.018; the lookup was rate-limited during this
  audit), PMID 34739369. Adding it as `type: primary_article` alongside
  `bakta-origin-docs` strengthens the annotation-provenance section.

- **Handoff tool citations in docs** — `docs/mobilome_evidence_reference.md`
  external-handoffs paragraph should name the primary descriptions now
  available in provenance: MOB-suite (robertson-2018-mob-suite),
  PlasmidFinder (carattoli-2014-plasmidfinder), AMRFinderPlus
  (feldgarden-2021-amrfinderplus), and VFDB (chen-2005-vfdb). External
  handoffs carry no source_ids in schema v1, so docs plus provenance is the
  v1-compatible home; a v2 schema could add handoff source_ids.

- **Extended plasmid mobility framing** — Garcillan-Barcia et al. 2025, "The
  extended mobility of plasmids" (NAR, DOI 10.1093/nar/gkaf652) reframes
  mobilization beyond canonical conjugation. Relevant when a future version
  revisits the helper-dependent wording, not required for v1.

## 6. Claim-wording guardrails for any expansion

- New supports/limitations text must keep passing the forbidden-phrase
  contract across markers.yaml, inference.yaml, and provenance.yaml (the
  scan in tests/test_mobilome_patch_01.py enforces this; the diffs in
  section 4 were verified against it locally).
- No source may be added without at least one stable identifier (DOI, PMID,
  PMCID, or URL) resolved during the same session that stages the change.
- Titles are copied verbatim from Crossref/PubMed, including slashes and
  hyphens ("toxin/antitoxin", "lysis-lysogeny"), not paraphrased.
- Renaming a source id requires updating every referencing marker/rule
  (grep first) plus goldens; the section 4 series did this for
  anjum-2014-pxo1 -> akhtar-2012-pxo1.

## 7. Handoff to the executing agent (Gemini)

Context: you are executing a pre-verified patch series in the
`WhyAdr/genbank-parser` repository. The five diffs in section 4 were
generated mechanically against baseline 48ac990 and passed `git apply
--check`, the full pytest suite (152 passed / 6 optional skips), scoped
Ruff, and golden regeneration in a clean clone. Your job is to re-verify,
apply, regenerate, and commit — not to redesign.

0. Environment and baseline:

~~~bash
cd <repo>
git status --short --branch        # expect clean at 48ac990
python -m pip install -e ".[test]"
~~~

1. Re-verify every citation independently before staging. Do not skip this
   even though the audit already did it; titles must come from these
   responses, not memory:

~~~bash
for doi in \
  10.1093/nar/gks231 \
  10.1093/nar/gkaa855 \
  10.1128/MMBR.00020-10 \
  10.1099/mgen.0.000206 \
  10.1016/j.plasmid.2011.12.012 \
  10.1038/nature21049 \
  10.1128/ecosalplus.esp-0003-2019 \
  10.1128/AAC.02412-14 \
  10.1093/nar/gki008 \
  10.1093/femsre/fuaf069 \
  10.1038/s41598-021-91456-0 \
  10.1111/j.1574-6976.2009.00168.x \
  10.1038/s41579-019-0299-x \
  10.1093/nar/gkj014 ; do
  curl -s "https://api.crossref.org/works/$doi" | \
    python3 -c "import json,sys; d=json.load(sys.stdin)['message']; print(d['title'][0])"
done
# PMID spot checks (expect matching titles):
#   22239982 -> "Two independent replicons can support replication..."
#   30052170 -> "MOB-suite: software tools..."
#   33045733 -> "Novel polyadenylylation-dependent neutralization..."
~~~

   Acceptance: each printed title matches the corresponding title in the
   section 4.1 diff, character for character except for en-dash rendering
   (Crossref prints Unicode en-dashes; the YAML stores ASCII hyphens — that
   substitution is intentional and already adjudicated).

2. Extract and apply the five diff blocks. Save each fenced block exactly as
   written (do not reflow), then:

~~~bash
for p in provenance.yaml markers.yaml inference.yaml \
         test_mobilome_database.py ci.yml ; do
  git apply --check "mobilome-patch-0.3-$p.diff" || exit 1
done
for p in provenance.yaml markers.yaml inference.yaml \
         test_mobilome_database.py ci.yml ; do
  git apply "mobilome-patch-0.3-$p.diff"
done
~~~

3. Regenerate the goldens (resource hashes and source_ids changed):

~~~bash
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
# expect: 152 passed, 6 skipped  (151 baseline + the new citation-integrity test)
python -m ruff check --no-cache src/genbank_parser/mobilome tests/test_mobilome_*.py
python -m ruff format --check --no-cache src/genbank_parser/mobilome tests/test_mobilome_*.py
python -m compileall -q src tests
python -m pip check
git diff --check
~~~

5. Stage exactly these paths and commit once (no git add . / -A):

~~~bash
git add -- \
  src/genbank_parser/data/mobilome/provenance.yaml \
  src/genbank_parser/data/mobilome/markers.yaml \
  src/genbank_parser/data/mobilome/inference.yaml \
  tests/test_mobilome_database.py \
  tests/golden/mobilome_v1.json \
  tests/golden/mobilome_v1.tsv \
  tests/golden/mobilome_v1.txt \
  .github/workflows/ci.yml \
  mobilome-patch-0.3.md
git diff --cached --name-only     # must equal the list above
git commit -m "fix(mobilome): correct and anchor source-science provenance at 1.1.0"
~~~

6. Post-commit smoke: `gbparse mobilome tests/fixtures/mobilome_evidence.gb
   --format json | head -8` must report catalog/inference 1.1.0, and
   `gbparse discover tests/fixtures/mobilome_evidence.gb --ruleset mobilome`
   must remain unchanged. Do not push; publication is separately authorized.

7. Optional follow-up (separate commit, only if time allows and after
   verifying the Bakta record): the docs handoff-citation paragraph from
   section 5, touching only docs/mobilome_evidence_reference.md and
   provenance.yaml (which would then move to 1.1.1 with golden regen).

## 8. Gate for this patch set

Full pytest (152 expected after the new citation-integrity test), scoped
Ruff check and format, compileall, pip check, git diff --check, regenerated
goldens, unchanged discover output, and version triplet 1.1.0/1.1.0/1.1.0
across catalog/provenance/inference with agreeing major versions. Catalog
and provenance content changed, so 1.0.x -> 1.1.0 minor bumps are the
correct semantic move (additions plus corrections, no removals of contracts);
a future marker-set change that deletes or reinterprets evidence would
require the v2 schema conversation instead.

## Execution record

Applied 2026-08-23 against baseline `48ac990` after independently verifying
the listed DOI metadata through Crossref and PubMed records, and checking the
HEPN/MNT and pXO1 claims against their primary full-text articles.

- Corrected source titles, PMIDs, and pXO1 attribution; added eight primary
  literature anchors and re-anchored the affected markers and disabled rule.
- Bumped catalog, provenance, and inference versions to `1.1.0` and updated
  the wheel CI assertions.
- Added citation-integrity regression coverage and regenerated JSON, TSV, and
  text goldens.
- Updated `changelogs.md` with the 0.6.2 source-science-audit entry.

Validation completed: targeted database tests pass (23 passed), the full
repository suite passes (158 passed), and scoped Ruff/format, mypy,
compileall, pip check, wheel build/install, version/reference closure,
discover, and diff checks pass. The existing Biopython malformed-LOCUS test
warning remains; no new warning was introduced.
