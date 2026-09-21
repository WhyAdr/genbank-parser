# Changelog: GenBank Parser & Annotation Engine

All notable changes to the `WhyAdr/genbank-parser` codebase are documented in this file.

---

## [0.9.1] - 2026-09-20

### Cohort hardening and provenance boundaries

- Made batch execution crash-resumable through a durable sibling in-progress
  tree, logical source-labelled replay argv, explicit required-output contracts,
  terminal outcome plus resume-action fields, and deterministic pruning/rekeying.
- Raised the persistent index to revision 3 with accepted-cohort key allocation,
  case-insensitive sample-key uniqueness, deterministic revision-2 migration,
  fail-closed structural validation, and side-effect-free WAL-preserving reads.
- Distinguished pre-commit rollback failure from post-commit backup-cleanup
  warnings while retaining recovery artifacts, and added complete Windows
  Python 3.12 release-suite and installed-schema instance gates.
- Preserved validation-threshold outcomes across verified batch resumes and
  validated manifest, environment, output, and stderr provenance before reuse.
- Keyed index sources and batch jobs by canonical resolved identity, with a
  revision-2 index migration and stable record-split filenames.
- Bound index and batch results to immutable source snapshots, implemented
  measured index parsing parallelism, and kept known-good indexes intact when
  no source was usable.
- Restored direct/index query type parity, complete canonical xref projections,
  bulk feature hydration, and cohort sample identity serialization.
- Guarded database/report and legacy batch-summary publication as coherent
  transactions, closed malformed zero-record and fail-policy paths, and
  expanded installed-wheel and cross-platform cohort gates.

---

## [0.9.0] - 2026-09-20

### Cohort records, persistent indexes, and resumable batch execution

- Added `gbparse records list|extract|filter|split` with schema-versioned
  record projections, exact ID/name selection, lossless GenBank round trips,
  gzip output, collision-safe split filenames, and an authoritative manifest.
- Added deterministic plain/gzip cohort discovery and raw-byte fingerprints,
  including symlink-safe duplicate suppression and collision-safe sample keys.
- Added normalized `gbparse.index.v1` SQLite build/update/query/status
  workflows with foreign-key/integrity validation, atomic replacement, and
  parameterized compilation of the existing bounded feature-query grammar.
- Added `gbparse batch` with explicit single-input adapters, isolated child
  argv execution, atomic per-job artifacts, `gbparse.batch.v1` manifests,
  hash-verified resume, and deterministic exit aggregation.
- Hardened `batch-summary` to share discovery, render in memory, publish
  atomically, report skipped inputs, and return nonzero for empty or wholly
  failed cohorts.
- Documented the v0.9.0 limits: no sequence archive, arbitrary SQL, shell
  execution, remote fetching, annotation reconciliation, or workflow-manager
  behavior.

---

## [0.8.5] - 2026-09-20

### CLI interoperability, safe queries, and complete evidence outputs

- Added source-aware plain/gzip/stdin parsing and shared atomic output
  publication with input/output collision protection and documented CLI exit
  statuses.
- Added safe declarative `query`, unified `export`, and structured,
  circular-aware `operons` commands. Query expressions use a bounded parser;
  they never evaluate Python, shell, or SQL expressions.
- Added canonical `gbparse.feature.v1` projections, JSONL, BED12, biological
  FAA/FFN/FNA, and candidate table2asn input exports. NCBI exports report
  unsupported qualifiers and retain the requirement for independent
  `table2asn` validation.
- Normalized summary, locus, CRISPR, phylogenomic, functional, and validation
  structured modes. `functional --format tsv` now emits real declared TSV
  rows, and `validate --fail-on` writes its report before returning status 1.
- Removed hidden machine-output truncation from annotation diffs and added
  typed compare evidence with count/presence/status modes.
- Preserved the existing MEOR, mobilome, neighborhood, region, and GFF3
  contracts. Scope exclusions are documented in
  `gbparse-v0.8.5-deferred-scope.md`.

## [0.8.2] - 2026-09-19

### TA/retron inference hardening and calibration reproducibility

- **Explicit spatial inference contract**: marker-cluster rules now declare
  `inference_mode: spatial_marker_cluster`; component rules cannot silently
  acquire marker-only fields, and legacy `max_circular_gap_bp` remains only a
  topology-neutral compatibility alias for one minor release cycle.
- **Injective compact-module selection**: tripartite and larger marker rules
  require deterministic marker-role-to-feature assignments, a globally
  coherent linear or circular covering extent, and deterministic
  non-overlapping compact modules. Valid submodules survive unrelated
  single-linkage bridges, while two-marker tandem arrays retain their
  documented aggregate behavior.
- **Structural policy framework**: custom catalogs can fail closed on strand,
  biological order, and typed intervening-feature policies; passing structural
  evidence is retained in hypothesis limitations. No speculative structural
  policy is activated in the packaged catalog.
- **Retron semantic split**: RT-msr/msd evidence is emitted as the tentative
  `retron_core_annotation_candidate` with `kind: retron`; the future retron TA
  aggregate remains disabled until a curated effector-correspondence contract
  exists. This replaces the 0.8.1
  `retron_rt_msdna_module_candidate` identifier without emitting a compatibility
  duplicate.
- **Calibration provenance**: the optional real-genome check is user-facing as
  a real-genome calibration, verifies the local fixture against a committed
  SHA-256 manifest, and explicitly does not claim a biological truth-set
  benchmark. The manifest distinguishes the local record ID from a public
  accession and marks the unreconstructable transformation provenance as
  incomplete. Mobilome v2 golden files are LF-stable.
- **Residual hardening**: compact-module assignment and packing guards now
  apply independently within connected marker clusters and surface skipped
  clusters as record-level spatial limitations. Circular structural policies
  follow the selected minimum covering arc without inspecting the complementary
  chromosome arc, and real v0.8.1 schema-1 custom inference catalogs migrate
  in memory during the v0.8.2 compatibility window.
- **Version bumps**: package 0.8.2 and mobilome inference resource 1.6.0;
  catalog and provenance resources remain 1.5.0.

## [0.8.1] - 2026-09-07

### Mobilome Spatial Clustering, Tripartite Modules, and repL Curation

- **C4 cluster engine**: Replaced the O(N x M) cross-product toxin-antitoxin
  pair emission with deterministic single-linkage marker clusters. A tandem
  array now collapses into one aggregate cluster hypothesis (with the observed
  cluster span and an explicit per-copy-pairing limitation) instead of one
  hypothesis per eligible pair; isolated toxins without a cognate still emit
  nothing; two-feature clusters keep the previous hypothesis ids, wording,
  and observed-gap limitations, so existing single-pair reports are unchanged.
- **C5 tripartite modules**: `required_marker_ids` now supports N-marker
  spatial clusters (N >= 2). Added rules for TAC modules
  (`higba_tac_module_candidate`: higb/higa/tac_chaperone), the
  omega-epsilon-zeta module (`omega_epsilon_zeta_module_candidate`), MqsRAC
  (`mqsrac_module_candidate`), and retron RT-msDNA modules
  (`retron_rt_msdna_module_candidate`), plus the Type III TenpIN and CptIN
  pair rules (`tenpn`/`tenpi`, `cptn`/`cpti`).
- **Marker expansion**: 15 new markers (tenpn, tenpi, cptn, cpti, zeta,
  epsilon, omega, mqsr, mqsa, mqsc, higb, higa, tac_chaperone, retron_rt,
  msdna) and the new `ta_accessory` facet for chaperone and regulator
  components.
- **repL dedicated curation (F5)**: `replication_initiation_candidate` now
  cites six verified RepL references (Catchpole 1992 PMID 1577254 and 1991
  PMID 1906970, Projan 1987 PMID 2822666, Khan and Novick 1982 PMID 7056699,
  Kwong 2017 PMID 29218034, Šprincová 2005 PMID 15907537); repL stays
  replication-mechanism-neutral.
- **Print-year update (F6)**: `christie-2025-t4ss` provenance carries the
  2026 print year (volume 50, online-first 2025-12-31) under its stable id.
- **Real-genome benchmark**: `scripts/benchmark_mobilome.py` and the opt-in
  `tests/test_mobilome_real_genome.py` check the integration-plan section
  16.9 calibration against a local, gitignored `PF_NNT_reoriented.gbff`,
  skipping clearly when absent and never staging the file. The standalone
  mobilome-parser prototype is documented as archived in its own repository.
- **Repository hygiene (F9)**: Removed the three paywalled full-text PDFs
  from version control; local copies now live in the gitignored
  `reference_pdfs/` directory and the catalog cites published records.
- **Fail-closed loader hardening**: toxin-antitoxin rules must declare at
  least two distinct `required_marker_ids` (duplicates already rejected).
- **Knowledge base**: 16 new verified provenance sources (six repL, ten
  tripartite-TA) with corrected identifiers: the deferred-review drafts
  mis-cited Bordes 2011 (correct: PNAS 108(20):8438-8443, PMID 21536872),
  Camacho 2002 (correct: Biol Chem 383(11):1701-1713, PMID 12530535), and
  Bobonis 2022 (correct PMID 35850148, not 35978194).
- **Version bumps**: package 0.8.1; mobilome catalog, provenance, and
  inference resources 1.5.0.

## [0.8.0] - 2026-09-06

### Mobilome Report Schema v2 and Hardening

- **Breaking report-contract repair**: Bumped the mobilome report contract to
  `gbparse.mobilome.v2`. The 0.7.0 release had introduced required fields
  (`handoffs[].source_ids`, `inventory[].spatial_limitations`) and a structured
  `missing_components` object array under the unchanged v1 schema version, so
  v1-labeled reports could no longer validate against the packaged v1 schema.
  Goldens were renamed to `mobilome_v2.*` and regenerated.
- **Text-report regression fix**: `RepliconInventory.spatial_limitations` are
  now printed in the text renderer (0.7.0 had silently dropped them there) and
  are carried in the TSV through a new `spatial_limitations_json` column
  (59-column header).
- **Fail-closed toxin-antitoxin gap policy**: the loader now rejects
  `toxin_antitoxin` rules that declare required markers without
  `max_circular_gap_bp`, preventing silently disabled pairing rules.
- **Knowledge-base addition**: added `makarova-2025-crispr-classification`
  (Nature Microbiology, PMID 41198952) anchoring subtype-suffixed Cas
  nomenclature, and extended the `crispr_cas_candidate` product matcher to
  recognize suffixed Cas names consistently with the gene matcher.
- **repL provenance scoping**: documented that the repX-specific citations
  anchor the repX pattern only; repL remains a local-policy candidate.
- **Test hardening**: regression tests for the linear-topology TA gap wording,
  scanner-level behavior of the expanded cas/virB/trb/repX regexes,
  retained-evidence facet/marker validation, handoff source resolution, and
  citation anchors for the 1.4.0 provenance additions.
- **Version bumps**: package 0.8.0; mobilome catalog, provenance, and
  inference resources 1.4.0.

## [0.7.0] - 2026-09-06

### Mobilome Catalog Hygiene, Marker Expansions, and Schema Decoupling

- **C2 Topology-aware TA gap wording**: Updated TA pair hypothesis limitation text in
  `inference.py` to reflect actual replicon topology (`Configured maximum gap ({inventory.topology}): ...`).
- **A4 Generic retained-evidence validator**: Removed hardcoded `pxo_numbered_product_annotation`
  whitelist in `database.py::_parse_disabled_rules`, allowing retained-evidence tokens from
  the union of `facet_ids | marker_ids`.
- **C7.2 Cas gene regex expansion**: Updated Cas regex in `markers.yaml` to Python-valid
  pattern `(?i)^cas(?:(?:[1-9]|1[0-4])[a-z]?|[A-F])$`, recognizing subtype-suffixed Cas names
  (e.g., `cas8f`, `cas12a`, `cas13d`).
- **C7.1 Conjugation MPF gene coverage**: Expanded `mpf_component` gene regex in `markers.yaml`
  with Agrobacterium Ti-plasmid `(?i)^virB(?:[1-9]|1[01])$` and IncP plasmid `(?i)^trb[B-L]$`.
- **C7.3 Replication initiation expansion**: Added `repX` and `repL` to
  `replication_initiation_candidate` in `markers.yaml`, accompanied by verified provenance
  entries for `anand-2008-repx` (PMID 18179418) and `tinsley-2006-repx` (PMID 16585744).
- **A3 External handoff source tracking**: Added `source_ids` to `ExternalHandoff` dataclass,
  JSON schema definition, `inference.yaml`, and JSON/TSV serializers.
- **A5 Structured missing components**: Standardized hypothesis `missing_components` with
  the `MissingComponent` dataclass, pairing component IDs with nested `missing_facets`
  across models, inference engine, reports, and schema.
- **Spatial limitations decoupling**: Isolated locus spatial-pairing warnings into
  `RepliconInventory.spatial_limitations`, keeping replicon-level classification limitations clean.
- **Catalog ID hygiene**: Disambiguated shadowed CRISPR marker IDs (`crispr_repeat_array` and
  `crispr_cas_candidate`) from their enclosing facet IDs and added validation in `database.py`
  to reject any marker ID that shadows an existing facet ID.
- **Version bumps and goldens**: Bumped package to 0.7.0, mobilome catalog, provenance, and
  inference resources to 1.3.0, updated CI wheel assertions, and regenerated versioned goldens.

## [0.6.5] - 2026-09-06

### Mobilome knowledge-base expansion and citation re-anchor

- Expanded the mobilome provenance knowledge base from 16 to 28 sources: Bakta
  (Schwengers 2021), plasmid extended mobility (Garcillan-Barcia 2025), ICEs
  (Johnson and Grossman 2015), integrons (Mazel 2006), toxin-antitoxin
  classification (Qiu 2022), fused MNT-HEPN regulation (Yao 2026), VFDB 2022
  (Liu 2022), MGE-AMR association (Partridge 2018), PHASTER, geNomad, ISEScan,
  and MobileElementFinder, all verified through Crossref, PubMed, and Unpaywall
  during the audit session.
- Re-anchored the AMRFinderPlus and VFDB citations that the 0.6.2 commit had
  landed on `phage_packaging_candidate` and `phage_structure_candidate`:
  `generic_amr_candidate` now carries `feldgarden-2021-amrfinderplus` and
  `partridge-2018-mge-amr`, and `generic_vf_candidate` now carries
  `chen-2005-vfdb` and `liu-2022-vfdb`; the two phage markers are back to
  local-policy-only sources, with a marker-level anchor regression test.
- Anchored the mobility and toxin-antitoxin markers and rules to the new
  sources and documented the modern Type VII classification of HEPN/MNT-type
  systems in the HEPN/MNT rule limitations without assigning a type number to
  any annotation match.
- Updated the mobilome evidence reference documentation with new source
  correspondence bullets and primary descriptions for the external handoff
  tools.
- Bumped the mobilome catalog, provenance, and inference resources to 1.2.0,
  moved the package version to 0.6.5 (resolving the drift between the 0.6.1 /
  0.6.2 changelog entries and the 0.6.0 pyproject version), updated the wheel
  CI assertions, regenerated the JSON, TSV, and text goldens, and extended the
  citation-integrity tests.

## [0.6.2] - 2026-08-23

### Mobilome source-science audit

- Corrected primary-article titles and the pXO1 provenance attribution, and
  added verified PMIDs for the existing mobility and arbitrium sources.
- Added primary literature anchors for relaxases, type IV secretion
  nomenclature, plasmid partition, CRISPR-Cas classification, insertion
  sequences, AMRFinderPlus, and VFDB.
- Bumped the mobilome catalog, provenance, and inference resources to 1.1.0;
  regenerated JSON, TSV, and text goldens and added citation-integrity tests.
- Kept the catalog annotation-supported and bounded: source evidence does not
  establish transfer, mechanism, phenotype, ancestry, or activity in a new
  genome.

## [0.6.1] - 2026-08-23

### Mobilome hardening

- Bumped the mobilome catalog and provenance resources to 1.0.1 after the
  patch-0.1 wording correction and MOB-suite provenance addition.
- Extended forbidden-claim enforcement through provenance resources and added
  database, CLI, toxin-antitoxin boundary, mobility-subset, and field-scoped
  scanner regressions.
- Recorded the patch-0.1 sequencing caveat: Phase 1 front-loaded the final
  catalog/provenance/inference content, Phase 5 pinned a hand-constructed
  empty-report golden, and the Phase 8 release included `report.py` and two
  test-file edits outside its nominal staging scope.

## [0.6.0] - 2026-08-23

### Typed mobilome evidence analysis

- Added native `gbparse mobilome` and `genbank_parser.mobilome`, reusing the canonical Biopython-backed `read_genbank()` parser for complete record inventory and typed feature identity.
- Added versioned marker, provenance, inference, and JSON Schema resources with exact-byte hashes; custom catalogs are loaded as one validated resource set.
- Added deterministic text, JSON, and normalized TSV reports, protected atomic output writing, record-level provenance, and explicit not-run handoffs for external typing workflows.
- Kept `gbparse discover --ruleset mobilome` unchanged as the annotation-island scanner; the new command is a separate replicon-centric evidence assessment.
- Added conservative component, toxin-antitoxin proximity, and cross-record helper-dependent hypotheses without mechanism, phenotype, transfer, ancestry, or element-identity conclusions.

## [0.5.0] - 2026-08-17

### Structured neighborhoods and optional visualization

- Added the pure `build_neighborhood()` API and schema-versioned text, TSV, and JSON outputs while retaining the printing compatibility wrapper.
- Hardened circular selection so oversized windows never duplicate a physical CDS; added origin-spanning compound-CDS unwrapping and circular last-to-first operon links.
- Added the optional `viz` extra and lazy dna_features_viewer renderer for headless SVG, PNG, and PDF output. The core wheel remains free of plotting dependencies.
- Added canonical mobilome/xenobiotics annotation-rule coloring, overlapping non-CDS context, and tested same-strand proximity overlays with evidence-calibrated legend wording.
- Added wheel-level visualization CI, a Python 3.10-3.13 core matrix, and parser-valid circular/context fixtures.
- Removed stale root ruleset copies and `requirements.txt`; `pyproject.toml` and packaged `src/genbank_parser/rulesets/` resources are authoritative.
- Deferred pyGenomeViz and comparative synteny until explicit alignment or orthology evidence is available.

## [0.4.1] - 2026-08-10

### MEOR identifier and reproducibility hardening
- Added independent MEOR catalog version `1.1` and emitted software/catalog versions in reports.
- Removed unresolved `K27540`, `K13060`, `K22363`, transferred EC `1.1.99.8`, and nondiscriminating `rhlRI` EC `2.3.1.184` from active evidence without guessing replacements.
- Downgraded wildcard ECs from structured Weight 3 to contextual Weight-1 note evidence.
- Replaced automatic audit `PASS` with `PRESENT_UNREVIEWED` and corrected marker-level KO/EC compatibility classification.
- Added a SHA-256 reference manifest and cross-platform deterministic audit paths.
- Moved the audit engine into `genbank_parser.meor.audit` with a thin script wrapper.
- Hardened MEOR clustering to require distinct physical genes on known strands and defined gap size as intervening bases.
- Added installed-wheel and package-data CI verification.

## [0.3.0] - 2026-08-09

### Breaking changes
- **Removed all legacy script wrappers** (`scripts/genbank_*.py`, `scripts/parse_bakta_summaries.py`, `scripts/genbank_parser.py`). The `gbparse` entry point is now the sole supported CLI. Use `pip install -e .` and invoke `gbparse <subcommand>` for all operations.
- Removed `scripts/genbank_parser.py` module-shadowing workaround (`sys.path` manipulation).
- Removed `test_legacy_scripts_execution` test.

## [0.2.1] - 2026-08-08

### Correctness and contract hardening
- Committed synthetic GenBank regression fixtures and raised the supported Python version to 3.10+.
- Locus lookup now prefers CDS features when paired `gene` and CDS records share a locus tag; neighborhood resolution no longer falls back to the first CDS.
- Region extraction always emits coordinates local to the returned sequence, preserves compound/reverse locations, records parent coordinates, and supports circular-origin windows.
- GFF3 export now respects `codon_start`, avoids duplicate hierarchy IDs, and marks origin-spanning compound parents.
- Compatibility strand access preserves `?` and `.`; coding density uses nonredundant CDS coverage.
- Validator reports total multi-record length, recognizes exceptional translation qualifiers, and reports unknown translation tables.
- Codon usage is translation-table aware and reports sense, stop, ambiguous, and excluded CDS counts.
- Discovery rules are package resources with functional text/JSON/TSV output and an observable `--operon-gap`; comparison, phylogenetic candidate, CRISPR/Cas, and diff contracts are safer and more explicit.

### Scope clarification
- Comparative matching remains annotation/xref based; phylogenetic output is candidate annotation matching; CRISPR output is an annotation scanner. Sequence similarity, HMM validation, and advanced evidence scoring remain deferred.

## [0.2.0] - 2026-08-08

### Major Architecture & Core Typed Data Model
- **Typed Model Layer (`src/genbank_parser/model.py`)**:
  - Implemented `@dataclass` models (`GenBankDocument`, `GenBankRecord`, `GenBankFeature`) preserving full Biopython `FeatureLocation` and `CompoundLocation` structures.
  - Preserved raw strand state ($+1, -1, 0, \text{None}$) with symbol mapping (`+`, `-`, `?`, `.`).
  - Added distinction between biological length (`len(location)`) and genomic span (`end - start + 1`).
  - Implemented biological sequence extraction via `feature.extract(record.seq)`.
  - Added backward-compatible dictionary-style item access on `GenBankFeature` (`feature['start']`, `feature['qualifiers']`, `feature['join_segments']`).
- **Canonical IO (`src/genbank_parser/io.py`)**:
  - Implemented `read_genbank()` as single source of truth for flatfile loading.
  - Enhanced `extract_xrefs()` with prioritized qualifier extraction for GO, COG, KEGG KO, Pfam, Rfam, and EC numbers (respecting INSDC `/EC_number` over Bakta `/note`).
- **Package Architecture (`pyproject.toml`)**:
  - Structured the codebase as a PEP 517/621 package under `src/genbank_parser/`.
  - Added `py.typed` PEP 561 marker.
  - Created entry point `gbparse = genbank_parser.cli:main`.
- **Legacy Script Compatibility** *(removed in 0.3.0)*:
  - Updated all 17 scripts in `scripts/` as lightweight CLI wrappers that imported from `genbank_parser`. These were removed in v0.3.0; use `gbparse` for all workflows.

---

### P0 Biological Correctness & Bug Fixes
- **GFF3 Exporter (`src/genbank_parser/gff.py`)**:
  - Fixed `join_segments` condition so ordinary single-segment features are never routed through compound parent/child splits.
  - Implemented accurate CDS translation phase calculation ($0, 1, 2$) across joined segments in 5' $\rightarrow$ 3' transcription order (properly handling reverse strand descending order).
  - Emitted full `##sequence-region` extents using true record lengths rather than min/max feature boundaries.
  - Added deterministic unique ID generation.
- **Validation Suite (`src/genbank_parser/validate.py`)**:
  - Replaced crude length estimation with biological translation checking against `transl_table` genetic codes and `codon_start` offsets.
  - Added pseudogene tolerance (`/pseudo`, `/pseudogene`) suppressing false missing translation errors.
  - Added structured severity levels (`ERROR`, `WARNING`, `INFO`) and machine-readable finding codes.
- **Sequence Extraction & Codon Bias (`src/genbank_parser/sequence.py`, `src/genbank_parser/codon.py`)**:
  - Sliced coding sequences biologically using `SeqFeature.extract()`, fixing inclusion of intervening non-coding gaps in joined/spliced CDSs.
  - Removed duplicate regex-based `ORIGIN` and `LOCUS` re-parsers across all modules.

---

### New Research Commands & Tools
- **Feature Query (`gbparse search` / `src/genbank_parser/query.py`)**:
  - Search annotations across genes, products, KOs, EC numbers, Pfam IDs, or regex patterns with text, TSV, CSV, or JSON output.
- **Region Extraction with Rebasing (`gbparse region` / `src/genbank_parser/region.py`)**:
  - Extract genomic sub-regions around locus tags with `--flank-genes`; Patch 2 makes returned feature coordinates local and retains `--rebase` as a compatibility flag.
- **Annotation Diff (`gbparse diff` / `src/genbank_parser/diff.py`)**:
  - Compare two annotation versions of the same genome (Bakta vs Prokka vs RefSeq), identifying identity-aware additions/removals, boundary shifts, product/gene renames, and KO/EC differences.
- **Declarative Rulesets (`src/genbank_parser/rulesets/*.yaml`)**:
  - Externalized mobilome and xenobiotic degradation discovery rules into YAML format.

---

### Testing & CI
- Created comprehensive `pytest` regression suite in `tests/` with 20 automated tests:
  - `tests/test_model.py`: Biological length vs span, strand symbols, dict compatibility.
  - `tests/test_parser.py`: Multi-record files, partial coordinates, xref parsing.
  - `tests/test_gff.py`: GFF3 compound vs ordinary feature export, CDS phases.
  - `tests/test_validate.py`: Translation verification, pseudogene tolerance, severity codes.
  - `tests/test_sequence_codon.py`: Biological extraction and codon usage calculation.
  - `tests/test_query_region_diff.py`: Search filters, region rebasing, annotation diffing.
  - `tests/test_cli.py`: Unified `gbparse` CLI subcommands.
- Added synthetic fixtures in `tests/fixtures/` (`simple_cds.gb`, `compound_joined.gb`, `special_cds.gb`, `multi_record_circular.gb`, `duplicate_locus.gb`).
- Added GitHub Actions CI workflow (`.github/workflows/ci.yml`) testing Python 3.10 through 3.13.
