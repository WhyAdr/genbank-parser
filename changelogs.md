# Changelog: GenBank Parser & Annotation Engine

All notable changes to the `WhyAdr/genbank-parser` codebase are documented in this file.

---

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
