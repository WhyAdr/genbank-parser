# Changelog: GenBank Parser & Annotation Engine

All notable changes to the `WhyAdr/genbank-parser` codebase are documented in this file.

---

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
- **Legacy Script Compatibility**:
  - Updated all 17 scripts in `scripts/` as lightweight CLI wrappers that import directly from `genbank_parser`, maintaining 100% backward compatibility for existing pipelines.

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
  - Extract genomic sub-regions around locus tags with `--flank-genes` and `--rebase` coordinate rewriting (for clinker, antiSMASH, and comparative loci).
- **Annotation Diff (`gbparse diff` / `src/genbank_parser/diff.py`)**:
  - Compare two annotation versions of the same genome (Bakta vs Prokka vs RefSeq), identifying added/removed CDSs, coordinate shifts, product/gene renames, and KO/EC differences.
- **Declarative Rulesets (`rulesets/*.yaml`)**:
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
  - `tests/test_cli.py`: Unified `gbparse` CLI subcommands and legacy script wrappers.
- Added synthetic fixtures in `tests/fixtures/` (`simple_cds.gb`, `compound_joined.gb`, `special_cds.gb`, `multi_record_circular.gb`, `duplicate_locus.gb`).
- Added GitHub Actions CI workflow (`.github/workflows/ci.yml`) testing Python 3.10 through 3.13.
