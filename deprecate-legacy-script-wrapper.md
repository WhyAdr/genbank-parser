# Remove Legacy Script Wrappers — Make `gbparse` the Sole CLI

All 17 files in `scripts/` are thin wrappers that do nothing but manipulate `sys.path` and delegate to `genbank_parser.*` modules. The unified `gbparse` entry point (installed via `pip install -e .`) already covers every subcommand with richer argument handling. This plan removes the wrappers, the module-shadowing workaround, and all documentation/test references to legacy script paths.

## Proposed Changes

### 1. Delete Legacy Script Wrappers

#### [DELETE] [scripts/genbank_codon.py](file:///d:/W/Skills%20Claude/genbank-feature-parser/scripts/genbank_codon.py)
#### [DELETE] [scripts/genbank_compare.py](file:///d:/W/Skills%20Claude/genbank-feature-parser/scripts/genbank_compare.py)
#### [DELETE] [scripts/genbank_crispr.py](file:///d:/W/Skills%20Claude/genbank-feature-parser/scripts/genbank_crispr.py)
#### [DELETE] [scripts/genbank_discover.py](file:///d:/W/Skills%20Claude/genbank-feature-parser/scripts/genbank_discover.py)
#### [DELETE] [scripts/genbank_extract.py](file:///d:/W/Skills%20Claude/genbank-feature-parser/scripts/genbank_extract.py)
#### [DELETE] [scripts/genbank_fasta.py](file:///d:/W/Skills%20Claude/genbank-feature-parser/scripts/genbank_fasta.py)
#### [DELETE] [scripts/genbank_functional.py](file:///d:/W/Skills%20Claude/genbank-feature-parser/scripts/genbank_functional.py)
#### [DELETE] [scripts/genbank_gff.py](file:///d:/W/Skills%20Claude/genbank-feature-parser/scripts/genbank_gff.py)
#### [DELETE] [scripts/genbank_locus.py](file:///d:/W/Skills%20Claude/genbank-feature-parser/scripts/genbank_locus.py)
#### [DELETE] [scripts/genbank_metadata.py](file:///d:/W/Skills%20Claude/genbank-feature-parser/scripts/genbank_metadata.py)
#### [DELETE] [scripts/genbank_neighborhood.py](file:///d:/W/Skills%20Claude/genbank-feature-parser/scripts/genbank_neighborhood.py)
#### [DELETE] [scripts/genbank_operons.py](file:///d:/W/Skills%20Claude/genbank-feature-parser/scripts/genbank_operons.py)
#### [DELETE] [scripts/genbank_parser.py](file:///d:/W/Skills%20Claude/genbank-feature-parser/scripts/genbank_parser.py)
#### [DELETE] [scripts/genbank_phylo.py](file:///d:/W/Skills%20Claude/genbank-feature-parser/scripts/genbank_phylo.py)
#### [DELETE] [scripts/genbank_sequence.py](file:///d:/W/Skills%20Claude/genbank-feature-parser/scripts/genbank_sequence.py)
#### [DELETE] [scripts/genbank_validate.py](file:///d:/W/Skills%20Claude/genbank-feature-parser/scripts/genbank_validate.py)
#### [DELETE] [scripts/parse_bakta_summaries.py](file:///d:/W/Skills%20Claude/genbank-feature-parser/scripts/parse_bakta_summaries.py)
#### [DELETE] [scripts/\_\_pycache\_\_/](file:///d:/W/Skills%20Claude/genbank-feature-parser/scripts/__pycache__/) (directory)

The `scripts/` directory itself will be deleted since no developer/build utilities remain in it.

---

### 2. Update Tests

#### [MODIFY] [test_cli.py](file:///d:/W/Skills%20Claude/genbank-feature-parser/tests/test_cli.py)

Remove `test_legacy_scripts_execution` and update the module docstring.

```diff
-"""Test unified gbparse CLI and backwards-compatible scripts."""
+"""Test unified gbparse CLI subcommands."""
 from pathlib import Path
-import subprocess
-import sys

 from genbank_parser.cli import main
```

```diff
-def test_legacy_scripts_execution(simple_cds_gbff: Path) -> None:
-    scripts_dir = Path(__file__).resolve().parent.parent / "scripts"
-    val_script = scripts_dir / "genbank_validate.py"
-
-    proc = subprocess.run([sys.executable, str(val_script), str(simple_cds_gbff)], capture_output=True, text=True)
-    assert proc.returncode == 0
-    assert "GENBANK FEATURE TABLE -- STRUCTURAL & BIOLOGICAL REPORT" in proc.stdout
```

---

### 3. Update SKILL.md

#### [MODIFY] [SKILL.md](file:///d:/W/Skills%20Claude/genbank-feature-parser/SKILL.md)

Remove all "Legacy Script" references and the dual-invocation table. Replace with a `gbparse`-only command reference.

```diff
-## Purpose
-
-Parse, validate, and analyze GenBank flatfiles (`.gb`, `.gbk`, `.gbff`, `.txt`) **computationally** — never loading the full file content into the AI context window. Powered by Biopython `Bio.SeqIO` and structured into an installable Python package (`genbank_parser`) with a unified CLI (`gbparse`) and backwards-compatible scripts under `scripts/`.
+## Purpose
+
+Parse, validate, and analyze GenBank flatfiles (`.gb`, `.gbk`, `.gbff`, `.txt`) **computationally** — never loading the full file content into the AI context window. Powered by Biopython `Bio.SeqIO` and structured into an installable Python package (`genbank_parser`) with a unified CLI (`gbparse`).
```

```diff
-## Unified CLI (`gbparse`) & Script Inventory
-
-You can invoke commands either via the unified `gbparse` CLI or through the legacy scripts under `scripts/`:
-
-| `gbparse` Subcommand | Legacy Script | Purpose | Typical Invocation |
-|---|---|---|---|
-| `gbparse validate` | `genbank_validate.py` | Biological translation QC & structure report | `gbparse validate INPUT.gbff [--json]` |
-| `gbparse summary` | `genbank_metadata.py` | LOCUS metadata, contigs, topologies & length | `gbparse summary INPUT.gbff` |
-| `gbparse extract` | `genbank_extract.py` | Tab-delimited annotation TSV export | `gbparse extract INPUT.gbff [output.tsv]` |
-| `gbparse search` | (new) | Search features by gene, product, KO, EC, Pfam | `gbparse search INPUT.gbff --gene ladA --format tsv` |
-| `gbparse locus` | `genbank_locus.py` | Single-locus qualifier deep-dive | `gbparse locus INPUT.gbff LOCUS_TAG` |
-| `gbparse neighborhood` | `genbank_neighborhood.py` | Circular-aware flanking gene viewer (+/- N) | `gbparse neighborhood INPUT.gbff LOCUS_TAG [window]` |
-| `gbparse region` | (new) | Sub-region extraction with valid local coordinates | `gbparse region INPUT.gbff --locus TAG --flank-genes 5 --output region.gbk` |
-| `gbparse fasta` | `genbank_fasta.py` | Export all CDS translations as protein FASTA | `gbparse fasta INPUT.gbff [proteins.faa]` |
-| `gbparse sequence` | `genbank_sequence.py` | Extract genome FASTA (.fna) & CDS (.ffn) | `gbparse sequence INPUT.gbff [--fna out.fna] [--ffn out.ffn]` |
-| `gbparse codon` | `genbank_codon.py` | Codon usage bias, RSCU & positional GC | `gbparse codon INPUT.gbff [--min-len 100]` |
-| `gbparse functional` | `genbank_functional.py` | COG distribution + metabolic completeness | `gbparse functional INPUT.gbff [--format json]` |
-| `gbparse discover` | `genbank_discover.py` | Scan annotation-supported mobilome/xenobiotic islands | `gbparse discover INPUT.gbff [--ruleset mobilome] [--format tsv]` |
-| `gbparse compare` | `genbank_compare.py` | Multi-genome marker presence/absence matrix | `gbparse compare genomes/ --targets "ladA,ssuD,K20938"` |
-| `gbparse diff` | (new) | Compare two annotation versions of a genome | `gbparse diff old.gbff new.gbff [--format json]` |
-| `gbparse phylo` | `genbank_phylo.py` | Annotation-based candidate phylogenetic markers | `gbparse phylo INPUT.gbff [--markers all] [--min-length 50]` |
-| `gbparse crispr` | `genbank_crispr.py` | CRISPR/Cas annotation scanner | `gbparse crispr INPUT.gbff [--window 15000]` |
-| `gbparse gff` | `genbank_gff.py` | Standard GFF3 export (with CDS phase & regions) | `gbparse gff INPUT.gbff [output.gff3] [--include-fasta]` |
-| `gbparse batch-summary` | `parse_bakta_summaries.py`| Bakta multi-isolate comparison tables | `gbparse batch-summary ./isolates/ --csv summary.csv` |
+## CLI Reference (`gbparse`)
+
+| Subcommand | Purpose | Typical Invocation |
+|---|---|---|
+| `validate` | Biological translation QC & structure report | `gbparse validate INPUT.gbff [--json]` |
+| `summary` | LOCUS metadata, contigs, topologies & length | `gbparse summary INPUT.gbff` |
+| `extract` | Tab-delimited annotation TSV export | `gbparse extract INPUT.gbff [output.tsv]` |
+| `search` | Search features by gene, product, KO, EC, Pfam | `gbparse search INPUT.gbff --gene ladA --format tsv` |
+| `locus` | Single-locus qualifier deep-dive | `gbparse locus INPUT.gbff LOCUS_TAG` |
+| `neighborhood` | Circular-aware flanking gene viewer (+/- N) | `gbparse neighborhood INPUT.gbff LOCUS_TAG [window]` |
+| `region` | Sub-region extraction with valid local coordinates | `gbparse region INPUT.gbff --locus TAG --flank-genes 5 --output region.gbk` |
+| `fasta` | Export all CDS translations as protein FASTA | `gbparse fasta INPUT.gbff [proteins.faa]` |
+| `sequence` | Extract genome FASTA (.fna) & CDS (.ffn) | `gbparse sequence INPUT.gbff [--fna out.fna] [--ffn out.ffn]` |
+| `codon` | Codon usage bias, RSCU & positional GC | `gbparse codon INPUT.gbff [--min-len 100]` |
+| `functional` | COG distribution + metabolic completeness | `gbparse functional INPUT.gbff [--format json]` |
+| `discover` | Scan annotation-supported mobilome/xenobiotic islands | `gbparse discover INPUT.gbff [--ruleset mobilome] [--format tsv]` |
+| `compare` | Multi-genome marker presence/absence matrix | `gbparse compare genomes/ --targets "ladA,ssuD,K20938"` |
+| `diff` | Compare two annotation versions of a genome | `gbparse diff old.gbff new.gbff [--format json]` |
+| `phylo` | Annotation-based candidate phylogenetic markers | `gbparse phylo INPUT.gbff [--markers all] [--min-length 50]` |
+| `crispr` | CRISPR/Cas annotation scanner | `gbparse crispr INPUT.gbff [--window 15000]` |
+| `gff` | Standard GFF3 export (with CDS phase & regions) | `gbparse gff INPUT.gbff [output.gff3] [--include-fasta]` |
+| `batch-summary` | Bakta multi-isolate comparison tables | `gbparse batch-summary ./isolates/ --csv summary.csv` |
```

---

### 4. Update changelogs.md

#### [MODIFY] [changelogs.md](file:///d:/W/Skills%20Claude/genbank-feature-parser/changelogs.md)

Add a new `[0.3.0]` section at the top. Update the v0.2.0 legacy script section to reflect removal.

```diff
 ---

+## [0.3.0] - 2026-08-09
+
+### Breaking changes
+- **Removed all legacy script wrappers** (`scripts/genbank_*.py`, `scripts/parse_bakta_summaries.py`, `scripts/genbank_parser.py`). The `gbparse` entry point is now the sole supported CLI. Use `pip install -e .` and invoke `gbparse <subcommand>` for all operations.
+- Removed `scripts/genbank_parser.py` module-shadowing workaround (`sys.path` manipulation).
+- Removed `test_legacy_scripts_execution` test.
+
 ## [0.2.1] - 2026-08-08
```

Update the v0.2.0 legacy section to reflect historical context:

```diff
-- **Legacy Script Compatibility**:
-  - Updated all 17 scripts in `scripts/` as lightweight CLI wrappers that import directly from `genbank_parser`. Legacy paths remain available, but some 0.2.0 CLI argument contracts changed; prefer `gbparse` for new workflows.
+- **Legacy Script Compatibility** *(removed in 0.3.0)*:
+  - Updated all 17 scripts in `scripts/` as lightweight CLI wrappers that imported from `genbank_parser`. These were removed in v0.3.0; use `gbparse` for all workflows.
```

Update the testing section:

```diff
-  - `tests/test_cli.py`: Unified `gbparse` CLI subcommands and legacy script wrappers.
+  - `tests/test_cli.py`: Unified `gbparse` CLI subcommands.
```

---

### 5. Update README.md

#### [MODIFY] [README.md](file:///d:/W/Skills%20Claude/genbank-feature-parser/README.md)

No `scripts/` references exist in README.md currently — it already documents only `gbparse`. No changes required.

---

### 6. Bump Version

#### [MODIFY] [pyproject.toml](file:///d:/W/Skills%20Claude/genbank-feature-parser/pyproject.toml)

```diff
-version = "0.2.1"
+version = "0.3.0"
```

---

## Open Questions

> [!IMPORTANT]
> The `scripts/` directory currently contains no developer/build/maintenance utilities — every file is a user-facing wrapper. After deletion, the `scripts/` directory itself will not exist. If you want to retain a `scripts/` directory for future dev tooling (e.g., a release script, fixture generator), we can create a placeholder. Otherwise the directory is simply removed.

> [!NOTE]
> The `model.py` method `to_dict()` has a docstring reading `"Convert to legacy feature dictionary format."` — this is unrelated to the script wrappers (it's the dict-style API compatibility layer for `GenBankFeature`). The `test_legacy_dict_preserves_unknown_strand_states` test in `test_patch2.py` tests this dict API, not the scripts. Both should be retained as-is.

## Verification Plan

### Automated Tests

```bash
# Confirm all remaining tests pass after script removal
pytest -v

# Confirm gbparse entry point still works
gbparse validate tests/fixtures/simple_cds.gb
gbparse --help
```

### Manual Verification

- Confirm `scripts/` directory no longer exists in the working tree.
- Confirm `git status` shows clean deletions staged together (per atomic restructuring rule).
- Grep the entire repo for stale `scripts/` or `legacy` references that should have been updated.
