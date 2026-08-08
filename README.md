# GenBank Feature Parser & Annotation Engine

A Biopython-powered genome-annotation query engine, validation suite, and CLI toolset (`gbparse`) for parsing, validating, analyzing, and mining GenBank flatfiles (`.gbff`, `.gbk`, `.gb`, `.txt`).

[![CI](https://github.com/WhyAdr/genbank-parser/actions/workflows/ci.yml/badge.svg)](https://github.com/WhyAdr/genbank-parser/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

---

## Key Capabilities & Biological Semantics

- **Standard-Compliant Core**: Built on **Biopython (`Bio.SeqIO`)** with a typed dataclass model (`GenBankDocument`, `GenBankRecord`, `GenBankFeature`).
- **Faithful Biological Locations**: Preserves `FeatureLocation` and `CompoundLocation` (`join`/`order`), calculating biological feature length (`len(location)`) and extracting biological coding sequences via `SeqFeature.extract()`.
- **GFF3 Export**: Correctly computes CDS translation phase ($0, 1, 2$) across multi-exon/joined segments in 5' $\rightarrow$ 3' transcription order, emits complete `##sequence-region` extents, and handles ordinary features without artificial parent splits.
- **QC & Semantic Validation**: Verifies translation integrity against genetic codes (`transl_table`) and `codon_start` offsets, with structured severity findings (`ERROR`, `WARNING`, `INFO`) and pseudogene tolerance.
- **18 Analysis Tools & Unified CLI**: Provides `gbparse` with subcommands covering feature search, sub-region extraction with coordinate rebasing (for clinker/antiSMASH), annotation diffing, codon bias (RSCU), phylogenomic markers, and mobilome discovery.

---

## Installation & Setup

Requires Python 3.8+ and `biopython>=1.80`.

```bash
# Clone the repository
git clone https://github.com/WhyAdr/genbank-parser.git
cd genbank-parser

# Install in editable mode
pip install -e .

# Or install with test dependencies
pip install -e .[test]
```

---

## Unified CLI (`gbparse`) Usage

```bash
# 1. Structural & biological validation
gbparse validate input.gbff [--json]

# 2. Record metadata and contig summary
gbparse summary input.gbff

# 3. Export tab-delimited annotation TSV
gbparse extract input.gbff annotations.tsv

# 4. Search features by gene, product, KO, EC, Pfam, or regex
gbparse search input.gbff --gene ladA --format tsv

# 5. Single-locus qualifier deep-dive
gbparse locus input.gbff LOCUS_TAG

# 6. View genomic neighborhood window (+/- N genes, circular-aware)
gbparse neighborhood input.gbff LOCUS_TAG 5

# 7. Extract genomic sub-region with coordinate rebasing
gbparse region input.gbff --locus LOCUS_TAG --flank-genes 5 --rebase --output region.gbk

# 8. Export protein FASTA from translations
gbparse fasta input.gbff proteins.faa

# 9. Extract genome FASTA (.fna) and CDS nucleotide slices (.ffn)
gbparse sequence input.gbff --fna genome.fna --ffn cds.ffn

# 10. Calculate codon usage, RSCU, and positional GC
gbparse codon input.gbff --output codon_usage.tsv

# 11. Functional profiling (COG distribution & pathway completeness)
gbparse functional input.gbff --format tsv

# 12. Mine mobilome islands, operons, and dark-matter clusters
gbparse discover input.gbff --cluster-gap 5000 --format text

# 13. Multi-genome comparative presence/absence matrix
gbparse compare ./genomes/ --targets "ladA,ssuD,K20938" --output matrix.tsv

# 14. Compare two annotation versions of the same genome
gbparse diff old.gbff new.gbff --format text

# 15. Extract phylogenomic core/housekeeping marker genes
gbparse phylo input.gbff --markers all --output-dir ./markers/

# 16. Detect CRISPR repeat arrays and Cas gene clusters
gbparse crispr input.gbff --window 15000

# 17. Convert GenBank to standard GFF3
gbparse gff input.gbff output.gff3 --include-fasta

# 18. Parse Bakta multi-isolate summary tables
gbparse batch-summary ./isolates/ --csv summary.csv
```

---

## Python API Usage

```python
from genbank_parser import read_genbank, extract_xrefs

# Load full document
doc = read_genbank("input.gbff")

# Iterate contigs/records
for rec in doc.records:
    print(f"Contig: {rec.id} ({rec.length:,} bp, topology: {rec.topology}, GC: {rec.gc_content:.1f}%)")
    
# Find a specific locus and extract biological sequence
match = doc.find_locus("ABC_00123")
if match:
    rec, feat = match
    print(feat.gene, feat.product, feat.length)
    nt_seq = feat.extract(rec.seq)
    xrefs = extract_xrefs(feat)
    print("KEGG KOs:", xrefs["kegg_kos"])
```

---

## Testing & Verification

Run the full pytest test suite:

```bash
pytest -v
```

---

## License

MIT License
