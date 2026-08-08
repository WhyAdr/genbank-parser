---
name: genbank-feature-parser
description: Parse GenBank flatfiles computationally. Validates structure, extracts biological sequences, searches annotations, and runs diagnostic bioinformatic analyses without loading entire files into the AI context window. Use this skill whenever the user asks to parse, validate, summarize, compare, diff, or analyze .gbff, .gbk, .gb, or .txt GenBank files; when they want gene neighborhoods, operon candidates, sub-regions with coordinate rebasing, COG/GO/KEGG/EC tables, codon usage (RSCU), phylogenomic markers, CRISPR detection, or GFF3 export; or when they mention tools like Bakta, Prokka, tbl2asn, or NCBI GenBank.
---

# GenBank Feature Parser & Annotation Engine

## Purpose

Parse, validate, and analyze GenBank flatfiles (`.gb`, `.gbk`, `.gbff`, `.txt`) **computationally** — never loading the full file content into the AI context window. Powered by Biopython `Bio.SeqIO` and structured into an installable Python package (`genbank_parser`) with a unified CLI (`gbparse`) and backwards-compatible scripts under `scripts/`.

## Epistemic Discipline

1. **Scripts parse; the AI reads summaries.** Never `view_file` on large GenBank files — invoke the CLI or script and read its compact summary.
2. **One canonical parser, reused everywhere.** All commands and scripts import from `genbank_parser` (powered by `read_genbank` and `GenBankDocument`). Never reconstruct custom parsers.
3. **Compact structured output only.** Tools emit diagnostic summaries, tabular TSVs, or schema-versioned JSON.

---

## Unified CLI (`gbparse`) & Script Inventory

You can invoke commands either via the unified `gbparse` CLI or through the legacy scripts under `scripts/`:

| `gbparse` Subcommand | Legacy Script | Purpose | Typical Invocation |
|---|---|---|---|
| `gbparse validate` | `genbank_validate.py` | Biological translation QC & structure report | `gbparse validate INPUT.gbff [--json]` |
| `gbparse summary` | `genbank_metadata.py` | LOCUS metadata, contigs, topologies & length | `gbparse summary INPUT.gbff` |
| `gbparse extract` | `genbank_extract.py` | Tab-delimited annotation TSV export | `gbparse extract INPUT.gbff [output.tsv]` |
| `gbparse search` | (new) | Search features by gene, product, KO, EC, Pfam | `gbparse search INPUT.gbff --gene ladA --format tsv` |
| `gbparse locus` | `genbank_locus.py` | Single-locus qualifier deep-dive | `gbparse locus INPUT.gbff LOCUS_TAG` |
| `gbparse neighborhood` | `genbank_neighborhood.py` | Circular-aware flanking gene viewer (+/- N) | `gbparse neighborhood INPUT.gbff LOCUS_TAG [window]` |
| `gbparse region` | (new) | Sub-region extraction with coordinate rebasing | `gbparse region INPUT.gbff --locus TAG --flank-genes 5 --rebase --output region.gbk` |
| `gbparse fasta` | `genbank_fasta.py` | Export all CDS translations as protein FASTA | `gbparse fasta INPUT.gbff [proteins.faa]` |
| `gbparse sequence` | `genbank_sequence.py` | Extract genome FASTA (.fna) & CDS (.ffn) | `gbparse sequence INPUT.gbff [--fna out.fna] [--ffn out.ffn]` |
| `gbparse codon` | `genbank_codon.py` | Codon usage bias, RSCU & positional GC | `gbparse codon INPUT.gbff [--min-len 100]` |
| `gbparse functional` | `genbank_functional.py` | COG distribution + metabolic completeness | `gbparse functional INPUT.gbff [--format json]` |
| `gbparse discover` | `genbank_discover.py` | Mine mobilome islands & dark-matter clusters | `gbparse discover INPUT.gbff [--cluster-gap 5000] [--rules rules.yaml]` |
| `gbparse compare` | `genbank_compare.py` | Multi-genome marker presence/absence matrix | `gbparse compare genomes/ --targets "ladA,ssuD,K20938"` |
| `gbparse diff` | (new) | Compare two annotation versions of a genome | `gbparse diff old.gbff new.gbff [--format json]` |
| `gbparse phylo` | `genbank_phylo.py` | Phylogenomic marker extraction (core/housekeeping) | `gbparse phylo INPUT.gbff [--markers all] [--output-dir ./markers/]` |
| `gbparse crispr` | `genbank_crispr.py` | CRISPR array + Cas gene spatial locator | `gbparse crispr INPUT.gbff [--window 15000]` |
| `gbparse gff` | `genbank_gff.py` | Standard GFF3 export (with CDS phase & regions) | `gbparse gff INPUT.gbff [output.gff3] [--include-fasta]` |
| `gbparse batch-summary` | `parse_bakta_summaries.py`| Bakta multi-isolate comparison tables | `gbparse batch-summary ./isolates/ --csv summary.csv` |

---

## Python API Usage

```python
from genbank_parser import read_genbank, extract_xrefs

# Load full typed document
doc = read_genbank("genome.gbff")

# Iterate contigs/records
for rec in doc.records:
    print(rec.id, rec.length, rec.topology, rec.gc_content)
    
# Find feature and extract biological sequence
match = doc.find_locus("ABC_00123")
if match:
    rec, feat = match
    print(feat.gene, feat.product, feat.length)
    nt_seq = feat.extract(rec.seq)
    xrefs = extract_xrefs(feat)
    print("KEGG KOs:", xrefs['kegg_kos'])
```

---

## Biological Semantics Handled

- **True Biological Extraction**: Slices joined/spliced features via `SeqFeature.extract()`, ignoring non-coding gaps.
- **Compound Locations**: Preserves `join()` and `order()` sub-segments, phases in GFF3, and biological lengths.
- **Strand Semantics**: Preserves $+1$, $-1$, $0$, and `None` states (rendered as `+`, `-`, `?`, `.`).
- **Circular Topology**: Handles circular contigs/plasmids and origin-spanning neighborhoods.
- **Semantic Cross-References**: Maps INSDC and Bakta `/db_xref`, `/note`, and `/EC_number` to typed identifiers (`go_terms`, `cog_ids`, `kegg_kos`, `pfam`, `rfam`, `ec_numbers`).
