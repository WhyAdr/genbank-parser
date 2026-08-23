---
name: genbank-feature-parser
description: Parse GenBank flatfiles computationally. Validates structure, extracts biological sequences, searches annotations, and runs diagnostic bioinformatic analyses without loading entire files into the AI context window. Use this skill whenever the user asks to parse, validate, summarize, compare, diff, or analyze .gbff, .gbk, .gb, or .txt GenBank files; when they want gene neighborhoods, operon candidates, local-coordinate sub-regions, COG/GO/KEGG/EC tables, codon usage (RSCU), candidate phylogenetic markers, CRISPR/Cas annotation scanning, GFF3 export, MEOR or petroleum-microbiology potential, hydrocarbon degradation, biosurfactants or bio-emulsifiers, alkB/ladA/almA, rhamnolipids, BTEX/PAH degradation, or anaerobic fumarate-addition markers; or when they mention tools like Bakta, Prokka, tbl2asn, or NCBI GenBank.
---

# GenBank Feature Parser & Annotation Engine

## Purpose

Parse, validate, and analyze GenBank flatfiles (`.gb`, `.gbk`, `.gbff`, `.txt`) **computationally** — never loading the full file content into the AI context window. Powered by Biopython `Bio.SeqIO` and structured into an installable Python package (`genbank_parser`) with a unified CLI (`gbparse`).

## Epistemic Discipline

1. **The CLI parses; the AI reads summaries.** Never `view_file` on large GenBank files — invoke the `gbparse` CLI and read its compact summary.
2. **One canonical parser, reused everywhere.** All commands import from `genbank_parser` (powered by `read_genbank` and `GenBankDocument`). Never reconstruct custom parsers.
3. **Compact structured output only.** Tools emit diagnostic summaries, tabular TSVs, or schema-versioned JSON.
4. **MEOR genomic potential is not phenotype.** `gbparse meor` detects annotation-supported candidates. Do not present hits as proof of expression, enzyme activity, hydrocarbon turnover, biosurfactant production, or field-scale oil recovery.
5. **Mobilome evidence is not transfer or identity.** `gbparse mobilome` retains annotation-supported record, feature, and component evidence; do not turn it into claims of plasmid autonomy, mechanism, transfer, co-transfer, expression, phenotype, AMR, virulence, ancestry, or phagemid identity.

---

## CLI Reference (`gbparse`)

| Subcommand | Purpose | Typical Invocation |
|---|---|---|
| `validate` | Biological translation QC & structure report | `gbparse validate INPUT.gbff [--json]` |
| `summary` | LOCUS metadata, contigs, topologies & length | `gbparse summary INPUT.gbff` |
| `extract` | Tab-delimited annotation TSV export | `gbparse extract INPUT.gbff [output.tsv]` |
| `search` | Search features by gene, product, KO, EC, Pfam | `gbparse search INPUT.gbff --gene ladA --format tsv` |
| `locus` | Single-locus qualifier deep-dive | `gbparse locus INPUT.gbff LOCUS_TAG` |
| `neighborhood` | Structured circular-aware flanking gene data and optional static figure | `gbparse neighborhood INPUT.gbff LOCUS_TAG [window] [--format json] [--visualize]` |
| `region` | Sub-region extraction with valid local coordinates | `gbparse region INPUT.gbff --locus TAG --flank-genes 5 --output region.gbk` |
| `fasta` | Export all CDS translations as protein FASTA | `gbparse fasta INPUT.gbff [proteins.faa]` |
| `sequence` | Extract genome FASTA (.fna) & CDS (.ffn) | `gbparse sequence INPUT.gbff [--fna out.fna] [--ffn out.ffn]` |
| `codon` | Codon usage bias, RSCU & positional GC | `gbparse codon INPUT.gbff [--min-len 100]` |
| `functional` | COG distribution + metabolic completeness | `gbparse functional INPUT.gbff [--format json]` |
| `discover` | Scan annotation-supported mobilome/xenobiotic islands | `gbparse discover INPUT.gbff [--ruleset mobilome] [--format tsv]` |
| `compare` | Multi-genome marker presence/absence matrix | `gbparse compare genomes/ --targets "ladA,ssuD,K20938"` |
| `diff` | Compare two annotation versions of a genome | `gbparse diff old.gbff new.gbff [--format json]` |
| `phylo` | Annotation-based candidate phylogenetic markers | `gbparse phylo INPUT.gbff [--markers all] [--min-length 50]` |
| `crispr` | CRISPR/Cas annotation scanner | `gbparse crispr INPUT.gbff [--window 15000]` |
| `gff` | Standard GFF3 export (with CDS phase & regions) | `gbparse gff INPUT.gbff [output.gff3] [--include-fasta]` |
| `batch-summary` | Bakta multi-isolate comparison tables | `gbparse batch-summary ./isolates/ --csv summary.csv` |
| `meor` | MEOR, hydrocarbon-degradation, biosurfactant and bio-emulsifier genomic-potential analysis | `gbparse meor INPUT.gbff --min-weight 2 --format json` |
| `mobilome` | Replicon-centric, annotation-supported mobilome evidence report | `gbparse mobilome INPUT.gbff --format json --min-evidence 2` |

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

Mobilome analysis deliberately uses a separate subpackage:

```python
from genbank_parser.mobilome import analyze_mobilome

report = analyze_mobilome("genome.gbff", include="all", min_evidence=2)
```

---

## Biological Semantics Handled

- **True Biological Extraction**: Slices joined/spliced features via `SeqFeature.extract()`, ignoring non-coding gaps.
- **Compound Locations**: Preserves `join()` and `order()` sub-segments, phases in GFF3, and biological lengths.
- **Strand Semantics**: Preserves $+1$, $-1$, $0$, and `None` states (rendered as `+`, `-`, `?`, `.`).
- **Circular Topology**: Handles circular contigs/plasmids and origin-spanning neighborhoods.
- **Neighborhood Context Overlays**: Optional non-CDS context, canonical mobilome/xenobiotics annotation-rule matches, and tested same-strand proximity links feed the structured result before dna_features_viewer renders it. These overlays are annotation-supported context, not evidence of horizontal transfer, co-transcription, or phenotype.
- **Semantic Cross-References**: Maps INSDC and Bakta `/db_xref`, `/note`, and `/EC_number` to typed identifiers (`go_terms`, `cog_ids`, `kegg_kos`, `pfam`, `rfam`, `ec_numbers`).
- **Source-Aware MEOR Evidence**: Active structured KO, fully specified structured EC, and gene evidence is Weight 3; product evidence is Weight 2; free-text `/note` evidence is Weight 1. Wildcard ECs and note-only identifiers are never promoted to structured evidence.
- **MEOR Scope**: Pathway completeness is genome-level annotation completeness, while clusters require distinct physical genes on the same contig and known common strand with at most N intervening bases. Neither is sequence-family confirmation or phenotype.
- **Mobilome Scope**: Record classification reflects declared source/record metadata, not topology or marker prediction. JSON and TSV retain all evidence reasons and disabled aggregate policies; only configured, cautious hypotheses use eligible non-pseudo/complete components. External typing, AMR, virulence, and boundary callers are reported as not-run handoffs with required provenance.
