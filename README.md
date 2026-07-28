# GenBank Feature Parser Suite

A Biopython-powered Python toolset and AI agent skill for parsing, validating, analyzing, and mining GenBank annotation files (`.gbff`, `.gbk`, `.gb`, `.txt`).

## Architecture & Key Principles

- **Standard-Compliant Parsing**: Built on **Biopython (`Bio.SeqIO`)** for robust flat-file parsing, multi-line qualifier continuation, quote stripping, and partial coordinate handling (`BeforePosition`/`AfterPosition`).
- **Contig-Boundary Aware**: All spatial operations (operon prediction, neighborhood windowing, mobilome island discovery) strictly partition features by contig before coordinate sorting.
- **Rich Semantic Cross-References**: Automatically extracts GO terms, COG IDs, KEGG KOs, Pfam, Rfam, and EC numbers from both INSDC `/EC_number` and Bakta/Prokka `/db_xref` and `/note` qualifiers.
- **16 Tool Modules**: Comprehensive CLI coverage spanning sequence extraction, locus indexing, operon prediction, CRISPR array detection, GFF3 conversion, codon usage bias (RSCU), and phylogenomic marker extraction.

---

## Installation & Setup

Requires Python 3.8+ and `biopython>=1.80`.

```bash
pip install -r requirements.txt
```

---

## Directory Layout

```text
genbank-feature-parser/
├── SKILL.md                        # Agent skill prompt and instructions
├── README.md                       # Documentation & usage guide
├── requirements.txt                # Python dependencies (biopython)
├── .gitignore                      # Git exclusion rules
├── pf_isolate_discovery_analysis.md # Example isolate discovery report
├── pf_supplementary_analysis.md    # Example genomic analysis report
└── scripts/
    ├── genbank_parser.py           # Core canonical parser library (Bio.SeqIO backend)
    ├── genbank_meor.py             # MEOR & biosurfactant biosynthesis discovery engine
    ├── genbank_discover.py         # Genomic island & operon discovery engine
    ├── genbank_functional.py       # Functional annotation (COG/GO/KEGG/EC) tables
    ├── genbank_gff.py              # GenBank to GFF3 converter & validator
    ├── genbank_codon.py            # Codon usage bias & RSCU calculator
    ├── genbank_crispr.py           # CRISPR array & repeat detector
    ├── genbank_phylo.py            # Phylogenomic marker extraction
    ├── genbank_compare.py          # Comparative synteny & ortholog analyzer
    ├── genbank_operons.py          # Intergenic distance operon predictor
    ├── genbank_neighborhood.py     # Genomic locus neighborhood window extractor
    ├── genbank_extract.py          # Feature & sequence extraction tool
    ├── genbank_fasta.py            # FASTA sequence exporter (CDS/gene/genome)
    ├── genbank_locus.py            # Locus tag indexer & feature search
    ├── genbank_metadata.py         # Record metadata & summary stats
    ├── genbank_sequence.py         # Nucleotide & translation fetcher
    └── genbank_validate.py         # Structural syntax & format validator
```

---

## Quick Start & Tool Usage Examples

### 1. Structural Validation
```bash
python scripts/genbank_validate.py input.gbff
```

### 2. Extract Functional Profiling (GO, COG, KEGG, Pfam, EC)
```bash
python scripts/genbank_functional.py input.gbff --format tsv
```

### 3. Convert GenBank to GFF3
```bash
python scripts/genbank_gff.py input.gbff output.gff3
```

### 4. Mine Dark-Matter Operons & Mobilome Islands
```bash
python scripts/genbank_discover.py input.gbff --format json
```

### 5. Calculate Codon Usage & RSCU
```bash
python scripts/genbank_codon.py input.gbff
```

### 6. Detect CRISPR Arrays & Cas Genes
```bash
python scripts/genbank_crispr.py input.gbff
```

### 7. Export Protein FASTA
```bash
python scripts/genbank_fasta.py input.gbff output_proteins.faa
```

---

## License

MIT License
