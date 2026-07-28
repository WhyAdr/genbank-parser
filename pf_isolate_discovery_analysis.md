# PF Isolate — Genome Discovery Analysis

**Genome**: *Bacillus pacificus* PF_NNT (reoriented assembly)
**Pipeline**: `genbank_discover.py` v3 (word-boundary regex, exclude-lists, interval distance)

---

## Genome Overview

| Metric | Value |
|---|---|
| Total CDS | 5,615 |
| Hypothetical / unknown | 552 (9.8%) |
| Contigs | 5 |
| Largest contig | contig_1: 5,205 kb (5,246 CDS) |

> [!NOTE]
> The hypothetical fraction dropped from **25% to 9.8%** after fixing `is_hypothetical()`. The previous version was misclassifying well-characterised domain-containing proteins (e.g., "HTH arsR-type domain-containing protein") as hypothetical, grossly inflating the dark matter count.

### Contig Summary

| Contig | Span | CDS | Hits | Density |
|---|---|---|---|---|
| contig_1 | 5,205.0 kb | 5,246 | 239 | 4.6% |
| contig_2 | 246.7 kb | 245 | 30 | **12.2%** |
| contig_3 | 81.7 kb | 93 | 10 | **10.8%** |
| contig_4 | 13.1 kb | 27 | 1 | 3.7% |
| contig_5 | 3.2 kb | 4 | 0 | 0.0% |

> [!IMPORTANT]
> **Contig 2** (246.7 kb) has 3× higher hit density than the chromosome. This is consistent with a large plasmid or megaplasmid carrying a disproportionate load of mobile elements, resistance genes, and HGT markers. **Contig 3** (81.7 kb) similarly elevated — likely a second replicon or large integrated genomic island.

### Karyogram

```
contig_1  |.  ++.#.+. ..+. ..+++# +++ ..+.+.. .++.+.. .#.+++++...+.+#.+..+.+ .....+. +.. +#++.. +..+...+..+..+  +...|
contig_2  |###.#|
contig_3  |#.|
contig_4  |#|
contig_5  | |
```

The karyogram reveals that hits on contig_1 are **not uniformly distributed** — there are clear hotspots around windows 8–9, 22–24, 41–47, 57–58, and 83–84 (each window = 50 kb). These correspond to genomic island-like clusters where mobile elements, resistance genes, and secondary metabolite loci co-localize.

---

## Category Breakdown

### 1. Heavy Metal & Stress Resistance — 52 hits (43 high-confidence)

This is the most biologically significant category for the PF isolate. The 52 hits (down from the artifactual 165 in the old script) represent genuine metal homeostasis and stress response genes.

#### Key findings:

**Arsenic resistance is dominant.** The genome carries a dense arsenic resistance arsenal:
- **14 ArsR-family regulators** scattered across the genome — an unusually high copy number. ArsR is the canonical arsenic-responsive transcriptional repressor. Multiple copies suggest either (a) recent gene duplication/HGT amplification under arsenic selection pressure, or (b) co-option for regulating other stress pathways (ArsR/SmtB family members also respond to cadmium, bismuth, and antimony).
- **3 ArsB membrane pumps** (PFNNT_016525, PFNNT_024155, PFNNT_025685) — the arsenical efflux system.
- **1 ArsC arsenate reductase** (PFNNT_024540) — reduces arsenate to arsenite for efflux.

**Mercury regulators without a complete mer operon.** There are **12 MerR-family transcriptional regulators** but no `merA` (mercuric reductase), `merB` (organomercurial lyase), or `merT/merC` (transport). This is significant: MerR-type HTH domains are one of the most versatile regulatory folds in bacteria and frequently regulate genes *other* than mercury resistance — including oxidative stress (SoxR), multidrug efflux (BmrR), and copper homeostasis. The PF isolate likely uses these as general stress sensors, not for mercury detoxification per se.

**Other metals:**
- **CzcD** (PFNNT_008680) — cobalt/zinc/cadmium efflux transporter
- **2 cadmium resistance genes** — CadA-type P-type ATPase (PFNNT_002475) and dedicated transporter (PFNNT_013070)
- **CopC** copper resistance (PFNNT_009495)
- **Tellurite resistance** (PFNNT_002445)

**Stress response:**
- SOD battery: 3 Mn-SODs + 1 Cu/Zn-SOD — excellent oxidative stress defense
- 2 catalases (KatB, KatE)
- 2 cold-shock proteins, 1 universal stress protein

#### Biological interpretation:

The PF isolate has a **generalist heavy metal and oxidative stress resistance profile** dominated by arsenic efflux. The high ArsR copy number and scattered distribution suggest this organism has experienced sustained selective pressure from arsenic-contaminated environments — consistent with soil or rhizosphere habitats where arsenic is naturally present. The lack of a complete *mer* operon indicates mercury resistance is not a primary adaptive trait.

---

### 2. Xenobiotic Degradation & Bioremediation — 27 hits (7 high-confidence)

#### Key findings:

**Ring-cleaving dioxygenases** — The MhqA cluster (PFNNT_005365/005370) on contig_1 is a **duplicated pair** of ring-cleaving dioxygenases. MhqA is specifically involved in methylhydroquinone detoxification — an oxidative stress response to plant-derived quinones. Both genes are co-directional (forward strand), suggesting a recent tandem duplication.

**Catechol dioxygenases** — Two copies of CatE (PFNNT_016780, PFNNT_021300) for meta-cleavage of catechol. This is a key step in aromatic compound degradation and could indicate capacity for lignin degradation intermediates.

**Haloacid dehalogenases** — Three hits (PFNNT_004580, PFNNT_015270 *menH*, PFNNT_020205 *ycjU*). PFNNT_004580 is the strongest candidate for genuine xenobiotic dehalogenation; the other two may be HAD-family hydrolases with broader substrate specificity.

**Alkanesulfonate monooxygenase** (PFNNT_014190, *ssuD*) — FMNH2-dependent enzyme for sulfonate catabolism. Important for sulfur scavenging in sulfur-limited environments.

**Anthranilate/hydroxyanthranilate degradation** (Cluster 3: PFNNT_010215 + PFNNT_010230) — This cluster suggests capacity for tryptophan-derived aromatic catabolism.

#### Biological interpretation:

The xenobiotic profile is characteristic of a **soil-dwelling saprophyte** — ring cleavage, catechol degradation, and sulfonate utilization are hallmarks of organisms that process plant-derived aromatics and compete for alternative sulfur sources. This is not a dedicated hydrocarbon degrader (no alkane hydroxylases, no *alkB*, no *ladA*), but it can process the aromatic intermediates that arise from lignin and humus decomposition.

---

### 3. Secondary Metabolites & BGCs — 14 hits (8 high-confidence)

#### Key findings:

**Bacillibactin siderophore cluster (Cluster 1)** — The dhbA/dhbE/dhbF cluster (PFNNT_011405–011425, spanning 11.7 kb) is a **canonical bacillibactin biosynthesis operon**:
- *dhbA*: 2,3-dihydroxybenzoate dehydrogenase
- *dhbE*: 2,3-dihydroxybenzoate adenylase
- *dhbF*: Non-ribosomal peptide synthetase

This is the primary iron acquisition system in *Bacillus* species. Bacillibactin is a catecholate-type siderophore and one of the strongest iron chelators known (stability constant log K > 48). Its presence confirms the isolate can compete for iron in iron-limited environments.

> [!NOTE]
> In the old script, `dhbA` and `dhbE` would have been falsely classified under "Xenobiotic Degradation" because their product annotations contain "dihydroxybenzoate" which matched the old `("benzoate", 2)` keyword. The v3 exclude-list correctly routes them only to Secondary Metabolites.

**Phenazine biosynthesis** — Two PhzF-family proteins (PFNNT_011490, PFNNT_014960). Phenazines are redox-active antibiotics — potent electron shuttles that generate reactive oxygen species in competing organisms. The presence of *two separate* phzF loci suggests either an intact phenazine BGC or remnants of one.

**Polyketide synthases** — PpsC (PFNNT_005305, phthiocerol synthesis PKS type I) and CurC (PFNNT_014460). PpsC is particularly interesting as phthiocerol-related PKS systems in Bacillus may produce lipopeptide antibiotics.

**Bacteriocins** — Two loci: a thiazole-containing bacteriocin maturation protein (PFNNT_006485) and a standalone bacteriocin (PFNNT_028230, on contig_2). The contig_2 bacteriocin is notable because it's within 5 kb of a transposase (PFNNT_028185), suggesting possible mobilization.

**Enterobactin esterase** (PFNNT_006350) — This doesn't synthesize enterobactin; it *degrades* it. This is a **siderophore piracy** gene — the isolate can steal iron from enterobactin-producing competitors by hydrolyzing their siderophore and releasing the chelated iron.

#### Biological interpretation:

The PF isolate has a **competitive secondary metabolite arsenal**: bacillibactin for iron acquisition, phenazines for antibiosis, PKS-derived lipopeptides, and enterobactin piracy. This is consistent with a rhizosphere or soil ecology where iron competition and antibiotic warfare are intense.

---

### 4. Mobilome & HGT Markers — 187 hits (96 high-confidence)

#### Key findings:

**IS element diversity** — The genome carries transposases from at least 6 IS families:
- **IS3 family** (dominant): ISBce19, ISBce15, ISBt2 — 15+ copies
- **IS1182 family**: ISBth7 — 5 copies
- **IS4 family**: 3 copies
- **IS200/IS605 family**: ISBce3, plus TnpB accessory transposase — high-priority as TnpB is the evolutionary ancestor of CRISPR-Cas12
- **IS110 family**: ISBth13 — 2 copies
- **IS256 family**: 1 copy
- **ISL3 family**: 2 copies

**Prophage remnant** — Cluster 30 (PFNNT_020445–020545) is a **7-gene phage cluster** spanning 17.7 kb on the minus strand. It contains tail proteins, head closure, portal, terminase, and an integrase — this is a largely intact **defective prophage** with structural genes but likely missing replication/lysis modules. The adjacent dark matter cluster (Cluster 8) with HNH nuclease and phage protein genes may be the remnant lysis cassette.

**Toxin-antitoxin systems** — Three distinct TA systems:
- **HicAB** (contig_2, PFNNT_027190/027195) — Type II TA, mRNA interferase
- **MntA/HepT** (contig_2, PFNNT_028110/028115) — Type VII TA
- **NdoA** (contig_1, PFNNT_001440) — Type II endoribonuclease

**Type IV secretion** (PFNNT_027655, contig_2) — A T4SS pilin on contig_2 suggests this replicon has (or had) conjugative transfer capability, reinforcing the interpretation that contig_2 is a conjugative plasmid.

**Restriction-modification** — Two restriction endonuclease pairs (PFNNT_001000/001010, PFNNT_001605) — defense against phage DNA.

#### Biological interpretation:

The mobilome is **highly active and diverse**. The IS3 family dominance with 15+ copies indicates ongoing genome plasticity. The contig_2 plasmid carries a particularly dense mobilome payload (TA systems, T4SS, transposases, recombinases) — this is a **conjugative resistance plasmid** that appears to be a hotspot for gene acquisition.

---

## Dark Matter Operons

### Pure Dark Matter

**Cluster 1 (contig_1, 243,842..244,704)** — A tandem array of **6 DUF3948 domain-containing proteins** on the forward strand. This is remarkable: DUF3948 is an uncharacterized domain found primarily in Firmicutes. A 6-copy tandem array suggests either (a) a specialized structural/surface protein with repeat domains, or (b) a recently expanded gene family under positive selection. **High priority for structural prediction (AlphaFold) and expression analysis.**

### Mixed Dark Matter — Notable Clusters

| Cluster | Location | Size | Anchor Gene(s) | Notes |
|---|---|---|---|---|
| 6 & 7 | ~3.18 Mb / ~3.43 Mb | 3 genes each | PadR regulator | **Duplicated pair** — DUF4097 + DUF1700 + PadR. Identical structure at two loci suggests recent duplication. PadR regulates phenolic acid responses. |
| 8 | ~3.97 Mb | 5 genes | HNH nuclease, Phage protein | Directly adjacent to the prophage (Cluster 30 in Mobilome). This is likely the **phage lysis/modification cassette**. |
| 13 | ~4.52 Mb | 5 genes | Phage shock protein A (PspA) | PspA responds to membrane stress from phage infection. The 3 DUFs may be phage-defense or membrane integrity genes. |
| 15 | contig_4, 7.5–13.1 kb | **17 genes** (71% hyp) | Zona occludens toxin (Zot) | **Highest priority.** A 17-gene cluster with Zot domain suggests a **phage-like element or pathogenicity island**. Zot was originally identified in CTXφ of *V. cholerae* and disrupts intestinal tight junctions. In *Bacillus*, Zot-like domains may have alternative functions (phage assembly or secretion). |

> [!WARNING]
> **Cluster 15 on contig_4** deserves immediate follow-up. A 17-gene mostly-hypothetical cluster carrying a Zona occludens toxin domain on the smallest contig (13.1 kb) is a strong candidate for a mobile pathogenicity island or cryptic prophage. This is the kind of locus that antiSMASH/PHASTER would flag.

---

## Mobilome-Adjacent Functional Islands

The most biologically interesting co-localizations:

### High-Priority Islands

| Mobilome Marker | Adjacent Functional Gene | Distance | Interpretation |
|---|---|---|---|
| PFNNT_001325/001330 (IS3 transposases) | hppD + hmgA (dioxygenases) | +3–5 kb | Aromatic catabolism genes flanked by IS3 elements — **possible acquired catabolic island** |
| PFNNT_009470/009475 (IS3 transposases) | CopC copper resistance | +4 kb | Copper resistance gene adjacent to IS3 — candidate for **mobilized copper resistance** |
| PFNNT_014510/014515 (transposases) | ArsR + HmoA + PKS CurC | 1–8 kb | Three functional genes near transposases — a **resistance/metabolism island** |
| PFNNT_025650/025655/025665 (IS3 + IS200/IS605 TnpB) | ArsB arsenical pump | +3–7 kb | Arsenical efflux pump flanked by two IS families — **classic mobilized resistance determinant** |
| PFNNT_028185 (transposase, contig_2) | Bacteriocin | +5 kb | Bacteriocin gene on the putative conjugative plasmid next to a transposase — **mobilized bacteriocin** |
| PFNNT_027380–027420 (transposase cluster, contig_2) | ArsR regulator (PFNNT_027425) | 0–9 kb | 5 different mobilome markers all adjacent to an ArsR — **genomic island carrying arsenic regulation on the plasmid** |

---

## Summary & Recommendations

### Ecological Profile

The PF isolate presents as a **rhizosphere-adapted, metal-tolerant soil bacterium** with:
- Dominant arsenic resistance (14 ArsR + 3 ArsB + ArsC)
- Competitive iron acquisition (bacillibactin + enterobactin piracy)
- Aromatic catabolism capability (ring cleavage, catechol degradation)
- Antibiotic production potential (phenazines, PKS, bacteriocins)
- A conjugative plasmid (contig_2) carrying mobilized resistance and bacteriocin genes

### Recommended Follow-Up

1. **Contig_4 Cluster 15** (17-gene Zot island) — Submit to PHASTER for prophage classification; run AlphaFold on the 12 hypothetical proteins
2. **DUF3948 6-copy array** — Structural prediction to determine if this is a surface protein or novel function
3. **PadR duplicated clusters** (6 & 7) — Expression analysis under phenolic acid stress to determine if both copies are functional
4. **Contig_2 plasmid** — Run MOB-suite or PlasmidFinder to classify the replicon type and confirm conjugative transfer capability
5. **antiSMASH** — Submit for comprehensive BGC detection, specifically to validate the phenazine and PKS clusters
6. **ArsR copy number** — Comparative genomics against other *B. pacificus* strains to determine if the 14-copy expansion is species-typical or strain-specific
