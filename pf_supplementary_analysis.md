# PF Isolate — Supplementary: Bacteriocins, IAA Synthesis & Biomass Utilization

**Context**: The PF isolate was obtained from **oil palm dried-ground palm kernel waste** — a lignocellulosic substrate rich in residual palm oil (triacylglycerols), hemicellulose (mannans, xylans), cellulose, lignin-derived phenolics, and chitin from fungal colonizers. This ecological origin reframes the entire genomic interpretation.

---

## 1. Bacteriocins — Three Distinct Systems

The discovery script flagged only two generic "bacteriocin" hits. The targeted search reveals **three structurally distinct bacteriocin systems** — each with a different mode of action:

### 1a. Circularin A / Uberolysin family (PFNNT_025145)

**Annotation**: `Circularin A/uberolysin family circular bacteriocin`
**Location**: contig_1, 4,806,921–4,807,142 (forward strand)

This is a **circular bacteriocin** — a class of ribosomally synthesized antimicrobial peptides where the N- and C-termini are covalently joined into a continuous backbone, making them exceptionally resistant to proteolytic degradation. Circularin A was originally characterized in *Clostridium beijerinckii* ATCC 25752, where it forms voltage-dependent pores in target cell membranes.

**Genomic context** is telling:

| Offset | Locus | Gene | Product | Interpretation |
|---|---|---|---|---|
| -3.2 kb | PFNNT_025135 | *acrA* | Membrane fusion protein | **Immunity/export** — efflux pump to protect self |
| -3.9 kb | PFNNT_025130 | *lolD* | ABC transporter | **Secretion** — ATP-driven export of the mature bacteriocin |
| -5.1 kb | PFNNT_025125 | *salY* | Macrolide ABC transporter permease | **Secretion** — transmembrane component |
| -1.6 kb | PFNNT_025140 | — | IS110 ISBth13 transposase | **Mobilization** — IS element flanking the cluster |
| +0.3 kb | PFNNT_025150 | — | ABC transporter permease | **Dedicated export** |
| +2.0 kb | PFNNT_025155 | — | Stage II sporulation protein M | Sporulation linkage |
| +2.6 kb | PFNNT_025160 | — | ABC transporter domain protein | Additional export |
| +9.2 kb | PFNNT_025195 | *lrgB* | Holin-like protein CidB | **Programmed cell death** |
| +9.8 kb | PFNNT_025200 | *yohJ* | Holin-like protein | **Programmed cell death** |

> [!IMPORTANT]
> The Circularin A gene sits in a **complete biosynthetic neighbourhood**: the bacteriocin structural gene is flanked upstream by ABC transporters and a membrane fusion protein (self-immunity + export), and downstream by additional ABC exporters. The IS110 transposase 1.6 kb upstream suggests this entire island may have been acquired via HGT. The nearby holin-like proteins (CidB, YohJ) are particularly interesting — in *Bacillus*, CidAB/LrgAB holins regulate programmed cell death and biofilm formation. Their proximity to a circular bacteriocin suggests coordinated lysis-mediated bacteriocin release during stationary phase or sporulation.

### 1b. Sonorensin-family bacteriocins (×2 copies)

**Locus 1**: PFNNT_012715 — `Heterocycloanthracin/sonorensin family bacteriocin` (contig_1, 2,487,228–2,487,515)
**Locus 2**: PFNNT_013050 — `Heterocycloanthracin/sonorensin family bacteriocin` (contig_1, 2,553,313–2,553,558)

Sonorensin is a **heterocycloanthracin (HCA)** — a class of ribosomally synthesized and post-translationally modified peptides (RiPPs) originally characterized in *Bacillus sonorensis* MT93. HCAs contain thiazole and oxazole heterocycles formed by cyclodehydration of Cys, Ser, and Thr residues. They have broad-spectrum activity against both Gram-positive and Gram-negative bacteria, and notably against *Listeria monocytogenes*.

**Sonorensin copy 1** (PFNNT_012715) neighbourhood:

| Offset | Locus | Product | Significance |
|---|---|---|---|
| -7.5 kb | PFNNT_012680 | **Triacylglycerol lipase** | Palm oil degradation (see §3) |
| -6.8 kb | PFNNT_012685 | IS200/IS605 TnpB transposase | Mobilization — TnpB is Cas12 ancestor |
| +5.6 kb | PFNNT_012740 | DUF3959 domain protein | Dark matter cluster 3 (mixed) |
| +9.7 kb | PFNNT_012770 | **Cytochrome P450** (*cypX*) | Oxidative xenobiotic metabolism |

**Sonorensin copy 2** (PFNNT_013050) neighbourhood:

| Offset | Locus | Product | Significance |
|---|---|---|---|
| -1.0 kb | PFNNT_013040 | CutA divalent cation tolerance | Metal homeostasis |
| +3.3 kb | PFNNT_013070 | **Cadmium resistance transporter** | Heavy metal efflux |
| +9.2 kb | PFNNT_013110 | LAGLIDADG homing endonuclease | Mobile element marker |

> [!NOTE]
> **Two sonorensin copies** at different chromosomal loci (~66 kb apart) is unusual. The first copy is flanked by a TnpB transposase and a triacylglycerol lipase — suggesting it sits in a mobile element that was acquired together with lipid-degradation capability. The second copy neighbours a cadmium resistance transporter. Both copies may be independently mobilized.

### 1c. Thiazole-containing bacteriocin maturation protein (PFNNT_006485)

This is not a structural bacteriocin gene but a **maturation enzyme** for thiazole-containing bacteriocins. It likely processes a precursor peptide (possibly encoded nearby) into the active form via thiazole ring formation. The maturation machinery is characteristic of the linear azol(in)e-containing peptide (LAP) class of RiPPs.

### Bacteriocin Summary

The PF isolate carries **three mechanistically distinct antimicrobial peptide systems**:

| System | Class | Mode of Action | Self-Immunity |
|---|---|---|---|
| Circularin A | Circular bacteriocin (Class IIc) | Membrane pore formation | ABC transporters + MFP |
| Sonorensin (×2) | HCA / RiPP | Heterocycle-mediated membrane disruption | Unknown (likely nearby) |
| LAP maturation | Linear azoline peptide | Thiazole-dependent | Unknown |

This is a **formidable competitive arsenal** for a palm kernel waste colonizer — the substrate is a nutrient-rich, microbially contested environment where antibiotic warfare determines community dominance.

---

## 2. Indole-3-Pyruvate Decarboxylase — IAA Biosynthesis

**Locus**: PFNNT_011935 — `Indole-3-pyruvate decarboxylase` (contig_1, 2,317,914–2,319,590, **reverse strand**)

This enzyme catalyses the **penultimate step in the indole-3-pyruvic acid (IPyA) pathway** of auxin (indole-3-acetic acid, IAA) biosynthesis:

```
L-Tryptophan → Indole-3-pyruvate → Indole-3-acetaldehyde → IAA
                    (IpdC/IPDC)         (aldehyde oxidase)
```

The IPyA pathway is the **primary IAA biosynthesis route** in plant growth-promoting rhizobacteria (PGPR). IAA produced by soil bacteria stimulates root elongation, lateral root formation, and root hair development in plants — enhancing nutrient uptake and establishing rhizosphere colonization advantage.

**Genomic context**:

| Offset | Locus | Gene | Product | Connection |
|---|---|---|---|---|
| -1.0 kb | PFNNT_011930 | — | Methyltransferase type 11 | SAM-dependent modification |
| -2.5 kb | PFNNT_011920 | *ubiG* | SAM-dependent methyltransferase | Ubiquinone pathway |
| +1.8 kb | PFNNT_011940 | *marR* | MarR-type regulator | **Regulatory** — MarR regulators often respond to aromatic signals (salicylate, phenolics) |

**Tryptophan supply**: The genome also carries a complete tryptophan biosynthesis operon at a separate locus:
- PFNNT_006400 *trpC* — indole-3-glycerol phosphate synthase
- PFNNT_006410 *trpB* — tryptophan synthase β
- PFNNT_006415 *trpA* — tryptophan synthase α

Plus tryptophan catabolism:
- PFNNT_013375 *kynA* — tryptophan 2,3-dioxygenase (kynurenine pathway)
- PFNNT_015205 *tspO* — tryptophan-rich protein

> [!IMPORTANT]
> **The combination of IpdC + complete Trp biosynthesis + MarR regulation is the signature of a PGPR strain.** In an oil palm waste context, this means the PF isolate isn't just a saprophyte decomposing waste — it has the genetic toolkit to **promote plant growth** if applied as a biofertilizer. The MarR regulator adjacent to IpdC suggests IAA production may be **induced by plant-derived phenolic signals** (ferulic acid, coumaric acid) released during lignin decomposition, creating a positive feedback loop: waste degradation releases signals that upregulate auxin production.

---

## 3. LigB and the Lignin Degradation Pathway

**Locus**: PFNNT_009180 *ligB* — `Extradiol ring-cleavage dioxygenase class III enzyme subunit B domain-containing protein` (contig_1, 1,774,367–1,775,128, **reverse strand**)

LigB is the β-subunit of **protocatechuate 4,5-dioxygenase** — the key ring-cleavage enzyme in the **protocatechuate branch of bacterial lignin degradation**. It catalyses the extradiol (meta) cleavage of protocatechuate, which is the central aromatic intermediate produced from lignin monomers (guaiacyl, syringyl, and hydroxyphenyl units):

```
Lignin polymer
    ↓ (depolymerization by laccases/peroxidases)
Ferulic acid / Coniferyl alcohol / Vanillin
    ↓ (demethylation, oxidation)
Protocatechuate
    ↓ (LigB — extradiol ring cleavage)
4-Carboxy-2-hydroxymuconate-6-semialdehyde
    ↓ (further catabolism)
TCA cycle intermediates
```

**Supporting genes across the genome** that complete this pathway:

| Locus | Gene | Product | Pathway Role |
|---|---|---|---|
| PFNNT_009180 | *ligB* | Protocatechuate extradiol dioxygenase | **Central ring cleavage** |
| PFNNT_016780 | *catE* | Catechol 2,3-dioxygenase | Meta-cleavage of catechol branch |
| PFNNT_021300 | *catE* (copy 2) | Catechol 2,3-dioxygenase | Redundant catechol processing |
| PFNNT_005365/70 | *mhqA* (×2) | Ring-cleaving dioxygenase | Methylhydroquinone detox (quinone intermediates) |
| PFNNT_017405 | *gloA* | Ring-cleaving dioxygenase | Additional ring cleavage |
| PFNNT_016675 | — | 2,4-dichlorophenol 6-monooxygenase | **Chlorophenol/phenolic hydroxylation** |

**LigB neighbourhood** (±15 kb):

The ligB locus sits in a metabolically dense region:
- +1.0 kb: *pldB* lysophospholipase — lipid degradation
- +3.4 kb: Acetyl-CoA hydrolase/transferase — CoA-mediated aromatic catabolism
- -1.3 kb: *araJ* MFS transporter — aromatic compound uptake

> [!NOTE]
> **Palm kernel waste is rich in lignin-derived phenolics.** During composting and microbial decomposition, the lignocellulosic matrix releases ferulic acid, p-coumaric acid, vanillin, and syringaldehyde. The PF isolate's LigB provides the enzymatic capability to cleave these aromatics and funnel them into central metabolism. This isn't just detoxification — it's **carbon acquisition from lignin**, giving the isolate a competitive advantage in a lignin-rich waste environment.

---

## 4. Biomass-Degrading Hydrolase Arsenal

The initial analysis completely overlooked the **polysaccharide and lipid degradation toolkit** — which is arguably the most ecologically relevant gene set for a palm kernel waste isolate. Here's the full inventory:

### 4a. Polysaccharide Hydrolases

| Locus | Gene | Product | Substrate | PKW Relevance |
|---|---|---|---|---|
| PFNNT_005990 | *amyA* | **Neopullulanase** | α-1,4/α-1,6 glucan | Starch/pullulan in PKW |
| PFNNT_013255 | *pulA* | **Type I pullulanase** | α-1,6 glucosidic bonds | Debranching starch |
| PFNNT_023055 | *pulA* (copy 2) | **Type I pullulanase** | α-1,6 glucosidic bonds | Redundant — high-demand enzyme |
| PFNNT_022360 | *frvX* | **Cellulase** | β-1,4 glucan (cellulose) | Cellulose in PKW fibre |
| PFNNT_023575 | — | **Endoglucanase** | β-1,4 glucan (cellulose) | Internal cellulose chain cleavage |
| PFNNT_017775 | — | **Chitinase** | β-1,4 GlcNAc (chitin) | Fungal cell wall degradation |
| PFNNT_012605 | — | **L-arabinolactonase** | Arabinoxylan intermediates | Hemicellulose processing |

**Two pullulanases** (PFNNT_013255 and PFNNT_023055) is notable. Pullulanases are debranching enzymes that hydrolyse α-1,6 glycosidic bonds in starch, amylopectin, and pullulan. Palm kernel expeller retains significant residual starch, and carrying two copies suggests this is a **high-priority catabolic function** under selection.

**The cellulase** (*frvX*, PFNNT_022360) and **endoglucanase** (PFNNT_023575) together provide cellulose degradation capability. *frvX* is annotated as a GH family cellulase — it likely performs exo-type hydrolysis, while the endoglucanase provides internal chain cleavage. Together they constitute a **minimal cellulolytic system**.

**The chitinase** (PFNNT_017775) is ecologically significant in PKW — the waste substrate is rapidly colonized by filamentous fungi (*Aspergillus*, *Trichoderma*, *Rhizopus*), and chitinase enables the PF isolate to **lyse fungal competitors** by degrading their cell walls. The chitinase neighbourhood includes copper chaperones (*copZ*) and an ArsR regulator, suggesting co-regulation with metal stress.

### 4b. Lipases and Esterases

Palm kernel waste retains 5–8% residual palm kernel oil (predominantly lauric acid C12:0 and myristic acid C14:0 triacylglycerols). The isolate carries a comprehensive lipolytic arsenal:

| Locus | Gene | Product | Substrate |
|---|---|---|---|
| PFNNT_005385 | — | **Lipase** | Triacylglycerols |
| PFNNT_010615 | *estA* | **Lipase** | Triacylglycerols |
| PFNNT_012680 | — | **Triacylglycerol lipase** | TAG (palm kernel oil) |
| PFNNT_023895 | — | **Lipase family protein** | Broad lipids |
| PFNNT_027350 | — | **Triacylglycerol lipase** (contig_2) | TAG |
| PFNNT_003570 | *cerA* | **Phospholipase CerA** | Phospholipids |
| PFNNT_010430 | — | Phospholipase C/D domain protein | Phospholipids |
| PFNNT_004675 | — | SGNH hydrolase-type esterase | Short-chain esters |
| PFNNT_012505 | *estB* | Esterase EstB | Carboxylesters |
| PFNNT_012965 | — | Carboxylesterase nap | Carboxylesters |
| PFNNT_015515 | *nlhH* | Carboxylesterase NlhH | Carboxylesters |

> [!IMPORTANT]
> **Five lipases** (including two specifically annotated as triacylglycerol lipases) is an exceptionally rich lipolytic toolkit. The triacylglycerol lipase on **contig_2** (PFNNT_027350, the putative conjugative plasmid) is particularly interesting — plasmid-borne lipases are common in industrial *Bacillus* strains and may have been acquired specifically for lipid-rich substrate exploitation. The PFNNT_012680 triacylglycerol lipase is only **700 bp from a TnpB transposase** (IS200/IS605), and **7.5 kb from sonorensin copy 1** — this entire region appears to be a mobile island carrying both lipid degradation and antimicrobial capabilities.

### 4c. Proteases

Two copies of **minor extracellular protease VpR** (PFNNT_012080 and PFNNT_023935, both >4 kb) — large secreted serine proteases for extracellular protein degradation. Palm kernel meal contains ~16% protein, and these proteases enable nitrogen acquisition from the substrate.

---

## 5. Revised Ecological Model

The initial analysis characterized the PF isolate as a "rhizosphere-adapted metal-tolerant soil bacterium." With the palm kernel waste context and the newly highlighted genes, the picture is fundamentally different and far more specific:

### The PF Isolate as a Palm Kernel Waste Specialist

```
                    Palm Kernel Waste
                          |
        +-----------------+-----------------+
        |                 |                 |
   Residual Oil     Lignocellulose      Protein
   (5-8% TAG)      (cellulose,         (~16%)
        |           hemicellulose,          |
        |           lignin)                 |
        v                 |                 v
  5 Lipases              |          2 VpR proteases
  (TAG → FFA)            |          (protein → amino acids)
        |                 |                 |
        |     +-----------+-----------+     |
        |     |           |           |     |
        |  Cellulase  2×Pullulanase  LigB   |
        |  Endogluc.  Neopullul.   2×CatE   |
        |     |        Arabinolact. MhqA×2  |
        |     v           v           v     |
        |   Glucose    Maltose    Protocat.  |
        |     |           |        → TCA     |
        +-----+-----------+---------+-------+
                          |
                    Central Metabolism
                          |
              +-----------+-----------+
              |           |           |
         Bacillibactin  IAA prod.  Bacteriocins
         (Fe capture)  (plant      (competitor
              |        promotion)   killing)
              v           v           v
         Iron in      Root growth  Community
         Fe-limited   stimulation  dominance
         substrate
```

### Why This Matters

1. **The isolate isn't just surviving in PKW — it's adapted to exploit it.** Five lipases for the residual oil, two pullulanases + neopullulanase for the starch, cellulase + endoglucanase for the fibre, and LigB + catechol dioxygenases for the lignin-derived phenolics.

2. **IAA production links waste degradation to plant growth promotion.** The IpdC enzyme converts tryptophan (from protein degradation) into auxin. If this isolate is applied as a biofertilizer alongside composted PKW, it could simultaneously (a) accelerate waste decomposition and (b) stimulate plant root growth via IAA secretion. The MarR regulator near IpdC may tune IAA production in response to phenolic signals released during lignin breakdown.

3. **The bacteriocin trio (Circularin A + 2× sonorensin + LAP maturation) provides competitive exclusion** against other microorganisms colonizing the nutrient-rich waste. This is essential in PKW, which is a non-sterile, high-nutrient substrate quickly colonized by diverse soil microbiota and fungi.

4. **Chitinase is anti-fungal weaponry.** PKW is rapidly colonized by *Aspergillus* and *Rhizopus* species. The chitinase enables the PF isolate to attack fungal competitors directly — their cell wall chitin is the substrate.

5. **The sonorensin-lipase-TnpB island** (~2.48 Mb) represents a potential **acquired metabolic-antimicrobial cassette**: the ability to degrade palm oil AND kill competitors was potentially co-mobilized on a single transposable element.

### Biotechnological Potential

| Application | Genes | Rationale |
|---|---|---|
| **Biofertilizer from PKW compost** | IpdC, dhbA/E/F, enterobactin esterase | IAA production + iron acquisition → plant growth promotion |
| **PKW bioconversion** | 5 lipases, 2 pullulanases, cellulase, endoglucanase | Complete substrate utilization for biorefinery feedstock |
| **Biocontrol agent** | Circularin A, 2× sonorensin, chitinase | Anti-bacterial + anti-fungal activity |
| **Phenolic waste detoxification** | LigB, CatE (×2), MhqA (×2), 2,4-DCP monooxygenase | Degrades toxic lignin phenolics in POME/PKW effluent |
| **Enzyme production** | Lipases (esp. contig_2 TAG lipase), pullulanases | Industrial enzyme candidates for palm oil processing |
