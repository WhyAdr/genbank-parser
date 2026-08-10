# MEOR & Hydrocarbon Biosynthesis Marker Reference

This reference documents MEOR catalog `1.1`: 48 marker definitions across 9 functional categories and 7 pathway-completeness models implemented by `gbparse meor` in `genbank-parser` v0.4.1. The catalog was initially migrated for behavioral parity from `WhyAdr/genbank-meor` commit `42d09105ca63a47dac78e0990767124ccde2dab9` (v1.1.0 behavior), then identifier-audited and corrected on 2026-08-10.

## Source Databases

| Source | Release / Version | Stable URL | Access Date |
|---|---|---|---|
| **CANT-HYD** | Release 1.0 (37 HMM profiles) | https://canthyd.ucalgary.ca/ | 2026-08-04 |
| **HMDB** | Hydrocarbon Monooxygenase DB | http://www.hydrocarbonmonooxygenase.com/ | 2026-08-04 |
| **HADEG** | Hydrocarbon Anaerobic & Degrading DB | https://hadeg.tu-braunschweig.de/ | 2026-08-04 |

The identifier catalog was audited offline on 2026-08-10 against the supplied
KEGG KO snapshot (`kegg_ko_list_2026-08-02.tsv`), COG-2024 definitions,
ExplorEnz EC data, and Pfam-A headers. Wildcard ECs are preserved as contextual
catalog metadata but cannot produce structured Weight-3 matches. The catalog
defines no COG or Pfam mappings; their reference parsers are audit support only,
and the runtime MEOR model does not consume them.

The audit removed direct definition mismatches (including stale AlkB,
lipopeptide, anaerobic-activation, glycolipid, and transfer-system KOs),
assigned KOs `K15666`-`K15668` to the fengycin/iturin marker, and corrected the
protocatechuate and dibenzothiophene EC numbers. Unresolved `ncrA` KO `K27540`,
transferred `alkJ` EC `1.1.99.8`, LasI KO `K13060` plus nondiscriminating EC
`2.3.1.184` on `rhlRI`, and inconsistent `etnE` KO `K22363` were removed from
active matching without guessed replacements. `rhlRI` retains structured KO
support for RhlI (`K13061`); RhlR is supported by gene/product annotations only.

### Audit status semantics

- `PRESENT_UNREVIEWED`: present in the supplied reference and not disproved by an automated rule; not a curated biological validation.
- `NONSPECIFIC`: broad identifier retained as context but excluded from structured Weight-3 matching.
- `REVIEW`: internally inconsistent or otherwise requiring curator judgment; disabled from active evidence in catalog `1.1`.
- `REPLACE`: deleted or transferred identifier; disabled until a unique successor is curated.
- `MISSING_FROM_REFERENCE`: absent from the supplied snapshot; disabled from active evidence.

Exact snapshot sizes and SHA-256 values are recorded in [`audit/reference_manifest.tsv`](../audit/reference_manifest.tsv).

---

## 48-Marker Scientific Provenance & Traceability Matrix

| Marker ID | Marker Name | Primary Source Database | Primary Annotation & Functional Interpretation Note |
|---|---|---|---|
| `pmoABC` | Particulate Methane Monooxygenase | CANT-HYD, HMDB | Short-chain alkane C1 methane oxidation (pmoA/B/C) |
| `mmoXYZ` | Soluble Methane Monooxygenase | CANT-HYD, HMDB | Soluble methane oxidation hydroxylase (mmoX/Y/Z/B/C/D) |
| `prmABCD` | Propane Monooxygenase | CANT-HYD, HMDB | C3 propane monooxygenation (prmA/B/C/D / bmo) |
| `pbmoABC` | Membrane Alkane/Butane Monooxygenase | HMDB | Membrane-bound butane & alkane monooxygenase |
| `sbmoXYZ` | Soluble Alkane Monooxygenase | HMDB | Soluble alkane monooxygenase subunits |
| `ahyA` | Molybdopterin Alkane C2 Hydroxylase | CANT-HYD | Anaerobic/aerobic C2 methylene hydroxylase |
| `alkB` | Alkane 1-monooxygenase | CANT-HYD, HMDB | C5-C16 medium alkane hydroxylation (alkB1/alkB2) |
| `CYP153` | Cytochrome P450 Alkane Hydroxylase | CANT-HYD, HMDB | Medium alkane hydroxylation (CYP153A family) |
| `rubAB` | Rubredoxin & Rubredoxin Reductase | CANT-HYD | Electron transfer system for AlkB hydroxylase |
| `alkJ` | Alcohol Dehydrogenase | CANT-HYD | Fatty alcohol to fatty aldehyde oxidation |
| `alkH` | Fatty Aldehyde Dehydrogenase | CANT-HYD | Fatty aldehyde to fatty acid oxidation |
| `alkK` | Fatty Acyl-CoA Synthetase | CANT-HYD | Fatty acid activation to acyl-CoA |
| `alkL` | Outer Membrane Alkane Import | CANT-HYD | Hydrocarbon outer membrane import channel |
| `ladA` | Long-chain Alkane Monooxygenase | CANT-HYD, HMDB | C15-C36 heavy paraffin wax oxidation (flavoprotein LadA/LadB) |
| `almA` | Long-chain Alkane Hydroxylase | CANT-HYD | C20-C32 paraffin wax degradation (flavin-binding AlmA) |
| `tmo_tod` | Toluene / BTEX Monooxygenase | CANT-HYD, HMDB | BTEX aromatic ring activation (tmoABCDEF / todC1C2BA) |
| `ndo_nah` | Naphthalene / PAH Dioxygenase | CANT-HYD | Polycyclic aromatic hydrocarbon initial ring dioxygenation |
| `catA_xylE` | Central Ring-Cleavage Dioxygenase | CANT-HYD | Catechol intradiol (catA) / extradiol (xylE) cleavage |
| `mahAB` | Benzene / Toluene Dioxygenase | CANT-HYD | Monocyclic aromatic ring dioxygenase (mahA/mahB) |
| `dmpO_tomA` | Phenol / Toluene Monooxygenase | CANT-HYD | Substituted phenol & toluene monooxygenation |
| `dszC` | Dibenzothiophene Monooxygenase | HADEG | Dibenzothiophene monooxygenase (4S biodesulfurization) |
| `pcaGH` | Protocatechuate 3,4-Dioxygenase | CANT-HYD | Protocatechuate ortho-cleavage (pcaG/pcaH) |
| `assA_masD` | Alkylsuccinate Synthase | CANT-HYD, HADEG | Anaerobic alkane fumarate addition catalytic subunit |
| `bssABC` | Benzylsuccinate Synthase | CANT-HYD, HADEG | Anaerobic toluene fumarate addition (bssA/bssB/bssC) |
| `nmsA` | Naphthylmethylsuccinate Synthase | CANT-HYD, HADEG | Anaerobic 2-methylnaphthalene fumarate addition |
| `cmdA_ebdA` | Ethylbenzene Dehydrogenase | CANT-HYD, HADEG | Anaerobic ethylbenzene alpha-subunit dehydrogenase |
| `abcA` | Anaerobic Benzene Carboxylase | HADEG | Anaerobic benzene carboxylase subunit A (abcA1/A2) |
| `ncrA` | Naphthalene Carboxylase | HADEG | Anaerobic naphthalene carboxylase (ncrA) |
| `rhlA` | HAA Synthase | HADEG | Rhamnolipid precursor HAA (3-hydroxyacyl-ACP) dimer synthase |
| `rhlB` | Rhamnosyltransferase I | HADEG | Mono-rhamnolipid glucosyltransferase |
| `rhlC` | Rhamnosyltransferase II | HADEG | Di-rhamnolipid glucosyltransferase |
| `rhlRI` | Rhamnolipid Quorum Sensing | HADEG | Rhamnolipid quorum sensing regulatory system (rhlR/rhlI) |
| `trehalolipid` | Trehalolipid Biosynthesis | HADEG | Trehalose glycolipid biosynthesis (sdtA/tshA/otsA/otsB) |
| `emt1` | MEL Biosynthesis | HADEG | Mannosylerythritol lipid erythritol transferase |
| `sble` | Sophorolipid Biosynthesis | HADEG | Sophorolipid glucosyltransferase & esterase |
| `srfA` | Surfactin Synthetase NRPS | HADEG | Surfactin lipopeptide biosurfactant NRPS (srfA-A/B/C/D) |
| `licA_D` | Lichenysin Synthetase NRPS | HADEG | Lichenysin lipopeptide biosurfactant NRPS (licA/B/C/D) |
| `fengycin_iturin` | Fengycin / Iturin NRPS | HADEG | Lipopeptide biosurfactant & antifungal NRPS modules |
| `mycA_C` | Mycosubtilin Synthetase | HADEG | Mycosubtilin lipopeptide NRPS (mycA/B/C) |
| `arfA_C` | Arthrofactin Synthetase | HADEG | Arthrofactin lipopeptide NRPS (arfA/B/C) |
| `pswP` | Serrawettin Synthetase | HADEG | Serrawettin lipopeptide synthetase |
| `sfp` | 4'-Phosphopantetheinyl Transferase | HADEG | NRPS phosphopantetheinyl transferase activator |
| `wza_wzb_wzc` | Polymeric EPS Secretion | HADEG | Bio-emulsifier (emulsan/alasan) exopolysaccharide export |
| `meor_acid_gas` | Organic Acid & Gas Drivers | HADEG, Literature | Formate lyase (pflB), acetate kinase (ackA), carbonic anhydrase |
| `isoA_I` | Isoprene Monooxygenase | HMDB | Gaseous isoprene alkene monooxygenase (isoA-I) |
| `xamo_xec` | Alkene / Ethene Monooxygenase | HMDB | Gaseous alkene & ethene monooxygenase (xamoA-F/xecA-E) |
| `etnE` | Epoxyalkane:CoM Transferase | HMDB | Epoxyalkane coenzyme M conjugation transferase |
| `mpdBC` | 2-Methylpropene Degradation | HMDB | Branched alkene degradation enzymes (mpdB/mpdC) |

---

## 1. Functional Categories & Marker Definitions (48 Markers)

### Category 1: Short-Chain Alkanes & Gaseous Hydrocarbons (C1–C4)
*Role: Methane, ethane, propane, and butane oxidation; reservoir gas utilization.*

1. **pmoABC** — Particulate Methane Monooxygenase (pmoA/pmoB/pmoC)
   - EC: `1.14.18.3` | KO: `K10944`, `K10945`, `K10946`
2. **mmoXYZ** — Soluble Methane Monooxygenase (mmoX/mmoY/mmoZ/mmoB/mmoC/mmoD)
   - EC: `1.14.13.25` | KO: `K16157`, `K16158`, `K16159`, `K16160`, `K16161`, `K16162`
3. **prmABCD** — Propane Monooxygenase (prmA/prmB/prmC/prmD / bmo)
   - EC: `1.14.13.-`
4. **pbmoABC** — Membrane-bound Alkane/Butane Monooxygenase (pbmoA/B/C)
   - EC: `1.14.13.-` | KO: `K21320`, `K21321`
5. **sbmoXYZ** — Soluble Alkane Monooxygenase Subunits (sbmoX/Y/Z)
   - EC: `1.14.13.-`
6. **ahyA** — Molybdopterin-family Alkane C2 Methylene Hydroxylase (ahyA)
   - EC: `1.17.99.-`

---

### Category 2: Medium-Chain n-Alkanes (C5–C16)
*Role: Liquid alkane degradation and crude oil viscosity reduction.*

7. **alkB** — Alkane 1-monooxygenase (alkB / alkB1 / alkB2)
   - EC: `1.14.15.3` | KO: `K00496`
8. **CYP153** — Cytochrome P450 Alkane Hydroxylase (CYP153A)
   - EC: `1.14.15.-`
9. **rubAB** — Rubredoxin & Rubredoxin Reductase (rubA/rubB)
   - EC: `1.18.1.1`
10. **alkJ** — Alcohol Dehydrogenase (alkJ / Alkane pathway)
    - EC: `1.1.1.1` | KO: `K00001`, `K13953`
11. **alkH** — Fatty Aldehyde Dehydrogenase (alkH)
    - EC: `1.2.1.3` | KO: `K00128`
12. **alkK** — Fatty Acyl-CoA Synthetase (alkK)
    - EC: `6.2.1.3` | KO: `K01897`
13. **alkL** — Outer Membrane Alkane Import Protein (alkL)

---

### Category 3: Long-Chain n-Alkanes & Heavy Paraffins (C18–C36+)
*Role: Heavy paraffin wax degradation, pour point reduction, and wax dispersion.*

14. **ladA** — Long-chain Alkane Monooxygenase LadA/LadB (Flavoprotein C15-C36)
    - EC: `1.14.14.28` | KO: `K20938`
15. **almA** — Long-chain Alkane Hydroxylase AlmA (Flavin-binding C20-C32)
    - EC: `1.14.14.-`

---

### Category 4: Aromatic & Polycyclic Aromatic Hydrocarbons (BTEX & PAHs)
*Role: Aromatic fraction breakdown, heavy crude liquefaction, and biodesulfurization.*

16. **tmo_tod** — Toluene / BTEX Monooxygenase & Dioxygenase (tmoABCDEF / todC1C2BA / xylMA)
    - EC: `1.14.13.-` | KO: `K15760`, `K15761`, `K15762`, `K03268`, `K16268`, `K18089`, `K18090`, `K16269`, `K16270`
17. **ndo_nah** — Naphthalene / PAH Dioxygenase (ndoABC / nahAcAb / phnAc / bphA1-4)
    - EC: `1.14.12.12`, `1.14.12.18` | KO: `K14579`, `K14580`, `K14581`
18. **catA_xylE** — Central Ring-Cleavage Dioxygenases (catA / xylE)
    - EC: `1.13.11.1`, `1.13.11.2` | KO: `K03381`, `K00446`
19. **mahAB** — Benzene / Toluene / Naphthalene Dioxygenase (mahA / mahB)
    - EC: `1.14.12.-`
20. **dmpO_tomA** — Phenol / Toluene 2-Monooxygenase (dmpO / tomA1-A4)
    - EC: `1.14.13.7`
21. **dszC** — Dibenzothiophene Monooxygenase (dszC / Biodesulfurization)
    - EC: `1.14.14.21` | KO: `K22219`
22. **pcaGH** — Protocatechuate 3,4-Dioxygenase (pcaG / pcaH)
    - EC: `1.13.11.3` | KO: `K00448`, `K00449`

---

### Category 5: Anaerobic Hydrocarbon Activation
*Role: Anoxic in-situ hydrocarbon activation via fumarate addition under reservoir conditions.*

23. **assA_masD** — Alkylsuccinate Synthase catalytic alpha subunit (assA / masD)
    - EC: `4.1.99.16`
24. **bssABC** — Benzylsuccinate Synthase (bssA / bssB / bssC)
    - EC: `4.1.99.11` | KO: `K07540`
25. **nmsA** — Naphthylmethylsuccinate Synthase (nmsA)
    - EC: `4.1.99.-`
26. **cmdA_ebdA** — Ethylbenzene Dehydrogenase Subunit Alpha (cmdA / ebdA)
    - EC: `1.17.99.2`
27. **abcA** — Anaerobic Benzene Carboxylase Subunit A (abcA1 / abcA2)
    - EC: `4.1.1.-`
28. **ncrA** — Naphthalene Carboxylase (ncrA)
    - EC: `4.1.1.-` | KO: none currently active

---

### Category 6: Glycolipid Biosurfactants
*Role: Surface tension reduction, interfacial tension (IFT) drop, and oil bio-emulsification.*

29. **rhlA** — HAA Synthase (rhlA / Rhamnolipid fatty acid precursor)
    - EC: `2.3.1.266`
30. **rhlB** — Rhamnosyltransferase I (rhlB / Mono-rhamnolipid)
    - EC: `2.4.1.189`
31. **rhlC** — Rhamnosyltransferase II (rhlC / Di-rhamnolipid)
    - EC: `2.4.1.298`
32. **rhlRI** — Rhamnolipid Quorum Sensing Regulators (rhlR / rhlI)
    - EC: none currently active | KO: `K13061` (RhlI; RhlR has gene/product evidence only)
33. **trehalolipid** — Trehalolipid Biosynthesis (sdtA / tshA / otsA / otsB)
    - EC: `2.4.1.15`, `3.1.3.12` | KO: `K00697`, `K01087`
34. **emt1** — Mannosylerythritol Lipid (MEL) Biosynthesis (emt1)
    - EC: `2.4.1.-`
35. **sble** — Sophorolipid Biosynthesis (sble / glucosyltransferase)
    - EC: `2.4.1.-`, `2.3.1.-`

---

### Category 7: Lipopeptide Biosurfactant NRPS BGCs
*Role: Potent lipopeptide biosurfactants (surfactin, lichenysin) for ultra-low IFT reduction.*

36. **srfA** — Surfactin Synthetase NRPS (srfA-A / srfA-B / srfA-C / srfAD)
    - EC: `2.7.7.-`, `6.3.2.-`
37. **licA_D** — Lichenysin Synthetase NRPS (licA / licB / licC / licD)
    - EC: `2.7.7.-`
38. **fengycin_iturin** — Fengycin / Iturin / Viscosin / Plipastatin NRPS (fenA-E / ituA-C / vsnA-C / ppsA-E)
    - EC: `2.7.7.-` | KO: `K15666`, `K15667`, `K15668`
39. **mycA_C** — Mycosubtilin Lipopeptide Synthetase (mycA / mycB / mycC)
    - EC: `2.7.7.-` | KO: `K15661`
40. **arfA_C** — Arthrofactin Lipopeptide Synthetase (arfA / arfB / arfC)
    - EC: `2.7.7.-` | KO: `K15658`
41. **pswP** — Serrawettin Lipopeptide Synthetase (pswP)
    - EC: `2.7.7.-`
42. **sfp** — 4'-Phosphopantetheinyl Transferase (sfp / NRPS activator)
    - EC: `2.7.8.7`

---

### Category 8: Polymeric Bio-emulsifiers & Organic Acid / Gas Production
*Role: Emulsification, carbonate rock dissolution, acid flooding, and reservoir repressurization.*

43. **wza_wzb_wzc** — Polymeric Bio-emulsifier Secretion (wza / wzb / wzc / EPS / Emulsan / Alasan)
    - EC: `3.1.3.-` | KO: `K01991`
44. **meor_acid_gas** — Organic Acid & Gas Drivers (pflB / ackA / pta / ca)
    - EC: `2.3.1.54`, `2.7.2.1`, `2.3.1.8`, `4.2.1.1` | KO: `K00656`, `K00925`, `K00625`, `K01672`, `K01673`

---

### Category 9: Gaseous & Aliphatic Alkenes
*Role: Alkene monooxygenation, epoxide ring-opening, and isoprene degradation.*

45. **isoA_I** — Isoprene Monooxygenase (isoA-I)
    - EC: `1.14.13.-`
46. **xamo_xec** — Alkene / Ethene Monooxygenase (xamoA-F / xecA-E)
    - EC: `1.14.13.-`
47. **etnE** — Epoxyalkane:CoM Transferase (etnE)
    - EC: `2.5.1.56` | KO: none currently active
48. **mpdBC** — 2-Methylpropene Degradation (mpdB / mpdC)

---

## 2. Predefined Pathway Completeness Models

1. **Medium Alkane Degradation Chain (C5–C16)** (6 steps)
   - Step 1: Primary Hydroxylation (`alkB` / `CYP153`)
   - Step 2: Electron Transfer (`rubAB`)
   - Step 3: Alcohol Oxidation (`alkJ`)
   - Step 4: Aldehyde Oxidation (`alkH`)
   - Step 5: Acyl-CoA Activation (`alkK`)
   - Step 6: Outer Membrane Import (`alkL`)
2. **Long-Chain Paraffin Wax Degradation (C18+)** (2 steps)
   - Step 1: Monooxygenase (`ladA`)
   - Step 2: Hydroxylase (`almA`)
3. **Aromatic & PAH Degradation** (2 steps)
   - Step 1: Ring Activation (`tmo_tod`, `ndo_nah`, `mahAB`, `dmpO_tomA`)
   - Step 2: Ring Cleavage (`catA_xylE`)
4. **Rhamnolipid Biosurfactant Biosynthesis** (3 steps)
   - Step 1: Precursor HAA Synthase (`rhlA`)
   - Step 2: Mono-rhamnolipid Synthase (`rhlB`)
   - Step 3: Di-rhamnolipid Synthase (`rhlC`)
5. **Lipopeptide Biosurfactant BGC** (2 steps)
   - Step 1: NRPS Modules (`srfA`, `licA_D`)
   - Step 2: PPTase Activator (`sfp`)
6. **Anaerobic Alkane Activation (AssA)** (1 step)
   - Step 1: Catalytic Subunit (`assA_masD`)
7. **Gaseous & Aliphatic Alkene Degradation** (2 steps)
   - Step 1: Monooxygenation (`isoA_I`, `xamo_xec`)
   - Step 2: Epoxide Conjugation (`etnE`)

---

## 3. Evidence Weighting Schema

Matches from `gbparse meor` are assigned confidence weights from 1 to 3:
- **Weight 3 (High Confidence)**: Match on an active structured KEGG KO, fully specified INSDC EC number, or exact gene symbol regex.
- **Weight 2 (Medium Confidence)**: Match on specific product keyword regex.
- **Weight 1 (Low Confidence)**: Match on broad note text, including a KO or EC identifier found only in `/note`.

Wildcard ECs are contextual catalog metadata and never produce Weight-3 evidence.

These annotation-supported candidates indicate genomic encoding potential. They do not demonstrate transcription, enzyme activity, hydrocarbon turnover, biosurfactant production, or enhanced oil recovery performance. The pathway scores are genome-level marker completeness, not cluster-local completeness or biochemical flux.
