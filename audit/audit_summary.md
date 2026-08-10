# MEOR identifier audit

This report was generated offline from the supplied snapshots. The audit is read-only and does not modify the marker catalog.

## Inputs

- `markers`: `src/genbank_parser/data/meor/markers.yaml`
- `kegg`: `audit-reference/kegg_ko_list_2026-08-02.tsv`
- `cog`: `audit-reference/cog-24.def.tab`
- `ec`: `audit-reference/explorenz-enzyme-data.xml`
- `pfam`: `audit-reference/Pfam-A.hmm`
- `ec_source_used`: `audit-reference/explorenz-enzyme-data.xml`

## Coverage

| Identifier type | Catalog assignments | Reference entries indexed |
| --- | ---: | ---: |
| KO | 49 | 24,514 |
| COG | 0 | 5,061 |
| EC | 53 | 8,343 entries; 8,343 history rows |
| Pfam | 0 | 0 aliases |

## Status counts

| Status | Count |
| --- | ---: |
| `NONSPECIFIC` | 24 |
| `PRESENT_UNREVIEWED` | 93 |

## Detected problems and compatibility qualifications

| Marker | Type | Identifier | Status | Notes |
| --- | --- | --- | --- | --- |
| `prmABCD` | `EC` | `1.14.13.-` | `NONSPECIFIC` | wildcard EC is broad contextual evidence; shared across markers: prmABCD, pbmoABC, sbmoXYZ, tmo_tod, isoA_I, xamo_xec |
| `pbmoABC` | `EC` | `1.14.13.-` | `NONSPECIFIC` | wildcard EC is broad contextual evidence; shared across markers: prmABCD, pbmoABC, sbmoXYZ, tmo_tod, isoA_I, xamo_xec |
| `sbmoXYZ` | `EC` | `1.14.13.-` | `NONSPECIFIC` | wildcard EC is broad contextual evidence; shared across markers: prmABCD, pbmoABC, sbmoXYZ, tmo_tod, isoA_I, xamo_xec |
| `ahyA` | `EC` | `1.17.99.-` | `NONSPECIFIC` | wildcard EC is broad contextual evidence |
| `CYP153` | `EC` | `1.14.15.-` | `NONSPECIFIC` | wildcard EC is broad contextual evidence |
| `almA` | `EC` | `1.14.14.-` | `NONSPECIFIC` | wildcard EC is broad contextual evidence |
| `tmo_tod` | `EC` | `1.14.13.-` | `NONSPECIFIC` | wildcard EC is broad contextual evidence; shared across markers: prmABCD, pbmoABC, sbmoXYZ, tmo_tod, isoA_I, xamo_xec |
| `mahAB` | `EC` | `1.14.12.-` | `NONSPECIFIC` | wildcard EC is broad contextual evidence |
| `nmsA` | `EC` | `4.1.99.-` | `NONSPECIFIC` | wildcard EC is broad contextual evidence |
| `abcA` | `EC` | `4.1.1.-` | `NONSPECIFIC` | wildcard EC is broad contextual evidence; shared across markers: abcA, ncrA |
| `ncrA` | `EC` | `4.1.1.-` | `NONSPECIFIC` | wildcard EC is broad contextual evidence; shared across markers: abcA, ncrA |
| `emt1` | `EC` | `2.4.1.-` | `NONSPECIFIC` | wildcard EC is broad contextual evidence; shared across markers: emt1, sble |
| `sble` | `EC` | `2.4.1.-` | `NONSPECIFIC` | wildcard EC is broad contextual evidence; shared across markers: emt1, sble |
| `sble` | `EC` | `2.3.1.-` | `NONSPECIFIC` | wildcard EC is broad contextual evidence |
| `srfA` | `EC` | `2.7.7.-` | `NONSPECIFIC` | wildcard EC is broad contextual evidence; shared across markers: srfA, licA_D, fengycin_iturin, mycA_C, arfA_C, pswP |
| `srfA` | `EC` | `6.3.2.-` | `NONSPECIFIC` | wildcard EC is broad contextual evidence |
| `licA_D` | `EC` | `2.7.7.-` | `NONSPECIFIC` | wildcard EC is broad contextual evidence; shared across markers: srfA, licA_D, fengycin_iturin, mycA_C, arfA_C, pswP |
| `fengycin_iturin` | `EC` | `2.7.7.-` | `NONSPECIFIC` | wildcard EC is broad contextual evidence; shared across markers: srfA, licA_D, fengycin_iturin, mycA_C, arfA_C, pswP |
| `mycA_C` | `EC` | `2.7.7.-` | `NONSPECIFIC` | wildcard EC is broad contextual evidence; shared across markers: srfA, licA_D, fengycin_iturin, mycA_C, arfA_C, pswP |
| `arfA_C` | `EC` | `2.7.7.-` | `NONSPECIFIC` | wildcard EC is broad contextual evidence; shared across markers: srfA, licA_D, fengycin_iturin, mycA_C, arfA_C, pswP |
| `pswP` | `EC` | `2.7.7.-` | `NONSPECIFIC` | wildcard EC is broad contextual evidence; shared across markers: srfA, licA_D, fengycin_iturin, mycA_C, arfA_C, pswP |
| `wza_wzb_wzc` | `EC` | `3.1.3.-` | `NONSPECIFIC` | wildcard EC is broad contextual evidence |
| `isoA_I` | `EC` | `1.14.13.-` | `NONSPECIFIC` | wildcard EC is broad contextual evidence; shared across markers: prmABCD, pbmoABC, sbmoXYZ, tmo_tod, isoA_I, xamo_xec |
| `xamo_xec` | `EC` | `1.14.13.-` | `NONSPECIFIC` | wildcard EC is broad contextual evidence; shared across markers: prmABCD, pbmoABC, sbmoXYZ, tmo_tod, isoA_I, xamo_xec |
| `tmo_tod` | `CONSISTENCY` | `tmo_tod` | `PRESENT_UNREVIEWED` | at least one KO EC matches; additional EC scope is broader |
| `ndo_nah` | `CONSISTENCY` | `ndo_nah` | `PRESENT_UNREVIEWED` | at least one KO EC matches; additional EC scope is broader |
| `trehalolipid` | `CONSISTENCY` | `trehalolipid` | `PRESENT_UNREVIEWED` | at least one KO EC matches; additional EC scope is broader |

## Reference index counts

- KEGG definitions: 24,514 (embedded EC annotations are parsed from `[EC:...]` brackets).
- COG definitions: 5,061; categories are preserved as complete strings.
- EC entries: 8,343; history rows: 8,343; class headings: 405.
- Pfam aliases: 0; aliases include versioned and unversioned accessions.
- Pfam indexing retains only catalog-requested accessions; the current catalog has no Pfam assignments, so the supplied Pfam snapshot was validated for presence but no model aliases were retained.

## Notes

`PRESENT_UNREVIEWED` means that an identifier exists and no implemented automatic rule disproved it; it is not curator-approved semantic validation. `NONSPECIFIC` wildcard ECs are retained only as broad contextual metadata and are excluded from structured Weight-3 matching. `DUPLICATED` identifies reuse across markers; `REVIEW` requires curator judgment.

COG/Pfam reference parsing is implemented, but the current MEOR catalog defines zero COG or Pfam assignments. The runtime model does not consume COG/Pfam evidence.
