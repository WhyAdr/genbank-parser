# MEOR identifier audit

This report was generated offline from the supplied snapshots. The audit is read-only and does not modify the marker catalog.

## Inputs

- `markers`: `src\genbank_parser\data\meor\markers.yaml`
- `kegg`: `audit-reference\kegg_ko_list_2026-08-02.tsv`
- `cog`: `audit-reference\cog-24.def.tab`
- `ec_xml`: `audit-reference\explorenz-enzyme-data.xml`
- `pfam`: `audit-reference\Pfam-A.hmm`
- `ec_source_used`: `audit-reference\explorenz-enzyme-data.xml`

## Coverage

| Identifier type | Catalog assignments | Reference entries indexed |
| --- | ---: | ---: |
| KO | 52 | 24,514 |
| COG | 0 | 5,061 |
| EC | 55 | 8,343 entries; 8,343 history rows |
| Pfam | 0 | 0 aliases |

## Status counts

| Status | Count |
| --- | ---: |
| `MISSING_FROM_REFERENCE` | 1 |
| `NONSPECIFIC` | 24 |
| `PASS` | 71 |
| `REPLACE` | 1 |
| `REVIEW` | 14 |

## Findings requiring review or remediation

| Marker | Type | Identifier | Status | Notes |
| --- | --- | --- | --- | --- |
| `prmABCD` | `EC` | `1.14.13.-` | `NONSPECIFIC` | wildcard EC is broad contextual evidence; shared across markers: prmABCD, pbmoABC, sbmoXYZ, tmo_tod, isoA_I, xamo_xec |
| `pbmoABC` | `EC` | `1.14.13.-` | `NONSPECIFIC` | wildcard EC is broad contextual evidence; shared across markers: prmABCD, pbmoABC, sbmoXYZ, tmo_tod, isoA_I, xamo_xec |
| `sbmoXYZ` | `EC` | `1.14.13.-` | `NONSPECIFIC` | wildcard EC is broad contextual evidence; shared across markers: prmABCD, pbmoABC, sbmoXYZ, tmo_tod, isoA_I, xamo_xec |
| `ahyA` | `EC` | `1.17.99.-` | `NONSPECIFIC` | wildcard EC is broad contextual evidence |
| `CYP153` | `EC` | `1.14.15.-` | `NONSPECIFIC` | wildcard EC is broad contextual evidence |
| `alkJ` | `EC` | `1.1.99.8` | `REPLACE` | transferred EC; successor(s): 1.1.2.7, 1.1.2.8 |
| `almA` | `EC` | `1.14.14.-` | `NONSPECIFIC` | wildcard EC is broad contextual evidence |
| `tmo_tod` | `KO` | `K03268` | `REVIEW` | embedded ECs: 1.14.12.3, 1.14.12.11, 1.14.12.26; KO/EC contradiction |
| `tmo_tod` | `KO` | `K16268` | `REVIEW` | embedded ECs: 1.14.12.3, 1.14.12.11, 1.14.12.26; KO/EC contradiction |
| `tmo_tod` | `KO` | `K18090` | `REVIEW` | embedded ECs: 1.18.1.3, 1.18.1.-; KO/EC contradiction |
| `tmo_tod` | `KO` | `K16269` | `REVIEW` | embedded ECs: 1.3.1.19, 1.3.1.119; KO/EC contradiction |
| `tmo_tod` | `KO` | `K16270` | `REVIEW` | embedded ECs: 1.13.11.-; KO/EC contradiction |
| `tmo_tod` | `EC` | `1.14.13.-` | `NONSPECIFIC` | wildcard EC is broad contextual evidence; shared across markers: prmABCD, pbmoABC, sbmoXYZ, tmo_tod, isoA_I, xamo_xec |
| `ndo_nah` | `KO` | `K14579` | `REVIEW` | embedded ECs: 1.14.12.12, 1.14.12.23, 1.14.12.24; KO/EC contradiction |
| `ndo_nah` | `KO` | `K14580` | `REVIEW` | embedded ECs: 1.14.12.12, 1.14.12.23, 1.14.12.24; KO/EC contradiction |
| `ndo_nah` | `KO` | `K14581` | `REVIEW` | embedded ECs: 1.18.1.7; KO/EC contradiction |
| `mahAB` | `EC` | `1.14.12.-` | `NONSPECIFIC` | wildcard EC is broad contextual evidence |
| `nmsA` | `EC` | `4.1.99.-` | `NONSPECIFIC` | wildcard EC is broad contextual evidence |
| `abcA` | `EC` | `4.1.1.-` | `NONSPECIFIC` | wildcard EC is broad contextual evidence; shared across markers: abcA, ncrA |
| `ncrA` | `KO` | `K27540` | `MISSING_FROM_REFERENCE` | KO identifier not present in KEGG snapshot |
| `ncrA` | `EC` | `4.1.1.-` | `NONSPECIFIC` | wildcard EC is broad contextual evidence; shared across markers: abcA, ncrA |
| `trehalolipid` | `KO` | `K00697` | `REVIEW` | embedded ECs: 2.4.1.15, 2.4.1.347; KO/EC contradiction |
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
| `etnE` | `KO` | `K22363` | `REVIEW` | embedded ECs: 4.4.1.23; KO/EC contradiction |
| `tmo_tod` | `CONSISTENCY` | `tmo_tod` | `REVIEW` | KO definition embeds EC(s) not compatible with marker ecs |
| `ndo_nah` | `CONSISTENCY` | `ndo_nah` | `REVIEW` | KO definition embeds EC(s) not compatible with marker ecs |
| `trehalolipid` | `CONSISTENCY` | `trehalolipid` | `REVIEW` | KO definition embeds EC(s) not compatible with marker ecs |
| `etnE` | `CONSISTENCY` | `etnE` | `REVIEW` | KO definition embeds EC(s) not compatible with marker ecs |

## Reference index counts

- KEGG definitions: 24,514 (embedded EC annotations are parsed from `[EC:...]` brackets).
- COG definitions: 5,061; categories are preserved as complete strings.
- EC entries: 8,343; history rows: 8,343; class headings: 405.
- Pfam aliases: 0; aliases include versioned and unversioned accessions.
- Pfam indexing retains only catalog-requested accessions; the current catalog has no Pfam assignments, so the supplied Pfam snapshot was validated for presence but no model aliases were retained.

## Notes

`NONSPECIFIC` wildcard ECs are retained as broad contextual evidence. `DUPLICATED` identifies reuse across markers; a high-confidence incompatible definition is reported as `DEFINITION_MISMATCH`. `REVIEW` consistency findings require curator judgment.
