# Pipeline Overview

This repository contains two related pipelines for building allele-specific expression AnnData outputs from concatenated-genome mapping results.

## Founder pipeline

Input assumption:
- Reads are mapped to a concatenated reference including B6J and alternate founder alleles.

Main behavior:
1. Load founder tissue annotation file: `{tissue}_Founder_processed_annotated_obs.csv`.
2. Read h5ad files across one or more plate IDs for each founder strain.
3. Parse 24bp barcode into bc1, bc2, bc3 segments.
4. Map segments to wells using `refs/barcode_map_df.csv`.
5. Reconstruct `cellID` and merge with annotation metadata.
6. Keep genotype rows matching `B6J` plus founder short code (for example `AJ`, `CASTJ`).
7. Set allele metadata in var (`B6J_allele`, `ALT_allele`, `allele`).
8. Remove chrY, chrM, and random contigs defined in script.
9. Write `{strain}_{tissue}_Founder_ASE.h5ad`.

## F1 pipeline

Input assumption:
- Reads are mapped to a concatenated reference including B6J and alternate founder alleles.

Main behavior:
1. Load F1 tissue annotation file: `{tissue}_F1_processed_annotated_obs.csv`.
2. Read h5ad files for a single plate ID per tissue and per strain.
3. Parse and map barcodes, reconstruct `cellID`, and merge annotations.
4. Keep genotype rows matching `B6J` plus the expected F1 genotype for that founder pair (for example `B6CASTF1J`).
5. Keep allele metadata for both parental alleles in var.
6. Remove chrY, chrM, and random contigs defined in script.
7. Write `{strain}_{tissue}_F1_ASE.h5ad`.

## Workflow separation

- F1 wrappers: `workflows/f1/`
- Founder wrappers: `workflows/founder/`
- Original processing scripts (unchanged): `scripts/`
