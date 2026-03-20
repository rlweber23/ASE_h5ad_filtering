# Guide

## Purpose

This repository documents two scripts used in a publication pipeline to build allele-specific expression AnnData files from concatenated-genome single-nucleus RNA-seq results.

The goal is explanation and reference, not packaging a fully automated framework.

## What each script does

### Founder script

Script: `scripts/Founder_make_strain_ASE_adata.py`

1. Uses founder annotations (`{tissue}_Founder_processed_annotated_obs.csv`).
2. For each tissue, combines one or more plate directories.
3. Reads strain-level h5ad files.
4. Splits each barcode into three 8bp segments and maps each segment to well IDs.
5. Reconstructs `cellID` from wells, subpool name, and plate ID.
6. Merges annotation metadata onto cells.
7. Keeps only cells from `B6J` or the founder-specific short genotype code.
8. Assigns allele labels for B6J vs alternate strain genes.
9. Removes chrY, chrM, and listed random contigs.
10. Writes `{strain}_{tissue}_Founder_ASE.h5ad`.

### F1 script

Script: `scripts/F1_make_strain_ASE_adata.py`

1. Uses F1 annotations (`{tissue}_F1_processed_annotated_obs.csv`).
2. For each tissue, uses one plate directory.
3. Reads strain-level h5ad files.
4. Performs the same barcode parsing and `cellID` reconstruction.
5. Merges annotation metadata.
6. Keeps only cells from `B6J` or the expected F1 genotype for that strain pair.
7. Assigns allele labels for B6J vs alternate strain genes.
8. Removes chrY, chrM, and listed random contigs.
9. Writes `{strain}_{tissue}_F1_ASE.h5ad`.

## Required inputs

1. Base directory containing h5ad files organized by plate and strain.
2. Annotation directory containing per-tissue annotation CSVs.
3. Barcode map CSV with columns `bc`, `bc_sequence`, `bc_well`.
4. Gene biotype TSV (gene ID, gene name, biotype).
5. Gene chromosome map TSV (chromosome and gene ID with version).

## Direct run commands

Run F1 for one tissue:

```bash
python scripts/F1_make_strain_ASE_adata.py \
  --tissue Heart \
  --base-path /path/to/cellbender \
  --barcode-map refs/barcode_map_df.csv \
  --genes-with-biotype refs/igvf_gtf_gene_biotypes.tsv \
  --gene-id-chr-map refs/gencode_M36_gene_ID_chromosome_map.tsv \
  --annotation-dir /path/to/annotations \
  --out-dir /path/to/output
```

Run founder for one tissue:

```bash
python scripts/Founder_make_strain_ASE_adata.py \
  --tissue Heart \
  --base-path /path/to/cellbender \
  --barcode-map refs/barcode_map_df.csv \
  --genes-with-biotype refs/igvf_gtf_gene_biotypes.tsv \
  --gene-id-chr-map refs/gencode_M36_gene_ID_chromosome_map.tsv \
  --annotation-dir /path/to/annotations \
  --out-dir /path/to/output
```

## Outputs

Both scripts write strain-by-tissue h5ad files under the output directory, with names:

- Founder: `{strain}_{tissue}_Founder_ASE.h5ad`
- F1: `{strain}_{tissue}_F1_ASE.h5ad`

## Notes

1. The Python script internals are intentionally unchanged from analysis usage.
2. Both scripts contain environment-specific default paths; pass explicit CLI paths in practice.
3. Tissue names must match the hardcoded tissue keys in each script.
