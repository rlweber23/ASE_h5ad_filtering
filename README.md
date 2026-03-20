# ASE Filtering Scripts

This is a reference repository for two publication scripts that generate allele-specific expression (ASE) AnnData outputs from concatenated-genome snRNA-seq alignments.


## Simple layout

- `scripts/`: unchanged Python scripts used in the analysis.
- `refs/`: lightweight reference files used by both scripts.
- `docs/GUIDE.md`: concise step-by-step explanation of founder and F1 workflows.

## Founder vs F1

Founder workflow (`scripts/Founder_make_strain_ASE_adata.py`):
- Uses founder annotation files.
- Combines multiple plates per tissue.
- Retains B6J plus founder-specific genotype labels.
- Writes `{strain}_{tissue}_Founder_ASE.h5ad`.

F1 workflow (`scripts/F1_make_strain_ASE_adata.py`):
- Uses F1 annotation files.
- Uses one plate per tissue.
- Retains B6J plus expected F1 genotype labels.
- Writes `{strain}_{tissue}_F1_ASE.h5ad`.

Both scripts:
- Parse barcodes into bc1/bc2/bc3 segments.
- Reconstruct `cellID` and merge annotations.
- Label alleles in gene metadata.
- Remove chrY/chrM and listed random contigs.

## Run directly (no wrapper required)

F1 example:

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

Founder example:

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


