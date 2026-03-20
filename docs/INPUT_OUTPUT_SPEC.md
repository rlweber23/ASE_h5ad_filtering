# Input and Output Specification

## Shared required inputs

1. Base path containing plate and strain directories with h5ad files.
2. Barcode map CSV with columns: `bc`, `bc_sequence`, `bc_well`.
3. Gene biotype TSV (no header): `gene_id`, `gene_name`, `biotype`.
4. Gene-to-chromosome TSV (no header): `chromosome`, `gene_id_with_version`.
5. Annotation CSV with `cellID` column and genotype metadata.

## h5ad assumptions

- `obs` must include `barcode`.
- `var` should include `gene_id`; if missing, scripts derive from `var_names`.

## Annotation naming conventions

- F1: `{tissue}_F1_processed_annotated_obs.csv`
- Founder: `{tissue}_Founder_processed_annotated_obs.csv`

## Primary generated outputs

- F1: `{out_dir}/{tissue}/{strain}_{tissue}_F1_ASE.h5ad`
- Founder: `{out_dir}/{tissue}/{strain}_{tissue}_Founder_ASE.h5ad`

## Added/derived columns

Added in `obs`:
- `bc1_sequence`, `bc2_sequence`, `bc3_sequence`
- `bc1_well`, `bc2_well`, `bc3_well`
- `subpool`
- `cellID`

Added in `var`:
- `gene_id_unmodified`
- `gene_name_unmodified`
- `allele`
- `B6J_allele`
- `ALT_allele`
- `chromosome`
