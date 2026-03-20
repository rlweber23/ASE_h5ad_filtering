# Hardcoded Paths and Portability Notes

## Where hardcoded defaults exist

Path defaults are embedded in argparse definitions in:
- `scripts/F1_make_strain_ASE_adata.py`
- `scripts/Founder_make_strain_ASE_adata.py`

These defaults point to environment-specific locations and should not be used as-is in a public reference setup.

## Portability strategy in this repository

Without editing Python internals, portability is handled by wrappers that always pass explicit CLI values.

Required environment variables:
- `ASE_BASE_PATH`
- `ASE_ANNOTATION_DIR`
- `ASE_OUT_DIR`

Optional environment variables:
- `ASE_BARCODE_MAP`
- `ASE_GENES_WITH_BIOTYPE`
- `ASE_GENE_ID_CHR_MAP`

Set these once per environment and run wrappers from `workflows/`.
