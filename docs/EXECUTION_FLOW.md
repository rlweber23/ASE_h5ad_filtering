# Execution Flow

## How shell wrappers call Python scripts

All wrappers pass explicit CLI arguments into the Python processing scripts.

F1 local wrapper:
- `workflows/f1/run_f1_local.sh <TISSUE>`
- Calls `scripts/F1_make_strain_ASE_adata.py`

F1 SLURM wrapper:
- `workflows/f1/run_f1_slurm.sh`
- Uses a tissue array and `SLURM_ARRAY_TASK_ID`
- Calls `scripts/F1_make_strain_ASE_adata.py`

Founder local wrapper:
- `workflows/founder/run_founder_local.sh <TISSUE>`
- Calls `scripts/Founder_make_strain_ASE_adata.py`

Founder SLURM wrapper:
- `workflows/founder/run_founder_slurm.sh`
- Uses a tissue array and `SLURM_ARRAY_TASK_ID`
- Calls `scripts/Founder_make_strain_ASE_adata.py`

## Runtime configuration

Wrappers read environment variables and pass values through Python CLI options:
- `ASE_BASE_PATH`
- `ASE_ANNOTATION_DIR`
- `ASE_OUT_DIR`
- Optional reference overrides:
  - `ASE_BARCODE_MAP`
  - `ASE_GENES_WITH_BIOTYPE`
  - `ASE_GENE_ID_CHR_MAP`

This avoids relying on internal hardcoded defaults and keeps Python logic unchanged.
