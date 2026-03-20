# ASE Filtering Scripts (Publication Reference Repository)

This repository contains scripts and wrappers used to build allele-specific single-nucleus RNA-seq AnnData outputs for founder and F1 mouse workflows.

## Repository layout

- `scripts/`: original Python processing scripts (logic preserved, no refactoring).
- `workflows/f1/`: F1 local and SLURM wrappers.
- `workflows/founder/`: founder local and SLURM wrappers.
- `refs/`: lightweight reference tables used by wrappers by default.
- `config/`: environment template for runtime paths.
- `docs/`: pipeline and I/O documentation.

## Pipelines

## Founder processing

Input context:
- Founder/inbred samples mapped to concatenated allele-aware references.

Core behavior:
1. Restrict processing to strain-relevant mappings represented in strain directories.
2. Build and align cell IDs from barcode segments.
3. Join founder annotation metadata.
4. Retain B6J plus founder-specific genotype labels.
5. Generate tissue- and strain-specific founder ASE h5ad outputs.

Entry points:
- Local: `workflows/founder/run_founder_local.sh`
- SLURM: `workflows/founder/run_founder_slurm.sh`

## F1 processing

Input context:
- F1 samples mapped to concatenated allele-aware references.

Core behavior:
1. Restrict to strain-relevant mappings.
2. Keep allele-specific information for both parental alleles.
3. Build and align cell IDs from barcode segments.
4. Join F1 annotation metadata.
5. Retain B6J plus strain-specific F1 genotype labels.
6. Generate tissue- and strain-specific F1 ASE h5ad outputs.

Entry points:
- Local: `workflows/f1/run_f1_local.sh`
- SLURM: `workflows/f1/run_f1_slurm.sh`

## Setup

1. Ensure Python environment has required packages:
   - `numpy`
   - `pandas`
   - `scanpy`
   - `anndata`
2. Set runtime paths using environment variables.

Example:

```bash
source config/env.example
# Edit values in your shell session or in a private env file.
```

Required variables:
- `ASE_BASE_PATH`
- `ASE_ANNOTATION_DIR`
- `ASE_OUT_DIR`

Optional overrides:
- `ASE_BARCODE_MAP` (default `refs/barcode_map_df.csv`)
- `ASE_GENES_WITH_BIOTYPE` (default `refs/igvf_gtf_gene_biotypes.tsv`)
- `ASE_GENE_ID_CHR_MAP` (default `refs/gencode_M36_gene_ID_chromosome_map.tsv`)

## Example usage

Run F1 for one tissue:

```bash
workflows/f1/run_f1_local.sh Heart
```

Run founder for one tissue:

```bash
workflows/founder/run_founder_local.sh Heart
```

Run on SLURM:

```bash
sbatch workflows/f1/run_f1_slurm.sh
sbatch workflows/founder/run_founder_slurm.sh
```

Adjust tissue arrays inside SLURM scripts as needed.

## Execution flow

Wrappers in `workflows/` call Python scripts in `scripts/` and pass explicit CLI path arguments. This removes dependence on embedded default paths and keeps code reproducible without modifying internal Python logic.

Details:
- `docs/EXECUTION_FLOW.md`
- `docs/INPUT_OUTPUT_SPEC.md`
- `docs/PIPELINE_OVERVIEW.md`
- `docs/HARDCODED_PATHS.md`

## Notes for reproducibility

- Python script internals are unchanged.
- Hardcoded default paths still exist inside Python argparse blocks, but wrappers override them with explicit values.
- Large data outputs (`*.h5ad`) are excluded via `.gitignore`.
