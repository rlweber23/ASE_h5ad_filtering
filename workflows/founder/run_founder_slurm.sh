#!/usr/bin/env bash
#SBATCH --job-name=founder_ase_proc
#SBATCH -A SEYEDAM_LAB
#SBATCH --nodes=1
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --mem=128G
#SBATCH --time=06:00:00
#SBATCH --array=1

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

tissues=("Heart")
tissue="${tissues[$SLURM_ARRAY_TASK_ID-1]}"

: "${ASE_BASE_PATH:?Set ASE_BASE_PATH to the directory containing igvf plate folders}"
: "${ASE_ANNOTATION_DIR:?Set ASE_ANNOTATION_DIR to annotation CSV directory}"
: "${ASE_OUT_DIR:?Set ASE_OUT_DIR to output directory for generated h5ad files}"

ASE_BARCODE_MAP="${ASE_BARCODE_MAP:-refs/barcode_map_df.csv}"
ASE_GENES_WITH_BIOTYPE="${ASE_GENES_WITH_BIOTYPE:-refs/igvf_gtf_gene_biotypes.tsv}"
ASE_GENE_ID_CHR_MAP="${ASE_GENE_ID_CHR_MAP:-refs/gencode_M36_gene_ID_chromosome_map.tsv}"

if [[ "${ASE_BARCODE_MAP}" != /* ]]; then ASE_BARCODE_MAP="${REPO_ROOT}/${ASE_BARCODE_MAP}"; fi
if [[ "${ASE_GENES_WITH_BIOTYPE}" != /* ]]; then ASE_GENES_WITH_BIOTYPE="${REPO_ROOT}/${ASE_GENES_WITH_BIOTYPE}"; fi
if [[ "${ASE_GENE_ID_CHR_MAP}" != /* ]]; then ASE_GENE_ID_CHR_MAP="${REPO_ROOT}/${ASE_GENE_ID_CHR_MAP}"; fi

PY_SCRIPT="${REPO_ROOT}/scripts/Founder_make_strain_ASE_adata.py"

python "${PY_SCRIPT}" -t "${tissue}" \
  --base-path "${ASE_BASE_PATH}" \
  --barcode-map "${ASE_BARCODE_MAP}" \
  --genes-with-biotype "${ASE_GENES_WITH_BIOTYPE}" \
  --gene-id-chr-map "${ASE_GENE_ID_CHR_MAP}" \
  --annotation-dir "${ASE_ANNOTATION_DIR}" \
  --out-dir "${ASE_OUT_DIR}"
