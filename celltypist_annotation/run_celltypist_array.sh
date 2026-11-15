#!/bin/bash
#SBATCH --job-name=celltypist
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail

# Expected TSV columns, no header:
# sample_id    h5_path    qc_barcodes_path
#
# Submit:
#   sbatch --array=1-N run_celltypist_array.sh samples.tsv

MANIFEST="${1:?Usage: sbatch --array=1-N run_celltypist_array.sh samples.tsv}"

mkdir -p logs results

line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MANIFEST")
sample_id=$(printf '%s\n' "$line" | cut -f1)
h5_path=$(printf '%s\n' "$line" | cut -f2)
barcode_path=$(printf '%s\n' "$line" | cut -f3)


export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

python sub_celltypist.py \
  --sample "$sample_id" \
  --h5 "$h5_path" \
  --barcodes "$barcode_path" \
  --output "results/${sample_id}.celltypist.tsv" \
  --model "Immune_All_Low.pkl"
