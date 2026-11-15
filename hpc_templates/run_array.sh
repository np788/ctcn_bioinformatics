#!/bin/bash
#SBATCH --job-name=analysis
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail

MANIFEST="${1:?Usage: sbatch --array=1-N run_array.sh samples.tsv}"

mkdir -p logs results

line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MANIFEST")
sample_id=$(printf '%s\n' "$line" | cut -f1)
input_path=$(printf '%s\n' "$line" | cut -f2-)

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

outdir="results/${sample_id}"
mkdir -p "$outdir"

echo "[$(date)] sample=${sample_id}"
echo "[$(date)] input=${input_path}"
echo "[$(date)] output=${outdir}"

# Replace
# python analysis.py --input "$input_path" --output "$outdir"

echo "[$(date)] completed ${sample_id}"
