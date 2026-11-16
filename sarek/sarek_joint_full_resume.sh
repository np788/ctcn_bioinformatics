#!/bin/bash
#SBATCH --job-name=sarek_joint_resume
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=6-00:00:00
#SBATCH --output=logs/sarek_joint_resume_%j.out
#SBATCH --error=logs/sarek_joint_resume_%j.err

set -euo pipefail

BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p \
    "$BASE/logs" \
    "$BASE/work/joint_full" \
    "$BASE/results_joint" \
    "$BASE/singularity"

module load Singularity/4.2.1

export NXF_SINGULARITY_CACHEDIR="$BASE/singularity"
export NXF_OPTS='-Xms1g -Xmx4g'

PARAMS_FILE="${PARAMS_FILE:-$BASE/config/sarek_joint_params.yaml}"

nextflow run nf-core/sarek \
    -r 3.10.0 \
    -profile singularity \
    -c "$BASE/config/slurm.config" \
    -params-file "$PARAMS_FILE" \
    -work-dir "$BASE/work/joint_full" \
    -resume
