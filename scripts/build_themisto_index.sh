#!/bin/bash
#SBATCH --job-name=themisto_idx
#SBATCH --partition=cpu
#SBATCH --qos=long
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --time=1-00:00:00
#SBATCH --output=themisto_idx_%j.log

set -euo pipefail

source $(conda info --base)/etc/profile.d/conda.sh
conda activate abaum

PROJDIR=/vol/projects/jjung/acinetobacter_mgx_2026
SUBSAMPLE_DIR=${PROJDIR}/analysis/subsample
OUTDIR=${PROJDIR}/analysis/themisto
THREADS=${SLURM_CPUS_PER_TASK:-16}

mkdir -p ${OUTDIR}

# --- 1. Build Themisto index ---
echo "Building Themisto index from $(wc -l < ${SUBSAMPLE_DIR}/themisto_input.tsv) genomes..."
themisto build \
    --input-file ${SUBSAMPLE_DIR}/themisto_input.tsv \
    --index-prefix ${OUTDIR}/ab_index \
    --temp-dir ${OUTDIR}/tmp \
    --n-threads ${THREADS}

echo "Index built."

# --- 2. Verify index ---
echo "Verifying index files..."
ls -lh ${OUTDIR}/ab_index*

echo "Done. Index prefix: ${OUTDIR}/ab_index"
