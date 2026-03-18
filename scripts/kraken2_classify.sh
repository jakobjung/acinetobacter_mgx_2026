#!/bin/bash
#SBATCH --job-name=kraken2_ab
#SBATCH --partition=cpu
#SBATCH --qos=long
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --time=2-00:00:00
#SBATCH --output=kraken2_%j.log

set -euo pipefail

source $(conda info --base)/etc/profile.d/conda.sh
conda activate abaum

PROJDIR=/vol/projects/jjung/acinetobacter_mgx_2026
FASTQDIR=${PROJDIR}/data/fastq
KRAKENDB=${PROJDIR}/data/kraken2_standard_db
OUTDIR=${PROJDIR}/analysis/kraken2
THREADS=${SLURM_CPUS_PER_TASK:-16}

mkdir -p ${KRAKENDB} ${OUTDIR}

# --- 1. Download pre-built standard Kraken2 database ---
if [ ! -f "${KRAKENDB}/hash.k2d" ]; then
    echo "Downloading Kraken2 standard database..."
    wget -q https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20240904.tar.gz \
        -O ${KRAKENDB}/k2_standard.tar.gz
    tar -xzf ${KRAKENDB}/k2_standard.tar.gz -C ${KRAKENDB}
    rm ${KRAKENDB}/k2_standard.tar.gz
    echo "Database downloaded."
else
    echo "Kraken2 database already exists."
fi

# --- 2. Run Kraken2 on all samples ---
echo "Running Kraken2 classification..."
for CRA in CRA010962 CRA010735 CRA010639; do
    for fq in ${FASTQDIR}/${CRA}/*/*.fastq.gz; do
        sample=$(basename ${fq} .fastq.gz)
        echo "Processing ${CRA}/${sample}..."

        kraken2 \
            --db ${KRAKENDB} \
            --threads ${THREADS} \
            --output ${OUTDIR}/${sample}.kraken.out \
            --report ${OUTDIR}/${sample}.kraken.report \
            ${fq}
    done
done

# --- 3. Summary: top species across all samples ---
echo ""
echo "=== Top species across all samples ==="
for report in ${OUTDIR}/*.kraken.report; do
    sample=$(basename ${report} .kraken.report)
    top=$(awk -F'\t' '$4 == "S" {gsub(/^ +/, "", $6); print $1"% "$6}' ${report} | head -5)
    echo "${sample}: ${top}"
done

echo "Done. Reports in ${OUTDIR}/"
