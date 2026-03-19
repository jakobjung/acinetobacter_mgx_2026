#!/bin/bash
#SBATCH --job-name=msweep_mgems
#SBATCH --partition=cpu
#SBATCH --qos=long
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --time=2-00:00:00
#SBATCH --output=msweep_mgems_%j.log

set -euo pipefail

source $(conda info --base)/etc/profile.d/conda.sh
conda activate abaum

PROJDIR=/vol/projects/jjung/acinetobacter_mgx_2026
FASTQDIR=${PROJDIR}/data/fastq_clean
THEMISTO_IDX=${PROJDIR}/analysis/themisto/ab_index
LABELS=${PROJDIR}/analysis/subsample/themisto_labels.tsv
OUTDIR=${PROJDIR}/analysis/msweep_mgems
THREADS=${SLURM_CPUS_PER_TASK:-16}

mkdir -p ${OUTDIR}

# --- 1. Pseudoalign, estimate abundances, and bin reads per sample ---
echo "Running Themisto + mSWEEP + mGEMS pipeline..."
for CRA in CRA010962 CRA010735 CRA010639; do
    for fq in ${FASTQDIR}/${CRA}/*.{fastq.gz,fq.gz}; do
        [ -f "${fq}" ] || continue
        sample=$(basename ${fq} | sed 's/\.fastq\.gz$//; s/\.fq\.gz$//')
        SAMPLE_DIR=${OUTDIR}/${sample}
        mkdir -p ${SAMPLE_DIR}/bins

        if [ -f "${SAMPLE_DIR}/msweep_abundances.txt" ]; then
            echo "Skipping ${sample} (already processed)"
            continue
        fi

        echo "=== Processing ${sample} ==="

        # --- Themisto pseudoalign ---
        echo "  Pseudoaligning..."
        themisto pseudoalign \
            --index-prefix ${THEMISTO_IDX} \
            --query-file ${fq} \
            --rc \
            --temp-dir ${SAMPLE_DIR}/tmp \
            --n-threads ${THREADS} \
            --out-file ${SAMPLE_DIR}/pseudoalignments.txt

        # --- mSWEEP abundance estimation ---
        echo "  Estimating abundances with mSWEEP..."
        mSWEEP \
            --themisto-1 ${SAMPLE_DIR}/pseudoalignments.txt \
            -i ${THEMISTO_IDX} \
            --grouping-list ${LABELS} \
            -o ${SAMPLE_DIR}/msweep \
            -t ${THREADS}

        # rename output
        mv ${SAMPLE_DIR}/msweep_abundances.txt ${SAMPLE_DIR}/msweep_abundances.txt 2>/dev/null || true

        # --- mGEMS read binning ---
        echo "  Binning reads with mGEMS..."
        mGEMS extract \
            --themisto-1 ${SAMPLE_DIR}/pseudoalignments.txt \
            -r ${fq} \
            --index ${THEMISTO_IDX} \
            --groups ${LABELS} \
            --probs ${SAMPLE_DIR}/msweep_probs.txt \
            -o ${SAMPLE_DIR}/bins \
            --min-abundance 0.01

        # --- Summary for this sample ---
        echo "  Top SCs in ${sample}:"
        sort -t$'\t' -k2 -rn ${SAMPLE_DIR}/msweep_abundances.txt | head -5

        # cleanup temp
        rm -rf ${SAMPLE_DIR}/tmp

        echo ""
    done
done

# --- 2. Summary across all samples ---
echo "=== Pipeline complete ==="
echo "Results in ${OUTDIR}/"
echo ""
echo "Samples with detected A. baumannii SCs (abundance >= 1%):"
for dir in ${OUTDIR}/*/; do
    sample=$(basename ${dir})
    if [ -f "${dir}/msweep_abundances.txt" ]; then
        detected=$(awk -F'\t' '$2 >= 0.01 {print $1}' ${dir}/msweep_abundances.txt | tr '\n' ', ')
        if [ -n "${detected}" ]; then
            echo "  ${sample}: ${detected}"
        fi
    fi
done

echo "Done."
