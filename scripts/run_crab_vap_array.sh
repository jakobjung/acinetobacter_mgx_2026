#!/bin/bash
#SBATCH --job-name=crab_vap
#SBATCH --partition=cpu
#SBATCH --qos=long
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --time=1-00:00:00
#SBATCH --output=crab_vap_%A_%a.log
#SBATCH --array=1-49

set -euo pipefail

source $(conda info --base)/etc/profile.d/conda.sh
conda activate abaum

PROJDIR=/vol/projects/jjung/acinetobacter_mgx_2026
RAWDIR=${PROJDIR}/data/crab_vap_test
CLEANDIR=${PROJDIR}/data/crab_vap_clean
HG38_IDX=${PROJDIR}/data/reference_sequences/hg38/GRCh38_noalt_as/GRCh38_noalt_as
THEMISTO_IDX=${PROJDIR}/analysis/themisto/ab_index
LABELS=${PROJDIR}/analysis/subsample/themisto_labels.tsv
OUTDIR=${PROJDIR}/analysis/crab_vap_results
THREADS=${SLURM_CPUS_PER_TASK:-8}

mkdir -p ${CLEANDIR} ${OUTDIR}

# get the sample for this array task
R1=$(ls ${RAWDIR}/*_1.fastq.gz | sed -n "${SLURM_ARRAY_TASK_ID}p")
if [ -z "${R1}" ]; then
    echo "No sample for task ${SLURM_ARRAY_TASK_ID}, exiting."
    exit 0
fi

sample=$(basename ${R1} _1.fastq.gz)
R2=${RAWDIR}/${sample}_2.fastq.gz
SAMPLE_DIR=${OUTDIR}/${sample}

# skip if already done
if [ -f "${SAMPLE_DIR}/msweep_abundances.txt" ]; then
    echo "Skipping ${sample} (already processed)"
    exit 0
fi

echo "=== Processing ${sample} (task ${SLURM_ARRAY_TASK_ID}) ==="
mkdir -p ${SAMPLE_DIR}/tmp

# --- 1. Human read removal ---
if [ ! -f "${CLEANDIR}/${sample}_1.fastq.gz" ]; then
    echo "  [1/3] Removing human reads..."
    bowtie2 \
        -x ${HG38_IDX} \
        -1 ${R1} \
        -2 ${R2} \
        --threads ${THREADS} \
        --sensitive \
        --un-conc-gz ${CLEANDIR}/${sample}_%.fastq.gz \
        -S /dev/null \
        2> ${CLEANDIR}/${sample}.bowtie2.log
    echo "  Done."
else
    echo "  [1/3] Clean reads exist."
fi

# --- 2. Themisto pseudoalign ---
echo "  [2/3] Pseudoaligning..."
themisto pseudoalign \
    -i ${THEMISTO_IDX} \
    -q ${CLEANDIR}/${sample}_1.fastq.gz \
    --rc \
    --temp-dir ${SAMPLE_DIR}/tmp \
    --n-threads ${THREADS} \
    -o ${SAMPLE_DIR}/pseudoalignments_R1.txt

themisto pseudoalign \
    -i ${THEMISTO_IDX} \
    -q ${CLEANDIR}/${sample}_2.fastq.gz \
    --rc \
    --temp-dir ${SAMPLE_DIR}/tmp \
    --n-threads ${THREADS} \
    -o ${SAMPLE_DIR}/pseudoalignments_R2.txt

# --- 3. mSWEEP ---
echo "  [3/3] Running mSWEEP..."
mSWEEP \
    -i ${LABELS} \
    --themisto-1 ${SAMPLE_DIR}/pseudoalignments_R1.txt \
    --themisto-2 ${SAMPLE_DIR}/pseudoalignments_R2.txt \
    -o ${SAMPLE_DIR}/msweep \
    -t ${THREADS} \
    --bin-reads \
    --min-abundance 0.01 \
    --write-probs

# Summary
echo "  Top SCs (>1%):"
awk -F'\t' '!/^#/ && $2 > 0.01 {printf "    %s\t%.4f\n", $1, $2}' ${SAMPLE_DIR}/msweep_abundances.txt | sort -t$'\t' -k2 -rn

# Cleanup
rm -rf ${SAMPLE_DIR}/tmp
rm -f ${SAMPLE_DIR}/pseudoalignments_R1.txt ${SAMPLE_DIR}/pseudoalignments_R2.txt

echo "Done: ${sample}"
