#!/bin/bash
#SBATCH --job-name=rm_human
#SBATCH --partition=cpu
#SBATCH --qos=long
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --time=2-00:00:00
#SBATCH --output=rm_human_%j.log

set -euo pipefail

source $(conda info --base)/etc/profile.d/conda.sh
conda activate abaum

PROJDIR=/vol/projects/jjung/acinetobacter_mgx_2026
FASTQDIR=${PROJDIR}/data/fastq
OUTDIR=${PROJDIR}/data/fastq_clean
HG38_IDX=${PROJDIR}/data/reference_sequences/hg38/hg38
THREADS=${SLURM_CPUS_PER_TASK:-8}

mkdir -p ${OUTDIR} ${PROJDIR}/data/reference_sequences/hg38

# --- 1. Download and index hg38 if not present ---
if [ ! -f "${HG38_IDX}.1.bt2" ]; then
    echo "Downloading hg38 bowtie2 index..."
    wget -q https://genome-idx.s3.amazonaws.com/bt2/GRCh38_noalt_as.zip \
        -O ${PROJDIR}/data/reference_sequences/hg38/hg38.zip
    unzip -o ${PROJDIR}/data/reference_sequences/hg38/hg38.zip \
        -d ${PROJDIR}/data/reference_sequences/hg38/
    rm ${PROJDIR}/data/reference_sequences/hg38/hg38.zip
    # rename to hg38 prefix
    for f in ${PROJDIR}/data/reference_sequences/hg38/GRCh38_noalt_as.*; do
        mv "$f" "${f/GRCh38_noalt_as/hg38}"
    done
    echo "hg38 index ready."
else
    echo "hg38 index already exists."
fi

# --- 2. Remove human reads from each sample ---
echo "Removing human reads..."
for CRA in CRA010962 CRA010735 CRA010639; do
    mkdir -p ${OUTDIR}/${CRA}
    for fq in ${FASTQDIR}/${CRA}/*/*.{fastq.gz,fq.gz}; do
        [ -f "${fq}" ] || continue
        sample=$(basename ${fq} | sed 's/\.fastq\.gz$//; s/\.fq\.gz$//')

        if [ -f "${OUTDIR}/${CRA}/${sample}.fastq.gz" ]; then
            echo "Skipping ${sample} (already processed)"
            continue
        fi

        echo "Processing ${CRA}/${sample}..."
        bowtie2 \
            -x ${HG38_IDX} \
            -U ${fq} \
            --threads ${THREADS} \
            --very-sensitive \
            --un-gz ${OUTDIR}/${CRA}/${sample}.fastq.gz \
            -S /dev/null \
            2>> ${OUTDIR}/${CRA}/${sample}.bowtie2.log

        # count reads before/after
        echo "${sample}: $(zcat ${fq} | awk 'END{print NR/4}') raw, $(zcat ${OUTDIR}/${CRA}/${sample}.fastq.gz | awk 'END{print NR/4}') after human removal"
    done
done

echo "Done. Clean reads in ${OUTDIR}/"
