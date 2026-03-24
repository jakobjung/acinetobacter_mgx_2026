#!/bin/bash
#SBATCH --job-name=crab_vap_test
#SBATCH --partition=cpu
#SBATCH --qos=long
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --time=2-00:00:00
#SBATCH --output=crab_vap_test_%j.log

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
THREADS=${SLURM_CPUS_PER_TASK:-16}

mkdir -p ${CLEANDIR} ${OUTDIR}

for R1 in ${RAWDIR}/*_1.fastq.gz; do
    sample=$(basename ${R1} _1.fastq.gz)
    R2=${RAWDIR}/${sample}_2.fastq.gz
    SAMPLE_DIR=${OUTDIR}/${sample}

    # skip if already done
    if [ -f "${SAMPLE_DIR}/msweep_abundances.txt" ]; then
        echo "Skipping ${sample} (already processed)"
        continue
    fi

    echo "========================================"
    echo "=== Processing ${sample} ==="
    echo "========================================"

    mkdir -p ${SAMPLE_DIR}/tmp

    # --- 1. Human read removal (PE mode) ---
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

        # bowtie2 --un-conc-gz names files with _1 and _2 suffix via %
        echo "  Human removal done:"
        raw=$(zcat ${R1} | awk 'END{print NR/4}')
        clean=$(zcat ${CLEANDIR}/${sample}_1.fastq.gz | awk 'END{print NR/4}')
        echo "    ${sample}: ${raw} raw PE reads -> ${clean} after human removal"
    else
        echo "  [1/3] Clean reads exist, skipping human removal"
    fi

    # --- 2. Themisto pseudoalign (PE: separate R1 and R2) ---
    echo "  [2/3] Pseudoaligning..."
    themisto pseudoalign \
        -i ${THEMISTO_IDX} \
        -q ${CLEANDIR}/${sample}_1.fastq.gz \
        --rc \
        --temp-dir ${SAMPLE_DIR}/tmp \
        --n-threads ${THREADS} \
        --sort-output \
        -o ${SAMPLE_DIR}/pseudoalignments_R1.txt

    themisto pseudoalign \
        -i ${THEMISTO_IDX} \
        -q ${CLEANDIR}/${sample}_2.fastq.gz \
        --rc \
        --temp-dir ${SAMPLE_DIR}/tmp \
        --n-threads ${THREADS} \
        --sort-output \
        -o ${SAMPLE_DIR}/pseudoalignments_R2.txt

    # --- 3. mSWEEP abundance estimation + read binning ---
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

    # --- Summary ---
    echo ""
    echo "  Results for ${sample}:"
    echo "  Top SCs (abundance > 1%):"
    awk -F'\t' '!/^#/ && $2 > 0.01 {printf "    %s\t%.4f\n", $1, $2}' ${SAMPLE_DIR}/msweep_abundances.txt | sort -t$'\t' -k2 -rn
    echo ""

    # cleanup temp
    rm -rf ${SAMPLE_DIR}/tmp
done

# --- Final summary ---
echo "========================================"
echo "=== SUMMARY: All samples ==="
echo "========================================"
for sample_dir in ${OUTDIR}/*/; do
    sample=$(basename ${sample_dir})
    if [ -f "${sample_dir}/msweep_abundances.txt" ]; then
        num_sc=$(awk -F'\t' '!/^#/ && $2 > 0.01 {count++} END {print count+0}' ${sample_dir}/msweep_abundances.txt)
        top=$(awk -F'\t' '!/^#/ && $2 > 0.01 {printf "%s(%.1f%%) ", $1, $2*100}' ${sample_dir}/msweep_abundances.txt | head -c 200)
        aligned=$(grep "num_aligned" ${sample_dir}/msweep_abundances.txt | cut -f2)
        total=$(grep "num_reads" ${sample_dir}/msweep_abundances.txt | cut -f2)
        echo "${sample}: ${num_sc} SCs >1% | aligned ${aligned}/${total} | ${top}"
    fi
done

echo ""
echo "Done. Results in ${OUTDIR}/"
