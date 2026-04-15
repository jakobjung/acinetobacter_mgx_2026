#!/bin/bash
#SBATCH --job-name=crab_v3
#SBATCH --partition=cpu
#SBATCH --qos=long
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --time=5-00:00:00
#SBATCH --output=crab_vap_v3_%j.log

set -euo pipefail

source $(conda info --base)/etc/profile.d/conda.sh
conda activate abaum

PROJDIR=/vol/projects/jjung/acinetobacter_mgx_2026
CLEANDIR=${PROJDIR}/data/crab_vap_clean
KRAKENDB=${PROJDIR}/data/kraken2_standard_db
THEMISTO_IDX=${PROJDIR}/analysis/themisto_v2/ab_index
LABELS=${PROJDIR}/analysis/subsample_v2/themisto_labels.tsv
OUTDIR=${PROJDIR}/analysis/crab_vap_results_v3
THREADS=${SLURM_CPUS_PER_TASK:-16}

mkdir -p ${OUTDIR}

TOTAL=$(ls ${CLEANDIR}/*_1.fastq.gz 2>/dev/null | wc -l)
COUNT=0

for R1 in ${CLEANDIR}/*_1.fastq.gz; do
    sample=$(basename ${R1} _1.fastq.gz)
    R2=${CLEANDIR}/${sample}_2.fastq.gz
    SAMPLE_DIR=${OUTDIR}/${sample}
    COUNT=$((COUNT + 1))

    # skip if already done
    if [ -f "${SAMPLE_DIR}/msweep_abundances.txt" ]; then
        echo "[${COUNT}/${TOTAL}] Skipping ${sample} (already processed)"
        continue
    fi

    echo "========================================"
    echo "[${COUNT}/${TOTAL}] Processing ${sample}"
    echo "========================================"

    mkdir -p ${SAMPLE_DIR}/tmp

    # --- Round 1: Kraken2 species-level classification ---
    if [ ! -f "${SAMPLE_DIR}/${sample}_abau_R1.fastq.gz" ]; then
        echo "  [1/4] Kraken2 species classification..."
        kraken2 \
            --db ${KRAKENDB} \
            --threads ${THREADS} \
            --paired \
            --output ${SAMPLE_DIR}/kraken2.out \
            --report ${SAMPLE_DIR}/kraken2.report \
            ${R1} ${R2}

        # count A. baumannii reads
        abau_reads=$(awk '$1 == "C" && $3 == 470' ${SAMPLE_DIR}/kraken2.out | wc -l)
        total_reads=$(wc -l < ${SAMPLE_DIR}/kraken2.out)
        echo "  Kraken2: ${abau_reads} / ${total_reads} reads classified as A. baumannii"

        # --- Extract A. baumannii reads (taxid 470 + children) ---
        echo "  Extracting A. baumannii reads..."
        # get read IDs classified as A. baumannii (taxid 470)
        awk '$1 == "C" && $3 == 470 {print $2}' ${SAMPLE_DIR}/kraken2.out > ${SAMPLE_DIR}/abau_readids.txt

        # extract PE reads using seqtk
        seqtk subseq ${R1} ${SAMPLE_DIR}/abau_readids.txt | gzip > ${SAMPLE_DIR}/${sample}_abau_R1.fastq.gz
        seqtk subseq ${R2} ${SAMPLE_DIR}/abau_readids.txt | gzip > ${SAMPLE_DIR}/${sample}_abau_R2.fastq.gz

        # cleanup large kraken output
        rm -f ${SAMPLE_DIR}/kraken2.out ${SAMPLE_DIR}/abau_readids.txt

        abau_extracted=$(zcat ${SAMPLE_DIR}/${sample}_abau_R1.fastq.gz | awk 'END{print NR/4}')
        echo "  Extracted ${abau_extracted} A. baumannii PE reads"
    else
        echo "  [1/4] A. baumannii reads already extracted"
    fi

    # skip if no A. baumannii reads
    abau_count=$(zcat ${SAMPLE_DIR}/${sample}_abau_R1.fastq.gz | awk 'END{print NR/4}')
    if [ "${abau_count}" -lt 100 ]; then
        echo "  Too few A. baumannii reads (${abau_count}), skipping mSWEEP"
        # write empty abundances file
        echo -e "#mSWEEP_version:\t\n#num_reads:\t${abau_count}\n#num_aligned:\t0\n#c_id\tmean_theta" > ${SAMPLE_DIR}/msweep_abundances.txt
        continue
    fi

    # --- Round 2: Themisto pseudoalign (SC-level, A. baumannii reads only) ---
    echo "  [2/4] Pseudoaligning A. baumannii reads..."
    themisto pseudoalign \
        -i ${THEMISTO_IDX} \
        -q ${SAMPLE_DIR}/${sample}_abau_R1.fastq.gz \
        --rc \
        --temp-dir ${SAMPLE_DIR}/tmp \
        --n-threads ${THREADS} \
        --sort-output \
        -o ${SAMPLE_DIR}/pseudoalignments_R1.txt

    themisto pseudoalign \
        -i ${THEMISTO_IDX} \
        -q ${SAMPLE_DIR}/${sample}_abau_R2.fastq.gz \
        --rc \
        --temp-dir ${SAMPLE_DIR}/tmp \
        --n-threads ${THREADS} \
        --sort-output \
        -o ${SAMPLE_DIR}/pseudoalignments_R2.txt

    # --- Round 2: mSWEEP abundance estimation ---
    echo "  [3/4] Running mSWEEP..."
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
    echo "  [4/4] Results for ${sample}:"
    awk -F'\t' '!/^#/ && $2 > 0.01 {printf "    %s\t%.4f\n", $1, $2}' ${SAMPLE_DIR}/msweep_abundances.txt | sort -t$'\t' -k2 -rn
    echo ""

    # cleanup
    rm -rf ${SAMPLE_DIR}/tmp
    rm -f ${SAMPLE_DIR}/pseudoalignments_R1.txt ${SAMPLE_DIR}/pseudoalignments_R2.txt
done

# --- Final summary ---
echo "========================================"
echo "=== FINAL SUMMARY (v3 — species-filtered) ==="
echo "========================================"
for sample_dir in ${OUTDIR}/*/; do
    sample=$(basename ${sample_dir})
    if [ -f "${sample_dir}/msweep_abundances.txt" ]; then
        aligned=$(grep "num_aligned" ${sample_dir}/msweep_abundances.txt | cut -f2)
        total=$(grep "num_reads" ${sample_dir}/msweep_abundances.txt | cut -f2)
        top=$(awk -F'\t' '!/^#/ && $2 > 0.01 {printf "%s(%.1f%%) ", $1, $2*100}' ${sample_dir}/msweep_abundances.txt)
        echo "${sample}	aligned=${aligned}/${total}	${top}"
    fi
done | tee ${OUTDIR}/summary_v3.tsv

echo ""
echo "Done."
