#!/bin/bash
#SBATCH --job-name=validate_crab
#SBATCH --partition=cpu
#SBATCH --qos=long
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --time=2-00:00:00
#SBATCH --output=validate_crab_%j.log

set -euo pipefail

source $(conda info --base)/etc/profile.d/conda.sh
conda activate abaum

PROJDIR=/vol/projects/jjung/acinetobacter_mgx_2026
MSWEEP_DIR=${PROJDIR}/analysis/crab_vap_results
SUBSAMPLE_DIR=${PROJDIR}/analysis/subsample
FASTQDIR=${PROJDIR}/data/crab_vap_clean
OUTDIR=${PROJDIR}/analysis/crab_vap_validation
THREADS=${SLURM_CPUS_PER_TASK:-8}

mkdir -p ${OUTDIR}/ref_sketches ${OUTDIR}/bin_sketches

# --- 1. Sketch reference genomes per SC (reuse if already done) ---
PREV_REF=${PROJDIR}/analysis/validation/ref_sketches
if [ -d "${PREV_REF}" ] && [ "$(ls ${PREV_REF}/*.msh 2>/dev/null | wc -l)" -gt 0 ]; then
    echo "Reusing existing reference sketches from ${PREV_REF}"
    ln -sf ${PREV_REF} ${OUTDIR}/ref_sketches_link
    REF_SKETCH_DIR=${PREV_REF}
else
    echo "Sketching reference genomes per SC..."
    REF_SKETCH_DIR=${OUTDIR}/ref_sketches

    python3 - ${SUBSAMPLE_DIR} ${REF_SKETCH_DIR} << 'PYEOF'
import sys
from collections import defaultdict

subsample_dir = sys.argv[1]
ref_dir = sys.argv[2]

cluster_to_label = {}
with open(f"{subsample_dir}/cluster_labels.tsv") as f:
    next(f)
    for line in f:
        parts = line.strip().split("\t")
        cluster_to_label[parts[0]] = f"SC_{parts[0]}_ST{parts[1]}"

sc_fastas = defaultdict(list)
with open(f"{subsample_dir}/subsampled_genomes.tsv") as f:
    for line in f:
        genome, cluster, fasta = line.strip().split("\t")
        label = cluster_to_label.get(cluster, f"SC_{cluster}")
        sc_fastas[label].append(f"{subsample_dir}/fastas/{genome}.fna")

for sc, fastas in sc_fastas.items():
    with open(f"{ref_dir}/{sc}_refs.txt", "w") as out:
        for f in fastas:
            out.write(f"{f}\n")

print(f"Created ref lists for {len(sc_fastas)} SCs")
PYEOF

    for ref_list in ${REF_SKETCH_DIR}/*_refs.txt; do
        sc=$(basename ${ref_list} _refs.txt)
        if [ ! -f "${REF_SKETCH_DIR}/${sc}.msh" ]; then
            mash sketch -l ${ref_list} -o ${REF_SKETCH_DIR}/${sc} -s 10000 2>/dev/null
        fi
    done
    echo "Reference sketching done."
fi

# --- 2. Extract binned reads, sketch, and validate ---
echo "Extracting and validating binned reads..."
echo ""
printf "sample\tSC\tabundance\tmash_dist\tverdict\n" > ${OUTDIR}/validation_summary.tsv

for sample_dir in ${MSWEEP_DIR}/*/; do
    sample=$(basename ${sample_dir})

    for bin_file in ${sample_dir}/*.bin; do
        [ -f "${bin_file}" ] || continue
        sc=$(basename ${bin_file} .bin)

        # get abundance
        abun=$(awk -F'\t' -v sc="${sc}" '!/^#/ && $1==sc {print $2}' ${sample_dir}/msweep_abundances.txt)

        # skip low abundance
        if (( $(echo "${abun} < 0.01" | bc -l) )); then
            continue
        fi

        # skip if already done
        if [ -f "${OUTDIR}/bin_sketches/${sample}_${sc}.msh" ]; then
            result=$(mash dist ${REF_SKETCH_DIR}/${sc}.msh ${OUTDIR}/bin_sketches/${sample}_${sc}.msh 2>/dev/null | sort -t$'\t' -k3 -n | head -1)
            dist=$(echo "${result}" | cut -f3)
        else
            # find clean fastq (PE: use R1)
            fq=${FASTQDIR}/${sample}_1.fastq.gz
            [ -f "${fq}" ] || continue

            # extract reads by index
            python3 - ${fq} ${bin_file} ${OUTDIR}/bin_sketches/${sample}_${sc}.fastq << 'PYEOF'
import sys, gzip

fq_path = sys.argv[1]
bin_path = sys.argv[2]
out_path = sys.argv[3]

indices = set()
with open(bin_path) as f:
    for line in f:
        indices.add(int(line.strip()))

opener = gzip.open if fq_path.endswith('.gz') else open
with opener(fq_path, 'rt') as fq, open(out_path, 'w') as out:
    read_num = 0
    for i, line in enumerate(fq):
        if i % 4 == 0:
            read_num += 1
            keep = read_num in indices
        if keep:
            out.write(line)
PYEOF

            num_reads=$(awk 'END{print NR/4}' ${OUTDIR}/bin_sketches/${sample}_${sc}.fastq)

            if [ "${num_reads}" -gt 10 ]; then
                mash sketch ${OUTDIR}/bin_sketches/${sample}_${sc}.fastq \
                    -o ${OUTDIR}/bin_sketches/${sample}_${sc} -s 10000 -r 2>/dev/null
            fi

            rm -f ${OUTDIR}/bin_sketches/${sample}_${sc}.fastq

            # compute distance
            ref_sketch=${REF_SKETCH_DIR}/${sc}.msh
            if [ ! -f "${ref_sketch}" ] || [ ! -f "${OUTDIR}/bin_sketches/${sample}_${sc}.msh" ]; then
                printf "${sample}\t${sc}\t${abun}\tNA\tno_sketch\n" >> ${OUTDIR}/validation_summary.tsv
                continue
            fi

            result=$(mash dist ${ref_sketch} ${OUTDIR}/bin_sketches/${sample}_${sc}.msh 2>/dev/null | sort -t$'\t' -k3 -n | head -1)
            dist=$(echo "${result}" | cut -f3)
        fi

        # verdict
        if [ -z "${dist}" ]; then
            verdict="no_dist"
        elif (( $(echo "${dist} < 0.05" | bc -l) )); then
            verdict="PASS"
        elif (( $(echo "${dist} < 0.1" | bc -l) )); then
            verdict="MODERATE"
        else
            verdict="FAIL"
        fi

        printf "${sample}\t${sc}\t${abun}\t${dist}\t${verdict}\n" >> ${OUTDIR}/validation_summary.tsv
    done
done

# --- 3. Summary ---
echo ""
echo "=== VALIDATION SUMMARY ==="
echo ""
echo "--- PASS (high confidence, dist < 0.05) ---"
awk -F'\t' '$5=="PASS"' ${OUTDIR}/validation_summary.tsv | sort -t$'\t' -k3 -rn
echo ""
echo "--- MODERATE (0.05 < dist < 0.1) ---"
awk -F'\t' '$5=="MODERATE"' ${OUTDIR}/validation_summary.tsv | sort -t$'\t' -k3 -rn
echo ""
echo "--- FAIL (dist > 0.1, likely false positive) ---"
awk -F'\t' '$5=="FAIL"' ${OUTDIR}/validation_summary.tsv | wc -l
echo " entries (not shown)"
echo ""
echo "Done. Full results: ${OUTDIR}/validation_summary.tsv"
