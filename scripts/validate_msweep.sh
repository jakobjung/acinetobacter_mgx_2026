#!/bin/bash
#SBATCH --job-name=validate_msweep
#SBATCH --partition=cpu
#SBATCH --qos=long
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --time=1-00:00:00
#SBATCH --output=validate_msweep_%j.log

set -euo pipefail

source $(conda info --base)/etc/profile.d/conda.sh
conda activate abaum

PROJDIR=/vol/projects/jjung/acinetobacter_mgx_2026
MSWEEP_DIR=${PROJDIR}/analysis/msweep_mgems
SUBSAMPLE_DIR=${PROJDIR}/analysis/subsample
FASTQDIR=${PROJDIR}/data/fastq_clean
OUTDIR=${PROJDIR}/analysis/validation
THREADS=${SLURM_CPUS_PER_TASK:-8}

mkdir -p ${OUTDIR}/ref_sketches ${OUTDIR}/bin_sketches

# --- 1. Sketch reference genomes per SC ---
echo "Sketching reference genomes per SC..."
python3 - ${SUBSAMPLE_DIR} ${OUTDIR}/ref_sketches << 'PYEOF'
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

for ref_list in ${OUTDIR}/ref_sketches/*_refs.txt; do
    sc=$(basename ${ref_list} _refs.txt)
    if [ ! -f "${OUTDIR}/ref_sketches/${sc}.msh" ]; then
        mash sketch -l ${ref_list} -o ${OUTDIR}/ref_sketches/${sc} -s 10000 2>/dev/null
    fi
done
echo "Reference sketching done."

# --- 2. Extract binned reads and sketch them ---
echo "Extracting and sketching binned reads..."
for sample_dir in ${MSWEEP_DIR}/*/; do
    sample=$(basename ${sample_dir})

    for bin_file in ${sample_dir}/*.bin; do
        [ -f "${bin_file}" ] || continue
        sc=$(basename ${bin_file} .bin)

        # skip if already sketched
        if [ -f "${OUTDIR}/bin_sketches/${sample}_${sc}.msh" ]; then
            continue
        fi

        # find the clean fastq
        fq=$(find ${FASTQDIR} -name "${sample}.fastq.gz" -o -name "${sample}.fq.gz" 2>/dev/null | head -1)
        [ -z "${fq}" ] && continue

        # extract reads by index using python (fast)
        python3 - ${fq} ${bin_file} ${OUTDIR}/bin_sketches/${sample}_${sc}.fastq << 'PYEOF'
import sys, gzip

fq_path = sys.argv[1]
bin_path = sys.argv[2]
out_path = sys.argv[3]

# load indices (assume 1-based from mSWEEP)
indices = set()
with open(bin_path) as f:
    for line in f:
        indices.add(int(line.strip()))

# extract reads
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

        # sketch the extracted reads
        num_reads=$(wc -l < ${OUTDIR}/bin_sketches/${sample}_${sc}.fastq)
        num_reads=$((num_reads / 4))

        if [ "${num_reads}" -gt 10 ]; then
            mash sketch ${OUTDIR}/bin_sketches/${sample}_${sc}.fastq \
                -o ${OUTDIR}/bin_sketches/${sample}_${sc} -s 10000 -r 2>/dev/null
        fi

        # cleanup fastq to save space
        rm -f ${OUTDIR}/bin_sketches/${sample}_${sc}.fastq
    done
done
echo "Binned read sketching done."

# --- 3. Compare binned reads to their SC references ---
echo ""
echo "=== Validation Results ==="
printf "sample\tSC\tabundance\tmash_dist\tp_value\tverdict\n" | tee ${OUTDIR}/validation_summary.tsv

for sketch in ${OUTDIR}/bin_sketches/*.msh; do
    [ -f "${sketch}" ] || continue

    filename=$(basename ${sketch} .msh)
    sample=$(echo ${filename} | sed 's/_SC_.*//')
    sc=$(echo ${filename} | sed "s/^${sample}_//")

    # get abundance from mSWEEP
    abun=$(awk -F'\t' -v sc="${sc}" '!/^#/ && $1==sc {print $2}' ${MSWEEP_DIR}/${sample}/msweep_abundances.txt)

    # find matching reference sketch
    ref_sketch=${OUTDIR}/ref_sketches/${sc}.msh
    if [ ! -f "${ref_sketch}" ]; then
        printf "${sample}\t${sc}\t${abun}\tNA\tNA\tno_ref\n" | tee -a ${OUTDIR}/validation_summary.tsv
        continue
    fi

    # compute mash distance (take closest reference)
    result=$(mash dist ${ref_sketch} ${sketch} 2>/dev/null | sort -t$'\t' -k3 -n | head -1)
    dist=$(echo "${result}" | cut -f3)
    pval=$(echo "${result}" | cut -f4)

    # verdict
    if (( $(echo "${dist} < 0.05" | bc -l) )); then
        verdict="PASS"
    elif (( $(echo "${dist} < 0.1" | bc -l) )); then
        verdict="MODERATE"
    else
        verdict="FAIL"
    fi

    printf "${sample}\t${sc}\t${abun}\t${dist}\t${pval}\t${verdict}\n" | tee -a ${OUTDIR}/validation_summary.tsv
done

echo ""
echo "=== Summary ==="
echo "PASS (high confidence):"
grep "PASS" ${OUTDIR}/validation_summary.tsv | sort -t$'\t' -k3 -rn
echo ""
echo "FAIL (likely false positive):"
grep "FAIL" ${OUTDIR}/validation_summary.tsv | sort -t$'\t' -k3 -rn
echo ""
echo "Done."
