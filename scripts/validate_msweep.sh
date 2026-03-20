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

mkdir -p ${OUTDIR}

# --- 1. Extract binned reads from index files ---
echo "Extracting binned reads..."
for sample_dir in ${MSWEEP_DIR}/*/; do
    sample=$(basename ${sample_dir})
    for bin_file in ${sample_dir}/*.bin; do
        [ -f "${bin_file}" ] || continue
        sc=$(basename ${bin_file} .bin)
        outfq=${OUTDIR}/${sample}_${sc}.fastq.gz

        if [ -f "${outfq}" ]; then
            continue
        fi

        # find the clean fastq for this sample
        fq=$(find ${FASTQDIR} -name "${sample}.fastq.gz" -o -name "${sample}.fq.gz" 2>/dev/null | head -1)
        if [ -z "${fq}" ]; then
            continue
        fi

        # extract reads using seqtk (bin file has 0-based read indices)
        seqtk subseq ${fq} ${bin_file} | gzip > ${outfq}
    done
done
echo "Extraction done."

# --- 2. Sketch reference genomes per SC with mash ---
echo "Sketching reference genomes..."
MASH_REF_DIR=${OUTDIR}/ref_sketches
mkdir -p ${MASH_REF_DIR}

# build per-SC combined reference fastas
python3 - ${SUBSAMPLE_DIR} ${MASH_REF_DIR} << 'PYEOF'
import sys
from collections import defaultdict

subsample_dir = sys.argv[1]
mash_ref_dir = sys.argv[2]

# load cluster labels
cluster_to_label = {}
with open(f"{subsample_dir}/cluster_labels.tsv") as f:
    next(f)
    for line in f:
        parts = line.strip().split("\t")
        cluster_to_label[parts[0]] = f"SC_{parts[0]}_ST{parts[1]}"

# load subsampled genomes
sc_fastas = defaultdict(list)
with open(f"{subsample_dir}/subsampled_genomes.tsv") as f:
    for line in f:
        genome, cluster, fasta = line.strip().split("\t")
        label = cluster_to_label.get(cluster, f"SC_{cluster}")
        sc_fastas[label].append(f"{subsample_dir}/fastas/{genome}.fna")

# write per-SC fasta lists
for sc, fastas in sc_fastas.items():
    with open(f"{mash_ref_dir}/{sc}_refs.txt", "w") as out:
        for f in fastas:
            out.write(f"{f}\n")

print(f"Created ref lists for {len(sc_fastas)} SCs")
PYEOF

# sketch each SC's references
for ref_list in ${MASH_REF_DIR}/*_refs.txt; do
    sc=$(basename ${ref_list} _refs.txt)
    if [ ! -f "${MASH_REF_DIR}/${sc}.msh" ]; then
        mash sketch -l ${ref_list} -o ${MASH_REF_DIR}/${sc} -s 10000 2>/dev/null
    fi
done
echo "Reference sketching done."

# --- 3. Compare binned reads to their assigned SC references ---
echo ""
echo "=== Validation Results ==="
echo "sample	SC	mash_dist	p_value	num_reads	verdict"

for binned_fq in ${OUTDIR}/*_SC_*.fastq.gz; do
    [ -f "${binned_fq}" ] || continue

    filename=$(basename ${binned_fq} .fastq.gz)
    # parse sample and SC from filename (e.g., CRR766456_SC_1_ST2)
    sample=$(echo ${filename} | sed 's/_SC_.*//')
    sc=$(echo ${filename} | sed 's/^[^_]*_//' | sed 's/^//')

    # check read count
    num_reads=$(zcat ${binned_fq} | awk 'END{print NR/4}')
    if [ "${num_reads}" -lt 10 ]; then
        echo "${sample}	${sc}	NA	NA	${num_reads}	too_few_reads"
        continue
    fi

    # sketch the binned reads
    mash sketch ${binned_fq} -o ${OUTDIR}/${filename} -s 10000 -r 2>/dev/null

    # find matching reference sketch
    ref_sketch=${MASH_REF_DIR}/${sc}.msh
    if [ ! -f "${ref_sketch}" ]; then
        echo "${sample}	${sc}	NA	NA	${num_reads}	no_ref_sketch"
        continue
    fi

    # compute mash distance
    result=$(mash dist ${ref_sketch} ${OUTDIR}/${filename}.msh 2>/dev/null | sort -t$'\t' -k3 -n | head -1)
    dist=$(echo "${result}" | cut -f3)
    pval=$(echo "${result}" | cut -f4)

    # verdict: dist < 0.05 = high confidence, 0.05-0.1 = moderate, > 0.1 = likely false positive
    if (( $(echo "${dist} < 0.05" | bc -l) )); then
        verdict="PASS"
    elif (( $(echo "${dist} < 0.1" | bc -l) )); then
        verdict="MODERATE"
    else
        verdict="FAIL"
    fi

    echo "${sample}	${sc}	${dist}	${pval}	${num_reads}	${verdict}"
done | tee ${OUTDIR}/validation_summary.tsv

echo ""
echo "Done. Results in ${OUTDIR}/validation_summary.tsv"
