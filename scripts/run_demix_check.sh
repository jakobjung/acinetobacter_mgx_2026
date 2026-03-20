#!/bin/bash
#SBATCH --job-name=demix_check
#SBATCH --partition=cpu
#SBATCH --qos=long
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --time=2-00:00:00
#SBATCH --output=demix_check_%j.log

set -euo pipefail

source $(conda info --base)/etc/profile.d/conda.sh
conda activate abaum

PROJDIR=/vol/projects/jjung/acinetobacter_mgx_2026
SUBSAMPLE_DIR=${PROJDIR}/analysis/subsample
MSWEEP_DIR=${PROJDIR}/analysis/msweep_mgems
OUTDIR=${PROJDIR}/analysis/demix_check
DEMIX_CHECK=${PROJDIR}/tools/demix_check/demix_check.py
THREADS=${SLURM_CPUS_PER_TASK:-16}

mkdir -p ${OUTDIR}/ref_set

# --- 1. Create ref_info.tsv for demix_check ---
echo "Creating reference info file..."
python3 - ${SUBSAMPLE_DIR} ${OUTDIR}/ref_set << 'PYEOF'
import sys

subsample_dir = sys.argv[1]
ref_set_dir = sys.argv[2]

# load subsampled genomes (genome, cluster, fasta_path)
genomes = []
with open(f"{subsample_dir}/subsampled_genomes.tsv") as f:
    for line in f:
        genome, cluster, fasta = line.strip().split("\t")
        genomes.append((genome, cluster, fasta))

# load cluster labels
cluster_to_label = {}
with open(f"{subsample_dir}/cluster_labels.tsv") as f:
    next(f)
    for line in f:
        parts = line.strip().split("\t")
        cluster = parts[0]
        st = parts[1]
        cluster_to_label[cluster] = f"SC_{cluster}_ST{st}"

# write ref_info.tsv
with open(f"{ref_set_dir}/ref_info.tsv", "w") as out:
    out.write("id\tcluster\tassembly\n")
    for genome, cluster, fasta in genomes:
        label = cluster_to_label.get(cluster, f"SC_{cluster}")
        out.write(f"{genome}\t{label}\t{subsample_dir}/fastas/{genome}.fna\n")

print(f"Wrote ref_info.tsv with {len(genomes)} genomes")
PYEOF

# --- 2. Run demix_check setup (compute mash distances and thresholds) ---
echo "Running demix_check setup..."
python3 ${DEMIX_CHECK} \
    --mode_setup \
    --ref ${OUTDIR}/ref_set \
    --threads ${THREADS} \
    --sketch_size 10000

# --- 3. Run demix_check on each sample ---
echo "Validating mSWEEP assignments..."
for sample_dir in ${MSWEEP_DIR}/*/; do
    sample=$(basename ${sample_dir})

    # skip if no bins
    if [ ! -d "${sample_dir}/bins" ] || [ -z "$(ls -A ${sample_dir}/bins/ 2>/dev/null)" ]; then
        continue
    fi

    # skip if already checked
    if [ -f "${OUTDIR}/${sample}/clu_score.tsv" ]; then
        echo "Skipping ${sample} (already checked)"
        continue
    fi

    echo "Checking ${sample}..."
    mkdir -p ${OUTDIR}/${sample}

    python3 ${DEMIX_CHECK} \
        --mode_check \
        --ref ${OUTDIR}/ref_set \
        --binned_reads_dir ${sample_dir}/bins \
        --msweep_abun ${sample_dir}/msweep_abundances.txt \
        --out_dir ${OUTDIR}/${sample} \
        --threads ${THREADS} \
        --min_abun 0.01 \
        --sketch_size 10000 \
        2>> ${OUTDIR}/${sample}/demix_check.log || echo "  ${sample} failed"
done

# --- 4. Summary ---
echo ""
echo "=== demix_check results ==="
echo "Score 1 = high confidence, Score 4 = likely false positive"
echo ""
for score_file in ${OUTDIR}/*/clu_score.tsv; do
    sample=$(basename $(dirname ${score_file}))
    echo "=== ${sample} ==="
    cat ${score_file}
    echo ""
done | grep -E "===|SC_.*ST[0-9]" | head -100

echo "Done. Results in ${OUTDIR}/"
