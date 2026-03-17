#!/bin/bash
#SBATCH --job-name=poppunk_ab
#SBATCH --partition=cpu
#SBATCH --qos=long
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --time=2-00:00:00
#SBATCH --output=poppunk_ab_%j.log

set -euo pipefail

# activate poppunk environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate poppunk

PROJDIR=/vol/projects/jjung/acinetobacter_mgx_2026
GENOMEDIR=${PROJDIR}/data/reference_sequences/ab_genomes/ncbi_dataset
OUTDIR=${PROJDIR}/analysis/poppunk
THREADS=${SLURM_CPUS_PER_TASK:-16}

mkdir -p ${OUTDIR}

# --- 1. Create PopPUNK input file (tab-separated: name\tpath) ---
INPUT_LIST=${OUTDIR}/poppunk_input.tsv

if [ -f "${INPUT_LIST}" ]; then
    echo "Input list exists: $(wc -l < ${INPUT_LIST}) genomes"
else
    echo "Creating input list..."
    find ${GENOMEDIR} -name "*.fna" | awk -F'/' '{print $(NF-1)"\t"$0}' > ${INPUT_LIST}
    echo "Found $(wc -l < ${INPUT_LIST}) genomes"
fi

# --- 2. Fit model with DBSCAN (database already created) ---
echo "Fitting PopPUNK DBSCAN model..."
poppunk --fit-model dbscan \
    --ref-db ${OUTDIR}/ab_db \
    --output ${OUTDIR}/ab_clusters_dbscan \
    --no-plot \
    --threads ${THREADS}

# --- 4. Summary ---
CLUSTER_FILE=${OUTDIR}/ab_clusters_dbscan/ab_clusters_dbscan_clusters.csv
echo "Clustering complete."
echo "Total genomes: $(tail -n +2 ${CLUSTER_FILE} | wc -l)"
echo "Total sequence clusters: $(tail -n +2 ${CLUSTER_FILE} | cut -d, -f2 | sort -u | wc -l)"

echo "Done. Cluster assignments: ${CLUSTER_FILE}"
