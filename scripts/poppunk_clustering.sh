#!/bin/bash
#SBATCH --job-name=poppunk_ab
#SBATCH --partition=cpu
#SBATCH --qos=long
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --time=2-00:00:00
#SBATCH --output=poppunk_ab_%j.log

set -euo pipefail

source $(conda info --base)/etc/profile.d/conda.sh
conda activate poppunk

PROJDIR=/vol/projects/jjung/acinetobacter_mgx_2026
OUTDIR=${PROJDIR}/analysis/poppunk
THREADS=${SLURM_CPUS_PER_TASK:-16}
DB_URL="http://ftp.ebi.ac.uk/pub/databases/pp_dbs/Acinetobacter_baumannii_v1_refs.tar.bz2"
DB_DIR=${OUTDIR}/Acinetobacter_baumannii_v1_refs

mkdir -p ${OUTDIR}

# --- 1. Download pre-built PopPUNK database ---
if [ -d "${DB_DIR}" ]; then
    echo "Pre-built database already exists."
else
    echo "Downloading pre-built A. baumannii PopPUNK database..."
    wget -q ${DB_URL} -O ${OUTDIR}/ab_refs.tar.bz2
    tar -xjf ${OUTDIR}/ab_refs.tar.bz2 -C ${OUTDIR}/
    rm ${OUTDIR}/ab_refs.tar.bz2
    echo "Database downloaded and extracted."
fi

# --- 2. Check for reference strains ---
echo "Checking for key reference strains..."
grep -c "GCA_000963815" ${DB_DIR}/*.csv 2>/dev/null && echo "GC1 (AB5075-UW) found" || echo "GC1 (AB5075-UW) NOT found"
grep -c "GCA_900088705" ${DB_DIR}/*.csv 2>/dev/null && echo "GC2 (BAL062) found" || echo "GC2 (BAL062) NOT found"

# --- 3. Summary of existing clusters ---
CLUSTER_FILE=$(find ${DB_DIR} -name "*_clusters.csv" | head -1)
if [ -n "${CLUSTER_FILE}" ]; then
    echo ""
    echo "Cluster summary:"
    echo "Total genomes: $(tail -n +2 ${CLUSTER_FILE} | wc -l)"
    echo "Total sequence clusters: $(tail -n +2 ${CLUSTER_FILE} | cut -d, -f2 | sort -u | wc -l)"
    echo ""
    echo "Top 20 largest clusters:"
    tail -n +2 ${CLUSTER_FILE} | cut -d, -f2 | sort | uniq -c | sort -rn | head -20
fi

echo "Done. Database: ${DB_DIR}"
