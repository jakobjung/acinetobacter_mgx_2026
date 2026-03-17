#!/bin/bash
#SBATCH --job-name=subsample_st
#SBATCH --partition=cpu
#SBATCH --qos=long
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --time=1-00:00:00
#SBATCH --output=subsample_st_%j.log

set -euo pipefail

source $(conda info --base)/etc/profile.d/conda.sh
conda activate abaum

PROJDIR=/vol/projects/jjung/acinetobacter_mgx_2026
POPPUNK_DIR=${PROJDIR}/analysis/poppunk/Acinetobacter_baumannii_v1_refs
GENOME_DIR=${PROJDIR}/data/reference_sequences/ab_genomes/ncbi_dataset/ncbi_dataset/data
OUTDIR=${PROJDIR}/analysis/subsample
THREADS=${SLURM_CPUS_PER_TASK:-8}

CLUSTER_CSV=${POPPUNK_DIR}/Acinetobacter_baumannii_v1_refs_clusters.csv
REFS_FILE=${POPPUNK_DIR}/Acinetobacter_baumannii_v1_refs.refs

mkdir -p ${OUTDIR}/fastas ${OUTDIR}/mlst

# --- 1. Map PopPUNK refs to available FASTAs ---
echo "Mapping PopPUNK refs to downloaded FASTAs..."
> ${OUTDIR}/ref_to_fasta.tsv
while read acc; do
    # find matching FASTA (accession may have version suffix in folder name)
    fasta=$(find ${GENOME_DIR}/${acc}* -name "*.fna" 2>/dev/null | head -1)
    if [ -n "${fasta}" ]; then
        echo -e "${acc}\t${fasta}" >> ${OUTDIR}/ref_to_fasta.tsv
    fi
done < ${REFS_FILE}
echo "Mapped $(wc -l < ${OUTDIR}/ref_to_fasta.tsv) / $(wc -l < ${REFS_FILE}) refs to FASTAs"

# --- 2. Subsample max 30 genomes per SC ---
echo "Subsampling max 30 per sequence cluster..."
> ${OUTDIR}/subsampled_genomes.tsv

# get cluster assignments for refs that have FASTAs
awk -F',' 'NR > 1' ${CLUSTER_CSV} | while IFS=',' read -r genome cluster; do
    fasta=$(grep "^${genome}	" ${OUTDIR}/ref_to_fasta.tsv | cut -f2)
    if [ -n "${fasta}" ]; then
        echo -e "${genome}\t${cluster}\t${fasta}"
    fi
done > ${OUTDIR}/genome_cluster_fasta.tsv

# subsample: take up to 30 per cluster
awk -F'\t' '{
    cluster=$2
    count[cluster]++
    if (count[cluster] <= 30) print
}' ${OUTDIR}/genome_cluster_fasta.tsv > ${OUTDIR}/subsampled_genomes.tsv

TOTAL=$(wc -l < ${OUTDIR}/subsampled_genomes.tsv)
CLUSTERS=$(cut -f2 ${OUTDIR}/subsampled_genomes.tsv | sort -u | wc -l)
echo "Subsampled: ${TOTAL} genomes across ${CLUSTERS} clusters"

# --- 3. Copy subsampled FASTAs ---
echo "Copying subsampled FASTAs..."
while IFS=$'\t' read -r genome cluster fasta; do
    cp ${fasta} ${OUTDIR}/fastas/${genome}.fna
done < ${OUTDIR}/subsampled_genomes.tsv

# --- 4. Run mlst to determine sequence types ---
echo "Running mlst on subsampled genomes..."
mlst --threads ${THREADS} ${OUTDIR}/fastas/*.fna > ${OUTDIR}/mlst/mlst_results.tsv
echo "MLST done."

# --- 5. Create SC labels (cluster -> dominant ST) ---
echo "Labelling clusters by dominant ST..."
python3 << 'PYEOF'
import csv
from collections import defaultdict, Counter

# load mlst results
mlst = {}
with open("${OUTDIR}/mlst/mlst_results.tsv") as f:
    for line in f:
        parts = line.strip().split("\t")
        genome = parts[0].split("/")[-1].replace(".fna", "")
        st = parts[2]  # ST column
        mlst[genome] = st

# load subsampled genomes with cluster assignments
cluster_sts = defaultdict(list)
with open("${OUTDIR}/subsampled_genomes.tsv") as f:
    for line in f:
        genome, cluster, fasta = line.strip().split("\t")
        if genome in mlst:
            cluster_sts[cluster].append(mlst[genome])

# determine dominant ST per cluster
with open("${OUTDIR}/cluster_labels.tsv", "w") as out:
    out.write("cluster\tdominant_ST\tcount\ttotal\n")
    for cluster, sts in sorted(cluster_sts.items()):
        counter = Counter(sts)
        dominant_st, count = counter.most_common(1)[0]
        label = f"SC_{cluster}_ST{dominant_st}"
        out.write(f"{cluster}\t{dominant_st}\t{count}\t{len(sts)}\n")

print("Cluster labels written.")
PYEOF

# --- 6. Create Themisto-compatible label file ---
echo "Creating Themisto label mapping..."
python3 << 'PYEOF'
import csv
from collections import defaultdict, Counter

# load cluster labels
cluster_to_label = {}
with open("${OUTDIR}/cluster_labels.tsv") as f:
    next(f)  # skip header
    for line in f:
        cluster, st, count, total = line.strip().split("\t")
        cluster_to_label[cluster] = f"SC_{cluster}_ST{st}"

# load subsampled genomes
with open("${OUTDIR}/themisto_input.tsv", "w") as fasta_list, \
     open("${OUTDIR}/themisto_labels.tsv", "w") as label_list:
    with open("${OUTDIR}/subsampled_genomes.tsv") as f:
        for line in f:
            genome, cluster, fasta_path = line.strip().split("\t")
            label = cluster_to_label.get(cluster, cluster)
            fasta_list.write(f"${OUTDIR}/fastas/{genome}.fna\n")
            label_list.write(f"{label}\n")

print("Themisto input files created.")
PYEOF

echo ""
echo "=== Summary ==="
echo "Subsampled genomes: $(wc -l < ${OUTDIR}/subsampled_genomes.tsv)"
echo "Sequence clusters: $(cut -f2 ${OUTDIR}/subsampled_genomes.tsv | sort -u | wc -l)"
echo "Cluster labels: ${OUTDIR}/cluster_labels.tsv"
echo "Themisto input: ${OUTDIR}/themisto_input.tsv"
echo "Themisto labels: ${OUTDIR}/themisto_labels.tsv"
echo "Done."
