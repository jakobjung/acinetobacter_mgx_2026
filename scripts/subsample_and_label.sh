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

# --- 1. Build lookup of all available FASTAs (fast: single find) ---
echo "Building FASTA lookup..."
find ${GENOME_DIR} -name "*.fna" > ${OUTDIR}/all_fastas.txt
# create mapping: accession_no_version -> fasta_path
awk -F'/' '{
    dir=$(NF-1)
    sub(/\.[0-9]+$/, "", dir)
    print dir"\t"$0
}' ${OUTDIR}/all_fastas.txt > ${OUTDIR}/acc_to_fasta.tsv
echo "Found $(wc -l < ${OUTDIR}/acc_to_fasta.tsv) FASTAs"

# --- 2. Join PopPUNK clusters with available FASTAs, subsample max 30/SC ---
echo "Joining clusters with FASTAs and subsampling..."
python3 - ${OUTDIR} ${CLUSTER_CSV} ${REFS_FILE} ${OUTDIR}/acc_to_fasta.tsv << 'PYEOF'
import sys
from collections import defaultdict

outdir = sys.argv[1]
cluster_csv = sys.argv[2]
refs_file = sys.argv[3]
acc_fasta_tsv = sys.argv[4]

# load fasta lookup: accession (no version) -> path
acc_to_fasta = {}
with open(acc_fasta_tsv) as f:
    for line in f:
        acc, fasta = line.strip().split("\t")
        acc_to_fasta[acc] = fasta

# load refs list
refs = set()
with open(refs_file) as f:
    for line in f:
        refs.add(line.strip())

# load cluster assignments for refs
cluster_genomes = defaultdict(list)
with open(cluster_csv) as f:
    next(f)  # skip header
    for line in f:
        genome, cluster = line.strip().split(",")
        if genome in refs and genome in acc_to_fasta:
            cluster_genomes[cluster].append((genome, acc_to_fasta[genome]))

# subsample max 30 per cluster
with open(f"{outdir}/subsampled_genomes.tsv", "w") as out:
    total = 0
    for cluster, genomes in sorted(cluster_genomes.items()):
        for genome, fasta in genomes[:30]:
            out.write(f"{genome}\t{cluster}\t{fasta}\n")
            total += 1

print(f"Subsampled: {total} genomes across {len(cluster_genomes)} clusters")
PYEOF

echo "$(wc -l < ${OUTDIR}/subsampled_genomes.tsv) genomes subsampled"

# --- 3. Copy subsampled FASTAs ---
echo "Copying subsampled FASTAs..."
cut -f1,3 ${OUTDIR}/subsampled_genomes.tsv | while IFS=$'\t' read -r genome fasta; do
    cp "${fasta}" "${OUTDIR}/fastas/${genome}.fna"
done
echo "Copied $(ls ${OUTDIR}/fastas/*.fna | wc -l) FASTAs"

# --- 4. Run mlst ---
echo "Running mlst on subsampled genomes..."
mlst --threads ${THREADS} ${OUTDIR}/fastas/*.fna > ${OUTDIR}/mlst/mlst_results.tsv
echo "MLST done: $(wc -l < ${OUTDIR}/mlst/mlst_results.tsv) results"

# --- 5. Label clusters by dominant ST and create Themisto files ---
echo "Labelling clusters and creating Themisto input..."
python3 - ${OUTDIR} << 'PYEOF'
import sys
from collections import defaultdict, Counter

outdir = sys.argv[1]

# load mlst results
mlst = {}
with open(f"{outdir}/mlst/mlst_results.tsv") as f:
    for line in f:
        parts = line.strip().split("\t")
        genome = parts[0].split("/")[-1].replace(".fna", "")
        st = parts[2]
        mlst[genome] = st

# load subsampled genomes
cluster_sts = defaultdict(list)
genomes_info = []
with open(f"{outdir}/subsampled_genomes.tsv") as f:
    for line in f:
        genome, cluster, fasta = line.strip().split("\t")
        genomes_info.append((genome, cluster))
        if genome in mlst:
            cluster_sts[cluster].append(mlst[genome])

# dominant ST per cluster
cluster_to_label = {}
with open(f"{outdir}/cluster_labels.tsv", "w") as out:
    out.write("cluster\tdominant_ST\tcount\ttotal\n")
    for cluster, sts in sorted(cluster_sts.items()):
        counter = Counter(sts)
        dominant_st, count = counter.most_common(1)[0]
        cluster_to_label[cluster] = f"SC_{cluster}_ST{dominant_st}"
        out.write(f"{cluster}\t{dominant_st}\t{count}\t{len(sts)}\n")

# Themisto input files
with open(f"{outdir}/themisto_input.tsv", "w") as fasta_list, \
     open(f"{outdir}/themisto_labels.tsv", "w") as label_list:
    for genome, cluster in genomes_info:
        label = cluster_to_label.get(cluster, cluster)
        fasta_list.write(f"{outdir}/fastas/{genome}.fna\n")
        label_list.write(f"{label}\n")

print(f"Created labels for {len(cluster_to_label)} clusters")
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
