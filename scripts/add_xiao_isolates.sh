#!/bin/bash
#SBATCH --job-name=add_xiao
#SBATCH --partition=cpu
#SBATCH --qos=long
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --time=1-00:00:00
#SBATCH --output=add_xiao_%j.log

set -euo pipefail

source $(conda info --base)/etc/profile.d/conda.sh
conda activate abaum

PROJDIR=/vol/projects/jjung/acinetobacter_mgx_2026
OUTDIR=${PROJDIR}/analysis/xiao_isolates
SUBSAMPLE_DIR=${PROJDIR}/analysis/subsample
THREADS=${SLURM_CPUS_PER_TASK:-8}

mkdir -p ${OUTDIR}

# --- 1. Download Xiao et al. isolate genomes (PRJNA679997) ---
if [ ! -d "${OUTDIR}/ncbi_dataset" ]; then
    echo "Downloading Xiao et al. isolate genomes..."
    datasets download genome accession PRJNA679997 \
        --include genome \
        --filename ${OUTDIR}/xiao_isolates.zip
    unzip -o ${OUTDIR}/xiao_isolates.zip -d ${OUTDIR}/ncbi_dataset
    rm ${OUTDIR}/xiao_isolates.zip
    echo "Downloaded $(find ${OUTDIR}/ncbi_dataset -name '*.fna' | wc -l) genomes"
else
    echo "Xiao isolate genomes already downloaded."
fi

# --- 2. MLST already run manually, skip ---
echo "MLST results already available."
echo "Pasteur: $(wc -l < ${OUTDIR}/mlst_pasteur.tsv) isolates"
echo "Oxford: $(wc -l < ${OUTDIR}/mlst_oxford.tsv) isolates"

# --- 3. Create combined reference panel ---
echo ""
echo "Building updated reference panel..."
NEW_SUBSAMPLE=${PROJDIR}/analysis/subsample_v2
mkdir -p ${NEW_SUBSAMPLE}/fastas ${NEW_SUBSAMPLE}/mlst

# copy existing subsampled genomes
cp ${SUBSAMPLE_DIR}/fastas/*.fna ${NEW_SUBSAMPLE}/fastas/
cp ${SUBSAMPLE_DIR}/subsampled_genomes.tsv ${NEW_SUBSAMPLE}/subsampled_genomes.tsv
cp ${SUBSAMPLE_DIR}/cluster_labels.tsv ${NEW_SUBSAMPLE}/cluster_labels.tsv

# add Xiao isolates as new SCs based on gyrB+gpi allelic profile
python3 - ${OUTDIR} ${NEW_SUBSAMPLE} << 'PYEOF'
import sys, re, glob, shutil
from collections import defaultdict

outdir = sys.argv[1]
new_subsample = sys.argv[2]

# load Oxford MLST results — group by gyrB + gpi profile
isolate_profiles = {}
with open(f"{outdir}/mlst_oxford.tsv") as f:
    for line in f:
        parts = line.strip().split("\t")
        genome_path = parts[0]
        # extract accession from directory path (e.g. .../GCA_022815885.1/...)
        path_parts = genome_path.split("/")
        accession = [p for p in path_parts if p.startswith("GCA_")]
        if not accession:
            continue
        accession = accession[0]
        # extract gyrB and gpi alleles
        gyrB = re.search(r'Oxf_gyrB\((\d+)\)', line)
        gpi = re.search(r'Oxf_gpi\((\d+)\)', line)
        if gyrB and gpi:
            profile = f"gyrB{gyrB.group(1)}_gpi{gpi.group(1)}"
        else:
            profile = "unknown"
        isolate_profiles[accession] = profile

# group by profile
profile_groups = defaultdict(list)
for acc, profile in isolate_profiles.items():
    profile_groups[profile].append(acc)

# find existing max cluster number
max_cluster = 0
with open(f"{new_subsample}/cluster_labels.tsv") as f:
    next(f)
    for line in f:
        parts = line.strip().split("\t")
        try:
            c = int(parts[0])
            if c > max_cluster:
                max_cluster = c
        except ValueError:
            pass

# add new SCs
new_clusters = []
with open(f"{new_subsample}/subsampled_genomes.tsv", "a") as sg:
    for profile, accessions in sorted(profile_groups.items()):
        max_cluster += 1
        cluster_id = max_cluster

        for acc in accessions[:30]:
            fasta_paths = glob.glob(f"{outdir}/ncbi_dataset/ncbi_dataset/data/{acc}/*.fna")
            if not fasta_paths:
                continue
            fasta = fasta_paths[0]
            dest = f"{new_subsample}/fastas/XIAO_{acc}.fna"
            shutil.copy(fasta, dest)
            sg.write(f"XIAO_{acc}\t{cluster_id}\t{dest}\n")

        new_clusters.append((cluster_id, profile, len(accessions)))

with open(f"{new_subsample}/cluster_labels.tsv", "a") as cl:
    for cluster_id, profile, count in new_clusters:
        cl.write(f"{cluster_id}\tXIAO_{profile}\t{min(count,30)}\t{min(count,30)}\n")

print(f"Added {len(new_clusters)} new SCs from Xiao et al. isolates")
for c, profile, n in new_clusters:
    print(f"  SC_{c}: {profile} ({n} isolates, max 30 used)")
PYEOF

# --- 4. Create new Themisto input files ---
echo "Creating Themisto input files..."
python3 - ${NEW_SUBSAMPLE} << 'PYEOF'
import sys

subsample_dir = sys.argv[1]

cluster_to_label = {}
with open(f"{subsample_dir}/cluster_labels.tsv") as f:
    next(f)
    for line in f:
        parts = line.strip().split("\t")
        cluster_to_label[parts[0]] = f"SC_{parts[0]}_ST{parts[1]}"

with open(f"{subsample_dir}/themisto_input.txt", "w") as fi, \
     open(f"{subsample_dir}/themisto_labels.tsv", "w") as fl:
    with open(f"{subsample_dir}/subsampled_genomes.tsv") as f:
        for line in f:
            genome, cluster, fasta = line.strip().split("\t")
            label = cluster_to_label.get(cluster, f"SC_{cluster}")
            fi.write(f"{subsample_dir}/fastas/{genome}.fna\n")
            fl.write(f"{label}\n")

total = sum(1 for _ in open(f"{subsample_dir}/themisto_input.txt"))
print(f"Themisto input: {total} genomes")
PYEOF

# --- 5. Build new Themisto index ---
echo "Building new Themisto index..."
THEMISTO_DIR=${PROJDIR}/analysis/themisto_v2
mkdir -p ${THEMISTO_DIR}/tmp

themisto build \
    -k 31 \
    --input-file ${NEW_SUBSAMPLE}/themisto_input.txt \
    --index-prefix ${THEMISTO_DIR}/ab_index \
    --temp-dir ${THEMISTO_DIR}/tmp \
    --n-threads ${THREADS} \
    --file-colors

echo ""
echo "=== Done ==="
echo "New reference panel: ${NEW_SUBSAMPLE}"
echo "New Themisto index: ${THEMISTO_DIR}/ab_index"
echo "Total genomes: $(wc -l < ${NEW_SUBSAMPLE}/themisto_input.txt)"
echo ""
echo "Next: re-run mSWEEP with the new index"
