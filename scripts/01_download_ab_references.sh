#!/bin/bash
#SBATCH --job-name=download_ab_refs
#SBATCH --partition=cpu
#SBATCH --qos=long
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --output=download_ab_refs_%j.log

OUTDIR=/vol/projects/jjung/acinetobacter_mgx_2026/data/reference_sequences/ab_genomes

mkdir -p ${OUTDIR}

# download all A. baumannii assemblies (dehydrated = metadata + stubs first)
datasets download genome taxon "Acinetobacter baumannii" \
    --assembly-level complete,chromosome,scaffold \
    --include genome \
    --dehydrated \
    --filename ${OUTDIR}/ab_genomes.zip

# unzip and rehydrate (actual sequence download)
unzip -o ${OUTDIR}/ab_genomes.zip -d ${OUTDIR}/ncbi_dataset
datasets rehydrate --directory ${OUTDIR}/ncbi_dataset/

echo "Download complete. Genome count:"
find ${OUTDIR}/ncbi_dataset -name "*.fna" | wc -l
