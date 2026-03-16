#!/bin/bash
#SBATCH --job-name=download_mNGS
#SBATCH --partition=cpu
#SBATCH --time=72:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --output=download_%j.log

mkdir -p /vol/projects/jjung/acinetobacter_mgx_2026/data/fastq

for CRA in CRA010962 CRA010735 CRA010639; do
    wget -r -np -nd --continue -P /vol/projects/jjung/acinetobacter_mgx_2026/data/fastq/${CRA} \
         ftp://download.big.ac.cn/gsa2/${CRA}/
done

# done

