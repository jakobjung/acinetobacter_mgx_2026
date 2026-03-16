#!/bin/bash
#SBATCH --job-name=download_mNGS
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --output=download_%j.log

mkdir -p /home/jjung/project/data/fastq

for CRA in CRA010962 CRA010735 CRA010639; do
    wget -r -np -nd -P /home/jjung/project/data/fastq/${CRA} \
         ftp://download.big.ac.cn/gsa/${CRA}/
done


