#!/bin/bash
#SBATCH --job-name=dl_crab_vap
#SBATCH --partition=cpu
#SBATCH --qos=long
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --time=1-00:00:00
#SBATCH --output=dl_crab_vap_%j.log

set -euo pipefail

PROJDIR=/vol/projects/jjung/acinetobacter_mgx_2026
OUTDIR=${PROJDIR}/data/crab_vap_test

mkdir -p ${OUTDIR}

# Download 5 tracheal metagenome samples from PRJNA681291 (CRAB VAP study)
# Using exact ENA FTP URLs
declare -A URLS
URLS[SRR13160906]="ftp.sra.ebi.ac.uk/vol1/fastq/SRR131/006/SRR13160906"
URLS[SRR13160908]="ftp.sra.ebi.ac.uk/vol1/fastq/SRR131/008/SRR13160908"
URLS[SRR13160909]="ftp.sra.ebi.ac.uk/vol1/fastq/SRR131/009/SRR13160909"
URLS[SRR13160924]="ftp.sra.ebi.ac.uk/vol1/fastq/SRR131/024/SRR13160924"
URLS[SRR13160943]="ftp.sra.ebi.ac.uk/vol1/fastq/SRR131/043/SRR13160943"

for SRR in "${!URLS[@]}"; do
    if [ -f "${OUTDIR}/${SRR}_1.fastq.gz" ]; then
        echo "Skipping ${SRR} (already downloaded)"
        continue
    fi
    echo "Downloading ${SRR}..."
    wget -c -P ${OUTDIR} \
        ftp://${URLS[$SRR]}/${SRR}_1.fastq.gz \
        ftp://${URLS[$SRR]}/${SRR}_2.fastq.gz
    echo "${SRR} done."
done

echo "Downloaded test samples to ${OUTDIR}/"
ls -lh ${OUTDIR}/
