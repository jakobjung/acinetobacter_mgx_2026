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
# PE WGS from endotracheal aspirates, NovaSeq
for SRR in SRR13160906 SRR13160908 SRR13160909 SRR13160924 SRR13160943; do
    if [ -f "${OUTDIR}/${SRR}_1.fastq.gz" ]; then
        echo "Skipping ${SRR} (already downloaded)"
        continue
    fi
    echo "Downloading ${SRR}..."
    wget -q -P ${OUTDIR} \
        ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR131/${SRR: -3}/${SRR}/${SRR}_1.fastq.gz \
        ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR131/${SRR: -3}/${SRR}/${SRR}_2.fastq.gz
    echo "${SRR} done."
done

echo "Downloaded test samples to ${OUTDIR}/"
ls -lh ${OUTDIR}/
