#!/bin/bash
#SBATCH --job-name=dl_crab_all
#SBATCH --partition=cpu
#SBATCH --qos=long
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --time=2-00:00:00
#SBATCH --output=dl_crab_all_%j.log

set -euo pipefail

PROJDIR=/vol/projects/jjung/acinetobacter_mgx_2026
OUTDIR=${PROJDIR}/data/crab_vap_test
URLFILE=${PROJDIR}/data/crab_vap_urls.tsv

mkdir -p ${OUTDIR}

while IFS=$'\t' read -r srr urls; do
    if [ -f "${OUTDIR}/${srr}_1.fastq.gz" ] && [ -f "${OUTDIR}/${srr}_2.fastq.gz" ]; then
        echo "Skipping ${srr} (already downloaded)"
        continue
    fi

    r1=$(echo ${urls} | cut -d';' -f1)
    r2=$(echo ${urls} | cut -d';' -f2)

    echo "Downloading ${srr}..."
    wget -c -q -P ${OUTDIR} ftp://${r1} ftp://${r2}
    echo "${srr} done."
done < ${URLFILE}

echo "All downloads complete."
ls -lh ${OUTDIR}/ | tail -5
echo "Total files: $(ls ${OUTDIR}/*.fastq.gz | wc -l)"
