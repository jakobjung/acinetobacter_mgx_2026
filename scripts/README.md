# Pipeline scripts

All scripts are SLURM batch scripts designed for an HPC cluster. Run with `sbatch <script.sh>`.

## Reference database construction

| Script | Description |
|--------|-------------|
| `download_ab_references.sh` | Downloads all _A. baumannii_ genome assemblies from NCBI using `datasets` (13,565 genomes) |
| `poppunk_clustering.sh` | Downloads pre-built BacPop PopPUNK _A. baumannii_ database (8,263 isolates, 420 SCs) from EBI FTP |
| `subsample_and_label.sh` | Maps PopPUNK reference accessions to downloaded FASTAs, subsamples max 30/SC, runs Pasteur MLST, labels SCs by dominant ST, creates Themisto input files |
| `build_themisto_index.sh` | Builds Themisto colored de Bruijn graph index (k=31) from subsampled reference panel |
| `add_xiao_isolates.sh` | Downloads 46 Xiao et al. CRAB isolate genomes (PRJNA679997), groups by Oxford MLST gyrB+gpi profile, adds as 13 new SCs, rebuilds Themisto index (v2, 715 genomes) |

## Data download

| Script | Description |
|--------|-------------|
| `download_crab_vap.sh` | Downloads 5 test samples from PRJNA681291 (Xiao et al. tracheal aspirate metagenomics) |
| `download_all_crab_vap.sh` | Downloads all 49 WGS paired-end metagenome samples from PRJNA681291 |

## Analysis pipeline

| Script | Description |
|--------|-------------|
| `run_crab_vap_test.sh` | Full pipeline (human removal + Themisto + mSWEEP) on 5 test samples with v1 reference panel |
| `run_crab_vap_all.sh` | Full pipeline on all 49 samples with v1 reference panel (669 genomes) |
| `run_crab_vap_v2.sh` | Full pipeline on all 49 samples with v2 reference panel (715 genomes, includes Xiao isolates) |
| `run_crab_vap_array.sh` | SLURM array version for parallel processing (not used in final analysis) |

## Validation

| Script | Description |
|--------|-------------|
| `validate_crab_vap.sh` | Mash distance validation of mSWEEP SC assignments — extracts binned reads, sketches with mash, compares to SC reference genomes |
