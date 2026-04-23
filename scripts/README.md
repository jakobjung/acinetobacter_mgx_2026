# Pipeline scripts

All `.sh` scripts are SLURM batch scripts designed for an HPC cluster. Run with `sbatch <script.sh>` in chronological order.

| # | Script | Description |
|---|--------|-------------|
| 01 | `01_download_ab_references.sh` | Downloads all _A. baumannii_ genome assemblies from NCBI (13,565 genomes) |
| 02 | `02_poppunk_clustering.sh` | Downloads pre-built BacPop PopPUNK _A. baumannii_ database (8,263 isolates, 420 SCs) |
| 03 | `03_subsample_and_label.sh` | Maps PopPUNK references to downloaded FASTAs, subsamples max 30/SC, runs Pasteur MLST, labels SCs, creates Themisto input files (669 genomes) |
| 04 | `04_build_themisto_index.sh` | Builds initial Themisto index (k=31) from subsampled reference panel |
| 05 | `05_download_crab_vap.sh` | Downloads all 49 WGS PE metagenome samples from PRJNA681291 (Xiao et al.) |
| 06 | `06_add_xiao_isolates.sh` | Downloads Xiao et al. isolate genomes (PRJNA679997), groups by Oxford MLST gyrB+gpi profile, adds 13 new SCs, rebuilds Themisto index v2 (715 genomes) |
| 07 | `07_run_crab_vap.sh` | Full pipeline: Kraken2 species filtering → extract _A. baumannii_ reads → Themisto pseudoalignment → mSWEEP abundance estimation |
| 08 | `08_plot_results.R` | Generates all figures: SC abundance bar chart, read count plot, SC count, heatmap, Kraken2 species composition, combined figures |
