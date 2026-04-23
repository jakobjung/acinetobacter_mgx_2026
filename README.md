# Metagenomics-based detection of _Acinetobacter baumannii_ poly-strain co-colonisation

Detection of within-GC2 sub-lineage diversity in carbapenem-resistant _A. baumannii_ (CRAB) from tracheal aspirate metagenomics, using the Themisto/mSWEEP pipeline with Kraken2 species-level pre-filtering.

## Dataset

Re-analysis of tracheal aspirate shotgun metagenomics from Xiao et al. 2022 (PMID: 35308401), a CRAB VAP study from a Chinese ICU (PRJNA681291). 49 samples: 24 CRAB-VAP (CRAB-I), 22 CRAB-colonisation (CRAB-C), 6 CRAB-negative (CRAB-N).

- **Sequencing**: Illumina NovaSeq 6000, PE 150bp
- **Sample type**: Endotracheal deep aspirates from mechanically ventilated ICU patients
- **Metagenomics accession**: PRJNA681291
- **WGS isolates accession**: PRJNA679997 (46 CRAB isolates)

## Pipeline

1. **Human read removal** — Bowtie2 against GRCh38 (PE sensitive mode)
2. **Species-level filtering** — Kraken2 classification, extract _A. baumannii_ reads only (taxid 470)
3. **Reference panel** — BacPop PopPUNK global panel (669 genomes) + 46 Xiao et al. isolates grouped by Oxford MLST (715 genomes total, 424 SCs)
4. **Themisto pseudoalignment** — Coloured de Bruijn graph index (k=31), PE pseudoalignment
5. **mSWEEP abundance estimation** — SC-level relative abundances with read binning (--min-abundance 0.01)

## Key findings

- CRAB-negative controls have essentially zero _A. baumannii_ reads (2-10 per sample)
- ~11 of 43 CRAB-positive samples show genuine multi-strain co-colonisation (2+ Oxford STs at >5% abundance)
- The dominant Oxford ST per sample matches culture-based WGS typing in nearly all cases
- In several samples, culture detected a minor sub-lineage while metagenomics reveals the true dominant population plus additional sub-lineages
- Sample S3 contains a non-GC2 lineage (Oxford ST1333, Pasteur ST46) at 71% alongside GC2 ST208 at 25% — a genuine cross-lineage co-infection missed by culture
- CRAB-I (VAP) samples show more co-colonisation complexity than CRAB-C samples

## Methodology reference

Pipeline approach based on Thorpe et al. 2024, _The Lancet Microbe_. DOI: 10.1016/S2666-5247(24)00113-7

## Repository structure

- `scripts/` — SLURM pipeline scripts and R plotting code
- `analysis/` — Summary output files, metadata, and methods description
- `data/` — Input data references and download metadata
