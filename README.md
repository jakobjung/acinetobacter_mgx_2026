# Metagenomics-based detection of _Acinetobacter baumannii_ poly-strain co-colonisation

Detection of within-GC2 sub-lineage diversity in carbapenem-resistant _A. baumannii_ (CRAB) from tracheal aspirate metagenomics, using the Themisto/mSWEEP pipeline.

## Dataset

Re-analysis of tracheal aspirate shotgun metagenomics from Xiao et al. 2022 (PMID: 35308401), a CRAB VAP study from a Chinese ICU (PRJNA681291). 49 samples: 24 CRAB-VAP (CRAB-I), 22 CRAB-colonisation (CRAB-C), 6 CRAB-negative (CRAB-N).

- **Sequencing**: Illumina NovaSeq 6000, PE 150bp
- **Sample type**: Endotracheal deep aspirates from mechanically ventilated ICU patients
- **Metagenomics accession**: PRJNA681291 (49 WGS paired-end samples)
- **WGS isolates accession**: PRJNA679997 (46 CRAB isolates from the same patients)

## Pipeline

1. **Human read removal** — bowtie2 (sensitive mode) against GRCh38, PE mode
2. **Reference database construction** — BacPop PopPUNK _A. baumannii_ database (8,263 isolates, 420 SCs), subsampled to max 30 genomes/SC (669 genomes), labelled by Pasteur MLST dominant ST
3. **Reference panel augmentation (v2)** — Added 46 Xiao et al. isolate genomes grouped into 13 new SCs by Oxford MLST gyrB+gpi allelic profile (715 genomes total)
4. **Themisto index** — Colored de Bruijn graph (k=31) built from the 715-genome reference panel
5. **Themisto pseudoalignment** — PE reads pseudoaligned against index (separate R1/R2, sorted output)
6. **mSWEEP** — Probabilistic abundance estimation per SC, with read binning (--min-abundance 0.01)
7. **Mash validation** — Mash distance comparison of binned reads against SC reference genomes

## Key findings

- SC_327_ST255 (Pasteur) is the dominant background strain across all samples, including CRAB-negative controls — likely represents the local ICU outbreak clone (Oxford ST208)
- CRAB-positive samples show additional GC2 sub-lineages (gyrB3_gpi97/ST208, gyrB38_gpi110/ST938, gyrB3_gpi157, gyrB3_gpi61, etc.) at 1-88% abundance
- ~35/43 CRAB-positive samples have 2+ sub-lineages at >5%, indicating widespread within-GC2 poly-strain co-colonisation
- CRAB-negative controls show only the background strain with no sub-lineage diversity
- One sample (S3/SRR13160953) contains SC_26_ST46 at 32% — a non-GC2 lineage, representing a genuine cross-lineage co-infection
- Culture-based typing (one colony per patient) reported only one Oxford ST per patient; metagenomics reveals additional sub-lineages missed by culture

## Methodology reference

Pipeline approach based on Thorpe et al. 2024, "Pan-pathogen deep sequencing of nosocomial bacterial pathogens in Italy", _The Lancet Microbe_. DOI: 10.1016/S2666-5247(24)00113-7

## Repository structure

- `scripts/` — SLURM pipeline scripts
- `analysis/` — Summary output files
- `data/` — Input data references and download metadata
