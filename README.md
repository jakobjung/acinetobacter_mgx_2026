# _A. baumannii_ poly-strain bloodstream infection analysis

- Data: blood mNGS from ICU patients (PRJCA016199, 3 CRA 
  accessions: CRA010962, CRA010735, CRA010639) downloading 
  to /home/jjung/project/data/fastq/
- Goal: identify GC1 (AB5075-UW, https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000963815.1/)/GC2 (BAL062, https://www.ncbi.nlm.nih.gov/nuccore/LT594095) co-infections using Themisto + mGEMS
- Pipeline: human read removal (bowtie2/hg38) → optional 
  Kraken2 A.bau filter → Themisto pseudoalign → mSWEEP 
  abundance → mGEMS read binning → clinical outcomes analysis
- Reference DB: A. baumannii global collection, PopPUNK SCs, 
  max 30 genomes/SC, labelled by dominant ST
- Single-end 75bp Illumina reads
- Following approach from Lancet Microbe 2024 (mSWEEP-mGEMS 
  pipeline) and reanalysing Frontiers FCIMB 2023 dataset