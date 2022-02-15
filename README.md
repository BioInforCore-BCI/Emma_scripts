# Emma_scripts

## AGeNT Pipeline

Scripts to generate array jobs. 

1. [AGeNT_Trimmer.sh](#agent_trimmersh)

1. [BWA_aligner.sh](#bwa_alignersh)


***


### AGeNT_Trimmer.sh

- Trims adaptor sequences and extracts molecular barcodes for SureSelect XT HS2. \
- Assumes a project structure as below. Concatenates multiple lanes into one for R1/R2. \
- Outputs a trimmed R1 and R2 file, and a single MBC file, into SCC_Trial/FASTQ_Trim.

Project root | raw | sample | .fastq.gz
--- | --- | --- | ---
SCC_Trial | FASTQ_Raw | 1 | 1_L001_R1.fastq.gz
|  |  | 1 | 1_L001_R2.fastq.gz
| |  |  | 1_L002_R1.fastq.gz
|  |  |  | 1_L002_R2.fastq.gz
| |  | 2 | 2_L001_R1.fastq.gz
|  |  |  | 2_L001_R2.fastq.gz
| |  |  | 2_L002_R1.fastq.gz
|  |  |  | 2_L002_R2.fastq.gz



### BWA_aligner.sh

- Aligns trimmed FASTQ files to reference genome GRCh38 using BWA-mem. \
- Uses samtools to convert SAM to BAM. \
- Outputs a bam file which retains MBC flags from FASTQ headers. \
- Reference genome path hardcoded and pre-indexed using BWA index. 
