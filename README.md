
# Emma's Scripts

Scripts to generate array jobs.

### AGeNT Pipeline:

1. [AGeNT_Trimmer.sh](#agent_trimmersh)

1. [BWA_aligner.sh](#bwa_alignersh)

1. [AGeNT_LocatIt.sh](#agent_locatitsh)

### Quality Checks:

1. [FASTQC.sh](#fastqcsh)

***

## AGeNT Pipeline:

#### AGeNT_Trimmer.sh

- Trims adaptor sequences and extracts molecular barcodes for SureSelect XT HS2. 
- Assumes a project structure as below. Concatenates two lanes into one for R1/R2. 
- Outputs a trimmed R1 and R2 file, and a single MBC file, into SCC_Trial/FASTQ_Trim.

Project root | raw | sample | .fastq.gz
--- | --- | --- | ---
SCC_Trial | FASTQ_Raw | 1 | 1_L001_R1.fastq.gz
|  |  |  | 1_L001_R2.fastq.gz
| |  |  | 1_L002_R1.fastq.gz
|  |  |  | 1_L002_R2.fastq.gz
| |  | 2 | 2_L001_R1.fastq.gz
|  |  |  | 2_L001_R2.fastq.gz
| |  |  | 2_L002_R1.fastq.gz
|  |  |  | 2_L002_R2.fastq.gz



#### BWA_aligner.sh

- Aligns trimmed FASTQ files to reference genome GRCh38 using BWA-mem. 
- Uses samtools to convert SAM to BAM. 
- Outputs a BAM file which retains MBC flags from FASTQ headers. 
- Reference genome path hardcoded and pre-indexed using BWA index. 


#### AGeNT_LocatIt.sh

- Tags read pairs in a BAM file with their MBC sequences and mark or merge MBC duplicates. 
- Requires BAM file annotated with MBC tags, produced using AGeNT_Trimmer.sh and BWA_aligner.sh. 
- Outputs a sorted BAM file. 
- Run in v2Duplex mode, i.e. 'duplex consensus mode'. 
- Uses a 'covered.bed' file from Agilent's SureDesign, specific to the panel. Path is currently hardcoded.   


***

## Quality Checks:

#### FASTQC.sh

-Runs FastQC on trimmed reads (R1/R2 files) produced using AGeNT_trimmer.sh. 
