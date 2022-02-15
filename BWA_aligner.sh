#!/bin/bash
# Script to generate job array for aligning trimmed reads to human genome GRCh38
# SureSelect XT HS2 Pipeline

today=`date +%Y-%m-%d`
DIR=$PWD
jobOutputDir=$DIR
jobName=BWA-$(basename $DIR)
AUTOSTART=0 #Default autostart off

## Process Arguments
while [ "$1" != "" ]; do
        case $1 in
		-a | --autostart )	AUTOSTART=1
					;;
		-n | --name )		shift
					jobName=$1
					;;
		-d | --directory )	shift
					if [[ -d $1 ]]; then
					        DIR=$1
					        echo "Will run on files in $DIR"
					else
					        echo "Specified directory $1 doesn't exist"
					        exit 1
					fi
					;;
		-h | --help )		echo "\
-a | --autostart	Automatically start the jobs, holding jobs so they run in the correct order
-n | --name		Sets the job name (default - UMI-VCF-$PWD)
-d | --directory	Root directory for the project
-h | --help		Display this message"
					exit 1
					;;
	esac
	shift
done

# Job script files
bwaJob=$jobOutputDir/$jobName\BWA_Job.$today\.sh
# Calculate number of samples to process (number of folders in FASTQ_Raw)
MAX=$(ls -d FASTQ_Trim/* | wc -l)

##
# Generate scripts for 'MAX' number of trimming jobs. Assumes FASTQ_Trim contains
# the output from previous trimming step: Trimmed R1 and R2 files, and the MBC file.
##

echo "
##!/bin/sh
#$ -wd $DIR		# use current working directory
#$ -o /data/scratch/$USER/
#$ -j y			# and put all output (inc errors) into it
#$ -m a			# Email on abort
#$ -pe smp 8		# Request 8 CPU cores
#$ -l h_rt=240:0:0	# Request 240 hour runtime
#$ -l h_vmem=10G		# Request 10G RAM / Core
#$ -t 1-$MAX		# run an array job of all the samples listed in FASTQ_Raw
#$ -tc 2		# only two jobs can run at the same time as the sams are massive
#$ -N $jobName-BWA

DIR=$DIR
" > $bwaJob

echo '

DIR=/data/BCI-DigitalPath/Elements/nourse/CNourse_Trial


module load bwa
module load java

## Get all the sample names from FASTQ_Trim
Samples=(ls FASTQ_Trim/*)

## Extract the file name at the position of the array job task ID
Sample=$(basename ${Samples[${SGE_TASK_ID}]})


#Move into sample directory
cd FASTQ_Trim/$Sample

#Extract file names for each FASTQ file: R1 and R2
R1=$(find -name "*R1*")
R2=$(find -name "*R2*")

#Reference genome GRCh38, pre-indexed using BWA index. 
REF=/data/BCI-DigitalPath/Genome/GRCh38_latest_genomic.fna

#Generate output file names
SAM_file=$(echo $R1 | sed -r 's/_R1.*//g' | sed -r 's/$/.sam/g')
BAM_file=$(echo $R1 | sed -r 's/_R1.*//g' | sed -r 's/$/.bam/g')


#Run BWA to align reads to GRCh38
bwa mem -C -t 8 $REF $R1 $R2 > $SAM_file

#Convert sam to bam using samtools
module load samtools

samtools view -S -b $SAM_file > $BAM_file

	
echo "Job complete"
        

' >> $bwaJob
