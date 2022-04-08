#!/bin/bash
# Script to generate job array for indexing bam files
# SureSelect XT HS2 Pipeline

today=`date +%Y-%m-%d`
DIR=$PWD
jobOutputDir=$DIR
jobName=Index-$(basename $DIR)
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
IndexJob=$jobOutputDir/$jobName\_Job.$today\.sh

# Calculate number of samples to process (number of folders in FASTQ_Trim)
MAX=$(ls -d FASTQ_Trim/* | wc -l)

##
# Generate scripts for 'MAX' number of trimming jobs. Assumes FASTQ_Trim/$Sample contains
# the output from previous step: Bam file of R1/R2 aligned to genome (GRCh38) with read groups
# and mate pairs fixed. 
##

echo "
##!/bin/sh
#$ -wd $DIR		# use current working directory
#$ -o /data/scratch/$USER/
#$ -j y			# and put all output (inc errors) into it
#$ -m a			# Email on abort
#$ -pe smp 1		# Request 1 CPU cores
#$ -l h_rt=1:0:0	# Request 1 hour runtime
#$ -l h_vmem=20G		# Request 20G RAM / Core
#$ -t 1-$MAX		# run an array job of all the samples listed in FASTQ_Trim
#$ -N $jobName-Index

DIR=$DIR
" > $IndexJob

echo '

DIR=/data/BCI-DigitalPath/Elements/nourse/CNourse_Trial

## Get all the sample names from FASTQ_Trim
Samples=(ls FASTQ_Trim/*)

## Extract the file name at the position of the array job task ID
Sample=$(basename ${Samples[${SGE_TASK_ID}]})

#Move into sample directory
cd FASTQ_Trim/$Sample

#Extract file name for bam file
BAM=$(basename $(find -maxdepth 1 -name "RG_fixed_mate_*locatIt_filtered.bam"))

echo " Working in FASTQ_Trim/${Sample}, indexing $BAM"

module load gatk

gatk BuildBamIndex I=$BAM

echo "Job complete"
        

' >> $IndexJob
