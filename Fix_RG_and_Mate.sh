#!/bin/bash
# Script to generate job array for fixing read groups and mate pairs. 
# SureSelect XT HS2 Pipeline

today=`date +%Y-%m-%d`
DIR=$PWD
jobOutputDir=$DIR
jobName=Fix_read_group-$(basename $DIR)
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
ValidateJob=$jobOutputDir/$jobName\_Job.$today\.sh
# Calculate number of samples to process (number of folders in FASTQ_Trim)
MAX=$(ls -d FASTQ_Trim/* | wc -l)

##
# Generate scripts for 'MAX' number of trimming jobs. Assumes FASTQ_Trim/$Sample contains
# the output from previous step: Bam file of R1/R2 aligned to genome (GRCh38).
##

echo "
##!/bin/sh
#$ -wd $DIR		# use current working directory
#$ -o /data/scratch/$USER/
#$ -j y			# and put all output (inc errors) into it
#$ -m a			# Email on abort
#$ -pe smp 1		# Request 1 CPU cores
#$ -l h_rt=1:0:0	# Request 72 hour runtime
#$ -l h_vmem=20G		# Request 12G RAM / Core
#$ -t 1-$MAX		# run an array job of all the samples listed in FASTQ_Trim
#$ -N $jobName-Fix_read_groups

DIR=$DIR
" > $ValidateJob

echo '

DIR=/data/BCI-DigitalPath/Elements/nourse/CNourse_Trial

module load gatk


## Get all the sample names from FASTQ_Trim
Samples=(ls FASTQ_Trim/*)

## Extract the file name at the position of the array job task ID
Sample=$(basename ${Samples[${SGE_TASK_ID}]})

#Move into sample directory
cd FASTQ_Trim/$Sample

echo "Working in FASTQ_Trim/${Sample} "


#Extract file name for bam file
BAM=$(basename $(find -maxdepth 1 -name "*locatIt_filtered.bam"))

echo "Fixing mate pairs for file: $BAM"

gatk FixMateInformation I=$BAM O=fixed_mate_${BAM}


#Extract file name for fixed mate pair bam file
BAM=$(basename $(find -maxdepth 1 -name "fixed_mate_*locatIt_filtered.bam"))

echo "Fixing read groups for file: $BAM"

#Store variables for command
OUTPUT=RG_${BAM}
LIB=lib_${Sample}
ID=sample_${Sample}

echo "Using ID: $ID and library: $LIB "

#The ID/lib used do not really matter, this is just to make bam compatible with later pipeline steps
gatk AddOrReplaceReadGroups \
I=$BAM \
O=$OUTPUT \
RGID=$ID \
RGLB=$LIB \
RGPL=ILLUMINA \
RGPU=$ID \
RGSM=$ID

echo "Job complete"
        

' >> $ValidateJob
