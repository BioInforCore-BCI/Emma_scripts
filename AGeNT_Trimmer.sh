#!/bin/bash
# Script to generate job array for trimming using AGeNT_2.0.5
# This removes adaptor sequences and extracts molecular barcodes for SureSelect XT HS2

today=`date +%Y-%m-%d`
DIR=$PWD
jobOutputDir=$DIR
jobName=Trim-$(basename $DIR)
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
trimJob=$jobOutputDir/$jobName\Trim_Job.$today\.sh
# Calculate number of samples to process (number of folders in FASTQ_Raw)
MAX=$(ls -d FASTQ_Raw/* | wc -l)

##
# Generate scripts for 'MAX' number of trimming jobs. Assumes FASTQ_Raw contains a 
# subdirectory for each patient sample, and contained within each of these is two sets 
# (lane 1 and lane 2) of forward and reverse reads. 
##

echo "
##!/bin/sh
#$ -wd $DIR		# use current working directory
#$ -o /data/scratch/$USER/
#$ -j y			# and put all output (inc errors) into it
#$ -m a			# Email on abort
#$ -pe smp 1		# Request 1 CPU cores
#$ -l h_rt=72:0:0	# Request 72 hour runtime
#$ -l h_vmem=5G		# Request 5G RAM / Core
#$ -t 1-$MAX		# run an array job of all the samples listed in FASTQ_Raw
#$ -N $jobName-Trim_Job

DIR=$DIR
" > $trimJob

echo '

DIR=/data/BCI-DigitalPath/Elements/nourse/CNourse_Trial

module load java

## Get all the sample names from FASTQ_Raw
Samples=(ls FASTQ_Raw/*)

## Extract the file name at the position of the array job task ID
Sample=$(basename ${Samples[${SGE_TASK_ID}]})

## Make directory for output
mkdir FASTQ_Trim/$Sample

#Move into sample directory
cd FASTQ_Raw/$Sample

#Extract file names for each FASTQ file, lanes 1 and 2, R1 and R2
R1L1=$(find -name "*L001*R1*")
R1L2=$(find -name "*L002*R1*")
R2L1=$(find -name "*L001*R2*")
R2L2=$(find -name "*L002*R2*")

#Generate output file names
R1=$(echo $R1L1 | sed -r 's/_L00.*//g' | sed -r 's/$/_R1.fastq.gz/g')
R2=$(echo $R2L1 | sed -r 's/_L00.*//g' | sed -r 's/$/_R2.fastq.gz/g')

#Concatenate lanes 1 and 2
cat $R1L1 $R1L2 > "$R1" 
cat $R2L1 $R2L2 > "$R2"


# Trim adapters using AGeNT
# AGeNT will output results in a single file for each read. 
# Also outputs a single merged MBC sequence file for R1/R2
java -jar /data/home/hfy041/AGeNT_2.0.5/agent/lib/trimmer-2.0.3.jar \
	-v2 \
    -fq1 $R1 \
	-fq2 $R2 \
	-out_loc $DIR/FASTQ_Trim/$Sample

	
echo "Job complete"
        

' >> $trimJob
