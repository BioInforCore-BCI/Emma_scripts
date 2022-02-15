#!/bin/bash
# Script to generate job array for LocatIt
# SureSelect XT HS2 Pipeline

today=`date +%Y-%m-%d`
DIR=$PWD
jobOutputDir=$DIR
jobName=LocatIt-$(basename $DIR)
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
LocatItJob=$jobOutputDir/$jobName\LocatIt_Job.$today\.sh
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
#$ -l h_rt=240:0:0	# Request 72 hour runtime
#$ -l h_vmem=12G		# Request 12G RAM / Core
#$ -t 1-$MAX		# run an array job of all the samples listed in FASTQ_Trim
#$ -N $jobName-LocatIt

DIR=$DIR
" > $LocatItJob

echo '

DIR=/data/BCI-DigitalPath/Elements/nourse/CNourse_Trial

module load java

#Extract file name for bam file
BAM=$(find -name "*.bam")

#Bed file, hard-coded currently
BED=/data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/SCC009_1_Covered.bed                      

#Get output name
OUTPUT=$(basename $BAM | sed -r 's/.bam/_locatIt.bam/g')

java -Xmx12G -jar /data/home/hfy041/AGeNT_2.0.5/agent/lib/locatit-2.0.5.jar -S -IB -v2Duplex -l $BED -o $OUTPUT $BAM

	
echo "Job complete"
        

' >> $LocatItJob