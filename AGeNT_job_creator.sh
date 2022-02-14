#!/bin/bash
# Script to generate job array for trimming using AGeNT_2.0.5

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
                -f | --fastq-suffix )   shift
                                        fastqSuffix=$1
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
# Generate scripts for 'MAX' number of trimming Job
##

echo "
##!/bin/sh
#$ -wd $DIR		# use current working directory
#$ -o /data/scratch/$USER/
#$ -j y			# and put all output (inc errors) into it
#$ -m a			# Email on abort
#$ -pe smp 1		# Request 1 CPU cores
#$ -l h_rt=24:0:0	# Request 24 hour runtime
#$ -l h_vmem=4G		# Request 4G RAM / Core
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

R1L1=$(find FASTQ_Raw/$Sample/ -name "*L001*R1*")
R1L2=$(find FASTQ_Raw/$Sample/ -name "*L002*R1*")
R2L1=$(find FASTQ_Raw/$Sample/ -name "*L001*R2*")
R2L2=$(find FASTQ_Raw/$Sample/ -name "*L002*R2*")


# Trim adapters using AGeNT
# /data/home/hfy041/AGeNT_2.0.5/agent/agent.sh 

cd FASTQ_Trim/$Sample


java -jar /data/home/hfy041/AGeNT_2.0.5/agent/lib/trimmer-2.0.3.jar \
	-v2 \
    -fq1 $DIR/$R1L1,$DIR/$R1L2 \
	-fq2 $DIR/$R2L1,$DIR/$R2L2 
        

' >> $trimJob
