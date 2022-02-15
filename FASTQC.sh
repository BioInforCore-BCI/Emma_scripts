#!/bin/bash
# Script to generate job array for running FastQC on reads trimmed using AGeNT Trimmer
# SureSelect XT HS2 Pipeline

today=`date +%Y-%m-%d`
DIR=$PWD
jobOutputDir=$DIR
jobName=FASTQC-$(basename $DIR)

## Process Arguments
while [ "$1" != "" ]; do
        case $1 in
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
-n | --name		Sets the job name (default - UMI-VCF-$PWD)
-d | --directory	Root directory for the project
-h | --help		Display this message"
					exit 1
					;;
	esac
	shift
done

# Job script files
FASTQC_Job=$jobOutputDir/$jobName\FASTQC_Job.$today\.sh
# Calculate number of samples to process (number of folders in FASTQ_Raw)
MAX=$(ls -d FASTQ_Trim/* | wc -l)


echo "
#!/bin/sh
#$ -wd $DIR		# use current working directory
#$ -o /data/scratch/$USER/
#$ -j y			# and put all output (inc errors) into it
#$ -m e			# Email on abort
#$ -pe smp 1		# Request 1 CPU cores
#$ -l h_rt=48:0:0	# Request 48 hour runtime (This shouldn't last more than a few minutes but in the case of large fastq might take longer)
#$ -l h_vmem=4G		# Request 4G RAM / Core
#$ -t 1-$MAX		# run an array job of all the samples listed in FASTQ_Raw
#$ -N FASTQC
FILEDIR=$FILEDIR" > $FASTQC_Job

echo '

## Get all the sample names from FASTQ_Trim
Samples=(ls FASTQ_Trim/*)

## Extract the file name at the position of the array job task ID
Sample=$(basename ${Samples[${SGE_TASK_ID}]})


#Move into sample directory
cd FASTQ_Trim/$Sample

#Create QC directory if does not exist
if ! [[ -f QC/ ]]; then mkdir QC; fi

module load fastqc

#Extract file names for each FASTQ file: R1 and R2
R1=$(find -name "*R1*")
R2=$(find -name "*R2*")

fastqc -o QC/ $R1
fastqc -o QC/ $R2
' >> $FASTQC_Job
