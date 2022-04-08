#!/bin/bash

today=`date +%Y-%m-%d`
DIR=$PWD
jobOutputDir=$DIR
jobName=Mutect2-$(basename $DIR)

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

# Get max number of files. 
MAX=$(ls -d /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/SH* | wc -l)

# Job script files
MUTECT2JOB=$jobOutputDir/$jobName\_Job.$today\.sh

##
# Generate scripts for 'MAX' number of trimming jobs. Assumes FASTQ_Trim/$Sample contains
# the output from previous step: Bam file of R1/R2 aligned to genome (GRCh38) with mate pairs
# and read groups fixed, and index files generated.
##

echo "
##!/bin/sh
#$ -cwd		# use current working directory
#$ -o /data/scratch/$USER/
#$ -j y			# and put all output (inc errors) into it
#$ -m a			# Email on abort
#$ -pe smp 1		# Request 1 CPU cores
#$ -l h_rt=240:0:0	# Request 240 hour runtime
#$ -l h_vmem=20G		# Request 20 RAM / Core
#$ -t 1-$MAX		# run an array job of all the samples listed in FASTQ_Trim
#$ -N $jobName-Mutect
" > $MUTECT2JOB

echo '

##########################################################################################################
# IMPORTANT: before running any variant calling, make a directory called "Patients" with a subdirectory  #
# for each patient (i.e. Patients/SH03). Then depending on what data is available for each patient       #
# add further subdirectories called "M", "T" and "N" for Metastatic, Tumour and Normal.                  #
#                                                                                                        #
# For example, if patient SH03 has all sample types and patient SH07 has only Normal and Tumour samples  #
# then the directory structure should be as follows:                                                     #
# -Patients                                                                                              #
#   --SH03                                                                                               #
#     ---N                                                                                               #
#     ---M                                                                                               #
#     ---T                                                                                               #
#   --SH07                                                                                               #
#     ---N                                                                                               #
#     ---T                                                                                               #
###########################################################################################################

Patients=($(ls -d /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/SH*))

## Extract the file name at the position of the array job task ID
Patient=$(basename ${Patients[${SGE_TASK_ID}-1]})

#Exit job if no patient ID for some reason
if [ "$Patient" = "" ]; then
	echo "No patient ID found, exiting..."
    exit 1
fi

echo "Running Mutect2 for sample: $Patient"

#Move into sample directory
cd /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/$Patient

if ! [[ -s N ]] && ! [[ -s T ]] && ! [[ -s M ]]; then
  	echo "Before running script, generate a directory called `Patients` with a subdirectory for each patient (i.e. Patients/SH03). 
  	Then depending on what data is available for each patient add further subdirectories called `M`, `T` and `N` for Metastatic, Tumour and Normal.
  	A detailed explanation is included at beginning of this script."
    exit 1
fi
 
if ! [[ -s N ]]; then
	echo "No somatic file found, run HaplotypeCaller instead."
    exit 1
fi


if [[ -s N ]]; then
	Normal_BAM=$(find /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/FASTQ_Trim/ -name "RG_fixed_mate_*_N_*locatIt_filtered.bam" -name "*${Patient}*")
	echo "Somatic/normal file found: $Normal_BAM"
fi


if [[ -s T ]]; then
	Tumour_BAM=$(find /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/FASTQ_Trim/  -name "RG_fixed_mate_*_T_*locatIt_filtered.bam" -name "*${Patient}*")
	echo "Tumour file found: $Tumour_BAM"
fi

	
if [[ -s M ]]; then
	Metastatic_BAM=$(find /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/FASTQ_Trim/  -name "RG_fixed_mate_*M1*locatIt_filtered.bam" -name "*${Patient}*")
	echo "Metastatic sample found: $Metastatic_BAM"
fi	
	

module load gatk/4.1.6.0
module load bcftools

REF=/data/BCI-DigitalPath/Genome/GRCh38.p12.genome.fa

BED=/data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/SCC009_1_Covered.bed

mkdir -p Mutect2/

if [[ -s N ]] && [[ -s T ]]; then
	echo "Running normal versus tumour using $Tumour_BAM $Normal_BAM and $REF"
	
 	#Get name of sample to provide normal sample name to Mutect
 	gatk GetSampleName -I $Normal_BAM -O normal_name.txt
 	NORMAL_SAMPLE_NAME=$(cat normal_name.txt)
 	gatk GetSampleName -I $Tumour_BAM -O tumour_name.txt
 	TUMOUR_SAMPLE_NAME=$(cat tumour_name.txt)
	
	#Run Mutect
 	gatk --java-options "-Xmx4g" Mutect2 --reference $REF --intervals $BED --input $Tumour_BAM --input $Normal_BAM --normal $NORMAL_SAMPLE_NAME --output Mutect2/${Patient}_tumour.vcf   
	
	#Filter mutect output
 		variant_file=Mutect2/${Patient}_tumour.vcf
  	output_variant_file=Mutect2/${Patient}_tumour.filter.vcf
  	echo "Running FilterMutectCalls using $variant_file and saving to $output_variant_file" 
  	gatk --java-options "-Xmx4g" FilterMutectCalls -O $output_variant_file -V $variant_file -R $REF
  	echo "Done running FilterMutectCalls"
   
  	#Mutect2 will output vcf files with sample name 
  	#To change this so TUMOUR and NORMAL is used in the vcf, we can reheader the vcf file
   
  	#Create text file with sample names corresponding to NORMAL or TUMOUR, needed by bcftools.
  	echo "$NORMAL_SAMPLE_NAME NORMAL 
$TUMOUR_SAMPLE_NAME TUMOR" > reheader_tumour.txt
	
	HEADER_FIXED=Reheader_${Patient}_mutect_tumour.final.vcf.gz
    
	bcftools reheader -s reheader_tumour.txt $output_variant_file > $HEADER_FIXED
		
fi


if [[ -s N ]] && [[ -s M ]]; then
	echo "Running normal versus metastatic"
	
 	#Get name of sample to provide normal sample name to Mutect
 	gatk GetSampleName -I $Normal_BAM -O normal_name.txt
 	NORMAL_SAMPLE_NAME=$(cat normal_name.txt)
 	gatk GetSampleName -I $Metastatic_BAM -O metastatic_name.txt
 	METASTATIC_SAMPLE_NAME=$(cat metastatic_name.txt)
	
	#Run Mutect
 	gatk --java-options "-Xmx4g" Mutect2 --reference $REF --intervals $BED --input $Metastatic_BAM --input $Normal_BAM --normal $NORMAL_SAMPLE_NAME --output Mutect2/${Patient}_metastatic.vcf
	
	#Filter mutect output
  	variant_file=Mutect2/${Patient}_metastatic.vcf
  	output_variant_file=Mutect2/${Patient}_metastatic.filter.vcf
  	echo "Running FilterMutectCalls using $variant_file and saving to $output_variant_file" 
  	gatk --java-options "-Xmx4g" FilterMutectCalls -O $output_variant_file -V $variant_file -R $REF
  	echo "Done running FilterMutectCalls"
    
  	#Mutect2 will output vcf files with sample name 
  	#To change this so TUMOUR and NORMAL is used in the vcf, we can reheader the vcf file
   
  	#Create text file with sample names corresponding to NORMAL or TUMOUR, needed by bcftools.
    
	echo "$NORMAL_SAMPLE_NAME NORMAL 
$METASTATIC_SAMPLE_NAME TUMOR" > reheader_metastatic.txt
	
	M_HEADER_FIXED=Reheader_${Patient}_mutect_metastatic.final.vcf.gz

	bcftools reheader -s reheader_metastatic.txt $M_output_variant_file > $M_HEADER_FIXED
	
fi





' >> $MUTECT2JOB
