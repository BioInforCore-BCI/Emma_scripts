#!/bin/bash

today=`date +%Y-%m-%d`
DIR=$PWD
jobOutputDir=$DIR
jobName=Strelka-$(basename $DIR)

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
STRELKAJOB=$jobOutputDir/$jobName\_Job.$today\.sh


echo "
##!/bin/sh
#$ -cwd		# use current working directory
#$ -o /data/scratch/$USER/
#$ -j y			# and put all output (inc errors) into it
#$ -m a			# Email on abort
#$ -pe smp 8		# Request 8 CPU cores
#$ -l h_rt=240:0:0	# Request 240 hour runtime
#$ -l h_vmem=5G		# Request 5G RAM / Core
#$ -t 1-$MAX		# run an array job of all the samples listed in FASTQ_Trim
#$ -l highmem
#$ -N Strelka_array
" > $STRELKAJOB

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

if [ "$Patient" = "" ]; then
	echo "No patient ID found, exiting..."
    exit 1
fi

echo "Running Strelka for sample: $Patient"

#Move into sample directory
cd /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/$Patient

mkdir -p Strelka/

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
module load annovar
module load bcftools/1.13 
module load python/3.6.3	
module load samtools

# Dependency: Add_GT_vcf.py script


REF=/data/BCI-DigitalPath/Genome/GRCh38.p12.genome.fa

STRELKA=/data/home/hfy041/Strelka-2.9.2/bin/configureStrelkaSomaticWorkflow.py

BED=/data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/SCC009_1_Covered.bed.bgz 


if [[ -s N ]] && [[ -s T ]]; then
	echo "Running normal versus tumour using $Tumour_BAM $Normal_BAM and $REF"
 	mkdir -p Strelka/Tumour
 	$STRELKA \
     	--normalBam $Normal_BAM \
         --tumorBam $Tumour_BAM \
         --ref $REF  \
         --callRegions $BED \
         --exome \
         --runDir Strelka/Tumour
     Strelka/Tumour/runWorkflow.py -m local -j 8
	
	
 	#Strelka vcfs
 	STRELKA_T_indels=/data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/${Patient}/Strelka/Tumour/results/variants/somatic.indels.vcf.gz
 	STRELKA_T_snvs=/data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/${Patient}/Strelka/Tumour/results/variants/somatic.snvs.vcf.gz
 	STRELKA_T_VCF=/data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/${Patient}/Strelka/Tumour/results/variants/${Patient}_somatic.all.vcf.gz
 	VCF_T_UNZIP=/data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/${Patient}/Strelka/Tumour/results/variants/${Patient}_temp.vcf
	VCF_T_FIXED=/data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/${Patient}/Strelka/Tumour/results/variants/${Patient}_GT_fixed.vcf
	STRELKA_T_FILT_VCF=/data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/${Patient}/${Patient}_strelka_tumour.final.vcf

	
 	#merge strelka vcfs with bcftools
 	bcftools concat -a $STRELKA_T_indels $STRELKA_T_snvs -o $STRELKA_T_VCF
 	bgzip -c -d $STRELKA_T_VCF > $VCF_UNZIP
 
  #Strelka has no GT column in output which causes issues with bcftools in later steps.
  #Adding GT column in here
 	python /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Add_GT_vcf.py $VCF_T_UNZIP $VCF_T_FIXED

	echo "Filtering $VCF_T_FIXED and saving to $STRELKA_T_FILT_VCF"
 
' >> $STRELKAJOB

echo "
	bcftools view -i \"%FILTER='PASS'\" $VCF_T_FIXED > $STRELKA_T_FILT_VCF
" >> $STRELKAJOB

echo '

	echo "Done running normal versus tumour"

fi


if [[ -s N ]] && [[ -s M ]]; then
echo "Running normal versus metastatic"
 	mkdir -p Strelka/Metastatic
	$STRELKA \
    	--normalBam $Normal_BAM \
        --tumorBam $Metastatic_BAM \
        --ref $REF  \
        --callRegions $BED \
        --exome \
        --runDir Strelka/Metastatic
    Strelka/Metastatic/runWorkflow.py -m local -j 8
    
    
  #Strelka vcfs
 	STRELKA_M_indels=/data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/${Patient}/Strelka/Metastatic/results/variants/somatic.indels.vcf.gz
 	STRELKA_M_snvs=/data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/${Patient}/Strelka/Metastatic/results/variants/somatic.snvs.vcf.gz
 	STRELKA_M_VCF=/data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/${Patient}/Strelka/Metastatic/results/variants/${Patient}_somatic.all.vcf.gz
 	VCF_UNZIP=/data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/${Patient}/Strelka/Metastatic/results/variants/${Patient}_temp.vcf
	VCF_FIXED=/data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/${Patient}/Strelka/Metastatic/results/variants/${Patient}_GT_fixed.vcf
	STRELKA_M_FILT_VCF=/data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/${Patient}/${Patient}_strelka_metastatic.final.vcf

	#merge strelka vcfs with bcftools
	bcftools concat -a $STRELKA_M_indels $STRELKA_M_snvs -o $STRELKA_M_VCF

	bgzip -c -d $STRELKA_M_VCF > $VCF_UNZIP

  #Strelka has no GT column in output which causes issues with bcftools in later steps.
  #Adding GT column in here
	python /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Add_GT_vcf.py $VCF_UNZIP $VCF_FIXED

	echo "Filtering $VCF_FIXED and saving to $STRELKA_M_FILT_VCF"
' >> $STRELKAJOB

echo "
	bcftools view -i \"%FILTER='PASS'\" $VCF_FIXED > $STRELKA_M_FILT_VCF
" >> $STRELKAJOB
  
echo '

 	echo "Done running normal versus metastatic"
	
fi

if ! [[ -s N ]] && [[ -s T ]]; then
 	echo "Running germline tumour using $Tumour_BAM and $REF"
	mkdir -p Strelka/Tumour
	$STRELKA \
        --bam $Tumour_BAM \
        --ref $REF  \
        --callRegions $BED \
        --exome \
        --runDir Strelka/Tumour
    Strelka/Tumour/runWorkflow.py -m local -j 8

  T_VCF=/data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/${Patient}/Strelka/Tumour/results/variants/variants.vcf.gz
  VCF_T_FIXED=/data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/${Patient}/Strelka/Tumour/results/variants/${Patient}_GT_fixed.vcf
  T_FILT_VCF=/data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/${Patient}/${Patient}_strelka_germline_tumour.final.vcf

  #Strelka has no GT column in output which causes issues with bcftools in later steps.
  #Adding GT column in here
  python /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Add_GT_vcf.py $T_VCF $VCF_T_FIXED
  
  echo "Filtering $VCF_T_FIXED and saving to $T_FILT_VCF"

' >> $STRELKAJOB

echo "
	bcftools view -i \"%FILTER='PASS'\" $T_VCF > $T_FILT_VCF
" >> $STRELKAJOB

echo '

	echo "Done running normal versus tumour"
fi


if ! [[ -s N ]] && [[ -s M ]]; then
	echo "Running germline metastatic using $Metastatic_BAM and $REF"
	mkdir -p Strelka/Metastatic
	$STRELKA \
        --bam $Metastatic_BAM \
        --ref $REF  \
        --callRegions $BED \
        --exome \
        --runDir Strelka/Metastatic
    Strelka/Metastatic/runWorkflow.py -m local -j 8
    
    M_VCF=/data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/${Patient}/Strelka/Metastatic/results/variants/variants.vcf.gz
    VCF_M_FIXED=/data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/${Patient}/Strelka/Metastatic/results/variants/${Patient}_GT_fixed.vcf
    M_FILT_VCF=/data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/${Patient}/${Patient}_strelka_germline_metastatic.final.vcf
    
    
    #Strelka has no GT column in output which causes issues with bcftools in later steps.
    #Adding GT column in here
   	python /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Add_GT_vcf.py $M_VCF $VCF_M_FIXED
    
    echo "Filtering $VCF_M_FIXED and saving to $M_FILT_VCF"

' >> $STRELKAJOB

echo " 
  bcftools view -i \"%FILTER='PASS'\" $M_VCF > $M_FILT_VCF 
" >> $STRELKAJOB

echo '
	
	echo "Done running normal versus metastatic"
	
fi



' >> $STRELKAJOB
