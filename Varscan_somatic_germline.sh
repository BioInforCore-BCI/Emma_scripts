#!/bin/bash

today=`date +%Y-%m-%d`
DIR=$PWD
jobOutputDir=$DIR
jobName=Varscan-$(basename $DIR)

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
VARSCANJOB=$jobOutputDir/$jobName\_Job.$today\.sh


echo "
##!/bin/sh
#$ -cwd		# use current working directory
#$ -o /data/scratch/$USER/
#$ -j y			# and put all output (inc errors) into it
#$ -m a			# Email on abort
#$ -pe smp 1		# Request 1 CPU cores
#$ -l h_rt=240:0:0	# Request 240 hour runtime
#$ -l h_vmem=50G		# Request 50G RAM / Core
#$ -t 1-$MAX		# run an array job of all the samples listed in FASTQ_Trim
#$ -l highmem
#$ -N Varscan_array
" > $VARSCANJOB

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
##########################################################################################################

Patients=($(ls -d /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/SH*))

## Extract the file name at the position of the array job task ID
Patient=$(basename ${Patients[${SGE_TASK_ID}-1]})

echo "Running Varscan for sample: $Patient"

#Move into sample directory
cd /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/$Patient

mkdir -p Varscan/

if [[ -s N ]]; then
	Normal_BAM=$(find /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/FASTQ_Trim/ -name "RG_fixed_mate_*_N_*locatIt_filtered.bam" -name "*${Patient}*")
	echo "Somatic/normal file found: $Normal_BAM
	"
fi

if [[ -s T ]]; then
	Tumour_BAM=$(find /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/FASTQ_Trim/  -name "RG_fixed_mate_*_T_*locatIt_filtered.bam" -name "*${Patient}*")
	echo "Tumour file found: $Tumour_BAM
	"
fi
	
if [[ -s M ]]; then
	Metastatic_BAM=$(find /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/FASTQ_Trim/  -name "RG_fixed_mate_*M1*locatIt_filtered.bam" -name "*${Patient}*")
	echo "Metastatic sample found: $Metastatic_BAM
	"
fi	
	

REF=/data/BCI-DigitalPath/Genome/GRCh38.p12.genome.fa

VARSCAN=/data/home/hfy041/VarScan.v2.3.9.jar

module load samtools
module load annovar
module load java
module load bcftools



if [[ -s N ]] && [[ -s T ]]; then
	
  cd /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/$Patient
	mkdir -p Varscan/Tumour
	cd Varscan/Tumour/
	
	echo "Running mpileup of normal versus tumour using $Tumour_BAM $Normal_BAM and $REF"
	samtools mpileup -f $REF -q 1 -B $Normal_BAM $Tumour_BAM > ${Patient}_normal-tumor.mpileup 
	
	echo "Running somatic varscan of normal versus tumour"
	java -jar /data/home/hfy041/Varscan-2.3.9/VarScan.v2.3.9.jar somatic ${Patient}_normal-tumor.mpileup ${Patient}_normal-tumour \
    --mpileup 1 \
    --min-coverage 10 \
    --min-avg-qual 20 \
    --min-read2 4 \
    --p-value 1 \
    --min-var-freq 0.01 \
    --strand-filter 1 \
    --output-vcf 1
    
    
	echo "Processing VarScan outputs by somatic status"
  java -jar /data/home/hfy041/Varscan-2.3.9/VarScan.v2.3.9.jar processSomatic *.snp.vcf --max-normal-freq 0.01
	java -jar /data/home/hfy041/Varscan-2.3.9/VarScan.v2.3.9.jar processSomatic *.indel.vcf --max-normal-freq 0.01
    
	bgzip ${Patient}_normal-tumour.snp.Somatic.vcf
	bgzip ${Patient}_normal-tumour.indel.Somatic.vcf 
    
	VARSCAN_T_SNPS=${Patient}_normal-tumour.snp.Somatic.vcf.gz
	VARSCAN_T_INDELS=${Patient}_normal-tumour.indel.Somatic.vcf.gz 
	VARSCAN_T_ALL=${Patient}_normal-tumour.all.Somatic.vcf.gz

	#merge varscan vcfs with bcftools
	bcftools index $VARSCAN_T_SNPS
	bcftools index $VARSCAN_T_INDELS
	bcftools concat -a $VARSCAN_T_SNPS $VARSCAN_T_INDELS -o $VARSCAN_T_ALL

	VARSCAN_T_ALL_FILTERED=/data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/${Patient}/${Patient}_varscan_tumour.final.vcf
' > $VARSCANJOB

echo "
  bcftools view -i \"%FILTER='PASS'\" $VARSCAN_T_ALL > $VARSCAN_T_ALL_FILTERED
" > $VARSCANJOB

echo '
      
	echo "Done running normal versus tumour"
fi



if [[ -s N ]] && [[ -s M ]]; then
	
  cd /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/$Patient
	mkdir -p Varscan/Metastatic
	cd Varscan/Metastatic/
	
	echo "Running mpileup of normal versus metastatic using $Metastatic_BAM $Normal_BAM and $REF"
	samtools mpileup -f $REF -q 1 -B $Normal_BAM $Metastatic_BAM > ${Patient}_normal-metastatic.mpileup 
	
	
	echo "Running somatic varscan of normal versus metastatic"
	java -jar /data/home/hfy041/Varscan-2.3.9/VarScan.v2.3.9.jar somatic ${Patient}_normal-metastatic.mpileup ${Patient}_normal-metastatic \
    --mpileup 1 \
    --min-coverage 10 \
    --min-avg-qual 20 \
    --min-read2 4 \
    --p-value 1 \
    --min-var-freq 0.01 \
    --strand-filter 1 \
    --output-vcf 1
    
    
	echo "Processing VarScan outputs by somatic status and confidence"
  java -jar /data/home/hfy041/Varscan-2.3.9/VarScan.v2.3.9.jar processSomatic *.snp.vcf --max-normal-freq 0.01
  java -jar /data/home/hfy041/Varscan-2.3.9/VarScan.v2.3.9.jar processSomatic *.indel.vcf --max-normal-freq 0.01
    
	bgzip ${Patient}_normal-metastatic.snp.Somatic.vcf
	bgzip ${Patient}_normal-metastatic.indel.Somatic.vcf   
  
  VARSCAN_M_SNPS=${Patient}_normal-metastatic.snp.Somatic.vcf.gz
	VARSCAN_M_INDELS=${Patient}_normal-metastatic.indel.Somatic.vcf.gz 
	VARSCAN_M_ALL=${Patient}_normal-metastatic.all.Somatic.vcf.gz


	#merge varscan vcfs with bcftools
	bcftools index $VARSCAN_M_SNPS
	bcftools index $VARSCAN_M_INDELS
	bcftools concat -a $VARSCAN_M_SNPS $VARSCAN_M_INDELS -o $VARSCAN_M_ALL


	VARSCAN_M_ALL_FILTERED=/data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/${Patient}/${Patient}_varscan_metastatic.final.vcf
' >> $VARSCANJOB

echo "
	bcftools view -i "%FILTER="PASS"" $VARSCAN_M_ALL > $VARSCAN_M_ALL_FILTERED
" >> $VARSCANJOB
    
echo '
	
	echo "Done running normal versus metastatic"
	
fi



if ! [[ -s N ]] && [[ -s T ]]; then
	
  cd /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/$Patient
	mkdir -p Varscan/Tumour
	cd Varscan/Tumour/
	
	echo "Running mpileup of germline tumour using $Tumour_BAM and $REF"
	samtools mpileup -f $REF -q 1 -B $Tumour_BAM > ${Patient}_germline-tumor.mpileup 
	
	echo "Running germline varscan of normal versus tumour"
	java -jar /data/home/hfy041/Varscan-2.3.9/VarScan.v2.3.9.jar mpileup2cns ${Patient}_germline-tumor.mpileup \
    --variants 1 \
    --mpileup 1 \
    --min-coverage 20 \
    --min-avg-qual 20 \
    --min-read2 4 \
    --p-value 1 \
    --min-var-freq 0.01 \
    --strand-filter 1 \
    --output-vcf 1 > ${Patient}_germline-tumor.vcf


	VARSCAN_T_ALL_FILTERED=/data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/${Patient}/${Patient}_varscan_germline-tumor.final.vcf
' >> $VARSCANJOB

echo "
  bcftools view -i \"%FILTER='PASS'\" ${Patient}_germline-tumor.vcf > $VARSCAN_T_ALL_FILTERED
" >> $VARSCANJOB

echo '
      
	echo "Done running germline tumour"
fi



if ! [[ -s N ]] && [[ -s M ]]; then
	
  cd /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/$Patient
	mkdir -p Varscan/Metastatic
	cd Varscan/Metastatic/
	
	echo "Running mpileup of germline metastatic using $Metastatic_BAM and $REF"
	samtools mpileup -f $REF -q 1 -B $Metastatic_BAM > ${Patient}_germline-metastatic.mpileup
	
	echo "Running germline varscan of normal versus tumour"
	java -jar /data/home/hfy041/Varscan-2.3.9/VarScan.v2.3.9.jar mpileup2cns ${Patient}_germline-metastatic.mpileup \
    --variants 1 \
    --mpileup 1 \
    --min-coverage 20 \
    --min-avg-qual 20 \
    --min-read2 4 \
    --p-value 1 \
    --min-var-freq 0.01 \
    --strand-filter 1 \
    --output-vcf 1 > ${Patient}_germline-metastatic.vcf


	VARSCAN_M_ALL_FILTERED=/data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/${Patient}/${Patient}_varscan_germline-metastatic.final.vcf
' >> $VARSCANJOB

echo "
	bcftools view -i "%FILTER='PASS'" ${Patient}_germline-metastatic.vcf > $VARSCAN_M_ALL_FILTERED
" >> $VARSCANJOB
      
	echo "Done running germline tumour
	"
fi



' >> $VARSCANJOB
