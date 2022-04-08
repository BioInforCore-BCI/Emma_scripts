#!/bin/bash

today=`date +%Y-%m-%d`
DIR=$PWD
jobOutputDir=$DIR
jobName=HaplotypeCaller-$(basename $DIR)

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
HaploCallerJOB=$jobOutputDir/$jobName\_Job.$today\.sh


echo "
##!/bin/sh
#$ -cwd		# use current working directory
#$ -o /data/scratch/$USER/
#$ -j y			# and put all output (inc errors) into it
#$ -m a			# Email on abort
#$ -pe smp 1		# Request 1 CPU cores
#$ -l h_rt=240:0:0	# Request 240 hour runtime
#$ -l h_vmem=20G		# Request 20G RAM / Core
#$ -t 1-$MAX		# run an array job of all the samples listed in FASTQ_Trim
#$ -N $jobName
" > $HaploCallerJOB

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

echo "Running HaplotypeCaller for sample: $Patient"

#Move into sample directory
cd /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/$Patient


if [[ -s N ]]; then
	echo "Somatic/normal file found, run Mutect2 instead"
  exit 1
fi

if ! [[ -s N ]]; then
	echo "No somatic/normal file found, running HaplotypeCaller for germline mode"
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
module load bcftools

REF=/data/BCI-DigitalPath/Genome/GRCh38.p12.genome.fa

BED=/data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/SCC009_1_Covered.bed

mkdir -p HaplotypeCaller/



if ! [[ -s N ]] && [[ -s T ]]; then

	#Run Mutect
	cd /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/${Patient}/HaplotypeCaller/
	
	mkdir -p Tumour/
	echo "Running germline versus tumour"	

 	#Run HaplotypeCaller and filter variants 
	gatk --java-options "-Xmx4g" GenotypeGVCFs -R $REF -V Tumour/$Patient.vcf -O Tumour/$Patient.Geno.vcf
       gatk VariantFiltration -V Metastatic/$Patient.Geno.vcf -O Metastatic/$Patient.filter.vcf
	
	OUTFILE=/data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/${Patient}/$Patient_Haplotype_caller_tumour_germline.vcf
	
' >> $HaploCallerJOB
 
echo "
	bcftools view -i \"%FILTER='PASS'\" Metastatic/$Patient.filter.vcf > $OUTFILE
" >> $HaploCallerJOB
	
echo '
 	
 	gatk GetSampleName -I $Tumour_BAM -O tumour_name.txt
 	TUMOUR_SAMPLE_NAME=$(cat tumour_name.txt)
 	
 	echo "$TUMOUR_SAMPLE_NAME TUMOUR" > reheader_tumour.txt
	
	T_HEADER_FIXED=Reheader_${Patient}_Haplotype_caller_tumour_germline.vcf.gz 

	bcftools reheader -s reheader_tumour.txt $T_output_variant_file > $T_HEADER_FIXED
  
  echo "Done running germline versus tumour"
	
fi


if ! [[ -s N ]] && [[ -s M ]]; then
	#Run Mutect
	cd /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/${Patient}/HaplotypeCaller/
	
	mkdir -p Metastatic/
 	echo "Running germline versus metastatic"
 	gatk --java-options "-Xmx4g" HaplotypeCaller -R $REF --intervals $BED -I $Metastatic_BAM -O Metastatic/$Patient.vcf -ERC GVCF
 	
 	gatk --java-options "-Xmx4g" GenotypeGVCFs -R $REF -V Metastatic/$Patient.vcf -O Metastatic/$Patient.Geno.vcf
      gatk VariantFiltration -V Metastatic/$Patient.Geno.vcf -O Metastatic/$Patient.filter.vcf
	
	OUTFILE=/data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/${Patient}/$Patient_Haplotype_caller_metastatic_germline.vcf
	
 ' >> $HaploCallerJOB
 
echo "
	bcftools view -i \"%FILTER='PASS'\" Metastatic/$Patient.filter.vcf > $OUTFILE
" >> $HaploCallerJOB
	
echo '

  	#Mutect2 will output vcf files with sample name 
  	#To change this so TUMOUR and NORMAL is used in the vcf, we can reheader the vcf file
   
  	#Create text file with sample names corresponding to NORMAL or TUMOUR, needed by bcftools.
 	gatk GetSampleName -I $Metastatic_BAM -O metastatic_name.txt
 	METASTATIC_SAMPLE_NAME=$(cat metastatic_name.txt)
 	
 	echo "$METASTATIC_SAMPLE_NAME TUMOR" > reheader_metastatic.txt

	M_HEADER_FIXED=Reheader_${Patient}_Haplotype_caller_metastatic_germline.vcf 

	bcftools reheader -s reheader_metastatic.txt $M_output_variant_file > $M_HEADER_FIXED

	     echo "Done running germline versus metastatic"
		
fi



' >> $HaploCallerJOB
