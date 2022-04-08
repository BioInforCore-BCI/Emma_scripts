#!/bin/bash

today=`date +%Y-%m-%d`
DIR=$PWD
jobOutputDir=$DIR
jobName=Annovar-$(basename $DIR)

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
ANNOVARJOB=$jobOutputDir/$jobName\_Job.$today\.sh



echo "
##!/bin/sh
#$ -cwd		# use current working directory
#$ -o /data/scratch/$USER/
#$ -j y			# and put all output (inc errors) into it
#$ -m a			# Email on abort
#$ -pe smp 1		# Request 1 CPU cores
#$ -l h_rt=1:0:0	# Request 72 hour runtime
#$ -l h_vmem=5G		# Request 12G RAM / Core
#$ -t 1-$MAX		# run an array job of all the samples listed in FASTQ_Trim
#$ -N $jobName-Annovar
" > $ANNOVARJOB

echo '

module load gatk/4.1.6.0
module load annovar
module load bcftools
module load samtools

###################Download databases if not already done################################
# 
# #Check what databases are available
# annotate_variation.pl -webfrom annovar -downdb avdblist -buildver hg38 .
# 
# #FASTA sequences for all annotated transcripts in RefSeq Gene
# annotate_variation.pl -downdb -buildver hg38 -webfrom annovar refGene humandb/
# 
# #clinical interpretation of missense variants (indels not supported)
# annotate_variation.pl -buildver hg38 -downdb -webfrom annovar intervar_20180118 humandb/ 
# 
# #International Cancer Genome Consortium version 28
# annotate_variation.pl -buildver hg38 -downdb -webfrom annovar icgc28 humandb/ 
# 
# #CLINVAR database with Variant Clinical Significance and Variant disease name
# annotate_variation.pl -buildver hg38 -downdb -webfrom annovar clinvar_20220320 humandb/ 
# 
# #Database of Functional Predictions and Annotations for Human Nonsynonymous and Splice-Site SNVs
# annotate_variation.pl -buildver hg38 -downdb -webfrom annovar dbnsfp30a humandb/
#
# #dbSNP150 with allelic splitting and left-normalization
# annotate_variation.pl -buildver hg38 -downdb -webfrom annovar avsnp150 humandb/ 
#########################################################################################



Patients=($(ls -d /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/SH*))

## Extract the file name at the position of the array job task ID
Patient=$(basename ${Patients[${SGE_TASK_ID}-1]})

if [ "$Patient" = "" ]; then
	echo "No patient ID found, exiting..."
    exit 1
fi

echo "Running Annovar for sample: $Patient"


REF=/data/BCI-DigitalPath/Genome/GRCh38.p12.genome.fa


cd /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/
mkdir -p Final_vcfs
mkdir -p Final_vcfs/Results

#Move into sample directory
cd /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/$Patient


if [[ -s N ]] && [[ -s T ]]; then
	
	echo "Cleaning up tumour vs normal vcf for $Patient"
	bcftools annotate -x INFO,^FORMAT/GT,FORMAT/PL /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/${Patient}/${Patient}_Candidate_tumour_variants.vcf.gz > /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/${Patient}/Tidy_${Patient}_Candidate_tumour_variants.vcf
	
	gatk LeftAlignAndTrimVariants -V /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/${Patient}/Tidy_${Patient}_Candidate_tumour_variants.vcf -O /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/Final_vcfs/${Patient}_Candidate_tumour_variants.vcf -R $REF
	
	echo "Running bgzip on tumour vs normal vcf for $Patient"
	bgzip -f /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/Final_vcfs/${Patient}_Candidate_tumour_variants.vcf
	
	
	echo "Running annovar on tumour vs normal vcf for $Patient"
	table_annovar.pl /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/Final_vcfs/${Patient}_Candidate_tumour_variants.vcf.gz /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/humandb/ -buildver hg38 -out /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/Final_vcfs/Results/${Patient}_tumour_annovar -remove -protocol refGene,icgc28,clinvar_20220320,intervar_20180118,dbnsfp30a,avsnp150 -operation g,f,f,f,f,f -nastring . -vcfinput -polish
	
	
fi




if [[ -s N ]] && [[ -s M ]]; then

	echo "Cleaning up metastatic vs normal vcf for $Patient"
	bcftools annotate -x INFO,^FORMAT/GT,FORMAT/PL /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/${Patient}/${Patient}_Candidate_metastatic_variants.vcf.gz > /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/${Patient}/Tidy_${Patient}_Candidate_metastatic_variants.vcf
	
	gatk LeftAlignAndTrimVariants -V /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/${Patient}/Tidy_${Patient}_Candidate_metastatic_variants.vcf -O /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/Final_vcfs/${Patient}_Candidate_metastatic_variants.vcf -R $REF
	
	echo "Running bgzip on metastatic vs normal vcf for $Patient"
	bgzip -f /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/Final_vcfs/${Patient}_Candidate_metastatic_variants.vcf
	
	echo "Running annovar on metastatic vs normal vcf for $Patient"
	table_annovar.pl /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/Final_vcfs/${Patient}_Candidate_metastatic_variants.vcf.gz /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/humandb/ -buildver hg38 -out /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/Final_vcfs/Results/${Patient}_metastatic_annovar -remove -protocol refGene,icgc28,clinvar_20220320,intervar_20180118,dbnsfp30a,avsnp150 -operation g,f,f,f,f,f -nastring . -vcfinput -polish

fi



if ! [[ -s N ]] && [[ -s M ]]; then

	echo "Cleaning up germline metastatic vcf for $Patient"
	bcftools annotate -x INFO,^FORMAT/GT,FORMAT/PL /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/${Patient}/${Patient}_Candidate_metastatic_variants.vcf.gz > /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/${Patient}/Tidy_${Patient}_Candidate_metastatic_variants.vcf
	
	gatk LeftAlignAndTrimVariants -V /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/${Patient}/Tidy_${Patient}_Candidate_metastatic_variants.vcf -O /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/Final_vcfs/${Patient}_Candidate_metastatic_variants.vcf -R $REF
	
	echo "Running bgzip on germline metastatic vcf for $Patient"
	bgzip -f /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/Final_vcfs/${Patient}_Candidate_metastatic_variants.vcf
	
	echo "Running annovar on germline metastatic vcf for $Patient"
	table_annovar.pl /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/Final_vcfs/${Patient}_Candidate_metastatic_variants.vcf.gz /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/humandb/ -buildver hg38 -out /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/Final_vcfs/Results/${Patient}_metastatic_annovar -remove -protocol refGene,icgc28,clinvar_20220320,intervar_20180118,dbnsfp30a,avsnp150 -operation g,f,f,f,f,f -nastring . -vcfinput -polish

fi



if ! [[ -s N ]] && [[ -s T ]]; then

	echo "Cleaning up germline tumour vcf for $Patient"
	bcftools annotate -x INFO,^FORMAT/GT,FORMAT/PL /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/${Patient}/${Patient}_Candidate_tumour_variants.vcf.gz > /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/${Patient}/Tidy_${Patient}_Candidate_tumour_variants.vcf
	
	gatk LeftAlignAndTrimVariants -V /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/${Patient}/Tidy_${Patient}_Candidate_tumour_variants.vcf -O /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/Final_vcfs/${Patient}_Candidate_tumour_variants.vcf -R $REF
	
	echo "Running bgzip on germline tumour vcf for $Patient"
	bgzip -f /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/Final_vcfs/${Patient}_Candidate_tumour_variants.vcf
	
	echo "Running annovar on germline tumour vcf for $Patient"
	table_annovar.pl /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/Final_vcfs/${Patient}_Candidate_tumour_variants.vcf.gz /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/humandb/ -buildver hg38 -out /data/BCI-DigitalPath/Elements/nourse/CNourse_Trial/Patients/Final_vcfs/Results/${Patient}_tumour_annovar -remove -protocol refGene,icgc28,clinvar_20220320,intervar_20180118,dbnsfp30a,avsnp150 -operation g,f,f,f,f,f -nastring . -vcfinput -polish
	
fi





' >> $ANNOVARJOB
