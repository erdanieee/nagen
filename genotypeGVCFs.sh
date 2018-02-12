#!/bin/bash
#$ -V
#$ -terse
#$ -m ea
#$ -M alberto.labarga.gutierrez@navarra.es
#$ -wd /datos/nagen
#$ -o /home/alabarga/jobs/opencga.stdout
#$ -e /home/alabarga/jobs/opencga.stderr
#$ -w e
##$ -l h_vmem=2G
#$ -l m_core=1

# Merge gVCFs using GATK's genotypeGVCFs tool
# Example of use sbatch genotypeGVCFs.sh E0001 /mnt/lustre/scratch/CBRA/projects/NAGEN/analysis/ /mnt/lustre/scratch/CBRA/data/indexed_genomes/bwa/hs37d5/hs37d5.fa 

hostname >&2
echo $0 >&2
printf START >&2; uptime >&2
date >&2

familyId=$1
familyInfo=$2
outputDir=$3
referenceGenome=/datos/nagen/reference_genome/hs37d5.fa


variantArguments=''

for sample in $(cat $familyInfo);do
variantArguments=$variantArguments' --variant '$outputDir/$sample'.g.vcf.gz ';done

##for sample in $(ls $outputDir/AC*.vcf.gz);do
##variantArguments=$variantArguments' --variant '$sample' ';done

# CombineGVCFs
/opt/gatk-4.0.0.0/gatk CombineGVCFs \
$variantArguments \
-O $outputDir/$familyId.vcf.gz \
-R $referenceGenome

# GenotypeGVCFs
/opt/gatk-4.0.0.0/gatk GenotypeGVCFs \
-V $outputDir/$familyId.vcf.gz \
-O $outputDir/$familyId.2.vcf.gz \
-R $referenceGenome

printf END >&2; uptime >&2
