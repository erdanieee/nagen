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

hostname >&2
echo $0 >&2
printf START >&2; uptime >&2
date >&2

familyId=$1
outputDir=$2
referenceGenome=/datos/nagen/reference_genome/hs37d5.fa

# GenotypeGVCFs
/opt/gatk-4.0.0.0/gatk GenotypeGVCFs \
-V $outputDir/$familyId.g.vcf.gz \
-O $outputDir/$familyId.vcf.gz \
-R $referenceGenome

printf END >&2; uptime >&2
