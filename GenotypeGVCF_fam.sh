#!/bin/bash
#$ -V
#$ -m ea
#$ -M alberto.labarga.gutierrez@navarra.es
#$ -o /home/alabarga/jobs/genotype_fam.stdout
#$ -e /home/alabarga/jobs/genotype_fam.stderr
#$ -w e
##$ -l h_vmem=2G
##$ -l m_core=1

# Declare directories
GENOMES_DIR=/datos/nagen
NAGEN_SCRIPTS=/datos/nagen/nagen_scripts
NAGEN_SAMPLES=/datos/nagen/nagen_samples
REF_GENOME=/datos/nagen/reference_genome
GATK_RESOURCES=/datos/nagen/gatk_resources
NAGEN_FAMILIES=/datos/nagen/nagen_families

# Merge gVCFs using GATK's genotypeGVCFs tool
# Example of use sbatch genotypeGVCFs.sh E0001 /mnt/lustre/scratch/CBRA/projects/NAGEN/analysis/ /mnt/lustre/scratch/CBRA/data/indexed_genomes/bwa/hs37d5/hs37d5.fa 

hostname >&2
echo $0 >&2
printf START >&2; uptime >&2
date >&2

familyId=NAGEN_02
familyInfo=sample.list

referenceGenome=/datos/nagen/reference_genome/hs37d5.fa

variantArguments=''

for sample in $(cat $familyInfo);do
 
 variantArguments=${variantArguments}' --variant '${NAGEN_SAMPLES}/${sample}/gvcf/${sample}.g.vcf.gz

done

# CombineGVCFs
/opt/gatk-4.0.0.0/gatk CombineGVCFs \
$variantArguments \
-O $familyId.g.vcf.gz \
-R $referenceGenome

# GenotypeGVCFs
/opt/gatk-4.0.0.0/gatk GenotypeGVCFs \
-V $familyId.g.vcf.gz \
-O $familyId.vcf.gz \
-R $referenceGenome

printf END >&2; uptime >&2
