#!/bin/bash
#$ -V
#$ -m ea
#$ -M alberto.labarga.gutierrez@navarra.es
#$ -wd /datos/nagen/nagen_samples
#$ -o /home/alabarga/jobs/opencga.stdout
#$ -e /home/alabarga/jobs/opencga.stderr
#$ -w e
##$ -l h_vmem=2G
#$ -l m_core=1

### Cat variants for the set of g.VCFs for a given sample
## Example: sbatch cat_gVCFs.sh AC5421 E0020 /mnt/lustre/scratch/CBRA/projects/NAGEN/NAGEN_01/gVCF/ /mnt/lustre/scratch/CBRA/projects/NAGEN/analysis/ /mnt/lustre/scratch/CBRA/data/indexed_genomes/bwa/hs37d5/hs37d5.fa 

sampleId=$1
inputDir=gvcf
outputDir=gvcf
referenceGenome=/datos/nagen/reference_genome/hs37d5.fa


# Cat gVCFs
/opt/gatk-4.0.0.0/gatk MergeVcfs \
-R $referenceGenome \
-I $inputDir/$sampleId'.1.g.vcf.gz' \
-I $inputDir/$sampleId'.2.g.vcf.gz' \
-I $inputDir/$sampleId'.3.g.vcf.gz' \
-I $inputDir/$sampleId'.4.g.vcf.gz' \
-I $inputDir/$sampleId'.5.g.vcf.gz' \
-I $inputDir/$sampleId'.6.g.vcf.gz' \
-I $inputDir/$sampleId'.7.g.vcf.gz' \
-I $inputDir/$sampleId'.8.g.vcf.gz' \
-I $inputDir/$sampleId'.9.g.vcf.gz' \
-I $inputDir/$sampleId'.10.g.vcf.gz' \
-I $inputDir/$sampleId'.11.g.vcf.gz' \
-I $inputDir/$sampleId'.12.g.vcf.gz' \
-I $inputDir/$sampleId'.13.g.vcf.gz' \
-I $inputDir/$sampleId'.14.g.vcf.gz' \
-I $inputDir/$sampleId'.15.g.vcf.gz' \
-I $inputDir/$sampleId'.16.g.vcf.gz' \
-I $inputDir/$sampleId'.17.g.vcf.gz' \
-I $inputDir/$sampleId'.18.g.vcf.gz' \
-I $inputDir/$sampleId'.19.g.vcf.gz' \
-I $inputDir/$sampleId'.20.g.vcf.gz' \
-I $inputDir/$sampleId'.21.g.vcf.gz' \
-I $inputDir/$sampleId'.22.g.vcf.gz' \
-I $inputDir/$sampleId'.X.g.vcf.gz' \
-I $inputDir/$sampleId'.Y.g.vcf.gz' \
-O $outputDir/$sampleId'.g.vcf.gz'
