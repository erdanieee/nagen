#!/bin/bash

##################################
# NAGEN VARIANT CALLING PIPELINE #
##################################

# Declare directories
NAGEN_SCRIPTS=/datos/nagen/nagen_scripts
REF_GENOME=/datos/nagen/reference_genome
GATK_RESOURCES=/datos/nagen/gatk_resources
GENOMES_DIR=/datos/nagen
BASEDIR=/datos/nagen/nagen_samples/

# The pipeline requires name of genome to analyse (SAMPLE_ID), 
# output directory (OUTPUT_DIR) and reference genome (REFERENCE)

# Check input parameters
if [ $# -ne 2 ]
	then
	echo "NaGen mapping reads and variant calling pipeline"
	echo "USAGE: nagen.sh <SAMPLE_ID> <REFERENCE>"
	exit 1
fi

# Input parameters
SAMPLE_ID=$1
REFERENCE=$2
WORKDIR=${BASEDIR}/${SAMPLE_ID}

cromosomas=$(cat /datos/nagen/nagen_scripts/chromosomes.txt)

######################
# BASE RECALIBRATION #
######################

#Run Recalibration of base quality scores
for chr in ${cromosomas}
do
    echo bqsr_${SAMPLE_ID}_$chr
	qsub -N bqsr_${SAMPLE_ID}_$chr -hold_jid RmDup_${SAMPLE_ID} $NAGEN_SCRIPTS/BaseRecal.sh $SAMPLE_ID $chr $REFERENCE
done


###########################
# STATISTICS ON BAM FILES #
###########################

# Run some summary statistics on the bam files
for chr in ${cromosomas}
do
    echo bamstat_${SAMPLE_ID}_$chr 
	qsub -wd $WORKDIR -hold_jid bqsr_${SAMPLE_ID}_$chr -N bamstat_${SAMPLE_ID}_$chr $NAGEN_SCRIPTS/BamStats.sh $SAMPLE_ID $chr $REFERENCE
done


################################################
# VARIANT CALLING: HAPLOTYPE CALLER, GVCF FILE #
################################################

# Run HaplotypeCaller

for chr in ${cromosomas}
do
    #qsub -wd $BASEDIR -N RenameSample_${GENOME}_${chr} $NAGEN_SCRIPTS/RenameSample.sh $GENOME \
    #recalibration/${GENOME}.sorted.rmdup.filtered.recalibrated.$chr.bam \
    #recalibration/${GENOME}.$chr.bam
    echo HaplotypeCaller_${SAMPLE_ID}_${chr}
	qsub -wd $WORKDIR -hold_jid bqsr_${SAMPLE_ID}_$chr -N HaplotypeCaller_${SAMPLE_ID}_${chr} $NAGEN_SCRIPTS/HaplotypeCaller.sh $SAMPLE_ID $chr $REFERENCE
done

##############
# GENOTYPING #
##############

# Run GenotypeGVCFs

JOB_IDS=()
for chr in ${cromosomas}
do	
	JOB_IDS+=(Genotype_${chr})
    qsub -wd $WORKDIR -hold_jid HaplotypeCaller_${SAMPLE_ID}_${chr} -N Genotype_${chr} $NAGEN_SCRIPTS/GenotypeGVCF.sh $SAMPLE_ID $chr
done

job_ids=$(IFS=","; echo "${JOB_IDS[*]}" )

###################
# MERGE VCF files #
###################

# Run MergeVcfs
qsub -wd $WORKDIR -hold_jid $job_ids -N MergeVcfs_${SAMPLE_ID} $NAGEN_SCRIPTS/MergeVCFs.sh $SAMPLE_ID

# Run cat_gVCFs
# qsub -wd $WORKDIR -hold_jid $job_ids -N cat_gVCFs_${SAMPLE_ID} $NAGEN_SCRIPTS/cat_gVCFs.sh $SAMPLE_ID

########################
# APPLY HARD FILTERING #
########################

# Run variant selection and variant filtration
qsub -wd $WORKDIR -hold_jid MergeVcfs_${SAMPLE_ID},cat_gVCFs_${SAMPLE_ID} -N HardFilters_${SAMPLE_ID} $NAGEN_SCRIPTS/HardFilters.sh $SAMPLE_ID

# Run variant recalibration
qsub -wd $WORKDIR -hold_jid MergeVcfs_${SAMPLE_ID},cat_gVCFs_${SAMPLE_ID} -N VariantRecalibration_${SAMPLE_ID} $NAGEN_SCRIPTS/VariantRecalibration.sh $SAMPLE_ID
