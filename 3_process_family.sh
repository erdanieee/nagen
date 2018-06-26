#!/bin/bash

# Declare directories
GENOMES_DIR=/datos/nagen
NAGEN_SCRIPTS=/datos/nagen/nagen_scripts
NAGEN_SAMPLES=/datos/nagen/nagen_samples
REF_GENOME=/datos/nagen/reference_genome
GATK_RESOURCES=/datos/nagen/gatk_resources
NAGEN_FAMILIES=/datos/nagen/nagen_families

REFERENCE=hs37d5

FAMILY_ID=$1

BASEDIR=${NAGEN_FAMILIES}/${FAMILY_ID}

/opt/opencga/bin/opencga.sh users login -u alabarga

qsub -wd ${BASEDIR} -N Genotype_Fam_${FAMILY_ID} ${NAGEN_SCRIPTS}/GenotypeGVCF_fam.sh ${FAMILY_ID}
qsub -wd $BASEDIR -hold_jid Genotype_Fam_${FAMILY_ID} -N FilterFam_${FAMILY_ID} $NAGEN_SCRIPTS/HardFilters_fam.sh $FAMILY_ID
qsub -wd $BASEDIR -hold_jid FilterFam_$FAMILY_ID -N LoadOpenCGA_${FAMILY_ID} $NAGEN_SCRIPTS/LoadOpenCGA.sh $FAMILY_ID

