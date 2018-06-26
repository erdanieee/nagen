# Declare directories

GENOMES_DIR=/datos/nagen
NAGEN_SCRIPTS=/datos/nagen/nagen_scripts
NAGEN_SAMPLES=/datos/nagen/nagen_samples
GATK_RESOURCES=/datos/nagen/gatk_resources
GENOME_REFERENCE=/datos/nagen/reference_genome

# Genome reference
REFERENCE=hs37d5
REF_GENOME=${GENOME_REFERENCE}/${REFERENCE}.fa
REF_GENOME_INDEX=${GENOME_REFERENCE}/${REFERENCE}.fa.fai

# Index the reference for mapping with GEM3
if [ ! -e ${REF_GENOME_INDEX} ]
    then
    echo "Reference ${REF_GENOME} is not available."
    qsub -wd ${GENOME_REFERENCE} -N RefGen_${REFERENCE} $NAGEN_SCRIPTS/RefGenome.sh $REFERENCE
    qrsh -hold_jid RefGen_${REFERENCE} echo 
else
    echo "Reference $REFERENCE is ready to be used"
fi
