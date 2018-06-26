#!/bin/bash

##################################
# NAGEN VARIANT CALLING PIPELINE #
##################################

# Declare directories
NAGEN_SCRIPTS=/datos/nagen/nagen_scripts
REF_GENOME=/datos/nagen/reference_genome
GATK_RESOURCES=/datos/nagen/gatk_resources
GENOMES_DIR=/datos/nagen


# The pipeline requires name of genome to analyse (GENOME), 
# output directory (OUTPUT_DIR) and reference genome (REFERENCE)

# Check input parameters
if [ $# -ne 3 ]
	then
	echo "NaGen mapping reads and variant calling pipeline"
	echo "USAGE: nagen.sh <GENOME> <OUTPUT_DIR> <REFERENCE>"
	exit 1
fi

# Input parameters
GENOME=$1
echo "Input GENOME name: $GENOME"

OUTPUT_DIR=$2
if [ ! -d "$OUTPUT_DIR" ]; then
 	mkdir $OUTPUT_DIR
	mkdir $OUTPUT_DIR/tmp
	mkdir $OUTPUT_DIR/reads
	ln -s $GENOMES_DIR/original_reads/$GENOME.*fastq $OUTPUT_DIR/reads/
fi
echo "Output directory: $OUTPUT_DIR"

REFERENCE=$3
echo "Reference genome: $REFERENCE"
echo


############################
# PREPARE REFERENCE GENOME #
############################

# Index the reference for mapping with GEM3
if [ ! -e "$REF_GENOME/"$REFERENCE".fa.fai" ]
	then
	echo $REFERENCE" Reference genome is not available."
	cd $REF_GENOME
	qsub -N 'RefGen_'$REFERENCE $NAGEN_SCRIPTS/RefGenome.sh
	qrsh -hold_jid 'RefGen_'$REFERENCE echo 
else
	echo "Reference genome is ready to be used"
fi


# Move to output directory
cd $OUTPUT_DIR


#########################
# PREPARE READS FOLDER  #
#########################

# Make a symbolic link of reads
mkdir reads
ln -s $GENOMES_DIR/original_reads/$GENOME.*.fastq reads/


###################
# QUALITY CONTROL #
###################

# Perform quality control of reads
mkdir fastqc
qsub -N QC_$GENOME $NAGEN_SCRIPTS/QualControl.sh $GENOME
qrsh -hold_jid QC_$GENOME echo "QC is done"


##############
# TRIM READS #
##############

# Remove adapters and trim reads with low quality
mkdir trimmed_reads
qsub -N trim_$GENOME $NAGEN_SCRIPTS/TrimReads.sh $GENOME
qrsh -hold_jid trim_$GENOME echo "Trimming is done"
mv *trimmed.fastq.gz trimmed_reads/


#####################
# MAPPING USING BWA #
#####################

# Map paired-end and single-end reads and aad read group info.
qsub -N bwa_$GENOME $NAGEN_SCRIPTS/MappingBWA.sh $GENOME $REFERENCE
qrsh -hold_jid bwa_$GENOME echo "Mapping is done"


################################################
# SAM TO BAM, SORT BAM AND SAMTOOLS STATISTICS #
################################################

# Merge SAM files into a single sorted BAM file and produce statistics
qsub -N sam2bam_$GENOME $NAGEN_SCRIPTS/Sam2Bam.sh $GENOME
qrsh -hold_jid sam2bam_$GENOME echo "Sorted bam files produced"

rm $GENOME*.sam

####################################
# MARKING DUPLICATES AND FILTERING #
####################################
 
# Marking and removing duplicates
qsub -N RmDup_$GENOME $NAGEN_SCRIPTS/BamRmDup.sh $GENOME
qrsh -hold_jid RmDup_$GENOME echo "Duplicates marked and removed"
mkdir bam
mv $GENOME.sorted.bam $GENOME.sorted.bam.bai $GENOME.sorted.flagstat.txt bam/

######################
# BASE RECALIBRATION #
######################

#Run Recalibration of base quality scores
for chr in `cat /datos/nagen/nagen_scripts/chromosomes.txt`
do
	qsub -N bqsr_$chr $NAGEN_SCRIPTS/BaseRecal.sh $GENOME $chr
done


###########################
# STATISTICS ON BAM FILES #
###########################

# Run some summary statistics on the bam files
for chr in `cat /datos/nagen/nagen_scripts/chromosomes.txt`
do
	qsub -N bamstat_$chr $NAGEN_SCRIPTS/BamStats.sh $GENOME $chr $REFERENCE
done
qsub -hold_jid


################################################
# VARIANT CALLING: HAPLOTYPE CALLER, GVCF FILE #
################################################

# Run HaplotypeCaller
for chr in `cat /datos/nagen/nagen_scripts/chromosomes.txt`
do
	qsub $NAGEN_SCRIPTS/HaplotypeCaller.sh $GENOME $chr
done
mkdir gvcf


##############
# GENOTYPING #
##############

# Run GenotypeGVCFs
for chr in `cat ../nagen_scripts/chromosomes.txt `
do
    qsub $NAGEN_SCRIPTS/GenotypeGVCFInd.sh $GENOME $chr
done


###################
# MERGE VCF files #
###################

# Run MergeVcfs
qsub $NAGEN_SCRIPTS/MergeVCFs.sh $GENOME


########################
# APPLY HARD FILTERING #
########################

# Run variant selection and variant filtration
qsub $NAGEN_SCRIPTS/HardFilters.sh $GENOME

