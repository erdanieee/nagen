#!/bin/bash

###################################
## NaGen Variant Calling Pipeline #
###################################
##
## This bash pipeline runs FastQC, Trimmomatic, BWA mapping, picard, samtools
## & various GATK4 tools on a single sample, scattered across chromosomes.
## Note: Pipeline NagenFamily should be used for analysing family members.
##
## Requirements/expectations :
## - Name of sample to analyse (SampleID) = Two ready fastq files
## - Output directory (OutputDir) for writing the results
## - Reference genome (Reference)
##
## Outputs :
## - Quality Control reports
## - Trimmed reads in fastq format
## - Bam processed file
## - vcf file for each sample
## - Statistics reports
##
## Runtime parameters should be optimized for Nasertic Platform implementation.#

# Declare directories
NaGen_Scripts=/datos/nagen/nagen_scripts
RefGenome=/datos/nagen/reference_genome
GATKResources=/datos/nagen/GATKResources
NaGenWD=/datos/nagen

# Check 3 input parameters have been parsed
if [ $# -ne 3 ]
	then
	echo "NaGen mapping reads and variant calling pipeline"
	echo "USAGE: nagen.sh <SampleID> <OutputDir> <Reference>"
	exit 1
fi

# Input parameters
# Take sampleID for fastq reads
SampleID=$1
echo "Input SampleID name: $SampleID"

# Take output directory
OutputDir=$2
if [ ! -d "$OutputDir" ]; then
 	mkdir $OutputDir
	mkdir $OutputDir/tmp
fi
echo "Output directory: $OutputDir"

# Take reference genome name
Reference=$3
echo "Reference genome: $Reference"
echo


############################
# Prepare Reference Genone #
############################

# Index the reference for mapping with BWA
# If reference genome index does not exist, produce it
if [ ! -e "$RefGenome/"$Reference".fa.fai" ]
	then
	echo $Reference" Reference genome is not available."
	cd $RefGenome
	qsub -N 'RefGen_'$Reference $NaGen_Scripts/RefGenome.sh
	qrsh -hold_jid 'RefGen_'$Reference echo 
else
	echo "Reference genome is ready to be used"
fi


# Move to output directory
cd $OutputDir


########################
# Prepare reads folder #
########################

# Make symbolic links for sequencing reads
echo "Preparing reads folder - Making symbolic links for reads" 
mkdir reads
ln -s $NaGenWD/original_reads/$SampleID.*.fastq reads/


###################
# Quality control #
###################

# Perform quality control of reads
echo "Performing quality control using FastQC"
mkdir fastqc
qsub -N 'QC_'$SampleID $NaGen_Scripts/QualControl.sh $SampleID
qrsh -hold_jid 'QC_'$SampleID echo "QC is done"
echo


##############
# Trim reads #
##############

# Remove Illumina adapters and trim reads with quality below threshold defined
echo "Trimming Illumina adapters and reads will low quality"
mkdir trimmed_reads
qsub -N 'trim_'$SampleID $NaGen_Scripts/TrimReads.sh $SampleID
qrsh -hold_jid 'trim_'$SampleID echo "Trimming is done"
mv *trimmed.fastq.gz trimmed_reads/
echo


#####################
# Mapping using BWA #
#####################

# Map paired-end and single-end reads and add read group info
echo "Mapping paired-end and single-end reads"
qsub -N 'bwa_'$SampleID $NaGen_Scripts/MappingBWA.sh $SampleID $Reference
qrsh -hold_jid 'bwa_'$SampleID echo "Mapping is done"
echo


###############################################
# SAM to BAM, sort BAM and obtains statistics #
###############################################

# Merge SAM files into a single sorted BAM file and produce statistics
echo "Merge SAM files and sort by coordinate into a BAM file"
qsub -N 'sam2bam_'$SampleID $NaGen_Scripts/Sam2Bam.sh $SampleID
qrsh -hold_jid 'sam2bam_'$SampleID echo "Sorted bam files produced"
echo
rm $SampleID*.sam


#####################################
# Removing duplicates and filtering #
#####################################

# Marking and removing duplicates
echo "Removing duplicates and applying further filtering"
qsub -N 'RmDup_'$SampleID $NaGen_Scripts/BamRmDup.sh $SampleID
qrsh -hold_jid 'RmDup_'$SampleID echo "Duplicates marked and removed"
mkdir bam
mv $SampleID.sorted.bam $SampleID.sorted.bam.bai $SampleID.sorted.flagstat.txt bam/
rm $SampleID.sorted.rmdup.bam $SampleID.sorted.rmdup.bai $SampleID.sorted.rmdup.bam.bai
echo


######################
# Base Recalibration #
######################

#Run Recalibration of base quality scores
echo "Running recalibration of base quality scores"
for chr in `cat /datos/nagen/nagen_scripts/chromosomes.txt`
do
	qsub -N 'bqsr_'$chr $NaGen_Scripts/BaseRecal.sh $SampleID $chr $Reference
done
Wait for chromosomes 1 and 4 to finish and continue
qrsh -hold_jid 'bqsr_1' echo "Recalibration is done"
qrsh -hold_jid 'bqsr_4' echo "Recalibration is done"
#Remove intermediate files
rm $SampleID.sorted.rmdup.filtered.bam $SampleID.sorted.rmdup.filtered.bam.bai
echo


###########################
# Statistics on BAM files #
###########################

# Run some summary statistics on the bam files
echo "Collecting alignment summary metrics for BAM files"
for chr in `cat /datos/nagen/nagen_scripts/chromosomes.txt`
do
	qsub -N 'bamstat_'$chr $NaGen_Scripts/BamStats.sh $SampleID $chr $Reference
done
qrsh -hold_jid 'bamstat_22' echo "Statistics have been obtained"
echo


#####################################
# Variant Calling: Haplotype caller #
#####################################

# Run HaplotypeCaller
echo "Running HaplotypeCaller on chromosomes"
for chr in `cat /datos/nagen/nagen_scripts/chromosomes.txt`
do
	qsub -N 'HapCall'_$chr $NaGen_Scripts/HaplotypeCaller.sh $SampleID $chr $Reference
done
qrsh -hold_jid 'HapCall_1' echo "Done"
qrsh -hold_jid 'HapCall_2' echo "Done"
echo
# Organize files
mkdir stats
mv *data.table *metrics.txt stats/
mv $SampleID'.sorted.rmdup.filtered.recalibrated'* bam/
mkdir vcfs


##############
# Genotyping #
##############

# Run GenotypeGVCFs
echo "Running Genotyping"
for chr in `cat /datos/nagen/nagen_scripts/chromosomes.txt`
do
    qsub -N 'Genotype_'$chr $NaGen_Scripts/GenotypeGVCFInd.sh $SampleID $chr
done
qrsh -hold_jid 'Genotype_1' echo "Done"
mkdir vcfs/gvcfs
mv *.g.vcf.gz *.g.vcf.gz.tbi vcfs/gvcfs/
echo


###################
# Merge VCF files #
###################

# Run MergeVcfs
echo "Merge the vcfs obtained for individual chromosomes"
qsub -N 'MergeVCFs_'$SampleID $NaGen_Scripts/MergeVCFs.sh $SampleID
qrsh -hold_jid 'MergeVCFs_'$SampleID echo "Done"
mkdir vcfs/vcfs
mv *.vcf.gz *.vcf.gz.tbi vcfs/vcfs
mv vcfs/vcfs/$SampleID.vcf.gz vcfs/vcfs/$SampleID.vcf.gz.tbi .
echo


########################
# APPLY HARD FILTERING #
########################

# Run variant selection and variant filtration
echo "Applying hard filterings"
qsub -N 'HardFilt_'$SampleID $NaGen_Scripts/HardFilters.sh $SampleID
qrsh -hold_jid 'HardFilt_'$SampleID echo "Done"
echo


#################################
# Merge hard filtered VCF files #
#################################

# Run MergeVcfs
echo "Merging the vcfs for SNPs and InDels"
qsub -N 'MergeHFVCFs_'$SampleID $NaGen_Scripts/MergeHFVCFs.sh $SampleID
qrsh -hold_jid 'MergeHFVCFs_'$SampleID echo "Done"
echo
mv *.vcf.gz *.vcf.gz.tbi vcfs/

rm -r tmp/
