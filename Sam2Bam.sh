#!/bin/bash
#$ -V
#$ -N Sam2Bam
#$ -m ea
#$ -o /datos/nagen/jobs/sam2bam.stdout
#$ -e /datos/nagen/jobs/sam2bam.stderr
#$ -w e
##$ -l h_vmem=2G
#$ -l m_core=4

# Sam to bam, sort bam and produce samtools statistics

abort()
{
    echo >&2 '
***************
*** ABORTED ***
***************
'
    echo "An error occurred. Exiting..." >&2
    exit 1
}

set -e
trap 'abort' 0

# Get sample name
SAMPLE_NAME=$1
SampleID=$2
REFERENCE=$3

hostname >&2
echo $0 >&2
printf START >&2; uptime >&2
date >&2

if [ ! -d "bam" ]; then
    mkdir bam
fi

# Merge SAM files into a single coordinate-sorted BAM file
/opt/jdk1.8.0_161/bin/java -Xmx2g -Djava.io.tmpdir=`pwd`/tmp -jar /opt/picard/build/libs/picard.jar \
		MergeSamFiles \
		I=sam/${SampleID}_pe.sam \
		I=sam/${SampleID}_se1.sam \
		I=sam/${SampleID}_se2.sam \
		O=bam/${SampleID}.sorted.bam \
		SO=coordinate \
		TMP_DIR=`pwd`/tmp

# Index BAM file
/opt/samtools-1.6/bin/samtools index bam/${SampleID}.sorted.bam

# Run samtools flagstats
/opt/samtools-1.6/bin/samtools flagstat bam/${SampleID}.sorted.bam > bam/${SampleID}.sorted.flagstat.txt

printf END >&2; uptime >&2

trap : 0

echo >&2 '
************
*** DONE *** 
************
'
