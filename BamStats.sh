#!/bin/bash
#$ -V
#$ -N BamStats
#$ -m ea
#$ -cwd
#$ -o /datos/nagen/jobs/bamstats.stdout
#$ -e /datos/nagen/jobs/bamstats.stderr
#$ -w e
##$ -l h_vmem=2G
#$ -l m_core=4

# Statistics on pre-processed bam files

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

# Get name of sample, chromosome and refence genome
SampleID=$1
Chromosome=$2
Reference=$3

hostname >&2
echo $0 >&2
printf START >&2; uptime >&2
date >&2

# Run some statistics
/opt/jdk1.8.0_161/bin/java -Xmx18000m -Djava.io.tmpdir=`pwd`/tmp -jar /opt/picard/build/libs/picard.jar \
		CollectAlignmentSummaryMetrics \
		I=recalibration/$SampleID.sorted.rmdup.filtered.recalibrated.$Chromosome.bam \
		O=stats/$SampleID.sorted.rmdup.filtered.recalibrated.$Chromosome.metrics.txt \
		REFERENCE_SEQUENCE=/datos/nagen/reference_genome/$Reference.fa \
		VALIDATION_STRINGENCY=LENIENT

# Index BAM file
/opt/samtools-1.6/bin/samtools flagstat recalibration/$SampleID.sorted.rmdup.filtered.recalibrated.$Chromosome.bam > stats/$SampleID.sorted.rmdup.filtered.recalibrated.$Chromosome.flagstat.txt

printf END >&2; uptime >&2

trap : 0

echo >&2 '
************
*** DONE ***
************
'
