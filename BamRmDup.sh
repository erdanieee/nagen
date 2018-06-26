#!/bin/bash
#$ -V
#$ -N BamRmDup
#$ -m ea
#$ -cwd
#$ -o /datos/nagen/jobs/bamrmdup.stdout
#$ -e /datos/nagen/jobs/bamrmdup.stderr
#$ -w e
#$ -l h_vmem=8G
#$ -l m_core=10

# Marking and removing duplicates

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

# Get genome name
SampleID=$1
#INPUT_FILES=$2

hostname >&2
echo $0 >&2
printf START >&2; uptime >&2
date >&2

bam_array=$(ls bam/*.sorted.bam | awk '{print "INPUT="$1}')

# Mark and remove duplicates
/opt/jdk1.8.0_161/bin/java -jar -Xmx2g -Djava.io.tmpdir=`pwd`/tmp -jar /opt/picard/build/libs/picard.jar \
		MarkDuplicates \
		$(echo ${bam_array[@]}) \
		OUTPUT=bam/${SampleID}.sorted.rmdup.bam \
		METRICS_FILE=bam/${SampleID}.sorted.rmdup.metrics.txt \
		REMOVE_DUPLICATES=true \
		ASSUME_SORTED=true \
		CREATE_INDEX=true

# Index sorted rmdup BAM file
/opt/samtools-1.6/bin/samtools index bam/${SampleID}.sorted.rmdup.bam

# Run some BAM filtering (reads unmapped and mapq quality below 10)
/opt/samtools-1.6/bin/samtools view -q 10 -F 4 bam/${SampleID}.sorted.rmdup.bam -o bam/${SampleID}.sorted.rmdup.filtered.bam

# Index filtered bam file
/opt/samtools-1.6/bin/samtools index bam/${SampleID}.sorted.rmdup.filtered.bam

printf END >&2; uptime >&2

trap : 0

echo >&2 '
************
*** DONE *** 
************
'
