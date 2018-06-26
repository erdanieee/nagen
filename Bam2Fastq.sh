#!/bin/bash
#$ -V
#$ -N bam2fastq
#$ -m ea
#$ -wd /datos/nagen
#$ -o /datos/nagen/jobs/bam2fastq.stdout
#$ -e /datos/nagen/jobs/bam2fastq.stderr
#$ -w e
##$ -l h_vmem=2G
#$ -l m_core=2

# Convert bam to fastq

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

# Get sequencing reads in bam format and the name of genome/sample
SeqBam=$1
SampleID=$2

hostname >&2
echo $0 >&2
printf START >&2; uptime >&2
date >&2

# First sort the bam file by coordinate
/opt/samtools-1.6/bin/samtools sort -n $SeqBam $SampleID.tmp.qsort.bam

# Use bamtofastq to convert the bam file into two fastq files (paired-end reads)
/opt/bib/active/bin/bamToFastq -i $SampleID.tmp.qsort.bam -fq $SampleID.1.fastq -fq2 $SampleID.2.fastq

rm $SampleID $SeqBam.qsort.bam

printf END >&2; uptime >&2

trap : 0

echo >&2 '
************
*** DONE ***
************
'
