#!/bin/bash
#$ -V
#$ -N TrimReads
#$ -m ea
#$ -o /datos/nagen/jobs/trimmomatic.stdout
#$ -e /datos/nagen/jobs/trimmomatic.stderr
#$ -w e
##$ -l h_vmem=2G
#$ -l m_core=20

# Trim reads and remove adapters with Trimmomatic

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
SampleID=$1

hostname >&2
echo $0 >&2
printf START >&2; uptime >&2
date >&2

# Run Trimmomatic for paired-end reads (PE)
/opt/jdk1.8.0_161/bin/java -jar /opt/trimmomatic/classes/trimmomatic.jar PE -threads 20 -phred33 reads/${SampleID}_1.fastq.gz reads/${SampleID}_2.fastq.gz trimmed_reads/${SampleID}_1_trimmed.fastq.gz trimmed_reads/${SampleID}_U1_trimmed.fastq.gz trimmed_reads/${SampleID}_2_trimmed.fastq.gz trimmed_reads/${SampleID}_U2_trimmed.fastq.gz ILLUMINACLIP:/opt/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

printf END >&2; uptime >&2

trap : 0

echo >&2 '
************
*** DONE *** 
************
'

