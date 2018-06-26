#!/bin/bash
#$ -V
#$ -N QualControl
#$ -m ea
#$ -o /datos/nagen/jobs/fastqc.stdout
#$ -e /datos/nagen/jobs/fastqc.stderr
#$ -w e
##$ -l h_vmem=2G
#$ -l m_core=10

# Perform quality control of sequencing reads

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

hostname >&2
echo $0 >&2
printf START >&2; uptime >&2
date >&2

if [ ! -d "fastqc" ]; then
    mkdir fastqc
fi

# Run quality control using FastQC
/opt/FastQC/fastqc -o fastqc/ -t 10 reads/${SampleID}_1.fastq.gz reads/${SampleID}_2.fastq.gz

printf END >&2; uptime >&2

trap : 0

echo >&2 '
************
*** DONE *** 
************
'
