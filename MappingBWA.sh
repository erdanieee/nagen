#!/bin/bash
#$ -V
#$ -N BWAmapping
#$ -m ea
#$ -o /datos/nagen/jobs/bwamapping.stdout
#$ -e /datos/nagen/jobs/bwamapping.stderr
#$ -w e
##$ -l h_vmem=2G
#$ -l m_core=20

# Map reads to reference genome using BWA

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

# Get reference genome and sample ID

SAMPLE_NAME=$1
SEQ_UNIT=$2
Reference=$3

hostname >&2
echo $0 >&2
printf START >&2; uptime >&2
date >&2

if [ ! -d "sam" ]; then
    mkdir sam
fi

# Map paired-end reads and add read group info.

/opt/bwa/bwa mem -t 20 /datos/nagen/reference_genome/$Reference.fa  \
                       trimmed_reads/${SEQ_UNIT}_1_trimmed.fastq.gz \
                       trimmed_reads/${SEQ_UNIT}_2_trimmed.fastq.gz \
                       -R '@RG\tID:'${SAMPLE_NAME}'\tSM:'${SAMPLE_NAME}'\tPL:Illumina\tCN:CBRA\tLB:Fragment' \
                       > sam/${SEQ_UNIT}_pe.sam

# Map single-end reads and add read group info.

/opt/bwa/bwa mem -t 10 /datos/nagen/reference_genome/$Reference.fa \
                       trimmed_reads/${SEQ_UNIT}_U1_trimmed.fastq.gz \
                       -R '@RG\tID:${SAMPLE_NAME}\tSM:${SAMPLE_NAME}\tPL:Illumina\tCN:CBRA\tLB:Fragment\tPU:${SAMPLE_NAME}' \
                       > sam/${SEQ_UNIT}_se1.sam


# /opt/bwa/bwa mem -t 10 /datos/nagen/reference_genome/hs37d5.fa trimmed_reads/JC5FBBXX_2_27nf_U1_trimmed.fastq.gz -R '@RG\tID:AC5444\tSM:AC5444\tPL:Illumina\tCN:CBRA\tLB:Fragment'  > sam/JC5FBBXX_2_27nf_se1.sam

/opt/bwa/bwa mem -t 10 /datos/nagen/reference_genome/$Reference.fa \
                       trimmed_reads/${SEQ_UNIT}_U2_trimmed.fastq.gz \
                       -R '@RG\tID:'${SAMPLE_NAME}'\tSM:'${SAMPLE_NAME}'\tPL:Illumina\tCN:CBRA\tLB:Fragment' \
                       > sam/${SEQ_UNIT}_se2.sam                      

#-R '@RG\tID:'$SampleID'\tSM:'$SampleID'\tPL:Illumina\tCN:CBRA\tLB:Fragment' > $SampleID'_se2'.sam

printf END >&2; uptime >&2

trap : 0

echo >&2 '
************************
*** DONE             *** 
************************
'

