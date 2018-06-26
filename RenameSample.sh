#!/bin/bash
#$ -V
#$ -N BamRmDup
#$ -m ea
#$ -cwd
#$ -o /datos/nagen/jobs/haplotypecaller.stdout
#$ -e /datos/nagen/jobs/haplotypecaller.stderr
#$ -w e
#$ -l h_vmem=16G
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
SAMPLE_NAME=$1
INPUT=$2
OUTPUT=$3

hostname >&2
echo $0 >&2
printf START >&2; uptime >&2
date >&2


/opt/jdk1.8.0_161/bin/java -jar -Xmx2g -Djava.io.tmpdir=`pwd`/tmp -jar /opt/picard/build/libs/picard.jar \
        AddOrReplaceReadGroups \
        I=${INPUT} ID=${SAMPLE_NAME} SM=${SAMPLE_NAME} PL=${SAMPLE_NAME} LB=Fragment PU=${SAMPLE_NAME} \
        O=${OUTPUT}


printf END >&2; uptime >&2

trap : 0

echo >&2 '
************
*** DONE *** 
************
'
