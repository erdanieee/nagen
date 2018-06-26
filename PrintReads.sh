#!/bin/bash
#$ -V
#$ -N BamRmDup
#$ -m ea
#$ -cwd
#$ -o /datos/nagen/jobs/bamrmdup.stdout
#$ -e /datos/nagen/jobs/bamrmdup.stderr
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
SampleID=$1
Reference=$2

hostname >&2
echo $0 >&2
printf START >&2; uptime >&2
date >&2

bam_array=$(ls recalibration/*.bam | awk '{print "-I "$1}')
echo ${bam_array[@]}

/opt/jdk1.8.0_161/bin/java -d64 -Xms4000m -Xmx4000m -Djava.io.tmpdir=`pwd`/tmp -jar /opt/gatk-4.0.0.0/gatk-package-4.0.0.0-local.jar \
     PrintReads \
     -R /datos/nagen/reference_genome/$Reference.fa \
     $(echo ${bam_array[@]}) \
     -O bam/${SampleID}.bam

printf END >&2; uptime >&2

trap : 0

echo >&2 '
************
*** DONE *** 
************
'
