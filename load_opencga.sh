#!/bin/bash
#$ -V
#$ -terse
#$ -m ea
#$ -M alberto.labarga.gutierrez@navarra.es
#$ -wd /home/alabarga/opencga_tutorial
#$ -o /home/alabarga/jobs/opencga.stdout
#$ -e /home/alabarga/jobs/opencga.stderr
#$ -w e
##$ -l h_vmem=16G
##$ -l m_core=4

# Load to OpenCGA

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

fileName=$1
outDir=$2 

hostname >&2
echo $0 >&2
printf START >&2; uptime >&2
date >&2

/opt/opencga/bin/opencga.sh  variant index --file $fileName  --calculate-stats --annotate -o $outDir

printf END >&2; uptime >&2
