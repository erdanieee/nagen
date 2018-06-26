#!/bin/bash
#$ -V
#$ -N alignproc
#$ -m ea
#$ -wd /datos/nagen/
#$ -o /datos/nagen/jobs/alignproc.stdout
#$ -e /datos/nagen/jobs/alignproc.stderr
#$ -w e
##$ -l h_vmem=2G
#$ -l m_core=20

hostname >&2
printf START >&2; uptime >&2
date >&2

# Take sam file as input
sam=$1

# Convert sam to bam
/opt/samtools-1.6/bin/samtools view -Sb $sam > "${sam%.sam}.bam"

printf END >&2; uptime >&2
