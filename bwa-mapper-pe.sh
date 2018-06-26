#!/bin/bash
#$ -V
#$ -m ea
#$ -cwd
#$ -o /datos/nagen/jobs/bwape.stdout
#$ -e /datos/nagen/jobs/bwape.stderr
#$ -w e
##$ -l h_vmem=2G
#$ -l m_core=20

hostname >&2
printf START >&2; uptime >&2
date >&2

# Read input files (genome name for sample and reference genome used)
GENOME=$1
REFERENCE=$2

# Run BWA on paired-end sequencing reads
/usr/local/bin/bwa mem -t 20 /datos/nagen/reference_genome/$REFERENCE /datos/nagen/P2026_101/trimmed_reads/$GENOME'_1_trimmed.fastq.gz' /datos/nagen/P2026_101/trimmed_reads/$GENOME'_2_trimmed.fastq.gz' -R '@RG\tID:'$GENOME'\tSM:'$GENOME'\tPL:Illumina\tCN:CBRA\tLB:Fragment' > $GENOME'.sam' # Paired-end

printf END >&2; uptime >&2
