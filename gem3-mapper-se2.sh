#!/bin/bash
#$ -V
#$ -m ea
#$ -M rodrigo.bacigalupe.perez@navarra.es
#$ -wd /datos/nagen/
#$ -o /datos/nagen/jobs/gem3mapper_se2.stdout
#$ -e /datos/nagen/jobs/gem3mapper_se2.stderr
#$ -w e
##$ -l h_vmem=2G
#$ -l m_core=20

hostname >&2
printf START >&2; uptime >&2
date >&2

GENOME=$1
REFERENCE=$2

# run gem-mapper for single end 2
/opt/gem3/bin/gem-mapper --threads 20 --max-reported-matches 100 --report-file $GENOME'_se2'.json -I /datos/nagen/reference_genome/$REFERENCE.fna.index.gem -i trimmed_reads/$GENOME'_U2_trimmed'.fastq.gz -o gem3_mapping/$GENOME.se2.sam

printf END >&2; uptime >&2
