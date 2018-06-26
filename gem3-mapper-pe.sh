#!/bin/bash
#$ -V
#$ -m ea
#$ -M rodrigo.bacigalupe.perez@navarra.es
#$ -wd /datos/nagen/
#$ -o /datos/nagen/jobs/gem3mapper.stdout
#$ -e /datos/nagen/jobs/gem3mapper.stderr
#$ -w e
##$ -l h_vmem=2G
#$ -l m_core=20

hostname >&2
printf START >&2; uptime >&2
date >&2

GENOME=$1
REFERENCE=$2

# run gem-mapper
/opt/gem3/bin/gem-mapper --threads 20 --max-reported-matches 100 \
--report-file P2026_101.json -I /datos/nagen/reference_genome/$REFERENCE.fna.index.gem -1 trimmed_reads/$GENOME'_1_trimmed'.fastq.gz -2 trimmed_reads/$GENOME'_2_trimmed'.fastq.gz -o gem3_mapping/$GENOME.pe.sam

printf END >&2; uptime >&2
