#!/bin/bash
#$ -V
#$ -N gem3mapper
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

# run gem-mapper 
/opt/gem3/bin/gem-mapper --threads 20 --max-reported-matches 100 --report-file /datos/nagen/gem3_mapping/P2026_101.json -I /datos/nagen/reference_genome/GRCh38.p12.index.gem -1 trimmed_reads/P2026_101_1_trimmed.fastq.gz -2 trimmed_reads/P2026_101_2_trimmed.fastq.gz -o gem3_mapping/P2026_101.pe.sam

printf END >&2; uptime >&2
