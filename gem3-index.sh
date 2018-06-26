#!/bin/bash
#$ -V
#$ -N gem3index
#$ -m ea
#$ -wd /datos/nagen/reference_genome
#$ -o /datos/nagen/jobs/gem3index.stdout
#$ -e /datos/nagen/jobs/gem3index.stderr
#$ -w e
##$ -l h_vmem=2G
#$ -l m_core=4

# Get reference genome
Reference=$1

# Index reference genome with gem-indexer
/opt/gem3/bin/gem-indexer -i $Reference.fna -o $Reference.fna.index
