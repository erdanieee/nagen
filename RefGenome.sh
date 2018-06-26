#!/bin/bash
#$ -V
#$ -m ea
#$ -wd /datos/nagen/reference_genome
#$ -o /datos/nagen/jobs/gem3index.stdout
#$ -e /datos/nagen/jobs/gem3index.stderr
#$ -w e
##$ -l h_vmem=2G
#$ -l m_core=4

# Prepare reference genome

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

hostname >&2
echo $0 >&2
printf START >&2; uptime >&2
date >&2

# Download genome from the internet
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
gunzip hs37d5.fa.gz

# Generate BWA index for mapping using BWA
/opt/bwa/bwa index -a bwtsw hs37d5.fa

# Generate the fasta file index
opt/samtools-1.6/bin/samtools faidx hs37d5.fa

# Create dict file
/opt/jdk1.8.0_161/bin/java -jar /opt/picard/build/libs/picard.jar CreateSequenceDictionary REFERENCE=hs37d5.fa OUTPUT=hs37d5.dict

printf END >&2; uptime >&2

trap : 0

echo >&2 '
************
*** DONE *** 
************
'
