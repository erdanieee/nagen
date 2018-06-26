#!/bin/bash
#$ -V
#$ -N GenotypeGVCF
#$ -m ea
#$ -M rodrigo.bacigalupe.perez@navarra.es
#$ -cwd
#$ -o /datos/nagen/jobs/genotype.stdout
#$ -e /datos/nagen/jobs/genotype.stderr
#$ -w e
##$ -l h_vmem=2G
#$ -l m_core=1

# Genotype GVCFs and produce VCF files

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
SAMPLE_ID=$1
chr=$2

hostname >&2
echo $0 >&2
printf START >&2; uptime >&2
date >&2

# GenotypeGVCFs
/opt/jdk1.8.0_161/bin/java -Djava.io.tmpdir=`pwd`/tmp -jar /opt/gatk-4.0.0.0/gatk-package-4.0.0.0-local.jar \
		   GenotypeGVCFs \
		-R /datos/nagen/reference_genome/hs37d5.fa \
		-V gvcf/${SAMPLE_ID}.${chr}.g.vcf.gz \
		-L ${chr} \
		-O vcf/${SAMPLE_ID}.${chr}.vcf.gz

printf END >&2; uptime >&2

trap : 0

echo >&2 '
************
*** DONE ***
************
'
