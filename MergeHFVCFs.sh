#!/bin/bash
#$ -V
#$ -N MergeHFVCFs
#$ -m ea
#$ -cwd
#$ -o /datos/nagen/jobs/mergehfvcfs.stdout
#$ -e /datos/nagen/jobs/mergehfvcfs.stderr
#$ -w e
##$ -l h_vmem=2G
#$ -l m_core=4

# Merge hard filetered vcf files

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

# Get sampleID
SampleID=$1

# Merge VCF files
/opt/jdk1.8.0_161/bin/java -Xmx18000m -Djava.io.tmpdir=`pwd`/tmp -jar /opt/gatk-4.0.0.0/gatk-package-4.0.0.0-local.jar \
	   MergeVcfs \
	-R /datos/nagen/reference_genome/hs37d5.fa \
	-I $SampleID.snvs.labeled.vcf.gz \
	-I $SampleID.indels.labeled.vcf.gz \
	-O $SampleID.labeled.vcf.gz

printf END >&2; uptime >&2

trap : 0

echo >&2 '
************
*** DONE ***
************
'
