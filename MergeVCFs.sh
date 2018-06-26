#!/bin/bash
#$ -V
#$ -N MergeVCFs
#$ -m ea
#$ -cwd
#$ -o /datos/nagen/jobs/mergevcfs.stdout
#$ -e /datos/nagen/jobs/mergevcfs.stderr
#$ -w e
##$ -l h_vmem=2G
#$ -l m_core=4

# Merge individual vcf files

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

# Get sample ID
SAMPLE_ID=$1

# Merge vcf files
/opt/jdk1.8.0_161/bin/java -Xmx18000m -Djava.io.tmpdir=`pwd`/tmp -jar /opt/gatk-4.0.0.0/gatk-package-4.0.0.0-local.jar \
	   MergeVcfs \
	-R /datos/nagen/reference_genome/hs37d5.fa \
	-I vcf/${SAMPLE_ID}.1.vcf.gz \
	-I vcf/${SAMPLE_ID}.2.vcf.gz \
	-I vcf/${SAMPLE_ID}.3.vcf.gz \
	-I vcf/${SAMPLE_ID}.4.vcf.gz \
	-I vcf/${SAMPLE_ID}.5.vcf.gz \
	-I vcf/${SAMPLE_ID}.6.vcf.gz \
	-I vcf/${SAMPLE_ID}.7.vcf.gz \
	-I vcf/${SAMPLE_ID}.8.vcf.gz \
	-I vcf/${SAMPLE_ID}.9.vcf.gz \
	-I vcf/${SAMPLE_ID}.10.vcf.gz \
	-I vcf/${SAMPLE_ID}.11.vcf.gz \
	-I vcf/${SAMPLE_ID}.12.vcf.gz \
	-I vcf/${SAMPLE_ID}.13.vcf.gz \
	-I vcf/${SAMPLE_ID}.14.vcf.gz \
	-I vcf/${SAMPLE_ID}.15.vcf.gz \
	-I vcf/${SAMPLE_ID}.16.vcf.gz \
	-I vcf/${SAMPLE_ID}.17.vcf.gz \
	-I vcf/${SAMPLE_ID}.18.vcf.gz \
	-I vcf/${SAMPLE_ID}.19.vcf.gz \
	-I vcf/${SAMPLE_ID}.20.vcf.gz \
	-I vcf/${SAMPLE_ID}.21.vcf.gz \
	-I vcf/${SAMPLE_ID}.22.vcf.gz \
	-I vcf/${SAMPLE_ID}.X.vcf.gz \
	-I vcf/${SAMPLE_ID}.Y.vcf.gz \
	-O vcf/${SAMPLE_ID}.vcf.gz

# Merge gvcf files
/opt/jdk1.8.0_161/bin/java -Xmx18000m -Djava.io.tmpdir=`pwd`/tmp -jar /opt/gatk-4.0.0.0/gatk-package-4.0.0.0-local.jar \
	   MergeVcfs \
	-R /datos/nagen/reference_genome/hs37d5.fa \
	-I gvcf/${SAMPLE_ID}.1.g.vcf.gz \
	-I gvcf/${SAMPLE_ID}.2.g.vcf.gz \
	-I gvcf/${SAMPLE_ID}.3.g.vcf.gz \
	-I gvcf/${SAMPLE_ID}.4.g.vcf.gz \
	-I gvcf/${SAMPLE_ID}.5.g.vcf.gz \
	-I gvcf/${SAMPLE_ID}.6.g.vcf.gz \
	-I gvcf/${SAMPLE_ID}.7.g.vcf.gz \
	-I gvcf/${SAMPLE_ID}.8.g.vcf.gz \
	-I gvcf/${SAMPLE_ID}.9.g.vcf.gz \
	-I gvcf/${SAMPLE_ID}.10.g.vcf.gz \
	-I gvcf/${SAMPLE_ID}.11.g.vcf.gz \
	-I gvcf/${SAMPLE_ID}.12.g.vcf.gz \
	-I gvcf/${SAMPLE_ID}.13.g.vcf.gz \
	-I gvcf/${SAMPLE_ID}.14.g.vcf.gz \
	-I gvcf/${SAMPLE_ID}.15.g.vcf.gz \
	-I gvcf/${SAMPLE_ID}.16.g.vcf.gz \
	-I gvcf/${SAMPLE_ID}.17.g.vcf.gz \
	-I gvcf/${SAMPLE_ID}.18.g.vcf.gz \
	-I gvcf/${SAMPLE_ID}.19.g.vcf.gz \
	-I gvcf/${SAMPLE_ID}.20.g.vcf.gz \
	-I gvcf/${SAMPLE_ID}.21.g.vcf.gz \
	-I gvcf/${SAMPLE_ID}.22.g.vcf.gz \
	-I gvcf/${SAMPLE_ID}.X.g.vcf.gz \
	-I gvcf/${SAMPLE_ID}.Y.g.vcf.gz \
	-O gvcf/${SAMPLE_ID}.g.vcf.gz

printf END >&2; uptime >&2

trap : 0

echo >&2 '
************
*** DONE *** 
************
'
