#!/bin/bash
#$ -V
#$ -N HaplotypeCaller
#$ -m ea
#$ -cwd
#$ -o /datos/nagen/jobs/haplotypecaller.stdout
#$ -e /datos/nagen/jobs/haplotypecaller.stderr
#$ -w e
##$ -l h_vmem=2G
#$ -l m_core=10

# Calling SNPs and indels for chromosomes

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

# Get sample name, chromosome and reference genome
SampleID=$1
Chromosome=$2
Reference=$3


hostname >&2
echo $0 >&2
printf START >&2; uptime >&2
date >&2

# Call SNPs and indels simultaneously using the Haplotypecaller
# -ERC GVCF \

/opt/jdk1.8.0_161/bin/java -Xmx18000m -Djava.io.tmpdir=`pwd`/tmp -jar /opt/gatk-4.0.0.0/gatk-package-4.0.0.0-local.jar \
	    HaplotypeCaller \
	 -R /datos/nagen/reference_genome/$Reference.fa \
	 -I recalibration/$SampleID.sorted.rmdup.filtered.recalibrated.$Chromosome.bam \
     -O gvcf/$SampleID.$Chromosome.g.vcf.gz \
     -ERC GVCF \
     -L $Chromosome \
     --sample-name $SampleID
	 
printf END >&2; uptime >&2

trap : 0

echo >&2 '
************
*** DONE ***
************
'

