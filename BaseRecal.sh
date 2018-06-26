#!/bin/bash
#$ -V
#$ -N BaseRecal
#$ -m ea
#$ -cwd
#$ -o /datos/nagen/jobs/baserecal.stdout
#$ -e /datos/nagen/jobs/baserecal.stderr
#$ -w e
##$ -l h_vmem=2G
#$ -l m_core=4

# Recalibration of base quality scores

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

# Get sample name and chromosome number
SampleID=$1
Chromosome=$2
Reference=$3

hostname >&2
echo $0 >&2
printf START >&2; uptime >&2
date >&2

# Run Recalibration of base quality scores previous to GATK calling
/opt/jdk1.8.0_161/bin/java -Xmx18000m -Djava.io.tmpdir=`pwd`/tmp -jar /opt/gatk-4.0.0.0/gatk-package-4.0.0.0-local.jar \
	    BaseRecalibrator \
	 -R /datos/nagen/reference_genome/$Reference.fa \
	 -I bam/$SampleID.sorted.rmdup.filtered.bam \
	 --use-original-qualities \
	 --known-sites /datos/nagen/gatk_resources/dbsnp_138.b37.vcf \
	 --known-sites /datos/nagen/gatk_resources/Mills_and_1000G_gold_standard.indels.b37.vcf \
	 --known-sites /datos/nagen/gatk_resources/1000G_phase1.indels.b37.vcf \
	 -L $Chromosome \
	 --output recalibration/$SampleID.sorted.rmdup.filtered.$Chromosome.recal_data.table

# Apply Recalibration
/opt/jdk1.8.0_161/bin/java -Xmx18000m -Djava.io.tmpdir=`pwd`/tmp -jar /opt/gatk-4.0.0.0/gatk-package-4.0.0.0-local.jar \
	    ApplyBQSR \
	 -R /datos/nagen/reference_genome/$Reference.fa \
	 -I bam/$SampleID.sorted.rmdup.filtered.bam \
	 -bqsr recalibration/$SampleID.sorted.rmdup.filtered.$Chromosome.recal_data.table \
	 --use-original-qualities \
	 -L $Chromosome \
	 --output recalibration/$SampleID.sorted.rmdup.filtered.recalibrated.$Chromosome.bam

printf END >&2; uptime >&2

trap : 0

echo >&2 '
************
*** DONE *** 
************
'
