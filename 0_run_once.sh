#!/bin/bash
#$ -V
#$ -N HaplotypeCaller
#$ -m ea
#$ -o /datos/nagen/jobs/haplotypecaller.stdout
#$ -e /datos/nagen/jobs/haplotypecaller.stderr
#$ -w e
##$ -l h_vmem=2G
#$ -l m_core=20

###
/opt/jdk1.8.0_161/bin/java -jar -Xmx2g -Djava.io.tmpdir=`pwd`/tmp -jar /opt/picard/build/libs/picard.jar AddOrReplaceReadGroups I=bam/AC5445.bam ID=AC5445 SM=AC5445 PL=ILLUMINA LB=Fragment PU=AC5445 O=bam/bamforHaplotypeCaller.bam
###

/opt/samtools-1.6/bin/samtools index /datos/nagen/nagen_samples/AC5445/bam/bamforHaplotypeCaller.bam

/opt/jdk1.8.0_161/bin/java -Xmx4g -Djava.io.tmpdir=`pwd`/tmp -jar /opt/gatk-4.0.0.0/gatk-package-4.0.0.0-local.jar \
        HaplotypeCaller \
     -R /datos/nagen/reference_genome/hs37d5.fa \
     -I bam/bamforHaplotypeCaller.bam \
     -ERC GVCF \
     --sample AC5445 \
     --output gcvf/AC5445.g.vcf.gz