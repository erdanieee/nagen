#!/bin/bash
#$ -V
#$ -N VariantRecalibration
#$ -m ea
#$ -cwd
#$ -o /datos/nagen/jobs/hardfilters.stdout
#$ -e /datos/nagen/jobs/hardfilters.stderr
#$ -w e
##$ -l h_vmem=16G
#$ -l m_core=4

# Apply variant recalibration to a VCF file
# https://software.broadinstitute.org/gatk/download/bundle
# https://software.broadinstitute.org/gatk/documentation/article.php?id=2805
# https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.3.0/org_broadinstitute_hellbender_tools_walkers_vqsr_VariantRecalibrator.php
# https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.3.0/org_broadinstitute_hellbender_tools_walkers_vqsr_ApplyVQSR.php

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

PATH=$PATH:/opt/R/3.4.3/bin/

# Get sample ID
SAMPLE_ID=$1
REF_GENOME=/datos/nagen/reference_genome/hs37d5.fa
VARIANCE_REFERENCE=/datos/nagen/gatk_resources

# -an InbreedingCoeff \
# --maxGaussians 4

/opt/gatk-4.0.4.0/gatk VariantRecalibrator \
   -R ${REF_GENOME} \
   -V vcf/${SAMPLE_ID}.vcf.gz \
   --resource hapmap,known=false,training=true,truth=true,prior=15.0:${VARIANCE_REFERENCE}/hapmap_3.3.b37.vcf \
   --resource omni,known=false,training=true,truth=false,prior=12.0:${VARIANCE_REFERENCE}/1000G_omni2.5.b37.vcf \
   --resource 1000G,known=false,training=true,truth=false,prior=10.0:${VARIANCE_REFERENCE}/1000G_phase1.snps.high_confidence.b37.vcf \
   --resource dbsnp,known=true,training=false,truth=false,prior=2.0:${VARIANCE_REFERENCE}/dbsnp_138.b37.vcf \
    -an DP \
    -an QD \
    -an FS \
    -an SOR \
    -an MQ \
    -an MQRankSum \
    -an ReadPosRankSum \
   -mode SNP \
   -O recalibration/output.snp.recal.vcf \
   --tranches-file recalibration/output.snp.tranches \
   --rscript-file recalibration/output.snp.plots.R
   
/opt/gatk-4.0.4.0/gatk  ApplyVQSR \
   -R ${REF_GENOME} \
   -V vcf/${SAMPLE_ID}.vcf.gz \
   -O vcf/${SAMPLE_ID}.recalibrated.snp.vcf.gz \
   -ts-filter-level 99.0 \
   --tranches-file recalibration/output.snp.tranches \
   --recal-file recalibration/output.snp.recal.vcf \
   -mode SNP

/opt/gatk-4.0.4.0/gatk VariantRecalibrator \
   -R ${REF_GENOME} \
   -V vcf/${SAMPLE_ID}.recalibrated.snp.vcf.gz \
   --resource hapmap,known=false,training=true,truth=true,prior=15.0:${VARIANCE_REFERENCE}/hapmap_3.3.b37.vcf \
   --resource mills,known=false,training=true,truth=false,prior=12.0:${VARIANCE_REFERENCE}/Mills_and_1000G_gold_standard.indels.b37.vcf \
   --resource dbsnp,known=true,training=false,truth=false,prior=2.0:${VARIANCE_REFERENCE}/dbsnp_138.b37.vcf \
    -an QD \
    -an DP \
    -an FS \
    -an SOR \
    -an MQRankSum \
    -an ReadPosRankSum \
   -mode INDEL \
   -O recalibration/output.indel.recal.vcf \
   --tranches-file recalibration/output.indel.tranches \
   --rscript-file recalibration/output.indel.plots.R
   
/opt/gatk-4.0.4.0/gatk  ApplyVQSR \
   -R ${REF_GENOME} \
   -V vcf/${SAMPLE_ID}.recalibrated.snp.vcf.gz \
   -O vcf/${SAMPLE_ID}.recalibrated.vcf.gz \
   -ts-filter-level 99.0 \
   --tranches-file recalibration/output.indel.tranches \
   --recal-file recalibration/output.indel.recal.vcf \
   -mode INDEL

echo >&2 '
************
*** DONE VariantRecalibration *** 
************
'
printf END >&2; uptime >&2
