#!/bin/bash
#$ -V
#$ -N hardfilters
#$ -m ea
#$ -cwd
#$ -o /datos/nagen/jobs/hardfilters.stdout
#$ -e /datos/nagen/jobs/hardfilters.stderr
#$ -w e
##$ -l h_vmem=16G
#$ -l m_core=4

# Apply hard filters to a VCF file

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

hostname >&2
echo $0 >&2
printf START >&2; uptime >&2
date >&2

## Apply hard filters.
# Extract SNPs from the call set 
/opt/gatk-4.0.0.0/gatk SelectVariants \
	-R /datos/nagen/reference_genome/hs37d5.fa \
	-V vcf/${SAMPLE_ID}.vcf.gz \
	--select-type-to-include SNP \
	-O vcf/${SAMPLE_ID}.snvs.vcf.gz

echo >&2 '
************
*** DONE SelectVariants *** 
************
'
# Apply hard filters for SNPs
# Following https://gatkforums.broadinstitute.org/gatk/
# discussion/2806/howto-apply-hard-filters-to-a-call-set)
/opt/gatk-4.0.0.0/gatk VariantFiltration \
	-R /datos/nagen/reference_genome/hs37d5.fa \
	-V vcf/${SAMPLE_ID}.snvs.vcf.gz \
	-filter "QD < 2.0" \
	--filter-name "QualByDepth" \
	-filter "FS > 60.0" \
	--filter-name "FisherStrand" \
	-filter "MQ < 40.0" \
	--filter-name "RMSMappingQuality" \
	-filter "MQRankSum < -12.5" \
	--filter-name "MappingQualityRankSumTest" \
	-filter "ReadPosRankSum < -8.0" \
	--filter-name "ReadPosRankSumTest" \
	-filter "SOR > 3.0" \
	--filter-name "StrandOddsRatio" \
	-O vcf/${SAMPLE_ID}.snvs.labeled.vcf.gz

echo >&2 '
************
*** DONE VariantFiltration ***
************
'

# Extract Indels from the call set
/opt/gatk-4.0.0.0/gatk SelectVariants \
	-R /datos/nagen/reference_genome/hs37d5.fa \
	-V vcf/${SAMPLE_ID}.vcf.gz \
	--select-type-to-include INDEL \
	-O vcf/${SAMPLE_ID}.indels.vcf.gz

echo >&2 '
************
*** DONE SelectVariants ***
************
'

# Apply hard filters for Indels
# See link above
/opt/gatk-4.0.0.0/gatk VariantFiltration \
	-R /datos/nagen/reference_genome/hs37d5.fa \
	-V vcf/${SAMPLE_ID}.indels.vcf.gz \
	-filter "QD < 2.0" \
	--filter-name "QualByDepth" \
	-filter "FS > 200.0" \
	--filter-name "FisherStrand" \
	-filter "ReadPosRankSum < -20.0" \
	--filter-name "ReadPosRankSumTest" \
	-filter "SOR > 10.0" \
	--filter-name "StrandOddsRatio" \
	-O vcf/${SAMPLE_ID}.indels.labeled.vcf.gz

echo >&2 '
************
*** DONE VariantFiltration ***
************
'

# Combine SNVs and Indels
/opt/gatk-4.0.0.0/gatk MergeVcfs \
	-R /datos/nagen/reference_genome/hs37d5.fa \
	-I vcf/${SAMPLE_ID}.snvs.labeled.vcf.gz \
	-I vcf/${SAMPLE_ID}.indels.labeled.vcf.gz \
	-O vcf/${SAMPLE_ID}.labeled.vcf.gz \

echo >&2 '
************
*** DONE MergeVcfs *** 
************
'
printf END >&2; uptime >&2
