#!/bin/bash
#$ -V
#$ -terse
#$ -m ea
#$ -M alberto.labarga.gutierrez@navarra.es
#$ -o /home/alabarga/jobs/hardfilters_fam.stdout
#$ -e /home/alabarga/jobs/hardfilters_fam.stderr
#$ -w e
#$ -l h_vmem=16G
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

familyId=$1
referenceGenome=/datos/nagen/reference_genome/hs37d5.fa

hostname >&2
echo $0 >&2
printf START >&2; uptime >&2
date >&2

## Apply hard filters.
# Extract SNPs from the call set 
/opt/gatk-4.0.0.0/gatk SelectVariants \
-R $referenceGenome \
-V $familyId'.vcf.gz' \
--select-type-to-include SNP \
-O $familyId'.snvs.vcf.gz'

echo >&2 '
************
*** DONE SelectVariants *** 
************
'
# Apply hard filters for SNPs
/opt/gatk-4.0.0.0/gatk VariantFiltration \
-R $referenceGenome \
-V $familyId'.snvs.vcf.gz' \
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
-O $familyId'.snvs.labeled.vcf.gz'

echo >&2 '
************
*** DONE VariantFiltration ***
************
'

# Extract Indels from the call set
/opt/gatk-4.0.0.0/gatk SelectVariants \
-R $referenceGenome \
-V $familyId'.vcf.gz' \
--select-type-to-include INDEL \
-O $familyId'.indels.vcf.gz'

echo >&2 '
************
*** DONE SelectVariants ***
************
'

# Apply hard filters for Indels
/opt/gatk-4.0.0.0/gatk VariantFiltration \
-R $referenceGenome \
-V $familyId'.indels.vcf.gz' \
-filter "QD < 2.0" \
--filter-name "QualByDepth" \
-filter "FS > 200.0" \
--filter-name "FisherStrand" \
-filter "ReadPosRankSum < -20.0" \
--filter-name "ReadPosRankSumTest" \
-filter "SOR > 10.0" \
--filter-name "StrandOddsRatio" \
-O $familyId'.indels.labeled.vcf.gz'

echo >&2 '
************
*** DONE VariantFiltration ***
************
'

# Combine SNVs and Indels
/opt/gatk-4.0.0.0/gatk MergeVcfs \
-R $referenceGenome \
-I $familyId'.snvs.labeled.vcf.gz' \
-I $familyId'.indels.labeled.vcf.gz' \
-O $familyId'.labeled.vcf.gz' \

echo >&2 '
************
*** DONE MergeVcfs *** 
************
'
printf END >&2; uptime >&2
