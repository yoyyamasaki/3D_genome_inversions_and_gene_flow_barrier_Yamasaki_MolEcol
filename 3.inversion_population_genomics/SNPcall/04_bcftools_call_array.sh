#!/bin/bash
#SNP call by bcftools. SNP call was conducted for each chromosome in parallel, then they were merged.
#Dependency: bcftools v 1.9, vcftools v0.1.16

PREFIX=Akkeshi_Quebec_WLD
BAM_LIST=./bam_list_akkeshi.txt
REF=~/db/PO_merged_maskXII_XXI_Y_20240517/PO_male_merged.v20240517.maskXII_XXI_Y.fa
VERSION=POmaskXIIXXIYv20240517
DIR=bcftools_call_akkeshi

mkdir -p ${DIR}

CHR=(XXX $(grep ">" ${REF}|perl -pe 's/>//g'))

bcftools mpileup -b ${BAM_LIST} -r ${CHR[$1]} -f $REF -a FORMAT/AD,FORMAT/DP,FORMAT/SP,INFO/AD|\
bcftools call -f GQ -M -m -V indels|\
bgzip > ./${DIR}/${PREFIX}_${CHR[$1]}_bcftools.${VERSION}.vcf.gz

tabix ./${DIR}/${PREFIX}_${CHR[$1]}_bcftools.${VERSION}.vcf.gz
