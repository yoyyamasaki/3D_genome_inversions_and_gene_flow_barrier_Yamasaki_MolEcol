#!/bin/bash
#SNP call by bcftools. SNP call was conducted for each chromosome in parallel, then they were merged.
#Dependency: bcftools v 1.9, vcftools v0.1.16

OUTPUT=Akkeshi_Quebec_WLD_bcftools.POmaskXIIXXIYv20240517.vcf.gz
DIR=bcftools_call_akkeshi

#SNP call in parallel
parallel -j 22 bash bcftools_call_array.sh ::: {1..22}

ls ./${DIR}/*_chr*.vcf.gz > vcf_list.txt

#MErge vcf files
vcf-concat -f vcf_list.txt |bgzip -c > ${DIR}/$OUTPUT

#Calculate depth for each site and individual
vcftools --gzvcf ${DIR}/$OUTPUT --site-mean-depth --out  ${DIR}/"$OUTPUT"_depth
vcftools --gzvcf ${DIR}/$OUTPUT --depth --out  ${DIR}/"$OUTPUT"_depth_ind

