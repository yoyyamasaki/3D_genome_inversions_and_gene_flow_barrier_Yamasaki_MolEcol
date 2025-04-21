#!/bin/sh
#Filtering of the vcf file for pixy and dsuite.
#Dependency: vcftools v0.1.16

INPUT=PO_JS_Quebec_WLD_bcftools.POmaskXIIXXIYv20240517.vcf.gz
MAXDP=45.09087
MINQ=20
MINGQ=20
MINDP=10
KEEP=PO_JS_Quebec.txt
INVERSION=syri20240917_PO_JS_asm5syri.INV.bed
SYNTENY_BED=syri20240917_PO_JS_asm5syri.SYNAL.bed
INV_SYNAL=syri20240917_PO_JS_asm5syri.INV_SYNAL.bed
OUT1=PO_JS_Quebec_WLD_bcftools.inINV.filtered.POmaskXIIXXIYv20240517.vcf.gz
OUT2=PO_JS_Quebec_WLD_bcftools.dsuite.inINV.filtered.rmsexchr.POmaskXIIXXIYv20240517.vcf.gz
OUT3=PO_JS_Quebec_WLD_bcftools.SYNAL.filtered.POmaskXIIXXIYv20240517.vcf.gz
OUT4=PO_JS_Quebec_WLD_bcftools.dsuite.SYNAL.filtered.rmsexchr.POmaskXIIXXIYv20240517.vcf.gz
WLD=WLD_bcftools.dsuite.filtered.POmaskXIIXXIYv20240517.vcf.gz

#For dsuite: Extract site which available in outgroup (G. wheatlandi).
vcftools --gzvcf ${INPUT} --recode -c --recode-INFO-all --keep WLD.txt --max-meanDP ${MAXDP} --remove-indels --minQ ${MINQ} --remove-filtered-all|\
vcftools --vcf - --recode -c --minDP ${MINDP} --minGQ ${MINGQ} --out ${WLD}| \
vcftools --vcf - --recode -c --max-missing 1|bgzip > ${WLD}

#Within inversion
##For pixy
vcftools --gzvcf ${INPUT} --recode -c --recode-INFO-all --keep ${KEEP} --max-meanDP ${MAXDP} --remove-indels --bed ${INVERSION} --max-alleles 2 --minQ ${MINQ} --remove-filtered-all|\
vcftools --vcf - --recode -c --exclude-bed ${INV_SYNAL}| \
vcftools --vcf - --recode -c --minDP ${MINDP} --minGQ ${MINGQ} \
vcftools --vcf - --recode -c --max-missing 0.8|bgzip > ${OUT1}
tabix ${OUT1}

##For dsuite:remove all sex chromosomes
vcftools --gzvcf ${INPUT} --recode -c --recode-INFO-all --max-meanDP ${MAXDP} --remove-indels --min-alleles 2 --max-alleles 2 --mac 1 --minQ ${MINQ} --remove-filtered-all --positions ${WLD} --not-chr chrXIX --not-chr chrY --not-chr chrIX|\
vcftools --vcf - --recode -c --exclude-bed ${INV_SYNAL}| \
vcftools --vcf - --recode -c --minDP ${MINDP} --minGQ ${MINGQ} --bed ${INVERSION}| \
vcftools --vcf - --recode -c --max-missing 0.8|bgzip > ${OUT2}
tabix ${OUT2}

#syntenic regions
##For pixy
vcftools --gzvcf ${INPUT} --recode -c --recode-INFO-all --keep ${KEEP} --max-meanDP ${MAXDP} --remove-indels --bed ${SYNTENY_BED} --max-alleles 2 --minQ ${MINQ} --remove-filtered-all|\
vcftools --vcf - --recode -c --exclude-bed ${INV_SYNAL}| \
vcftools --vcf - --recode -c --minDP ${MINDP} --minGQ ${MINGQ}| \
vcftools --vcf - --recode -c --max-missing 0.8|bgzip > ${OUT3}
tabix ${OUT3}

##For dsuite:remove all sex chromosomes
vcftools --gzvcf ${INPUT} --recode -c --recode-INFO-all --max-meanDP ${MAXDP} --remove-indels --min-alleles 2 --max-alleles 2 --mac 1 --minQ ${MINQ} --remove-filtered-all --positions ${WLD} --not-chr chrXIX --not-chr chrY --not-chr chrIX|\
vcftools --vcf - --recode -c --exclude-bed ${INV_SYNAL}| \
vcftools --vcf - --recode -c --minDP ${MINDP} --minGQ ${MINGQ} --bed ${SYNTENY_BED}| \
vcftools --vcf - --recode -c --max-missing 0.8|bgzip > ${OUT4}
tabix ${OUT4}

