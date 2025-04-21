#!/bin/bash
#Conduct D statistics analysis by dsuite.
#Dependency: dsuite v0.5 r52

DSUITE=~/software/Dsuitev05
VCF_INV=/home/yo/HDD/yo_analysis/HDD4_20230223/stickleback/variant_call/PO_merged_maskXII_XXI_Y_20240517/bcftools_call/PO_JS_Quebec_WLD_bcftools.dsuite.inINV.filtered.rmsexchr.POmaskXIIXXIYv20240517.vcf.gz
VCF_OUTINV=/home/yo/HDD/yo_analysis/HDD4_20230223/stickleback/variant_call/PO_merged_maskXII_XXI_Y_20240517/bcftools_call/PO_JS_Quebec_WLD_bcftools.dsuite.outINV.filtered.rmsexchr.POmaskXIIXXIYv20240517.vcf.gz
VCF_SYNAL=/home/yo/HDD/yo_analysis/HDD4_20230223/stickleback/variant_call/PO_merged_maskXII_XXI_Y_20240517/bcftools_call/PO_JS_Quebec_WLD_bcftools.dsuite.SYNAL.filtered.rmsexchr.POmaskXIIXXIYv20240517.vcf.gz
POPMAP=popmap_all.txt
OUTPREFIX_INV=dsuite_tree_Quebec_PO_JS_WLD.INV
OUTPREFIX_SYNAL=dsuite_tree_Quebec_PO_JS_WLD.SYNAL
TREE=tree_dsuite_tree_Quebec_PO_JS_WLD.tre

${DSUITE} Dtrios -t ${TREE} --o ${OUTPREFIX_INV} ${VCF_INV} ${POPMAP}

${DSUITE} Dtrios -t ${TREE} --o ${OUTPREFIX_SYNAL} ${VCF_SYNAL} ${POPMAP}

