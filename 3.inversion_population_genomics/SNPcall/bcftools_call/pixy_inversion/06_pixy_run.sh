#!/bin/bash
#Run pixy to calculate Fst and dxy.
#Dependency: pixy v1.2.7.beta1

VCF_inINV=/home/yo/HDD/yo_analysis/HDD4_20230223/stickleback/variant_call/PO_merged_maskXII_XXI_Y_20240517/bcftools_call/PO_JS_Quebec_WLD_bcftools.inINV.filtered.POmaskXIIXXIYv20240517.vcf.gz

VCF_SYNAL=/home/yo/HDD/yo_analysis/HDD4_20230223/stickleback/variant_call/PO_merged_maskXII_XXI_Y_20240517/bcftools_call/PO_JS_Quebec_WLD_bcftools.SYNAL.filtered.POmaskXIIXXIYv20240517.vcf.gz

OUT_HEADER_inINV=PO_JS_Quebec_WLD_bcftools.inINV.filtered.POmaskXIIXXIYv20240517
OUT_HEADER_outINV=PO_JS_Quebec_WLD_bcftools.outINV.filtered.POmaskXIIXXIYv20240517
OUT_HEADER_SYNAL=PO_JS_Quebec_WLD_bcftools.SYNAL.filtered.POmaskXIIXXIYv20240517

OUT_HEADER_inINV_female=PO_JS_Quebec_WLD_bcftools.inINV.female.filtered.POmaskXIIXXIYv20240517
OUT_HEADER_outINV_female=PO_JS_Quebec_WLD_bcftools.outINV.female.filtered.POmaskXIIXXIYv20240517
OUT_HEADER_SYNAL_female=PO_JS_Quebec_WLD_bcftools.SYNAL.female.filtered.POmaskXIIXXIYv20240517


pixy --stats pi fst dxy --vcf ${VCF_inINV} --population popmap.txt --window_size 25000 --n_cores 20 --fst_type wc --output_prefix ${OUT_HEADER_inINV}

pixy --stats pi fst dxy --vcf ${VCF_SYNAL} --population popmap.txt --window_size 25000 --n_cores 20 --fst_type wc --output_prefix ${OUT_HEADER_SYNAL}

